#!/bin/bash
#
# Dark Gene Pipeline - Sample Analysis Script
# Docker 이미지를 사용하여 샘플 분석 실행
#

set -e

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

pipeline_timestamp() { date '+%Y-%m-%d %H:%M:%S %z'; }

# 사용법 출력
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    -w, --work-dir YYMM        Work directory (e.g., 2601) [required]
    -s, --sample SAMPLE        Sample name [required]
    -d, --data-dir DIR         Data directory path (default: ./data)
    -r, --ref-dir DIR          Reference directory path (default: ./refs)
    -c, --cleanup              Clean up work directory after completion
    --skip-cnv                 Skip CNV analysis (single sample 시 필요, cohort mode는 2+ 샘플 필요)
    --aligner ALIGNER          Aligner to use: bwa-mem (default) or bwa-mem2
    --variant-caller CALLER    Variant caller: gatk (default), deepvariant, or strelka2
    --skip-vep                 Skip VEP annotation (use legacy snpEff mode)
    --shared-ref-dir DIR       Shared reference root (default: /data/reference)
    --fresh                    Delete sample work/ and .nextflow cache, then run (no -resume)
    -h, --help                 Show this help message

Environment:
    CHOWN_SPEC                 Numeric uid:gid for output ownership (default: \$(id -u):\$(id -g); e.g. ken:genolyx on host)
                               Task containers use NXF_DOCKER_TASK_USER. After each run, the whole order dir under analysis/output/log
                               gets chown + group write (g+rwX) + setgid on directories so genolyx (or CHOWN_SPEC gid) can add samples without 755-only-owner issues.

Example:
    $0 -w 2601 -s Sample_A10
    $0 --work-dir 2601 --sample Sample_A10 --cleanup
    $0 -w 2601 -s Sample_A10 --aligner bwa-mem2 --variant-caller deepvariant

Directory Structure:
    <data-dir>/
    ├── fastq/<work-dir>/<sample>/     # Input FASTQ files
    ├── analysis/<work-dir>/<sample>/  # Intermediate files
    ├── output/<work-dir>/<sample>/    # Final results
    └── log/<work-dir>/<sample>/       # Logs

EOF
    exit 1
}

# 파라미터 초기화
WORK_DIR=""
SAMPLE_NAME=""
DATA_DIR="$(pwd)/data"
REF_DIR="$(pwd)/refs"
CLEANUP=""
SKIP_CNV=""
ALIGNER="bwa-mem"
VARIANT_CALLER="gatk"
SKIP_VEP="true"
SHARED_REF_DIR="/data/reference"
FRESH=""

# 파라미터 파싱
while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -s|--sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -r|--ref-dir)
            REF_DIR="$2"
            shift 2
            ;;
        -c|--cleanup)
            CLEANUP="--cleanup"
            shift
            ;;
        --skip-cnv)
            SKIP_CNV="--skip_cnv"
            shift
            ;;
        --aligner)
            ALIGNER="$2"
            shift 2
            ;;
        --variant-caller)
            VARIANT_CALLER="$2"
            shift 2
            ;;
        --skip-vep)
            SKIP_VEP="true"
            shift
            ;;
        --no-skip-vep)
            SKIP_VEP="false"
            shift
            ;;
        --shared-ref-dir)
            SHARED_REF_DIR="$2"
            shift 2
            ;;
        --fresh)
            FRESH="1"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo -e "${RED}Error: Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# 필수 파라미터 확인
if [ -z "$WORK_DIR" ] || [ -z "$SAMPLE_NAME" ]; then
    echo -e "${RED}Error: Work directory and sample name are required${NC}"
    usage
fi

# 절대 경로로 변환
DATA_DIR=$(realpath "$DATA_DIR")
REF_DIR=$(realpath "$REF_DIR")

# 디렉토리 존재 확인
FASTQ_DIR="${DATA_DIR}/fastq/${WORK_DIR}/${SAMPLE_NAME}"
if [ ! -d "${DATA_DIR}/data/refs" ] || [ ! -d "${DATA_DIR}/data/bed" ]; then
    echo -e "${RED}Error: data/refs and data/bed required. Use -d <project_root> (e.g. -d .)${NC}"
    exit 1
fi
if [ ! -d "$FASTQ_DIR" ]; then
    echo -e "${RED}Error: FASTQ directory not found: ${FASTQ_DIR}${NC}"
    exit 1
fi

# R1/R2 파일 확인 (*_R1_*, *_R1.*, *_1.*, *_R2_* 등)
R1_COUNT=$(find "$FASTQ_DIR" \( -name "*_R1_*" -o -name "*_R1.*" -o -name "*_1.fq.gz" -o -name "*_1.fastq.gz" \) | wc -l)
R2_COUNT=$(find "$FASTQ_DIR" \( -name "*_R2_*" -o -name "*_R2.*" -o -name "*_2.fq.gz" -o -name "*_2.fastq.gz" \) | wc -l)

if [ "$R1_COUNT" -eq 0 ] || [ "$R2_COUNT" -eq 0 ]; then
    echo -e "${RED}Error: R1/R2 FASTQ files not found in ${FASTQ_DIR}${NC}"
    exit 1
fi

# 결과물 소유자: 기본은 실행 사용자(uid:gid)
CHOWN_SPEC="${CHOWN_SPEC:-$(id -u):$(id -g)}"

# analysis|output|log/<WORK_DIR> 전체: chown CHOWN_SPEC + 그룹 rwx + 디렉터리 setgid (동일 gid가 새 항목에 상속)
# $1 = quiet 일 때 진행 메시지 생략 (ensure_sample_output_dirs 직후 등)
repair_order_tree_permissions() {
    local quiet="${1:-}" base t any=0
    for base in analysis output log; do
        t="${DATA_DIR}/${base}/${WORK_DIR}"
        [ -d "$t" ] || continue
        any=1
    done
    [ "$any" -eq 1 ] || return 0

    if command -v docker >/dev/null 2>&1 && id -nG | grep -qw docker && docker info &>/dev/null; then
        if docker run --rm --platform linux/amd64 \
            --name "gx-exome-perm-${WORK_DIR}" \
            -e WORK_DIR="$WORK_DIR" \
            -e CHOWN_SPEC="$CHOWN_SPEC" \
            -v "${DATA_DIR}/analysis:/fa" \
            -v "${DATA_DIR}/output:/fo" \
            -v "${DATA_DIR}/log:/fl" \
            alpine:3.20 \
            sh -c 'for base in /fa /fo /fl; do
                t="$base/$WORK_DIR"
                [ -d "$t" ] || continue
                chown -R "$CHOWN_SPEC" "$t"
                chmod -R g+rwX "$t"
                find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
            done'; then
            [ "$quiet" = quiet ] || echo "  Order trees (${WORK_DIR}): ${CHOWN_SPEC}, g+rwX, setgid on dirs (docker)"
            return 0
        fi
    fi
    if sudo -n true 2>/dev/null; then
        for base in analysis output log; do
            t="${DATA_DIR}/${base}/${WORK_DIR}"
            [ -d "$t" ] || continue
            sudo chown -R "${CHOWN_SPEC}" "$t" 2>/dev/null || true
            sudo chmod -R g+rwX "$t" 2>/dev/null || true
            sudo find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
        done
        [ "$quiet" = quiet ] || echo "  Order trees (${WORK_DIR}): ${CHOWN_SPEC}, g+rwX, setgid on dirs (sudo)"
        return 0
    fi
    for base in analysis output log; do
        t="${DATA_DIR}/${base}/${WORK_DIR}"
        [ -d "$t" ] || continue
        chmod -R g+rwX "$t" 2>/dev/null || true
        find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true
    done
    [ "$quiet" = quiet ] || echo -e "${YELLOW}  Order trees: collaborative chmod only (for chown use docker group or passwordless sudo).${NC}"
}

# 출력 디렉토리 생성. 과거 Docker가 analysis/<work>/ 를 root로 만들면 일반 사용자 mkdir 불가 → docker로 보정
ensure_sample_output_dirs() {
    local a o l
    a="${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}"
    o="${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}"
    l="${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}"
    if mkdir -p "${a}/.nextflow" "$o" "$l" 2>/dev/null; then
        repair_order_tree_permissions quiet
        return 0
    fi
    if ! command -v docker >/dev/null 2>&1 || ! docker info &>/dev/null; then
        echo -e "${RED}Error: cannot create sample directories (permission denied).${NC}"
        echo "  e.g. sudo chown -R \"\$(id -u):\$(id -g)\" \"${DATA_DIR}/analysis/${WORK_DIR}\" \"${DATA_DIR}/output/${WORK_DIR}\" \"${DATA_DIR}/log/${WORK_DIR}\""
        exit 1
    fi
    echo -e "${YELLOW}Creating output dirs via docker (host mkdir failed, e.g. parent dir owned by root)...${NC}"
    docker run --rm --platform linux/amd64 \
        --name "gx-exome-mkdir-${WORK_DIR}" \
        -e WORK_DIR="$WORK_DIR" \
        -e SAMPLE_NAME="$SAMPLE_NAME" \
        -e CHOWN_SPEC="$CHOWN_SPEC" \
        -v "${DATA_DIR}/analysis:/fa" \
        -v "${DATA_DIR}/output:/fo" \
        -v "${DATA_DIR}/log:/fl" \
        alpine:3.20 \
        sh -c 'mkdir -p "/fa/${WORK_DIR}/${SAMPLE_NAME}/.nextflow" "/fo/${WORK_DIR}/${SAMPLE_NAME}" "/fl/${WORK_DIR}/${SAMPLE_NAME}" && \
            for base in /fa /fo /fl; do \
                t="$base/${WORK_DIR}"; \
                [ -d "$t" ] || continue; \
                chown -R "$CHOWN_SPEC" "$t"; \
                chmod -R g+rwX "$t"; \
                find "$t" -type d -exec chmod g+s {} + 2>/dev/null || true; \
            done'
}
ensure_sample_output_dirs

echo "======================================"
echo "Carrier Screening Pipeline - Sample Analysis"
echo "======================================"
echo ""
echo -e "${GREEN}Configuration:${NC}"
echo "  Work Directory: ${WORK_DIR}"
echo "  Sample Name: ${SAMPLE_NAME}"
echo "  FASTQ Directory: ${FASTQ_DIR}"
echo "  Data Directory: ${DATA_DIR}"
echo "  Reference Directory: ${REF_DIR}"
echo "  Cleanup: ${CLEANUP:-disabled}"
echo "  Skip CNV: $([ -n "$SKIP_CNV" ] && echo "enabled" || echo "disabled")"
echo "  Aligner: ${ALIGNER}"
echo "  Variant Caller: ${VARIANT_CALLER}"
echo "  VEP Annotation: $([ "$SKIP_VEP" = "true" ] && echo "disabled (snpEff)" || echo "enabled")"
echo "  Shared Ref Dir: ${SHARED_REF_DIR}"
echo "  Fresh run: $([ -n "$FRESH" ] && echo "yes (work + .nextflow cleared, no -resume)" || echo "no")"
echo "  File ownership (CHOWN_SPEC): ${CHOWN_SPEC}"
echo ""

ANALYSIS_SAMPLE_DIR="${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}"
if [ -n "$FRESH" ]; then
    echo -e "${YELLOW}--fresh: removing Nextflow work and session cache...${NC}"
    # Docker often leaves work/ as root:root on the host mount — plain rm fails
    rm -rf "${ANALYSIS_SAMPLE_DIR}/work" "${ANALYSIS_SAMPLE_DIR}/.nextflow" 2>/dev/null || true
    if [ -e "${ANALYSIS_SAMPLE_DIR}/work" ] || [ -e "${ANALYSIS_SAMPLE_DIR}/.nextflow" ]; then
        echo -e "${YELLOW}  Remaining paths not removable as current user (Docker root); trying sudo...${NC}"
        if ! sudo rm -rf "${ANALYSIS_SAMPLE_DIR}/work" "${ANALYSIS_SAMPLE_DIR}/.nextflow"; then
            echo -e "${RED}Error: could not remove work/.nextflow. Run:${NC}"
            echo "  sudo rm -rf \"${ANALYSIS_SAMPLE_DIR}/work\" \"${ANALYSIS_SAMPLE_DIR}/.nextflow\""
            exit 1
        fi
    fi
    echo "  Cleared: ${ANALYSIS_SAMPLE_DIR}/work"
    echo "  Cleared: ${ANALYSIS_SAMPLE_DIR}/.nextflow"
    echo ""
fi

NEXTFLOW_RESUME="-resume"
if [ -n "$FRESH" ]; then
    NEXTFLOW_RESUME=""
fi

# Paraphase/SMAca/ExpansionHunter 등 biocontainer 전용 모듈은 variant caller 종류와
# 무관하게 항상 Docker 컨테이너가 필요하므로 -profile docker 를 항상 활성화한다.
NXF_PROFILE="-profile docker"

# Docker 이미지 확인
if ! docker images | grep -q "gx-exome"; then
    echo -e "${RED}Error: Docker image 'gx-exome' not found${NC}"
    echo "Please build the image first:"
    echo "  docker-compose -f docker/docker-compose.yml build"
    exit 1
fi

# 분석 시작
echo -e "${YELLOW}Starting analysis...${NC}"
echo "  Started at: $(pipeline_timestamp)"
echo ""

# Docker 컨테이너 실행 (태스크는 NXF_DOCKER_TASK_USER, 종료 시 docker chown)
# nextflow.config가 projectDir/../data/refs, data/bed 경로 사용 → /app/data 마운트 필요
# /data/reference 마운트: data/refs 내 심볼릭 링크 대상 경로 접근을 위해 필요
# HOST_WORK_DIR: Nextflow이 태스크 컨테이너에 work dir를 마운트할 때 호스트 경로를 써야 합니다.
# 컨테이너 내부 경로(/data/analysis/...)를 쓰면 호스트 Docker 데몬이 경로를 찾지 못해
# DeepVariant 등 태스크 컨테이너가 작업 디렉터리를 마운트하지 못합니다.
HOST_WORK_DIR="${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}/work"

# Docker binary path on the host — bind-mounted into the container so
# Nextflow (docker.enabled=true) can spawn task containers via the host daemon.
DOCKER_BIN="$(which docker)"

# 비-root 로 Nextflow를 띄우면 work/ 해시 디렉터리도 동일 uid 로 생성되어
# NXF_DOCKER_TASK_USER 태스크 컨테이너와 쓰기 권한이 맞음. (root 메인 + 1000 태스크면 .command.trace Permission denied)
# passwd 에 없으면 HOME 이 비어 Nextflow가 잘못된 경로에 쓰므로 HOME/NXF_HOME 고정.
# docker.sock 이 root:docker 이면 컨테이너에 docker GID 보조 그룹 필요.
DOCKER_GROUP_ARGS=()
if DOCKER_SOCK_GID="$(getent group docker 2>/dev/null | cut -d: -f3)" && [ -n "${DOCKER_SOCK_GID}" ]; then
    DOCKER_GROUP_ARGS=(--group-add "${DOCKER_SOCK_GID}")
fi

# docker ps NAME: order id (-w) + sample (sanitize for Docker: [a-zA-Z0-9][a-zA-Z0-9_.-]*)
NF_DOCKER_NAME_RAW="gx-exome-${WORK_DIR}-${SAMPLE_NAME}"
NF_DOCKER_NAME="$(printf '%s' "$NF_DOCKER_NAME_RAW" | sed -e 's/[^a-zA-Z0-9_.-]/-/g' -e 's/^[-_.]*//' -e 's/^$/carrier-unknown/')"
NF_DOCKER_NAME="${NF_DOCKER_NAME:0:200}"

docker run --rm -t --name "$NF_DOCKER_NAME" \
    -u "${CHOWN_SPEC}" \
    "${DOCKER_GROUP_ARGS[@]}" \
    -e HOME=/tmp \
    -e NXF_HOME=/tmp/.nextflow \
    -v "${DATA_DIR}/fastq:/data/fastq:ro" \
    -v "${DATA_DIR}/analysis:/data/analysis" \
    -v "${DATA_DIR}/output:/data/output" \
    -v "${DATA_DIR}/log:/data/log" \
    -v "${DATA_DIR}/data:/app/data:ro" \
    -v "${DATA_DIR}/bin:/app/bin:ro" \
    -v "${DATA_DIR}/fastq:${DATA_DIR}/fastq:ro" \
    -v "${DATA_DIR}/analysis:${DATA_DIR}/analysis" \
    -v "${DATA_DIR}/output:${DATA_DIR}/output" \
    -v "${DATA_DIR}/log:${DATA_DIR}/log" \
    -v "${DATA_DIR}/data:${DATA_DIR}/data:ro" \
    -v /data/reference:/data/reference:ro \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v "${DOCKER_BIN}:/usr/local/bin/docker:ro" \
    -e NXF_OPTS="-Xms1g -Xmx4g" \
    -e NXF_CACHE_DIR="${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}/.nextflow" \
    -e NXF_DATA_DIR="${DATA_DIR}" \
    -e NXF_DOCKER_TASK_USER="${CHOWN_SPEC}" \
    -e HOST_WORK_DIR="${HOST_WORK_DIR}" \
    gx-exome:latest \
    bash -c "
        cd ${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME} && \
        nextflow -log ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/nextflow.log run /app/bin/main.nf \
            -ansi-log false ${NEXTFLOW_RESUME} ${NXF_PROFILE} \
            --fastq_dir ${DATA_DIR}/fastq/${WORK_DIR}/${SAMPLE_NAME} \
            --ref_fasta ${DATA_DIR}/data/refs/GRCh38.fasta \
            --ref_fai ${DATA_DIR}/data/refs/GRCh38.fasta.fai \
            --ref_dict ${DATA_DIR}/data/refs/GRCh38.dict \
            --ref_bwa_indices ${DATA_DIR}/data/refs/bwa_index \
            --backbone_bed ${DATA_DIR}/data/bed/Twist_Exome2.0_plus_Comprehensive_Exome_Spikein_targets_covered_annotated_hg38.bed \
            --backbone_bed_gz ${DATA_DIR}/data/bed/Twist_Exome2.0_plus_Comprehensive_Exome_Spikein_targets_covered_annotated_hg38.bed.gz \
            --backbone_bed_tbi ${DATA_DIR}/data/bed/Twist_Exome2.0_plus_Comprehensive_Exome_Spikein_targets_covered_annotated_hg38.bed.gz.tbi \
            --vep_cache_dir ${DATA_DIR}/data/refs/vep_cache \
            --dark_genes_plus_bed ${DATA_DIR}/data/bed/dark_genes_plus.bed \
            --hba_bed ${DATA_DIR}/data/bed/hba_targets.bed \
            --cyp21a2_bed ${DATA_DIR}/data/bed/cyp21a2_targets.bed \
            --eh_catalog ${DATA_DIR}/data/bed/variant_catalog_grch38.json \
            --aligner ${ALIGNER} \
            --variant_caller ${VARIANT_CALLER} \
            --skip_vep ${SKIP_VEP} \
            ${SKIP_CNV} \
            --outdir ${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME} \
            --output_dir ${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME} \
            --sample_name ${SAMPLE_NAME} \
            ${CLEANUP} \
            -work-dir \"\$HOST_WORK_DIR\" \
            -with-report ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/report.html \
            -with-trace ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/trace.txt \
            -with-timeline ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/timeline.html \
            -with-dag ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/dag.html
    "

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}✅ Analysis completed successfully!${NC}"
    echo ""
    repair_order_tree_permissions
    echo "Results:"
    echo "  Output: ${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}"
    echo "  Logs: ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}"
    echo ""

    # analysis.completed 마커 (FASTQ 가 ro/타 소유면 실패할 수 있음 — 파이프라인 성공과 무관)
    if touch "${FASTQ_DIR}/analysis.completed" 2>/dev/null; then
        echo -e "${GREEN}✅ Created analysis.completed marker${NC}"
    elif touch "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}/analysis.completed" 2>/dev/null; then
        echo -e "${GREEN}✅ Created analysis.completed in output dir${NC}"
    else
        echo -e "${YELLOW}⚠ Skipped analysis.completed (no write permission); pipeline succeeded.${NC}"
    fi
    echo "  Finished at: $(pipeline_timestamp)"

    exit 0
else
    echo -e "${RED}❌ Analysis failed with exit code: ${EXIT_CODE}${NC}"
    echo "  Stopped at: $(pipeline_timestamp)"
    echo ""
    repair_order_tree_permissions || true
    echo "Check logs:"
    echo "  ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/nextflow.log"
    echo ""
    exit $EXIT_CODE
fi
