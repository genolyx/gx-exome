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
    -h, --help                 Show this help message

Example:
    $0 -w 2601 -s Sample_A10
    $0 --work-dir 2601 --sample Sample_A10 --cleanup

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

# 출력 디렉토리 생성 (Nextflow .nextflow 캐시용)
mkdir -p "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}"
mkdir -p "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}/.nextflow"
mkdir -p "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}"
mkdir -p "${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}"

echo "======================================"
echo "Dark Gene Pipeline - Sample Analysis"
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
echo ""

# Docker 이미지 확인
if ! docker images | grep -q "carrier-screening"; then
    echo -e "${RED}Error: Docker image 'carrier-screening' not found${NC}"
    echo "Please build the image first:"
    echo "  docker-compose -f docker/docker-compose.yml build"
    exit 1
fi

# 분석 시작
echo -e "${YELLOW}Starting analysis...${NC}"
echo ""

# Docker 컨테이너 실행 (root로 실행 후 출력 파일 소유권 수정)
# nextflow.config가 projectDir/../data/refs, data/bed 경로 사용 → /app/data 마운트 필요
docker run --rm -t \
    -v "${DATA_DIR}/fastq:/data/fastq:ro" \
    -v "${DATA_DIR}/analysis:/data/analysis" \
    -v "${DATA_DIR}/output:/data/output" \
    -v "${DATA_DIR}/log:/data/log" \
    -v "${DATA_DIR}/data:/app/data:ro" \
    -e NXF_OPTS="-Xms1g -Xmx4g" \
    -e NXF_CACHE_DIR="/data/analysis/${WORK_DIR}/${SAMPLE_NAME}/.nextflow" \
    carrier-screening:latest \
    bash -c "
        cd /data/analysis/${WORK_DIR}/${SAMPLE_NAME} && \
        nextflow -log /data/log/${WORK_DIR}/${SAMPLE_NAME}/nextflow.log run /app/bin/main.nf \
            -ansi-log false \
            -resume \
            --fastq_dir /data/fastq/${WORK_DIR}/${SAMPLE_NAME} \
            ${SKIP_CNV} \
            --outdir /data/analysis/${WORK_DIR}/${SAMPLE_NAME} \
            --output_dir /data/output/${WORK_DIR}/${SAMPLE_NAME} \
            --sample_name ${SAMPLE_NAME} \
            ${CLEANUP} \
            -work-dir ./work \
            -with-report /data/log/${WORK_DIR}/${SAMPLE_NAME}/report.html \
            -with-trace /data/log/${WORK_DIR}/${SAMPLE_NAME}/trace.txt \
            -with-timeline /data/log/${WORK_DIR}/${SAMPLE_NAME}/timeline.html \
            -with-dag /data/log/${WORK_DIR}/${SAMPLE_NAME}/dag.html
    "

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}✅ Analysis completed successfully!${NC}"
    echo ""
    # Docker가 root로 생성한 파일 소유권 수정 (passwordless sudo 가능 시)
    if sudo -n true 2>/dev/null; then
        sudo chown -R "$(id -u):$(id -g)" "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE_NAME}" "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}" "${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}" 2>/dev/null || true
    fi
    echo "Results:"
    echo "  Output: ${DATA_DIR}/output/${WORK_DIR}/${SAMPLE_NAME}"
    echo "  Logs: ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}"
    echo ""
    
    # analysis.completed 마커 생성
    touch "${FASTQ_DIR}/analysis.completed"
    echo -e "${GREEN}✅ Created analysis.completed marker${NC}"
    
    exit 0
else
    echo -e "${RED}❌ Analysis failed with exit code: ${EXIT_CODE}${NC}"
    echo ""
    echo "Check logs:"
    echo "  ${DATA_DIR}/log/${WORK_DIR}/${SAMPLE_NAME}/nextflow.log"
    echo ""
    exit $EXIT_CODE
fi
