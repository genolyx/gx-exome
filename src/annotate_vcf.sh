#!/bin/bash
#
# Annotate one raw VCF with Ensembl VEP. Does not run the Nextflow pipeline.
# Static VEP flags come from bin/modules/vep_annotation.flags, the same file
# VEP_ANNOTATION reads. snpEff is not invoked here.
#

set -euo pipefail

RED='\033[0;31m'
NC='\033[0m'

VEP_IMAGE="ensemblorg/ensembl-vep:release_111.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FLAG_FILE="${REPO_DIR}/bin/modules/vep_annotation.flags"

usage() {
    cat << EOF
Usage: $0 --vcf VCF --assembly GRCh38|GRCh37 --out OUT.vcf.gz [OPTIONS]

Annotate a raw VCF (no CSQ / ANN required) with the same VEP command used by
VEP_ANNOTATION. Writes OUT.vcf.gz and OUT.vcf.gz.tbi.

Options:
    --vcf VCF            Input .vcf or .vcf.gz (required)
    --assembly NAME      GRCh38 or GRCh37 (required)
    --out PATH           Output .vcf.gz path (required). A .tbi is written beside it.
    --sample ID          Sample id for the stats file name (default: input filename stem)
    --data-dir DIR       Same default as run_analysis.sh: ./data
                         Cache:  \${DATA_DIR}/data/refs/vep_cache
                         FASTA:  \${DATA_DIR}/data/refs/\${assembly}.fasta (+ .fai)
    --fork N             VEP --fork (default: number of CPUs)
    -h, --help           Show this help message

EOF
    exit 1
}

die() {
    echo -e "${RED}Error: $*${NC}" >&2
    exit 1
}

VCF=""
ASSEMBLY=""
OUT=""
SAMPLE=""
DATA_DIR="$(pwd)/data"
FORK=""

require_arg() {
    [[ $# -ge 2 && -n "${2:-}" ]] || die "$1 requires a value"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)      require_arg "$1" "${2:-}"; VCF="$2"; shift 2 ;;
        --assembly) require_arg "$1" "${2:-}"; ASSEMBLY="$2"; shift 2 ;;
        --out)      require_arg "$1" "${2:-}"; OUT="$2"; shift 2 ;;
        --sample)   require_arg "$1" "${2:-}"; SAMPLE="$2"; shift 2 ;;
        --data-dir) require_arg "$1" "${2:-}"; DATA_DIR="$2"; shift 2 ;;
        --fork)     require_arg "$1" "${2:-}"; FORK="$2"; shift 2 ;;
        -h|--help)  usage ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "$VCF" ]] || die "--vcf is required"
[[ -n "$ASSEMBLY" ]] || die "--assembly is required"
[[ -n "$OUT" ]] || die "--out is required"
[[ "$ASSEMBLY" == "GRCh38" || "$ASSEMBLY" == "GRCh37" ]] || die "--assembly must be GRCh38 or GRCh37"
[[ "$OUT" == *.vcf.gz ]] || die "--out must be a .vcf.gz path"
[[ -f "$FLAG_FILE" ]] || die "Shared VEP flags not found: $FLAG_FILE"

if [[ -z "$FORK" ]]; then
    FORK="$(nproc)"
fi
[[ "$FORK" =~ ^[1-9][0-9]*$ ]] || die "--fork must be a positive integer"

if [[ ! -f "$VCF" ]]; then
    die "Input VCF not found: $VCF"
fi
if [[ ! -s "$VCF" ]]; then
    die "Input VCF is empty: $VCF"
fi
case "$VCF" in
    *.vcf.gz)
        gzip -t "$VCF" || die "Input VCF is not a valid gzip file: $VCF"
        gz_empty=1
        set +o pipefail
        gzip -dc "$VCF" | head -c 1 | grep -q . && gz_empty=0
        set -o pipefail
        if [[ "$gz_empty" -eq 1 ]]; then
            die "Input VCF is empty: $VCF"
        fi
        ;;
    *.vcf) ;;
    *) die "Input VCF must be .vcf or .vcf.gz: $VCF" ;;
esac

DATA_DIR="$(realpath "$DATA_DIR")"
cache_path="${DATA_DIR}/data/refs/vep_cache"
fasta_path="${DATA_DIR}/data/refs/${ASSEMBLY}.fasta"
fai_path="${fasta_path}.fai"

if [[ ! -d "$cache_path" ]]; then
    die "VEP cache directory not found: ${cache_path}"
fi
CACHE="$(realpath "$cache_path")"
shopt -s nullglob
assembly_caches=( "${CACHE}/homo_sapiens/"*_"${ASSEMBLY}" )
shopt -u nullglob
if [[ ${#assembly_caches[@]} -eq 0 ]]; then
    die "No ${ASSEMBLY} VEP cache under ${CACHE}/homo_sapiens (refusing to annotate with a different assembly cache)"
fi

if [[ ! -s "$fasta_path" ]]; then
    die "Reference FASTA not found: ${fasta_path}"
fi
if [[ ! -s "$fai_path" ]]; then
    die "Reference FASTA index not found: ${fai_path}"
fi
FASTA="$(realpath "$fasta_path")"
FAI="$(realpath "$fai_path")"

if [[ -z "$SAMPLE" ]]; then
    stem="$(basename "$VCF")"
    stem="${stem%.gz}"
    stem="${stem%.vcf}"
    SAMPLE="$stem"
fi
[[ -n "$SAMPLE" && "$SAMPLE" != */* ]] || die "--sample must be a single path segment"

if ! docker image inspect "$VEP_IMAGE" >/dev/null 2>&1; then
    die "Docker image ${VEP_IMAGE} is not present locally (this script does not build or pull it)"
fi

VCF="$(realpath "$VCF")"
OUT="$(realpath -m "$OUT")"
OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"

STATS="${OUT_DIR}/${SAMPLE}_vep_summary.html"
PARTIAL="${OUT_DIR}/.${SAMPLE}.vep.partial.vcf.gz"

VEP_FLAGS=()
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" || "$line" == \#* ]] && continue
    VEP_FLAGS+=("$line")
done < "$FLAG_FILE"
[[ ${#VEP_FLAGS[@]} -gt 0 ]] || die "Shared VEP flags file is empty: $FLAG_FILE"

success=0
cleanup() {
    if [[ "$success" -ne 1 ]]; then
        rm -f "$PARTIAL" "${PARTIAL}.tbi" "$STATS"
    fi
}
trap cleanup EXIT

FLAG_PAYLOAD="$(printf '%s\n' "${VEP_FLAGS[@]}")"

set +e
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e TMPDIR=/tmp \
    -e VCF="$VCF" \
    -e PARTIAL="$PARTIAL" \
    -e STATS="$STATS" \
    -e CACHE="$CACHE" \
    -e FASTA="$FASTA" \
    -e ASSEMBLY="$ASSEMBLY" \
    -e FORK="$FORK" \
    -e FLAG_PAYLOAD="$FLAG_PAYLOAD" \
    -v "${VCF}:${VCF}:ro" \
    -v "${FASTA}:${FASTA}:ro" \
    -v "${FAI}:${FAI}:ro" \
    -v "${CACHE}:${CACHE}:ro" \
    -v "${OUT_DIR}:${OUT_DIR}" \
    "$VEP_IMAGE" \
    bash -lc '
set -euo pipefail
cleanup_partial() {
    ec=$?
    if [[ $ec -ne 0 ]]; then
        rm -f "$PARTIAL" "${PARTIAL}.tbi" "$STATS"
    fi
}
trap cleanup_partial EXIT
mapfile -t FLAGS <<< "$FLAG_PAYLOAD"
vep \
    --input_file "$VCF" \
    --output_file "$PARTIAL" \
    --stats_file "$STATS" \
    "${FLAGS[@]}" \
    --dir_cache "$CACHE" \
    --assembly "$ASSEMBLY" \
    --fasta "$FASTA" \
    --fork "$FORK"
tabix -p vcf "$PARTIAL"
csq_ok=0
set +o pipefail
gzip -dc "$PARTIAL" | grep -m1 -q "^##INFO=<ID=CSQ," && csq_ok=1
set -o pipefail
if [[ "$csq_ok" -ne 1 ]]; then
    echo "VEP output is missing INFO/CSQ header: $PARTIAL" >&2
    exit 1
fi
'
ec=$?
set -e

if [[ "$ec" -ne 0 ]]; then
    exit "$ec"
fi

mv -f "$PARTIAL" "$OUT"
mv -f "${PARTIAL}.tbi" "${OUT}.tbi"
success=1

echo "Wrote ${OUT}"
echo "Wrote ${OUT}.tbi"
