#!/usr/bin/env bash
#
# Lightweight backup of comparable pipeline outputs (VCF, QC stats, optional raw VCF, summary).
# Usage: run from project root or any cwd; pass absolute or relative paths as needed.
#
set -euo pipefail

usage() {
    echo "Usage: $0 -w WORK_DIR -s SAMPLE [-d DATA_DIR] [--tag TAG] [--no-raw] [--no-summary]"
    echo "  -d DATA_DIR   Project/data root (default: current directory). Same as run_analysis.sh -d ."
    echo "  --tag TAG     Backup folder name under results_backup/ (default: date + bwa_gatk_minimal)"
    echo "  --no-raw      Skip *_raw.vcf.gz from analysis/variant (saves space)"
    echo "  --no-summary  Skip output summary *.txt"
    exit 1
}

WORK_DIR=""
SAMPLE=""
DATA_DIR="$(pwd)"
TAG=""
INCLUDE_RAW=1
INCLUDE_SUMMARY=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--work-dir) WORK_DIR="$2"; shift 2 ;;
        -s|--sample)   SAMPLE="$2"; shift 2 ;;
        -d|--data-dir) DATA_DIR="$2"; shift 2 ;;
        --tag)         TAG="$2"; shift 2 ;;
        --no-raw)      INCLUDE_RAW=0; shift ;;
        --no-summary)  INCLUDE_SUMMARY=0; shift ;;
        -h|--help)     usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$WORK_DIR" || -z "$SAMPLE" ]]; then
    usage
fi

DATA_DIR="$(realpath "$DATA_DIR")"

if [[ -z "$TAG" ]]; then
    TAG="$(date +%Y%m%d)_bwa_gatk_minimal"
fi

DST="${DATA_DIR}/results_backup/${TAG}"
mkdir -p "${DST}/vcf" "${DST}/qc" "${DST}/alignment" "${DST}/variant_raw" "${DST}/summary"

shopt -s nullglob

# Final filtered VCF (output)
for f in "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/vcf/"*.vcf.gz \
         "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/vcf/"*.vcf.gz.tbi; do
    [[ -e "$f" ]] && cp -a "$f" "${DST}/vcf/"
done

# QC in output
for f in "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/qc/"*.stats.txt \
         "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/qc/"*_duplicate_metrics.txt \
         "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/qc/"*_qc_metrics.txt \
         "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/qc/"*_target_coverage.txt; do
    [[ -e "$f" ]] && cp -a "$f" "${DST}/qc/"
done

# Duplicate metrics may live only under analysis/alignment
for f in "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE}/alignment/"*_duplicate_metrics.txt; do
    [[ -e "$f" ]] && cp -a "$f" "${DST}/alignment/"
done

if [[ "$INCLUDE_RAW" -eq 1 ]]; then
    for f in "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE}/variant/"*_raw.vcf.gz \
             "${DATA_DIR}/analysis/${WORK_DIR}/${SAMPLE}/variant/"*_raw.vcf.gz.tbi; do
        [[ -e "$f" ]] && cp -a "$f" "${DST}/variant_raw/"
    done
fi

if [[ "$INCLUDE_SUMMARY" -eq 1 ]]; then
    for f in "${DATA_DIR}/output/${WORK_DIR}/${SAMPLE}/summary/"*.txt; do
        [[ -e "$f" ]] && cp -a "$f" "${DST}/summary/"
    done
fi

shopt -u nullglob

{
    echo "created_at=$(date -Iseconds)"
    echo "data_dir=${DATA_DIR}"
    echo "work_dir=${WORK_DIR}"
    echo "sample=${SAMPLE}"
    echo "destination=${DST}"
    echo "include_raw=${INCLUDE_RAW}"
    echo "include_summary=${INCLUDE_SUMMARY}"
} >> "${DST}/MANIFEST.txt"

echo "Backup done: ${DST}"
