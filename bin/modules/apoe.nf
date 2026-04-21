// ============================================================
// apoe.nf — APOE ε2/ε3/ε4 from rs429358 + rs7412 (GRCh38 VCF)
//
// Complements pgx_custom (same SNPs) with diplotype + AD risk context.
// Uses python:3.12-bookworm + pysam like RUN_PGX_CUSTOM.
// ============================================================

process RUN_APOE {
    tag "$sample_id"
    label 'apoe'
    publishDir "${params.outdir}/apoe", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(bam), path(bai), path(apoe_py)

    output:
    path "apoe_result.json", emit: result_json
    path "apoe_summary.txt", emit: summary_txt
    path "apoe_review.txt",   emit: review_txt

    script:
    """
    set +e
    export HOME=\$PWD
    pip install --user pysam 2>&1 | tail -3
    export PATH="\$PWD/.local/bin:\$PATH"

    python3 ${apoe_py} \\
        --vcf "${vcf}" \\
        --bam "${bam}" \\
        --sample-id "${sample_id}" \\
        --json-out apoe_result.json \\
        --summary-out apoe_summary.txt \\
        --review-out apoe_review.txt

    exit 0
    """
}
