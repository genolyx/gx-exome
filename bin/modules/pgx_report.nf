// ============================================================
// pgx_report.nf — Combined PGx panel HTML report
//
// Merges PharmCAT (21 genes) + custom extended panel (30 genes)
// into a single interactive HTML table for portal display.
// ============================================================

process GENERATE_PGX_PANEL_REPORT {
    tag "$sample_id"
    label 'pgx_finalize'
    publishDir "${params.outdir}/pgx", mode: 'copy', pattern: 'pgx_panel_report.html'

    input:
    tuple val(sample_id), path(pgx_staging), path(custom_json), path(report_script)

    output:
    path "pgx_panel_report.html", emit: report_html

    script:
    """
    python3 ${report_script} \
        --sample-id "${sample_id}" \
        --phenotype-json "${pgx_staging}/${sample_id}_pgx.phenotype.json" \
        --report-json "${pgx_staging}/${sample_id}_pgx.report.json" \
        --custom-json "${custom_json}" \
        --output pgx_panel_report.html

    if [ ! -f pgx_panel_report.html ]; then
        echo '<html><body><h1>PGx Panel Report</h1><p>Error: report generation failed.</p></body></html>' > pgx_panel_report.html
    fi
    """
}
