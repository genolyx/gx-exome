// ============================================================
// pgx.nf — Optional PharmCAT pharmacogenomics (PGx) add-on
//
// Runs after the final sample VCF (annotated if VEP ran) is available.
// Does not change SNV/indel calling. Soft-fail: always exits 0; status in pgx_meta.json.
// Outputs under ${params.outdir}/pgx/ (and copied to output_dir/pgx/ in onComplete).
// ============================================================

process RUN_PGX_PHARMCAT {
    tag "$sample_id"
    label 'pgx'
    publishDir "${params.outdir}/pgx", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(pgx_finalize_py)
    val reference_genome

    output:
    path "pgx_meta.json",   emit: meta
    path "pgx_result.json", emit: result
    path "${sample_id}_pgx.match.json", optional: true
    path "${sample_id}_pgx.phenotype.json", optional: true
    path "${sample_id}_pgx.report.json", optional: true
    path "${sample_id}_pgx.report.html", optional: true
    path "pharmcat.stderr.log", optional: true
    path "pharmcat_version.txt", optional: true

    script:
    """
    set +e
    export TMPDIR=\$PWD
    export JAVA_MAX_HEAP=\${JAVA_MAX_HEAP:-6G}

    # Version string for metadata (best-effort before pipeline run)
    pharmcat_pipeline -V 2>/dev/null | head -n 1 > pharmcat_version.txt || echo "unknown" > pharmcat_version.txt
    TOOL_VER=\$(cat pharmcat_version.txt)

    pharmcat_pipeline -o . -bf ${sample_id}_pgx -reporterJson -reporterHtml "${vcf}" 2> pharmcat.stderr.log
    RC=\$?

    PYBIN=python3
    command -v python3 >/dev/null 2>&1 || PYBIN=python
    \$PYBIN ${pgx_finalize_py} \\
        --sample-id "${sample_id}" \\
        --reference-genome "${reference_genome}" \\
        --image "${params.pgx_container}" \\
        --tool-version "\$TOOL_VER" \\
        --exit-code "\$RC" \\
        --stderr pharmcat.stderr.log \\
        --work-dir .

    PY=\$?
    if [ "\$PY" -ne 0 ]; then
        TS=\$(date -Iseconds 2>/dev/null || date)
        cat > pgx_meta.json << EOFM
{
  "schema": "gx_exome_pgx_meta_v1",
  "tool": "PharmCAT",
  "tool_version": "\$TOOL_VER",
  "container_image": "${params.pgx_container}",
  "reference_genome": "${reference_genome}",
  "resource_versions": {
    "pharmcat_docker_image": "${params.pgx_container}",
    "note": "CPIC/PharmGKB data bundled with PharmCAT release"
  },
  "sample_id": "${sample_id}",
  "timestamp": "\$TS",
  "exit_status": "error",
  "message": "pgx_finalize.py failed (python3 missing or script error); see stderr"
}
EOFM
        echo '{"schema":"gx_exome_pgx_result_v1","sample_id":"${sample_id}","error":"pgx_finalize_failed"}' > pgx_result.json
    fi

    exit 0
    """
}
