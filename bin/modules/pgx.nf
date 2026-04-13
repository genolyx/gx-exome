// ============================================================
// pgx.nf — PharmCAT + Python finalizer (two processes)
//
// PharmCAT images do not always expose python3 on PATH; JSON contract is
// built in a dedicated Python container (label pgx_finalize).
// ============================================================

process RUN_PGX_PHARMCAT {
    tag "$sample_id"
    label 'pgx'

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(outside_call)

    output:
    tuple val(sample_id), path("pgx_staging"), emit: staged

    script:
    """
    set +e
    export TMPDIR=\$PWD
    export JAVA_MAX_HEAP=\${JAVA_MAX_HEAP:-6G}

    mkdir -p pgx_staging

    pharmcat_pipeline -V 2>/dev/null | head -n 1 > pgx_staging/pharmcat_version.txt || echo "unknown" > pgx_staging/pharmcat_version.txt

    pharmcat_pipeline -0 -o pgx_staging -bf ${sample_id}_pgx -reporterJson "${vcf}" 2> pgx_staging/pharmcat.stderr.log
    echo \$? > pgx_staging/pharmcat_exit.code

    # If Aldy produced a CYP2D6 outside call, re-run Phenotyper + Reporter
    # with the outside call merged in. pharmcat_pipeline does not support -po;
    # the Java tool (pharmcat.jar) does.
    if [ -s "${outside_call}" ]; then
        cp "${outside_call}" pgx_staging/aldy_cyp2d6_outside_call.tsv
        MATCH_JSON="pgx_staging/${sample_id}_pgx.match.json"
        PHARMCAT_JAR=\$(ls /pharmcat/pharmcat*.jar 2>/dev/null | head -1)
        if [ -f "\$MATCH_JSON" ] && [ -n "\$PHARMCAT_JAR" ]; then
            java -Xmx\${JAVA_MAX_HEAP} -jar "\$PHARMCAT_JAR" \\
                -phenotyper -pi "\$MATCH_JSON" \\
                -po "${outside_call}" \\
                -reporter -reporterJson \\
                -o pgx_staging -bf ${sample_id}_pgx \\
                2>> pgx_staging/pharmcat.stderr.log
            echo \$? > pgx_staging/pharmcat_po_exit.code
        fi
    fi
    exit 0
    """
}

process FINALIZE_PGX_JSON {
    tag "$sample_id"
    label 'pgx_finalize'
    publishDir "${params.outdir}/pgx", mode: 'copy'

    input:
    tuple val(sample_id), path(pgx_staging), path(pgx_finalize_py), val(reference_genome), val(pgx_image)

    output:
    path "pgx_meta.json",   emit: meta
    path "pgx_result.json", emit: result
    path "${sample_id}_pgx.match.json", optional: true
    path "${sample_id}_pgx.phenotype.json", optional: true
    path "${sample_id}_pgx.report.json", optional: true
    path "${sample_id}_pgx.report.html", optional: true
    path "pharmcat.stderr.log", optional: true
    path "pharmcat_version.txt", optional: true
    path "aldy_cyp2d6_outside_call.tsv", optional: true
    path "pgx_summary.txt", optional: true
    path "pgx_finalize.stderr", optional: true

    script:
    """
    set +e
    # v2: stageInMode=copy fix — ensures pgx_finalize.py is a real file, not a dangling symlink
    cp -a pgx_staging/. .
    TOOL_VER=\$(head -n1 pharmcat_version.txt 2>/dev/null || echo unknown)
    EC=\$(cat pharmcat_exit.code 2>/dev/null || echo 1)

    python3 ${pgx_finalize_py} \\
        --sample-id "${sample_id}" \\
        --reference-genome "${reference_genome}" \\
        --image "${pgx_image}" \\
        --tool-version "\$TOOL_VER" \\
        --exit-code "\$EC" \\
        --stderr pharmcat.stderr.log \\
        --work-dir . 2> pgx_finalize.stderr

    PY=\$?
    if [ "\$PY" -ne 0 ]; then
        TS=\$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || date -Iseconds 2>/dev/null || date)
        cat > pgx_meta.json << EOFM
{
  "schema": "gx_exome_pgx_meta_v1",
  "tool": "PharmCAT",
  "tool_version": "\$TOOL_VER",
  "container_image": "${pgx_image}",
  "reference_genome": "${reference_genome}",
  "resource_versions": {
    "pharmcat_docker_image": "${pgx_image}",
    "note": "CPIC/PharmGKB data bundled with PharmCAT release"
  },
  "sample_id": "${sample_id}",
  "timestamp": "\$TS",
  "exit_status": "error",
  "message": "FINALIZE_PGX_JSON: pgx_finalize.py failed — see pgx_finalize.stderr in pgx/"
}
EOFM
        export SAMPLE_ID="${sample_id}"
        python3 << 'PYFALLBACK'
import json
import os

sid = os.environ.get("SAMPLE_ID", "")
detail = ""
try:
    with open("pgx_finalize.stderr", "r", encoding="utf-8", errors="replace") as f:
        detail = f.read()[:4000]
except OSError:
    pass
out = {
    "schema": "gx_exome_pgx_result_v1",
    "sample_id": sid,
    "error": "pgx_finalize_failed",
}
if detail.strip():
    out["detail"] = detail
with open("pgx_result.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
    f.write("\\n")
PYFALLBACK
    fi

    # Build a concise actionable summary from the phenotype JSON
    export SAMPLE_ID="${sample_id}"
    python3 << 'PYSUMMARY' 2>/dev/null || true
import json
import os

sid = os.environ.get("SAMPLE_ID", "")
phen_file = f"{sid}_pgx.phenotype.json"

if not os.path.isfile(phen_file):
    with open("pgx_summary.txt", "w") as f:
        f.write(f"PGx Summary - {sid}\\nNo phenotype data available.\\n")
    raise SystemExit(0)

with open(phen_file) as f:
    data = json.load(f)

lines = []
lines.append(f"PGx Summary - {sid}")
lines.append("=" * 72)
lines.append("")

actionable = []
normal = []

RISK_FUNCTIONS = {"no function", "decreased function", "unfavorable response allele"}
SKIP_PHENOTYPES = {"no result", "n/a", ""}

seen_genes = set()
gene_reports = data.get("geneReports", {})
# v2: geneReports.CPIC.<gene> / geneReports.DPWG.<gene>
# v3: geneReports.<gene> (flat)
source_gene_pairs = []
if "CPIC" in gene_reports or "DPWG" in gene_reports:
    for source in ("CPIC", "DPWG"):
        reports = gene_reports.get(source, {})
        if isinstance(reports, dict):
            for gn in sorted(reports.keys()):
                source_gene_pairs.append((source, gn, reports[gn]))
else:
    for gn in sorted(gene_reports.keys()):
        gd = gene_reports[gn]
        if isinstance(gd, dict) and gd.get("recommendationDiplotypes"):
            source_gene_pairs.append(("CPIC", gn, gd))

for source, gene_name, gene_data in source_gene_pairs:
    if gene_name in seen_genes:
        continue
    if not isinstance(gene_data, dict):
        continue
    call_src = gene_data.get("callSource", "")
    rec_dips = gene_data.get("recommendationDiplotypes", [])
    if not rec_dips:
        continue

        for dip in rec_dips:
            a1 = dip.get("allele1") or {}
            a2 = dip.get("allele2") or {}
            n1 = a1.get("name", "")
            n2 = a2.get("name", "")
            fn1 = (a1.get("function") or "").strip()
            fn2 = (a2.get("function") or "").strip()
            phenotypes = dip.get("phenotypes", [])
            activity = dip.get("activityScore")
            diplotype = f"{n1}/{n2}" if n2 else n1
            phenotype_str = ", ".join(p for p in phenotypes if p) if phenotypes else ""

            if phenotype_str.lower() in SKIP_PHENOTYPES:
                continue

            functions = [f for f in (fn1, fn2) if f]
            has_risk = any(f.lower() in RISK_FUNCTIONS for f in functions)

            row = f"{gene_name:12s}  {diplotype:35s}  {phenotype_str}"
            if activity and str(activity) != "None":
                row += f"  (activity: {activity})"
            if call_src == "OUTSIDE":
                row += "  [Aldy]"
            if has_risk:
                row += "  <<"

            if has_risk:
                actionable.append((gene_name, row, functions, phenotype_str))
            else:
                normal.append(row)

            seen_genes.add(gene_name)
            break

if actionable:
    lines.append("*** ACTIONABLE — genes with reduced/no-function alleles ***")
    lines.append("-" * 72)
    for gene, row, fns, phen in actionable:
        lines.append(f"  {row}")
        fn_detail = " + ".join(fns)
        lines.append(f"{'':14s}  Functions: {fn_detail}")
        lines.append("")
    lines.append("  >> Check report.html for specific drug dose recommendations.")
    lines.append("")

if normal:
    lines.append("Normal / reference results:")
    lines.append("-" * 72)
    for row in normal:
        lines.append(f"  {row}")
    lines.append("")

lines.append("Full report: pgx/<sample>_pgx.report.html")

with open("pgx_summary.txt", "w") as f:
    f.write("\\n".join(lines) + "\\n")

PYSUMMARY

    exit 0
    """
}
