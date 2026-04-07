// ============================================================
// aldy.nf — Aldy CYP2D6 star-allele caller
//
// Aldy calls CYP2D6 diplotypes from the BAM (not VCF) using its
// integer-linear-programming model.  The exome profile assumes
// exactly 2 gene copies (no CNV/fusion detection possible from
// capture data).
//
// Output: PharmCAT "outside call" TSV so PharmCAT can translate
// the diplotype into phenotype + CPIC drug recommendations.
// ============================================================

process RUN_ALDY_CYP2D6 {
    tag "$sample_id"
    label 'aldy'
    publishDir "${params.outdir}/pgx", mode: 'copy', pattern: 'aldy_cyp2d6_result.json'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("aldy_cyp2d6_outside_call.tsv"), emit: outside_call
    path "${sample_id}.CYP2D6.aldy", optional: true, emit: raw
    path "aldy_cyp2d6.stderr", optional: true
    path "aldy_cyp2d6_result.json", emit: result_json

    script:
    """
    set +e

    export HOME=\$PWD
    pip install --user aldy 2>&1 | tail -5
    export PATH="\$PWD/.local/bin:\$PATH"

    aldy genotype \\
        -p exome \\
        -g CYP2D6 \\
        -o "${sample_id}.CYP2D6.aldy" \\
        "${bam}" 2> aldy_cyp2d6.stderr

    export ALDY_EC=\$?
    export ALDY_SAMPLE_ID="${sample_id}"

    python3 << 'PYPARSE'
import json
import os
import re

sample_id = os.environ["ALDY_SAMPLE_ID"]
aldy_ec = int(os.environ.get("ALDY_EC", "1"))
aldy_file = f"{sample_id}.CYP2D6.aldy"

stderr_text = ""
try:
    with open("aldy_cyp2d6.stderr", "r", encoding="utf-8", errors="replace") as f:
        stderr_text = f.read()
except OSError:
    pass

diplotype = None
minor_detail = None

if os.path.isfile(aldy_file):
    with open(aldy_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            cols = line.split("\\t")
            if len(cols) >= 4 and cols[1] == "CYP2D6":
                raw_call = cols[3]
                parts = raw_call.split("/")
                cleaned = []
                for p in parts:
                    p = p.strip()
                    m = re.match(r"\\*(\\d+)", p)
                    if m:
                        cleaned.append(f"*{m.group(1)}")
                    elif p.startswith("*"):
                        cleaned.append(p)
                if len(cleaned) >= 2:
                    diplotype = "/".join(cleaned[:2])
                    minor_detail = raw_call
                break

result = {
    "schema": "gx_exome_aldy_cyp2d6_v1",
    "sample_id": sample_id,
    "gene": "CYP2D6",
    "tool": "Aldy",
    "profile": "exome",
    "exit_code": aldy_ec,
}

if diplotype:
    result["diplotype"] = diplotype
    result["status"] = "success"
    if minor_detail and minor_detail != diplotype:
        result["minor_detail"] = minor_detail
    with open("aldy_cyp2d6_outside_call.tsv", "w") as f:
        f.write("CYP2D6\\t" + diplotype + "\\n")
else:
    result["status"] = "no_call"
    result["message"] = "Aldy did not produce a CYP2D6 diplotype call"
    if stderr_text.strip():
        result["stderr_tail"] = stderr_text[-2000:]
    with open("aldy_cyp2d6_outside_call.tsv", "w") as f:
        f.write("")

with open("aldy_cyp2d6_result.json", "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
    f.write("\\n")

PYPARSE
    exit 0
    """
}
