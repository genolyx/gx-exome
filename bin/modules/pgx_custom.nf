// ============================================================
// pgx_custom.nf — Extended PGx genotyping for non-PharmCAT genes
//
// Queries the annotated VCF at curated pharmacogenomic positions
// defined in pgx_custom_variants.tsv.  Produces a JSON result
// and a human-readable summary published alongside PharmCAT
// output under pgx/.
// ============================================================

process RUN_PGX_CUSTOM {
    tag "$sample_id"
    label 'pgx_custom'
    publishDir "${params.outdir}/pgx", mode: 'copy', pattern: 'pgx_custom_*.{json,txt}'

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(variants_tsv)

    output:
    path "pgx_custom_result.json", emit: result_json
    path "pgx_custom_summary.txt", emit: summary_txt

    script:
    """
    set +e
    export HOME=\$PWD
    pip install --user pysam 2>&1 | tail -3
    export PATH="\$PWD/.local/bin:\$PATH"

    export SAMPLE_ID="${sample_id}"
    export VCF_PATH="${vcf}"
    export VARIANTS_TSV="${variants_tsv}"

    python3 << 'PYGENOTYPE'
import json
import os
import sys

sample_id = os.environ["SAMPLE_ID"]
vcf_path = os.environ["VCF_PATH"]
tsv_path = os.environ["VARIANTS_TSV"]

try:
    import pysam
except ImportError:
    sys.stderr.write("pysam not available\\n")
    with open("pgx_custom_result.json", "w") as f:
        json.dump({"schema": "gx_exome_pgx_custom_v1", "sample_id": sample_id,
                    "error": "pysam_not_installed"}, f, indent=2)
        f.write("\\n")
    with open("pgx_custom_summary.txt", "w") as f:
        f.write(f"Extended PGx — {sample_id}\\nError: pysam not available\\n")
    sys.exit(0)

targets = []
with open(tsv_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\\t")
        if len(cols) < 8:
            continue
        targets.append({
            "gene": cols[0],
            "rsid": cols[1],
            "chrom": cols[2],
            "pos": int(cols[3]),
            "ref": cols[4],
            "alt": cols[5],
            "variant_name": cols[6],
            "clinical_significance": cols[7],
            "drugs": cols[8] if len(cols) > 8 else "",
            "evidence_level": cols[9] if len(cols) > 9 else "",
            "clinpgx_url": cols[10] if len(cols) > 10 else "",
        })

vcf = pysam.VariantFile(vcf_path)
csq_fields = []
for rec in vcf.header.records:
    if rec.key == "INFO" and rec.get("ID") == "CSQ":
        desc = rec.get("Description", "")
        if "Format:" in desc:
            fmt = desc.split("Format:")[1].strip().rstrip('"')
            csq_fields = fmt.split("|")
        break

existing_var_idx = csq_fields.index("Existing_variation") if "Existing_variation" in csq_fields else -1
symbol_idx = csq_fields.index("SYMBOL") if "SYMBOL" in csq_fields else -1
consequence_idx = csq_fields.index("Consequence") if "Consequence" in csq_fields else -1
hgvsp_idx = csq_fields.index("HGVSp") if "HGVSp" in csq_fields else -1
hgvsc_idx = csq_fields.index("HGVSc") if "HGVSc" in csq_fields else -1
sift_idx = csq_fields.index("SIFT") if "SIFT" in csq_fields else -1
polyphen_idx = csq_fields.index("PolyPhen") if "PolyPhen" in csq_fields else -1
clin_sig_idx = csq_fields.index("CLIN_SIG") if "CLIN_SIG" in csq_fields else -1
gnomad_af_idx = csq_fields.index("gnomADe_AF") if "gnomADe_AF" in csq_fields else -1

results = []

for t in targets:
    chrom = t["chrom"]
    pos_0 = t["pos"] - 1
    pos_1 = t["pos"]

    entry = {
        "gene": t["gene"],
        "rsid": t["rsid"],
        "variant_name": t["variant_name"],
        "clinical_significance": t["clinical_significance"],
        "drugs": t["drugs"],
        "evidence_level": t["evidence_level"],
        "clinpgx_url": t["clinpgx_url"],
        "chrom": chrom,
        "pos": t["pos"],
        "expected_ref": t["ref"],
        "expected_alt": t["alt"],
    }

    found = False
    try:
        for rec in vcf.fetch(chrom, pos_0, pos_1):
            if rec.pos != t["pos"]:
                continue

            gt = None
            for sample_name in rec.samples:
                gt_tuple = rec.samples[sample_name]["GT"]
                if gt_tuple is not None:
                    alleles = [rec.alleles[i] if i is not None else "." for i in gt_tuple]
                    gt = "/".join(alleles)
                break

            csq_str = rec.info.get("CSQ")
            rsid_match = False
            vep_gene = ""
            vep_consequence = ""
            vep_hgvsp = ""
            vep_hgvsc = ""
            vep_sift = ""
            vep_polyphen = ""
            vep_clin_sig = ""
            vep_gnomad_af = ""

            if csq_str:
                csq_list = csq_str if isinstance(csq_str, (list, tuple)) else [csq_str]
                for csq_entry in csq_list:
                    parts = csq_entry.split("|")
                    if existing_var_idx >= 0 and existing_var_idx < len(parts):
                        ev = parts[existing_var_idx]
                        if t["rsid"] in ev.split("&"):
                            rsid_match = True
                    if symbol_idx >= 0 and symbol_idx < len(parts):
                        vep_gene = parts[symbol_idx]
                    if consequence_idx >= 0 and consequence_idx < len(parts):
                        vep_consequence = parts[consequence_idx]
                    if hgvsp_idx >= 0 and hgvsp_idx < len(parts):
                        vep_hgvsp = parts[hgvsp_idx]
                    if hgvsc_idx >= 0 and hgvsc_idx < len(parts):
                        vep_hgvsc = parts[hgvsc_idx]
                    if sift_idx >= 0 and sift_idx < len(parts):
                        vep_sift = parts[sift_idx]
                    if polyphen_idx >= 0 and polyphen_idx < len(parts):
                        vep_polyphen = parts[polyphen_idx]
                    if clin_sig_idx >= 0 and clin_sig_idx < len(parts):
                        vep_clin_sig = parts[clin_sig_idx]
                    if gnomad_af_idx >= 0 and gnomad_af_idx < len(parts):
                        vep_gnomad_af = parts[gnomad_af_idx]
                    if rsid_match:
                        break

            is_ref = gt in (None, f"{rec.ref}/{rec.ref}") if rec else True
            is_het = gt is not None and len(set(gt.split("/"))) > 1
            is_hom_alt = gt is not None and not is_ref and not is_het

            if is_ref:
                zygosity = "homozygous_ref"
            elif is_het:
                zygosity = "heterozygous"
            elif is_hom_alt:
                zygosity = "homozygous_alt"
            else:
                zygosity = "unknown"

            entry["genotype"] = gt
            entry["zygosity"] = zygosity
            entry["rsid_confirmed"] = rsid_match
            entry["status"] = "variant_found"
            if vep_gene:
                entry["vep_gene"] = vep_gene
            if vep_consequence:
                entry["vep_consequence"] = vep_consequence
            if vep_hgvsp:
                entry["vep_hgvsp"] = vep_hgvsp
            if vep_hgvsc:
                entry["vep_hgvsc"] = vep_hgvsc
            if vep_sift:
                entry["vep_sift"] = vep_sift
            if vep_polyphen:
                entry["vep_polyphen"] = vep_polyphen
            if vep_clin_sig:
                entry["vep_clin_sig"] = vep_clin_sig.replace("&", ", ")
            if vep_gnomad_af:
                entry["gnomad_af"] = vep_gnomad_af
            found = True
            break
    except Exception as e:
        entry["status"] = "query_error"
        entry["error"] = str(e)
        results.append(entry)
        continue

    if not found:
        entry["genotype"] = f"{t['ref']}/{t['ref']}"
        entry["zygosity"] = "homozygous_ref"
        entry["status"] = "no_variant"
    results.append(entry)

vcf.close()

gene_groups = {}
for r in results:
    g = r["gene"]
    if g not in gene_groups:
        gene_groups[g] = []
    gene_groups[g].append(r)

output = {
    "schema": "gx_exome_pgx_custom_v1",
    "sample_id": sample_id,
    "total_targets": len(targets),
    "variants_found": sum(1 for r in results if r["status"] == "variant_found"),
    "genes": gene_groups,
}

with open("pgx_custom_result.json", "w") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
    f.write("\\n")

lines = []
lines.append(f"Extended PGx Panel — {sample_id}")
lines.append("=" * 72)
lines.append(f"Targets queried: {len(targets)}   Variants detected: {output['variants_found']}")
lines.append("")

carrier = []
normal = []

for gene in sorted(gene_groups.keys()):
    variants = gene_groups[gene]
    has_variant = any(v["status"] == "variant_found" and v["zygosity"] != "homozygous_ref" for v in variants)
    if has_variant:
        for v in variants:
            if v["status"] == "variant_found" and v["zygosity"] != "homozygous_ref":
                zyg = "HET" if v["zygosity"] == "heterozygous" else "HOM"
                drugs = v.get("drugs", "")
                ev = v.get("evidence_level", "")
                url = v.get("clinpgx_url", "")
                carrier.append({
                    "line1": f"  {v['gene']:12s}  {v['rsid']:16s}  {v['variant_name']:20s}  {v['genotype']:12s} ({zyg})",
                    "line2": f"{'':14s}  Drugs: {drugs}" if drugs else "",
                    "line3": f"{'':14s}  Evidence: {ev}   ClinPGx: {url}" if ev or url else "",
                })
    else:
        rsids = ", ".join(v["rsid"] for v in variants)
        url = variants[0].get("clinpgx_url", "")
        normal.append(f"  {gene:12s}  {rsids:40s}  {url}")

if carrier:
    lines.append("*** VARIANTS DETECTED ***")
    lines.append("-" * 90)
    for item in carrier:
        lines.append(item["line1"])
        if item["line2"]:
            lines.append(item["line2"])
        if item["line3"]:
            lines.append(item["line3"])
        lines.append("")

if normal:
    lines.append("Homozygous reference (no variant):")
    lines.append("-" * 90)
    for row in normal:
        lines.append(row)
    lines.append("")

lines.append("Evidence levels: 1A/1B = CPIC guideline or FDA label;  2A/2B = moderate evidence;")
lines.append("                 3 = single study / weak evidence;  4 = case report / in vitro only")
lines.append("")
lines.append("Source: ClinPGx (clinpgx.org) variant annotations + VCF genotype")
lines.append("Full details: pgx/pgx_custom_result.json")

with open("pgx_custom_summary.txt", "w") as f:
    f.write("\\n".join(lines) + "\\n")

PYGENOTYPE

    if [ ! -f pgx_custom_result.json ]; then
        echo '{"schema":"gx_exome_pgx_custom_v1","sample_id":"${sample_id}","error":"script_failed"}' > pgx_custom_result.json
    fi
    if [ ! -f pgx_custom_summary.txt ]; then
        echo "Extended PGx — ${sample_id}" > pgx_custom_summary.txt
        echo "Error: genotyping script failed" >> pgx_custom_summary.txt
    fi
    exit 0
    """
}
