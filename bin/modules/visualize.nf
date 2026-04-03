process GENERATE_VISUAL_EVIDENCE {
    tag "$sample_id"
    label 'visualize'
    publishDir "${params.outdir}/snapshots", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(vcf_idx), path(eh_realigned_bam), path(eh_realigned_bai), path(eh_json), path(smn_unified_bam), path(smn_unified_bai)
    path eh_images
    path ref_fasta
    path ref_fai
    path viz_env
    path gtf
    path gtf_idx
    
    output:
    path "*_visual_report.html", emit: snapshots

    script:
    """
    # --- Shared Environment Setup ---
    # Use the pre-built environment passed from PREPARE_VIZ_RESOURCES
    export PATH=\$PWD/viz_env/bin:\$PATH

    # Symlink Indices
    if [ ! -f "${bam}.bai" ]; then ln -s ${bai} ${bam}.bai; fi
    if [ ! -f "${ref_fasta}.fai" ]; then ln -s ${ref_fai} ${ref_fasta}.fai; fi
    # Real files (not symlinks) so igv-reports embeds them; symlinks often yield an empty EH track in the HTML bundle
    cp -f ${eh_realigned_bam} FMR1_eh_realigned.bam
    cp -f ${eh_realigned_bai} FMR1_eh_realigned.bam.bai

    # Single SMA panel: SMN1 NM_000344 c.840 pileup (GRCh38 ref C vs alt T); ±120 bp around smn_c840_sm1_grch38_1bp
    printf "chr5\\t${params.smn_c840_sm1_grch38_1bp - 120}\\t${params.smn_c840_sm1_grch38_1bp + 120}\\tSMA_SMN1_c840\\n" > regions.bed
    printf "chr16\\t173300\\t173800\\tHb_Constant_Spring\\n" >> regions.bed

    printf '##gff-version 3\\n' > sma_landmarks.gff3
    printf 'chr5\\tgx_exome\\tregion\\t%d\\t%d\\t.\\t+\\t.\\tID=SMN1_c840;Name=NM_000344_c.840_refC_altT\\n' ${params.smn_c840_sm1_grch38_1bp} ${params.smn_c840_sm1_grch38_1bp} >> sma_landmarks.gff3

    # SMN unified BAM (same bwa-mem2 + mini-ref as SMACA_RUN) → project chr5_local onto hg38 chr5 for IGV
    export SMN_START=${params.smn_unified_region_start}
    export PRIMARY_BAM="${bam}"
    export UNIFIED_BAM="${smn_unified_bam}"
    python3 << 'PYSMN'
import os
import pysam

# Mini-ref slice starts at SMN_START (1-based g.); SAM POS is 1-based; BAM ref_start is 0-based.
off = int(os.environ["SMN_START"]) - 1
primary = os.environ["PRIMARY_BAM"]
unified = os.environ["UNIFIED_BAM"]

in_pri = pysam.AlignmentFile(primary, "rb")
out = pysam.AlignmentFile("smn_combined_hg38.bam", "wb", template=in_pri)
if out.get_tid("chr5") < 0:
    raise SystemExit("chr5 missing from primary BAM header")

in_u = pysam.AlignmentFile(unified, "rb")
hdr_u = in_u.header


def read_to_sam_line(r):
    # pysam: older builds use to_string(header); newer use to_string() only
    try:
        return r.to_string(hdr_u)
    except TypeError:
        return r.to_string()


for read in in_u:
    if read.is_unmapped:
        continue
    # Cannot mutate reference_id on reads from in_u — pysam validates against unified header (chr5_local only).
    # Round-trip through SAM with hg38 RNAME/POS, then parse with primary header.
    s = read_to_sam_line(read)
    parts = s.split("\t")
    if len(parts) < 11:
        continue
    if parts[2] == "chr5_local":
        parts[2] = "chr5"
    parts[3] = str(int(parts[3]) + off)
    if parts[6] == "chr5_local":
        parts[6] = "="
    if parts[6] == "=":
        parts[7] = str(int(parts[7]) + off)
    line = "\t".join(parts)
    try:
        lifted = pysam.AlignedSegment.fromstring(line, out.header)
    except AttributeError:
        lifted = pysam.AlignedSegment.from_string(line, out.header)
    except ValueError as e:
        raise SystemExit(f"SMN lift failed for read {parts[0]}: {e}") from e
    out.write(lifted)

in_u.close()
out.close()
in_pri.close()
PYSMN
    samtools sort -o smn_combined_final.bam smn_combined_hg38.bam
    samtools index smn_combined_final.bam

    # RefSeq GTF is already provided (staged as refGene.sorted.gtf.gz or similar)
    # Symlink to strict names for config
    ln -s ${gtf} refGene.gtf.gz
    ln -s ${gtf_idx} refGene.gtf.gz.tbi

    # --- Prepare Fragile X Image Track (Popup) ---
    # --- Prepare Fragile X Image Track (Popup) ---
    fx_track_logic=""
    eh_img=\$(find . -maxdepth 1 -name "${sample_id}*.svg" | head -n 1)
    
    if [ -n "\$eh_img" ]; then
        echo "Found Fragile X Image, embedding in GFF..."
        
        # 1. Create GFF with Button Trigger (Short string, safe)
        cat <<EOF > create_fx_track.py
with open('fx_image.gff3', 'w') as f:
    f.write('##gff-version 3\\n')
    import urllib.parse
    
    # JavaScript to show the hidden div and scroll to it
    js_action = "document.getElementById('fmr1_eh_panel').scrollIntoView({behavior: 'smooth'});"
    
    btn_html = f'<button onclick="{js_action}" style="background:orange; color:white; font-weight:bold; padding:5px; cursor:pointer;">JUMP TO FMR1 PANEL</button>'
    encoded_note = urllib.parse.quote(btn_html)
    
    f.write(f'chrX\\tREViewer\\tTrigger\\t147911600\\t147912450\\t.\\t.\\t.\\tID=FX1;Name=Click to View;Note={encoded_note}\\n')
EOF
        python3 create_fx_track.py
        
        if [ -f "fx_image.gff3" ]; then
            touch fx_track_ready
        fi
    fi

    # --- Build Track Config (Dynamic JSON) ---
    cat <<EOF > build_config.py
import json
import os

tracks = []

# 1. BAM Track
tracks.append({
    "type": "alignment",
    "format": "bam",
    "url": "${bam}",
    "indexURL": "${bam}.bai",
    "name": "${sample_id}",
    "displayMode": "EXPANDED", 
    "colorBy": "base",
    "height": 500,
    "minMappingQuality": 0,
    "samplingDepth": 1000000,
    "showSoftClips": True,
    "visibilityWindow": -1
})

# 1b. SMN1 c.840 — unified mini-ref (same as SMAca); use with Color by base to see T vs ref C at landmark
tracks.append({
    "type": "alignment",
    "format": "bam",
    "url": "smn_combined_final.bam",
    "indexURL": "smn_combined_final.bam.bai",
    "name": "SMN1 c.840 pileup (unified ref, SMAca-matched)",
    "displayMode": "EXPANDED",
    "colorBy": "base",
    "height": 550,
    "minMappingQuality": 0,
    "samplingDepth": 1000000,
    "showSoftClips": True,
    "visibilityWindow": -1
})

# 1c. SMN1 c.840 position (ref C vs alt T)
tracks.append({
    "type": "annotation",
    "format": "gff3",
    "url": "sma_landmarks.gff3",
    "name": "SMN1 c.840 (NM_000344)",
    "displayMode": "EXPANDED",
    "color": "darkgreen",
    "height": 60,
    "visibilityWindow": -1
})

# 1d. FMR1 — ExpansionHunter repeat-aware realignment
tracks.append({
    "type": "alignment",
    "format": "bam",
    "url": "FMR1_eh_realigned.bam",
    "indexURL": "FMR1_eh_realigned.bam.bai",
    "name": "FMR1 (ExpansionHunter realigned)",
    "displayMode": "EXPANDED",
    "colorBy": "base",
    "height": 500,
    "minMappingQuality": 0,
    "samplingDepth": 1000000,
    "showSoftClips": True,
    "visibilityWindow": -1
})

# 3. RefSeq Genes (Detailed)
tracks.append({
    "type": "annotation",
    "format": "gtf",
    "url": "refGene.gtf.gz",
    "indexURL": "refGene.gtf.gz.tbi",
    "name": "RefSeq Genes (Detail)",
    "displayMode": "EXPANDED",
    "visibilityWindow": -1
})

# 4. RefSeq Genes (Structure Only)
tracks.append({
    "type": "annotation",
    "format": "gtf",
    "url": "refGene.gtf.gz",
    "indexURL": "refGene.gtf.gz.tbi",
    "name": "RefSeq Genes (Structure)",
    "displayMode": "SQUISHED",
    "color": "darkblue",
    "visibilityWindow": -1
})

# 5. Fragile X Image Track (Conditional)
if os.path.exists('fx_track_ready') and os.path.exists('fx_image.gff3'):
    tracks.append({
        "type": "annotation",
        "format": "gff3",
        "url": "fx_image.gff3",
        "name": "Fragile X (Graph Trigger)",
        "displayMode": "EXPANDED",
        "color": "orange",
        "height": 40,
        "visibilityWindow": -1
    })


# 3. VCF Track
vcf_file = "${vcf}"
if vcf_file and vcf_file != "null" and os.path.exists(vcf_file):
    # Ensure index link
    if not os.path.exists(vcf_file + ".tbi") and os.path.exists("${vcf_idx}"):
         os.symlink("${vcf_idx}", vcf_file + ".tbi")
         
    tracks.append({
        "type": "variant",
        "format": "vcf",
        "url": vcf_file,
        "indexURL": vcf_file + ".tbi",
        "name": "Called Variants",
        "displayMode": "SQUISHED",
        "color": "purple",
        "visibilityWindow": -1
    })

with open("track_config.json", "w") as f:
    json.dump(tracks, f, indent=2)
EOF

    python3 build_config.py

    echo "Generating IGV Report..."
    
    create_report regions.bed \\
        --fasta ${ref_fasta} \\
        --genome hg38 \\
        --track-config track_config.json \\
        --flanking 1500 \\
        --output temp_report.html \\
        --title "Visual Evidence: ${sample_id}"

    # --- Alpha thalassemia review note (IGV report footer) ---
    python3 -c "
html = open('temp_report.html').read()
alpha_box = '''<div id='alpha_thal_review_note' style='margin-top:24px;padding:16px;border:1px solid #2c5282;border-radius:8px;background:#f0f7ff;max-width:900px;'>
<h3 style='margin:0 0 8px 0;'>Alpha thalassemia (review)</h3>
<p style='margin:0;line-height:1.5;color:#1a365d;'>The fallback analysis runs Freebayes on <code>chr16:172000–178000</code> with ploidy 4 mainly for Hb Constant Spring at a fixed position — it is not a full α-thal SNV/indels genotype caller for all carriers.</p>
</div>'''
html = html.replace('</body>', alpha_box + chr(10) + '</body>')
open('temp_report.html', 'w').write(html)
"

    # --- Fragile X: EH CGG summary + REViewer graph (visible; primary view for repeat count) ---
    export EH_JSON_PATH="${eh_json}"
    export EH_SVG_PATH="\${eh_img:-}"
    python3 << 'PYFMR1'
import html as html_mod
import json
import os
import re
from pathlib import Path

eh_path = Path(os.environ["EH_JSON_PATH"])
svg_path = (os.environ.get("EH_SVG_PATH") or "").strip()
html_path = Path("temp_report.html")
text = html_path.read_text(encoding="utf-8", errors="replace")


def format_genotype(gt):
    if gt is None:
        return "—"
    if isinstance(gt, (list, tuple)):
        return "/".join(str(x) for x in gt)
    if isinstance(gt, dict):
        return json.dumps(gt)
    return str(gt)


def extract_rows(locus):
    rows = []
    for v in locus.get("Variants") or []:
        if v.get("VariantType") != "Repeat":
            continue
        rows.append(
            {
                "unit": v.get("RepeatUnit") or "",
                "genotype": v.get("Genotype"),
                "ci": v.get("GenotypeConfidenceInterval"),
                "ref": v.get("ReferenceRegion"),
            }
        )
    return rows


def fmr1_block(data):
    lr = data.get("LocusResults") or {}
    for key, locus in lr.items():
        lid = str(locus.get("LocusId") or key)
        if "FMR1" not in lid and "FMR1" not in str(key):
            continue
        rows = extract_rows(locus)
        return {
            "locus_id": lid,
            "coverage": locus.get("Coverage"),
            "allele_count": locus.get("AlleleCount"),
            "rows": rows,
        }
    for key, locus in lr.items():
        rows = extract_rows(locus)
        if not rows:
            continue
        if any((r["unit"] or "").upper() == "CGG" for r in rows):
            lid = str(locus.get("LocusId") or key)
            return {
                "locus_id": lid,
                "coverage": locus.get("Coverage"),
                "allele_count": locus.get("AlleleCount"),
                "rows": rows,
            }
    return None

summary_inner = ""
if eh_path.exists():
    try:
        with open(eh_path, encoding="utf-8") as f:
            data = json.load(f)
        blk = fmr1_block(data)
        if blk and blk["rows"]:
            parts = [
                "<table style='border-collapse:collapse;margin-top:8px;font-size:14px;'>",
                "<tr><th style='text-align:left;padding:4px 12px 4px 0;border-bottom:1px solid #ccc;'>Motif</th>"
                "<th style='text-align:left;padding:4px 12px;border-bottom:1px solid #ccc;'>CGG repeats (genotype)</th>"
                "<th style='text-align:left;padding:4px 0;border-bottom:1px solid #ccc;'>95% CI</th></tr>",
            ]
            for r in blk["rows"]:
                gt_s = format_genotype(r["genotype"])
                ci = r["ci"]
                ci_s = html_mod.escape(str(ci)) if ci is not None else "—"
                parts.append(
                    "<tr><td style='padding:6px 12px 6px 0;'>"
                    + html_mod.escape(str(r["unit"]))
                    + "</td><td style='padding:6px 12px;font-weight:600;'>"
                    + html_mod.escape(gt_s)
                    + "</td><td style='padding:6px 0;'>"
                    + ci_s
                    + "</td></tr>"
                )
            parts.append("</table>")
            cov = blk["coverage"]
            extra = ""
            if cov is not None:
                extra += f"<p style='margin:8px 0 0 0;color:#4a5568;font-size:13px;'>Locus coverage (EH estimate): {html_mod.escape(str(cov))}</p>"
            summary_inner = "".join(parts) + extra
        elif blk:
            summary_inner = "<p style='margin:0;color:#744210;'>FMR1 found in JSON but no repeat variant rows.</p>"
        else:
            summary_inner = "<p style='margin:0;color:#744210;'>No FMR1 entry under <code>LocusResults</code> in ExpansionHunter JSON.</p>"
    except Exception as e:
        summary_inner = f"<p style='margin:0;color:#c53030;'>Could not parse ExpansionHunter JSON: {html_mod.escape(str(e))}</p>"
else:
    summary_inner = "<p style='margin:0;color:#744210;'>ExpansionHunter JSON not found.</p>"

svg_block = ""
if svg_path and Path(svg_path).exists():
    svg_block = (
        "<h3 style='margin:20px 0 10px 0;font-size:16px;'>REViewer (read graph)</h3>"
        "<div style='overflow-x:auto;border:1px solid #e2e8f0;border-radius:6px;padding:8px;background:#fafafa;'>"
        + Path(svg_path).read_text(encoding="utf-8", errors="replace")
        + "</div>"
    )
else:
    svg_block = "<p style='margin:12px 0 0 0;color:#4a5568;font-size:13px;'>REViewer SVG was not generated (e.g. FMR1 absent from EH VCF). Use the table above and the IGV panel <code>FragileX_FMR1</code> for the primary BAM.</p>"

# Groovy wraps this script in double-triple-quotes; Python panel string must use f''' not f-double-triple-quote
panel = f'''<div id='fmr1_eh_panel' style='margin:0 0 24px 0;padding:20px;border:2px solid #dd6b20;border-radius:10px;background:linear-gradient(180deg,#fffaf0 0%,#ffffff 40%);max-width:1100px;'>
<h2 style='margin:0 0 8px 0;color:#9c4221;font-size:20px;'>Fragile X — FMR1 (CGG)</h2>
<p style='margin:0 0 12px 0;line-height:1.5;color:#2d3748;font-size:14px;'>
<b>ExpansionHunter</b> estimates <b>CGG repeat copy number</b> below (not inferred from the primary exome BAM pileup). The read graph shows how EH/REViewer represent spanning and in-repeat reads.
</p>
{summary_inner}
{svg_block}
</div>'''

m = re.search(r"(<body[^>]*>)", text, re.I)
if m:
    text = text[: m.end()] + panel + text[m.end() :]
else:
    text = panel + text
html_path.write_text(text, encoding="utf-8")
PYFMR1

    mv temp_report.html ${sample_id}_visual_report.html

    # Verify output
    if [ ! -f "${sample_id}_visual_report.html" ]; then
        echo "Error: Report generation failed."
        exit 1
    fi
    """
}
