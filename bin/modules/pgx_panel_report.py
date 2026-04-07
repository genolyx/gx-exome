#!/usr/bin/env python3
"""
pgx_panel_report.py — Generate combined PGx panel HTML report.

Merges PharmCAT phenotype/report JSON (CPIC/DPWG genes) with the
extended-panel custom JSON into a single interactive HTML table.

Usage:
    python3 pgx_panel_report.py \
        --sample-id SAMPLE \
        --phenotype-json path/to/phenotype.json \
        --report-json path/to/report.json \
        --custom-json path/to/pgx_custom_result.json \
        --output pgx_panel_report.html
"""
import argparse
import json
import os
import html as html_mod
from datetime import datetime


RISK_FUNCTIONS = {"no function", "decreased function", "unfavorable response allele"}


def load_drug_map(report_path):
    drug_map = {}
    if not report_path or not os.path.isfile(report_path):
        return drug_map
    with open(report_path) as f:
        rdata = json.load(f)
    for src in ("CPIC", "DPWG"):
        for gn, gv in rdata.get("genes", {}).get(src, {}).items():
            if gn in drug_map:
                continue
            rd = gv.get("relatedDrugs", [])
            if isinstance(rd, list) and rd:
                if isinstance(rd[0], dict):
                    names = [d.get("name", "") for d in rd[:5]]
                else:
                    names = [str(d) for d in rd[:5]]
                if len(rd) > 5:
                    names.append(f"+{len(rd) - 5} more")
                drug_map[gn] = ", ".join(names)
    return drug_map


def load_pharmcat_rows(phenotype_path, drug_map):
    rows = []
    if not os.path.isfile(phenotype_path) or os.path.getsize(phenotype_path) == 0:
        return rows
    with open(phenotype_path) as f:
        pdata = json.load(f)
    for source in ("CPIC", "DPWG"):
        reports = pdata.get("geneReports", {}).get(source, {})
        if not isinstance(reports, dict):
            continue
        for gene_name in sorted(reports.keys()):
            if any(r["gene"] == gene_name for r in rows):
                continue
            gd = reports[gene_name]
            if not isinstance(gd, dict):
                continue
            call_src = gd.get("callSource", "")
            rec_dips = gd.get("recommendationDiplotypes", [])
            drugs_str = drug_map.get(gene_name, "\u2014")
            if not rec_dips:
                rows.append({
                    "gene": gene_name, "variant": "\u2014", "genotype": "No call",
                    "phenotype": "\u2014", "drugs": drugs_str, "evidence": "CPIC/DPWG",
                    "zygosity": "\u2014", "source": "PharmCAT", "status": "no_call",
                    "url": f"https://www.clinpgx.org/hgnc/{gene_name}",
                })
                continue
            dip = rec_dips[0]
            a1 = dip.get("allele1") or {}
            a2 = dip.get("allele2") or {}
            n1, n2 = a1.get("name", ""), a2.get("name", "")
            fn1 = (a1.get("function") or "").strip()
            fn2 = (a2.get("function") or "").strip()
            phenotypes = dip.get("phenotypes", [])
            activity = dip.get("activityScore")
            diplotype = f"{n1}/{n2}" if n2 else n1
            phen = ", ".join(p for p in phenotypes if p) if phenotypes else "\u2014"
            functions = [fn for fn in (fn1, fn2) if fn]
            has_risk = any(fn.lower() in RISK_FUNCTIONS for fn in functions)
            src_label = "Aldy (via PharmCAT)" if call_src == "OUTSIDE" else "PharmCAT"
            status = "actionable" if has_risk else ("normal" if phen != "\u2014" else "no_call")
            rows.append({
                "gene": gene_name, "variant": diplotype, "genotype": diplotype,
                "phenotype": phen, "drugs": drugs_str, "evidence": "CPIC/DPWG",
                "zygosity": "Diplotype", "source": src_label, "status": status,
                "url": f"https://www.clinpgx.org/hgnc/{gene_name}",
                "activity": str(activity) if activity and str(activity) != "None" else "",
            })
    return rows


def load_custom_rows(custom_path):
    rows = []
    if not os.path.isfile(custom_path) or os.path.getsize(custom_path) == 0:
        return rows
    with open(custom_path) as f:
        cdata = json.load(f)
    genes = cdata.get("genes", {})
    for gene_name in sorted(genes.keys()):
        for v in genes[gene_name]:
            is_variant = v.get("status") == "variant_found" and v.get("zygosity") != "homozygous_ref"
            zyg = v.get("zygosity", "homozygous_ref")
            if zyg == "heterozygous":
                zyg_label = "HET"
            elif zyg == "homozygous_alt":
                zyg_label = "HOM ALT"
            else:
                zyg_label = "REF"
            rows.append({
                "gene": v.get("gene", gene_name),
                "variant": v.get("variant_name", "\u2014"),
                "genotype": v.get("genotype", "\u2014"),
                "phenotype": v.get("clinical_significance", "\u2014"),
                "drugs": v.get("drugs", "\u2014"),
                "evidence": v.get("evidence_level", "\u2014"),
                "zygosity": zyg_label,
                "source": "Extended Panel",
                "status": "variant" if is_variant else "reference",
                "url": v.get("clinpgx_url", f"https://www.clinpgx.org/hgnc/{gene_name}"),
                "rsid": v.get("rsid", ""),
            })
    return rows


def esc(s):
    return html_mod.escape(str(s)) if s else ""


def status_badge(status):
    colors = {
        "actionable": ("#c53030", "#fff5f5", "Actionable"),
        "variant": ("#b7791f", "#fffff0", "Variant"),
        "normal": ("#276749", "#f0fff4", "Normal"),
        "reference": ("#718096", "#f7fafc", "Reference"),
        "no_call": ("#a0aec0", "#edf2f7", "No Call"),
    }
    color, bg, label = colors.get(status, ("#718096", "#f7fafc", status))
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:12px;'
        f'font-size:11px;font-weight:600;color:{color};background:{bg};'
        f'border:1px solid {color}30;">{label}</span>'
    )


def evidence_badge(ev):
    ev_str = str(ev).strip()
    if not ev_str or ev_str == "\u2014":
        return ""
    if "1A" in ev_str or "1B" in ev_str:
        cls = "ev-1A"
    elif "2A" in ev_str or "2B" in ev_str:
        cls = "ev-2A"
    elif "3" in ev_str:
        cls = "ev-3"
    elif "4" in ev_str:
        cls = "ev-4"
    else:
        cls = "ev-3"
    return f'<span class="evidence-badge {cls}">{esc(ev_str)}</span>'


def drug_cell(drugs):
    if not drugs or drugs == "\u2014":
        return '<span style="color:#a0aec0;">\u2014</span>'
    parts = drugs.split(", ")
    shown = ", ".join(parts[:3])
    if len(parts) > 3:
        shown += f' <span style="color:#718096;font-size:11px;">(+{len(parts) - 3})</span>'
    return shown


def render_table(rows, section_id, columns):
    parts = [f'<table id="{section_id}">', "<thead><tr>"]
    for idx, col in enumerate(columns):
        parts.append(f"<th onclick=\"sortTable('{section_id}', {idx})\">{col}</th>")
    parts.append("</tr></thead><tbody>")
    for r in rows:
        row_class = ""
        if r["status"] == "actionable":
            row_class = "actionable"
        elif r["status"] == "variant":
            row_class = "variant-row"
        parts.append(
            f'<tr class="{row_class} pgx-row" '
            f'data-status="{r["status"]}" data-gene="{esc(r["gene"])}">'
        )
        parts.append(
            f'<td><span class="gene-name">'
            f'<a href="{esc(r["url"])}" target="_blank">{esc(r["gene"])}</a></span>'
        )
        if r.get("rsid"):
            parts.append(f'<br><span class="rsid">{esc(r["rsid"])}</span>')
        parts.append("</td>")
        parts.append(f"<td>{esc(r['variant'])}</td>")
        parts.append(f"<td><code>{esc(r['genotype'])}</code></td>")
        parts.append(f"<td>{esc(r['zygosity'])}</td>")
        parts.append(f"<td>{esc(r['phenotype'])}</td>")
        parts.append(f"<td>{drug_cell(r['drugs'])}</td>")
        parts.append(f"<td>{evidence_badge(r['evidence'])}</td>")
        parts.append(f"<td>{status_badge(r['status'])}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return "\n".join(parts)


CSS = """\
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif; background: #f7fafc; color: #2d3748; }
.header { background: linear-gradient(135deg, #2c5282, #2b6cb0); color: white; padding: 24px 32px; }
.header h1 { font-size: 22px; font-weight: 700; margin-bottom: 4px; }
.header .subtitle { font-size: 14px; opacity: 0.85; }
.stats { display: flex; gap: 16px; margin-top: 16px; flex-wrap: wrap; }
.stat { background: rgba(255,255,255,0.15); padding: 8px 16px; border-radius: 8px; text-align: center; }
.stat .num { font-size: 24px; font-weight: 700; }
.stat .label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; opacity: 0.8; }
.container { max-width: 1400px; margin: 0 auto; padding: 20px; }
.controls { background: white; padding: 12px 16px; border-radius: 8px; margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
.controls input { padding: 8px 12px; border: 1px solid #e2e8f0; border-radius: 6px; font-size: 14px; width: 260px; }
.controls select { padding: 8px 12px; border: 1px solid #e2e8f0; border-radius: 6px; font-size: 14px; background: white; }
.controls label { font-size: 13px; color: #4a5568; display: flex; align-items: center; gap: 4px; cursor: pointer; }
.section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); margin-bottom: 20px; overflow: hidden; }
.section-header { padding: 14px 20px; font-size: 15px; font-weight: 700; border-bottom: 2px solid #e2e8f0; display: flex; justify-content: space-between; align-items: center; cursor: pointer; user-select: none; }
.section-header:hover { background: #f7fafc; }
.section-header .toggle { font-size: 12px; color: #718096; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
thead th { background: #edf2f7; padding: 10px 12px; text-align: left; font-weight: 600; color: #4a5568; border-bottom: 2px solid #e2e8f0; position: sticky; top: 0; cursor: pointer; white-space: nowrap; }
thead th:hover { background: #e2e8f0; }
tbody td { padding: 10px 12px; border-bottom: 1px solid #f0f0f0; vertical-align: top; }
tbody tr:hover { background: #f7fafc; }
tbody tr.actionable { background: #fff5f5; }
tbody tr.actionable:hover { background: #fee; }
tbody tr.variant-row { background: #fffff0; }
tbody tr.variant-row:hover { background: #fefcbf; }
a { color: #3182ce; text-decoration: none; }
a:hover { text-decoration: underline; }
.gene-name { font-weight: 700; font-size: 14px; }
.rsid { font-size: 11px; color: #718096; }
.evidence-badge { display: inline-block; padding: 1px 6px; border-radius: 4px; font-size: 11px; font-weight: 600; }
.ev-1A, .ev-1B { background: #c6f6d5; color: #22543d; }
.ev-2A, .ev-2B { background: #fefcbf; color: #744210; }
.ev-3 { background: #e2e8f0; color: #4a5568; }
.ev-4 { background: #f7fafc; color: #a0aec0; }
.footer { text-align: center; padding: 16px; font-size: 12px; color: #a0aec0; }
.hidden { display: none; }
"""

JS = """\
function toggleSection(id) {
  var el = document.getElementById(id);
  var toggle = document.getElementById(id + '-toggle');
  if (el.classList.contains('hidden')) {
    el.classList.remove('hidden');
    toggle.innerHTML = '&#9660;';
  } else {
    el.classList.add('hidden');
    toggle.innerHTML = '&#9654;';
  }
}
function filterTable() {
  var search = document.getElementById('searchBox').value.toLowerCase();
  var statusF = document.getElementById('statusFilter').value;
  var hideRef = document.getElementById('hideRef').checked;
  var rows = document.querySelectorAll('.pgx-row');
  for (var i = 0; i < rows.length; i++) {
    var row = rows[i];
    var text = row.textContent.toLowerCase();
    var status = row.getAttribute('data-status');
    var show = true;
    if (search && text.indexOf(search) === -1) show = false;
    if (statusF && status !== statusF) show = false;
    if (hideRef && (status === 'reference' || status === 'normal')) show = false;
    row.style.display = show ? '' : 'none';
  }
}
function sortTable(tableId, colIdx) {
  var table = document.getElementById(tableId);
  var tbody = table.querySelector('tbody');
  var rows = Array.prototype.slice.call(tbody.querySelectorAll('tr'));
  var dir = table.getAttribute('data-sort-dir') === 'asc' ? 'desc' : 'asc';
  table.setAttribute('data-sort-dir', dir);
  rows.sort(function(a, b) {
    var aText = (a.children[colIdx] ? a.children[colIdx].textContent : '').trim().toLowerCase();
    var bText = (b.children[colIdx] ? b.children[colIdx].textContent : '').trim().toLowerCase();
    return dir === 'asc' ? aText.localeCompare(bText) : bText.localeCompare(aText);
  });
  for (var i = 0; i < rows.length; i++) { tbody.appendChild(rows[i]); }
}
"""


def build_html(sample_id, pharmcat_rows, custom_rows):
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    total_genes = len(set(r["gene"] for r in pharmcat_rows)) + len(set(r["gene"] for r in custom_rows))
    actionable_count = sum(1 for r in pharmcat_rows if r["status"] == "actionable")
    variant_count = sum(1 for r in custom_rows if r["status"] == "variant")
    ext_genes = len(set(r["gene"] for r in custom_rows))

    pc_table = render_table(
        pharmcat_rows, "pharmcat-table",
        ["Gene", "Diplotype", "Genotype", "Type", "Phenotype", "Drugs", "Evidence", "Status"],
    )
    cx_table = render_table(
        custom_rows, "custom-table",
        ["Gene", "Variant", "Genotype", "Zygosity", "Significance", "Drugs", "Evidence", "Status"],
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PGx Panel Report \u2014 {esc(sample_id)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="header">
  <h1>Pharmacogenomics (PGx) Panel Report</h1>
  <div class="subtitle">{esc(sample_id)} &mdash; {now}</div>
  <div class="stats">
    <div class="stat"><div class="num">{total_genes}</div><div class="label">Genes Tested</div></div>
    <div class="stat"><div class="num">{actionable_count}</div><div class="label">Actionable (PharmCAT)</div></div>
    <div class="stat"><div class="num">{variant_count}</div><div class="label">Variants (Extended)</div></div>
    <div class="stat"><div class="num">{len(pharmcat_rows)}</div><div class="label">PharmCAT Genes</div></div>
    <div class="stat"><div class="num">{ext_genes}</div><div class="label">Extended Panel Genes</div></div>
  </div>
</div>
<div class="container">
<div class="controls">
  <input type="text" id="searchBox" placeholder="Search gene, drug, variant..." oninput="filterTable()">
  <select id="statusFilter" onchange="filterTable()">
    <option value="">All Results</option>
    <option value="actionable">Actionable Only</option>
    <option value="variant">Variants Only</option>
    <option value="normal">Normal Only</option>
    <option value="reference">Reference Only</option>
  </select>
  <label><input type="checkbox" id="hideRef" onchange="filterTable()"> Hide reference/normal</label>
</div>
<div class="section">
  <div class="section-header" onclick="toggleSection('pharmcat-body')">
    <span>PharmCAT \u2014 CPIC/DPWG Guided Genes ({len(pharmcat_rows)} genes)</span>
    <span class="toggle" id="pharmcat-body-toggle">&#9660;</span>
  </div>
  <div id="pharmcat-body">{pc_table}</div>
</div>
<div class="section">
  <div class="section-header" onclick="toggleSection('custom-body')">
    <span>Extended PGx Panel \u2014 Curated Variant Genotyping ({ext_genes} genes, {len(custom_rows)} variants)</span>
    <span class="toggle" id="custom-body-toggle">&#9660;</span>
  </div>
  <div id="custom-body">{cx_table}</div>
</div>
<div class="footer">
  Evidence: <span class="evidence-badge ev-1A">1A/1B</span> CPIC guideline / FDA label &nbsp;
  <span class="evidence-badge ev-2A">2A/2B</span> Moderate evidence &nbsp;
  <span class="evidence-badge ev-3">3</span> Limited evidence &nbsp;
  <span class="evidence-badge ev-4">4</span> Case report / in-vitro
  <br><br>
  Sources: <a href="https://www.clinpgx.org" target="_blank">ClinPGx</a> (PharmGKB + CPIC) &bull;
  <a href="https://pharmcat.clinpgx.org" target="_blank">PharmCAT</a>
  &bull; Generated {now}
</div>
</div>
<script>{JS}</script>
</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(description="Generate PGx panel HTML report")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--phenotype-json", required=True)
    parser.add_argument("--report-json", default="")
    parser.add_argument("--custom-json", required=True)
    parser.add_argument("--output", default="pgx_panel_report.html")
    args = parser.parse_args()

    drug_map = load_drug_map(args.report_json)
    pharmcat_rows = load_pharmcat_rows(args.phenotype_json, drug_map)
    custom_rows = load_custom_rows(args.custom_json)
    html_content = build_html(args.sample_id, pharmcat_rows, custom_rows)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(html_content)


if __name__ == "__main__":
    main()
