process GENERATE_SUMMARY_REPORT {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path manta_vcfs
    path gcnv_vcfs
    path paraphase_jsons
    path hba_results
    path cyp21a2_results
    path eh_vcfs
    path smaca_results
    path intron_depth_reports
    path snapshots
    path bed_file
    path hba_paralog_tsvs

    output:
    path "*_summary_report.txt", emit: summary
    path "*_detailed_report.txt", emit: detailed

    script:
    """
    # --- Dependency Setup ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem || true
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
            && chmod +x micromamba_bin || true
    fi
    
    [ -f micromamba_bin ] && ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge python=3.9 pysam pandas -y || true
    export PATH=\$PWD/env/bin:\$PATH

    # --- Python Summary Script ---
    cat <<EOF > summarize.py
import os
import glob
import json
import gzip
import re
import ast

SMN_CNV_TOTAL_ASSUMED = ${params.smn_cnv_est_total_copies}

def silent_carrier_dup_g27134_allele_counts(raw):
    # g.27134T>G DUP_MARK cell (silent carrier dup context); not c.840 exon 7.
    if not raw or raw == '-':
        return None
    m = re.search(r'\\[.*\\]', raw, re.DOTALL)
    if not m:
        return None
    try:
        arr = ast.literal_eval(m.group(0))
    except Exception:
        return None
    if not isinstance(arr, list) or len(arr) != 4:
        return None
    if not all(isinstance(r, list) for r in arr):
        return None
    A = int(sum(arr[0]))
    C = int(sum(arr[1]))
    G = int(sum(arr[2]))
    T = int(sum(arr[3]))
    tg = T + G
    ct = C + T
    t_frac_tg = f"{T / tg:.3f}" if tg > 0 else "?"
    c_frac_ct = f"{C / ct:.3f}" if ct > 0 else "?"
    return {
        'A': A, 'C': C, 'G': G, 'T': T,
        'T_frac_TG': t_frac_tg,
        'C_frac_CT': c_frac_ct,
    }

def silent_carrier_summary_label(raw_val, acgt):
    # g.27134T>G: True if SMAca G read count > 0; else legacy "G" in text before [[...]].
    if raw_val == '-':
        return "False"
    raw_head = raw_val.split('[', 1)[0]
    if acgt is not None:
        g = acgt['G']
        call = "True" if g > 0 else "False"
        return f"{call} (g27134_G={g}; Raw: {raw_head})"
    if raw_head.strip().startswith('G'):
        return f"True (g27134_G=NA; Raw: {raw_head})"
    return f"False (g27134_G=NA; Raw: {raw_head})"

# --- DMD Target Mapping ---
def load_dmd_targets(bed_path):
    intervals = []
    try:
        with open(bed_path, 'r') as f:
            for line in f:
                # Assuming standard bed: chrom start end name ...
                if line.startswith('chrX') and 'DMD' in line:
                    parts = line.strip().split('\\t')
                    if len(parts) >= 3:
                        intervals.append((int(parts[1]), int(parts[2])))
    except Exception as e:
        print(f"Warning: Could not read BED file {bed_path}: {e}")
        return []

    if not intervals:
        return []

    intervals.sort()
    
    # Cluster intervals (merge if < 500bp apart, serving as 'exons')
    clusters = []
    if intervals:
        curr_start, curr_end = intervals[0]
        for i in range(1, len(intervals)):
            next_start, next_end = intervals[i]
            if next_start - curr_end < 500:
                curr_end = max(curr_end, next_end)
            else:
                clusters.append((curr_start, curr_end))
                curr_start, curr_end = next_start, next_end
        clusters.append((curr_start, curr_end))
    
    # Map to Exon Numbers (Highest Coord = Exon 1 for Minus Strand gene)
    # Sort Exons by coordinate descending
    clusters.sort(key=lambda x: x[0], reverse=True)
    
    dmd_exons = []
    for i, (start, end) in enumerate(clusters):
        exon_num = i + 1
        dmd_exons.append({'exon': exon_num, 'start': start, 'end': end})
        
    print(f"Loaded {len(dmd_exons)} DMD Exon Targets from BED.")
    return dmd_exons

dmd_exons = load_dmd_targets("${bed_file}")

def get_dmd_affected_exons(chrom, start, end):
    if chrom != 'chrX': return []
    affected = []
    # Simple overlap check
    # Variant: [start, end]
    # Exon: [e_start, e_end]
    for exon in dmd_exons:
        # Overlap condition: not (end < e_start or start > e_end)
        # Note: input pos might be single point, SVs have SVLEN.
        # We handle this logic in the caller, here we expect start/end ints.
        if not (end < exon['start'] or start > exon['end']):
            affected.append(exon['exon'])
    
    affected.sort()
    return affected

def format_exon_range(exon_list):
    if not exon_list: return ""
    # Simplify list: 1, 2, 3, 5, 6 -> 1-3, 5-6
    ranges = []
    if not exon_list: return ""
    
    start = exon_list[0]
    prev = exon_list[0]
    
    for x in exon_list[1:]:
        if x == prev + 1:
            prev = x
        else:
            if start == prev: ranges.append(str(start))
            else: ranges.append(f"{start}-{prev}")
            start = x
            prev = x
    
    if start == prev: ranges.append(str(start))
    else: ranges.append(f"{start}-{prev}")
    
    return "Exons " + ", ".join(ranges)

def parse_manta(files):
    results = {}
    print(f"Parsing {len(files)} Manta VCFs...")
    for f in files:
        sample = os.path.basename(f).replace('_manta.vcf.gz', '').replace('_manta.vcf', '')
        calls = []
        dmd_calls = []
        try:
            with gzip.open(f, 'rt') as vcf:
                for line in vcf:
                    if line.startswith('#'): continue
                    parts = line.split('\\t')
                    chrom = parts[0]
                    pos = int(parts[1])
                    info = parts[7]
                    
                    if 'SVTYPE=DEL' in info or 'SVTYPE=DUP' in info:
                        svlen = 0
                        for field in info.split(';'):
                            if field.startswith('SVLEN='):
                                svlen = int(field.split('=')[1])
                                break
                        
                        end = pos + abs(svlen) 
                        
                        # CLN3 Check (approx chr16:28477000-28503000)
                        if chrom == 'chr16' and 28470000 <= pos <= 28510000:
                             calls.append(f"DEL(CLN3-Like):{pos}:{svlen}")
                        
                        # DMD Check (chrX:31115794-33357558)
                        elif chrom == 'chrX' and 31115794 <= pos <= 33357558:
                             sv_type = "DEL" if "SVTYPE=DEL" in info else "DUP"
                             
                             # Identify Exons
                             # Manta DEL SVLEN is negative, so end is pos + abs(SVLEN)
                             # Manta DUP SVLEN is positive? Check spec. Usually SVLEN is length of event.
                             # For DUP, pos is start.
                             
                             affected = get_dmd_affected_exons('chrX', pos, end)
                             exon_str = format_exon_range(affected)
                             
                             dmd_calls.append(f"{sv_type}(DMD):{pos}:{svlen} | {exon_str}")
                             
                        # General large deletions (>1kb)
                        elif abs(svlen) > 1000:
                             calls.append(f"LARGE_DEL:{chrom}:{pos}:{svlen}")
            
            res_str = ""
            if calls: res_str += "; ".join(calls)
            if dmd_calls: 
                if res_str: res_str += " | "
                # Clean marker for Summary
                res_str += "DMD_SV:" + "; ".join(dmd_calls)
                
            results[sample] = res_str
        except Exception as e:
            results[sample] = f"Error: {e}"
    return results

def parse_gcnv(files):
    results = {}
    print(f"Parsing {len(files)} GCNV VCFs...")
    for f in files:
        sample = os.path.basename(f).replace('_cnv.vcf.gz', '')
        calls = []
        dmd_calls = []
        try:
            with gzip.open(f, 'rt') as vcf:
                for line in vcf:
                    if line.startswith('#'): continue
                    parts = line.split('\\t')
                    chrom = parts[0]
                    pos = int(parts[1])
                    fmt_def = parts[8].split(':')
                    if 'CN' in fmt_def:
                        cn_idx = fmt_def.index('CN')
                        fmt_val = parts[9].split(':')
                        cn = fmt_val[cn_idx]
                        
                        if cn != '2':
                             # CLN3 Check
                             if chrom == 'chr16' and 28470000 <= pos <= 28510000:
                                  calls.append(f"CLN3_CN={cn}")
                             # DMD Check (chrX:31115794-33357558)
                             elif chrom == 'chrX' and 31115794 <= pos <= 33357558:
                                  # GCNV calls are per-interval (exon) usually
                                  affected = get_dmd_affected_exons('chrX', pos, pos+10) # Point check
                                  exon_str = ""
                                  if affected:
                                      exon_str = f"(Exon {affected[0]})"
                                      
                                  dmd_calls.append(f"DMD_CN={cn}:{pos}{exon_str}")
        
            res_str = ""
            if calls: res_str += "; ".join(calls)
            if dmd_calls:
                if res_str: res_str += " | "
                res_str += "DMD_CNV:" + "; ".join(dmd_calls)
                
            results[sample] = res_str
        except Exception as e:
            results[sample] = f"Error: {e}"
    return results

def parse_paraphase(files):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_paraphase.json', '')
        try:
            with open(f) as json_file:
                # Check if empty
                if os.stat(f).st_size == 0:
                    results[sample] = "Error: JSON Empty"
                    continue
                    
                data = json.load(json_file)
                if not data:
                    results[sample] = "No calls"
                    continue
                summary = []
                for gene, info in data.items():
                    # Paraphase may emit null for a gene (JSON null) or final_haplotypes: null
                    if not isinstance(info, dict):
                        continue
                    # Generic CN
                    if 'total_cn' in info:
                        summary.append(f"{gene}_CN={info['total_cn']}")
                    
                    # SMN specific
                    if 'smn1_cn' in info:
                        summary.append(f"SMN1={info['smn1_cn']}")
                    if 'smn2_cn' in info:
                        summary.append(f"SMN2={info['smn2_cn']}")
                    
                    # GBA specific (Parsing Haplotypes for Gene vs Pseudogene)
                    if gene == 'GBA' and 'final_haplotypes' in info:
                        haps = info['final_haplotypes']
                        if haps is not None:
                            gba_count = 0
                            gbap1_count = 0
                            other_count = 0
                            hap_names = []
                            for h in haps:
                                # Handle both string and dict formats
                                if isinstance(h, str):
                                    h_name = h
                                else:
                                    h_name = h.get('haplotype', 'unknown')
                                hap_names.append(h_name)
                                lower_name = h_name.lower()
                                if 'gbap1' in lower_name:
                                    gbap1_count += 1
                                elif 'gba' in lower_name:
                                    gba_count += 1
                                else:
                                    other_count += 1
                            summary.append(f"GBA_Gene_CN={gba_count}")
                            summary.append(f"GBAP1_CN={gbap1_count}")
                            if other_count > 0:
                                summary.append(f"Other_CN={other_count}")
                            summary.append(f"Haps: {', '.join(hap_names)}")

                results[sample] = "; ".join(summary)
        except json.JSONDecodeError:
            results[sample] = "Error: Invalid JSON"
        except Exception as e:
            results[sample] = f"Error: {str(e)}"
    return results

def parse_fallback(files, label):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace(f'_{label}_fallback.txt', '')
        try:
            with open(f) as txt:
                content = txt.read()
                ratio = "N/A"
                status = "Unknown"
                
                for line in content.splitlines():
                    if "Ratio" in line: ratio = line.split(':')[-1].strip()
                    if "Status:" in line: status = line.split('Status:')[-1].strip()
                
                # Estimate Copy Number from Ratio
                # HBA Logic: Ratio ~2.0 = 4 Copies -> CN = Ratio * 2
                # CYP21A2 Logic: Ratio ~1.0 = 2 Copies -> CN = Ratio * 2
                est_cn = "?"
                try:
                    r_val = float(ratio)
                    # Simple rounding
                    est_cn = round(r_val * 2)
                except:
                    pass

                # Combine into concise string
                results[sample] = f"Est_CN={est_cn}|Ratio={ratio}|{status}"
        except:
            results[sample] = "Error"
    return results

def parse_hba_paralog_trace(files):
    # Lines from hba_paralog_pileup.py: HBA_PARALOG_CN_RATIO / HBA_PARALOG_CN_TRACE (comment lines in TSV).
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_hba_paralog.tsv', '')
        lines = []
        try:
            with open(f) as fh:
                for line in fh:
                    if line.startswith('# ') and 'HBA_PARALOG' in line:
                        lines.append(line[2:].strip())
        except Exception:
            pass
        results[sample] = lines
    return results

def parse_eh(files):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_eh.vcf', '')
        try:
            summary = []
            with open(f) as vcf:
                for line in vcf:
                    if line.startswith('#'): continue
                    parts = line.split('\\t')
                    if len(parts) >= 10:
                        info = parts[7]
                        fmt_keys = parts[8].split(':')
                        fmt_vals = parts[9].strip().split(':')
                        
                        repid_match = re.search(r'REPID=([^;]+)', info)
                        if repid_match:
                            gene = repid_match.group(1)
                            if 'REPCN' in fmt_keys:
                                idx = fmt_keys.index('REPCN')
                                cn = fmt_vals[idx]
                                summary.append(f"{gene}:{cn}")
            results[sample] = " | ".join(summary) if summary else "Ref"
        except:
            results[sample] = "Error"
    return results

def filter_smaca_for_detailed_report(line):
    # Detailed report: keep CN/SMAca + HGVS norm (G→C, A→T); drop verbose per-base ACGT.
    if not line or line == "No Data" or "|" not in line:
        return line
    parts = [p.strip() for p in line.split("|") if p.strip()]
    has_norm = any("_hgvs_norm_" in p for p in parts)
    keep = []
    for p in parts:
        if "=" not in p:
            continue
        k = p.split("=", 1)[0]
        if k in (
            "SMN1_cov_frac",
            "C_Ratio",
            "Cov(1,2)",
            "SMN1_CN_est",
            "SMN2_CN_est",
            "CNV_est_total_assumed",
            "SilentCarrier",
            "c840_pileup_mode",
            "c840_SMN1_depth_frac_of_pair",
        ):
            keep.append(p)
            continue
        if k.startswith("exon7_") or k.startswith("silent_g27134_"):
            keep.append(p)
            continue
        if k.startswith("SMAca_cov_frac_minus_c840"):
            keep.append(p)
            continue
        if k.startswith("c840_SMN1_hgvs_norm_") or k.startswith("c840_SMN2_hgvs_norm_"):
            keep.append(p)
            continue
        if k.startswith("c840_merged_hgvs_norm_"):
            keep.append(p)
            continue
        if has_norm:
            continue
        if k in (
            "c840_SMN1_hgvs_C",
            "c840_SMN1_hgvs_T",
            "c840_SMN1_hgvs_depth",
            "c840_SMN2_hgvs_C",
            "c840_SMN2_hgvs_T",
            "c840_SMN2_hgvs_depth",
        ):
            keep.append(p)
            continue
        if k in (
            "c840_merged_hgvs_C",
            "c840_merged_hgvs_T",
            "c840_merged_hgvs_depth",
            "c840_merged_hgvs_C_frac_CplusT_coding",
        ):
            keep.append(p)
    return "|".join(keep) if keep else line

def parse_smaca(files):
    results = {}
    extras = {}
    for f in files:
        sample = os.path.basename(f).replace('_smaca.txt', '')
        extras[sample] = []
        try:
             with open(f) as txt:
                 content = txt.read()
             lines = content.splitlines()
             # Prefer stable block appended by SMACA_RUN (daemon-friendly)
             for i, line in enumerate(lines):
                 if line.strip() == '##DARK_GENE_PIPELINE_SUMMARY' and i + 1 < len(lines):
                     results[sample] = lines[i + 1].strip()
                     j = i + 2
                     while j < len(lines) and lines[j].strip().startswith('##C840'):
                         extras[sample].append(lines[j].strip())
                         j += 1
                     break
             else:
                 header = []
                 data = []
                 for line in lines:
                     line = line.strip()
                     if line.startswith('# id') or line.startswith('id'):
                         header = [h.strip() for h in line.replace('#', '').strip().split(',')]
                     elif not line.startswith('#') and len(line) > 5:
                         if '|' in line:
                             data = [d.strip() for d in line.split('|')]
                         elif ',' in line:
                             data = [d.strip() for d in line.split(',')]
                 if header and data:
                     min_len = min(len(header), len(data))
                     row = {header[i]: data[i] for i in range(min_len)}
                     s1 = row.get('avg_cov_SMN1', '?')
                     s2 = row.get('avg_cov_SMN2', '?')
                     variant_col = 'g.27134T>G'
                     raw_val = row.get(variant_col, '-')
                     acgt = silent_carrier_dup_g27134_allele_counts(raw_val) if raw_val != '-' else None
                     silent_carrier = silent_carrier_summary_label(raw_val, acgt)
                     s1_val = None
                     s2_val = None
                     try:
                        s1_val = float(s1)
                        s1 = f"{s1_val:.2f}"
                     except Exception:
                        pass
                     try:
                        s2_val = float(s2)
                        s2 = f"{s2_val:.2f}"
                     except Exception:
                        pass
                     c_ratio = "?"
                     try:
                         if isinstance(s1_val, float) and isinstance(s2_val, float):
                             total_cov = s1_val + s2_val
                             if total_cov > 0:
                                 c_ratio = f"{s1_val / total_cov:.3f}"
                             else:
                                 c_ratio = "0.00"
                     except Exception:
                        pass
                     cn_est1 = "?"
                     cn_est2 = "?"
                     try:
                         if isinstance(s1_val, float) and isinstance(s2_val, float):
                             tot = s1_val + s2_val
                             if tot > 0:
                                 sm1_frac = s1_val / tot
                                 cn_est1 = int(round(SMN_CNV_TOTAL_ASSUMED * sm1_frac))
                                 cn_est2 = SMN_CNV_TOTAL_ASSUMED - cn_est1
                     except Exception:
                         pass
                     pi_b = row.get('Pi_b', '?')
                     cov_s1b = row.get('cov_SMN1_b', '?')
                     cov_s2b = row.get('cov_SMN2_b', '?')
                     exon7 = f"|exon7_c840_Pi_b={pi_b}|cov_SMN1_b={cov_s1b}|cov_SMN2_b={cov_s2b}"
                     silent_allelic = ""
                     if acgt:
                         silent_allelic = (
                             f"|silent_g27134_A={acgt['A']}|silent_g27134_C={acgt['C']}|silent_g27134_G={acgt['G']}|silent_g27134_T={acgt['T']}"
                             f"|silent_g27134_T_frac_TplusG={acgt['T_frac_TG']}|silent_g27134_C_frac_CplusT={acgt['C_frac_CT']}"
                         )
                     else:
                         silent_allelic = "|silent_g27134_alleles=NA"
                     results[sample] = (
                         f"SMN1_cov_frac={c_ratio}|C_Ratio={c_ratio}"
                         f"|Cov(1,2)={s1},{s2}|SMN1_CN_est={cn_est1}|SMN2_CN_est={cn_est2}|CNV_est_total_assumed={SMN_CNV_TOTAL_ASSUMED}"
                         f"|SilentCarrier={silent_carrier}{exon7}{silent_allelic}"
                         f"|c840_pileup=NA|c840_SMN1_depth_frac_of_pair=NA|SMAca_cov_frac_minus_c840_depth_frac=NA"
                     )
                 else:
                     results[sample] = "Parsing Failed"
        except Exception as e:
            results[sample] = f"Error: {e}"
    return results, extras

def parse_intron_verify(files):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_intron_depth.txt', '')
        try:
             with open(f) as txt:
                 content = txt.read()
                 warnings = []
                 if "WARNING:" in content:
                     for line in content.splitlines():
                         if "WARNING:" in line:
                             warnings.append(line.strip())
                 
                 results[sample] = warnings
        except:
             results[sample] = []
    return results

# Main
manta_files = glob.glob("*_manta.vcf.gz")
gcnv_files = glob.glob("*_cnv.vcf.gz")
para_files = glob.glob("*_paraphase.json")
hba_files = glob.glob("*_hba_fallback.txt")
hba_paralog_files = glob.glob("*_hba_paralog.tsv")
cyp_files = glob.glob("*_cyp21a2_fallback.txt")
eh_files = glob.glob("*_eh.vcf")
sma_files = glob.glob("*_smaca.txt")
ver_files = glob.glob("*_intron_depth.txt")

m_res = parse_manta(manta_files)
g_res = parse_gcnv(gcnv_files)
p_res = parse_paraphase(para_files)
h_res = parse_fallback(hba_files, "hba")
hba_paralog_lines = parse_hba_paralog_trace(hba_paralog_files)
c_res = parse_fallback(cyp_files, "cyp21a2")
e_res = parse_eh(eh_files)
s_res, smaca_extras = parse_smaca(sma_files)
v_res = parse_intron_verify(ver_files)

all_samples = sorted(
    set(m_res.keys()) | set(g_res.keys()) | set(p_res.keys()) | set(hba_paralog_lines.keys())
)

# Generate Individual Reports
header = "Sample\\tParaphase(SMN/GBA)\\tSMAca\\tFragileX(FMR1)\\tHBA_Ratio\\tCYP21A2_Ratio\\tLarge_SVs(Manta/GCNV)\\tQC_Warnings"

for s in all_samples:
    # 1. Individual Summary Report
    with open(f"{s}_summary_report.txt", "w") as out:
        out.write(header + "\\n")
        
        sv_combined = []
        
        # Clean up Manta/GCNV for summary (remove detailed DMD tags, kept in Detailed Report)
        m_val = m_res.get(s, "")
        if "DMD_SV:" in m_val:
             pass
        if m_val: sv_combined.append(m_val)
        
        g_val = g_res.get(s, "")
        if g_val: sv_combined.append(g_val)
        
        sv_str = "; ".join(sv_combined) if sv_combined else "-"
        
        # Warnings
        warns = v_res.get(s, [])
        warn_str = " | ".join(warns) if warns else "PASS"

        line = f"{s}\\t{p_res.get(s, '-')}\\t{s_res.get(s, '-')}\\t{e_res.get(s, '-')}\\t{h_res.get(s, '-')}\\t{c_res.get(s, '-')}\\t{sv_str}\\t{warn_str}"
        out.write(line + "\\n")

    # 2. Individual Detailed Report
    with open(f"{s}_detailed_report.txt", "w") as out:
        out.write(f"SAMPLE: {s}\\n")
        out.write("="*40 + "\\n")
        
        # Verification Warnings
        warns = v_res.get(s, [])
        if warns:
             out.write("!!! QUALITY WARNINGS !!!\\n")
             for w in warns:
                 out.write(f"  {w}\\n")
             out.write("\\n")
        
        out.write("PARAPHASE RESULTS (SMN1/2, GBA - Breakdown):\\n")
        out.write(f"  {p_res.get(s, 'No Data')}\\n\\n")

        out.write("SMAca CHECK (exon7/c840: Pi_b cov_SMN*_b; silent carrier: g.27134 SilentCarrier):\\n")
        smaca_val = filter_smaca_for_detailed_report(s_res.get(s, 'No Data'))
        smaca_fmt = smaca_val.replace('|', '\\n  ')
        out.write(f"  {smaca_fmt}\\n")
        for _ex in smaca_extras.get(s, []):
            out.write(f"  {_ex}\\n")
        out.write("\\n")

        out.write("HBA ANALYSIS (Alpha Thalassemia - Dosage):\\n")
        hba_val = h_res.get(s, 'No Data')
        hba_fmt = hba_val.replace('|', '\\n  ')
        out.write(f"  {hba_fmt}\\n")
        for _hl in hba_paralog_lines.get(s, []):
            out.write(f"  {_hl}\\n")
        out.write("\\n")
        
        out.write("CYP21A2 ANALYSIS (CAH - Dosage):\\n")
        cyp_val = c_res.get(s, 'No Data')
        cyp_fmt = cyp_val.replace('|', '\\n  ')
        out.write(f"  {cyp_fmt}\\n\\n")

        out.write("EXPANSION HUNTER (Fragile X / FMR1):\\n")
        out.write(f"  {e_res.get(s, 'No Data')}\\n\\n")
        
        # DMD SECTION
        out.write("DMD ANALYSIS (ChrX:31.1M-33.3M):\\n")
        dmd_found = False
        
        m_val = m_res.get(s, "")
        if "DMD_SV:" in m_val:
            # Extract DMD part
            parts = m_val.split("DMD_SV:")
            if len(parts) > 1:
                dmd_hits = parts[1].split(" | ")[0] # Logic if it's in middle
                out.write(f"  [Manta] {dmd_hits}\\n")
                dmd_found = True
                
        g_val = g_res.get(s, "")
        if "DMD_CNV:" in g_val:
             parts = g_val.split("DMD_CNV:")
             if len(parts) > 1:
                 dmd_hits = parts[1].split(" | ")[0]
                 out.write(f"  [GCNV]  {dmd_hits}\\n")
                 dmd_found = True
        
        if not dmd_found:
            out.write("  No CNVs/SVs detected in DMD region.\\n")
        out.write("\\n")

        out.write("LARGE SVs (Manta/GCNV - Rest):\\n")
        sv_combined = []
        if m_res.get(s): sv_combined.append(f"Manta: {m_res[s]}")
        if g_res.get(s): sv_combined.append(f"GCNV: {g_res[s]}")
        
        if sv_combined:
            for item in sv_combined:
                out.write(f"  {item}\\n")
        else:
            out.write("  No large events detected.\\n")

        out.write("\\n" + "-"*40 + "\\n\\n")
        
print("Summary and Detailed reports generated.")
EOF

    python3 summarize.py
    """
}
