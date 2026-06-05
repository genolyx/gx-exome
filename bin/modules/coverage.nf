process DEPTH_ANALYSIS {
    tag "$sample_id"
    label 'paraphase' // Reuse existing env with samtools+mosdepth
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed
    path dark_genes_bed

    output:
    path "${sample_id}_target_coverage.txt", emit: cov
    tuple val(sample_id), path("${sample_id}_target_coverage.txt"), emit: cov_tuple // For dependency joining
    path "${sample_id}_qc_metrics.txt", emit: qc
    path "${sample_id}_intron_depth.txt", emit: intron_report

    script:
    """
    export TMPDIR=\$PWD
    export HOME=\$PWD

    if ! command -v samtools >/dev/null 2>&1 || ! command -v mosdepth >/dev/null 2>&1; then
        export PIP_CACHE_DIR=\$PWD/.cache/pip
        export XDG_CACHE_HOME=\$PWD/.cache/xdg
        export PATH=\$PWD:\$PATH
        export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
        export MAMBA_ROOT_PREFIX=\$PWD/micromamba
        mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX \$PIP_CACHE_DIR \$XDG_CACHE_HOME
        wget -q --no-check-certificate https://curl.se/ca/cacert.pem
        export SSL_CERT_FILE=\$PWD/cacert.pem
        export MAMBA_SSL_VERIFY=false
        if [ ! -f "micromamba_bin" ]; then
            wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
                && chmod +x micromamba_bin
        fi
        if [ -f micromamba_bin ] && [ ! -x ./env/bin/samtools ]; then
            ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 mosdepth=0.3.3 -y
        fi
        [ -x ./env/bin/samtools ] && [ -x ./env/bin/mosdepth ] || { echo "ERROR: samtools/mosdepth env missing"; exit 1; }
        export PATH=\$PWD/env/bin:\$PATH
    fi

    ST=\$(command -v samtools)
    MD=\$(command -v mosdepth)
    
    # --- Analysis ---

    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    \$ST bedcov ${bed} ${bam} > ${sample_id}_target_coverage.txt

    \$MD --by ${bed} --thresholds 20,50,100 ${sample_id} ${bam}

    # 3. Parse Summary for User
    echo "Sample: ${sample_id}" > ${sample_id}_qc_metrics.txt
    echo "QC Metrics (on targets)" >> ${sample_id}_qc_metrics.txt
    
    MEAN_DEPTH=\$(grep "total_region" ${sample_id}.mosdepth.summary.txt | awk '{print \$4}')
    echo "Mean Depth: \$MEAN_DEPTH" >> ${sample_id}_qc_metrics.txt

    PCT_20X=\$(grep "^total.*	20	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')
    PCT_50X=\$(grep "^total.*	50	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')
    PCT_100X=\$(grep "^total.*	100	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')

    echo "% Bases > 20x: \$PCT_20X" >> ${sample_id}_qc_metrics.txt
    echo "% Bases > 50x: \$PCT_50X" >> ${sample_id}_qc_metrics.txt
    echo "% Bases > 100x: \$PCT_100X" >> ${sample_id}_qc_metrics.txt
    
    
    # --- 4. Intron Depth Verification (Merged) ---
    echo "Region\\tMean_Depth\\tMax_Depth\\tPercent_Covered" > ${sample_id}_intron_depth.txt
    
    \$ST depth -a -b ${dark_genes_bed} ${bam} | \\
    awk -F'\\t' -v bed="${dark_genes_bed}" '
        BEGIN {
            while ((getline line < bed) > 0) {
                split(line, f, "\\t")
                key = f[1] ":" f[2] "-" f[3]
                name[key] = f[4]
                order[++n] = key
            }
        }
        {
            key = \$1 ":" (NR > 0 ? "" : "")
            region = \$1 ":" \$2
            for (k in name) {
                split(k, rf, /[:-]/)
                if (\$1 == rf[1] && \$2 >= rf[2]+1 && \$2 <= rf[3]) {
                    sum[k] += \$3
                    count[k]++
                    if (\$3 > 0) covered[k]++
                    if (\$3 > max[k]) max[k] = \$3
                    break
                }
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                k = order[i]
                if (count[k] > 0)
                    printf "%s\\t%.2f\\t%d\\t%.2f\\n", name[k], sum[k]/count[k], max[k], (covered[k]/count[k])*100
                else
                    printf "%s\\t0\\t0\\t0\\n", name[k]
            }
        }
    ' >> ${sample_id}_intron_depth.txt

    # Same windows as HBA1_HBA2_Region + SMN1_SMN2_Region rows above (dark_genes_plus.bed) — mean depth, not a separate hardcoded locus/median.
    echo "\\n[Dark-gene SMA / HBA depth QC]" >> ${sample_id}_intron_depth.txt
    ALPHA_MEAN=\$(awk -F'\\t' '\$1=="HBA1_HBA2_Region" {print \$2+0}' ${sample_id}_intron_depth.txt)
    SMA_MEAN=\$(awk -F'\\t' '\$1=="SMN1_SMN2_Region" {print \$2+0}' ${sample_id}_intron_depth.txt)
    echo "HBA1_HBA2_Region mean depth (alpha): \$ALPHA_MEAN" >> ${sample_id}_intron_depth.txt
    echo "SMN1_SMN2_Region mean depth (SMA/SMN): \$SMA_MEAN" >> ${sample_id}_intron_depth.txt

    WARN_MIN=${params.depth_warning_min_dark}
    INC_SMN=${params.depth_warning_include_smn}
    if [ "\$INC_SMN" != "true" ]; then
        echo "Note: depth QC threshold uses HBA only (depth_warning_include_smn=false); SMN/SMA rely on Paraphase/SMAca." >> ${sample_id}_intron_depth.txt
    fi
    if [ "\$WARN_MIN" = "0" ] || [ "\$WARN_MIN" = "0.0" ]; then
        echo "Note: depth QC warning disabled (depth_warning_min_dark=0)." >> ${sample_id}_intron_depth.txt
    else
        if [ "\$INC_SMN" = "true" ]; then
            TRIGGER=\$(awk -v a="\$ALPHA_MEAN" -v s="\$SMA_MEAN" -v t="\$WARN_MIN" 'BEGIN{print ((a+0 < t+0) || (s+0 < t+0)) ? 1 : 0}')
        else
            TRIGGER=\$(awk -v a="\$ALPHA_MEAN" -v t="\$WARN_MIN" 'BEGIN{print (a+0 < t+0) ? 1 : 0}')
        fi
        if [ "\$TRIGGER" = "1" ]; then
            if [ "\$INC_SMN" = "true" ]; then
                FAIL_WHY=\$(awk -v a="\$ALPHA_MEAN" -v s="\$SMA_MEAN" -v t="\$WARN_MIN" 'BEGIN{
                    fa=(a+0<t); fs=(s+0<t);
                    if (fa && fs) print "Both regions fail.";
                    else if (fa) print "HBA1_HBA2_Region fails only (SMN1_SMN2_Region is OK).";
                    else if (fs) print "SMN1_SMN2_Region fails only (HBA1_HBA2_Region is OK).";
                    else print "";
                }')
                echo "WARNING: Dark-gene mean depth below \$WARN_MIN x — \$FAIL_WHY Values: HBA1_HBA2_Region=\$ALPHA_MEAN, SMN1_SMN2_Region=\$SMA_MEAN (same rows as Region table; not mosdepth panel-wide mean)." >> ${sample_id}_intron_depth.txt
                echo "WARNING: SMN window often looks much lower than exome-wide depth (SMN1/SMN2 homology); use Paraphase/SMAca for SMA context." >> ${sample_id}_intron_depth.txt
            else
                echo "WARNING: HBA1_HBA2_Region mean depth below \$WARN_MIN x (value=\$ALPHA_MEAN). SMN depth not used for this gate (depth_warning_include_smn=false)." >> ${sample_id}_intron_depth.txt
            fi
        fi
    fi
    """
}

// -------------------------------------------------------
// Genome-wide mosdepth per-base (tabix) + optional gene-level TSV for service-daemon
// Publishes under qc/ (same tree as SAMTOOLS_BAM_STATS) for extract_qc_summary / gene coverage.
// -------------------------------------------------------
process MOSDEPTH_PER_BASE_QC {
    tag "$sample_id"
    label 'paraphase'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    // Note: path(..., mode: 'copy') is not valid inside tuple(); plain path() relies on NXF_DATA_DIR mount for BED.
    tuple val(sample_id), path(bam), path(bai), path(gene_bed), val(build_gene_tsv)

    output:
    path "${sample_id}.mosdepth.per-base.bed.gz", emit: per_base
    path "${sample_id}.mosdepth.per-base.bed.gz.tbi", emit: per_base_index
    path "${sample_id}_gene_coverage.tsv", emit: gene_tsv

    script:
    """
    set -e
    export TMPDIR=\$PWD
    export HOME=\$PWD
    export PIP_CACHE_DIR=\$PWD/.cache/pip
    export XDG_CACHE_HOME=\$PWD/.cache/xdg
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX \$PIP_CACHE_DIR \$XDG_CACHE_HOME

    wget -q --no-check-certificate https://curl.se/ca/cacert.pem || true
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \\
            && chmod +x micromamba_bin || true
    fi

    if [ -f micromamba_bin ] && [ ! -x ./env/bin/samtools ]; then
        ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 mosdepth=0.3.3 python=3.11 -y
    fi
    [ -x ./env/bin/samtools ] && [ -x ./env/bin/mosdepth ] && [ -x ./env/bin/tabix ] && [ -x ./env/bin/python3 ] || { echo "ERROR: samtools/mosdepth/tabix/python3 env missing"; exit 1; }

    export PATH=\$PWD/env/bin:\$PATH

    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    PREFIX="${sample_id}.mosdepth"
    mosdepth -t ${task.cpus} "\$PREFIX" ${bam}

    PB="\${PREFIX}.per-base.bed.gz"
    if [ ! -f "\$PB" ]; then
        echo "ERROR: expected \$PB from mosdepth"
        ls -la
        exit 1
    fi

    if tabix -p bed "\$PB" 2>/dev/null; then
        echo "tabix index OK"
    else
        echo "WARN: tabix failed; re-bgzip for block gzip"
        zcat "\$PB" | bgzip -c > "\$PB.tmp" && mv "\$PB.tmp" "\$PB"
        tabix -p bed "\$PB"
    fi

    # Decode embedded copy of modules/mosdepth_gene_coverage_tsv.py (avoids staging .py into paraphase container).
    if [ "${build_gene_tsv}" = "true" ]; then
        base64 -d <<'MGC_B64' > mosdepth_gene_coverage_tsv.py
IyEvdXNyL2Jpbi9lbnYgcHl0aG9uMwojIFNvdXJjZSBvZiB0cnV0aCBmb3IgdGhlIGJhc2U2NCBw
YXlsb2FkIGluIGNvdmVyYWdlLm5mIChwcm9jZXNzIE1PU0RFUFRIX1BFUl9CQVNFX1FDKTsga2Vl
cCB0aGVtIGluIHN5bmMuCiIiIgpCdWlsZCBwZXItZ2VuZSBjb3ZlcmFnZSBUU1YgZnJvbSBtb3Nk
ZXB0aCBwZXItYmFzZSAodGFiaXgtaW5kZXhlZCkgKyBwYW5lbCBCRUQuCgpCRUQgY29sdW1uIDQg
aXMgdHJlYXRlZCBhcyB0aGUgZ2VuZSBrZXkgKEhHTkMtc3R5bGUgdGFyZ2V0cykuIEludGVydmFs
cyBwZXIgZ2VuZSBhcmUKbWVyZ2VkIG9uIGVhY2ggY2hyb21vc29tZTsgcGVyLWJhc2Ugcm93cyBh
cmUgY2xpcHBlZCB0byB0aG9zZSBpbnRlcnZhbHMuCgpPdXRwdXQgY29sdW1ucyBtYXRjaCBkYWVt
b24tZnJpZW5kbHkgbmFtZXM6IGdlbmUsIG1lYW5fY292ZXJhZ2UsIHBjdF9iYXNlc18xMHgsIHBj
dF9iYXNlc18yMHguCiIiIgpmcm9tIF9fZnV0dXJlX18gaW1wb3J0IGFubm90YXRpb25zCgppbXBv
cnQgYXJncGFyc2UKaW1wb3J0IG9zCmltcG9ydCByZQppbXBvcnQgc3VicHJvY2VzcwppbXBvcnQg
c3lzCmZyb20gY29sbGVjdGlvbnMgaW1wb3J0IGRlZmF1bHRkaWN0CmZyb20gdHlwaW5nIGltcG9y
dCBEaWN0LCBMaXN0LCBUdXBsZQoKCkludGVydmFsID0gVHVwbGVbc3RyLCBpbnQsIGludF0gICMg
Y2hyb20sIHN0YXJ0LCBlbmQgKDAtYmFzZWQgaGFsZi1vcGVuKQoKIyBQcmVmZXIgdGFiaXggZnJv
bSB0aGUgc2FtZSBjb25kYSBlbnYgYXMgdGhpcyBpbnRlcnByZXRlciAoaGFuZGxlcyBjb250YWlu
ZXJzCiMgd2hlcmUgUEFUSCBtYXkgbm90IGluY2x1ZGUgdGhlIG1pY3JvbWFtYmEtY3JlYXRlZCBl
bnYvYmluIGF0IFB5dGhvbiBsYXVuY2ggdGltZSkuCl9UQUJJWF9DQU5ESURBVEUgPSBvcy5wYXRo
LmpvaW4ob3MucGF0aC5kaXJuYW1lKG9zLnBhdGguYWJzcGF0aChzeXMuZXhlY3V0YWJsZSkpLCAi
dGFiaXgiKQpfVEFCSVggPSBfVEFCSVhfQ0FORElEQVRFIGlmIG9zLnBhdGguaXNmaWxlKF9UQUJJ
WF9DQU5ESURBVEUpIGVsc2UgInRhYml4IgoKX0dFTkVfU1lNID0gcmUuY29tcGlsZShyIl5bQS1a
YS16XVtBLVphLXowLTktXXswLDI0fSQiKQpfU0tJUF9UT0tfUFJFRklYRVMgPSAoCiAgICAiRU5T
VCIsCiAgICAiTk1fIiwKICAgICJOUl8iLAogICAgIlhNXyIsCiAgICAiWFJfIiwKICAgICJFTlNH
IiwKICAgICJDQ0RTIiwKICAgICJDTElOSUQiLAogICAgIkxPQyIsCikKCgpkZWYgX2dlbmVfa2V5
c19mcm9tX2JlZF9jb2w0KG5hbWU6IHN0cikgLT4gTGlzdFtzdHJdOgogICAgIiIiCiAgICBIR05D
LWxpa2Ugc3ltYm9scyBpbiBjb2x1bW4gNC4gVHdpc3QvdmVuZG9yIEJFRHMgdXNlIGBgUEFIO05N
Xy4uLjtFTlNULi4uYGA7IHBsYWluIEJFRHMgdXNlIGBgUEFIYGAuCiAgICAiIiIKICAgIHJhdyA9
IChuYW1lIG9yICIiKS5zdHJpcCgpCiAgICBpZiBub3QgcmF3OgogICAgICAgIHJldHVybiBbXQog
ICAga2V5czogTGlzdFtzdHJdID0gW10KICAgIGZvciBjaHVuayBpbiByYXcucmVwbGFjZSgiLCIs
ICI7IikucmVwbGFjZSgifCIsICI7Iikuc3BsaXQoIjsiKToKICAgICAgICB0ID0gY2h1bmsuc3Ry
aXAoKQogICAgICAgIGlmIG5vdCB0OgogICAgICAgICAgICBjb250aW51ZQogICAgICAgIHVsID0g
dC51cHBlcigpCiAgICAgICAgaWYgYW55KHVsLnN0YXJ0c3dpdGgocCkgZm9yIHAgaW4gX1NLSVBf
VE9LX1BSRUZJWEVTKToKICAgICAgICAgICAgY29udGludWUKICAgICAgICBpZiBfR0VORV9TWU0u
bWF0Y2godCk6CiAgICAgICAgICAgIGtleXMuYXBwZW5kKHVsKQogICAgaWYga2V5czoKICAgICAg
ICByZXR1cm4ga2V5cwogICAgaWYgX0dFTkVfU1lNLm1hdGNoKHJhdyk6CiAgICAgICAgcmV0dXJu
IFtyYXcudXBwZXIoKV0KICAgIHJldHVybiBbXQoKCmRlZiBfcGFyc2VfYmVkKHBhdGg6IHN0cikg
LT4gRGljdFtzdHIsIExpc3RbSW50ZXJ2YWxdXToKICAgICIiIkdyb3VwIGludGVydmFscyBieSBI
R05DIHN5bWJvbCBkZXJpdmVkIGZyb20gY29sdW1uIDQuIiIiCiAgICBieV9nZW5lOiBEaWN0W3N0
ciwgTGlzdFtJbnRlcnZhbF1dID0gZGVmYXVsdGRpY3QobGlzdCkKICAgIHdpdGggb3BlbihwYXRo
LCBlbmNvZGluZz0idXRmLTgiLCBlcnJvcnM9InJlcGxhY2UiKSBhcyBmOgogICAgICAgIGZvciBs
aW5lIGluIGY6CiAgICAgICAgICAgIGxpbmUgPSBsaW5lLnN0cmlwKCkKICAgICAgICAgICAgaWYg
bm90IGxpbmUgb3IgbGluZS5zdGFydHN3aXRoKCIjIikgb3IgbGluZS5zdGFydHN3aXRoKCJ0cmFj
ayAiKToKICAgICAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgICAgIHBhcnRzID0gbGluZS5z
cGxpdCgiXHQiKQogICAgICAgICAgICBpZiBsZW4ocGFydHMpIDwgNDoKICAgICAgICAgICAgICAg
IGNvbnRpbnVlCiAgICAgICAgICAgIGNocm9tLCBzLCBlLCBuYW1lID0gcGFydHNbMF0sIHBhcnRz
WzFdLCBwYXJ0c1syXSwgcGFydHNbM10uc3RyaXAoKQogICAgICAgICAgICBpZiBub3QgbmFtZToK
ICAgICAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgICAgIHRyeToKICAgICAgICAgICAgICAg
IHN0YXJ0X2kgPSBpbnQocykKICAgICAgICAgICAgICAgIGVuZF9pID0gaW50KGUpCiAgICAgICAg
ICAgIGV4Y2VwdCBWYWx1ZUVycm9yOgogICAgICAgICAgICAgICAgY29udGludWUKICAgICAgICAg
ICAgaWYgZW5kX2kgPD0gc3RhcnRfaToKICAgICAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAg
ICAgIGZvciBna2V5IGluIF9nZW5lX2tleXNfZnJvbV9iZWRfY29sNChuYW1lKToKICAgICAgICAg
ICAgICAgIGJ5X2dlbmVbZ2tleV0uYXBwZW5kKChjaHJvbSwgc3RhcnRfaSwgZW5kX2kpKQogICAg
cmV0dXJuIGJ5X2dlbmUKCgpkZWYgX21lcmdlX2ludGVydmFscyhpbnRlcnZhbHM6IExpc3RbSW50
ZXJ2YWxdKSAtPiBMaXN0W0ludGVydmFsXToKICAgIGJ5X2Nocm9tOiBEaWN0W3N0ciwgTGlzdFtU
dXBsZVtpbnQsIGludF1dXSA9IGRlZmF1bHRkaWN0KGxpc3QpCiAgICBmb3IgY2hyb20sIHMsIGUg
aW4gaW50ZXJ2YWxzOgogICAgICAgIGJ5X2Nocm9tW2Nocm9tXS5hcHBlbmQoKHMsIGUpKQogICAg
bWVyZ2VkOiBMaXN0W0ludGVydmFsXSA9IFtdCiAgICBmb3IgY2hyb20sIGl2cyBpbiBzb3J0ZWQo
YnlfY2hyb20uaXRlbXMoKSwga2V5PWxhbWJkYSB4OiB4WzBdKToKICAgICAgICBpdnMuc29ydChr
ZXk9bGFtYmRhIHg6ICh4WzBdLCB4WzFdKSkKICAgICAgICBjdXJfcywgY3VyX2UgPSBpdnNbMF0K
ICAgICAgICBmb3IgcywgZSBpbiBpdnNbMTpdOgogICAgICAgICAgICBpZiBzIDw9IGN1cl9lOgog
ICAgICAgICAgICAgICAgY3VyX2UgPSBtYXgoY3VyX2UsIGUpCiAgICAgICAgICAgIGVsc2U6CiAg
ICAgICAgICAgICAgICBtZXJnZWQuYXBwZW5kKChjaHJvbSwgY3VyX3MsIGN1cl9lKSkKICAgICAg
ICAgICAgICAgIGN1cl9zLCBjdXJfZSA9IHMsIGUKICAgICAgICBtZXJnZWQuYXBwZW5kKChjaHJv
bSwgY3VyX3MsIGN1cl9lKSkKICAgIHJldHVybiBtZXJnZWQKCgpkZWYgX2FsdF9jaHJvbShjOiBz
dHIpIC0+IHN0cjoKICAgIGlmIGMuc3RhcnRzd2l0aCgiY2hyIik6CiAgICAgICAgcmV0dXJuIGNb
MzpdIGlmIGxlbihjKSA+IDMgZWxzZSBjCiAgICByZXR1cm4gImNociIgKyBjCgoKZGVmIF90YWJp
eF9saW5lcyhiZ3o6IHN0ciwgY2hyb206IHN0ciwgc3RhcnQ6IGludCwgZW5kOiBpbnQpIC0+IExp
c3Rbc3RyXToKICAgICIiIlF1ZXJ5IHBlci1iYXNlIEJFRDsgdXNlIDAtYmFzZWQgaGFsZi1vcGVu
IGludGVydmFscyAobWF0Y2hlcyBCRUQgc3RvcmFnZSkuIiIiCiAgICBmb3IgY2ggaW4gKGNocm9t
LCBfYWx0X2Nocm9tKGNocm9tKSk6CiAgICAgICAgciA9IGYie2NofTp7c3RhcnR9LXtlbmR9Igog
ICAgICAgIHRyeToKICAgICAgICAgICAgb3V0ID0gc3VicHJvY2Vzcy5jaGVja19vdXRwdXQoCiAg
ICAgICAgICAgICAgICBbX1RBQklYLCAiLS16ZXJvLWJhc2VkIiwgYmd6LCByXSwKICAgICAgICAg
ICAgICAgIHN0ZGVycj1zdWJwcm9jZXNzLkRFVk5VTEwsCiAgICAgICAgICAgICAgICB0ZXh0PVRy
dWUsCiAgICAgICAgICAgICkKICAgICAgICBleGNlcHQgc3VicHJvY2Vzcy5DYWxsZWRQcm9jZXNz
RXJyb3I6CiAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgaWYgb3V0LnN0cmlwKCk6CiAgICAg
ICAgICAgIHJldHVybiBvdXQuc3BsaXRsaW5lcygpCiAgICByZXR1cm4gW10KCgpkZWYgX2FjY3Vt
dWxhdGVfaW50ZXJ2YWwoCiAgICBiZ3o6IHN0ciwKICAgIGNocm9tOiBzdHIsCiAgICBnczogaW50
LAogICAgZ2U6IGludCwKICAgIHRvdGFsX2Jhc2VzOiBpbnQsCiAgICBzdW1fZGVwdGg6IGludCwK
ICAgIGJhc2VzX2dlXzEwOiBpbnQsCiAgICBiYXNlc19nZV8yMDogaW50LAopIC0+IFR1cGxlW2lu
dCwgaW50LCBpbnQsIGludF06CiAgICAiIiJBZ2dyZWdhdGUgZGVwdGggb3ZlciBbZ3MsIGdlKS4g
SW50ZXJ2YWwgbGVuZ3RoIGlzIGFsd2F5cyAoZ2UgLSBncyk7IGdhcHMgaW4gdGFiaXggPSBkZXB0
aCAwLiIiIgogICAgc3BhbiA9IGdlIC0gZ3MKICAgIGlmIHNwYW4gPD0gMDoKICAgICAgICByZXR1
cm4gdG90YWxfYmFzZXMsIHN1bV9kZXB0aCwgYmFzZXNfZ2VfMTAsIGJhc2VzX2dlXzIwCiAgICBj
b3ZlcmVkID0gMAogICAgZm9yIGxpbmUgaW4gX3RhYml4X2xpbmVzKGJneiwgY2hyb20sIGdzLCBn
ZSk6CiAgICAgICAgcGFydHMgPSBsaW5lLnJzdHJpcCgpLnNwbGl0KCJcdCIpCiAgICAgICAgaWYg
bGVuKHBhcnRzKSA8IDQ6CiAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgdHJ5OgogICAgICAg
ICAgICBycyA9IGludChwYXJ0c1sxXSkKICAgICAgICAgICAgcmVfID0gaW50KHBhcnRzWzJdKQog
ICAgICAgICAgICBkZXB0aCA9IGludChmbG9hdChwYXJ0c1szXSkpCiAgICAgICAgZXhjZXB0IChW
YWx1ZUVycm9yLCBJbmRleEVycm9yKToKICAgICAgICAgICAgY29udGludWUKICAgICAgICBvc18g
PSBtYXgocnMsIGdzKQogICAgICAgIG9lID0gbWluKHJlXywgZ2UpCiAgICAgICAgaWYgb2UgPD0g
b3NfOgogICAgICAgICAgICBjb250aW51ZQogICAgICAgIGxuID0gb2UgLSBvc18KICAgICAgICBj
b3ZlcmVkICs9IGxuCiAgICAgICAgc3VtX2RlcHRoICs9IGRlcHRoICogbG4KICAgICAgICBpZiBk
ZXB0aCA+PSAxMDoKICAgICAgICAgICAgYmFzZXNfZ2VfMTAgKz0gbG4KICAgICAgICBpZiBkZXB0
aCA+PSAyMDoKICAgICAgICAgICAgYmFzZXNfZ2VfMjAgKz0gbG4KICAgIHVuY292ZXJlZCA9IHNw
YW4gLSBtaW4oY292ZXJlZCwgc3BhbikKICAgIHRvdGFsX2Jhc2VzICs9IHNwYW4KICAgICMgdW5j
b3ZlcmVkIGJhc2VzIGNvbnRyaWJ1dGUgMCBkZXB0aCBhbmQgZG8gbm90IGNvdW50IHRvd2FyZCA+
PTEwIC8gPj0yMAogICAgcmV0dXJuIHRvdGFsX2Jhc2VzLCBzdW1fZGVwdGgsIGJhc2VzX2dlXzEw
LCBiYXNlc19nZV8yMAoKCmRlZiB3cml0ZV90c3YoYmd6OiBzdHIsIGJlZDogc3RyLCBvdXRfcGF0
aDogc3RyKSAtPiBOb25lOgogICAgZ2VuZXMgPSBfcGFyc2VfYmVkKGJlZCkKICAgIHJvd3M6IExp
c3RbVHVwbGVbc3RyLCBmbG9hdCwgZmxvYXQsIGZsb2F0XV0gPSBbXQogICAgZm9yIGdlbmUgaW4g
c29ydGVkKGdlbmVzLmtleXMoKSk6CiAgICAgICAgbWVyZ2VkID0gX21lcmdlX2ludGVydmFscyhn
ZW5lc1tnZW5lXSkKICAgICAgICB0YiA9IHNkID0gYjEwID0gYjIwID0gMAogICAgICAgIGZvciBj
aHJvbSwgZ3MsIGdlIGluIG1lcmdlZDoKICAgICAgICAgICAgdGIsIHNkLCBiMTAsIGIyMCA9IF9h
Y2N1bXVsYXRlX2ludGVydmFsKGJneiwgY2hyb20sIGdzLCBnZSwgdGIsIHNkLCBiMTAsIGIyMCkK
ICAgICAgICBpZiB0YiA8PSAwOgogICAgICAgICAgICBtZWFuX2MgPSAwLjAKICAgICAgICAgICAg
cDEwID0gcDIwID0gMC4wCiAgICAgICAgZWxzZToKICAgICAgICAgICAgbWVhbl9jID0gc2QgLyB0
YgogICAgICAgICAgICBwMTAgPSAxMDAuMCAqIGIxMCAvIHRiCiAgICAgICAgICAgIHAyMCA9IDEw
MC4wICogYjIwIC8gdGIKICAgICAgICByb3dzLmFwcGVuZCgoZ2VuZSwgbWVhbl9jLCBwMTAsIHAy
MCkpCgogICAgd2l0aCBvcGVuKG91dF9wYXRoLCAidyIsIGVuY29kaW5nPSJ1dGYtOCIpIGFzIGY6
CiAgICAgICAgZi53cml0ZSgiZ2VuZVx0bWVhbl9jb3ZlcmFnZVx0cGN0X2Jhc2VzXzEweFx0cGN0
X2Jhc2VzXzIweFxuIikKICAgICAgICBmb3IgZ2VuZSwgbWVhbl9jLCBwMTAsIHAyMCBpbiByb3dz
OgogICAgICAgICAgICBmLndyaXRlKGYie2dlbmV9XHR7bWVhbl9jOi40Zn1cdHtwMTA6LjRmfVx0
e3AyMDouNGZ9XG4iKQoKCmRlZiBtYWluKCkgLT4gaW50OgogICAgYXAgPSBhcmdwYXJzZS5Bcmd1
bWVudFBhcnNlcihkZXNjcmlwdGlvbj0iR2VuZS1sZXZlbCBjb3ZlcmFnZSBUU1YgZnJvbSBtb3Nk
ZXB0aCBwZXItYmFzZSArIEJFRC4iKQogICAgYXAuYWRkX2FyZ3VtZW50KCJwZXJfYmFzZV9iZ3oi
LCBoZWxwPSJtb3NkZXB0aCAqLnBlci1iYXNlLmJlZC5neiAodGFiaXggaW5kZXhlZCkiKQogICAg
YXAuYWRkX2FyZ3VtZW50KCJnZW5lX2JlZCIsIGhlbHA9IkJFRCB3aXRoIEhHTkMgKG9yIGxhYmVs
KSBpbiBjb2x1bW4gNCIpCiAgICBhcC5hZGRfYXJndW1lbnQoIm91dF90c3YiLCBoZWxwPSJPdXRw
dXQgVFNWIHBhdGgiKQogICAgYXJncyA9IGFwLnBhcnNlX2FyZ3MoKQogICAgdHJ5OgogICAgICAg
IHdyaXRlX3RzdihhcmdzLnBlcl9iYXNlX2JneiwgYXJncy5nZW5lX2JlZCwgYXJncy5vdXRfdHN2
KQogICAgZXhjZXB0IE9TRXJyb3IgYXMgZToKICAgICAgICBwcmludChmIkVSUk9SOiB7ZX0iLCBm
aWxlPXN5cy5zdGRlcnIpCiAgICAgICAgcmV0dXJuIDEKICAgIHJldHVybiAwCgoKaWYgX19uYW1l
X18gPT0gIl9fbWFpbl9fIjoKICAgIHN5cy5leGl0KG1haW4oKSkK
MGC_B64
        python3 mosdepth_gene_coverage_tsv.py "\$PB" ${gene_bed} ${sample_id}_gene_coverage.tsv
    else
        echo -e "gene\\tmean_coverage\\tpct_bases_10x\\tpct_bases_20x" > ${sample_id}_gene_coverage.tsv
    fi
    """
}