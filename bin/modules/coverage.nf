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
cnQgYXJncGFyc2UKaW1wb3J0IHJlCmltcG9ydCBzdWJwcm9jZXNzCmltcG9ydCBzeXMKZnJvbSBj
b2xsZWN0aW9ucyBpbXBvcnQgZGVmYXVsdGRpY3QKZnJvbSB0eXBpbmcgaW1wb3J0IERpY3QsIExp
c3QsIFR1cGxlCgoKSW50ZXJ2YWwgPSBUdXBsZVtzdHIsIGludCwgaW50XSAgIyBjaHJvbSwgc3Rh
cnQsIGVuZCAoMC1iYXNlZCBoYWxmLW9wZW4pCgpfR0VORV9TWU0gPSByZS5jb21waWxlKHIiXltB
LVphLXpdW0EtWmEtejAtOS1dezAsMjR9JCIpCl9TS0lQX1RPS19QUkVGSVhFUyA9ICgKICAgICJF
TlNUIiwKICAgICJOTV8iLAogICAgIk5SXyIsCiAgICAiWE1fIiwKICAgICJYUl8iLAogICAgIkVO
U0ciLAogICAgIkNDRFMiLAogICAgIkNMSU5JRCIsCiAgICAiTE9DIiwKKQoKCmRlZiBfZ2VuZV9r
ZXlzX2Zyb21fYmVkX2NvbDQobmFtZTogc3RyKSAtPiBMaXN0W3N0cl06CiAgICAiIiIKICAgIEhH
TkMtbGlrZSBzeW1ib2xzIGluIGNvbHVtbiA0LiBUd2lzdC92ZW5kb3IgQkVEcyB1c2UgYGBQQUg7
Tk1fLi4uO0VOU1QuLi5gYDsgcGxhaW4gQkVEcyB1c2UgYGBQQUhgYC4KICAgICIiIgogICAgcmF3
ID0gKG5hbWUgb3IgIiIpLnN0cmlwKCkKICAgIGlmIG5vdCByYXc6CiAgICAgICAgcmV0dXJuIFtd
CiAgICBrZXlzOiBMaXN0W3N0cl0gPSBbXQogICAgZm9yIGNodW5rIGluIHJhdy5yZXBsYWNlKCIs
IiwgIjsiKS5yZXBsYWNlKCJ8IiwgIjsiKS5zcGxpdCgiOyIpOgogICAgICAgIHQgPSBjaHVuay5z
dHJpcCgpCiAgICAgICAgaWYgbm90IHQ6CiAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgdWwg
PSB0LnVwcGVyKCkKICAgICAgICBpZiBhbnkodWwuc3RhcnRzd2l0aChwKSBmb3IgcCBpbiBfU0tJ
UF9UT0tfUFJFRklYRVMpOgogICAgICAgICAgICBjb250aW51ZQogICAgICAgIGlmIF9HRU5FX1NZ
TS5tYXRjaCh0KToKICAgICAgICAgICAga2V5cy5hcHBlbmQodWwpCiAgICBpZiBrZXlzOgogICAg
ICAgIHJldHVybiBrZXlzCiAgICBpZiBfR0VORV9TWU0ubWF0Y2gocmF3KToKICAgICAgICByZXR1
cm4gW3Jhdy51cHBlcigpXQogICAgcmV0dXJuIFtdCgoKZGVmIF9wYXJzZV9iZWQocGF0aDogc3Ry
KSAtPiBEaWN0W3N0ciwgTGlzdFtJbnRlcnZhbF1dOgogICAgIiIiR3JvdXAgaW50ZXJ2YWxzIGJ5
IEhHTkMgc3ltYm9sIGRlcml2ZWQgZnJvbSBjb2x1bW4gNC4iIiIKICAgIGJ5X2dlbmU6IERpY3Rb
c3RyLCBMaXN0W0ludGVydmFsXV0gPSBkZWZhdWx0ZGljdChsaXN0KQogICAgd2l0aCBvcGVuKHBh
dGgsIGVuY29kaW5nPSJ1dGYtOCIsIGVycm9ycz0icmVwbGFjZSIpIGFzIGY6CiAgICAgICAgZm9y
IGxpbmUgaW4gZjoKICAgICAgICAgICAgbGluZSA9IGxpbmUuc3RyaXAoKQogICAgICAgICAgICBp
ZiBub3QgbGluZSBvciBsaW5lLnN0YXJ0c3dpdGgoIiMiKSBvciBsaW5lLnN0YXJ0c3dpdGgoInRy
YWNrICIpOgogICAgICAgICAgICAgICAgY29udGludWUKICAgICAgICAgICAgcGFydHMgPSBsaW5l
LnNwbGl0KCJcdCIpCiAgICAgICAgICAgIGlmIGxlbihwYXJ0cykgPCA0OgogICAgICAgICAgICAg
ICAgY29udGludWUKICAgICAgICAgICAgY2hyb20sIHMsIGUsIG5hbWUgPSBwYXJ0c1swXSwgcGFy
dHNbMV0sIHBhcnRzWzJdLCBwYXJ0c1szXS5zdHJpcCgpCiAgICAgICAgICAgIGlmIG5vdCBuYW1l
OgogICAgICAgICAgICAgICAgY29udGludWUKICAgICAgICAgICAgdHJ5OgogICAgICAgICAgICAg
ICAgc3RhcnRfaSA9IGludChzKQogICAgICAgICAgICAgICAgZW5kX2kgPSBpbnQoZSkKICAgICAg
ICAgICAgZXhjZXB0IFZhbHVlRXJyb3I6CiAgICAgICAgICAgICAgICBjb250aW51ZQogICAgICAg
ICAgICBpZiBlbmRfaSA8PSBzdGFydF9pOgogICAgICAgICAgICAgICAgY29udGludWUKICAgICAg
ICAgICAgZm9yIGdrZXkgaW4gX2dlbmVfa2V5c19mcm9tX2JlZF9jb2w0KG5hbWUpOgogICAgICAg
ICAgICAgICAgYnlfZ2VuZVtna2V5XS5hcHBlbmQoKGNocm9tLCBzdGFydF9pLCBlbmRfaSkpCiAg
ICByZXR1cm4gYnlfZ2VuZQoKCmRlZiBfbWVyZ2VfaW50ZXJ2YWxzKGludGVydmFsczogTGlzdFtJ
bnRlcnZhbF0pIC0+IExpc3RbSW50ZXJ2YWxdOgogICAgYnlfY2hyb206IERpY3Rbc3RyLCBMaXN0
W1R1cGxlW2ludCwgaW50XV1dID0gZGVmYXVsdGRpY3QobGlzdCkKICAgIGZvciBjaHJvbSwgcywg
ZSBpbiBpbnRlcnZhbHM6CiAgICAgICAgYnlfY2hyb21bY2hyb21dLmFwcGVuZCgocywgZSkpCiAg
ICBtZXJnZWQ6IExpc3RbSW50ZXJ2YWxdID0gW10KICAgIGZvciBjaHJvbSwgaXZzIGluIHNvcnRl
ZChieV9jaHJvbS5pdGVtcygpLCBrZXk9bGFtYmRhIHg6IHhbMF0pOgogICAgICAgIGl2cy5zb3J0
KGtleT1sYW1iZGEgeDogKHhbMF0sIHhbMV0pKQogICAgICAgIGN1cl9zLCBjdXJfZSA9IGl2c1sw
XQogICAgICAgIGZvciBzLCBlIGluIGl2c1sxOl06CiAgICAgICAgICAgIGlmIHMgPD0gY3VyX2U6
CiAgICAgICAgICAgICAgICBjdXJfZSA9IG1heChjdXJfZSwgZSkKICAgICAgICAgICAgZWxzZToK
ICAgICAgICAgICAgICAgIG1lcmdlZC5hcHBlbmQoKGNocm9tLCBjdXJfcywgY3VyX2UpKQogICAg
ICAgICAgICAgICAgY3VyX3MsIGN1cl9lID0gcywgZQogICAgICAgIG1lcmdlZC5hcHBlbmQoKGNo
cm9tLCBjdXJfcywgY3VyX2UpKQogICAgcmV0dXJuIG1lcmdlZAoKCmRlZiBfYWx0X2Nocm9tKGM6
IHN0cikgLT4gc3RyOgogICAgaWYgYy5zdGFydHN3aXRoKCJjaHIiKToKICAgICAgICByZXR1cm4g
Y1szOl0gaWYgbGVuKGMpID4gMyBlbHNlIGMKICAgIHJldHVybiAiY2hyIiArIGMKCgpkZWYgX3Rh
Yml4X2xpbmVzKGJnejogc3RyLCBjaHJvbTogc3RyLCBzdGFydDogaW50LCBlbmQ6IGludCkgLT4g
TGlzdFtzdHJdOgogICAgIiIiUXVlcnkgcGVyLWJhc2UgQkVEOyB1c2UgMC1iYXNlZCBoYWxmLW9w
ZW4gaW50ZXJ2YWxzIChtYXRjaGVzIEJFRCBzdG9yYWdlKS4iIiIKICAgIGZvciBjaCBpbiAoY2hy
b20sIF9hbHRfY2hyb20oY2hyb20pKToKICAgICAgICByID0gZiJ7Y2h9OntzdGFydH0te2VuZH0i
CiAgICAgICAgdHJ5OgogICAgICAgICAgICBvdXQgPSBzdWJwcm9jZXNzLmNoZWNrX291dHB1dCgK
ICAgICAgICAgICAgICAgIFsidGFiaXgiLCAiLS16ZXJvLWJhc2VkIiwgYmd6LCByXSwKICAgICAg
ICAgICAgICAgIHN0ZGVycj1zdWJwcm9jZXNzLkRFVk5VTEwsCiAgICAgICAgICAgICAgICB0ZXh0
PVRydWUsCiAgICAgICAgICAgICkKICAgICAgICBleGNlcHQgc3VicHJvY2Vzcy5DYWxsZWRQcm9j
ZXNzRXJyb3I6CiAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgaWYgb3V0LnN0cmlwKCk6CiAg
ICAgICAgICAgIHJldHVybiBvdXQuc3BsaXRsaW5lcygpCiAgICByZXR1cm4gW10KCgpkZWYgX2Fj
Y3VtdWxhdGVfaW50ZXJ2YWwoCiAgICBiZ3o6IHN0ciwKICAgIGNocm9tOiBzdHIsCiAgICBnczog
aW50LAogICAgZ2U6IGludCwKICAgIHRvdGFsX2Jhc2VzOiBpbnQsCiAgICBzdW1fZGVwdGg6IGlu
dCwKICAgIGJhc2VzX2dlXzEwOiBpbnQsCiAgICBiYXNlc19nZV8yMDogaW50LAopIC0+IFR1cGxl
W2ludCwgaW50LCBpbnQsIGludF06CiAgICAiIiJBZ2dyZWdhdGUgZGVwdGggb3ZlciBbZ3MsIGdl
KS4gSW50ZXJ2YWwgbGVuZ3RoIGlzIGFsd2F5cyAoZ2UgLSBncyk7IGdhcHMgaW4gdGFiaXggPSBk
ZXB0aCAwLiIiIgogICAgc3BhbiA9IGdlIC0gZ3MKICAgIGlmIHNwYW4gPD0gMDoKICAgICAgICBy
ZXR1cm4gdG90YWxfYmFzZXMsIHN1bV9kZXB0aCwgYmFzZXNfZ2VfMTAsIGJhc2VzX2dlXzIwCiAg
ICBjb3ZlcmVkID0gMAogICAgZm9yIGxpbmUgaW4gX3RhYml4X2xpbmVzKGJneiwgY2hyb20sIGdz
LCBnZSk6CiAgICAgICAgcGFydHMgPSBsaW5lLnJzdHJpcCgpLnNwbGl0KCJcdCIpCiAgICAgICAg
aWYgbGVuKHBhcnRzKSA8IDQ6CiAgICAgICAgICAgIGNvbnRpbnVlCiAgICAgICAgdHJ5OgogICAg
ICAgICAgICBycyA9IGludChwYXJ0c1sxXSkKICAgICAgICAgICAgcmVfID0gaW50KHBhcnRzWzJd
KQogICAgICAgICAgICBkZXB0aCA9IGludChmbG9hdChwYXJ0c1szXSkpCiAgICAgICAgZXhjZXB0
IChWYWx1ZUVycm9yLCBJbmRleEVycm9yKToKICAgICAgICAgICAgY29udGludWUKICAgICAgICBv
c18gPSBtYXgocnMsIGdzKQogICAgICAgIG9lID0gbWluKHJlXywgZ2UpCiAgICAgICAgaWYgb2Ug
PD0gb3NfOgogICAgICAgICAgICBjb250aW51ZQogICAgICAgIGxuID0gb2UgLSBvc18KICAgICAg
ICBjb3ZlcmVkICs9IGxuCiAgICAgICAgc3VtX2RlcHRoICs9IGRlcHRoICogbG4KICAgICAgICBp
ZiBkZXB0aCA+PSAxMDoKICAgICAgICAgICAgYmFzZXNfZ2VfMTAgKz0gbG4KICAgICAgICBpZiBk
ZXB0aCA+PSAyMDoKICAgICAgICAgICAgYmFzZXNfZ2VfMjAgKz0gbG4KICAgIHVuY292ZXJlZCA9
IHNwYW4gLSBtaW4oY292ZXJlZCwgc3BhbikKICAgIHRvdGFsX2Jhc2VzICs9IHNwYW4KICAgICMg
dW5jb3ZlcmVkIGJhc2VzIGNvbnRyaWJ1dGUgMCBkZXB0aCBhbmQgZG8gbm90IGNvdW50IHRvd2Fy
ZCA+PTEwIC8gPj0yMAogICAgcmV0dXJuIHRvdGFsX2Jhc2VzLCBzdW1fZGVwdGgsIGJhc2VzX2dl
XzEwLCBiYXNlc19nZV8yMAoKCmRlZiB3cml0ZV90c3YoYmd6OiBzdHIsIGJlZDogc3RyLCBvdXRf
cGF0aDogc3RyKSAtPiBOb25lOgogICAgZ2VuZXMgPSBfcGFyc2VfYmVkKGJlZCkKICAgIHJvd3M6
IExpc3RbVHVwbGVbc3RyLCBmbG9hdCwgZmxvYXQsIGZsb2F0XV0gPSBbXQogICAgZm9yIGdlbmUg
aW4gc29ydGVkKGdlbmVzLmtleXMoKSk6CiAgICAgICAgbWVyZ2VkID0gX21lcmdlX2ludGVydmFs
cyhnZW5lc1tnZW5lXSkKICAgICAgICB0YiA9IHNkID0gYjEwID0gYjIwID0gMAogICAgICAgIGZv
ciBjaHJvbSwgZ3MsIGdlIGluIG1lcmdlZDoKICAgICAgICAgICAgdGIsIHNkLCBiMTAsIGIyMCA9
IF9hY2N1bXVsYXRlX2ludGVydmFsKGJneiwgY2hyb20sIGdzLCBnZSwgdGIsIHNkLCBiMTAsIGIy
MCkKICAgICAgICBpZiB0YiA8PSAwOgogICAgICAgICAgICBtZWFuX2MgPSAwLjAKICAgICAgICAg
ICAgcDEwID0gcDIwID0gMC4wCiAgICAgICAgZWxzZToKICAgICAgICAgICAgbWVhbl9jID0gc2Qg
LyB0YgogICAgICAgICAgICBwMTAgPSAxMDAuMCAqIGIxMCAvIHRiCiAgICAgICAgICAgIHAyMCA9
IDEwMC4wICogYjIwIC8gdGIKICAgICAgICByb3dzLmFwcGVuZCgoZ2VuZSwgbWVhbl9jLCBwMTAs
IHAyMCkpCgogICAgd2l0aCBvcGVuKG91dF9wYXRoLCAidyIsIGVuY29kaW5nPSJ1dGYtOCIpIGFz
IGY6CiAgICAgICAgZi53cml0ZSgiZ2VuZVx0bWVhbl9jb3ZlcmFnZVx0cGN0X2Jhc2VzXzEweFx0
cGN0X2Jhc2VzXzIweFxuIikKICAgICAgICBmb3IgZ2VuZSwgbWVhbl9jLCBwMTAsIHAyMCBpbiBy
b3dzOgogICAgICAgICAgICBmLndyaXRlKGYie2dlbmV9XHR7bWVhbl9jOi40Zn1cdHtwMTA6LjRm
fVx0e3AyMDouNGZ9XG4iKQoKCmRlZiBtYWluKCkgLT4gaW50OgogICAgYXAgPSBhcmdwYXJzZS5B
cmd1bWVudFBhcnNlcihkZXNjcmlwdGlvbj0iR2VuZS1sZXZlbCBjb3ZlcmFnZSBUU1YgZnJvbSBt
b3NkZXB0aCBwZXItYmFzZSArIEJFRC4iKQogICAgYXAuYWRkX2FyZ3VtZW50KCJwZXJfYmFzZV9i
Z3oiLCBoZWxwPSJtb3NkZXB0aCAqLnBlci1iYXNlLmJlZC5neiAodGFiaXggaW5kZXhlZCkiKQog
ICAgYXAuYWRkX2FyZ3VtZW50KCJnZW5lX2JlZCIsIGhlbHA9IkJFRCB3aXRoIEhHTkMgKG9yIGxh
YmVsKSBpbiBjb2x1bW4gNCIpCiAgICBhcC5hZGRfYXJndW1lbnQoIm91dF90c3YiLCBoZWxwPSJP
dXRwdXQgVFNWIHBhdGgiKQogICAgYXJncyA9IGFwLnBhcnNlX2FyZ3MoKQogICAgdHJ5OgogICAg
ICAgIHdyaXRlX3RzdihhcmdzLnBlcl9iYXNlX2JneiwgYXJncy5nZW5lX2JlZCwgYXJncy5vdXRf
dHN2KQogICAgZXhjZXB0IE9TRXJyb3IgYXMgZToKICAgICAgICBwcmludChmIkVSUk9SOiB7ZX0i
LCBmaWxlPXN5cy5zdGRlcnIpCiAgICAgICAgcmV0dXJuIDEKICAgIHJldHVybiAwCgoKaWYgX19u
YW1lX18gPT0gIl9fbWFpbl9fIjoKICAgIHN5cy5leGl0KG1haW4oKSkK
MGC_B64
        python3 mosdepth_gene_coverage_tsv.py "\$PB" ${gene_bed} ${sample_id}_gene_coverage.tsv
    else
        echo -e "gene\\tmean_coverage\\tpct_bases_10x\\tpct_bases_20x" > ${sample_id}_gene_coverage.tsv
    fi
    """
}