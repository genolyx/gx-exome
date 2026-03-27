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
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & Tools
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    # Create env with samtools and mosdepth
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 mosdepth=0.3.3 -y
    export PATH=\$PWD/env/bin:\$PATH
    
    # --- Analysis ---

    # Fix: Symlink index for mosdepth/samtools (expects .bam.bai)
    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    # 1. Basic Bedcov (Keep existing)
    samtools bedcov ${bed} ${bam} > ${sample_id}_target_coverage.txt

    # 2. Mosdepth QC (Targeted)
    mosdepth --by ${bed} --thresholds 20,50,100 ${sample_id} ${bam}

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
    
    while read -r chrom start end name; do
        # Calculate depth stats for this region
        samtools depth -r "\$chrom:\$start-\$end" -a ${bam} | \\
        awk -v name="\$name" '{sum+=\$3; if(\$3>0) covered++; if(\$3>max) max=\$3; count++} END {if (count>0) print name "\t" sum/count "\t" max "\t" (covered/count)*100; else print name "\t0\t0\t0"}' >> ${sample_id}_intron_depth.txt
    done < ${dark_genes_bed}

    # Specific Loci Checks (Alpha / SMA)
    echo "\\n[Specific Loci Verification]" >> ${sample_id}_intron_depth.txt
    
    # Calculate Median Depth using awk/sort
    ALPHA_MED=\$(samtools depth -r "chr16:164000-183000" -a ${bam} | cut -f3 | sort -n | awk '{a[i++]=\$1;} END {if (i==0) print 0; else print a[int(i/2)];}')
    SMA_MED=\$(samtools depth -r "chr5:70920000-71000000" -a ${bam} | cut -f3 | sort -n | awk '{a[i++]=\$1;} END {if (i==0) print 0; else print a[int(i/2)];}')
    
    echo "Alpha-Cluster Median Depth: \$ALPHA_MED" >> ${sample_id}_intron_depth.txt
    echo "SMA Locus Median Depth: \$SMA_MED" >> ${sample_id}_intron_depth.txt
    
    # Check Threshold (15x)
    FAILED=0
    if [ "\$ALPHA_MED" -lt 15 ]; then FAILED=1; fi
    if [ "\$SMA_MED" -lt 15 ]; then FAILED=1; fi
    
    if [ "\$FAILED" -eq 1 ]; then
        echo "WARNING: Median depth below 15x in Alpha-Cluster (\$ALPHA_MED) or SMA (\$SMA_MED)." >> ${sample_id}_intron_depth.txt
    fi
    """
}
