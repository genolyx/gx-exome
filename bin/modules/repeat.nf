process EXPANSION_HUNTER {
    tag "$sample_id"
    label 'expansionhunter'
    publishDir "${params.outdir}/repeat", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path variant_catalog // JSON specifying FMR1

    output:
    tuple val(sample_id), path("${sample_id}_eh.json"), path("${sample_id}_eh.vcf"), emit: results
    tuple val(sample_id), path("*.svg"), emit: images, optional: true
    tuple val(sample_id), path("${sample_id}_eh_realigned.sorted.bam"), path("${sample_id}_eh_realigned.sorted.bam.bai"), emit: eh_realigned

    script:
    """
    # Micromamba + NXF_DOCKER_TASK_USER: HOME 없으면 /.conda 실패, || true 는 samtools 미설치를 숨김
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
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
            && chmod +x micromamba_bin || true
    fi
    
    if [ -f micromamba_bin ] && [ ! -x ./env/bin/samtools ]; then
        ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge expansionhunter=4.0.2 samtools -y
    fi
    [ -x ./env/bin/samtools ] && [ -x ./env/bin/ExpansionHunter ] || { echo "ERROR: ExpansionHunter/samtools env missing"; exit 1; }

    ST=\$PWD/env/bin/samtools
    EH=\$PWD/env/bin/ExpansionHunter

    # Install REViewer manually (Standalone Binary)
    if [ ! -f "REViewer" ]; then
        # Check for cached copy first (optional optimization)
        wget -qO REViewer.gz https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz
        gzip -d REViewer.gz
        chmod +x REViewer
    fi

    # ExpansionHunter can be picky about index naming
    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    # --- Auto-Sex Detection ---
    # Heuristic: Ratio of aligned reads on chrY vs chrX
    \$ST idxstats ${bam} > idxstats.txt
    
    # Extract mapped read counts (Column 3) for chrX and chrY
    read_x=\$(grep -w "chrX" idxstats.txt | cut -f3)
    read_y=\$(grep -w "chrY" idxstats.txt | cut -f3)

    # Handle cases where chromosome might be missing or zero coverage
    if [ -z "\$read_x" ]; then read_x=1; fi
    if [ -z "\$read_y" ]; then read_y=0; fi

    # Calculate Ratio (Y / X)
    # Threshold: 0.05. Females (XX) have typically <0.01 Y reads (noise). Males (XY) have distinct Y mapping.
    sex="female"
    is_male=\$(awk -v y=\$read_y -v x=\$read_x 'BEGIN {print (x>0 && (y/x) > 0.05) ? 1 : 0}')

    if [ "\$is_male" -eq 1 ]; then
        sex="male"
    fi
    
    echo "Inferred Sex: \$sex (chrX: \$read_x, chrY: \$read_y, Y/X Ratio)"

    # 1. Run ExpansionHunter
    # v4.0.2 automatically generates realigned BAMs (no flag needed)
    \$EH \\
        --reads ${bam} \\
        --reference ${ref_fasta} \\
        --variant-catalog ${variant_catalog} \\
        --sex \$sex \\
        --output-prefix ${sample_id}_eh

    # EH v4 realigned BAM is unsorted; sort+index for REViewer and IGV visual evidence
    \$ST sort ${sample_id}_eh_realigned.bam -o ${sample_id}_eh_realigned.sorted.bam
    \$ST index ${sample_id}_eh_realigned.sorted.bam

    # 2. Run REViewer (Visualization) for FMR1
    # Check if FMR1 was called in the VCF to avoid errors if missing
    if grep -q "FMR1" ${sample_id}_eh.vcf; then
        echo "Generating REViewer visualization for FMR1..."

        ./REViewer \\
            --reads ${sample_id}_eh_realigned.sorted.bam \\
            --vcf ${sample_id}_eh.vcf \\
            --reference ${ref_fasta} \\
            --catalog ${variant_catalog} \\
            --locus FMR1 \\
            --output-prefix ${sample_id}_FMR1
    else
        echo "FMR1 not found in VCF, skipping visualization."
    fi

    # 3. CFTR intron 9 poly-TG / poly-T tract (carrier screening: 5T/7T/9T x TG11/12/13).
    # EH preserves per-allele order across both VariantIds in the CFTR_IVS9 locus,
    # so REPCN columns can be combined into haplotypes downstream in summary.nf.
    if grep -qE "CFTR_TG|CFTR_polyT" ${sample_id}_eh.vcf; then
        echo "Generating REViewer visualization for CFTR_IVS9..."
        ./REViewer \\
            --reads ${sample_id}_eh_realigned.sorted.bam \\
            --vcf ${sample_id}_eh.vcf \\
            --reference ${ref_fasta} \\
            --catalog ${variant_catalog} \\
            --locus CFTR_IVS9 \\
            --output-prefix ${sample_id}_CFTR_IVS9 || echo "WARN: REViewer failed on CFTR_IVS9 (continuing)."
    else
        echo "CFTR_IVS9 not found in VCF, skipping visualization."
    fi

    """
}
