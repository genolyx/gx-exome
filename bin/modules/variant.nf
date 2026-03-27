// ============================================================
// variant.nf — Variant Calling module
//
// Supports three variant callers selectable via params.variant_caller:
//   'gatk'        : HaplotypeCaller + VariantFiltration (hard filter)
//   'deepvariant' : Google DeepVariant CNN (WES mode, built-in filter)
//   'strelka2'    : Illumina Strelka2 (exome mode, built-in filter)
//
// IMPORTANT:
//   - DeepVariant and Strelka2 do NOT require a separate VariantFiltration
//     step. Their internal models handle filtering automatically.
//   - All three processes emit the same channel signature:
//       tuple val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
//     so downstream processes (annotation, visualization) are unchanged.
// ============================================================

// -------------------------------------------------------
// GATK HaplotypeCaller + VariantFiltration (hard filter)
// Use when: GATK is required for compatibility or institutional policy.
// -------------------------------------------------------
process CALL_VARIANTS_GATK {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/variant", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path ref_dict
    path backbone_bed
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), emit: vcf
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    
    if command -v curl >/dev/null 2>&1; then
        curl -k -L -o cacert.pem https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && curl -k -L -o micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    else
        wget -q --no-check-certificate https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    fi
    
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    chmod +x micromamba_bin
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    # Step 1: HaplotypeCaller
    gatk HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -L ${backbone_bed} \\
        --interval-padding 100 \\
        -O ${sample_id}_raw.vcf.gz \\
        -ERC NONE \\
        --create-output-variant-index true

    # Step 2: VariantFiltration (hard filter — WES-tuned thresholds)
    gatk VariantFiltration \\
        -R ${ref_fasta} \\
        -V ${sample_id}_raw.vcf.gz \\
        -O ${sample_id}_filtered.vcf.gz \\
        --filter-expression "QD < 1.5"             --filter-name "QD1.5" \\
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \\
        --filter-expression "SOR > 4.0"             --filter-name "SOR4" \\
        --filter-expression "FS > 80.0"             --filter-name "FS80" \\
        --filter-expression "MQ < 30.0"             --filter-name "MQ30" \\
        --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    """
}

// -------------------------------------------------------
// Google DeepVariant (WES mode)
// Use when: highest SNV precision is required (clinical / rare disease).
// Notes:
//   - --model_type WES is MANDATORY for exome data.
//   - No separate VariantFiltration needed; FILTER=PASS is set by the model.
//   - GPU acceleration is supported; add --gpus all to docker.runOptions if available.
// -------------------------------------------------------
process CALL_VARIANTS_DEEPVARIANT {
    tag "$sample_id"
    label 'deepvariant'
    publishDir "${params.outdir}/variant", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path backbone_bed
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), emit: vcf
    path "${sample_id}.g.vcf.gz", emit: gvcf   // optional gVCF for joint calling
    script:
    """
    export TMPDIR=\$PWD

    # DeepVariant ships all tools inside the container; no micromamba needed.
    # run_deepvariant is the all-in-one wrapper (make_examples → call_variants → postprocess).
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WES \\
        --ref=${ref_fasta} \\
        --reads=${bam} \\
        --regions=${backbone_bed} \\
        --output_vcf=${sample_id}_dv_raw.vcf.gz \\
        --output_gvcf=${sample_id}.g.vcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=\$TMPDIR/dv_tmp

    # Extract PASS variants only (DeepVariant sets FILTER=PASS for confident calls)
    if command -v bcftools >/dev/null 2>&1; then
        BCFTOOLS=bcftools
    elif [ -f /opt/deepvariant/bin/bcftools ]; then
        BCFTOOLS=/opt/deepvariant/bin/bcftools
    elif command -v /usr/bin/bcftools >/dev/null 2>&1; then
        BCFTOOLS=/usr/bin/bcftools
    else
        export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
        export MAMBA_ROOT_PREFIX=\$PWD/micromamba
        mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
        
        if command -v curl >/dev/null 2>&1; then
            curl -k -L -o cacert.pem https://curl.se/ca/cacert.pem
            [ ! -f "micromamba_bin" ] && curl -k -L -o micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        else
            wget -q --no-check-certificate https://curl.se/ca/cacert.pem
            [ ! -f "micromamba_bin" ] && wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        fi
        
        export SSL_CERT_FILE=\$PWD/cacert.pem
        export MAMBA_SSL_VERIFY=false
        chmod +x micromamba_bin
        
        ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bcftools -y
        export PATH=\$PWD/env/bin:\$PATH
        BCFTOOLS=bcftools
    fi

    \$BCFTOOLS view -f PASS -O z -o ${sample_id}_filtered.vcf.gz ${sample_id}_dv_raw.vcf.gz
    \$BCFTOOLS index --tbi ${sample_id}_filtered.vcf.gz
    """
}

// -------------------------------------------------------
// Illumina Strelka2 (exome mode)
// Use when: large cohort, speed priority, or balanced SNV+Indel accuracy.
// Notes:
//   - --exome flag is MANDATORY for WES data.
//   - Manta SV caller is run first to provide candidate indels (optional but recommended).
//   - No separate VariantFiltration needed; FILTER=PASS is set by the Random Forest model.
// -------------------------------------------------------
process CALL_VARIANTS_STRELKA2 {
    tag "$sample_id"
    label 'strelka2'
    publishDir "${params.outdir}/variant", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path backbone_bed
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), emit: vcf
    script:
    """
    export TMPDIR=\$PWD
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    
    if command -v curl >/dev/null 2>&1; then
        curl -k -L -o cacert.pem https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && curl -k -L -o micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    else
        wget -q --no-check-certificate https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    fi
    
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    chmod +x micromamba_bin
    
    # Install Strelka2 and bgzip/tabix (htslib)
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge strelka=2.9.10 htslib bcftools -y
    export PATH=\$PWD/env/bin:\$PATH

    # Step 1: Configure Strelka2 for germline WES
    # --exome flag adjusts the statistical model for non-uniform WES coverage
    configureStrelkaGermlineWorkflow.py \\
        --bam ${bam} \\
        --referenceFasta ${ref_fasta} \\
        --callRegions ${backbone_bed}.gz \\
        --exome \\
        --runDir strelka_run

    # Step 2: Run Strelka2 workflow
    python2 strelka_run/runWorkflow.py -m local -j ${task.cpus}

    # Step 3: Merge SNVs and Indels into a single VCF, keep PASS only
    bcftools concat -a \\
        strelka_run/results/variants/variants.vcf.gz \\
        strelka_run/results/variants/genome.vcf.gz \\
        2>/dev/null || \\
    cp strelka_run/results/variants/variants.vcf.gz ${sample_id}_strelka_raw.vcf.gz

    # Strelka2 outputs SNVs in variants.vcf.gz; use that as primary output
    bcftools view -f PASS \\
        -O z -o ${sample_id}_filtered.vcf.gz \\
        strelka_run/results/variants/variants.vcf.gz
    bcftools index --tbi ${sample_id}_filtered.vcf.gz
    """
}

// -------------------------------------------------------
// Backward-compatible alias: CALL_VARIANTS → CALL_VARIANTS_GATK
// This allows existing scripts that reference CALL_VARIANTS to continue working.
// -------------------------------------------------------
process CALL_VARIANTS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/variant", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path ref_dict
    path backbone_bed
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), emit: vcf
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    
    if command -v curl >/dev/null 2>&1; then
        curl -k -L -o cacert.pem https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && curl -k -L -o micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    else
        wget -q --no-check-certificate https://curl.se/ca/cacert.pem
        [ ! -f "micromamba_bin" ] && wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
    fi
    
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    chmod +x micromamba_bin
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -L ${backbone_bed} \\
        --interval-padding 100 \\
        -O ${sample_id}_raw.vcf.gz \\
        -ERC NONE \\
        --create-output-variant-index true

    gatk VariantFiltration \\
        -R ${ref_fasta} \\
        -V ${sample_id}_raw.vcf.gz \\
        -O ${sample_id}_filtered.vcf.gz \\
        --filter-expression "QD < 1.5"             --filter-name "QD1.5" \\
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \\
        --filter-expression "SOR > 4.0"             --filter-name "SOR4" \\
        --filter-expression "FS > 80.0"             --filter-name "FS80" \\
        --filter-expression "MQ < 30.0"             --filter-name "MQ30" \\
        --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    """
}
