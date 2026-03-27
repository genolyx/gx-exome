// ============================================================
// align.nf — Alignment module
//
// Supports two aligners selectable via params.aligner:
//   'bwa-mem'  : classic BWA-MEM (INDEX_BWA / ALIGN_AND_SORT)
//   'bwa-mem2' : BWA-MEM2 (INDEX_BWA_MEM2 / ALIGN_AND_SORT_BWA_MEM2)
//
// Common downstream processes (MARK_DUPLICATES, SAMTOOLS_BAM_STATS)
// are shared regardless of aligner choice.
// ============================================================

// -------------------------------------------------------
// BWA-MEM (classic) — Index
// -------------------------------------------------------
process INDEX_BWA {
    tag "bwa_index"
    label 'bwa'
    publishDir "${params.outdir}/bwa_index", mode: 'copy'
    input:
    path ref_fasta
    output:
    tuple path(ref_fasta), path("*.{amb,ann,bwt,pac,sa}"), emit: indices
    script:
    """
    export TMPDIR=\$PWD
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa=0.7.17 -y
    export PATH=\$PWD/env/bin:\$PATH
    bwa index ${ref_fasta}
    """
}

// -------------------------------------------------------
// BWA-MEM2 — Index
// Generates: .0123, .amb, .ann, .bwt.2bit.64, .pac
// NOTE: This is a one-time operation (~30 min for GRCh38).
//       Pre-building on the server is strongly recommended.
// -------------------------------------------------------
process INDEX_BWA_MEM2 {
    tag "bwa_mem2_index"
    label 'bwa_mem2'
    publishDir "${params.outdir}/bwa_mem2_index", mode: 'copy'
    input:
    path ref_fasta
    output:
    tuple path(ref_fasta), path("*.{0123,amb,ann,bwt.2bit.64,pac}"), emit: indices
    script:
    """
    export TMPDIR=\$PWD
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa-mem2=2.2.1 samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH
    bwa-mem2 index ${ref_fasta}
    """
}

// -------------------------------------------------------
// BWA-MEM (classic) — Align & Sort
// -------------------------------------------------------
process ALIGN_AND_SORT {
    tag "$sample_id"
    label 'bwa'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path ref_fai
    path bwa_indices
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa=0.7.17 samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH

    # Pipe bwa mem directly into samtools sort to avoid large intermediate SAM
    bwa mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \\
        ${ref_fasta} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -m 768M -o ${sample_id}.bam -

    samtools index ${sample_id}.bam
    """
}

// -------------------------------------------------------
// BWA-MEM2 — Align & Sort
// Key differences from BWA-MEM:
//   1. Uses 'bwa-mem2 mem' command
//   2. Index files have different extensions (.0123, .bwt.2bit.64)
//   3. Requires ~64 GB RAM to load the index
// -------------------------------------------------------
process ALIGN_AND_SORT_BWA_MEM2 {
    tag "$sample_id"
    label 'bwa_mem2'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path ref_fai
    path bwa_mem2_indices
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa-mem2=2.2.1 samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH

    # bwa-mem2 uses the same CLI as bwa mem, but the binary is 'bwa-mem2'
    # Pipe directly into samtools sort to avoid large intermediate SAM
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \\
        ${ref_fasta} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -m 768M -o ${sample_id}.bam -

    samtools index ${sample_id}.bam
    """
}

// -------------------------------------------------------
// MarkDuplicates — shared for both aligners
// -------------------------------------------------------
process MARK_DUPLICATES {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: bam
    path "${sample_id}.duplicate_metrics.txt", emit: metrics
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$XDG_CACHE_HOME
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.md.bam \\
        -M ${sample_id}.duplicate_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT
    """
}

// -------------------------------------------------------
// BAM QC Stats — shared for both aligners
// -------------------------------------------------------
process SAMTOOLS_BAM_STATS {
    tag "$sample_id"
    label 'bwa'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    path "${sample_id}.stats.txt", emit: stats
    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH

    samtools stats -@ ${task.cpus} ${bam} > ${sample_id}.stats.txt
    """
}
