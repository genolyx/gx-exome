process MANTA_SV {
    tag "$sample_id"
    label 'manta'
    publishDir "${params.outdir}/sv", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}_manta.vcf.gz"), emit: vcf

    script:
    """
    export TMPDIR=\$PWD

    if ! command -v configManta.py >/dev/null 2>&1; then
        export PATH=\$PWD:\$PATH
        export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
        export MAMBA_ROOT_PREFIX=\$PWD/micromamba
        mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
        wget -q --no-check-certificate https://curl.se/ca/cacert.pem
        export SSL_CERT_FILE=\$PWD/cacert.pem
        export MAMBA_SSL_VERIFY=false
        if [ ! -f "micromamba_bin" ]; then
            wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
                && chmod +x micromamba_bin
        fi
        [ -f micromamba_bin ] && ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge manta python=2.7 -y
        export PATH=\$PWD/env/bin:\$PATH
    fi

    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${ref_fasta} \\
        --runDir manta_run

    python2 manta_run/runWorkflow.py -m local -j ${task.cpus}

    mv manta_run/results/variants/diploidSV.vcf.gz ${sample_id}_manta.vcf.gz
    """
}
