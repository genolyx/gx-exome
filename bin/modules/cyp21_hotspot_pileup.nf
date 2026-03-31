// BAM pileup at the seven NM_000500.9 sites in cyp21a2_hotspots.tsv (cross-check vs SNV VCF).
// Imports cyp21_paralog_pileup.py from the same task work dir.

process CYP21_HOTSPOT_PILEUP {
    tag "$sample_id"
    publishDir "${params.outdir}/fallback", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(hotspots_tsv), path(hotspot_py), path(paralog_pileup_py)

    output:
    path "${sample_id}_cyp21_hotspot_pileup.tsv", emit: tsv

    script:
    """
    export TMPDIR=\$PWD
    export HOME=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX

    wget -q --no-check-certificate https://curl.se/ca/cacert.pem || true
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \\
            && chmod +x micromamba_bin || true
    fi

    [ -f micromamba_bin ] && ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge python=3.9 pysam -y || true
    export PATH=\$PWD/env/bin:\$PATH

    # cyp21_paralog_pileup.py is staged in the work dir (Nextflow input); hotspot script imports it.

    python3 ${hotspot_py} \\
        --bam ${bam} \\
        --sites ${hotspots_tsv} \\
        -o ${sample_id}_cyp21_hotspot_pileup.tsv
    """
}
