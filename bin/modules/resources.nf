process PREPARE_VIZ_RESOURCES {
    tag "setup_resources"
    label 'local'
    publishDir "${params.outdir}/resources", mode: 'copy'

    input:
    path ref_fasta
    path ref_fai

    output:
    path "viz_env", emit: env
    path "refGene.sorted.gtf.gz", emit: gtf
    path "refGene.sorted.gtf.gz.tbi", emit: gtf_index
    path "smn_ref_wide.fa", emit: smn_ref
    path "smn_ref_wide.fa.*", emit: smn_index

    script:
    """
    # --- 1. Setup Shared Environment ---
    # We build the environment ONCE here.
    
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
    
    echo "Creating Shared Visualization Environment..."
    # Install: igv-reports, htslib (tabix), bwa, samtools
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./viz_env -c bioconda -c conda-forge python=3.9 igv-reports pysam htslib bwa samtools -y
    
    # Activate for use in this script
    export PATH=\$PWD/viz_env/bin:\$PATH

    # --- 2. Download and Process RefSeq GTF ---
    echo "Downloading RefSeq GTF..."
    wget -q --no-check-certificate -O refGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
    
    # Sort and Index (using installed tabix)
    echo "Indexing GTF..."
    gzip -d -c refGene.gtf.gz | grep -v "^#" | sort -k1,1 -k4,4n | bgzip > refGene.sorted.gtf.gz
    tabix -p gff refGene.sorted.gtf.gz

    # --- 3. Prepare SMN Mini-Reference ---
    # Extract chr5:70920000-70960000 (Widened for Silent Carrier SNP g.27134)
    echo "Generating SMN Mini-Reference..."
    SMN_START=70920000
    SMN_END=70960000
    
    # Access ref using fai
    samtools faidx ${ref_fasta} chr5:\${SMN_START}-\${SMN_END} > smn_ref_wide.fa
    
    # Rename header for local alignment
    sed -i "s/>.*/ >chr5/" smn_ref_wide.fa
    
    # Index
    bwa index smn_ref_wide.fa
    """
}
