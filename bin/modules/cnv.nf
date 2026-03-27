process PREPROCESS_INTERVALS {
    tag "intervals"
    label 'gatk'
    
    input:
    path backbone_bed
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    path "preprocessed.interval_list", emit: interval_list

    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    
    # Keep runtime env isolated
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    
    # Ensure cache dir exists
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk PreprocessIntervals \\
        -L ${backbone_bed} \\
        -R ${ref_fasta} \\
        --bin-length 0 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O preprocessed.interval_list
    """
}

process COLLECT_READ_COUNTS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/cnv/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path interval_list
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.counts.tsv"), emit: counts // Request TSV for readability check
    tuple val(sample_id), path("${sample_id}.counts.hdf5"), emit: counts_hdf5

    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk CollectReadCounts \\
        -I ${bam} \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${ref_fasta} \\
        --format TSV \\
        -O ${sample_id}.counts.tsv

    gatk CollectReadCounts \\
        -I ${bam} \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${ref_fasta} \\
        -O ${sample_id}.counts.hdf5
    """
}

process GCNV_CLARITY {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/cnv", mode: 'copy'

    input:
    tuple val(sample_id), path(counts)
    path pon_tar
    path interval_list

    output:
    tuple val(sample_id), path("${sample_id}_gcnv_calls"), emit: gcnv_calls

    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 gcnvkernel fastprogress=1.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    tar -xf ${pon_tar}

    # 2. Determine Germline Contig Ploidy (using PoN model)
    gatk DetermineGermlineContigPloidy \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --input ${counts} \\
        --model ploidy-model \\
        --output . \\
        --output-prefix ${sample_id} \\
        --verbosity DEBUG

    # 3. Germline CNV Caller
    gatk GermlineCNVCaller \\
        --run-mode CASE \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --input ${counts} \\
        --contig-ploidy-calls ${sample_id}-calls \\
        --model gcnv-model \\
        --output ${sample_id}_gcnv_calls \\
        --output-prefix ${sample_id} \\
        --verbosity DEBUG
    """
}

process ANNOTATE_INTERVALS {
    tag "annotate"
    label 'gatk'
    
    input:
    path interval_list
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    path "annotated.interval_list", emit: annotated_intervals

    script:
    """
    # --- Dependency Setup ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk AnnotateIntervals \\
        -L ${interval_list} \\
        -R ${ref_fasta} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O annotated.interval_list
    """
}

process GCNV_COHORT_RUN {
    tag "cohort_N\${counts_hdf5.size()}"
    label 'gatk'
    publishDir "${params.outdir}/cnv/cohort", mode: 'copy'

    input:
    path counts_hdf5 // List of all .hdf5 files
    path annotated_intervals
    path interval_list

    output:
    path "gcnv-calls/*", emit: calls
    path "gcnv-model", emit: model
    path "ploidy-model", emit: ploidy_model
    path "ploidy-calls", emit: ploidy_calls

    script:
    """
    # --- Dependency Setup ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    
    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 gcnvkernel fastprogress=1.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH
    
    INPUT_ARGS=""
    for f in *.hdf5; do
        INPUT_ARGS="\$INPUT_ARGS -I \$f"
    done
    
    python3 -c "
import sys

def get_priors(contig):
    c = contig.lower()
    if 'chrx' in c or (c == 'x'):
        return '0.01\\t0.49\\t0.49\\t0.01'
    elif 'chry' in c or (c == 'y'):
        return '0.49\\t0.49\\t0.01\\t0.01'
    elif 'chrm' in c or (c == 'mt') or (c == 'm'):
        return '0.01\\t0.01\\t0.97\\t0.01'
    else:
        return '0.01\\t0.01\\t0.97\\t0.01'

with open('${interval_list}', 'r') as f_in, open('contig_ploidy_priors.tsv', 'w') as f_out:
    f_out.write('CONTIG_NAME\\tPLOIDY_PRIOR_0\\tPLOIDY_PRIOR_1\\tPLOIDY_PRIOR_2\\tPLOIDY_PRIOR_3\\n')
    for line in f_in:
        if line.startswith('@SQ'):
            parts = line.strip().split()
            for p in parts:
                if p.startswith('SN:'):
                    contig = p.split(':')[1]
                    priors = get_priors(contig)
                    f_out.write(f'{contig}\\t{priors}\\n')
"

    # 1. Determine Cohort Ploidy
    mkdir -p ploidy-model
    gatk DetermineGermlineContigPloidy \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        \$INPUT_ARGS \\
        --contig-ploidy-priors contig_ploidy_priors.tsv \\
        --mapping-error-rate 0.2 \\
        --output . \\
        --output-prefix ploidy \\
        --verbosity DEBUG
        
    # 2. Cohort CNV Calling
    mkdir -p gcnv-model
    mkdir -p cohort-calls
    
    gatk GermlineCNVCaller \\
        --run-mode COHORT \\
        -L ${interval_list} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        \$INPUT_ARGS \\
        --contig-ploidy-calls ploidy-calls \\
        --annotated-intervals ${annotated_intervals} \\
        --output . \\
        --output-prefix gcnv \\
        --verbosity DEBUG
        
    mv gcnv-calls/* gcnv-calls/ 2>/dev/null || true
    """
}

process GCNV_CASE_RUN {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/cnv/case/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(counts_hdf5)
    path annotated_intervals
    path interval_list
    path gcnv_model
    path ploidy_model // This is the ploidy MODEL directory from cohort run

    output:
    tuple val(sample_id), path("${sample_id}-calls"), emit: calls
    tuple val(sample_id), path("${sample_id}-ploidy-calls"), emit: ploidy_calls

    script:
    """
    # --- Dependency Setup ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 gcnvkernel fastprogress=1.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    mkdir -p ${sample_id}-ploidy-calls
    mkdir -p ${sample_id}-calls

    # 1. Determine Sample Ploidy (Case Mode)
    gatk DetermineGermlineContigPloidy \\
        --model ${ploidy_model} \\
        -I ${counts_hdf5} \\
        -O ${sample_id}-ploidy-calls \\
        --output-prefix ploidy \\
        --verbosity DEBUG

    # Fix nesting caused by --output-prefix:
    # Moves ${sample_id}-ploidy-calls/ploidy-calls/SAMPLE_0 -> ${sample_id}-ploidy-calls/SAMPLE_0
    mv ${sample_id}-ploidy-calls/ploidy-calls/* ${sample_id}-ploidy-calls/
    rmdir ${sample_id}-ploidy-calls/ploidy-calls

    # 2. Germline CNV Caller (Case Mode)
    gatk GermlineCNVCaller \\
        --run-mode CASE \\
        --model ${gcnv_model} \\
        -I ${counts_hdf5} \\
        --contig-ploidy-calls ${sample_id}-ploidy-calls \\
        -O ${sample_id}-calls \\
        --output-prefix ${sample_id} \\
        --verbosity DEBUG
        
    # Fix nesting caused by --output-prefix in GermlineCNVCaller
    if [ -d "${sample_id}-calls/${sample_id}-calls" ]; then
        mv ${sample_id}-calls/${sample_id}-calls/* ${sample_id}-calls/
        rmdir ${sample_id}-calls/${sample_id}-calls
    fi

    # Copy interval_list.tsv for POSTPROCESS consistency
    cp ${gcnv_model}/interval_list.tsv ${sample_id}-calls/interval_list.tsv
    """
}

process POSTPROCESS_GCNV {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/cnv/segments", mode: 'copy'

    input:
    tuple val(sample_id), path(calls_path), path(ploidy_calls_path)
    path model_path
    path ref_fasta
    path ref_fai
    path ref_dict
    path interval_list
    path gcnv_configs

    output:
    tuple val(sample_id), path("${sample_id}_cnv.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}_cnv.vcf.gz.tbi"), emit: index
    path "${sample_id}_segments.vcf.gz", emit: segments

    script:
    """
    # --- Dependency Setup ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 gcnvkernel fastprogress=1.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    # Prepare working directory base
    mkdir -p work_calls
    mkdir -p work_ploidy

    # Copy provided configs from Cohort run
    # gcnv_configs is passed as a space-separated list of paths by Nextflow
    if [ -n "${gcnv_configs}" ]; then
       # Use a loop or simple cp, filtering for existence
       for f in ${gcnv_configs}; do
         if [ -f "\$f" ]; then
           cp "\$f" work_calls/
         fi
       done
    fi

    # 1. Setup Calls
    echo "DEBUG: Searching for SAMPLE_0 in ${calls_path}..."
    
    if [ -d "${calls_path}/SAMPLE_0" ]; then
        echo "DEBUG: Found SAMPLE_0 at root."
        SHARD_INDEX=0
        mkdir -p work_calls/SAMPLE_0
        cp -r ${calls_path}/SAMPLE_0/* work_calls/SAMPLE_0/
        cp ${calls_path}/*.json work_calls/ 2>/dev/null || true
        
    else
        echo "DEBUG: Validating nested search..."
        CASE_NESTED_PATH=\$(find ${calls_path} -name "SAMPLE_0" -type d | head -n 1)
        
        if [ -n "\$CASE_NESTED_PATH" ]; then
            echo "DEBUG: Found nested SAMPLE_0 at \$CASE_NESTED_PATH"
            SHARD_INDEX=0
            mkdir -p work_calls/SAMPLE_0
            cp -r \$CASE_NESTED_PATH/* work_calls/SAMPLE_0/
            
            CASE_PARENT=\$(dirname \$CASE_NESTED_PATH)
            cp \$CASE_PARENT/*.json work_calls/ 2>/dev/null || true
            
        elif [ -f "${calls_path}/sample_name.txt" ]; then
            # COHORT MODE
            DIR_NAME=\$(basename ${calls_path})
            SHARD_INDEX=\${DIR_NAME#SAMPLE_}
            
            mkdir -p work_calls/SAMPLE_\${SHARD_INDEX}
            cp -r ${calls_path}/* work_calls/SAMPLE_\${SHARD_INDEX}/
        else
            echo "Error: Could not locate valid calls structure in ${calls_path}"
            ls -R ${calls_path}
            exit 1
        fi
    fi
    
    # 2. Setup Ploidy (Must match theCalls Index)
    mkdir -p work_ploidy/SAMPLE_\${SHARD_INDEX}
    
    # Try to find the matching SAMPLE_X in the ploidy path
    if [ -d "${ploidy_calls_path}/SAMPLE_\${SHARD_INDEX}" ]; then
         cp -r ${ploidy_calls_path}/SAMPLE_\${SHARD_INDEX}/* work_ploidy/SAMPLE_\${SHARD_INDEX}/
    elif [ "\${SHARD_INDEX}" == "0" ] && [ -d "${ploidy_calls_path}/SAMPLE_0" ]; then
         # Fallback for Case Mode if path structure is slightly different (e.g. nested)
         cp -r ${ploidy_calls_path}/SAMPLE_0/* work_ploidy/SAMPLE_0/
    else
         echo "Error: Could not find matching ploidy calls for Index \${SHARD_INDEX} in ${ploidy_calls_path}"
         ls -R ${ploidy_calls_path}
         exit 1
    fi

    # Copy interval_list.tsv from the MODEL directory to the SHARD directory root
    # This guarantees it matches exactly what the model was built with
    cp ${model_path}/interval_list.tsv work_calls/interval_list.tsv

    gatk PostprocessGermlineCNVCalls \\
        --model-shard-path ${model_path} \\
        --calls-shard-path work_calls \\
        --allosomal-contig chrX --allosomal-contig chrY \\
        --contig-ploidy-calls work_ploidy \\
        --sample-index \${SHARD_INDEX} \\
        --output-genotyped-intervals ${sample_id}_cnv.vcf.gz \\
        --output-genotyped-segments ${sample_id}_segments.vcf.gz \\
        --output-denoised-copy-ratios ${sample_id}_denoised_copy_ratios.tsv \\
        --sequence-dictionary ${ref_dict}
    """
}
