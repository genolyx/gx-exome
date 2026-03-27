process FALLBACK_ANALYSIS {
    tag "$sample_id"
    label 'paraphase'
    publishDir "${params.outdir}/fallback", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path hba_bed
    path cyp21a2_bed
    path backbone_bed
    path ref_fasta
    path ref_index

    output:
    path "${sample_id}_hba_fallback.txt", emit: hba_report
    path "${sample_id}_cyp21a2_fallback.txt", emit: cyp21a2_report

    script:
    """
    # --- Dependency Setup (Micromamba) ---
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
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 freebayes=1.3.6 -y
    export PATH=\$PWD/env/bin:\$PATH
    
    # --- HBA Analysis ---
    
    awk '\$1=="chr16" || \$1=="16"' ${backbone_bed} > chr16_targets.bed
    
    samtools view -b -F 1024 ${bam} > no_dup.bam
    samtools index no_dup.bam
    
    samtools bedcov chr16_targets.bed no_dup.bam > chr16_cov.txt
    samtools bedcov ${hba_bed} no_dup.bam > hba_cov.txt

    freebayes \\
        --fasta-reference ${ref_fasta} \\
        --region chr16:172000-178000 \\
        --min-mapping-quality 0 \\
        --ploidy 4 \\
        --min-alternate-fraction 0.1 \\
        --min-alternate-count 3 \\
        no_dup.bam > ${sample_id}_hba_freebayes.vcf

    # Python HBA Analysis
    python3 -c "
import statistics
import sys

def parse_bedcov(filename):
    depths = []
    total_bases = 0
    total_len = 0
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            start, end, read_sum = int(parts[1]), int(parts[2]), int(parts[-1])
            length = end - start
            if length > 0:
                depth = read_sum / length
                depths.append(depth)
                total_bases += read_sum
                total_len += length
    return depths, total_bases, total_len

# Coverage Analysis
chr16_depths, _, _ = parse_bedcov('chr16_cov.txt')
median_cov = statistics.median(chr16_depths) if chr16_depths else 0

_, hba_sum, hba_len = parse_bedcov('hba_cov.txt')
hba_mean = hba_sum / hba_len if hba_len > 0 else 0

ratio = hba_mean / median_cov if median_cov > 0 else 0

status = 'PASS'
if ratio > 1.8:
    status = 'PASS (Normal 4-Copy Dosage: Negative for --SEA/--MED/--FIL/--THAI)'
elif ratio > 1.3:
    status = 'WARNING: Dosage suggests 3 Copies (Likely Silent Carrier -a/aa, e.g. -a3.7 or -a4.2)'
elif ratio > 0.8:
    status = 'WARNING: Dosage suggests 2 Copies (Likely Alpha Thal Trait --/aa or -a/-a)'
else:
    status = 'CRITICAL: Dosage suggests <=1 Copy (HbH Disease Risk --/-a)'

detected_variants = []
try:
    with open('${sample_id}_hba_freebayes.vcf', 'r') as matches:
        for line in matches:
            if line.startswith('#'): continue
            cols = line.strip().split('\\t')
            # CHROM POS ID REF ALT QUAL FILTER INFO ...
            chrom, pos, ref, alt, qual = cols[0], cols[1], cols[3], cols[4], float(cols[5])
            if int(pos) == 173598:
                 detected_variants.append(f'{chrom}:{pos} {ref}>{alt} (Qual:{qual:.1f}) - **POTENTIAL HB CONSTANT SPRING**')

except Exception as e:
    detected_variants.append(f'Error reading VCF: {e}')

with open('${sample_id}_hba_fallback.txt', 'w') as out:
    out.write(f'Sample: ${sample_id}\\n')
    out.write(f'HBA Interval Mean Depth: {hba_mean:.2f}\\n')
    out.write(f'Chr16 Median Target Depth: {median_cov:.2f}\\n')
    out.write(f'Ratio (HBA/Median): {ratio:.4f}\\n')
    out.write(f'Status: {status}\\n')
    out.write('\\n--- Targeted SNV Check (MAPQ 0 Enabled) ---\\n')
    if detected_variants:
        out.write('Hb Constant Spring Variant Detected:\\n')
        for v in detected_variants:
            out.write(f'  {v}\\n')
    else:
        out.write('No Hb Constant Spring variant (chr16:173598) detected.\\n')
    "

    # --- CYP21A2 Analysis ---
    
    awk '\$1=="chr6" || \$1=="6"' ${backbone_bed} > chr6_targets.bed

    samtools bedcov chr6_targets.bed no_dup.bam > chr6_cov.txt
    samtools bedcov ${cyp21a2_bed} no_dup.bam > cyp21a2_cov.txt

    # Python CYP21A2 Analysis
    python3 -c "
import statistics

def parse_bedcov(filename):
    depths = []
    total_bases = 0
    total_len = 0
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            start, end, read_sum = int(parts[1]), int(parts[2]), int(parts[-1])
            length = end - start
            if length > 0:
                depth = read_sum / length
                depths.append(depth)
                total_bases += read_sum
                total_len += length
    return depths, total_bases, total_len

# Chr6 Median
chr6_depths, _, _ = parse_bedcov('chr6_cov.txt')
if not chr6_depths:
    median_cov = 0
else:
    median_cov = statistics.median(chr6_depths)

# CYP21A2 Mean
_, gene_sum, gene_len = parse_bedcov('cyp21a2_cov.txt')
if gene_len > 0:
    gene_mean = gene_sum / gene_len
else:
    gene_mean = 0

# Ratio
if median_cov > 0:
    ratio = gene_mean / median_cov
else:
    ratio = 0

status = 'PASS'
if ratio < 0.7:
    status = 'WARNING: Possible Deletion Detected (Depth Proxy)'

with open('${sample_id}_cyp21a2_fallback.txt', 'w') as out:
    out.write(f'Sample: ${sample_id}\\n')
    out.write(f'CYP21A2 Interval Mean Depth: {gene_mean:.2f}\\n')
    out.write(f'Chr6 Median Target Depth: {median_cov:.2f}\\n')
    out.write(f'Ratio (CYP21A2/Median Chr6): {ratio:.4f}\\n')
    out.write(f'Status: {status}\\n')
    "
    """
}
