process PARAPHASE_RUN {
    tag "$sample_id"
    label 'paraphase'
    publishDir "${params.outdir}/pseudogene/paraphase", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}_paraphase.json"), emit: json
    tuple val(sample_id), path("*.vcf.gz"), emit: vcf, optional: true

    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    
    # Keep runtime env isolated
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    # Fix: Create directories to ensure they exist and are writable
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # Fix: Install SSL certs for Micromamba
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem

    # Fix: Use Micromamba for robust dependency resolution
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    # Fix: Disable SSL verify to bypass missing CA certs
    export MAMBA_SSL_VERIFY=false
    
    # Create temp environment with dependencies from Bioconda
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge python=3.9 pip samtools=1.16.1 minimap2=2.24 -y
    
    export PATH=\$PWD/env/bin:\$PATH

    # Install latest Paraphase via pip (supports --targeted)
    pip install paraphase --upgrade

    # Fix: Symlink BAM and COPY index to ensure perfect timestamp relationship
    # Using local names avoids confusion
    ln -s ${bam} input.bam
    cp ${bai} input.bam.bai

    # Debug: Check pysam access and list supported genes
    python3 -c "
import pysam
import os
import paraphase
from pathlib import Path
from collections import Counter

print('DEBUG: Paraphase version:', paraphase.__version__)
print('DEBUG: Paraphase path:', paraphase.__file__)

# Check genes in build 38
data_dir = Path(paraphase.__file__).parent / 'data' / '38'
if data_dir.exists():
    genes = [p.name for p in data_dir.iterdir() if p.is_dir()]
    print('DEBUG: Supported genes in data/38:', sorted(genes))
else:
    print('DEBUG: data/38 not found')

try:
    # Check BAM access
    bam = pysam.AlignmentFile('input.bam', 'rb')
    print('DEBUG: BAM has index?', bam.has_index())
    
    # Analyze GBA reads
    print('DEBUG: Analyzing GBA reads (chr1:155230000-155242500)...')
    mapqs = []
    flags = []
    count = 0
    for read in bam.fetch('chr1', 155230000, 155242500):
        count += 1
        mapqs.append(read.mapping_quality)
        flags.append(read.flag)
        if count < 5:
            print(f'DEBUG: Sample read: {read.query_name}, MAPQ: {read.mapping_quality}, Flag: {read.flag}')
            
    print(f'DEBUG: Total GBA reads: {count}')
    if count > 0:
        print(f'DEBUG: MAPQ distribution: {Counter(mapqs)}')
    
    # Analyze SMN reads
    print('DEBUG: Analyzing SMN reads (chr5:70925000-70954000)...')
    smn_count = 0
    smn_mapqs = []
    for read in bam.fetch('chr5', 70925000, 70954000):
        smn_count += 1
        smn_mapqs.append(read.mapping_quality)
    print(f'DEBUG: Total SMN reads: {smn_count}')
    if smn_count > 0:
        print(f'DEBUG: SMN MAPQ distribution: {Counter(smn_mapqs)}')

except Exception as e:
    print(f'DEBUG: Pysam error: {e}')
    "

    paraphase \
        -b input.bam \
        -r ${ref_fasta} \
        -g smn1,gba,pms2,strc,rccx,cyp21a2,vWF,hba1 \
        --targeted \
        --min-variant-frequency 0.05 \
        --write-nocalls-in-vcf \
        -p ${sample_id} \
        -o .
    
    echo "DEBUG: Finding all JSON files"
    find . -name "*.json"

    # Rename output output handling
    # If explicit name fails, try finding any json produced by paraphase
    if [ -f "${sample_id}.json" ]; then
        mv "${sample_id}.json" "${sample_id}_paraphase.json"
    else
        # Try to find any JSON that looks like output (excluding internal configs if any)
        # Using simple find default
        FOUND_JSON=\$(find . -maxdepth 1 -name "*.json" | head -n 1)
        if [ ! -z "\$FOUND_JSON" ]; then
            echo "DEBUG: Found alternative JSON: \$FOUND_JSON"
            mv "\$FOUND_JSON" "${sample_id}_paraphase.json"
        fi
    fi
    """
}

process SMACA_RUN {
    tag "$sample_id"
    label 'smaca'
    publishDir "${params.outdir}/pseudogene/smaca", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    // SMAca might need reference or other configs depending on version

    output:
    tuple val(sample_id), path("${sample_id}_smaca.txt"), emit: txt

    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    
    # Keep runtime env isolated
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & SMAca
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge smaca -y
    export PATH=\$PWD/env/bin:\$PATH

    # SMAca v1.2.3 usage
    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    smaca \\
        --reference hg38 \\
        --output ${sample_id}_smaca.txt \\
        ${bam}
    """
}

process PARAPHASE_RESCUE {
    tag "$sample_id"
    label 'paraphase'
    publishDir "${params.outdir}/rescue/paraphase", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path dark_genes_bed

    output:
    path "${sample_id}_paraphase.json", emit: json
    path "${sample_id}_paraphase.vcfs/*", emit: vcfs, optional: true
    path "${sample_id}_paraphase.bam", emit: bam, optional: true

    script:
    """
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    
    # Use SHARED cache for packages (speed up)
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    
    # Keep runtime env isolated
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    # Fix: Create directories to ensure they exist and are writable
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # Fix: Install SSL certs for Micromamba
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem

    # Fix: Use Micromamba for robust dependency resolution
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    # Fix: Disable SSL verify to bypass missing CA certs
    export MAMBA_SSL_VERIFY=false
    
    # Create temp environment with dependencies from Bioconda
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge python=3.9 pip samtools=1.16.1 minimap2=2.24 -y
    
    export PATH=\$PWD/env/bin:\$PATH

    # Install latest Paraphase via pip
    pip install paraphase --upgrade

    # --- PARAPHASE PATCHING FOR RESCUE MODE ---
    # Paraphase v3.4.0 hardcodes min_mapq=50 for realignment and min_mapq=5 for phasing.
    # To enable rescue of multi-mapping reads (MAPQ 0) in dark regions, we must patch the source.
    
    SITE_PACKAGES="env/lib/python3.9/site-packages"
    echo "Patching Paraphase at \$SITE_PACKAGES..."
    
    # Patch 1: Allow MAPQ 0 in BamRealigner AND Reduce min_aln for Short Reads
    if [ -f "\$SITE_PACKAGES/paraphase/prepare_bam_and_vcf.py" ]; then
        sed -i 's/min_mapq = 50/min_mapq = 0/g' \$SITE_PACKAGES/paraphase/prepare_bam_and_vcf.py
        sed -i 's/min_aln = 800/min_aln = 50/g' \$SITE_PACKAGES/paraphase/prepare_bam_and_vcf.py
        echo "Patched prepare_bam_and_vcf.py (min_mapq=50 -> 0, min_aln=800 -> 50)"
    else
        echo "WARNING: prepare_bam_and_vcf.py not found"
    fi
    
    # Patch 2: Allow MAPQ 0 in Phaser AND Disable Coverage Check
    if [ -f "\$SITE_PACKAGES/paraphase/phaser.py" ]; then
        # MapQ default: 5 -> 0
        sed -i 's/min_mapq=5/min_mapq=0/g' \$SITE_PACKAGES/paraphase/phaser.py
        # Coverage check: median <= 8 -> median <= -1 (disable check)
        sed -i 's/self.region_avg_depth.median <= 8/self.region_avg_depth.median <= -1/g' \$SITE_PACKAGES/paraphase/phaser.py
        
        
        # Patch 3: Force min_mapping_quality=0 in pileup calls using Python for safety
        cat <<EOF > patch_phaser.py
import os
import sys


site_packages = "\$SITE_PACKAGES"
paraphase_dir = os.path.join(site_packages, "paraphase")

if not os.path.exists(paraphase_dir):
    print(f"Error: {paraphase_dir} not found")
    sys.exit(1)

files_to_patch = []
for root, dirs, files in os.walk(paraphase_dir):
    for file in files:
        if file.endswith(".py"):
            files_to_patch.append(os.path.join(root, file))

for file_path in files_to_patch:
    print(f"Patching {file_path}...")
    with open(file_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    modified = False
    for i, line in enumerate(lines):
        clean_line = line.strip()
        
        # Check for 'truncate=True,'
        if "truncate=True," in clean_line:
            # Check if next line has min_base_quality
            next_has_base_qual = False
            if i + 1 < len(lines):
                if "min_base_quality" in lines[i+1]:
                    next_has_base_qual = True
            
            # If NO min_base_quality follows, we must inject min_mapping_quality here
            if not next_has_base_qual:
                if "min_mapping_quality" not in line:
                    new_lines.append(line.replace("truncate=True,", "truncate=True, min_mapping_quality=0,"))
                    modified = True
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)
                 
        # Check for 'min_base_quality='
        elif "min_base_quality=" in clean_line:
            # Inject min_mapping_quality here
            if "min_mapping_quality" not in line:
                new_lines.append(line.replace("min_base_quality=self.MEAN_BASE_QUAL,", "min_base_quality=self.MEAN_BASE_QUAL, min_mapping_quality=0,"))
                modified = True
            else:
                new_lines.append(line)
                
        else:
            new_lines.append(line)
    
    if modified:
        with open(file_path, "w") as f:
            f.writelines(new_lines)
        print(f"  - Patched {file_path}")


print("Successfully patched phaser.py with safe Python script")
EOF
        python3 patch_phaser.py
        
        echo "Patched phaser.py (min_mapq=5 -> 0, disabled coverage check, forced pileup min_mapq=0)"
    else
         echo "WARNING: phaser.py not found"
    fi

    # --- CONFIG GENERATION ---
    # Create a custom config file that maps the Dark Genes BED regions to Paraphase gene keys.
    
    cat <<EOF > generate_config.py
import sys
import yaml
import os
import paraphase

# Load default config from installed package
data_path = os.path.join(os.path.dirname(paraphase.__file__), "data", "38")
try:
    with open(os.path.join(data_path, "config.yaml"), "r") as f:
        config = yaml.safe_load(f)
except Exception as e:
    print(f"Error loading default config: {e}")
    sys.exit(1)

# Mapping from BED name column to Paraphase Config Keys
mapping = {
    'SMN1_SMN2_Region': 'smn1',
    'GBA_GBAP1_Region': 'GBA',
    'CYP21A2_RCCX_Region': 'rccx',
    'HBA1_HBA2_Region': 'hba',
    'CYP2D6_Region': 'CYP2D6'
}

# Parse BED file
new_regions = {}
try:
    with open("${dark_genes_bed}", "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.strip().split()
            if len(parts) >= 4:
                chrom, start, end, name = parts[0], parts[1], parts[2], parts[3]
                if name in mapping:
                    key = mapping[name]
                    new_regions[key] = f"{chrom}:{start}-{end}"
except Exception as e:
    print(f"Error reading BED file: {e}")
    sys.exit(1)

# Create Rescue Config
rescue_config = {}
for key, region in new_regions.items():
    if key in config:
        # Copy existing settings
        rescue_config[key] = config[key]
        # OVERRIDE regions with our continuous blocks
        rescue_config[key]['realign_region'] = region
        rescue_config[key]['extract_regions'] = region
        
        # Remove gene2_region to ensure we only use our defined continuous block (EXCEPT for GBA)
        # GBA needs gene2_region (GBAP1) for gene conversion analysis
        if 'gene2_region' in rescue_config[key]:
             if key == 'GBA':
                 print(f"Preserved gene2_region for {key} (L444P analysis)")
                 # Ensure extract_regions covers both? 
                 # We set extract_regions = region (which is the BIG block from BED).
                 # Paraphase might append gene2_region to extraction list if it exists?
                 # Or does it trust 'extract_regions' we set?
                 # We set extract_regions = region (155210000-155245000). This covers BOTH.
                 # So preserving gene2_region logic allows internal mapping.
             else:
                 del rescue_config[key]['gene2_region']
                 print(f"Configured {key}: {region} (removed gene2_region)")
        else:
             print(f"Configured {key}: {region}")
             
        # Remove explicit boundaries to allow them to default to the new realign_region
        # This prevents 0 coverage due to wide default boundaries vs narrow target regions
        for b_key in ['left_boundary', 'right_boundary', 'gene_start', 'gene_end']:
            if b_key in rescue_config[key]:
                del rescue_config[key][b_key]
                print(f"Removed {b_key} from {key}")

    else:
        print(f"Warning: {key} not found in default config keys")

with open("rescue_config.yaml", "w") as f:
    yaml.dump(rescue_config, f)
EOF

    # Run config generation
    python3 generate_config.py

    # Determine genes to run from the generated config
    GENES_TO_RUN=\$(python3 -c "import yaml; print(','.join(yaml.safe_load(open('rescue_config.yaml')).keys()))")
    
    echo "Starting Paraphase Rescue for genes: \$GENES_TO_RUN"
    
    # Create output directory
    mkdir -p ${sample_id}_paraphase.vcfs

    # Check for BAM index
    # Paraphase strictly requires bam.bai
    if [ ! -f "${bam}.bai" ]; then
        cp ${bai} ${bam}.bai
    fi

    # Run Paraphase with Custom Config and Patched Code
    paraphase \
        -b ${bam} \
        -r ${ref_fasta} \
        -o . \
        -c rescue_config.yaml \
        -g \$GENES_TO_RUN \
        --targeted \
        --min-variant-frequency 0.05 \
        --write-nocalls-in-vcf

    # Rename output for consistency
    
    # JSON
    FOUND_JSON=\$(find . -maxdepth 1 -name "*paraphase.json" | head -n 1)
    if [ ! -z "\$FOUND_JSON" ]; then
        mv "\$FOUND_JSON" "${sample_id}_paraphase.json"
    fi

    # BAM (and index)
    FOUND_BAM=\$(find . -maxdepth 1 -name "*paraphase.bam" | head -n 1)
    if [ ! -z "\$FOUND_BAM" ]; then
        mv "\$FOUND_BAM" "${sample_id}_paraphase.bam"
        # Check for index (standard naming or with .bai extension)
        if [ -f "\$FOUND_BAM.bai" ]; then
             mv "\$FOUND_BAM.bai" "${sample_id}_paraphase.bam.bai"
        elif [ -f "\${FOUND_BAM%.bam}.bai" ]; then
             mv "\${FOUND_BAM%.bam}.bai" "${sample_id}_paraphase.bam.bai"
        fi
    fi
    
    # VCFs
    # Identify generated dir
    FOUND_VCF_DIR=\$(find . -maxdepth 1 -type d -name "*_paraphase_vcfs" | head -n 1)
    TARGET_VCF_DIR="${sample_id}_paraphase.vcfs"
    
    if [ ! -z "\$FOUND_VCF_DIR" ] && [ "\$FOUND_VCF_DIR" != "./\$TARGET_VCF_DIR" ] && [ "\$FOUND_VCF_DIR" != "\$TARGET_VCF_DIR" ]; then
         # If target exists and is empty, remove it (we created it with mkdir earlier)
         rmdir \$TARGET_VCF_DIR 2>/dev/null || true
         # Move if still different
         if [ -d "\$FOUND_VCF_DIR" ]; then
            mv "\$FOUND_VCF_DIR" "\$TARGET_VCF_DIR"
         fi
    fi
    """
}



process CONSOLIDATE_RESCUE {
    publishDir "${params.outdir}/rescue", mode: 'copy'
    
    input:
    path paraphase_jsons
    path hba_results
    path cyp21a2_results
    path eh_vcfs

    output:
    path "summary_report.txt"
    path "detailed_report.txt"

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

    # Install Micromamba & Python
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge python=3.9 -y
    export PATH=\$PWD/env/bin:\$PATH

    # --- Python Aggregation Script ---
    cat <<EOF > consolidate.py
import json
import os
import glob
import re

def parse_paraphase(files):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_paraphase.json', '')
        try:
            with open(f) as json_file:
                data = json.load(json_file)
                # Structure: { "gene": { "total_cn": X, "final_haplotypes": [...] } }
                results[sample] = {}
                if not data:
                    results[sample] = "No Calls (Low Quality/Coverage)"
                    continue

                for gene, info in data.items():
                    summary = []
                    # Copy Numbers
                    if 'total_cn' in info:
                        summary.append(f"{gene}_CN={info['total_cn']}")
                    if 'smn1_cn' in info:
                        summary.append(f"SMN1={info['smn1_cn']}")
                    if 'smn2_cn' in info:
                        summary.append(f"SMN2={info['smn2_cn']}")
                    
                    # Haplotypes details
                    if 'final_haplotypes' in info:
                         haps = info['final_haplotypes']
                         hap_names = [h['haplotype'] for h in haps if 'haplotype' in h]
                         summary.append(f"Haps: {', '.join(hap_names)}")
                    
                    results[sample][gene] = "; ".join(summary)
        except Exception as e:
            results[sample] = f"Error: {e}"
    return results

def parse_fallback(files, type_label):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace(f'_{type_label}_fallback.txt', '')
        try:
            with open(f) as txt_file:
                content = txt_file.readlines()
                # Parse the key metrics directly from file content
                # Format: 
                # CYP21A2 Interval Mean Depth: 10.28
                # Chr6 Median Target Depth: 55.00
                # Ratio (CYP21A2/Median Chr6): 0.1870
                # Status: ...
                
                details = []
                status = "Unknown"
                ratio = "N/A"
                
                for line in content:
                    line = line.strip()
                    if "Ratio" in line:
                        ratio = line.split(':')[-1].strip()
                    if "Status:" in line:
                        status = line.split('Status:')[-1].strip()
                    if "CN:" in line or "Carrier" in line or "Detected" in line:
                         details.append(line)

                # Construct a more detailed string
                if type_label == 'cyp21a2':
                    results[sample] = f"Ratio={ratio} | {status}"
                else:
                    # Generic fallback (HBA)
                    if not details: 
                        results[sample] = "No variants detected"
                    else:
                        results[sample] = "; ".join(details)
                        
        except Exception as e:
             results[sample] = f"Error: {e}"
    return results

def parse_eh(files):
    results = {}
    for f in files:
        sample = os.path.basename(f).replace('_eh.vcf', '')
        try:
            summary = []
            with open(f) as vcf:
                genes_found = False
                for line in vcf:
                    if line.startswith('#'): continue
                    parts = line.split('\\t')
                    if len(parts) >= 10:
                        info = parts[7]
                        fmt_keys = parts[8].split(':')
                        fmt_vals = parts[9].strip().split(':')
                        
                        # Extract Gene Name from INFO (REPID=FMR1)
                        repid_match = re.search(r'REPID=([^;]+)', info)
                        if repid_match:
                            gene = repid_match.group(1)
                            genes_found = True
                            
                            # Extract Repeat Counts (REPCN) and CI (REPCI)
                            # FORMAT: GT:SO:REPCN:REPCI...
                            if 'REPCN' in fmt_keys:
                                idx = fmt_keys.index('REPCN')
                                cn = fmt_vals[idx]
                                
                                # Optional: Add Confidence Interval if available
                                ci_str = ""
                                if 'REPCI' in fmt_keys:
                                    ci_idx = fmt_keys.index('REPCI')
                                    ci_str = f" (CI: {fmt_vals[ci_idx]})"
                                    
                                summary.append(f"{gene} Repeats: {cn}{ci_str}")
            
            if not summary:
                 if genes_found: results[sample] = "No expansions called (Ref)"
                 else: results[sample] = "No catalog matches"
            else:
                 results[sample] = " | ".join(summary)
                 
        except Exception as e:
            results[sample] = f"Error: {e}"
    return results

# --- Main Execution ---

# 1. Gather all data
paraphase_files = glob.glob("*_paraphase.json")
hba_files = glob.glob("*_hba_fallback.txt")
cyp21a2_files = glob.glob("*_cyp21a2_fallback.txt")
eh_files = glob.glob("*_eh.vcf")

p_results = parse_paraphase(paraphase_files)
h_results = parse_fallback(hba_files, "hba")
c_results = parse_fallback(cyp21a2_files, "cyp21a2")
e_results = parse_eh(eh_files)

# 2. Get all samples
all_samples = set(p_results.keys()) | set(h_results.keys()) | set(c_results.keys()) | set(e_results.keys())

# 3. Generating Summary Report (Table)
with open("summary_report.txt", "w") as out:
    header = "Sample\\tParaphase(SMN/GBA)\\tHBA_Analysis\\tCYP21A2_Analysis\\tFragileX(FMR1)"
    out.write(header + "\\n")
    print(header)
    
    for sample in sorted(all_samples):
        # Format columns - simplified for table
        p_str = str(p_results.get(sample, "N/A"))[:50] + "..." if len(str(p_results.get(sample, "N/A"))) > 50 else str(p_results.get(sample, "N/A"))
        h_str = h_results.get(sample, "N/A")
        c_str = c_results.get(sample, "N/A")
        e_str = e_results.get(sample, "N/A")
        
        line = f"{sample}\\t{p_str}\\t{h_str}\\t{c_str}\\t{e_str}"
        out.write(line + "\\n")
        print(line)

# 4. Detailed Report (Verbose)
with open("detailed_report.txt", "w") as out:
    for sample in sorted(all_samples):
        out.write(f"SAMPLE: {sample}\\n")
        out.write("="*40 + "\\n")
        
        out.write("PARAPHASE RESULTS (SMN1/2, GBA):\\n")
        # Pretty print the dictionary or string
        p_val = p_results.get(sample, "No Data")
        if isinstance(p_val, dict):
            for g, res in p_val.items():
                out.write(f"  - {g}: {res}\\n")
        else:
            out.write(f"  {p_val}\\n")
        out.write("\\n")

        out.write("HBA ANALYSIS (Alpha Thalassemia):\\n")
        out.write(f"  {h_results.get(sample, 'No Data')}\\n\\n")
        
        out.write("CYP21A2 ANALYSIS (CAH):\\n")
        out.write(f"  {c_results.get(sample, 'No Data')}\\n\\n")

        out.write("EXPANSION HUNTER (Fragile X / FMR1):\\n")
        out.write(f"  {e_results.get(sample, 'No Data')}\\n")
        
        out.write("\\n" + "-"*40 + "\\n\\n")
EOF

    python3 consolidate.py
    """
}
