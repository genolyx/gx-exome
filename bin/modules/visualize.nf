process GENERATE_VISUAL_EVIDENCE {
    tag "$sample_id"
    label 'visualize'
    publishDir "${params.outdir}/snapshots", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(vcf_idx)
    path eh_images
    path ref_fasta
    path ref_fai
    path viz_env
    path gtf
    path gtf_idx
    path smn_ref
    path smn_indices
    
    output:
    path "*_visual_report.html", emit: snapshots

    script:
    """
    # --- Shared Environment Setup ---
    # Use the pre-built environment passed from PREPARE_VIZ_RESOURCES
    export PATH=\$PWD/viz_env/bin:\$PATH

    # Symlink Indices
    if [ ! -f "${bam}.bai" ]; then ln -s ${bai} ${bam}.bai; fi
    if [ ! -f "${ref_fasta}.fai" ]; then ln -s ${ref_fai} ${ref_fasta}.fai; fi

    # First region = default IGV panel: c.840 discriminative base (C vs T on unified + hg38 tracks)
    printf "chr5\\t70951920\\t70951970\\tSMA_c840_SNP_Zoom\\n" > regions.bed
    printf "chr5\\t70952050\\t70952100\\tSMA_Silent_Carrier_SNP\\n" >> regions.bed
    printf "chr16\\t173300\\t173800\\tHb_Constant_Spring\\n" >> regions.bed
    printf "chr16\\t160000\\t190000\\tAlpha_Thal_Structural\\n" >> regions.bed
    printf "chr1\\t155232540\\t155240432\\tGBA_Locus\\n" >> regions.bed
    printf "chr6\\t32038000\\t32045000\\tCYP21A2_Locus\\n" >> regions.bed
    printf "chrX\\t147911919\\t147912110\\tFragileX_FMR1\\n" >> regions.bed

    # --- SMN Combined Pileup Strategy ---
    # Resources (smn_ref.fa and indices) provided by PREPARE_VIZ_RESOURCES
    SMN_START=70920000
    
    # 2. Extract Reads from BOTH SMN1 and SMN2 loci
    # Widened windows to include g.27134 (Silent Carrier SNP)
    samtools view -b ${bam} "chr5:70920000-70960000" "chr5:70020000-70080000" > smn_raw_reads.bam
    samtools fastq smn_raw_reads.bam > smn_reads.fq
    
    # 3. Re-align to SMN1 ref
    bwa mem -t 2 ${smn_ref} smn_reads.fq | samtools view -b - > smn_realigned_local.bam
    
    # 4. Shift Coordinates Back to hg38
    python3 -c "
import pysam
infile = pysam.AlignmentFile('smn_realigned_local.bam', 'rb')
outfile = pysam.AlignmentFile('smn_combined_hg38.bam', 'wb', template=infile)
offset = \$SMN_START - 1 

for read in infile:
    if not read.is_unmapped:
        read.reference_start += offset
        if not read.mate_is_unmapped and read.next_reference_name == 'chr5':
             read.next_reference_start += offset
    outfile.write(read)
infile.close()
outfile.close()
"
    samtools sort -o smn_combined_final.bam smn_combined_hg38.bam
    samtools index smn_combined_final.bam

    # RefSeq GTF is already provided (staged as refGene.sorted.gtf.gz or similar)
    # Symlink to strict names for config
    ln -s ${gtf} refGene.gtf.gz
    ln -s ${gtf_idx} refGene.gtf.gz.tbi

    # --- Prepare Fragile X Image Track (Popup) ---
    # --- Prepare Fragile X Image Track (Popup) ---
    fx_track_logic=""
    eh_img=\$(find . -maxdepth 1 -name "${sample_id}*.svg" | head -n 1)
    
    if [ -n "\$eh_img" ]; then
        echo "Found Fragile X Image, embedding in GFF..."
        
        # 1. Create GFF with Button Trigger (Short string, safe)
        cat <<EOF > create_fx_track.py
with open('fx_image.gff3', 'w') as f:
    f.write('##gff-version 3\\n')
    import urllib.parse
    
    # JavaScript to show the hidden div and scroll to it
    js_action = "document.getElementById('fx_graph_container').style.display='block'; document.getElementById('fx_graph_container').scrollIntoView({behavior: 'smooth'});"
    
    btn_html = f'<button onclick="{js_action}" style="background:orange; color:white; font-weight:bold; padding:5px; cursor:pointer;">SHOW GRAPH AT BOTTOM</button>'
    encoded_note = urllib.parse.quote(btn_html)
    
    f.write(f'chrX\\tREViewer\\tTrigger\\t147911919\\t147912110\\t.\\t.\\t.\\tID=FX1;Name=Click to View;Note={encoded_note}\\n')
EOF
        python3 create_fx_track.py
        
        if [ -f "fx_image.gff3" ]; then
            touch fx_track_ready
        fi
    fi

    # --- Build Track Config (Dynamic JSON) ---
    cat <<EOF > build_config.py
import json
import os

tracks = []

# 1. BAM Track
tracks.append({
    "type": "alignment",
    "format": "bam",
    "url": "${bam}",
    "indexURL": "${bam}.bai",
    "name": "${sample_id}",
    "displayMode": "EXPANDED", 
    "colorBy": "base",
    "height": 500,
    "minMappingQuality": 0,
    "samplingDepth": 1000000,
    "showSoftClips": True,
    "visibilityWindow": -1
})

# 2. SMN Combined Pileup (Unified Alignment)
tracks.append({
    "type": "alignment",
    "format": "bam",
    "url": "smn_combined_final.bam",
    "indexURL": "smn_combined_final.bam.bai",
    "name": "SMN1+SMN2 Pileup (Unified)",
    "displayMode": "EXPANDED", 
    "colorBy": "base",
    "height": 500,
    "minMappingQuality": 0,
    "samplingDepth": 1000000,
    "showSoftClips": True,
    "visibilityWindow": -1
})

# 3. RefSeq Genes (Detailed)
tracks.append({
    "type": "annotation",
    "format": "gtf",
    "url": "refGene.gtf.gz",
    "indexURL": "refGene.gtf.gz.tbi",
    "name": "RefSeq Genes (Detail)",
    "displayMode": "EXPANDED",
    "visibilityWindow": -1
})

# 4. RefSeq Genes (Structure Only)
tracks.append({
    "type": "annotation",
    "format": "gtf",
    "url": "refGene.gtf.gz",
    "indexURL": "refGene.gtf.gz.tbi",
    "name": "RefSeq Genes (Structure)",
    "displayMode": "SQUISHED",
    "color": "darkblue",
    "visibilityWindow": -1
})

# 5. Fragile X Image Track (Conditional)
if os.path.exists('fx_track_ready') and os.path.exists('fx_image.gff3'):
    tracks.append({
        "type": "annotation",
        "format": "gff3",
        "url": "fx_image.gff3",
        "name": "Fragile X (Graph Trigger)",
        "displayMode": "EXPANDED",
        "color": "orange",
        "height": 40,
        "visibilityWindow": -1
    })


# 3. VCF Track
vcf_file = "${vcf}"
if vcf_file and vcf_file != "null" and os.path.exists(vcf_file):
    # Ensure index link
    if not os.path.exists(vcf_file + ".tbi") and os.path.exists("${vcf_idx}"):
         os.symlink("${vcf_idx}", vcf_file + ".tbi")
         
    tracks.append({
        "type": "variant",
        "format": "vcf",
        "url": vcf_file,
        "indexURL": vcf_file + ".tbi",
        "name": "Called Variants",
        "displayMode": "SQUISHED",
        "color": "purple",
        "visibilityWindow": -1
    })

with open("track_config.json", "w") as f:
    json.dump(tracks, f, indent=2)
EOF

    python3 build_config.py

    echo "Generating IGV Report..."
    
    create_report regions.bed \\
        --fasta ${ref_fasta} \\
        --genome hg38 \\
        --track-config track_config.json \\
        --flanking 1000 \\
        --output temp_report.html \\
        --title "Visual Evidence: ${sample_id}"

    # --- Alpha thalassemia review note (IGV report footer) ---
    python3 -c "
html = open('temp_report.html').read()
alpha_box = '''<div id='alpha_thal_review_note' style='margin-top:24px;padding:16px;border:1px solid #2c5282;border-radius:8px;background:#f0f7ff;max-width:900px;'>
<h3 style='margin:0 0 8px 0;'>Alpha thalassemia (review)</h3>
<p style='margin:0;line-height:1.5;color:#1a365d;'>The fallback analysis runs Freebayes on <code>chr16:172000–178000</code> with ploidy 4 mainly for Hb Constant Spring at a fixed position — it is not a full α-thal SNV/indels genotype caller for all carriers.</p>
</div>'''
html = html.replace('</body>', alpha_box + chr(10) + '</body>')
open('temp_report.html', 'w').write(html)
"

    # --- 2. Inject Hidden SVG at Bottom ---
    if [ -n "\$eh_img" ]; then
        echo "Injecting hidden SVG container..."
        python3 -c "
html_content = open('temp_report.html').read()
svg_content = open('\$eh_img').read()

# Injection: Hidden Div by default.
injection = f'''
<div id='fx_graph_container' style='display:none; margin-top: 50px; padding: 20px; border-top: 4px solid orange; background:#fff;'>
    <h2>Fragile X (FMR1) Graph <button onclick=\"this.parentElement.style.display='none'\" style='float:right; color:red;'>Close</button></h2>
    <div style='overflow-x: auto;'>
        {svg_content}
    </div>
</div>
</body>
'''
new_html = html_content.replace('</body>', injection)
with open('${sample_id}_visual_report.html', 'w') as f:
    f.write(new_html)
"
    else
        mv temp_report.html ${sample_id}_visual_report.html
    fi

    # Verify output
    if [ ! -f "${sample_id}_visual_report.html" ]; then
        echo "Error: Report generation failed."
        exit 1
    fi
    """
}
