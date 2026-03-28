// ============================================================
// annotation.nf — Post-variant-calling Annotation module
//
// Runs Ensembl VEP (Variant Effect Predictor) on the filtered VCF
// produced by any of the three variant callers (GATK / DeepVariant / Strelka2).
//
// VEP adds the following information to the VCF INFO field (CSQ tag):
//   - Gene symbol, Ensembl Gene ID
//   - Transcript (MANE Select preferred)
//   - HGVSc, HGVSp (HGVS notation)
//   - Consequence (effect: missense_variant, stop_gained, etc.)
//   - gnomAD exomes/genomes allele frequency
//   - dbSNP rsID
//   - SIFT, PolyPhen-2 prediction scores
//   - ClinVar significance (for cross-reference; primary ClinVar logic stays in daemon)
//
// The output VCF (*_{gatk,deepvariant,strelka2}_annotated.vcf.gz) is the primary input to service-daemon.
// service-daemon parses the CSQ field and no longer needs to query
// gnomAD/dbSNP VCF files directly, reducing daemon memory usage significantly.
//
// VEP cache must be pre-installed on the server. See SERVER_SETUP.md.
// ============================================================

process VEP_ANNOTATION {
    tag "$sample_id"
    label 'vep'
    publishDir "${params.outdir}/variant", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path vep_cache_dir
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample_id),
          path("${sample_id}_${params.variant_caller}_annotated.vcf.gz"),
          path("${sample_id}_${params.variant_caller}_annotated.vcf.gz.tbi"),  emit: vcf
    path "${sample_id}_${params.variant_caller}_vep_summary.html",             emit: summary

    script:
    // VEP fork count: use all available CPUs
    def vep_forks = task.cpus
    """
    export TMPDIR=\$PWD

    # ── VEP annotation ──────────────────────────────────────
    # --offline        : use local cache (no internet required at runtime)
    # --cache          : enable cache mode
    # --dir_cache      : path to pre-downloaded VEP cache
    # --assembly       : GRCh38 (must match cache version)
    # --format vcf     : input is VCF
    # --vcf            : output in VCF format (adds CSQ INFO field)
    # --compress_output bgzip : bgzip the output for tabix indexing
    # --hgvs           : add HGVSc and HGVSp
    # --hgvsg          : add HGVSg (genomic HGVS)
    # --symbol         : add gene symbol (HGNC)
    # --mane           : annotate with MANE Select transcript
    # --canonical      : mark canonical transcripts
    # --af_gnomad      : add gnomAD exomes + genomes AF
    # --af_1kg         : add 1000 Genomes AF
    # --sift b         : add SIFT score and prediction
    # --polyphen b     : add PolyPhen-2 score and prediction
    # --check_existing : check against known variants (dbSNP rsID)
    # --variant_class  : add SO variant class
    # --numbers        : add exon/intron numbers
    # --gene_phenotype : add gene-phenotype associations
    # --pick            : pick one consequence per variant (MANE > canonical > most severe)
    # --pick_order     : priority for --pick
    # --fork           : parallel processing
    vep \
        --input_file ${vcf} \
        --output_file ${sample_id}_${params.variant_caller}_annotated.vcf.gz \
        --stats_file ${sample_id}_${params.variant_caller}_vep_summary.html \
        --format vcf \
        --vcf \
        --compress_output bgzip \
        --offline \
        --cache \
        --dir_cache ${vep_cache_dir} \
        --assembly GRCh38 \
        --fasta ${ref_fasta} \
        --hgvs \
        --hgvsg \
        --symbol \
        --mane \
        --canonical \
        --af_gnomad \
        --af_1kg \
        --sift b \
        --polyphen b \
        --check_existing \
        --variant_class \
        --numbers \
        --gene_phenotype \
        --pick \
        --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \
        --fork ${vep_forks} \
        --force_overwrite \
        --no_progress

    # ── Index the output VCF ────────────────────────────────
    tabix -p vcf ${sample_id}_${params.variant_caller}_annotated.vcf.gz
    """
}

// ============================================================
// VEP_ANNOTATION_DOCKER
//
// Alternative process using the official Ensembl VEP Docker image.
// Use this if VEP is not installed natively on the server.
// The VEP cache directory must still be pre-downloaded on the host
// and mounted into the container via Docker volume.
// ============================================================
process VEP_ANNOTATION_DOCKER {
    tag "$sample_id"
    label 'vep_docker'
    publishDir "${params.outdir}/variant", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path vep_cache_dir
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample_id),
          path("${sample_id}_${params.variant_caller}_annotated.vcf.gz"),
          path("${sample_id}_${params.variant_caller}_annotated.vcf.gz.tbi"),  emit: vcf
    path "${sample_id}_${params.variant_caller}_vep_summary.html",             emit: summary

    script:
    def vep_forks = task.cpus
    """
    export TMPDIR=\$PWD

    vep \
        --input_file ${vcf} \
        --output_file ${sample_id}_${params.variant_caller}_annotated.vcf.gz \
        --stats_file ${sample_id}_${params.variant_caller}_vep_summary.html \
        --format vcf \
        --vcf \
        --compress_output bgzip \
        --offline \
        --cache \
        --dir_cache ${vep_cache_dir} \
        --assembly GRCh38 \
        --fasta ${ref_fasta} \
        --hgvs \
        --hgvsg \
        --symbol \
        --mane \
        --canonical \
        --af_gnomad \
        --af_1kg \
        --sift b \
        --polyphen b \
        --check_existing \
        --variant_class \
        --numbers \
        --gene_phenotype \
        --pick \
        --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \
        --fork ${vep_forks} \
        --force_overwrite \
        --no_progress

    tabix -p vcf ${sample_id}_${params.variant_caller}_annotated.vcf.gz
    """
}
