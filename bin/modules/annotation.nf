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
//
// Static VEP flags live in modules/vep_annotation.flags. src/annotate_vcf.sh
// reads that same file. Do not paste a second copy into these process scripts.
// ============================================================

def vepAnnotationStaticFlags() {
    def flagFile = file("${projectDir}/modules/vep_annotation.flags")
    if (!flagFile.exists()) {
        throw new IllegalStateException("Missing shared VEP flags: ${flagFile}")
    }
    def tokens = flagFile.readLines()
        .collect { it.trim() }
        .findAll { it && !it.startsWith('#') }
    if (tokens.isEmpty()) {
        throw new IllegalStateException("Shared VEP flags file is empty: ${flagFile}")
    }
    return tokens.join(' ')
}

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
    // Static flags: modules/vep_annotation.flags (shared with src/annotate_vcf.sh).
    def vep_flags = vepAnnotationStaticFlags()
    def vep_forks = task.cpus
    """
    export TMPDIR=\$PWD

    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}_${params.variant_caller}_annotated.vcf.gz \\
        --stats_file ${sample_id}_${params.variant_caller}_vep_summary.html \\
        ${vep_flags} \\
        --dir_cache ${vep_cache_dir} \\
        --assembly GRCh38 \\
        --fasta ${ref_fasta} \\
        --fork ${vep_forks}

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
    def vep_flags = vepAnnotationStaticFlags()
    def vep_forks = task.cpus
    """
    export TMPDIR=\$PWD

    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}_${params.variant_caller}_annotated.vcf.gz \\
        --stats_file ${sample_id}_${params.variant_caller}_vep_summary.html \\
        ${vep_flags} \\
        --dir_cache ${vep_cache_dir} \\
        --assembly GRCh38 \\
        --fasta ${ref_fasta} \\
        --fork ${vep_forks}

    tabix -p vcf ${sample_id}_${params.variant_caller}_annotated.vcf.gz
    """
}
