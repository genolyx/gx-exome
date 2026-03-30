nextflow.enable.dsl=2

// Set defaults to avoid warnings
params.hba_bed = null
params.cyp21a2_bed = null
params.cleanup = false
params.output_dir = null
params.sample_name = null
params.skip_vep = false
params.vep_cache_dir = "${projectDir}/../data/refs/vep_cache"

// -------------------------------------------------------
// Module imports
// -------------------------------------------------------
// Alignment — both BWA-MEM and BWA-MEM2 processes
include { INDEX_BWA; ALIGN_AND_SORT;
          INDEX_BWA_MEM2; ALIGN_AND_SORT_BWA_MEM2;
          MARK_DUPLICATES; SAMTOOLS_BAM_STATS } from './modules/align'

// Variant Calling — three callers
include { CALL_VARIANTS_GATK;
          CALL_VARIANTS_DEEPVARIANT;
          CALL_VARIANTS_STRELKA2 } from './modules/variant'

// VEP Annotation — runs after variant calling
include { VEP_ANNOTATION } from './modules/annotation'

// CNV
include { GCNV_CLARITY; PREPROCESS_INTERVALS; COLLECT_READ_COUNTS;
          ANNOTATE_INTERVALS; GCNV_COHORT_RUN; GCNV_CASE_RUN;
          POSTPROCESS_GCNV } from './modules/cnv'

// Pseudogene
include { PARAPHASE_RUN; SMN_UNIFIED_C840_BAM; SMACA_RUN; PARAPHASE_RESCUE } from './modules/pseudogene'

// Other analyses
include { GENERATE_SUMMARY_REPORT } from './modules/summary'
include { EXPANSION_HUNTER }        from './modules/repeat'
include { MANTA_SV }                from './modules/sv'
include { DEPTH_ANALYSIS }          from './modules/coverage'
include { FALLBACK_ANALYSIS }       from './modules/fallback'
include { HBA_PARALOG_PILEUP }      from './modules/hba_paralog'
include { GENERATE_VISUAL_EVIDENCE } from './modules/visualize'
include { PREPARE_VIZ_RESOURCES }   from './modules/resources'

// -------------------------------------------------------
// Workflow
// -------------------------------------------------------
workflow {
    if (params.fastq_dir == null) {
        error "Please provide --fastq_dir"
    }

    // Validate aligner parameter
    def validAligners = ['bwa-mem', 'bwa-mem2']
    if (!validAligners.contains(params.aligner)) {
        error "Invalid aligner '${params.aligner}'. Choose one of: ${validAligners.join(', ')}"
    }

    // Validate variant_caller parameter
    def validCallers = ['gatk', 'deepvariant', 'strelka2']
    if (!validCallers.contains(params.variant_caller)) {
        error "Invalid variant_caller '${params.variant_caller}'. Choose one of: ${validCallers.join(', ')}"
    }

    println "=" * 60
    println "Carrier Screening Pipeline"
    println "  Aligner        : ${params.aligner}"
    println "  Variant Caller : ${params.variant_caller}"
    println "  VEP Annotation : ${params.skip_vep ? 'SKIPPED' : 'ENABLED'}"
    println "=" * 60

    // Channel for FASTQ pairs
    Channel
        .fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2}*.{fastq,fq}.gz")
        .set { fastq_ch }

    ref_fasta = file(params.ref_fasta)
    ref_fai   = file(params.ref_fai)

    // -------------------------------------------------------
    // 0. Prepare Shared Visualization Resources
    // -------------------------------------------------------
    PREPARE_VIZ_RESOURCES(ref_fasta, ref_fai)

    // -------------------------------------------------------
    // 1. Alignment — BWA-MEM or BWA-MEM2
    // -------------------------------------------------------
    if (params.aligner == 'bwa-mem2') {
        // BWA-MEM2: use pre-built index if available, otherwise build it
        if (params.ref_bwa_mem2_indices) {
            bwa_mem2_indices = Channel.fromPath(
                "${params.ref_bwa_mem2_indices}/*.{0123,amb,ann,bwt.2bit.64,pac}"
            ).collect()
        } else {
            INDEX_BWA_MEM2(ref_fasta)
            bwa_mem2_indices = INDEX_BWA_MEM2.out.indices.map { it[1] }.collect()
        }
        ALIGN_AND_SORT_BWA_MEM2(fastq_ch, ref_fasta, ref_fai, bwa_mem2_indices)
        raw_bam_ch = ALIGN_AND_SORT_BWA_MEM2.out.bam

    } else {
        // BWA-MEM (classic): use pre-built index if available, otherwise build it
        if (params.ref_bwa_indices) {
            bwa_indices = Channel.fromPath(
                "${params.ref_bwa_indices}/*.{amb,ann,bwt,pac,sa}"
            ).collect()
        } else {
            INDEX_BWA(ref_fasta)
            bwa_indices = INDEX_BWA.out.indices.map { it[1] }.collect()
        }
        ALIGN_AND_SORT(fastq_ch, ref_fasta, ref_fai, bwa_indices)
        raw_bam_ch = ALIGN_AND_SORT.out.bam
    }

    // MarkDuplicates — shared regardless of aligner
    MARK_DUPLICATES(raw_bam_ch)
    SAMTOOLS_BAM_STATS(MARK_DUPLICATES.out.bam)

    // Use MarkDuplicated BAM for all downstream processes
    bam_ch = MARK_DUPLICATES.out.bam

    // -------------------------------------------------------
    // 2. Track 1: Global CNV & SV
    // -------------------------------------------------------
    interval_list_ch = Channel.empty()

    if (!params.skip_cnv) {
        if (params.interval_list) {
            interval_list_ch = Channel.value(file(params.interval_list))
        } else if (params.backbone_bed) {
            PREPROCESS_INTERVALS(file(params.backbone_bed), ref_fasta, ref_fai, file(params.ref_dict))
            interval_list_ch = PREPROCESS_INTERVALS.out.interval_list
        }

        if (params.backbone_bed || params.interval_list) {
            COLLECT_READ_COUNTS(bam_ch, interval_list_ch, ref_fasta, ref_fai, file(params.ref_dict))
        }
    }

    if (!params.skip_cnv && (params.backbone_bed || params.interval_list)) {
        if (params.pon_tar) {
            GCNV_CLARITY(COLLECT_READ_COUNTS.out.counts_hdf5, file(params.pon_tar), interval_list_ch)
        } else {
            ANNOTATE_INTERVALS(interval_list_ch, ref_fasta, ref_fai, file(params.ref_dict))

            if (params.gcnv_model) {
                model_ch        = Channel.value(file("${params.gcnv_model}/gcnv-model"))
                ploidy_model_ch = Channel.value(file("${params.gcnv_model}/ploidy-model"))

                GCNV_CASE_RUN(
                    COLLECT_READ_COUNTS.out.counts_hdf5,
                    ANNOTATE_INTERVALS.out.annotated_intervals,
                    interval_list_ch,
                    model_ch,
                    ploidy_model_ch
                )

                postprocess_inputs  = GCNV_CASE_RUN.out.calls.join(GCNV_CASE_RUN.out.ploidy_calls)
                postprocess_model   = model_ch
                postprocess_configs = Channel.value([])

            } else {
                all_counts_hdf5 = COLLECT_READ_COUNTS.out.counts_hdf5.map { it[1] }.collect()
                GCNV_COHORT_RUN(
                    all_counts_hdf5,
                    ANNOTATE_INTERVALS.out.annotated_intervals,
                    interval_list_ch
                )

                GCNV_COHORT_RUN.out.calls.flatten()
                    .branch {
                        configs: it.name.endsWith('.json')
                        dirs:    it.isDirectory() && it.name != "gcnv-model" && it.name != "ploidy-model"
                        other:   true
                    }
                    .set { cohort_output_parts }

                postprocess_configs = cohort_output_parts.configs.collect()

                calls_ch = cohort_output_parts.dirs
                    .map { dir ->
                        def sample_name_file = dir.resolve("sample_name.txt")
                        if (sample_name_file.exists()) {
                            def sample_id = sample_name_file.text.trim()
                            return [sample_id, dir]
                        }
                        return null
                    }
                    .filter { it != null }

                postprocess_inputs = calls_ch.combine(GCNV_COHORT_RUN.out.ploidy_calls)
                postprocess_model  = GCNV_COHORT_RUN.out.model
            }

            POSTPROCESS_GCNV(
                postprocess_inputs,
                postprocess_model,
                ref_fasta,
                ref_fai,
                file(params.ref_dict),
                interval_list_ch,
                postprocess_configs
            )
        }
    }

    MANTA_SV(bam_ch, ref_fasta, ref_fai)

    // 2b. Target Coverage & Intron Verification
    if (params.backbone_bed) {
        dark_genes_plus = file(params.dark_genes_plus_bed)
        DEPTH_ANALYSIS(bam_ch, file(params.backbone_bed), dark_genes_plus)
    }

    // 2c. Fallback Analysis (HBA / CYP21A2)
    if (params.backbone_bed) {
        FALLBACK_ANALYSIS(
            bam_ch,
            file(params.hba_bed),
            file(params.cyp21a2_bed),
            file(params.backbone_bed),
            ref_fasta,
            ref_fai
        )
    }

    hba_paralog_py = Channel.fromPath("${projectDir}/modules/hba_paralog_pileup.py", checkIfExists: true)
    hba_paralog_in = bam_ch
        .combine(Channel.value(file(params.hba_paralog_sites, checkIfExists: true)))
        .combine(hba_paralog_py)
        .map { sid, bam, bai, sites, py -> tuple(sid, bam, bai, sites, py) }
    HBA_PARALOG_PILEUP(hba_paralog_in)

    // -------------------------------------------------------
    // 3. Track 2: Pseudogene Resolution
    // -------------------------------------------------------
    // Channel.fromPath ensures the helper script is staged into each task work dir (Docker smaca image).
    smaca_append_ch = Channel.fromPath("${projectDir}/modules/smaca_append_summary.py", checkIfExists: true)
    PARAPHASE_RUN(bam_ch, ref_fasta, ref_fai)
    // Re-align SMN1+SMN2 region reads to one SMN1 haplotype slice; pileup c.840 in smaca_append_summary.py
    smn_unified_in = bam_ch.combine(Channel.value(ref_fasta)).combine(Channel.value(ref_fai))
    SMN_UNIFIED_C840_BAM(smn_unified_in)
    smaca_joined = bam_ch.combine(smaca_append_ch).map { sid, bam, bai, py -> tuple(sid, bam, bai, py) }
        .join(SMN_UNIFIED_C840_BAM.out.unified.map { sid, ubam, ubai -> tuple(sid, ubam, ubai) }, by: 0)
    SMACA_RUN(smaca_joined)
    PARAPHASE_RESCUE(bam_ch, ref_fasta, file(params.dark_genes_plus_bed))

    // -------------------------------------------------------
    // 4. Track 3: Repeat Expansion
    // -------------------------------------------------------
    EXPANSION_HUNTER(bam_ch, ref_fasta, ref_fai, file(params.eh_catalog ?: params.ref_fasta))

    // -------------------------------------------------------
    // 5. Track 4: Variant Calling — selectable via params.variant_caller
    // -------------------------------------------------------
    vcf_ch = Channel.empty()

    if (params.backbone_bed) {
        if (params.variant_caller == 'deepvariant') {
            // DeepVariant: built-in CNN filter, no VariantFiltration needed
            CALL_VARIANTS_DEEPVARIANT(
                bam_ch,
                ref_fasta,
                ref_fai,
                file(params.backbone_bed)
            )
            vcf_ch = CALL_VARIANTS_DEEPVARIANT.out.vcf

        } else if (params.variant_caller == 'strelka2') {
            // Strelka2: built-in Random Forest filter, no VariantFiltration needed
            // --callRegions requires a bgzip+tabix .bed.gz file
            CALL_VARIANTS_STRELKA2(
                bam_ch,
                ref_fasta,
                ref_fai,
                file(params.backbone_bed_gz),
                file(params.backbone_bed_tbi)
            )
            vcf_ch = CALL_VARIANTS_STRELKA2.out.vcf

        } else {
            // GATK: HaplotypeCaller + VariantFiltration (hard filter)
            CALL_VARIANTS_GATK(
                bam_ch,
                ref_fasta,
                ref_fai,
                file(params.ref_dict),
                file(params.backbone_bed)
            )
            vcf_ch = CALL_VARIANTS_GATK.out.vcf
        }
    }

    // -------------------------------------------------------
    // 6. VEP Annotation — runs on the variant-called VCF
    //    Adds gene, transcript, HGVSc/p, gnomAD AF, dbSNP rsID,
    //    SIFT, PolyPhen-2 to the INFO/CSQ field.
    //    service-daemon parses CSQ directly; no longer needs to
    //    query gnomAD/dbSNP VCF files at runtime.
    //    Skip with --skip_vep true (service-daemon falls back to snpEff).
    // -------------------------------------------------------
    annotated_vcf_ch = Channel.empty()

    if (!params.skip_vep && params.backbone_bed) {
        vep_cache = file(params.vep_cache_dir)
        VEP_ANNOTATION(vcf_ch, vep_cache, ref_fasta, ref_fai)
        annotated_vcf_ch = VEP_ANNOTATION.out.vcf
    } else {
        // Pass through the raw VCF if VEP is skipped
        // service-daemon will use snpEff as fallback annotation
        annotated_vcf_ch = vcf_ch
    }

    // -------------------------------------------------------
    // 1b. Visual Evidence (IGV Snapshots)
    // -------------------------------------------------------
    viz_input    = bam_ch.join(annotated_vcf_ch, remainder: true)
    eh_images_all = EXPANSION_HUNTER.out.images.map { it[1] }.collect()

    GENERATE_VISUAL_EVIDENCE(
        viz_input,
        eh_images_all,
        ref_fasta,
        ref_fai,
        PREPARE_VIZ_RESOURCES.out.env,
        PREPARE_VIZ_RESOURCES.out.gtf,
        PREPARE_VIZ_RESOURCES.out.gtf_index,
        PREPARE_VIZ_RESOURCES.out.smn_ref,
        PREPARE_VIZ_RESOURCES.out.smn_index
    )

    // -------------------------------------------------------
    // Final Consolidation
    // -------------------------------------------------------
    manta_vcf_ch = MANTA_SV.out.vcf.map { it[1] }.collect()
    gcnv_vcf_ch  = params.skip_cnv ? Channel.value([]) : POSTPROCESS_GCNV.out.vcf.map { it[1] }.collect()

    GENERATE_SUMMARY_REPORT(
        manta_vcf_ch,
        gcnv_vcf_ch,
        PARAPHASE_RESCUE.out.json.collect(),
        FALLBACK_ANALYSIS.out.hba_report.collect(),
        FALLBACK_ANALYSIS.out.cyp21a2_report.collect(),
        EXPANSION_HUNTER.out.results.map { it[2] }.collect(),
        SMACA_RUN.out.txt.map { it[1] }.collect(),
        DEPTH_ANALYSIS.out.intron_report.collect(),
        GENERATE_VISUAL_EVIDENCE.out.snapshots.collect(),
        file(params.backbone_bed),
        HBA_PARALOG_PILEUP.out.tsv.collect()
    )
}

// -------------------------------------------------------
// Completion handler
// -------------------------------------------------------
workflow.onComplete {
    def workDir = workflow.workDir
    def runName = workflow.runName

    println "\n" + "-" * 60
    println "Pipeline Execution Summary"
    println "-" * 60
    println "Run Name       : ${runName}"
    println "Status         : ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    println "Duration       : ${workflow.duration}"
    println "Aligner        : ${params.aligner}"
    println "Variant Caller : ${params.variant_caller}"
    println "VEP Annotation : ${params.skip_vep ? 'SKIPPED' : 'ENABLED'}"
    println "Work Dir       : ${workDir}"
    println "Analysis Dir   : ${params.outdir}"
    if (params.output_dir) {
        println "Output Dir     : ${params.output_dir}"
    }
    println "-" * 60

    if (workflow.success) {
        if (params.output_dir) {
            println "Copying results to output directory..."

            def outputDir = file(params.output_dir)
            outputDir.mkdirs()

            def summaryOut     = file("${params.output_dir}/summary")
            def snapshotsOut   = file("${params.output_dir}/snapshots")
            def vcfOut         = file("${params.output_dir}/vcf")
            def cnvOut         = file("${params.output_dir}/cnv")
            def svOut          = file("${params.output_dir}/sv")
            def repeatOut      = file("${params.output_dir}/repeat")
            def pseudogeneOut  = file("${params.output_dir}/pseudogene")
            def qcOut          = file("${params.output_dir}/qc")
            def pipelineInfoOut = file("${params.output_dir}/pipeline_info")

            summaryOut.mkdirs()
            snapshotsOut.mkdirs()
            vcfOut.mkdirs()
            cnvOut.mkdirs()
            svOut.mkdirs()
            repeatOut.mkdirs()
            pseudogeneOut.mkdirs()
            qcOut.mkdirs()
            pipelineInfoOut.mkdirs()

            def copyScript = """
            #!/bin/bash
            set -e

            if [ -d "${params.outdir}/summary" ]; then
                cp -r ${params.outdir}/summary/* ${summaryOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/snapshots" ]; then
                cp -r ${params.outdir}/snapshots/* ${snapshotsOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/variant" ]; then
                # Copy annotated + filtered VCFs (separate cp: both can coexist after caller-specific naming)
                cp ${params.outdir}/variant/*_annotated.vcf.gz* ${vcfOut}/ 2>/dev/null || true
                cp ${params.outdir}/variant/*_filtered.vcf.gz* ${vcfOut}/ 2>/dev/null || true
                # VEP HTML summaries (e.g. *_deepvariant_vep_summary.html)
                cp ${params.outdir}/variant/*_vep_summary.html ${vcfOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/cnv/segments" ]; then
                cp ${params.outdir}/cnv/segments/*.vcf.gz* ${cnvOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/sv" ]; then
                cp ${params.outdir}/sv/*.vcf.gz* ${svOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/repeat" ]; then
                cp ${params.outdir}/repeat/*.vcf ${repeatOut}/ 2>/dev/null || true
                cp ${params.outdir}/repeat/*.json ${repeatOut}/ 2>/dev/null || true
                cp ${params.outdir}/repeat/*.svg ${repeatOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/pseudogene" ]; then
                cp -r ${params.outdir}/pseudogene/* ${pseudogeneOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/coverage" ]; then
                cp ${params.outdir}/coverage/*_qc_metrics.txt ${qcOut}/ 2>/dev/null || true
                cp ${params.outdir}/coverage/*_target_coverage.txt ${qcOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/alignment" ]; then
                cp ${params.outdir}/alignment/*_duplicate_metrics.txt ${qcOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/qc" ]; then
                cp ${params.outdir}/qc/*.stats.txt ${qcOut}/ 2>/dev/null || true
                cp ${params.outdir}/qc/*.bam.stats ${qcOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/pipeline_info" ]; then
                cp ${params.outdir}/pipeline_info/trace.txt ${pipelineInfoOut}/ 2>/dev/null || true
                cp ${params.outdir}/pipeline_info/timeline.html ${pipelineInfoOut}/ 2>/dev/null || true
                cp ${params.outdir}/pipeline_info/report.html ${pipelineInfoOut}/ 2>/dev/null || true
            fi

            # Pipeline completion marker (for service-daemon)
            cat > ${params.output_dir}/pipeline_complete.json << MARKER
            {
                "status": "SUCCESS",
                "run_name": "${runName}",
                "duration": "${workflow.duration}",
                "sample_name": "${params.sample_name ?: 'unknown'}",
                "aligner": "${params.aligner}",
                "variant_caller": "${params.variant_caller}",
                "vep_annotation": "${params.skip_vep ? 'skipped' : 'enabled'}",
                "completed_at": "\$(date -Iseconds)",
                "output_dir": "${params.output_dir}",
                "analysis_dir": "${params.outdir}"
            }
MARKER

            echo "Results copied to ${params.output_dir}"
            """

            def copyProc = ["bash", "-c", copyScript].execute()
            copyProc.waitFor()

            if (copyProc.exitValue() == 0) {
                println "✓ Results copied to output directory successfully"
            } else {
                println "⚠ Warning: Some files may not have been copied"
            }
        }

        if (params.cleanup) {
            println "Cleanup Flag Detected. Removing intermediate files..."
            println "Deleting Work Directory: ${workDir}"
            ["rm", "-rf", "${workDir}"].execute().waitFor()
            println "Automatic Cleanup Complete."
        } else {
            println "TIP: Intermediate files kept in analysis directory"
            println "     To cleanup: rm -rf ${params.outdir}"
        }
        println "-" * 60 + "\n"
    }
}
