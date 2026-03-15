nextflow.enable.dsl=2

// Set defaults to avoid warnings
params.hba_bed = null
params.cyp21a2_bed = null
params.cleanup = false // Default: Manual cleanup
params.output_dir = null // Output directory for portal files
params.sample_name = null // Sample name for tracking

include { INDEX_BWA; ALIGN_AND_SORT; MARK_DUPLICATES } from './modules/align'
include { GCNV_CLARITY; PREPROCESS_INTERVALS; COLLECT_READ_COUNTS; ANNOTATE_INTERVALS; GCNV_COHORT_RUN; GCNV_CASE_RUN; POSTPROCESS_GCNV } from './modules/cnv'
include { PARAPHASE_RUN; SMACA_RUN; PARAPHASE_RESCUE } from './modules/pseudogene'
include { GENERATE_SUMMARY_REPORT } from './modules/summary'
include { EXPANSION_HUNTER } from './modules/repeat'
include { MANTA_SV } from './modules/sv' 
include { DEPTH_ANALYSIS } from './modules/coverage' 
include { FALLBACK_ANALYSIS } from './modules/fallback' 
include { CALL_VARIANTS } from './modules/variant' 
include { GENERATE_VISUAL_EVIDENCE } from './modules/visualize'
include { PREPARE_VIZ_RESOURCES } from './modules/resources'


workflow {
    if (params.fastq_dir == null) {
        error "Please provide --fastq_dir"
    }
    
    // Channel for FASTQ pairs
    Channel
        .fromFilePairs("${params.fastq_dir}/*_{1,2,R1,R2}*.{fastq,fq}.gz")
        .set { fastq_ch }

    ref_fasta = file(params.ref_fasta)
    ref_fai = file(params.ref_fai)

    // Handling Indices for BWA
    if (params.ref_bwa_indices) {
        // Collect all BWA index files explicitly
        bwa_indices = Channel.fromPath("${params.ref_bwa_indices}/*.{amb,ann,bwt,pac,sa}").collect()
    } else {
        INDEX_BWA(ref_fasta)
        bwa_indices = INDEX_BWA.out.indices.map{ it[1] }.collect()
    }

    // 0. Prepare Shared Resources (Env, GTF, SMN Ref)
    PREPARE_VIZ_RESOURCES(ref_fasta, ref_fai)


    // 1. Alignment & Sorting
    ALIGN_AND_SORT(fastq_ch, ref_fasta, ref_fai, bwa_indices)
    
    // Mark Duplicates (Consuming Sorted BAM)
    MARK_DUPLICATES(ALIGN_AND_SORT.out.bam)
    
    // Use Marked BAM for downstream
    bam_ch = MARK_DUPLICATES.out.bam



    // 2. Track 1: Global CNV & SV
    // Preprocess Intervals (One-time, only if bed provided AND not supplied externally)
    interval_list_ch = Channel.empty()
    
    if (params.interval_list) {
        interval_list_ch = Channel.value(file(params.interval_list))
    } else if (params.backbone_bed) {
        PREPROCESS_INTERVALS(file(params.backbone_bed), ref_fasta, ref_fai, file(params.ref_dict))
        interval_list_ch = PREPROCESS_INTERVALS.out.interval_list
    }
    
    if (params.backbone_bed || params.interval_list) {
        // Collect Read Counts (Per sample)
        COLLECT_READ_COUNTS(bam_ch, interval_list_ch, ref_fasta, ref_fai, file(params.ref_dict))
    }

    if (!params.skip_cnv && (params.backbone_bed || params.interval_list)) { 
        if (params.pon_tar) {
             // CASE MODE (PoN provided)
             GCNV_CLARITY(COLLECT_READ_COUNTS.out.counts_hdf5, file(params.pon_tar), interval_list_ch)
        } else {
             // COHORT MODE (No PoN, learn from batch)
             ANNOTATE_INTERVALS(interval_list_ch, ref_fasta, ref_fai, file(params.ref_dict))
             
             // --- GCNV Mode Selection ---
        if (params.gcnv_model) {
            // CASE MODE (Recycle existing model)
            // Assuming params.gcnv_model is the PARENT directory containing "gcnv-model" and "ploidy-model".
            
            model_ch = Channel.value(file("${params.gcnv_model}/gcnv-model"))
            ploidy_model_ch = Channel.value(file("${params.gcnv_model}/ploidy-model"))
            
            GCNV_CASE_RUN(
                COLLECT_READ_COUNTS.out.counts_hdf5,
                ANNOTATE_INTERVALS.out.annotated_intervals,
                interval_list_ch, // Use interval_list_ch which is either from params.interval_list or PREPROCESS_INTERVALS
                model_ch,
                ploidy_model_ch
            )
            
            // Prepare channels for POSTPROCESS
            // Calls are per sample. Ploidy calls are per sample.
            // Join them by sample_id.
            postprocess_inputs = GCNV_CASE_RUN.out.calls
                .join(GCNV_CASE_RUN.out.ploidy_calls)
            
            // We pass the SAME model to postprocess
            postprocess_model = model_ch
            
            // Case mode handles configs internally (inside calls_path)
            postprocess_configs = Channel.value([])
            
        } else {
            // COHORT MODE (Train new model)
            // Collect all HDF5 files for cohort run
            all_counts_hdf5 = COLLECT_READ_COUNTS.out.counts_hdf5.map{it[1]}.collect()

            GCNV_COHORT_RUN(
                all_counts_hdf5,
                ANNOTATE_INTERVALS.out.annotated_intervals,
                interval_list_ch // Use interval_list_ch
            )
            
            // Prepare channels for POSTPROCESS
            // GCNV_COHORT_RUN emits "gcnv-calls/*". We need to separate configs (json) from sample dirs.
            
            GCNV_COHORT_RUN.out.calls.flatten()
                .branch {
                    configs: it.name.endsWith('.json')
                    dirs: it.isDirectory() && it.name != "gcnv-model" && it.name != "ploidy-model"
                    other: true
                }
                .set { cohort_output_parts }
            
            // Collect configs into a list to pass to every postprocess instance
            postprocess_configs = cohort_output_parts.configs.collect()
            
            // Flatten calls to [sample_id, calls_path]
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

            // Ploidy calls are one directory. We need to pass this SAME directory to every postprocess instance?
            // OR does POSTPROCESS accept the cohort ploidy directory?
            // Yes, "NOTE: The user must specify the ... contig-ploidy calls directory ... produced by DetermineGermlineContigPloidy."
            // So we pass the cohort ploidy calls directory.
            
            // We need to combine calls_ch with the single ploidy_calls object
            postprocess_inputs = calls_ch.combine(GCNV_COHORT_RUN.out.ploidy_calls)
            
            postprocess_model = GCNV_COHORT_RUN.out.model
        }

        POSTPROCESS_GCNV(
             postprocess_inputs, // [sample_id, calls_path, ploidy_calls_path]
             postprocess_model,
             ref_fasta, 
             ref_fai, 
             file(params.ref_dict),
             interval_list_ch, // Use interval_list_ch
             postprocess_configs
         )
        }
    }
    MANTA_SV(bam_ch, ref_fasta, ref_fai)

    // 2b. Track 1 Extension: Target Coverage & Intron Verification (Consolidated)
    if (params.backbone_bed) {
        dark_genes_plus = file("${projectDir}/../data/bed/dark_genes_plus.bed")
        DEPTH_ANALYSIS(bam_ch, file(params.backbone_bed), dark_genes_plus)
    }

    // 2c. Track 1 Extension: Fallback Analysis (HBA/CYP21A2)
    if (params.backbone_bed) {
         FALLBACK_ANALYSIS(bam_ch, 
             file(params.hba_bed ?: "${projectDir}/../data/bed/hba_targets.bed"), 
             file(params.cyp21a2_bed ?: "${projectDir}/../data/bed/cyp21a2_targets.bed"), 
             file(params.backbone_bed), 
             ref_fasta, 
             ref_fai
         )
    }

    // 3. Track 2: Pseudogene Resolution
    PARAPHASE_RUN(bam_ch, ref_fasta, ref_fai)
    SMACA_RUN(bam_ch)

    // 3b. Track 2b: Rescue Track (Dark Genes Plus)
    dark_genes_plus = file("${projectDir}/../data/bed/dark_genes_plus.bed")
    PARAPHASE_RESCUE(bam_ch, ref_fasta, dark_genes_plus)

    // VERIFY_INTRON_DEPTH Logic is now merged into DEPTH_ANALYSIS above studio

    // 4. Track 3: Repeat Expansion
    EXPANSION_HUNTER(bam_ch, ref_fasta, ref_fai, file(params.eh_catalog ?: params.ref_fasta)) 
    
    // 5. Track 4: Variant Calling (SNV/Indel) - Optional but recommended
    vcf_ch = Channel.empty()
    if (params.backbone_bed) {
         CALL_VARIANTS(bam_ch, ref_fasta, ref_fai, file(params.ref_dict), file(params.backbone_bed))
         vcf_ch = CALL_VARIANTS.out.vcf // Use filtered VCF
    }

    // 1b. Visual Evidence (IGV Snapshots + Fragile X Visuals)
    // Join BAM with VCF (optional) to show variant calls
    // Using remainder: true to keep samples even if VCF is missing
    
    viz_input = bam_ch.join(vcf_ch, remainder: true)
    
    // Pass ALL EH images; module will find the matching one
    eh_images_all = EXPANSION_HUNTER.out.images.map{it[1]}.collect()
    
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

    // Final Consolidation
    // Collect specific outputs (handling optionals safely if channels exist)
    
    // Safety: Define empty channels if processes didn't run
    manta_vcf_ch = MANTA_SV.out.vcf.map{it[1]}.collect()
    
    // GCNV VCFs might be empty if skip_cnv is true
    gcnv_vcf_ch = params.skip_cnv ? Channel.value([]) : POSTPROCESS_GCNV.out.vcf.map{it[1]}.collect()
    
    GENERATE_SUMMARY_REPORT(
        manta_vcf_ch,
        gcnv_vcf_ch,
        PARAPHASE_RESCUE.out.json.collect(), 
        FALLBACK_ANALYSIS.out.hba_report.collect(), 
        FALLBACK_ANALYSIS.out.cyp21a2_report.collect(), 
        EXPANSION_HUNTER.out.results.map{ it[2] }.collect(),
        SMACA_RUN.out.txt.map{ it[1] }.collect(),
        DEPTH_ANALYSIS.out.intron_report.collect(),
        GENERATE_VISUAL_EVIDENCE.out.snapshots.collect(), // Enforce dependency
        file(params.backbone_bed ?: "${projectDir}/../data/bed/Twist_Exome2.0_plus_Comprehensive_Exome_Spikein_targets_covered_annotated_hg38.bed")
    )
}

workflow.onComplete {
    def workDir = workflow.workDir
    def runName = workflow.runName
    
    println "\n" + "-"*60
    println "Pipeline Execution Summary"
    println "-"*60
    println "Run Name    : ${runName}"
    println "Status      : ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    println "Duration    : ${workflow.duration}"
    println "Work Dir    : ${workDir}"
    println "Analysis Dir: ${params.outdir}"
    if (params.output_dir) {
        println "Output Dir  : ${params.output_dir}"
    }
    println "-"*60
    
    if (workflow.success) {
        // Copy critical results to output directory
        if (params.output_dir) {
            println "Copying results to output directory..."
            
            def outputDir = file(params.output_dir)
            outputDir.mkdirs()
            
            // Create subdirectories
            def summaryOut = file("${params.output_dir}/summary")
            def snapshotsOut = file("${params.output_dir}/snapshots")
            def vcfOut = file("${params.output_dir}/vcf")
            def cnvOut = file("${params.output_dir}/cnv")
            def svOut = file("${params.output_dir}/sv")
            def repeatOut = file("${params.output_dir}/repeat")
            def pseudogeneOut = file("${params.output_dir}/pseudogene")
            def qcOut = file("${params.output_dir}/qc")
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
            
            // Copy files using shell commands
            def copyScript = """
            #!/bin/bash
            set -e
            
            # Summary reports
            if [ -d "${params.outdir}/summary" ]; then
                cp -r ${params.outdir}/summary/* ${summaryOut}/ 2>/dev/null || true
            fi
            
            # Snapshots
            if [ -d "${params.outdir}/snapshots" ]; then
                cp -r ${params.outdir}/snapshots/* ${snapshotsOut}/ 2>/dev/null || true
            fi
            
            # VCF files (filtered only)
            if [ -d "${params.outdir}/variant" ]; then
                cp ${params.outdir}/variant/*_filtered.vcf.gz* ${vcfOut}/ 2>/dev/null || true
            fi
            
            # CNV segments
            if [ -d "${params.outdir}/cnv/segments" ]; then
                cp ${params.outdir}/cnv/segments/*.vcf.gz* ${cnvOut}/ 2>/dev/null || true
            fi
            
            # SV calls
            if [ -d "${params.outdir}/sv" ]; then
                cp ${params.outdir}/sv/*.vcf.gz* ${svOut}/ 2>/dev/null || true
            fi
            
            # Repeat expansion results
            if [ -d "${params.outdir}/repeat" ]; then
                cp ${params.outdir}/repeat/*.vcf ${repeatOut}/ 2>/dev/null || true
                cp ${params.outdir}/repeat/*.json ${repeatOut}/ 2>/dev/null || true
                cp ${params.outdir}/repeat/*.svg ${repeatOut}/ 2>/dev/null || true
            fi
            
            # Pseudogene results
            if [ -d "${params.outdir}/pseudogene" ]; then
                cp -r ${params.outdir}/pseudogene/* ${pseudogeneOut}/ 2>/dev/null || true
            fi
            
            # QC metrics (coverage, duplicate metrics)
            if [ -d "${params.outdir}/coverage" ]; then
                cp ${params.outdir}/coverage/*_qc_metrics.txt ${qcOut}/ 2>/dev/null || true
                cp ${params.outdir}/coverage/*_target_coverage.txt ${qcOut}/ 2>/dev/null || true
            fi
            if [ -d "${params.outdir}/alignment" ]; then
                cp ${params.outdir}/alignment/*_duplicate_metrics.txt ${qcOut}/ 2>/dev/null || true
            fi
            
            # Pipeline info (trace, timeline, report)
            if [ -d "${params.outdir}/pipeline_info" ]; then
                cp ${params.outdir}/pipeline_info/trace.txt ${pipelineInfoOut}/ 2>/dev/null || true
                cp ${params.outdir}/pipeline_info/timeline.html ${pipelineInfoOut}/ 2>/dev/null || true
                cp ${params.outdir}/pipeline_info/report.html ${pipelineInfoOut}/ 2>/dev/null || true
            fi
            
            # Pipeline completion marker with metadata (for service-daemon)
            cat > ${params.output_dir}/pipeline_complete.json << MARKER
            {
                "status": "SUCCESS",
                "run_name": "${runName}",
                "duration": "${workflow.duration}",
                "sample_name": "${params.sample_name ?: 'unknown'}",
                "completed_at": "$(date -Iseconds)",
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
        println "-"*60 + "\n"
    }
}
