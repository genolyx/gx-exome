// ============================================================
// scratch.nf — SSD Scratch lifecycle management
//
// Strategy overview:
//   Alignment (BWA + samtools sort + MarkDuplicates) is the most
//   I/O-intensive phase of the exome pipeline. By directing the
//   Nextflow workDir to a local SSD the pipeline achieves:
//     HDD 180 MB/s → SSD 3,000 MB/s (up to ~16× faster BAM I/O)
//
//   Flow:
//     1. SCRATCH_SETUP     — validate SSD dir, check free space,
//                            clean up any leftovers from prior runs
//     2. (alignment runs here with workDir = SSD)
//     3. SCRATCH_MOVE_FINAL — relay final .md.bam downstream;
//                             log SSD usage post-alignment
//     4. Remaining pipeline continues on SSD workDir; publishDir
//        copies final results to HDD analysis/output dirs
//     5. On success: cleanup=true in ssd_scratch profile removes
//        the SSD workDir automatically
//     6. On failure: onComplete handler + run_analysis.sh exit
//        trap remove the SSD scratch dir
//
// Activation: --use-ssd flag in run_analysis.sh (adds
//   -profile docker,ssd_scratch and redirects HOST_WORK_DIR)
// ============================================================

// -------------------------------------------------------
// SCRATCH_SETUP — validate SSD directory before alignment
// Emits scratch_dir val as a barrier channel: alignment
// processes only start after this process succeeds.
// -------------------------------------------------------
process SCRATCH_SETUP {
    tag "${sample_name}"
    label 'scratch'

    input:
    val sample_name
    val scratch_dir

    output:
    val scratch_dir, emit: scratch_dir

    script:
    """
    set -euo pipefail

    mkdir -p "${scratch_dir}"

    # Remove leftover BAMs/BAIs from any previous failed run
    find "${scratch_dir}" -maxdepth 3 \\( -name '*.bam' -o -name '*.bai' \\) \\
        -delete 2>/dev/null || true

    # Check available space and warn if below safe threshold
    AVAIL_KB=\$(df -Pk "${scratch_dir}" | awk 'NR==2{print \$4}')
    AVAIL_GB=\$(( AVAIL_KB / 1048576 ))
    TOTAL_KB=\$(df -Pk "${scratch_dir}" | awk 'NR==2{print \$2}')
    TOTAL_GB=\$(( TOTAL_KB / 1048576 ))

    echo "[SSD] Initialized: ${scratch_dir}"
    echo "[SSD] Available: \${AVAIL_GB} GB / \${TOTAL_GB} GB total"

    if [ "\${AVAIL_GB}" -lt 150 ]; then
        echo "[SSD] WARNING: < 150 GB free — BWA-MEM2 sort + MarkDuplicates may need 80–120 GB per sample"
    fi

    # Confirm this mount is actually faster than HDD (ROTA=0 means SSD/NVMe)
    MOUNT_DEV=\$(df "${scratch_dir}" | awk 'NR==2{print \$1}' | sed 's|[0-9]*\$||; s|/dev/||')
    ROTA=\$(cat /sys/block/\${MOUNT_DEV}/queue/rotational 2>/dev/null || echo "unknown")
    if [ "\${ROTA}" = "1" ]; then
        echo "[SSD] WARNING: device \${MOUNT_DEV} appears to be a spinning HDD (ROTA=1) — SSD benefit may be limited"
    elif [ "\${ROTA}" = "0" ]; then
        echo "[SSD] Device \${MOUNT_DEV}: SSD/NVMe confirmed (ROTA=0)"
    else
        echo "[SSD] Device type unknown (rotational check not available in container)"
    fi
    """
}

// -------------------------------------------------------
// SCRATCH_MOVE_FINAL — relay final .md.bam through the
// SSD workDir, log post-alignment SSD usage, and hand the
// BAM to downstream processes via Nextflow channels.
//
// publishDir 'copy' in MARK_DUPLICATES already wrote the
// HDD copy under outdir/alignment/.  This process does NOT
// re-copy; it creates symlinks to ensure the output tuple
// is tracked correctly when workDir is on SSD.
// -------------------------------------------------------
process SCRATCH_MOVE_FINAL {
    tag "${sample_id}"
    label 'scratch'

    input:
    tuple val(sample_id), path(final_bam), path(final_bai)
    val  scratch_dir

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: bam

    script:
    """
    set -euo pipefail

    # Staged inputs already carry the correct filenames; verify they exist
    test -f "${final_bam}" || { echo "[SSD] ERROR: ${final_bam} not found"; exit 1; }
    test -f "${final_bai}" || { echo "[SSD] ERROR: ${final_bai} not found"; exit 1; }

    # Ensure output names match (staged names should already match, but be explicit)
    if [ ! -f "${sample_id}.md.bam" ]; then
        ln -sf "\$(realpath ${final_bam})" "${sample_id}.md.bam"
    fi
    if [ ! -f "${sample_id}.md.bai" ]; then
        ln -sf "\$(realpath ${final_bai})" "${sample_id}.md.bai"
    fi

    # Log SSD usage after the alignment steps (BWA + sort + MarkDup)
    AVAIL_KB=\$(df -Pk "${scratch_dir}" 2>/dev/null | awk 'NR==2{print \$4}' || echo 0)
    AVAIL_GB=\$(( AVAIL_KB / 1048576 ))
    USED=\$(du -sh "${scratch_dir}" 2>/dev/null | cut -f1 || echo "unknown")
    echo "[SSD] Post-alignment scratch usage: \${USED}  (free: \${AVAIL_GB} GB)"
    echo "[SSD] Final BAM handed to downstream: ${sample_id}.md.bam"
    """
}

// -------------------------------------------------------
// SCRATCH_CLEANUP_ON_FAILURE — remove the SSD scratch tree
// after a pipeline failure to reclaim SSD space.
// errorStrategy 'ignore' prevents a cleanup failure from
// masking the original pipeline error.
//
// Note: this process is invoked programmatically from the
// workflow.onComplete handler in main.nf (not via a
// standard process call in the workflow {} block).
// -------------------------------------------------------
process SCRATCH_CLEANUP_ON_FAILURE {
    tag "${sample_name}"
    label 'scratch'
    errorStrategy 'ignore'

    input:
    val sample_name
    val scratch_dir

    output:
    stdout

    script:
    """
    echo "[SSD] Cleaning up scratch after failure: ${scratch_dir}"
    rm -rf "${scratch_dir}" 2>/dev/null || true
    echo "[SSD] Done"
    """
}
