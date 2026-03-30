#!/usr/bin/env python3
"""Synthetic BAM: one read per (chrom, pos) covering ref, for hba_paralog_pileup demo."""
import csv
import sys
from pathlib import Path

import pysam

# demo reads: 81bp, pos at center (40 offset)
READ_LEN = 81
HALF = READ_LEN // 2


def main() -> None:
    ref = Path(__file__).parent / "chr16.fa"
    sites = Path(__file__).parent / "hba_sites.tsv"
    out = Path(__file__).parent / "demo_synthetic.bam"
    fa = pysam.FastaFile(str(ref))
    chrom = "chr16"
    ln = fa.get_reference_length(chrom)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": ln}],
    }

    positions: set[int] = set()
    with open(sites, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            positions.add(int(row["hba1_pos"]))
            positions.add(int(row["hba2_pos"]))

    out_bam = pysam.AlignmentFile(str(out), "wb", header=header)
    q = 0
    for pos in sorted(positions):
        s0 = max(0, pos - 1 - HALF)
        e0 = min(ln, pos - 1 + HALF + 1)
        seq = fa.fetch(chrom, s0, e0)
        if len(seq) < READ_LEN:
            continue
        # center read on pos (1-based)
        start0 = pos - 1 - HALF
        if start0 < 0:
            start0 = 0
        seq = fa.fetch(chrom, start0, start0 + READ_LEN)
        if len(seq) != READ_LEN:
            continue
        aln = pysam.AlignedSegment(header=out_bam.header)
        aln.query_name = f"demo_{pos}"
        aln.query_sequence = seq
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = start0
        aln.mapping_quality = 60
        aln.cigar = [(0, READ_LEN)]  # M
        aln.query_qualities = pysam.qualitystring_to_array("~" * READ_LEN)
        out_bam.write(aln)
        q += 1
    out_bam.close()
    fa.close()
    pysam.index(str(out))
    print(f"Wrote {q} reads to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
