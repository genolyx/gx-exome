#!/usr/bin/env python3
"""
Extract homologous single-base differences between HBA1 and HBA2 on GRCh38.

Workflow:
  1) Fetch gene-span sequences from a reference FASTA (default: Ensembl GRCh38 gene bounds).
  2) Global pairwise alignment → list mismatch columns (homologous sites with different ref bases).
  3) Optional: count k-mer occurrences in the genome (pysam) to flag non-unique windows.
  4) Optional: print Ensembl overlap URLs for manual dbSNP/gnomAD follow-up (no API calls by default).

Defaults (GRCh38, Ensembl gene coordinates, inclusive 1-based):
  HBA2: chr16:172876-173710
  HBA1: chr16:176680-177522

Does NOT substitute for clinical validation — freeze coordinates for your assembly build.

Next step: pileup each site in a BAM (SMN-style A/C/G/T counts):
  python3 scripts/hba_paralog_pileup.py --bam sample.bam --sites hba_paralog_mismatches.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# --- Default GRCh38 gene spans (Ensembl; includes introns within gene bounds) ---
DEFAULT_HBA2 = ("chr16", 172876, 173710)
DEFAULT_HBA1 = ("chr16", 176680, 177522)


def parse_region(s: str) -> tuple[str, int, int]:
    """Parse chr16:172876-173710 or 16:172876-173710 (1-based inclusive)."""
    m = re.match(r"^(chr)?(\d+):(\d+)-(\d+)$", s.strip(), re.I)
    if not m:
        raise ValueError(f"Bad region (expected chrN:start-end): {s!r}")
    chrom = f"chr{m.group(2)}"
    start, end = int(m.group(3)), int(m.group(4))
    if end < start:
        raise ValueError("end < start")
    return chrom, start, end


def fetch_fasta_region(fasta: str, chrom: str, start: int, end: int) -> str:
    """1-based inclusive coordinates; returns upper-case sequence."""
    import pysam

    fa = pysam.FastaFile(fasta)
    try:
        # pysam 0-based half-open
        seq = fa.fetch(chrom, start - 1, end)
    finally:
        fa.close()
    return seq.upper()


def align_pairwise(seq1: str, seq2: str):
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2
    alns = aligner.align(seq1, seq2)
    return alns[0]


def mismatches_from_alignment(
    alignment,
    chrom: str,
    start1: int,
    start2: int,
) -> list[tuple[int, str, int, str]]:
    """
    Map aligned columns to genomic coordinates.
    start1/start2: 1-based start of each extracted interval on chrom.
    """
    a = alignment[0]
    b = alignment[1]
    g1 = start1
    g2 = start2
    out: list[tuple[int, str, int, str]] = []
    for i in range(len(a)):
        ca, cb = a[i], b[i]
        if ca == "-":
            g2 += 1
        elif cb == "-":
            g1 += 1
        else:
            if ca != cb:
                out.append((g1, ca, g2, cb))
            g1 += 1
            g2 += 1
    return out


def revcomp(s: str) -> str:
    t = str.maketrans("ACGT", "TGCA")
    return s.translate(t)[::-1]


def count_kmer_occurrences(
    fasta: str, kmer: str, chrom: str | None = None
) -> int:
    """
    Count forward + reverse-complement occurrences.
    If chrom is set, only that chromosome (fast for WES-relevance on chr16).
    Else scan primary assembly chromosomes (slow on first use per k-mer).
    """
    import pysam

    kmer_u = kmer.upper()
    rc = revcomp(kmer_u)
    total = 0
    fa = pysam.FastaFile(fasta)
    try:
        if chrom:
            seq = fa.fetch(chrom).upper()
            total += seq.count(kmer_u)
            if rc != kmer_u:
                total += seq.count(rc)
        else:
            for ref in fa.references:
                if ref.startswith("chr") and ref not in (
                    "chrM",
                    "chrMT",
                ) and "_" not in ref:
                    seq = fa.fetch(ref).upper()
                    total += seq.count(kmer_u)
                    if rc != kmer_u:
                        total += seq.count(rc)
    finally:
        fa.close()
    return total


def kmer_at(fasta: str, chrom: str, pos: int, k: int) -> str:
    """Extract k-mer centered on pos (1-based); pad if near edge."""
    import pysam

    half = k // 2
    s = pos - half
    e = pos + half
    fa = pysam.FastaFile(fasta)
    try:
        seq = fa.fetch(chrom, s - 1, e).upper()
    finally:
        fa.close()
    return seq


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--fasta",
        required=True,
        help="Path to GRCh38 reference FASTA (indexed: .fai)",
    )
    ap.add_argument(
        "--hba2-region",
        default="chr16:172876-173710",
        help="HBA2 interval (default: Ensembl GRCh38 gene span)",
    )
    ap.add_argument(
        "--hba1-region",
        default="chr16:176680-177522",
        help="HBA1 interval (default: Ensembl GRCh38 gene span)",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write TSV path (default: stdout)",
    )
    ap.add_argument(
        "--kmer",
        type=int,
        metavar="K",
        help="If set, extract K-mer centered on each HBA1 position and count occurrences; adds columns",
    )
    ap.add_argument(
        "--kmer-full-genome",
        action="store_true",
        help="With --kmer, scan whole primary assembly (slow); default is chr16 only",
    )
    ap.add_argument(
        "--print-ensembl-links",
        action="store_true",
        help="Add column with Ensembl region overlap URL for manual variant DB lookup",
    )
    args = ap.parse_args()

    c2, s2, e2 = parse_region(args.hba2_region)
    c1, s1, e1 = parse_region(args.hba1_region)
    if c1 != c2:
        print("HBA1 and HBA2 must be on the same chromosome.", file=sys.stderr)
        sys.exit(1)

    seq_hba2 = fetch_fasta_region(args.fasta, c2, s2, e2)
    seq_hba1 = fetch_fasta_region(args.fasta, c1, s1, e1)

    aln = align_pairwise(seq_hba1, seq_hba2)
    mm = mismatches_from_alignment(aln, c1, s1, s2)

    fieldnames = [
        "hba1_chr",
        "hba1_pos",
        "hba1_ref",
        "hba2_chr",
        "hba2_pos",
        "hba2_ref",
    ]
    if args.kmer:
        fieldnames.extend(
            [
                f"kmer_{args.kmer}",
                "kmer_count_scope",
                f"kmer_count_rc",
            ]
        )
    if args.print_ensembl_links:
        fieldnames.append("ensembl_overlap_url")

    kmer_cache: dict[str, int] = {}
    rows = []
    for g1, b1, g2, b2 in mm:
        row = {
            "hba1_chr": c1,
            "hba1_pos": g1,
            "hba1_ref": b1,
            "hba2_chr": c2,
            "hba2_pos": g2,
            "hba2_ref": b2,
        }
        if args.kmer:
            km = kmer_at(args.fasta, c1, g1, args.kmer)
            row[f"kmer_{args.kmer}"] = km
            if len(km) != args.kmer:
                row["kmer_count_scope"] = ""
                row["kmer_count_rc"] = ""
            else:
                scope = "genome" if args.kmer_full_genome else c1
                row["kmer_count_scope"] = scope
                if km not in kmer_cache:
                    kmer_cache[km] = count_kmer_occurrences(
                        args.fasta,
                        km,
                        chrom=None if args.kmer_full_genome else c1,
                    )
                row["kmer_count_rc"] = kmer_cache[km]
        if args.print_ensembl_links:
            # 1bp region for overlap lookup (manual gnomAD/dbSNP follow-up)
            row["ensembl_overlap_url"] = (
                f"https://www.ensembl.org/Homo_sapiens/Location/View?r={c1}:{g1}-{g1}"
            )
        rows.append(row)

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    finally:
        if args.output:
            out.close()

    print(
        f"# mismatches: {len(mm)} (gene-span alignment; verify against your assembly)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
