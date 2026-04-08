#!/usr/bin/env python3
"""
Pileup A/C/G/T at CYP21A2 vs CYP21A1P paralog mismatch sites (from extract_cyp21_paralog_sites.py).

Same reference-forward logic as hba_paralog_pileup.py. CYP21A2 on chr6 is forward (+) in GRCh38.

Example:
  python3 extract_cyp21_paralog_sites.py --fasta GRCh38.fa -o sites.tsv
  python3 cyp21_paralog_pileup.py --bam sample.bam --sites sites.tsv

Summarizes across all sites:
  ratio_cyp21a2_ref_mass = sum(reads matching CYP21A2 ref) / (sum CYP21A2 ref + sum CYP21A1P ref reads).
  est_copies_cyp21a2 = total_diploid_copies * ratio (default total_diploid_copies=2).
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

import pysam


def pick_chr6(bam: pysam.AlignmentFile) -> str | None:
    for c in ("chr6", "6"):
        if c in bam.references:
            return c
    return None


def normalize_contig(refs: set[str], chrom: str) -> str | None:
    """Map chr6/6 to whichever exists in the BAM header."""
    if chrom in refs:
        return chrom
    if chrom == "chr6" and "6" in refs:
        return "6"
    if chrom == "6" and "chr6" in refs:
        return "chr6"
    return None


def complement_dna(b: str) -> str:
    m = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return m.get(b.upper(), "N")


def read_base_on_reference_forward(aln: pysam.AlignedSegment, query_pos: int) -> str:
    b = aln.query_sequence[query_pos].upper()
    if aln.is_reverse:
        return complement_dna(b)
    return b


def base_in_hgvs_coding_space(ref_forward_base: str, minus_strand_gene: bool) -> str:
    if not minus_strand_gene:
        return ref_forward_base
    return complement_dna(ref_forward_base)


def pileup_acgt_at_1bp(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos_1based: int,
    *,
    hgvs_coding_minus_strand: bool = False,
):
    """Count A/C/G/T/N at one reference position (1-based). Same semantics as SMN pileup."""
    start0 = pos_1based - 1
    end0 = pos_1based
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    try:
        for col in bam.pileup(
            contig=chrom,
            start=start0,
            stop=end0,
            truncate=True,
            min_mapping_quality=0,
            min_base_quality=0,
            max_depth=1_000_000,
        ):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                if pr.query_position is None:
                    continue
                aln = pr.alignment
                qp = pr.query_position
                rf = read_base_on_reference_forward(aln, qp)
                b = base_in_hgvs_coding_space(rf, hgvs_coding_minus_strand)
                if b in counts:
                    counts[b] += 1
                else:
                    counts["N"] += 1
    except (ValueError, OSError) as e:
        print(f"WARNING: pileup failed at {chrom}:{pos_1based}: {e}", file=sys.stderr)
        return None
    return counts


def depth(ct: dict | None) -> int:
    if not ct:
        return 0
    return int(sum(ct.values()))


def frac_ref(ct: dict | None, ref_base: str) -> str:
    f = frac_ref_float(ct, ref_base)
    if f is None:
        return "NA"
    return f"{f:.4f}"


def frac_ref_float(ct: dict | None, ref_base: str) -> float | None:
    if not ct:
        return None
    r = ref_base.upper()
    if r not in ct:
        return None
    d = sum(ct.values())
    if d == 0:
        return None
    return ct[r] / d


def fmt_acgt(ct: dict | None) -> str:
    if not ct:
        return "A=NA C=NA G=NA T=NA depth=NA"
    d = depth(ct)
    return (
        f"A={ct['A']} C={ct['C']} G={ct['G']} T={ct['T']} "
        f"N={ct['N']} depth={d}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True, help="Indexed WGS/WES BAM")
    ap.add_argument(
        "--sites",
        required=True,
        type=Path,
        help="TSV from extract_cyp21_paralog_sites.py (tab-separated, header row)",
    )
    ap.add_argument(
        "--hgvs-minus-strand",
        action="store_true",
        help="Tally bases in HGVS coding sense for minus-strand genes (CYP21A2 is + strand; leave off)",
    )
    ap.add_argument(
        "--max-sites",
        type=int,
        default=0,
        help="If >0, only process first N data rows (debug)",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write TSV (default: stdout)",
    )
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print ##CYP21_PARALOG_PILEUP lines to stderr (CN summary still printed)",
    )
    ap.add_argument(
        "--no-cn-summary",
        action="store_true",
        help="Do not print ##HBA_PARALOG_CN_RATIO line or append summary to -o file",
    )
    ap.add_argument(
        "--total-diploid-copies",
        type=int,
        default=2,
        help="Scale factor for est_copies from ref-mass ratio (default 2 = diploid CYP21A2 alleles; heuristic vs pseudogene arm)",
    )
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_refs = set(bam.references)
    default_chrom = pick_chr6(bam)
    if not default_chrom:
        bam.close()
        print("BAM has no chr6 or 6 contig.", file=sys.stderr)
        sys.exit(1)

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    cn_summary_for_file: str | None = None
    try:
        fieldnames = [
            "idx",
            "chrom",
            "cyp21a2_pos",
            "cyp21a2_ref",
            "cyp21a2_ref_frac",
            "cyp21a2_pileup",
            "cyp21a1p_pos",
            "cyp21a1p_ref",
            "cyp21a1p_ref_frac",
            "cyp21a1p_pileup",
        ]
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()

        sum_cyp21a2_ref_reads = 0
        sum_cyp21a1p_ref_reads = 0
        sum_cyp21a2_depth = 0
        sum_cyp21a1p_depth = 0
        n_rows = 0
        site_frac_cyp21a2: list[float] = []
        site_frac_cyp21a1p: list[float] = []

        with open(args.sites, newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for idx, row in enumerate(rdr, start=1):
                if args.max_sites and idx > args.max_sites:
                    break
                try:
                    p1 = int(row["cyp21a2_pos"])
                    p2 = int(row["cyp21a1p_pos"])
                    r1 = row["cyp21a2_ref"].strip().upper()
                    r2 = row["cyp21a1p_ref"].strip().upper()
                except (KeyError, ValueError) as e:
                    print(f"Skip row {idx}: {e}", file=sys.stderr)
                    continue

                chrom = (row.get("cyp21a2_chr") or row.get("cyp21a1p_chr") or "").strip()
                if not chrom:
                    chrom = default_chrom
                chrom = normalize_contig(bam_refs, chrom) or default_chrom

                c1 = pileup_acgt_at_1bp(
                    bam,
                    chrom,
                    p1,
                    hgvs_coding_minus_strand=args.hgvs_minus_strand,
                )
                c2 = pileup_acgt_at_1bp(
                    bam,
                    chrom,
                    p2,
                    hgvs_coding_minus_strand=args.hgvs_minus_strand,
                )

                if c1 and r1 in c1:
                    sum_cyp21a2_ref_reads += int(c1[r1])
                if c2 and r2 in c2:
                    sum_cyp21a1p_ref_reads += int(c2[r2])
                sum_cyp21a2_depth += depth(c1)
                sum_cyp21a1p_depth += depth(c2)
                n_rows += 1
                f1 = frac_ref_float(c1, r1)
                f2 = frac_ref_float(c2, r2)
                if f1 is not None:
                    site_frac_cyp21a2.append(f1)
                if f2 is not None:
                    site_frac_cyp21a1p.append(f2)

                w.writerow(
                    {
                        "idx": idx,
                        "chrom": chrom,
                        "cyp21a2_pos": p1,
                        "cyp21a2_ref": r1,
                        "cyp21a2_ref_frac": frac_ref(c1, r1),
                        "cyp21a2_pileup": fmt_acgt(c1),
                        "cyp21a1p_pos": p2,
                        "cyp21a1p_ref": r2,
                        "cyp21a1p_ref_frac": frac_ref(c2, r2),
                        "cyp21a1p_pileup": fmt_acgt(c2),
                    }
                )

                if not args.quiet:
                    print(
                        f"##CYP21_PARALOG_PILEUP idx={idx} chrom={chrom} "
                        f"cyp21a2_pos={p1} expected_ref={r1} | cyp21a1p_pos={p2} expected_ref={r2}",
                        file=sys.stderr,
                    )
                    print(
                        f"  CYP21A2_ref {fmt_acgt(c1)} | CYP21A1P_ref {fmt_acgt(c2)}",
                        file=sys.stderr,
                    )

        if not args.no_cn_summary:
            tot_ref = sum_cyp21a2_ref_reads + sum_cyp21a1p_ref_reads
            if tot_ref > 0:
                ratio_a2 = sum_cyp21a2_ref_reads / tot_ref
                tc = args.total_diploid_copies
                est_a2 = (tc * ratio_a2) if tc else float("nan")
                est_a1p = tc - est_a2 if tc else float("nan")
                ratio_s = f"{ratio_a2:.6f}"
                est1_s = f"{est_a2:.4f}"
                est2_s = f"{est_a1p:.4f}"
            else:
                ratio_a2 = None
                ratio_s = "NA"
                est1_s = "NA"
                est2_s = "NA"

            summary_body = (
                "CYP21_PARALOG_CN_RATIO "
                f"cyp21a2_ref_reads_total={sum_cyp21a2_ref_reads} "
                f"cyp21a1p_ref_reads_total={sum_cyp21a1p_ref_reads} "
                f"ratio_cyp21a2_ref_mass={ratio_s} "
                f"(CYP21A2_ref/(CYP21A2_ref+CYP21A1P_ref) over {n_rows} sites; "
                f"est_copies_cyp21a2={est1_s} est_copies_cyp21a1p_arm={est2_s} "
                f"scaled_from_total_diploid_copies={args.total_diploid_copies}) "
                f"sum_depth_cyp21a2={sum_cyp21a2_depth} sum_depth_cyp21a1p={sum_cyp21a1p_depth}"
            )

            def _fmt_stats(xs: list[float]) -> str:
                if not xs:
                    return "min=NA mean=NA max=NA n=0"
                return (
                    f"min={min(xs):.4f} mean={statistics.mean(xs):.4f} "
                    f"max={max(xs):.4f} n={len(xs)}"
                )

            tot_depth = sum_cyp21a2_depth + sum_cyp21a1p_depth
            depth_ratio_a2 = (
                (sum_cyp21a2_depth / tot_depth) if tot_depth > 0 else float("nan")
            )
            pooled1 = (
                (sum_cyp21a2_ref_reads / sum_cyp21a2_depth)
                if sum_cyp21a2_depth > 0
                else None
            )
            pooled2 = (
                (sum_cyp21a1p_ref_reads / sum_cyp21a1p_depth)
                if sum_cyp21a1p_depth > 0
                else None
            )
            frac_a1p_mass = (
                (sum_cyp21a1p_ref_reads / tot_ref) if tot_ref > 0 else float("nan")
            )
            p1s = f"{pooled1:.6f}" if pooled1 is not None else "NA"
            p2s = f"{pooled2:.6f}" if pooled2 is not None else "NA"
            drs = f"{depth_ratio_a2:.6f}" if tot_depth > 0 else "NA"
            f2ms = f"{frac_a1p_mass:.6f}" if tot_ref > 0 else "NA"

            r1_over_tot = (
                f"{sum_cyp21a2_ref_reads}/{tot_ref}" if tot_ref > 0 else "NA/NA"
            )
            r2_over_tot = (
                f"{sum_cyp21a1p_ref_reads}/{tot_ref}" if tot_ref > 0 else "NA/NA"
            )
            trace_body = (
                "CYP21_PARALOG_CN_TRACE "
                f"tot_ref_reads_denominator={tot_ref} "
                f"(cyp21a2_ref_reads_total+cyp21a1p_ref_reads_total) "
                f"frac_cyp21a2_ref_mass={r1_over_tot}={ratio_s} "
                f"frac_cyp21a1p_ref_mass={r2_over_tot}={f2ms} "
                f"formula_est_cyp21a2={args.total_diploid_copies}*"
                f"{sum_cyp21a2_ref_reads}/({sum_cyp21a2_ref_reads}+{sum_cyp21a1p_ref_reads})={est1_s} "
                f"formula_est_cyp21a1p_arm={args.total_diploid_copies}*"
                f"{sum_cyp21a1p_ref_reads}/({sum_cyp21a2_ref_reads}+{sum_cyp21a1p_ref_reads})={est2_s} "
                f"sum_depth_total={tot_depth} "
                f"depth_mass_ratio_cyp21a2={drs} "
                f"(sum_depth_cyp21a2/sum_depth_total; compare to ratio_cyp21a2_ref_mass={ratio_s}) "
                f"pooled_ref_reads_over_depth_cyp21a2={p1s} "
                f"pooled_ref_reads_over_depth_cyp21a1p={p2s} "
                f"per_site_ref_frac_cyp21a2[{_fmt_stats(site_frac_cyp21a2)}] "
                f"per_site_ref_frac_cyp21a1p[{_fmt_stats(site_frac_cyp21a1p)}]"
            )

            print("##" + summary_body, file=sys.stderr)
            print("##" + trace_body, file=sys.stderr)
            cn_summary_for_file = "# " + summary_body + "\n# " + trace_body
        else:
            cn_summary_for_file = None

    finally:
        bam.close()
        if args.output:
            out.close()
        if not args.no_cn_summary and cn_summary_for_file is not None:
            if args.output:
                with open(args.output, "a") as fa:
                    for line in cn_summary_for_file.split("\n"):
                        fa.write(line + "\n")
            else:
                for line in cn_summary_for_file.split("\n"):
                    print(line, file=sys.stdout)


if __name__ == "__main__":
    main()
