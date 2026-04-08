#!/usr/bin/env python3
"""
Pileup A/C/G/T at HBA1 vs HBA2 paralog mismatch sites (from extract_hba_paralog_sites.py).

Uses the same reference-forward / HGVS coding logic as smaca_append_summary.py pileup
for SMN c.840: bases are counted in transcript (coding) sense when the gene is on
the minus strand. HBA1/HBA2 on chr16 are forward (+) in GRCh38 — default is
hgvs_coding_minus_strand=False (reference forward = coding).

Example:
  python3 extract_hba_paralog_sites.py --fasta GRCh38.fa -o sites.tsv
  python3 hba_paralog_pileup.py --bam sample.bam --sites sites.tsv

Summarizes across all sites:
  ratio_hba1_ref_mass = sum(reads matching HBA1 ref) / (sum HBA1 ref + sum HBA2 ref reads).
  est_copies_hba1 = total_alpha_copies * ratio (default total_alpha_copies=4 → ~2 vs ~2 when balanced).
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

import pysam


def pick_chr16(bam: pysam.AlignmentFile) -> str | None:
    for c in ("chr16", "16"):
        if c in bam.references:
            return c
    return None


def normalize_contig(refs: set[str], chrom: str) -> str | None:
    """Map chr16/16 to whichever exists in the BAM header."""
    if chrom in refs:
        return chrom
    if chrom == "chr16" and "16" in refs:
        return "16"
    if chrom == "16" and "chr16" in refs:
        return "chr16"
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
        help="TSV from extract_hba_paralog_sites.py (tab-separated, header row)",
    )
    ap.add_argument(
        "--hgvs-minus-strand",
        action="store_true",
        help="Tally bases in HGVS coding sense for minus-strand genes (HBA is + strand; leave off)",
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
        help="Do not print ##HBA_PARALOG_PILEUP lines to stderr (CN summary still printed)",
    )
    ap.add_argument(
        "--no-cn-summary",
        action="store_true",
        help="Do not print ##HBA_PARALOG_CN_RATIO line or append summary to -o file",
    )
    ap.add_argument(
        "--total-alpha-copies",
        type=int,
        default=4,
        help="Assumed total alpha-globin gene copies in diploid ref (default 4 = 2xHBA1 + 2xHBA2) for scaling ratio→copies",
    )
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_refs = set(bam.references)
    default_chrom = pick_chr16(bam)
    if not default_chrom:
        bam.close()
        print("BAM has no chr16 or 16 contig.", file=sys.stderr)
        sys.exit(1)

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    cn_summary_for_file: str | None = None
    try:
        fieldnames = [
            "idx",
            "chrom",
            "hba1_pos",
            "hba1_ref",
            "hba1_ref_frac",
            "hba1_pileup",
            "hba2_pos",
            "hba2_ref",
            "hba2_ref_frac",
            "hba2_pileup",
        ]
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()

        sum_hba1_ref_reads = 0
        sum_hba2_ref_reads = 0
        sum_hba1_depth = 0
        sum_hba2_depth = 0
        n_rows = 0
        site_frac_hba1: list[float] = []
        site_frac_hba2: list[float] = []

        with open(args.sites, newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for idx, row in enumerate(rdr, start=1):
                if args.max_sites and idx > args.max_sites:
                    break
                try:
                    p1 = int(row["hba1_pos"])
                    p2 = int(row["hba2_pos"])
                    r1 = row["hba1_ref"].strip().upper()
                    r2 = row["hba2_ref"].strip().upper()
                except (KeyError, ValueError) as e:
                    print(f"Skip row {idx}: {e}", file=sys.stderr)
                    continue

                chrom = (row.get("hba1_chr") or row.get("hba2_chr") or "").strip()
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
                    sum_hba1_ref_reads += int(c1[r1])
                if c2 and r2 in c2:
                    sum_hba2_ref_reads += int(c2[r2])
                sum_hba1_depth += depth(c1)
                sum_hba2_depth += depth(c2)
                n_rows += 1
                f1 = frac_ref_float(c1, r1)
                f2 = frac_ref_float(c2, r2)
                if f1 is not None:
                    site_frac_hba1.append(f1)
                if f2 is not None:
                    site_frac_hba2.append(f2)

                w.writerow(
                    {
                        "idx": idx,
                        "chrom": chrom,
                        "hba1_pos": p1,
                        "hba1_ref": r1,
                        "hba1_ref_frac": frac_ref(c1, r1),
                        "hba1_pileup": fmt_acgt(c1),
                        "hba2_pos": p2,
                        "hba2_ref": r2,
                        "hba2_ref_frac": frac_ref(c2, r2),
                        "hba2_pileup": fmt_acgt(c2),
                    }
                )

                if not args.quiet:
                    # Grep-friendly summary (same spirit as ##C840 lines)
                    print(
                        f"##HBA_PARALOG_PILEUP idx={idx} chrom={chrom} "
                        f"hba1_pos={p1} expected_ref={r1} | hba2_pos={p2} expected_ref={r2}",
                        file=sys.stderr,
                    )
                    print(
                        f"  HBA1_coding {fmt_acgt(c1)} | HBA2_coding {fmt_acgt(c2)}",
                        file=sys.stderr,
                    )

        if not args.no_cn_summary:
            tot_ref = sum_hba1_ref_reads + sum_hba2_ref_reads
            if tot_ref > 0:
                ratio_hba1 = sum_hba1_ref_reads / tot_ref
                tc = args.total_alpha_copies
                est_hba1 = (tc * ratio_hba1) if tc else float("nan")
                est_hba2 = tc - est_hba1 if tc else float("nan")
                ratio_s = f"{ratio_hba1:.6f}"
                est1_s = f"{est_hba1:.4f}"
                est2_s = f"{est_hba2:.4f}"
            else:
                ratio_hba1 = None
                ratio_s = "NA"
                est1_s = "NA"
                est2_s = "NA"

            summary_body = (
                "HBA_PARALOG_CN_RATIO "
                f"hba1_ref_reads_total={sum_hba1_ref_reads} "
                f"hba2_ref_reads_total={sum_hba2_ref_reads} "
                f"ratio_hba1_ref_mass={ratio_s} "
                f"(HBA1_ref/(HBA1_ref+HBA2_ref) over {n_rows} sites; "
                f"est_copies_hba1={est1_s} est_copies_hba2={est2_s} "
                f"scaled_from_total_alpha_copies={args.total_alpha_copies}) "
                f"sum_depth_hba1={sum_hba1_depth} sum_depth_hba2={sum_hba2_depth}"
            )

            def _fmt_stats(xs: list[float]) -> str:
                if not xs:
                    return "min=NA mean=NA max=NA n=0"
                return (
                    f"min={min(xs):.4f} mean={statistics.mean(xs):.4f} "
                    f"max={max(xs):.4f} n={len(xs)}"
                )

            tot_depth = sum_hba1_depth + sum_hba2_depth
            depth_ratio_hba1 = (
                (sum_hba1_depth / tot_depth) if tot_depth > 0 else float("nan")
            )
            pooled1 = (
                (sum_hba1_ref_reads / sum_hba1_depth) if sum_hba1_depth > 0 else None
            )
            pooled2 = (
                (sum_hba2_ref_reads / sum_hba2_depth) if sum_hba2_depth > 0 else None
            )
            frac_hba2_mass = (
                (sum_hba2_ref_reads / tot_ref) if tot_ref > 0 else float("nan")
            )
            p1s = f"{pooled1:.6f}" if pooled1 is not None else "NA"
            p2s = f"{pooled2:.6f}" if pooled2 is not None else "NA"
            drs = f"{depth_ratio_hba1:.6f}" if tot_depth > 0 else "NA"
            f2ms = f"{frac_hba2_mass:.6f}" if tot_ref > 0 else "NA"

            r1_over_tot = (
                f"{sum_hba1_ref_reads}/{tot_ref}" if tot_ref > 0 else "NA/NA"
            )
            r2_over_tot = (
                f"{sum_hba2_ref_reads}/{tot_ref}" if tot_ref > 0 else "NA/NA"
            )
            trace_body = (
                "HBA_PARALOG_CN_TRACE "
                f"tot_ref_reads_denominator={tot_ref} "
                f"(hba1_ref_reads_total+hba2_ref_reads_total) "
                f"frac_hba1_ref_mass={r1_over_tot}={ratio_s} "
                f"frac_hba2_ref_mass={r2_over_tot}={f2ms} "
                f"formula_est_hba1={args.total_alpha_copies}*"
                f"{sum_hba1_ref_reads}/({sum_hba1_ref_reads}+{sum_hba2_ref_reads})={est1_s} "
                f"formula_est_hba2={args.total_alpha_copies}*"
                f"{sum_hba2_ref_reads}/({sum_hba1_ref_reads}+{sum_hba2_ref_reads})={est2_s} "
                f"sum_depth_total={tot_depth} "
                f"depth_mass_ratio_hba1={drs} "
                f"(sum_depth_hba1/sum_depth_total; compare to ratio_hba1_ref_mass={ratio_s}) "
                f"pooled_ref_reads_over_depth_hba1={p1s} "
                f"pooled_ref_reads_over_depth_hba2={p2s} "
                f"per_site_ref_frac_hba1[{_fmt_stats(site_frac_hba1)}] "
                f"per_site_ref_frac_hba2[{_fmt_stats(site_frac_hba2)}]"
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
