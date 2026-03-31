#!/usr/bin/env python3
"""
BAM pileup at the seven NM_000500.9 / GRCh38 chr6 sites in cyp21a2_hotspots.tsv.

Uses the same reference-forward pileup semantics as cyp21_paralog_pileup.py (CYP21A2 is + strand).
Outputs per-site A/C/G/T counts, depth, ref/expected-alt fractions for cross-check vs variant VCF calls.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pysam

# Same directory as this script (Nextflow stages both .py files in task work dir)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cyp21_paralog_pileup import (  # noqa: E402
    depth,
    fmt_acgt,
    normalize_contig,
    pick_chr6,
    pileup_acgt_at_1bp,
)


def load_hotspot_rows(tsv_path: Path) -> list[dict]:
    rows: list[dict] = []
    with open(tsv_path, newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6 or parts[0] == "id":
                continue
            hid, label, chrom, pos_s, ref, alts = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
                parts[4],
                parts[5],
            )
            note = parts[6] if len(parts) > 6 else ""
            alt_list = [a.strip().upper() for a in alts.split(",") if a.strip()]
            rows.append(
                {
                    "id": hid,
                    "label": label,
                    "chrom": chrom,
                    "pos": int(pos_s),
                    "ref": ref.upper(),
                    "alts": alt_list,
                    "note": note,
                }
            )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True, help="Indexed BAM aligned to GRCh38")
    ap.add_argument(
        "--sites",
        required=True,
        type=Path,
        help="cyp21a2_hotspots.tsv (tab-separated)",
    )
    ap.add_argument("-o", "--output", required=True, type=Path, help="Output TSV path")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    default_chrom = pick_chr6(bam)
    bam.close()
    if not default_chrom:
        print("BAM has no chr6 or 6 contig.", file=sys.stderr)
        sys.exit(1)

    sites = load_hotspot_rows(args.sites)
    if not sites:
        print("No hotspot rows loaded from --sites.", file=sys.stderr)
        sys.exit(1)

    fieldnames = [
        "id",
        "label",
        "chrom",
        "pos",
        "ref",
        "expected_alts",
        "A",
        "C",
        "G",
        "T",
        "N",
        "depth",
        "ref_count",
        "alt_expected_count",
        "ref_frac",
        "alt_frac",
        "pileup",
    ]

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()

        for h in sites:
            chrom = h["chrom"]
            chrom = normalize_contig(args.bam, chrom) or default_chrom
            pos = h["pos"]
            ref = h["ref"]
            alts = h["alts"]

            ct = pileup_acgt_at_1bp(
                args.bam,
                chrom,
                pos,
                hgvs_coding_minus_strand=False,
            )
            d = depth(ct) if ct else 0
            ref_count = int(ct[ref]) if ct and ref in ct else 0
            alt_expected_count = 0
            if ct:
                for a in alts:
                    if a in ct:
                        alt_expected_count += int(ct[a])

            ref_frac = f"{ref_count / d:.6f}" if d else "NA"
            alt_frac = f"{alt_expected_count / d:.6f}" if d else "NA"

            w.writerow(
                {
                    "id": h["id"],
                    "label": h["label"],
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "expected_alts": ",".join(alts),
                    "A": ct["A"] if ct else 0,
                    "C": ct["C"] if ct else 0,
                    "G": ct["G"] if ct else 0,
                    "T": ct["T"] if ct else 0,
                    "N": ct["N"] if ct else 0,
                    "depth": d,
                    "ref_count": ref_count,
                    "alt_expected_count": alt_expected_count,
                    "ref_frac": ref_frac,
                    "alt_frac": alt_frac,
                    "pileup": fmt_acgt(ct),
                }
            )

        out.write(
            "# ##CYP21_HOTSPOT_PILEUP_GRCh38 chr6 reference coordinates; "
            "A/C/G/T in reference forward strand (matches cyp21_paralog_pileup). "
            "expected_alts = catalogued hotspot alternate(s) from cyp21a2_hotspots.tsv (not inferred from data). "
            "alt_frac = reads matching those alts / depth. "
            "Compare alt_frac to variant caller GT at the same positions.\n"
        )


if __name__ == "__main__":
    main()
