#!/usr/bin/env python3
# Source of truth for the base64 payload in coverage.nf (process MOSDEPTH_PER_BASE_QC); keep them in sync.
"""
Build per-gene coverage TSV from mosdepth per-base (tabix-indexed) + panel BED.

BED column 4 is treated as the gene key (HGNC-style targets). Intervals per gene are
merged on each chromosome; per-base rows are clipped to those intervals.

Output columns match daemon-friendly names: gene, mean_coverage, pct_bases_10x, pct_bases_20x.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


Interval = Tuple[str, int, int]  # chrom, start, end (0-based half-open)

_GENE_SYM = re.compile(r"^[A-Za-z][A-Za-z0-9-]{0,24}$")
_SKIP_TOK_PREFIXES = (
    "ENST",
    "NM_",
    "NR_",
    "XM_",
    "XR_",
    "ENSG",
    "CCDS",
    "CLINID",
    "LOC",
)


def _gene_keys_from_bed_col4(name: str) -> List[str]:
    """
    HGNC-like symbols in column 4. Twist/vendor BEDs use ``PAH;NM_...;ENST...``; plain BEDs use ``PAH``.
    """
    raw = (name or "").strip()
    if not raw:
        return []
    keys: List[str] = []
    for chunk in raw.replace(",", ";").replace("|", ";").split(";"):
        t = chunk.strip()
        if not t:
            continue
        ul = t.upper()
        if any(ul.startswith(p) for p in _SKIP_TOK_PREFIXES):
            continue
        if _GENE_SYM.match(t):
            keys.append(ul)
    if keys:
        return keys
    if _GENE_SYM.match(raw):
        return [raw.upper()]
    return []


def _parse_bed(path: str) -> Dict[str, List[Interval]]:
    """Group intervals by HGNC symbol derived from column 4."""
    by_gene: Dict[str, List[Interval]] = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track "):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, s, e, name = parts[0], parts[1], parts[2], parts[3].strip()
            if not name:
                continue
            try:
                start_i = int(s)
                end_i = int(e)
            except ValueError:
                continue
            if end_i <= start_i:
                continue
            for gkey in _gene_keys_from_bed_col4(name):
                by_gene[gkey].append((chrom, start_i, end_i))
    return by_gene


def _merge_intervals(intervals: List[Interval]) -> List[Interval]:
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, s, e in intervals:
        by_chrom[chrom].append((s, e))
    merged: List[Interval] = []
    for chrom, ivs in sorted(by_chrom.items(), key=lambda x: x[0]):
        ivs.sort(key=lambda x: (x[0], x[1]))
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e))
    return merged


def _alt_chrom(c: str) -> str:
    if c.startswith("chr"):
        return c[3:] if len(c) > 3 else c
    return "chr" + c


def _tabix_lines(bgz: str, chrom: str, start: int, end: int) -> List[str]:
    """Query per-base BED; use 0-based half-open intervals (matches BED storage)."""
    for ch in (chrom, _alt_chrom(chrom)):
        r = f"{ch}:{start}-{end}"
        try:
            out = subprocess.check_output(
                ["tabix", "--zero-based", bgz, r],
                stderr=subprocess.DEVNULL,
                text=True,
            )
        except subprocess.CalledProcessError:
            continue
        if out.strip():
            return out.splitlines()
    return []


def _accumulate_interval(
    bgz: str,
    chrom: str,
    gs: int,
    ge: int,
    total_bases: int,
    sum_depth: int,
    bases_ge_10: int,
    bases_ge_20: int,
) -> Tuple[int, int, int, int]:
    """Aggregate depth over [gs, ge). Interval length is always (ge - gs); gaps in tabix = depth 0."""
    span = ge - gs
    if span <= 0:
        return total_bases, sum_depth, bases_ge_10, bases_ge_20
    covered = 0
    for line in _tabix_lines(bgz, chrom, gs, ge):
        parts = line.rstrip().split("\t")
        if len(parts) < 4:
            continue
        try:
            rs = int(parts[1])
            re_ = int(parts[2])
            depth = int(float(parts[3]))
        except (ValueError, IndexError):
            continue
        os_ = max(rs, gs)
        oe = min(re_, ge)
        if oe <= os_:
            continue
        ln = oe - os_
        covered += ln
        sum_depth += depth * ln
        if depth >= 10:
            bases_ge_10 += ln
        if depth >= 20:
            bases_ge_20 += ln
    uncovered = span - min(covered, span)
    total_bases += span
    # uncovered bases contribute 0 depth and do not count toward >=10 / >=20
    return total_bases, sum_depth, bases_ge_10, bases_ge_20


def write_tsv(bgz: str, bed: str, out_path: str) -> None:
    genes = _parse_bed(bed)
    rows: List[Tuple[str, float, float, float]] = []
    for gene in sorted(genes.keys()):
        merged = _merge_intervals(genes[gene])
        tb = sd = b10 = b20 = 0
        for chrom, gs, ge in merged:
            tb, sd, b10, b20 = _accumulate_interval(bgz, chrom, gs, ge, tb, sd, b10, b20)
        if tb <= 0:
            mean_c = 0.0
            p10 = p20 = 0.0
        else:
            mean_c = sd / tb
            p10 = 100.0 * b10 / tb
            p20 = 100.0 * b20 / tb
        rows.append((gene, mean_c, p10, p20))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("gene\tmean_coverage\tpct_bases_10x\tpct_bases_20x\n")
        for gene, mean_c, p10, p20 in rows:
            f.write(f"{gene}\t{mean_c:.4f}\t{p10:.4f}\t{p20:.4f}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Gene-level coverage TSV from mosdepth per-base + BED.")
    ap.add_argument("per_base_bgz", help="mosdepth *.per-base.bed.gz (tabix indexed)")
    ap.add_argument("gene_bed", help="BED with HGNC (or label) in column 4")
    ap.add_argument("out_tsv", help="Output TSV path")
    args = ap.parse_args()
    try:
        write_tsv(args.per_base_bgz, args.gene_bed, args.out_tsv)
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
