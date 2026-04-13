#!/usr/bin/env python3
"""
Infer APOE ε2 / ε3 / ε4 isoforms from the two defining SNPs (GRCh38):
  rs429358 — chr19:44908684 T>C
  rs7412  — chr19:44908822 C>T

Cis haplotypes (chromosome):
  ε4: C at rs429358, C at rs7412
  ε3: T at rs429358, C at rs7412
  ε2: T at rs429358, T at rs7412

Unphased genotypes can match more than one diplotype (notably ε3/ε3 vs ε2/ε4).
Output is informational; not a clinical diagnosis — see JSON disclaimer.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from typing import Any

# GRCh38 (1-based), must match pgx_custom_variants.tsv
RS429358 = ("chr19", 44908684, "T", "C")
RS7412 = ("chr19", 44908822, "C", "T")

ISO_TO_BASES = {
    "ε2": ("T", "T"),
    "ε3": ("T", "C"),
    "ε4": ("C", "C"),
}


def _norm_chrom(c: str) -> str:
    c = c.strip()
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def _iso_sort_key(x: str) -> int:
    return ("ε2", "ε3", "ε4").index(x)


def _get_genotype_at(
    vcf_path: str,
    chrom: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
) -> tuple[list[str] | None, bool | None, str | None]:
    """Allele bases [a1, a2], phased flag, error."""
    try:
        import pysam
    except ImportError:
        return None, None, "pysam_not_installed"

    chrom = _norm_chrom(chrom)
    vcf = pysam.VariantFile(vcf_path)
    try:
        for rec in vcf.fetch(chrom, pos - 1, pos):
            if rec.pos != pos:
                continue
            for sample_name in rec.samples:
                s = rec.samples[sample_name]
                gt = s.get("GT")
                if gt is None or any(x is None for x in gt):
                    return None, None, "missing_gt"
                phased = bool(s.get("phased", False))
                alleles = [rec.alleles[i] for i in gt]
                if len(alleles) != 2:
                    return None, None, "invalid_gt"
                for a in alleles:
                    if len(a) != 1:
                        return None, None, "indel_at_site"
                    au = a.upper()
                    if au not in (ref_allele.upper(), alt_allele.upper()):
                        return None, None, f"unexpected_allele:{au}"
                return [alleles[0].upper(), alleles[1].upper()], phased, None
        # No VCF row: homozygous reference
        return [ref_allele.upper(), ref_allele.upper()], False, None
    finally:
        vcf.close()


def _haplotype_from_chr(b429: str, b7412: str) -> str | None:
    for iso, (a429, a7412) in ISO_TO_BASES.items():
        if b429 == a429 and b7412 == a7412:
            return iso
    return None


def _nt_counts(bases: list[str], ref: str, alt: str) -> tuple[int, int] | None:
    n_ref = sum(1 for x in bases if x == ref)
    n_alt = sum(1 for x in bases if x == alt)
    if n_ref + n_alt != 2:
        return None
    return n_ref, n_alt


def _enumerate_diplotypes(
    nC_429: int,
    nT_429: int,
    nC_7412: int,
    nT_7412: int,
) -> list[tuple[str, str]]:
    """All unordered ε/ε pairs consistent with observed base counts."""
    alleles = ("ε2", "ε3", "ε4")
    valid: set[tuple[str, str]] = set()
    for h1, h2 in product(alleles, repeat=2):
        a429_1, a7412_1 = ISO_TO_BASES[h1]
        a429_2, a7412_2 = ISO_TO_BASES[h2]
        c429 = (1 if a429_1 == "C" else 0) + (1 if a429_2 == "C" else 0)
        t429 = 2 - c429
        c7412 = (1 if a7412_1 == "C" else 0) + (1 if a7412_2 == "C" else 0)
        t7412 = 2 - c7412
        if (c429, t429, c7412, t7412) == (nC_429, nT_429, nC_7412, nT_7412):
            pair = tuple(sorted([h1, h2], key=_iso_sort_key))
            valid.add(pair)
    return sorted(valid, key=lambda p: (_iso_sort_key(p[0]), _iso_sort_key(p[1])))


def _risk_note(diplotypes: list[tuple[str, str]], resolved: bool) -> str:
    if not diplotypes:
        return "Could not resolve APOE isoforms from genotypes."
    eps4_max = 0
    for a, b in diplotypes:
        eps4_max = max(eps4_max, (a == "ε4") + (b == "ε4"))
    lines = [
        "APOE ε4 is associated with higher population risk of late-onset Alzheimer disease (dose-related); "
        "ε2 may be associated with lower risk relative to ε3. This is not a diagnosis or deterministic prediction.",
    ]
    if not resolved:
        lines.append(
            "Genotypes did not fully phase: more than one ε diplotype may fit the observed SNPs — "
            "interpret with caution or use phased data or orthogonal assay if clinically required.",
        )
    if eps4_max == 0:
        lines.append("Resolved possibility(ies) carry no ε4 allele (relative to common ε3 risk).")
    elif eps4_max == 1:
        lines.append("At least one resolved diplotype includes one ε4 allele.")
    else:
        lines.append("At least one resolved diplotype includes two ε4 alleles.")
    return " ".join(lines)


def run(vcf_path: str, sample_id: str) -> dict[str, Any]:
    c1, p1, r1, a1 = RS429358
    c2, p2, r2, a2 = RS7412

    g429, ph429, e1 = _get_genotype_at(vcf_path, c1, p1, r1, a1)
    g7412, ph7412, e2 = _get_genotype_at(vcf_path, c2, p2, r2, a2)

    out: dict[str, Any] = {
        "schema": "gx_exome_apoe_v1",
        "sample_id": sample_id,
        "reference_genome": "GRCh38",
        "snps": {
            "rs429358": {
                "chrom": c1,
                "pos": p1,
                "genotype_bases": g429,
                "phased": ph429,
                "error": e1,
            },
            "rs7412": {
                "chrom": c2,
                "pos": p2,
                "genotype_bases": g7412,
                "phased": ph7412,
                "error": e2,
            },
        },
        "disclaimer": (
            "Informational research use. APOE genotyping from exome VCF does not replace "
            "clinical-grade assays; phasing may be ambiguous. Not for diagnosis or medical decisions without a qualified professional."
        ),
    }

    if e1 or e2 or g429 is None or g7412 is None:
        out["status"] = "error"
        out["error"] = e1 or e2 or "genotype_missing"
        return out

    # Fully phased at both loci: cis haplotypes from allele order
    if ph429 and ph7412:
        h1 = _haplotype_from_chr(g429[0], g7412[0])
        h2 = _haplotype_from_chr(g429[1], g7412[1])
        if h1 and h2:
            pair = tuple(sorted([h1, h2], key=_iso_sort_key))
            out["status"] = "ok"
            out["phasing"] = "phased"
            out["diplotypes_possible"] = [list(pair)]
            out["diplotype_display"] = f"{pair[0]}/{pair[1]}"
            out["risk_context"] = _risk_note([pair], resolved=True)
            return out
        out["phasing_note"] = (
            "VCF marked phased but rs429358/rs7412 pair does not match canonical ε2/ε3/ε4; "
            "falling back to unphased enumeration."
        )

    # rs429358: ref T alt C — count C vs T
    ct429 = _nt_counts(g429, "T", "C")
    ct7412 = _nt_counts(g7412, "C", "T")
    if ct429 is None or ct7412 is None:
        out["status"] = "error"
        out["error"] = "invalid_base_counts"
        return out

    nT_429, nC_429 = ct429[0], ct429[1]  # ref T, alt C
    nC_7412, nT_7412 = ct7412[0], ct7412[1]  # ref C, alt T

    diplotypes = _enumerate_diplotypes(nC_429, nT_429, nC_7412, nT_7412)
    if not diplotypes:
        out["status"] = "error"
        out["error"] = "no_valid_cis_haplotype_pair"
        out["note"] = (
            "Observed genotypes do not match any combination of canonical ε2/ε3/ε4 haplotypes "
            "(rare recombinant or genotyping error)."
        )
        return out

    resolved = len(diplotypes) == 1
    displays = [f"{a}/{b}" for a, b in diplotypes]
    out["status"] = "ok"
    out["phasing"] = "unphased" if not resolved else "inferred_unique"
    out["diplotypes_possible"] = [list(p) for p in diplotypes]
    out["diplotype_display"] = displays[0] if resolved else " OR ".join(displays)
    out["risk_context"] = _risk_note(diplotypes, resolved=resolved)
    return out


def write_summary(path: str, data: dict[str, Any]) -> None:
    lines = [
        f"APOE isoforms (Alzheimer risk context) — {data.get('sample_id', '')}",
        "=" * 72,
    ]
    snps = data.get("snps", {})
    for label in ("rs429358", "rs7412"):
        s = snps.get(label, {})
        g = s.get("genotype_bases")
        lines.append(f"{label}: {g if g else 'N/A'}  (GRCh38 {s.get('chrom')}:{s.get('pos')})")
    lines.append("")
    if data.get("status") != "ok":
        lines.append(f"Status: {data.get('status', 'unknown')}")
        if data.get("error"):
            lines.append(f"Error: {data['error']}")
        if data.get("note"):
            lines.append(data["note"])
    else:
        lines.append(f"Diplotype: {data.get('diplotype_display', '')}")
        lines.append(f"Phasing: {data.get('phasing', '')}")
        if len(data.get("diplotypes_possible") or []) > 1:
            lines.append("Possible diplotypes (unphased):")
            for p in data["diplotypes_possible"]:
                lines.append(f"  — {'/'.join(p)}")
        lines.append("")
        lines.append(data.get("risk_context", ""))
    lines.append("")
    lines.append(data.get("disclaimer", ""))
    lines.append("")
    lines.append("Full JSON: apoe/apoe_result.json")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _fmt_gt_bases(g: Any) -> str:
    if isinstance(g, list) and len(g) == 2:
        return f"{g[0]}/{g[1]}"
    return str(g) if g else "N/A"


def write_review(path: str, data: dict[str, Any]) -> None:
    """Standalone review document for clinical/lab review folders (also embedded in detailed report)."""
    sid = data.get("sample_id", "")
    lines: list[str] = [
        "=" * 78,
        "APOE REVIEW — Alzheimer disease risk context (isoforms ε2, ε3, ε4)",
        "=" * 78,
        "",
        f"Sample: {sid}",
        "Reference: GRCh38",
        "Defining variants: rs429358 (chr19:44908684), rs7412 (chr19:44908822)",
        "",
        "—" * 78,
        "1. GENOTYPES",
        "—" * 78,
    ]
    snps = data.get("snps", {})
    for label in ("rs429358", "rs7412"):
        s = snps.get(label, {})
        g = s.get("genotype_bases")
        lines.append(
            f"  {label}  {s.get('chrom', '')}:{s.get('pos', '')}  →  {_fmt_gt_bases(g)}"
        )
    lines.extend(["", "—" * 78, "2. ISOFORM / DIPLOTYPE", "—" * 78])

    st = data.get("status")
    if st == "skipped":
        lines.append(f"  {data.get('message', 'APOE not run.')}")
    elif st != "ok":
        lines.append(f"  Status: {st}")
        if data.get("error"):
            lines.append(f"  Error: {data['error']}")
        if data.get("note"):
            lines.append(f"  Note: {data['note']}")
    else:
        lines.append(f"  Diplotype: {data.get('diplotype_display', '')}")
        lines.append(f"  Phasing: {data.get('phasing', '')}")
        if len(data.get("diplotypes_possible") or []) > 1:
            lines.append("  Possible diplotypes (if unphased):")
            for p in data["diplotypes_possible"]:
                lines.append(f"    • {'/'.join(p)}")
        if data.get("phasing_note"):
            lines.append(f"  Note: {data['phasing_note']}")

    lines.extend(["", "—" * 78, "3. INTERPRETATION (population risk, not diagnosis)", "—" * 78])
    if st == "ok":
        lines.append("  " + data.get("risk_context", "").replace("\n", "\n  "))
    else:
        lines.append("  (No isoform interpretation — genotyping incomplete or skipped.)")

    lines.extend(["", "—" * 78, "4. LIMITATIONS", "—" * 78])
    lines.append("  " + data.get("disclaimer", "").replace("\n", "\n  "))
    lines.extend([
        "",
        "—" * 78,
        "5. METHOD",
        "—" * 78,
        "  Isoforms inferred from exome VCF genotypes at rs429358 and rs7412 (cis rules for ε2/ε3/ε4).",
        "  Not a substitute for clinical-grade assays when required.",
        "",
        f"Machine-readable: apoe/apoe_result.json",
        "=" * 78,
        "",
    ])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="APOE ε2/ε3/ε4 from VCF (GRCh38)")
    p.add_argument("--vcf", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--json-out", default="apoe_result.json")
    p.add_argument("--summary-out", default="apoe_summary.txt")
    p.add_argument("--review-out", default="apoe_review.txt")
    args = p.parse_args()

    data = run(args.vcf, args.sample_id)
    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    write_summary(args.summary_out, data)
    write_review(args.review_out, data)


if __name__ == "__main__":
    main()
