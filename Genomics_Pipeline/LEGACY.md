# LEGACY: Genomics Pipeline Reference Code

This directory contains the original Genomics Pipeline scripts that were used as reference for implementing the Carrier Screening service modules in `genolyx/service-daemon`.

## Migration Status

The core logic from each subdirectory has been integrated into the service-daemon's carrier_screening plugin:

| Original Directory | Migrated To (service-daemon) | Status |
|---|---|---|
| `phenotype_portal/main.py` | `app/services/carrier_screening/vcf_parser.py` | Migrated (VCF parsing, HPO/AF/ClinVar filtering) |
| `phenotype_portal/main.py` | `app/services/carrier_screening/annotator.py` | Migrated (ClinVar, gnomAD, snpEff annotation) |
| `Carrier_result/generate.py` | `app/services/carrier_screening/report.py` | Migrated (Jinja2 HTML template, PDF generation) |
| `Carrier_result/generate.py` | `app/services/carrier_screening/review.py` | Migrated (disease-gene mapping, result.json) |
| `genetic_reporter/acmg_classifier.py` | `app/services/carrier_screening/acmg.py` | Migrated (ACMG classification) |

## Key Differences

The service-daemon implementation improves upon the original in several ways. First, it uses a modular class-based architecture instead of monolithic Flask routes. Second, it supports multiple annotation sources (ClinVar, gnomAD, ClinGen, MANE, HPO, curated DB, HGMD) through a unified `VariantAnnotator` facade. Third, it provides configurable filtering through `VariantFilterConfig` with HPO-based gene filtering, max AF thresholds, and ClinVar significance filters. Fourth, it generates multi-language reports (EN/CN/KO) with couple report support. Finally, it produces a structured `result.json` for Portal review with disease grouping and carrier status determination.

## Do Not Modify

These files are kept for historical reference. All active development should be done in the `genolyx/service-daemon` repository.
