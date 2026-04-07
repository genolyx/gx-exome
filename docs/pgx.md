# Pharmacogenomics (PGx) add-on (PharmCAT)

PGx runs **by default** on every exome workflow (reflex from exome). It does **not** replace alignment or variant calling; it only consumes the filtered (and VEP-annotated, if enabled) VCF. Turn it off with **`--skip_pgx true`** when needed.

## When it runs

- **Default: on** for every exome run (reflex from exome). Controlled by Nextflow **`--skip_pgx`** (default: `false`).
- To disable from `run_analysis.sh`: **`--skip-pgx`** (passes `--skip_pgx true` to Nextflow).
- Requires **`--backbone_bed`** (capture BED) so variant calling ran; PGx uses the same final VCF channel as the rest of the workflow.

## Tooling (pinned stack)

| Item | Value |
|------|--------|
| Tool | [PharmCAT](https://pharmcat.org/) (`pharmcat_pipeline`) |
| Default container | `pgkb/pharmcat:2.15.5` |
| Override | `--pgx_container '<image:tag>'` |
| Reference label | `--pgx_reference_genome` (default `GRCh38`) — must match the VCF / GRCh38 alignment |
| CPIC / PharmGKB data | Bundled with the PharmCAT release (see `pgx_meta.json` / PharmCAT docs) |

PharmCAT applies its own VCF preprocessor for PGx positions; genome build must be consistent with the exome pipeline (this repo assumes **GRCh38**).

## Process

- **Nextflow process**: `RUN_PGX_PHARMCAT` in `bin/modules/pgx.nf`
- **Label**: `pgx` (Docker image from `params.pgx_container`)
- **Failure behavior**: The process **always exits 0** (soft fail). If PharmCAT fails, `pgx/pgx_meta.json` has `"exit_status": "error"` and a message; the main exome workflow still completes.

## Output contract (for service-daemon)

All artifacts are written under a single directory:

| Path (under `outdir` and mirrored to `output_dir`) | Description |
|-----------------------------------------------------|-------------|
| `pgx/pgx_meta.json` | Tool name, version string, container image, reference genome, resource notes, sample id, ISO timestamp, `exit_status` (`success` / `error`), optional `message` |
| `pgx/pgx_result.json` | Normalized summary: schema id, sample id, artifact filenames, and embedded **named allele matcher** / **phenotype** (and **reporter** JSON when present) from PharmCAT |
| `pgx/<sample>_pgx.match.json` | PharmCAT named allele matcher (when success) |
| `pgx/<sample>_pgx.phenotype.json` | Phenotyper output (when success) |
| `pgx/<sample>_pgx.report.json` | Reporter JSON (when `-reporterJson` produces it) |
| `pgx/<sample>_pgx.report.html` | HTML report (optional for portals) |
| `pgx/pharmcat.stderr.log` | PharmCAT stderr (useful when `exit_status` is `error`) |
| `pgx/pharmcat_version.txt` | `pharmcat_pipeline -V` line (best-effort) |

Downstream can glob `**/pgx/pgx_meta.json` or `analysis_dir/pgx/` relative to the job root.

## CYP2D6 — Aldy integration

PharmCAT **does not call CYP2D6 from VCF** in its default mode because CYP2D6 phenotype prediction requires structural-variant / copy-number awareness that standard VCFs cannot represent. The pipeline uses **Aldy** as a dedicated CYP2D6 star-allele caller.

### How it works

1. **`RUN_ALDY_CYP2D6`** (`bin/modules/aldy.nf`, label `aldy`) runs on the mark-duplicated BAM (in parallel with variant calling / VEP).
2. Aldy calls CYP2D6 diplotypes using the **exome** profile (`-p exome -g CYP2D6`).
3. The diplotype is written as a **PharmCAT outside-call TSV** (`CYP2D6\t*1/*4`).
4. **`RUN_PGX_PHARMCAT`** receives the TSV via the `-po` flag, so PharmCAT translates the CYP2D6 diplotype into phenotype + CPIC drug recommendations alongside all other PGx genes.
5. If Aldy produces no call (e.g., insufficient coverage), the TSV is empty and PharmCAT runs without CYP2D6 outside calls.

### Exome limitations

With exome data, Aldy **assumes exactly 2 gene copies** — it cannot detect:
- Gene deletions (`*5`)
- Gene duplications (ultrarapid metabolizer configurations)
- CYP2D6/CYP2D7 hybrid alleles

SNP/indel-based star alleles (`*3`, `*4`, `*6`, `*9`, `*10`, `*17`, `*41`, etc.) are still called and cover the most clinically relevant variants.

### Container

| Item | Value |
|------|-------|
| Tool | [Aldy](https://github.com/0xTCG/aldy) (pip-installed at runtime) |
| Container | `python:3.12-bookworm` (shared with `pgx_finalize`) |
| License | Free for academic / non-commercial use |

### Output

| Path (under `outdir`) | Description |
|------------------------|-------------|
| `pgx/aldy_cyp2d6_outside_call.tsv` | PharmCAT outside-call TSV (copied into `pgx_staging` by `RUN_PGX_PHARMCAT`) |
| Work dir only: `aldy_cyp2d6_result.json` | Aldy result summary (schema `gx_exome_aldy_cyp2d6_v1`) |
| Work dir only: `<sample>.CYP2D6.aldy` | Raw Aldy output with minor-allele decomposition |

## Implementation notes

- The pipeline uses **two PharmCAT processes** (`RUN_PGX_PHARMCAT` + `FINALIZE_PGX_JSON`) and one Aldy process (`RUN_ALDY_CYP2D6`).
- `RUN_PGX_PHARMCAT` runs inside the PharmCAT image (label `pgx`); `FINALIZE_PGX_JSON` runs in `python:3.12-bookworm` (label `pgx_finalize`) because the PharmCAT image may not expose `python3` on `PATH`.
- `pgx_finalize.py` builds `pgx_meta.json` and `pgx_result.json` and embeds PharmCAT JSON outputs when present.
- **Nextflow 25+**: `Channel.into` was removed; the workflow uses **`multiMap`** to fan out the annotated VCF to PGx, IGV, and summary.

### Overriding container images

PharmCAT and Aldy images are pinned in `nextflow.config`. To override, add a small config file or extend the main config:

```groovy
process {
  withLabel: 'pgx' {
    container = 'pgkb/pharmcat:3.2.0'
  }
  withLabel: 'aldy' {
    container = 'gx-exome/aldy:4.8.3'
  }
}
```

`params.pgx_container` in `main.nf` is still passed into the PGx process script for metadata in `pgx_meta.json`; keep it in sync with the image you run.
