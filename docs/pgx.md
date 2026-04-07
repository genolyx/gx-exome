# Pharmacogenomics (PGx) add-on (PharmCAT)

The exome pipeline can run an **optional** PGx interpretation step after the final sample VCF is produced. It does **not** replace alignment or variant calling; it only consumes the filtered (and VEP-annotated, if enabled) VCF.

## When it runs

- Controlled by Nextflow parameter **`--skip_pgx`** (default: `true` → PGx **off**).
- To enable from `run_analysis.sh`: **`--enable-pgx`** (passes `--skip_pgx false` to Nextflow).
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

## Implementation notes

- Wrapper script: `bin/modules/pgx_finalize.py` builds `pgx_meta.json` and `pgx_result.json` and embeds PharmCAT JSON outputs when present.
- The PharmCAT image must include **Python 3** for the finalizer (included in `pgkb/pharmcat` images used by the team).
