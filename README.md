# Carrier Screening Pipeline

Nextflow 기반 Carrier Screening 분석 파이프라인입니다. WES(Whole Exome Sequencing) FASTQ 입력으로부터 변이 호출, CNV/SV 검출, 반복 확장 분석, 유사유전자 해석까지 수행하며, **service-daemon**과 연동하여 annotation, review, report 생성을 지원합니다.

## Architecture

```
┌─────────────┐     Submit      ┌──────────────────┐     Docker Run     ┌─────────────────────┐
│   Portal     │ ──────────────→│  service-daemon   │ ─────────────────→│  Pipeline Container  │
│  (Web UI)    │                │  (FastAPI)        │                   │  (Nextflow)          │
│              │     result     │                   │     output_dir    │                      │
│  - Review    │ ←──────────────│  - Queue Manager  │ ←─────────────────│  - Align & Sort      │
│  - Report    │     .json      │  - Annotator      │  vcf/qc/summary   │  - Variant Calling   │
│              │                │  - VCF Parser     │                   │  - CNV/SV Detection  │
│              │   Create       │  - Review Builder │                   │  - Repeat Expansion  │
│              │   Report       │  - Report Gen     │                   │  - Pseudogene Res.   │
│              │ ──────────────→│  - PDF Generator  │                   │  - Visual Evidence   │
│              │     .pdf       │                   │                   │  - Summary Report    │
│              │ ←──────────────│                   │                   │                      │
└─────────────┘                └──────────────────┘                   └─────────────────────┘
```

## Pipeline Workflow

파이프라인은 5개의 분석 트랙으로 구성됩니다.

| Track | Module | Description | Output |
|-------|--------|-------------|--------|
| Alignment | `align.nf` | BWA-MEM2 정렬, MarkDuplicates | BAM, duplicate_metrics |
| Variant Calling | `variant.nf` | GATK HaplotypeCaller + VariantFiltration | Filtered VCF |
| CNV/SV | `cnv.nf`, `sv.nf` | GATK gCNV (Cohort/Case), Manta | CNV VCF, SV VCF |
| Repeat | `repeat.nf` | ExpansionHunter (FMR1 등) | EH VCF, JSON, SVG |
| Pseudogene | `pseudogene.nf` | Paraphase (SMN1/2, GBA), SMAca | JSON, TXT |

추가적으로 `coverage.nf`에서 QC 메트릭을, `visualize.nf`에서 IGV 스냅샷을, `summary.nf`에서 통합 리포트를 생성합니다.

## Output Directory Structure

파이프라인 완료 시 `--output_dir`에 다음 구조로 결과가 복사됩니다.

```
output_dir/
├── vcf/                    # Filtered VCF (SNV/Indel)
├── cnv/                    # CNV segment VCFs
├── sv/                     # Structural variant VCFs
├── repeat/                 # Repeat expansion results (VCF, JSON, SVG)
├── pseudogene/             # Paraphase, SMAca results
├── summary/                # Summary & detailed reports (TXT)
├── snapshots/              # IGV snapshots (HTML)
├── qc/                     # QC metrics, coverage, duplicate metrics
├── pipeline_info/          # Nextflow trace, timeline, report
└── pipeline_complete.json  # Completion marker (for service-daemon)
```

`pipeline_complete.json`은 service-daemon이 파이프라인 완료를 감지하는 데 사용됩니다.

## Quick Start

### Prerequisites

- Nextflow >= 23.04
- Docker or Singularity
- Reference genome (GRCh38)
- Target BED file

### Local Execution

```bash
# 1. Clone repository
git clone https://github.com/genolyx/carrier-screening.git
cd carrier-screening

# 2. Run pipeline
nextflow run bin/main.nf \
  --fastq_dir /path/to/fastq/sample \
  --outdir results \
  --output_dir output \
  --sample_name Sample_001 \
  -profile docker
```

### Docker Compose (with service-daemon)

```bash
# 1. Configure environment
cp .env.example .env
# Edit .env with actual paths

# 2. Start all services
cd docker
docker compose up -d

# 3. Access services
# - Dashboard: http://localhost:5000
# - service-daemon API: http://localhost:8000
# - Test Portal: open service-daemon/portal/index.html
```

## service-daemon Integration

service-daemon은 파이프라인 결과를 받아 annotation, review, report 생성을 수행합니다. 자세한 내용은 [genolyx/service-daemon](https://github.com/genolyx/service-daemon) 레포를 참조하세요.

### Workflow

1. **Portal Submit** → service-daemon이 주문을 받아 Nextflow 파이프라인 실행
2. **Pipeline Complete** → `pipeline_complete.json` 생성, service-daemon이 감지
3. **Process Results** → VCF에 ClinVar/gnomAD/HPO annotation 적용, ACMG 분류, `result.json` 생성
4. **Portal Review** → 리뷰어가 변이 선택, 코멘트 작성
5. **Create Report** → 선택된 변이로 다국어 PDF 리포트 생성

### API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/order/carrier_screening/submit` | 분석 주문 제출 |
| GET | `/order/{id}/status` | 주문 상태 조회 |
| GET | `/order/{id}/result` | result.json 조회 (리뷰용) |
| POST | `/order/{id}/report` | 리포트 생성 |
| GET | `/orders` | 주문 목록 조회 |

## Directory Structure

```
carrier-screening/
├── bin/                    # Nextflow pipeline
│   ├── main.nf            # Main workflow
│   ├── nextflow.config    # Configuration
│   └── modules/           # Process modules (10 files)
├── docker/                # Docker configuration
│   ├── Dockerfile         # Pipeline container
│   ├── docker-compose.yml # Full stack (pipeline + service-daemon)
│   ├── entrypoint.sh      # Container entrypoint
│   └── supervisord.conf   # Process management
├── dashboard/             # Local analysis dashboard (Flask)
├── daemon/                # [DEPRECATED] Legacy daemon → use service-daemon
├── Genomics_Pipeline/     # [LEGACY] Reference code → migrated to service-daemon
├── data/                  # Reference data, BED files
├── docs/                  # Documentation
├── src/                   # Utility scripts
└── .env.example           # Environment configuration template
```

## Configuration

주요 Nextflow 파라미터는 `bin/nextflow.config`에서 설정합니다. Docker Compose 환경 변수는 `.env.example`을 참조하세요.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fastq_dir` | (required) | FASTQ 파일 디렉토리 |
| `ref_fasta` | `data/refs/GRCh38.fasta` | Reference genome |
| `backbone_bed` | Twist Exome 2.0 | Target BED file |
| `outdir` | `results` | Analysis directory (intermediate) |
| `output_dir` | (optional) | Output directory (final results) |
| `sample_name` | (optional) | Sample name for tracking |
| `gcnv_model` | (optional) | Pre-trained gCNV model path |
| `skip_cnv` | `false` | Skip CNV analysis |
