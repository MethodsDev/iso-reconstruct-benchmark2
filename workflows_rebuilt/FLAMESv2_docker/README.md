# FLAMESv2 Docker Integration

This directory contains a fresh FLAMES v2 integration for the isoform benchmark workflows.

## Scope

- Uses upstream FLAMES v2 (`mritchielab/FLAMES`) as new software, not an in-place upgrade of legacy FLAMES.
- Keeps a benchmark-compatible interface and outputs.
- Pairs with workflow file: `../FLAMESv2.wdl`.

## Components

- `Dockerfile`: runtime image for FLAMES v2 execution.
- `FLAMESv2-runner.Rscript`: CLI wrapper used by WDL/task execution.
- `build_docker.sh`: builds versioned and `latest` images and runs a help check.
- `push_docker.sh`: pushes versioned and `latest` images.
- `testing/Makefile`: local smoke test command.

## Current Design

### Input interface (benchmark compatibility)

The wrapper accepts BAM-first inputs:

- `--bam`
- `--genome`
- `--gtf`
- `--data_type`
- `--ncpu`
- `--quant_only` (optional)

Even though FLAMES `BulkPipeline()` is FASTQ-first, the runner converts BAM to FASTQ internally via:

- `samtools bam2fq`

### FLAMES pipeline path

The runner creates a FLAMES config and runs:

- `BulkPipeline(...)`
- `run_FLAMES(...)`

Key config decisions currently enforced:

- `pipeline_parameters.oarfish_quantification = FALSE`
  - Forces non-Oarfish FLAMES quantification (`transcript_count.csv.gz`) for stable benchmark output mapping.
- `pipeline_parameters.do_isoform_identification = !quant_only`
  - `quant_only=false`: reference-guided isoform ID + quantification.
  - `quant_only=true`: skip isoform discovery, quantify against reference-guided transcript model.

### Output contract

Runner emits benchmark-standard files:

- `<sample>.FLAMES.gff3`
- `<sample>.FLAMES.counts.tsv`

Output mapping behavior:

- GFF3:
  - prefers FLAMES-produced `isoform_annotated*.gff3/gtf`.
  - in `quant_only` mode, falls back to input GTF if no isoform annotation is generated.
- Counts:
  - primary path: parse `transcript_count.csv.gz` and collapse sample columns to `gene_id, transcript_id, count`.
  - fallback path: parse `*.quant` if needed.

## Modes

### 1) Reference-guided isoform identification + quantification

- WDL: `quant_only = false`
- Runner: no `--quant_only`

This runs genome alignment, isoform identification, read realignment, transcript quantification.

### 2) Quant-only

- WDL: `quant_only = true`
- Runner: `--quant_only`

This skips isoform identification and performs transcript quantification from reference annotation.

## WDL integration

`../FLAMESv2.wdl` includes:

- task input `Boolean quant_only = false`
- command flag wiring to runner `--quant_only` when enabled
- output files:
  - `~{sample_id}.FLAMES.gff3`
  - `~{sample_id}.FLAMES.counts.tsv`

## Build and test

From `FLAMESv2_docker`:

```bash
./build_docker.sh
```

Smoke test:

```bash
cd testing
make test
```

## Known pitfalls / future work

1. FLAMES config template uses `type = "sc_3end"` as a base because FLAMES currently lacks a dedicated `type = "bulk"` preset in `create_config()`.
   - We override key fields for bulk benchmarking behavior.
2. If upstream FLAMES changes output file names or config structure, update `FLAMESv2-runner.Rscript` mapping logic first.
3. If future workflows want Oarfish quantification, add an explicit switch and separate output adapter path.
4. If multi-sample bulk is required later, update runner to pass multiple FASTQs/sample names and revise count-collapsing logic.

## Versioning notes

- Image name: `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2`
- Version source: `VERSION.txt`
- Current scaffold version: `0.1.0`
