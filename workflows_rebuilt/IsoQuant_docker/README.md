# IsoQuant Docker Workflow

This directory contains the Docker image, runner script, and tests for running IsoQuant in the iso-reconstruct benchmark workflows.

## Version

- IsoQuant software version: `v3.13.0`
- Docker version token (`VERSION.txt`): `3.13.0-BI`
- Runner output version token default: `v3.13.0`

## Key Files

- `Dockerfile`: builds the IsoQuant container
- `IsoQuant-runner.py`: wrapper that runs IsoQuant and standardizes outputs
- `build_docker.sh`: builds `:3.13.0` and `:latest` image tags and verifies runner invocation
- `push_docker.sh`: pushes Docker tags
- `testing/Makefile`: mode-level test runs (`quant_only`, `ref_guided`, `ref_free`)

Related WDL update:
- `../IsoQuant.wdl`

## Build

From this directory:

```bash
./build_docker.sh
```

This builds:
- `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:3.13.0`
- `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest`

## Test

From `testing/`:

```bash
make clean
make test
```

Targets:
- `quant_only`
- `ref_guided`
- `ref_free`

## Output Naming

Runner outputs include a version token in the filename:

- counts: `<sample_id>.isoquant-v3.13.0.IsoQuant.counts.tsv`
- gtf (when produced): `<sample_id>.isoquant-v3.13.0.IsoQuant.gtf`

The version token can be overridden via:

```bash
--isoquant_version_tag <token>
```

## Quant File Selection by Mode (Validated)

IsoQuant `v3.13.0` output naming differs by mode, so `IsoQuant-runner.py` selects counts files with mode-aware fallbacks.

### 1) `--quant_only`
- Expected source file copied to final output:
  - `OUT/OUT.transcript_counts.tsv`
- GTF output: not produced (expected)

### 2) ref-guided (`--gtf` provided, no `--quant_only`)
- Observed source file copied to final output:
  - `OUT/OUT.discovered_transcript_counts.tsv`
- GTF output: produced from
  - `OUT/OUT.transcript_models.gtf`

### 3) ref-free (no `--gtf`, no `--quant_only`)
- Observed source file copied to final output:
  - `OUT/OUT.discovered_transcript_counts.tsv`
- GTF output: produced from
  - `OUT/OUT.transcript_models.gtf`

## Runner Counts Fallback Logic (non-quant-only)

For non-quant-only runs, the runner checks these candidates in order and uses the first existing file:

1. `OUT/OUT.transcript_model_counts.tsv`
2. `OUT/OUT.discovered_transcript_counts.tsv`
3. `OUT/OUT.transcript_counts.tsv`

If none exist, the runner raises an error.

## Notes

- The testing Makefile uses `docker run -i` (not `-it`) to work in non-interactive environments.
- The testing Makefile references the explicit image tag `:latest`.
