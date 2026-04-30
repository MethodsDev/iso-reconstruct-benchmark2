# Flair Docker Notes: `bam2Bed12` Regression and Fix

## Summary

When running Flair on Terra/Cromwell, the workflow failed at:

```text
bam2Bed12 ... 
ModuleNotFoundError: No module named 'flair.bam2Bed12'
```

This is not a Terra-specific issue. It is an upstream packaging regression in Flair `v3.0.0` (and current `master` at the time of investigation).

## Root Cause

- Flair ships a `bam2Bed12` executable/entrypoint.
- That entrypoint expects Python module `flair.bam2Bed12`.
- In Flair `v3.0.0`, `src/flair/bam2Bed12.py` is missing.
- So the entrypoint is present but points to code that is no longer included.

During source inspection, older tags (for example `v2.2.0`) still contained:

- `src/flair/bam2Bed12.py`
- `src/flair/samJuncs.py`

Those files are removed in `v3.0.0`, while `bin/bam2Bed12` remains.

## Fix Implemented Here

We fixed this in the Docker image directly, without relying on bedtools fallback:

1. Install Flair from source tag `v3.0.0`.
2. Replace `/opt/conda/envs/flair_env/bin/bam2Bed12` with a self-contained implementation that:
   - parses `-i/--input_bam`
   - reads BAM via `pysam`
   - infers strand from `ts` tag / orientation logic
   - emits BED12 records expected by `flair correct`
3. Keep required tools present in the environment (`samtools`, `minimap2`, etc.).
4. Runner now calls `bam2Bed12 -i ...` directly (no fallback path).

## Validation

Validated with local testing:

- `make test` in `workflows_rebuilt/Flair_docker/testing`
  - `ref_guided`: pass
  - `quant_only`: pass
- `make wdl_test` in `workflows_rebuilt/Flair_docker/testing`
  - `wdl_ref_guided`: pass
  - `wdl_quant_only`: pass (with optional GTF output behavior unchanged)

## Practical Impact

This restores expected Flair behavior in environments like Terra where the workflow depends on `bam2Bed12` for generating BED12 input to `flair correct`.
