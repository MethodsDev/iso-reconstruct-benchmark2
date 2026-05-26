# ESPRESSO integration

This directory contains the Docker packaging and Python wrapper used by
`../ESPRESSO.wdl` to run ESPRESSO for isoform reconstruction.

The container installs ESPRESSO 1.6.0 from Bioconda and copies
`ESPRESSO-runner.py` into `/usr/local/bin`. `VERSION.txt` is the wrapper image
version used by `build_docker.sh` and `push_docker.sh`; it is not the upstream
ESPRESSO version.

## Files

- `Dockerfile`: builds the runtime image with the `espresso_env` conda
  environment, ESPRESSO 1.6.0, `samtools`, and the local wrapper.
- `ESPRESSO-runner.py`: normalizes inputs and runs the ESPRESSO S/C/Q steps.
- `../ESPRESSO.wdl`: exposes the wrapper as a WDL task/workflow.
- `testing/`: small SIRV test data, miniwdl input JSON files, and Makefile
  targets for ref-guided and annotation-free runs.

## Runner behavior

`ESPRESSO-runner.py` accepts:

```bash
ESPRESSO-runner.py \
  --output_prefix sample \
  --genome genome.fa \
  --bam alignments.bam \
  [--gtf annotation.gtf] \
  [--ncpu 4] \
  [--sort_buffer_memGB 10]
```

The wrapper performs these steps:

1. Converts the input BAM to `<output_prefix>.sam` with `samtools view`.
2. Writes `espresso_samples.tsv` with one sample named `espresso`.
3. Runs `ESPRESSO_S.pl`.
4. Runs `ESPRESSO_C.pl`.
5. Runs `ESPRESSO_Q.pl`.
6. Copies ESPRESSO's final outputs to stable benchmark filenames:
   - `<output_prefix>.espresso.v1.6.0.gtf`
   - `<output_prefix>.espresso.v1.6.0.counts.tsv`

## Execution modes

The same runner supports two modes.

### Ref-guided

When `--gtf` is provided, the wrapper passes `-A <annotation.gtf>` to both
`ESPRESSO_S.pl` and `ESPRESSO_Q.pl`.

Example:

```bash
ESPRESSO-runner.py \
  --output_prefix sample.ref_guided \
  --genome genome.fa \
  --gtf annotation.gtf \
  --bam alignments.bam
```

### Annotation-free / denovo

When `--gtf` is omitted, the wrapper does not pass `-A` to ESPRESSO. For
`ESPRESSO_S.pl`, it also adds `--alignment_read_groups`, which ESPRESSO 1.6.0
requires when no annotation is supplied.

Example:

```bash
ESPRESSO-runner.py \
  --output_prefix sample.denovo \
  --genome genome.fa \
  --bam alignments.bam
```

In denovo mode, ESPRESSO output transcript IDs are still generated, but gene
assignments that would normally come from the annotation may be `NA`.

## WDL integration

`../ESPRESSO.wdl` defines `espressoWorkflow`, which wraps `espressoTask`.

Required workflow inputs:

- `sample_id`
- `inputBAM`
- `inputBAMIndex`
- `referenceGenomeFasta`
- `referenceGenomeIndex`

Optional workflow inputs:

- `referenceAnnotationGTF`: enables ref-guided mode when supplied; omit it for
  denovo mode.
- `docker`: Docker image tag to run. Defaults to
  `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest`.

Workflow outputs:

- `espresso_gtf`
- `espresso_counts`

The WDL command conditionally includes the annotation:

```wdl
~{"--gtf " + referenceAnnotationGTF}
```

If `referenceAnnotationGTF` is undefined, miniwdl/Cromwell render this as an
empty string and the runner executes in denovo mode.

## Build and publish

From this directory:

```bash
./build_docker.sh
./push_docker.sh
```

The scripts tag and push both:

- `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:<VERSION.txt>`
- `us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest`

After changes to `ESPRESSO-runner.py`, rebuild and push the image before
running WDLs that use the default published Docker tag.

## Tests

Run direct Docker smoke tests from `testing/`:

```bash
make ref_guided
make denovo
```

Run WDL smoke tests:

```bash
make wdl_ref_guided
make wdl_denovo
```

The WDL test inputs are:

- `testing/wdl_inputs/ref_guided.inputs.json`
- `testing/wdl_inputs/denovo.inputs.json`

Use `make clean` to remove local test outputs.

## Notes

- ESPRESSO itself is installed from Bioconda as version 1.6.0.
- The benchmark output filenames use the hard-coded token `v1.6.0` in
  `ESPRESSO-runner.py`.
- `inputBAMIndex` and `referenceGenomeIndex` are declared in the WDL to stage
  index files with the BAM and FASTA. The Python runner uses the BAM and FASTA
  paths directly.
