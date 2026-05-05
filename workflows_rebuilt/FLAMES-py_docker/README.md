# FLAMES-py Wrapper Notes

This wrapper runs the legacy `FLAMES-py` pipeline through `FLAMES-runner.py`.

## Strand-Specificity Handling

The wrapper currently switches between two bundled JSON config files:

- `flames.pbio.config.json`
- `flames.ont.config.json`

Despite those filenames, the only material difference between the two config files is the value of:

- `isoform_parameters.strand_specific`

Specifically:

- `flames.pbio.config.json` sets `"strand_specific": 1`
- `flames.ont.config.json` sets `"strand_specific": 0`

No other parameters differ between the two files.

## Practical Interpretation

Because the configs differ only by `strand_specific`, they should be selected according to the library strandedness you want FLAMES to assume, not according to sequencing technology name alone.

In other words:

- use the `strand_specific = 1` config for strand-specific libraries
- use the `strand_specific = 0` config for non-strand-specific libraries

The current filenames are historical/convenience labels, but they should not be interpreted as meaning that the wrapper is applying broader PacBio-specific or ONT-specific parameterization.

## Recommendation

If this wrapper is revised further, the config-selection logic should be documented and named in terms of strandedness rather than technology, since strandedness is the only actual distinction in the bundled configs.
