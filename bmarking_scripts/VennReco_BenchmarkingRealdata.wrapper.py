#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
import glob

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


TOOLDIR = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser(
        description="run benchmarking for reconstructions on real data examining known and novel isoforms reconstructed",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--ref_expressed_or_kept_gtf",
        type=str,
        required=True,
        help="reference isoforms to compare - you decide (if reduced reference, leverage the kept set for comparison of novel)",
    )

    parser.add_argument(
        "--dataset_name", type=str, required=True, help="name of dataset"
    )

    args = parser.parse_args()

    ref_gtf = os.path.abspath(args.ref_expressed_or_kept_gtf)
    dataset_name = args.dataset_name

    tool_gtfs = list()
    for filename in glob.glob("raw_prog_results/*gtf") + glob.glob(
        "raw_prog_results/*gff"
    ):
        if re.search("bambu.gtf$", filename) is not None:
            tool_gtfs.append(["Bambu.gtf", filename])
        elif re.search("Mandalor.*gtf$", filename) is not None:
            tool_gtfs.append(["Mandalorion.gtf", filename])
        elif re.search("IsoQuant.gtf$", filename) is not None:
            tool_gtfs.append(["IsoQuant.gtf", filename])
        elif re.search("LRAA.gtf$", filename) is not None:
            tool_gtfs.append(["LRAA.gtf", filename])
        elif re.search("stringtie.gtf$", filename) is not None:
            tool_gtfs.append(["StringTie.gtf", filename])
        elif re.search("IsoSeq.*gff$", filename) is not None:
            tool_gtfs.append(["IsoSeq.gtf", filename])
        else:
            raise RuntimeError("Error, not recognizing gtf file: {}".format(filename))

    # make symlinks
    for target_fname, fname in tool_gtfs:
        subprocess.check_call("ln -sf {} {}".format(fname, target_fname), shell=True)

    # just get the target names
    tool_gtfs = [x[0] for x in tool_gtfs]

    # run analysis
    cmd = " ".join(
        [
            "python3",
            os.path.join(TOOLDIR, "VennReco_BenchmarkingRealdata.py"),
            "--ref_expressed_or_kept_gtf " + ref_gtf,
            "--reco_gtfs " + " ".join(tool_gtfs),
            "--dataset_name " + dataset_name,
        ]
    )
    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
    sys.exit(0)


if __name__ == "__main__":
    main()
