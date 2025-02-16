#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import csv

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="extract TP and FP from gtf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--progname", type=str, required=True, help="name of the program"
    )

    parser.add_argument(
        "--class_assignments",
        type=str,
        required=True,
        help="file full_transcriptome_TPR_FDR_F1.class_assignments.tsv",
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="prog gtf filename",
    )

    args = parser.parse_args()

    progname = args.progname
    class_assignments_file = args.class_assignments
    gtf_filename = args.gtf

    logger.info(
        "Getting transcripts corresponding to TP and FP categories for {}".format(
            progname
        )
    )

    transcript_id_to_class = dict()
    with open(class_assignments_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["progname"] == progname:
                if row["class"] in ("TP", "FP"):
                    transcript_id_to_class[row["transcript_ids"]] = row["class"]

    ofh_TPs = open("LRAA-TPs.gtf", "wt")
    ofh_FPs = open("LRAA-FPs.gtf", "wt")

    with open(gtf_filename) as fh:
        for line in fh:
            if line[0] == "\n" or line[0] == "#":
                continue

            line = line.rstrip()
            m = re.search('transcript_id "([^"]+)"', line)
            if m is not None:
                transcript_id = m.group(1)
                if transcript_id in transcript_id_to_class:
                    class_type = transcript_id_to_class[transcript_id]
                    ofh = ofh_TPs if class_type == "TP" else ofh_FPs
                    print(line, file=ofh)
            else:
                raise RuntimeError("No transcript_id parsed from line: " + line)

    print("Done")
    sys.exit(0)


if __name__ == "__main__":
    main()
