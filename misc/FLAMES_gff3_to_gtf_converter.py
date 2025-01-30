#!/usr/bin/env python3

import sys, os, re
from collections import defaultdict
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    usage = "\n\tusage: {} flames.gff3 > flames.gtf\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    flames_gff3_file = sys.argv[1]

    transcript_info = defaultdict(dict)

    transcript_id_to_gene_id = dict()

    logger.info(f"-parsing {flames_gff3_file}")

    with open(flames_gff3_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 9:
                continue

            chrom = vals[0]
            feature_type = vals[2]

            if feature_type not in ("exon", "transcript"):
                continue

            lend = vals[3]
            rend = vals[4]

            strand = vals[6]

            info = vals[8]

            if feature_type == "transcript":
                m = re.search("ID=transcript:([^;]+);.*Parent=gene:([^;]+);", info)
                if m is not None:
                    transcript_id = m.group(1)
                    gene_id = m.group(2)
                    transcript_id_to_gene_id[transcript_id] = gene_id
                else:
                    raise RuntimeError("couldn't extract gene_id from {}".format(info))

                continue

            transcript_id = None

            m = re.search("Parent=transcript:([^;]+);", info)
            if m is not None:
                transcript_id = m.group(1)

            else:
                raise RuntimeError(
                    "couldn't extract transcript_id from {}".format(info)
                )

            if transcript_id not in transcript_info:

                transcript_info[transcript_id]["chrom"] = chrom
                transcript_info[transcript_id]["strand"] = strand
                transcript_info[transcript_id]["exons"] = list()

            transcript_info[transcript_id]["exons"].append([lend, rend])

    # write gtf
    logger.info("-writing gtf output")

    for transcript_id, info_dict in transcript_info.items():

        chrom = info_dict["chrom"]
        strand = info_dict["strand"]
        exons = info_dict["exons"]

        gene_id = transcript_id_to_gene_id[transcript_id]

        for exon in exons:
            lend, rend = exon

            print(
                "\t".join(
                    [
                        chrom,
                        "flames",
                        "exon",
                        lend,
                        rend,
                        ".",
                        strand,
                        ".",
                        f'transcript_id "{transcript_id}"; gene_id "{gene_id}";',
                    ]
                )
            )

    logger.info("Done")

    sys.exit(0)


if __name__ == "__main__":
    main()
