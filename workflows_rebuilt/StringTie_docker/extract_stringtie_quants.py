#!/usr/bin/env python3

import sys, os, re

usage = "usage: stringtie.gtf output_file.quants.tsv\n\n"

if len(sys.argv) < 3:
    exit(usage)

stringtie_gtf = sys.argv[1]
quants_output_file = sys.argv[2]


def get_key_val_pairs(info):
    key_val_dict = dict()

    key_vals = info.split(";")
    for key_val in key_vals:
        if key_val == "":
            continue
        key_val = key_val.strip()
        key_val = key_val.replace('"', "")
        key, val = key_val.split(" ")
        key_val_dict[key] = val

    return key_val_dict


# print header
print("\t".join(["gene_id", "transcript_id", "cov", "FPKM", "TPM"]))

with open(stringtie_gtf, "rt") as fh:
    with open(quants_output_file, "wt") as ofh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 9:
                continue

            if vals[2] != "transcript":
                continue

            info_line = vals[8]
            info_dict = get_key_val_pairs(info_line)

            print(
                "\t".join(
                    [
                        info_dict[x]
                        for x in ("gene_id", "transcript_id", "cov", "FPKM", "TPM")
                    ]
                ),
                file=ofh,
            )


exit(0)
