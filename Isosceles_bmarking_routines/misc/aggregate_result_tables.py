#!/usr/bin/env python3

import sys, os, re
import pandas as pd


def main():

    usage = "\n\n\tusage: {} tsv.files.list\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        exit(usage)

    all_data = None

    tsv_file_list = sys.argv[1]
    with open(tsv_file_list, "rt") as fh:
        for tsv_filename in fh:
            tsv_filename = tsv_filename.rstrip()
            if tsv_filename.split(".")[-1] == "csv":
                df = pd.read_csv(tsv_filename, sep=",")
            else:
                df = pd.read_csv(tsv_filename, sep="\t")
            data_name = tsv_filename.split("/")[-2]
            df[["data_name"]] = data_name
            all_data = pd.concat([all_data, df])

    all_data.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
