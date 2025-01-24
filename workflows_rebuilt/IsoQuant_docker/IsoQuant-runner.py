#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="isoquant runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output_prefix", required=True, help="output file prefix")
    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")
    parser.add_argument(
        "--gtf", type=str, default=None, required=False, help="input gtf file"
    )
    parser.add_argument(
        "--bam", type=str, required=True, help="input bam alignment file"
    )
    parser.add_argument(
        "--data_type",
        type=str,
        required=True,
        choices=["assembly", "pacbio_ccs", "nanopore"],
    )
    parser.add_argument(
        "--ncpu", type=int, required=False, default=4, help="num threads"
    )
    parser.add_argument(
        "--quant_only",
        action="store_true",
        default=False,
        help="only perform quantification, no isoform discovery.",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    num_threads = args.ncpu
    quant_only_flag = args.quant_only

    if quant_only_flag and gtf_file is None:
        raise RuntimeError("Error, must specify --gtf if --quant_only set")

    ## begin

    cmd = " ".join(
        [
            "/usr/local/src/IsoQuant/isoquant.py",
            f"--reference {genome_fasta}",
            f"--bam {input_bam}",
            f"--data_type {data_type}",
            f"--threads {num_threads}",
        ]
    )

    if gtf_file is not None:
        cmd += f" --genedb {gtf_file} "

    if quant_only_flag:
        cmd += f" --no_model_reconstruction"

    run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
