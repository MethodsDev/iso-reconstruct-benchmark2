#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="flair runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output_prefix", required=True, help="output file prefix")
    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")
    parser.add_argument("--gtf", type=str, required=False, help="input gtf file")
    parser.add_argument(
        "--bam", type=str, required=True, help="input bam alignment file"
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
        raise RuntimeError("need gtf file if quant_only mode")

    ## begin

    cmd = f"stringtie {bam_file} -L --ref {genome_fasta} -p {num_threads} -o {output_prefix}.stringtie.gtf"

    if gtf_file:
        cmd += f" -G {gtf_file}"

    if quant_only_flag:
        cmd += " -e"

    run_cmd(cmd)

    # get quants
    cmd = f"extract_stringtie_quants.py {output_prefix}.stringtie.gtf {output_prefix}.stringtie.quant.tsv"
    run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
