#!/usr/bin/env python

import sys, os, re
import subprocess
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="__add_descr__",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output_prefix", required=True, help="output file prefix")
    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")
    parser.add_argument("--gtf", type=str, required=True, help="input gtf file")
    parser.add_argument(
        "--bam", type=str, required=True, help="input bam alignment file"
    )
    parser.add_argument(
        "--ncpu", type=int, required=False, default=4, help="num threads"
    )
    parser.add_argument(
        "--sort_buffer_memGB",
        type=int,
        required=False,
        default=10,
        help="sort memory in GB",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    sort_buffer_memGB = args.sort_buffer_memGB
    num_threads = args.ncpu

    # prep inputs
    run_cmd(f"samtools view -o {output_prefix}.sam {bam_file}")

    samples_file = f"espresso_samples.tsv"
    with open(samples_file, "wt") as ofh:
        print("\t".join([f"{output_prefix}.sam", "espresso"]), file=ofh)

    # run espresso S
    cmd = " ".join(
        [
            "perl /opt/conda/envs/espresso_env/bin/ESPRESSO_S.pl",
            "--sort_buffer_size {}".format(sort_buffer_memGB),
            f"-L {samples_file}",
            f"-F {genome_fasta}",
            f"-A {gtf_file}",
            "-O . ",
            f"-T {num_threads}",
        ]
    )
    run_cmd(cmd)

    # run espresso C
    cmd = " ".join(
        [
            "perl /opt/conda/envs/espresso_env/bin/ESPRESSO_C.pl",
            f"--sort_buffer_size {sort_buffer_memGB}",
            "-I .",
            f"-F {genome_fasta}",
            "-X 0",
            f"-T {num_threads}",
        ]
    )
    run_cmd(cmd)

    # run espresso Q
    cmd = " ".join(
        [
            "perl /opt/conda/envs/espresso_env/bin/ESPRESSO_Q.pl",
            f"-L {samples_file}.updated",
            f"-A {gtf_file}",
            f"-T {num_threads}",
        ]
    )
    run_cmd(cmd)

    run_cmd(f"cp espresso_samples_N2_R0_updated.gtf {output_prefix}.espresso.gtf")
    run_cmd(
        f"cp espresso_samples_N2_R0_abundance.esp {output_prefix}.espresso.counts.tsv"
    )

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
