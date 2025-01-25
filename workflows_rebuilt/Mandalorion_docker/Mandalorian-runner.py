#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import shutil


def main():

    parser = argparse.ArgumentParser(
        description="Mandalorian runner",
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

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    num_threads = args.ncpu

    cmd = f"samtools bam2fq {bam_file} > temp.fastq"
    run_cmd(cmd)

    tmpdir = "__mandotmpdir"
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)

    os.makedirs(tmpdir)

    cmd = " ".join(
        [
            "python3 /usr/local/src/Mandalorion/Mando.py",
            f"-G {genome_fasta}",
            f"-f temp.fastq",
            f"-p {tmpdir}",
            f"-t {num_threads}",
        ]
    )

    if gtf_file is not None:
        cmd += f" -g {gtf_file}"

    run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
