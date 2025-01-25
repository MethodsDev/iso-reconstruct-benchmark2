#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import shutil


def main():

    parser = argparse.ArgumentParser(
        description="Oarfish runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output_prefix", required=True, help="output file prefix")
    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")
    parser.add_argument("--gtf", type=str, required=True, help="input gtf file")
    parser.add_argument(
        "--fastq", type=str, required=False, help="input reads in fastq format"
    )

    parser.add_argument(
        "--ncpu", type=int, required=False, default=4, help="num threads"
    )

    parser.add_argument(
        "--mode",
        choices=["byAlignment", "byReads"],
        required=True,
        help="quant method to use",
    )

    parser.add_argument(
        "--seq_tech",
        choices=["ont-cdna", "ont-drna", "pac-bio", "pac-bio-hifi"],
        default="pac-bio",
        help="seq tech mode for byReads mode",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    fastq_file = args.fastq
    oarfish_mode = args.mode

    num_threads = args.ncpu

    seq_tech = args.seq_tech

    run_cmd(f"gffread {gtf_file} -g {genome_fasta} -w transcriptome.fa")

    if oarfish_mode == "byAlignment":

        cmd = " ".join(
            [
                "minimap2 -a ",
                f"-t {num_threads}",
                "transcriptome.fa",
                "tmp.fastq",
                f" | samtools sort -n --threads {num_threads} -o mapped.namesorted.bam",
            ]
        )
        run_cmd(cmd)

        cmd = " ".join(
            [
                "oarfish",
                "--alignments mapped.namesorted.bam",
                f"-o {output_prefix}.Oarfish.byAlignment",
                "--model-coverage",
                f"--threads {num_threads}",
            ]
        )
        run_cmd(cmd)

    else:

        # by reads directly
        cmd = " ".join(
            [
                "oarfish",
                f"--reads {fastq_file}",
                f"-o {output_prefix}.Oarfish.byReads",
                "--model-coverage",
                f"--threads {num_threads}",
                "--reference  transcriptome.fa",
                f"--seq-tech {seq_tech}",
            ]
        )

        run_cmd(cmd)

    print("done")

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
