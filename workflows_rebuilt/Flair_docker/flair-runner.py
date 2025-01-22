#!/usr/bin/env python

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
    parser.add_argument("--gtf", type=str, required=True, help="input gtf file")
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

    ## begin

    run_cmd(f"samtools bam2fq {bam_file} > temp.fastq")
    run_cmd(f"bam2Bed12 -i {bam_file} > {output_prefix}-reads.bed")

    if quant_only_flag:
        gtf_for_quant = gtf_file
    else:
        # perform isoform ID ref-guided

        cmd = " ".join(
            [
                "flair correct",
                f"-q {output_prefix}-reads.bed",
                f"--genome {genome_fasta}",
                f"--gtf {gtf_file}",
                f"-o {output_prefix}",
                f"-t {num_threads}",
            ]
        )
        run_cmd(cmd)

        cmd = " ".join(
            [
                "flair collapse",
                f"-g {genome_fasta}",
                f"-f {gtf_file}",
                "-r temp.fastq",
                f"-q {output_prefix}_all_corrected.bed",
                f"-o {output_prefix}",
                f"-t {num_threads}",
            ]
        )
        # include opts?
        # --stringent --check_splice --generate_map --annotation_reliant generate

        run_cmd(cmd)

        run_cmd(f"cp {output_prefix}.isoforms.gtf {output_prefix}.flair.gtf")

        gtf_for_quant = f"{output_prefix}.flair.gtf"

    # run quant
    run_cmd(f"gffread {gtf_for_quant} -g {genome_fasta} -w gtf_transcriptome.fa")

    manifest_filename = "flair_manifest.tsv"
    with open(manifest_filename, "wt") as ofh:
        print("\t".join(["flair", "condition1", "batch1", "temp.fastq"]), file=ofh)

    cmd = " ".join(
        [
            "flair quantify",
            f"-r {manifest_filename}",
            "-i gtf_transcriptome.fa",
            f"-o {output_prefix}.flair",
            f"-t {num_threads}",
        ]
    )
    run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
