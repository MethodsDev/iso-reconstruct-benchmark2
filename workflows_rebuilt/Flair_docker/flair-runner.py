#!/usr/bin/env python

import sys, os, re, shutil
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
    parser.add_argument(
        "--flair_version_tag",
        type=str,
        required=False,
        default="v3.0.0",
        help="version token appended to final output filenames",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    num_threads = args.ncpu
    quant_only_flag = args.quant_only
    flair_version_tag = args.flair_version_tag
    versioned_output_prefix = f"{output_prefix}.flair-{flair_version_tag}"

    ## begin
    flair_cmd = get_flair_cmd()

    run_cmd(f"samtools bam2fq {bam_file} > temp.fastq")
    generate_reads_bed(bam_file=bam_file, output_prefix=output_prefix)

    if quant_only_flag:
        gtf_for_quant = gtf_file
    else:
        # perform isoform ID ref-guided

        cmd = " ".join(
            [
                f"{flair_cmd} correct",
                f"-q {output_prefix}-reads.bed",
                f"-f {gtf_file}",
                f"-o {output_prefix}",
                f"-t {num_threads}",
            ]
        )
        run_cmd(cmd)

        cmd = " ".join(
            [
                f"{flair_cmd} collapse",
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

        run_cmd(
            f"cp {output_prefix}.isoforms.gtf {versioned_output_prefix}.flair.gtf"
        )

        gtf_for_quant = f"{versioned_output_prefix}.flair.gtf"

    # run quant
    run_cmd(f"gffread {gtf_for_quant} -g {genome_fasta} -w gtf_transcriptome.fa")

    manifest_filename = "flair_manifest.tsv"
    with open(manifest_filename, "wt") as ofh:
        print("\t".join(["flair", "condition1", "batch1", "temp.fastq"]), file=ofh)

    cmd = " ".join(
        [
            f"{flair_cmd} quantify",
            f"-r {manifest_filename}",
            "-i gtf_transcriptome.fa",
            f"-o {versioned_output_prefix}",
            f"-t {num_threads}",
        ]
    )
    run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


def generate_reads_bed(bam_file, output_prefix):
    reads_bed_file = f"{output_prefix}-reads.bed"
    try:
        run_cmd(f"bam2Bed12 -i {bam_file} > {reads_bed_file}")
        return
    except subprocess.CalledProcessError:
        # FLAIR 3.0.0 conda packaging can ship a broken bam2Bed12 entrypoint.
        # Fallback to bedtools to generate BED12 from BAM.
        print(
            "WARNING: bam2Bed12 failed; falling back to bedtools bamtobed -bed12",
            file=sys.stderr,
        )
        run_cmd(f"bedtools bamtobed -bed12 -i {bam_file} > {reads_bed_file}")
        return


def get_flair_cmd():
    for cmd in ("flair", "flair.py"):
        if shutil.which(cmd):
            return cmd
    raise RuntimeError("Could not find FLAIR executable (`flair` or `flair.py`) in PATH")


if __name__ == "__main__":
    main()
