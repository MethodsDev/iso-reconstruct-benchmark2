#!/usr/bin/env python

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

    run_cmd(f"samtools bam2fq ~{bam_file} > temp.fastq")

    cmd = " ".join(
        [
            "pbmm2 align",
            f"--num-threads {num_threads}",
            "--preset ISOSEQ",
            "--sort",
            f"{referenceGenomeFasta}",
            "temp.fastq",
            "pbmm_aligned.bam",
        ]
    )

    # Ref-free isoform ID (needed for ref-guided too)
    cmd = " ".join(
        [
            "isoseq3 collapse",
            "--do-not-collapse-extra-5exons",
            "pbmm_aligned.bam",
            f"{output_prefix}.IsoSeq.ref-free.ID.gff",
        ]
    )

    run_cmd(cmd)

    if gtf_file is not None:
        # ref-filtered isoform ID
        cmd = "pigeon prepare {output_prefix}.IsoSeq.ref-free.ID.gff"
        run_cmd(cmd)

        cmd = "pigeon prepare {gtf_file} {genome_fasta}"
        run_cmd(cmd)

        cmd = " ".join(
            [
                "pigeon classify",
                "~{output_prefix}.IsoSeq.ref-free.ID.gff",
                f"{gtf_file}",
                f"{referenceGenomeFasta}",
                "--fl pbmm_aligned.flnc_count.txt",
                "-d .",
            ]
        )
        run_cmd(cmd)

        sorted_gff_file = gff_file
        sorted_gff_file = sorted_gff_file.replace(".gff", ".sorted.gff")

        cmd = " ".join(
            [
                "pigeon classify",
                "{output_prefix}.IsoSeq.ref-free.ID.sorted.gff",
                f"{sorted_gff_file}",
                f"{referenceGenomeFasta}",
                "--fl pbmm_aligned.flnc_count.txt",
                "-d .",
            ]
        )
        run_cmd(cmd)

        cmd = " ".join(
            [
                "pigeon filter",
                "pbmm_aligned_classification.txt",
                f"--isoforms {output_prefix}.IsoSeq.ref-free.ID.sorted.gff",
            ]
        )
        run_cmd(cmd)

        run_cmd(
            f"cp {output_prefix}.IsoSeq.ref-free.ID.sorted.filtered_lite.gff ~{output_prefix}.IsoSeq.ref-filtered.ID.gff"
        )

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
