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
        "--ncpu", type=int, required=False, default=4, help="num threads"
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = os.path.abspath(args.genome)
    ref_gtf_file = os.path.abspath(args.gtf) if args.gtf is not None else None
    bam_file = os.path.abspath(args.bam)
    num_threads = args.ncpu

    ## begin

    os.chdir(os.path.dirname(output_prefix))

    run_cmd(f"samtools bam2fq {bam_file} > temp.fastq")

    # run pbmm2
    cmd = " ".join(
        [
            "pbmm2 align",
            f"--num-threads {num_threads}",
            "--preset ISOSEQ",
            "--sort",
            f"{genome_fasta}",
            "temp.fastq",
            "pbmm_aligned.bam",
        ]
    )
    run_cmd(cmd)

    # Ref-free isoform ID (needed for ref-guided too)
    denovo_gff_file = f"{output_prefix}.IsoSeq.ref-free.ID.gff"
    cmd = " ".join(
        [
            "isoseq3 collapse",
            "--do-not-collapse-extra-5exons",
            "pbmm_aligned.bam",
            denovo_gff_file,
        ]
    )

    run_cmd(cmd)

    if ref_gtf_file is not None:

        #########################
        # ref-filtered isoform ID

        # prep de novo gff file
        cmd = f"pigeon prepare {denovo_gff_file}"
        run_cmd(cmd)

        sorted_denovo_gff_file = denovo_gff_file
        sorted_denovo_gff_file = sorted_denovo_gff_file.replace(".gff", ".sorted.gff")

        # prep ref gtf file
        cmd = f"pigeon prepare {ref_gtf_file} {genome_fasta}"
        run_cmd(cmd)

        sorted_ref_gtf_file = ref_gtf_file
        sorted_ref_gtf_file = sorted_ref_gtf_file.replace(".gtf", ".sorted.gtf")

        flnc_count_file = denovo_gff_file
        flnc_count_file = flnc_count_file.replace(".gff", ".flnc_count.txt")

        cmd = " ".join(
            [
                "pigeon classify",
                sorted_denovo_gff_file,
                sorted_ref_gtf_file,
                genome_fasta,
                f"--fl {flnc_count_file}",
                f"-o {os.path.basename(output_prefix)}",
            ]
        )
        run_cmd(cmd)

        cmd = " ".join(
            [
                "pigeon filter",
                f"{output_prefix}_classification.txt",
                f"--isoforms {sorted_denovo_gff_file}",
            ]
        )
        run_cmd(cmd)

        sorted_filtered_denovo_gff_file = sorted_denovo_gff_file
        sorted_filtered_denovo_gff_file = sorted_filtered_denovo_gff_file.replace(
            ".sorted.gff", ".sorted.filtered_lite.gff"
        )

        run_cmd(
            f"cp {sorted_filtered_denovo_gff_file} {output_prefix}.IsoSeq.ref-filtered.ID.gff"
        )

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
