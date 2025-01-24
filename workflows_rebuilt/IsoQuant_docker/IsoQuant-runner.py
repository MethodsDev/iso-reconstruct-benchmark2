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
    data_type = args.data_type

    if quant_only_flag and gtf_file is None:
        raise RuntimeError("Error, must specify --gtf if --quant_only set")

    ## begin

    output_dir = f"{output_prefix}.isoquant.outdir"

    cmd = " ".join(
        [
            "/usr/local/src/IsoQuant/isoquant.py",
            f"--reference {genome_fasta}",
            f"--bam {bam_file}",
            f"--data_type {data_type}",
            f"--threads {num_threads}",
            f"--output {output_dir}",
        ]
    )

    if gtf_file is not None:
        cmd += f" --genedb {gtf_file} "

    if quant_only_flag:
        cmd += f" --no_model_construction"

    run_cmd(cmd)

    isoquant_models_gtf = f"{output_dir}/OUT/OUT.transcript_models.gtf"
    if os.path.exists(isoquant_models_gtf):
        run_cmd(f"cp {isoquant_models_gtf} {output_prefix}.IsoQuant.gtf")

    if quant_only_flag:
        counts_file = f"{output_dir}/OUT/OUT.transcript_counts.tsv"
    else:
        counts_file = f"{output_dir}/OUT/OUT.transcript_model_counts.tsv"

    run_cmd(f"cp {counts_file} {output_prefix}.IsoQuant.counts.tsv")

    run_cmd(f"ln -sf {output_dir}  isoquant_output_dir")

    print("done")

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
