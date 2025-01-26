#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="FLAMES runner",
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
        "--output_dir", required=False, default=None, help="directory for output files"
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    num_threads = args.ncpu
    output_dir = args.output_dir

    if output_dir is not None:
        if not (os.path.exists):
            os.makedirs(output_dir)

    datasetName = os.path.basename(output_prefix)
    talonPrefix = datasetName
    dataType = "pacbio"

    cmd = f"talon_label_reads --f {bam_file} --t 1 --o {talonPrefix} --g {genome_fasta}"
    run_cmd(cmd)

    cmd = f"samtools calmd -@ {num_threads} --reference {genome_fasta} {talonPrefix}_labeled.sam > {talonPrefix}_labeled.md.sam"
    run_cmd(cmd)

    cmd = f"talon_initialize_database --f {gtf_file} --g {datasetName} --a {datasetName} --o {datasetName}"
    run_cmd(cmd)

    with open(f"{talonPrefix}.csv", "wt") as ofh:
        print(
            f"{datasetName},{datasetName},{dataType},{talonPrefix}_labeled.md.sam",
            file=ofh,
        )

    cmd = f"talon --build {datasetName} --db {datasetName}.db --o {talonPrefix}_raw --f {talonPrefix}.csv --threads {num_threads}"
    run_cmd(cmd)

    cmd = f"talon_filter_transcripts --db {datasetName}.db -a {datasetName} --datasets {datasetName} --o {talonPrefix}_filter"
    run_cmd(cmd)

    cmd = f"talon_create_GTF --build {datasetName} --db {datasetName}.db -a {datasetName} --o {talonPrefix} --whitelist {talonPrefix}_filter"
    run_cmd(cmd)

    cmd = f"talon_abundance --db {datasetName}.db --whitelist {talonPrefix}_filter --o {talonPrefix}_Quant --build {datasetName} -a {datasetName}"
    run_cmd(cmd)

    if output_dir is not None:
        cmd = f"cp {datasetName}* {output_dir}/"
        run_cmd(cmd)

    sys.exit(0)


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
