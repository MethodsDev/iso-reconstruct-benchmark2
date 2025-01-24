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
        "--data_type",
        type=str,
        required=True,
        choices=["pacbio_ccs", "nanopore"],
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    genome_fasta = args.genome
    gtf_file = args.gtf
    bam_file = args.bam
    num_threads = args.ncpu
    data_type = args.data_type

    json_dir = os.path.dirname(os.path.abspath(__file__))
    config_json = (
        os.path.join(json_dir, "flames.pbio.config.json")
        if data_type == "pacbio_ccs"
        else os.path.join(json_dir, "flames.ont.config.json")
    )

    os.makedirs("temp_fq_dir")
    cmd = f"samtools bam2fq {bam_file} > temp_fq_dir/temp.fastq"
    run_cmd(cmd)

    cmd = " ".join(
        [
            "python3 /usr/local/src/FLAMES/python/bulk_long_pipeline.py",
            f"--gff3 {gtf_file}",
            f"--genomefa {genome_fasta}",
            "--fq_dir temp_fq_dir",
            f"--inbam {bam_file}",
            "--outdir . ",
            f"--config_file {config_json}",
        ]
    )

    run_cmd(cmd)

    output_gff3_filename = output_prefix + ".FLAMES.gff3"
    run_cmd(f"cp isoform_annotated.gff3 {output_gff3_filename}")

    output_counts_filename = output_prefix + ".FLAMES.counts.tsv"

    # make a counts file.
    with open(output_gff3_filename) as fh:
        with open(output_counts_filename, "wt") as ofh:
            print("\t".join(["gene_id", "transcript_id", "count"]), file=ofh)
            for line in fh:
                line = line.rstrip()
                vals = line.split("\t")
                if len(vals) < 9:
                    continue
                if vals[2] == "transcript":
                    # ID=transcript:SIRV706;transcript_id=SIRV706;Parent=gene:SIRV7;support_count=6;source=known
                    info = vals[8]
                    info_dict = get_key_val_dict(info)
                    gene_id, transcript_id, count = (
                        info_dict["Parent"],
                        info_dict["transcript_id"],
                        info_dict["support_count"],
                    )
                    print("\t".join([gene_id, transcript_id, count]), file=ofh)

    sys.exit(0)


def get_key_val_dict(info):
    key_val_pairs = info.split(";")
    ret_dict = dict()

    for key_val_pair in key_val_pairs:
        key, val = key_val_pair.split("=")
        ret_dict[key] = val

    return ret_dict


def run_cmd(cmd):

    print("CMD: " + cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

    return


if __name__ == "__main__":
    main()
