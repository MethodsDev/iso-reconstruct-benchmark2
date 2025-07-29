#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
import glob
import pandas as pd


# Note, I suspect this is overcounting TPs and maybe FPs and also results in inconsistent truth sets - each method has a different effective truth set being applied based on the logic wihtin. Renameing as 'flawed.' until resolved.


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


TOOLDIR = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser(
        description="run benchmarking for reconstructions on real data examining known and novel isoforms reconstructed",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--truth_gtf",
        type=str,
        required=True,
        help="Isoforms defined as truth",
    )

    parser.add_argument(
        "--exclude_gtf",
        type=str,
        required=False,
        default=None,
        help="Isoforms to exclude from benchmarking (eg. those provided for ref-guided mode)",
    )

    parser.add_argument(
        "--dataset_name", type=str, required=True, help="name of dataset"
    )

    args = parser.parse_args()

    ref_gtf = os.path.abspath(args.truth_gtf)
    dataset_name = args.dataset_name
    exclude_gtf = args.exclude_gtf

    #

    tool_gtf_pairs = list()
    for filename in glob.glob("raw_prog_results/*gtf") + glob.glob(
        "raw_prog_results/*gff"
    ):

        if re.search("FLAMES.gff3", filename) is not None:
            # convert to gtf
            outputfile = "FLAMES.gtf"
            cmd = "/home/bhaas/projects/LRAA_project/LRAA_PAPER_Analyses/iso-reconstruct-benchmark2/misc/FLAMES_gff3_to_gtf_converter.py {} > {}".format(
                filename, outputfile
            )
            subprocess.check_cal(cmd, shell=True)
            tool_gtf_pairs.append([outputfile, None])

        elif re.search("bambu.gtf$", filename) is not None:
            tool_gtf_pairs.append(["Bambu.gtf", filename])
        elif re.search("Mandalor.*gtf$", filename) is not None:
            tool_gtf_pairs.append(["Mandalorion.gtf", filename])
        elif re.search("IsoQuant.gtf$", filename) is not None:
            tool_gtf_pairs.append(["IsoQuant.gtf", filename])
        elif re.search("LRAA.gtf$", filename) is not None:
            tool_gtf_pairs.append(["LRAA.gtf", filename])
        elif re.search("stringtie.gtf$", filename) is not None:
            tool_gtf_pairs.append(["StringTie.gtf", filename])
        elif re.search("IsoSeq.*gff$", filename) is not None:
            tool_gtf_pairs.append(["IsoSeq.gtf", filename])
        elif re.search("isosceles.*gtf$", filename) is not None:
            tool_gtf_pairs.append(["Isosceles.gtf", filename])
        elif re.search("flair.gtf$", filename) is not None:
            tool_gtf_pairs.append(["Flair.gtf", filename])
        elif re.search("talon.gtf$", filename) is not None:
            tool_gtf_pairs.append(["TALON.gtf", filename])
        elif re.search("espresso.gtf", filename) is not None:
            tool_gtf_pairs.append(["ESPRESSO.gtf", filename])

        else:
            raise RuntimeError("Error, not recognizing gtf file: {}".format(filename))

    # make symlinks
    tool_names = list()
    tool_gtfs = list()
    for target_fname, fname in tool_gtf_pairs:
        tool_gtfs.append(target_fname)
        if fname is not None:
            subprocess.check_call(
                "ln -sf {} {}".format(fname, target_fname), shell=True
            )
        tool_basename = re.sub("\\.(gtf|gff3|gff)$", "", target_fname, flags=re.I)
        tool_names.append(tool_basename)

    # run gffcompare
    cmd = f"gffcompare -o {dataset_name}.denovo {ref_gtf} " + " ".join(tool_gtfs)
    if exclude_gtf is not None:
        cmd += " " + exclude_gtf

    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

    tracking_file = f"{dataset_name}.denovo.tracking"
    if not os.path.exists(tracking_file):
        raise RuntimeError("Error, missing output file: " + tracking_file)

    generate_metrics(dataset_name, tracking_file, tool_names, exclude_gtf)

    sys.exit(0)


def generate_metrics(dataset, tracking_file, ordered_tools, exclude_gtf):

    metrics_by_tool = {
        tool: {"TP": 0, "FP": 0, "FN": 0, "Precision": 0, "Recall": 0, "F1": 0}
        for tool in ordered_tools
    }

    with open(tracking_file, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            ref = columns[4]
            tool_predictions = columns[5:]
            if exclude_gtf is not None:
                # skip entries that show up in the exclude category
                if tool_predictions[-1] != "-":
                    continue
                else:
                    # trim out the exclude column
                    tool_predictions = tool_predictions[:-1]

            non_dash_count = sum(
                1 for prediction in tool_predictions if prediction != "-"
            )

            for tool_index, prediction in enumerate(tool_predictions):
                tool_name = ordered_tools[tool_index]

                if ref != "-" and prediction != "-":
                    metrics_by_tool[tool_name]["TP"] += 1
                elif ref != "-" and non_dash_count >= 1 and prediction == "-":
                    metrics_by_tool[tool_name]["FN"] += 1
                elif ref == "-" and non_dash_count == 1 and prediction != "-":
                    metrics_by_tool[tool_name]["FP"] += 1

    for tool, metrics in metrics_by_tool.items():
        TP = metrics["TP"]
        FP = metrics["FP"]
        FN = metrics["FN"]
        Precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        Recall = TP / (TP + FN) if (TP + FN) > 0 else 0
        F1 = (
            2 * (Precision * Recall) / (Precision + Recall)
            if (Precision + Recall) > 0
            else 0
        )
        metrics["Precision"] = Precision
        metrics["Recall"] = Recall
        metrics["F1"] = F1

    results = []
    for tool_name, metrics in metrics_by_tool.items():
        results.append(
            {
                "Tool": tool_name,
                "TP": metrics["TP"],
                "FP": metrics["FP"],
                "FN": metrics["FN"],
                "Precision": metrics["Precision"],
                "Recall": metrics["Recall"],
                "F1": metrics["F1"],
            }
        )

    df = pd.DataFrame(results)

    df.to_csv(f"{dataset}.denovo_ID_summary.tsv", sep="\t", index=False)

    return


if __name__ == "__main__":
    main()
