#!/usr/bin/env python3

import sys
import os
import pandas as pd
import argparse


def generate_metrics(dataset, tracking_file, ordered_tools):

    metrics_by_tool = {
        tool: {"TP": 0, "FP": 0, "FN": 0, "Precision": 0, "Recall": 0, "F1": 0}
        for tool in ordered_tools
    }

    with open(tracking_file, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            ref = columns[4]
            tool_predictions = columns[5:]

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


def main():

    parser = argparse.ArgumentParser(
        description="eval TP, FP, FN from tracking file of de novo bmarking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--tracking", type=str, required=True, help="tracking file")

    parser.add_argument("--dataset", type=str, required=True, help="name of dataset")

    parser.add_argument(
        "--ordered_tools",
        type=str,
        required=True,
        nargs="+",
        help="tools ordered for comparison",
    )

    args = parser.parse_args()

    tracking_file = args.tracking
    ordered_tools = args.ordered_tools
    dataset = args.dataset

    generate_metrics(dataset, tracking_file, ordered_tools)

    print("Completed metrics generation.")

    sys.exit(0)


if __name__ == "__main__":
    main()
