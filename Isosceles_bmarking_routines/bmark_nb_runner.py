#!/usr/bin/env python3

import sys, os, re
import argparse
import logging
import papermill as pm
import glob
import shutil

PYLIB_DIR = os.path.join(os.path.dirname(__file__), "pylib")
sys.path.insert(0, PYLIB_DIR)
import QuantParser


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


notebook_template_dir = os.path.join(os.path.dirname(__file__), "template_notebooks")

analysisType_to_notebook = {
    "quant_only": os.path.join(notebook_template_dir, "template.quant-only.ipynb"),
    "ref_guided": os.path.join(notebook_template_dir, "template.ref_reduced.ipynb"),
    "ref_free": os.path.join(notebook_template_dir, "template.denovo.ipynb"),
}

QuantParser.FLAMES_gff3_converter = os.path.join(
    os.path.dirname(__file__), "../misc/FLAMES_gff3_to_gtf_converter.py"
)


def main():

    parser = argparse.ArgumentParser(
        description="run benchmarking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--analysisType",
        type=str,
        required=True,
        choices=["quant_only", "ref_guided", "ref_free"],
        help="type of benchmarking to run. Choices: [quant_only, ref_guided, ref_free]",
    )

    parser.add_argument(
        "--truth_gtf",
        required=True,
        type=str,
        help="truth gtf file",
    )

    parser.add_argument(
        "--truth_reduced_gtf",
        required=False,
        default=None,
        type=str,
        help="truth reduced gtf file - that gtf provided as a guide to ref-guided ID methods",
    )

    parser.add_argument(
        "--truth_quant", required=True, type=str, help="truth quant file"
    )

    parser.add_argument(
        "--Venn_mode",
        action="store_true",
        default=False,
        help="incorporate agreed-upon predicted isoforms into the truth set",
    )

    args = parser.parse_args()
    analysisType = args.analysisType
    truth_gtf = args.truth_gtf
    truth_reduced_gtf = args.truth_reduced_gtf
    truth_quant = args.truth_quant

    inputs_dict = prep_files()

    inputs_dict["PYLIB_DIR"] = PYLIB_DIR
    inputs_dict["REF_gtf_file"] = truth_gtf
    inputs_dict["REF_quant_file"] = truth_quant
    inputs_dict["REF_reduced_gtf_file"] = truth_reduced_gtf

    if args.Venn_mode:
        inputs_dict["IGNORE_NONUNIQUE_NONREF"] = True

    template_notebook = analysisType_to_notebook[analysisType]

    notebook_to_run = os.path.abspath("template." + analysisType + ".ipynb")
    shutil.copyfile(template_notebook, notebook_to_run)

    notebook_output = os.path.abspath(analysisType + ".ipynb")

    pm.execute_notebook(
        notebook_to_run,
        notebook_output,
        parameters=inputs_dict,
    )

    print("Done.")

    sys.exit(0)


def prep_files():

    method_type_to_filename = dict()

    filenames = glob.glob("raw_prog_results/*")

    if len(filenames) == 0:
        raise RuntimeError("Error, not finding result files at raw_prog_results/")

    for filename in filenames:
        token, filename = QuantParser.process_file(filename)
        if token is not None:
            method_type_to_filename[token] = filename

    return method_type_to_filename


if __name__ == "__main__":
    main()
