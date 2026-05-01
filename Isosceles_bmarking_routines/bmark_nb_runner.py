#!/usr/bin/env python3

import sys, os
import argparse
import logging
import glob
import shutil
import tempfile
import json
import subprocess

import yaml

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
    "quant_only_no_truthset": os.path.join(
        notebook_template_dir, "template.quant-only-no-truthset.ipynb"
    ),
}

QuantParser.FLAMES_gff3_converter = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "misc", "FLAMES_gff3_to_gtf_converter.py")
)


REGISTRY_FILENAME = "tool_registry.yaml"


def find_default_registry(start_dir: str) -> str:
    """Walk up from start_dir looking for tool_registry.yaml."""
    d = os.path.abspath(start_dir)
    while True:
        candidate = os.path.join(d, REGISTRY_FILENAME)
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(d)
        if parent == d:
            raise FileNotFoundError(
                f"Could not find {REGISTRY_FILENAME} walking up from {start_dir}. "
                "Pass --registry explicitly."
            )
        d = parent


def load_report_set(path: str) -> dict:
    """Optional report-set YAML: {include: [name, ...], overrides: {name: {display: bool, venn: bool}}}."""
    with open(path, "rt") as fh:
        data = yaml.safe_load(fh) or {}
    include = data.get("include")
    if include is not None and not isinstance(include, list):
        raise ValueError(f"{path}: 'include' must be a list of registry entry names")
    overrides = data.get("overrides") or {}
    if not isinstance(overrides, dict):
        raise ValueError(f"{path}: 'overrides' must be a mapping of name -> {{display, venn}}")
    return {"include": include, "overrides": overrides}


def execute_notebook_isolated(notebook_to_run, notebook_output, inputs_dict):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(inputs_dict, f)
        params_file = f.name

    with tempfile.TemporaryDirectory() as temp_dir:
        runtime_dir = os.path.join(temp_dir, "jupyter_runtime")
        os.makedirs(runtime_dir, exist_ok=True)

        try:
            cmd = [
                "papermill",
                "-f",
                params_file,
                "--kernel",
                "python3",
                "--progress-bar",
                "--log-level",
                "INFO",
                notebook_to_run,
                notebook_output,
            ]
            cmdstr = " ".join(cmd)
            logger.info(cmdstr)

            env = os.environ.copy()
            env["JUPYTER_RUNTIME_DIR"] = runtime_dir

            result = subprocess.run(
                cmdstr,
                text=True,
                env=env,
                timeout=3600,
                shell=True,
            )

            if result.returncode != 0:
                print(f"Papermill stderr: {result.stderr}")
                print(f"Papermill stdout: {result.stdout}")
                raise subprocess.CalledProcessError(
                    result.returncode, cmd, result.stdout, result.stderr
                )

            return result

        except subprocess.TimeoutExpired:
            print(f"Notebook execution timed out: {notebook_to_run}")
            raise
        except subprocess.CalledProcessError as e:
            print(f"Papermill execution failed: {e}")
            raise
        finally:
            if os.path.exists(params_file):
                os.unlink(params_file)


def main():

    parser = argparse.ArgumentParser(
        description="run benchmarking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--analysisType",
        type=str,
        required=True,
        choices=["quant_only", "ref_guided", "ref_free", "quant_only_no_truthset"],
        help="type of benchmarking to run.",
    )

    parser.add_argument("--truth_gtf", required=True, type=str, help="truth gtf file")
    parser.add_argument(
        "--truth_reduced_gtf",
        required=False,
        default=None,
        type=str,
        help="truth reduced gtf file - that gtf provided as a guide to ref-guided ID methods",
    )
    parser.add_argument(
        "--truth_quant", required=False, type=str, help="truth quant file"
    )
    parser.add_argument(
        "--Venn_mode",
        action="store_true",
        default=False,
        help="incorporate agreed-upon predicted isoforms into the truth set",
    )
    parser.add_argument(
        "--registry",
        required=False,
        default=None,
        type=str,
        help=(
            "Path to tool_registry.yaml. If omitted, walks up from the "
            "current directory looking for one."
        ),
    )
    parser.add_argument(
        "--report_set",
        required=False,
        default=None,
        type=str,
        help=(
            "Optional report-set YAML overriding which registry entries "
            "are included and how. Schema: {include: [names], "
            "overrides: {name: {display: bool, venn: bool}}}."
        ),
    )

    args = parser.parse_args()

    registry_path = args.registry or find_default_registry(os.getcwd())
    logger.info(f"Loading registry: {registry_path}")
    entries = QuantParser.load_registry(registry_path)

    if args.report_set:
        rs = load_report_set(args.report_set)
        active_entries = build_active_entries(entries, rs)
        logger.info(
            f"Applied report set {args.report_set}: "
            f"{len(active_entries)} entries active out of {len(entries)} registered"
        )
    else:
        active_entries = [e for e in entries if e["display"]]
        logger.info(
            f"No report set; {len(active_entries)} entries with display:true "
            f"out of {len(entries)} registered"
        )

    program_files = prep_files(active_entries)

    if not program_files:
        raise RuntimeError(
            "No tool outputs matched the active registry entries in raw_prog_results/"
        )

    inputs_dict = {
        "PYLIB_DIR": PYLIB_DIR,
        "REF_gtf_file": args.truth_gtf,
        "REF_quant_file": args.truth_quant,
        "program_files": program_files,
        "IGNORE_NONUNIQUE_NONREF": bool(args.Venn_mode),
    }
    if args.truth_reduced_gtf is not None:
        inputs_dict["REF_reduced_gtf_file"] = args.truth_reduced_gtf

    template_notebook = analysisType_to_notebook[args.analysisType]

    notebook_to_run = os.path.abspath("template." + args.analysisType + ".ipynb")
    shutil.copyfile(template_notebook, notebook_to_run)

    notebook_output = os.path.abspath(args.analysisType + ".ipynb")

    execute_notebook_isolated(notebook_to_run, notebook_output, inputs_dict)

    print("Done.")
    sys.exit(0)


def build_active_entries(entries, report_set):
    """Apply a report_set to the registry: keep only entries in `include`,
    apply per-name `overrides` to display/venn."""
    by_name = {e["name"]: e for e in entries}

    include = report_set["include"]
    if include is not None:
        unknown = [n for n in include if n not in by_name]
        if unknown:
            raise ValueError(
                f"report_set 'include' references unknown registry entries: {unknown}"
            )
        active = [dict(by_name[n]) for n in include]
    else:
        active = [dict(e) for e in entries if e["display"]]

    overrides = report_set["overrides"]
    if overrides:
        active_by_name = {e["name"]: e for e in active}
        for name, ov in overrides.items():
            if name not in active_by_name:
                raise ValueError(
                    f"report_set 'overrides' references entry '{name}' "
                    f"that is not in the active set"
                )
            for k, v in ov.items():
                if k not in ("display", "venn"):
                    raise ValueError(
                        f"report_set override for '{name}': only display/venn "
                        f"are settable, got '{k}'"
                    )
                active_by_name[name][k] = v

    return active


def prep_files(active_entries):
    """Glob raw_prog_results/* and route every file through QuantParser
    against active_entries. Build {name: {quant, gtf, color, family, venn}}."""

    filenames = sorted(glob.glob("raw_prog_results/*"))
    if not filenames:
        raise RuntimeError("Error, not finding result files at raw_prog_results/")

    program_files: dict = {}
    for filename in filenames:
        if not os.path.isfile(filename):
            continue
        kind, entry, processed_path = QuantParser.process_file(filename, active_entries)
        if kind is None:
            continue

        rec = program_files.setdefault(
            entry["name"],
            {
                "quant": None,
                "gtf": None,
                "color": entry["color"],
                "family": entry["family"],
                "venn": bool(entry["venn"]),
            },
        )
        if kind == "quant":
            rec["quant"] = processed_path
        elif kind == "gtf":
            # only propagate the GTF to the notebook if the entry asks
            # to use its own gtf for intron derivation
            if entry["gtf_source"] == "own":
                rec["gtf"] = processed_path
            # else: file is recognized (no deposit warning) but unused

    # Drop any entry that didn't yield a quant file -- a tool can't be
    # plotted without a quant.
    incomplete = [n for n, r in program_files.items() if r["quant"] is None]
    for n in incomplete:
        logger.warning(
            f"entry {n}: matched gtf but no quant in raw_prog_results/; dropping"
        )
        del program_files[n]

    # Warn about own-gtf entries that produced a quant but no gtf; those
    # will silently fall back to REF_gtf at the notebook layer, which is
    # probably wrong for tools with private transcript IDs.
    by_name = {e["name"]: e for e in active_entries}
    for n, r in program_files.items():
        if by_name[n]["gtf_source"] == "own" and r["gtf"] is None:
            logger.warning(
                f"entry {n}: gtf_source is 'own' but no GTF was matched; "
                f"REF_gtf will be substituted, transcript IDs may not align"
            )

    return program_files


if __name__ == "__main__":
    main()
