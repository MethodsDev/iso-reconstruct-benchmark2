#!/usr/bin/env python3

"""
Registry-driven quant-file parser.

Each entry in the tool registry (tool_registry.yaml) declares one
(tool, version) combination: a regex matched against files in
raw_prog_results/, plus the column indices needed to extract
(transcript_id, TPM) and the display config used downstream.

The parser walks entries in the order they appear in the registry and
returns on the first match -- so versioned patterns must be listed
before the corresponding less-specific -vNA pattern within a family.
"""

import os
import re
import subprocess
from typing import List, Optional, Tuple

import pandas as pd
import yaml


PROCESSED_DIR = "processed_prog_results"

# Set by bmark_nb_runner.py before process_file() is called for any
# entry with gtf_converter: flames_gff3.
FLAMES_gff3_converter: Optional[str] = None


_ENTRY_DEFAULTS = {
    "quant_pattern": None,
    "quant_id_col": None,
    "quant_tpm_col": None,
    "quant_skip_rows": 0,
    "quant_no_header": False,
    "gtf_pattern": None,
    "gtf_converter": None,
    "gtf_source": "ref",
    "color": "gray",
    "display": True,
    "venn": True,
}


def load_registry(path: str) -> List[dict]:
    with open(path, "rt") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict) or "tools" not in data:
        raise ValueError(f"{path}: expected a top-level 'tools' list")

    entries: List[dict] = []
    seen_names = set()
    for raw in data["tools"]:
        if "name" not in raw or "family" not in raw:
            raise ValueError(f"{path}: entry missing 'name' or 'family': {raw!r}")
        if raw["name"] in seen_names:
            raise ValueError(f"{path}: duplicate tool name '{raw['name']}'")
        if raw.get("quant_pattern") is None and raw.get("gtf_pattern") is None:
            raise ValueError(
                f"{path}: entry '{raw['name']}' has neither quant_pattern nor gtf_pattern"
            )
        if raw.get("quant_pattern") is not None and (
            raw.get("quant_id_col") is None or raw.get("quant_tpm_col") is None
        ):
            raise ValueError(
                f"{path}: entry '{raw['name']}' has quant_pattern "
                f"but missing quant_id_col / quant_tpm_col"
            )
        e = {**_ENTRY_DEFAULTS, **raw}
        seen_names.add(e["name"])
        entries.append(e)
    return entries


def process_file(
    input_filename: str, entries: List[dict]
) -> Tuple[Optional[str], Optional[dict], Optional[str]]:
    """
    Match input_filename against the registry. On a quant match,
    write the normalized TSV to processed_prog_results/<name>.tsv and
    return ("quant", entry, that_path). On a GTF match, either run the
    declared converter and return ("gtf", entry, processed_path) or
    pass the raw input through (return ("gtf", entry, input_filename)).
    Returns (None, None, None) on no match.

    First-match-wins. Quant patterns are checked before gtf patterns;
    within each kind, registry order decides.
    """

    bn = os.path.basename(input_filename)

    for e in entries:
        if e["quant_pattern"] is not None and re.search(e["quant_pattern"], bn):
            os.makedirs(PROCESSED_DIR, exist_ok=True)
            output_filename = os.path.join(PROCESSED_DIR, f"{e['name']}.tsv")
            make_tsv(
                input_filename,
                output_filename,
                "\t",
                e["quant_id_col"],
                e["quant_tpm_col"],
                no_header=e["quant_no_header"],
                skip_rows=e["quant_skip_rows"],
            )
            return ("quant", e, output_filename)

    for e in entries:
        if e["gtf_pattern"] is not None and re.search(e["gtf_pattern"], bn):
            converter = e["gtf_converter"]
            if converter is None:
                # raw GTF is consumed as-is
                return ("gtf", e, input_filename)
            if converter == "flames_gff3":
                if FLAMES_gff3_converter is None:
                    raise RuntimeError(
                        "FLAMES gff3 converter path not set "
                        "(QuantParser.FLAMES_gff3_converter)"
                    )
                os.makedirs(PROCESSED_DIR, exist_ok=True)
                out = os.path.join(PROCESSED_DIR, f"{e['name']}.gtf")
                cmd = f"{FLAMES_gff3_converter} {input_filename} > {out}"
                subprocess.check_call(cmd, shell=True)
                return ("gtf", e, out)
            raise ValueError(
                f"entry '{e['name']}': unknown gtf_converter '{converter}'"
            )

    return (None, None, None)


def make_tsv(
    input_filename,
    output_filename,
    delim,
    transcript_id_field,
    tpm_field,
    no_header=False,
    skip_rows=0,
):

    print("-processing {}".format(input_filename))

    tpm_val_dict = dict()
    with open(input_filename, "rt") as fh:
        if skip_rows > 0:
            for _ in range(skip_rows):
                next(fh)

        if not no_header:
            next(fh)
        for line in fh:
            line = line.rstrip()
            vals = line.split(delim)
            transcript_id, tpm = vals[transcript_id_field], vals[tpm_field]
            tpm_val_dict[transcript_id] = float(tpm)

    df = pd.DataFrame(list(tpm_val_dict.items()), columns=["transcript_id", "TPM"])
    df["TPM"] = df["TPM"] / df["TPM"].sum() * 1e6

    df.to_csv(output_filename, sep="\t", index=False)

    print("Done writing {}".format(output_filename))
