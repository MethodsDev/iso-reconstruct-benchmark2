#!/usr/bin/env python3

import sys, os, re
import subprocess
import pandas as pd

FLAMES_gff3_converter = None


def process_file(input_filename):

    # ensure output dir exists.
    output_dir = "processed_prog_results"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ## indicate which columns contain the transcript_id and the read count or TPM value.

    ############
    ## FLAMES ##
    ############

    if re.search("FLAMES.counts.tsv$", input_filename) is not None:
        # FLAMES quant
        output_filename = os.path.join(output_dir, "FLAMES.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            1,
            2,
        )
        return ("FLAMES_quant_file", output_filename)

    elif re.search("FLAMES.gff3$", input_filename) is not None:
        # FLAMES gtf file
        # convert to gtf
        flames_gtf_filename = os.path.join(output_dir, "FLAMES.gtf")
        if FLAMES_gff3_converter is None:
            raise RuntimeError("flames gff3 converter not set")
        cmd = f"{FLAMES_gff3_converter} {input_filename} > {flames_gtf_filename}"
        subprocess.check_call(cmd, shell=True)
        return ("FLAMES_gtf_file", flames_gtf_filename)

    ##############
    ## IsoQuant ##
    ##############

    elif re.search("IsoQuant.counts.tsv$", input_filename) is not None:
        # IsoQuant quant
        output_filename = os.path.join(output_dir, "IsoQuant.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            1,
        )
        return ("IsoQuant_quant_file", output_filename)

    elif re.search("IsoQuant.gtf$", input_filename) is not None:
        # IsoQuant gtf
        return ("IsoQuant_gtf_file", input_filename)

    ############
    ## IsoSeq ##
    ############

    elif re.search("IsoSeq.*\\.abundance.txt$", input_filename) is not None:
        # IsoSeq
        output_filename = os.path.join(output_dir, "IsoSeq.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            2,
            skip_rows=3,
        )
        return ("IsoSeq_quant_file", output_filename)

    elif re.search("IsoSeq.ref-filtered.ID.gff$", input_filename) is not None:
        # ref-guided provided if exists
        return ("IsoSeq_gtf_file", input_filename)
    elif re.search("IsoSeq.ref-free.ID.gff$", input_filename) is not None:
        # otherwise, ref-free
        return ("IsoSeq_gtf_file", input_filename)

    ##########
    ## LRAA ##
    ##########

    elif re.search("LRAA.*.quant.expr$", input_filename) is not None:
        # LRAA
        output_filename = os.path.join(output_dir, "LRAA.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            1,
            3,
        )
        return ("LRAA_quant_file", output_filename)

    elif re.search("LRAA.*.gtf$", input_filename) is not None:
        output_filename = os.path.join(output_dir, "LRAA.gtf")
        with open(output_filename, "wt") as ofh:
            with open(input_filename, "rt") as fh:
                for line in fh:
                    line = line.rstrip()
                    if line != "":
                        print(line, file=ofh)

        return ("LRAA_gtf_file", output_filename)

    #################
    ## Mandalorion ##
    #################

    elif (
        re.search("Mandalorian.Isoforms.filtered.clean.quant$", input_filename)
        is not None
    ):
        # Mandalorion
        output_filename = os.path.join(output_dir, "Mandalorion.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            2,
        )
        return ("Mandalorion_quant_file", output_filename)

    elif (
        re.search("Mandalorian.Isoforms.filtered.clean.gtf$", input_filename)
        is not None
    ):
        return ("Mandalorion_gtf_file", input_filename)

    ###########################
    ## Oarfish - byAlignment ##
    ###########################

    elif re.search("Oarfish.byAlignment.quant$", input_filename) is not None:
        # Oarfish - byAlignment
        output_filename = os.path.join(output_dir, "Oarfish_align.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            2,
        )
        return ("Oarfish_align_quant_file", output_filename)

    #######################
    ## Oarfish - byReads ##
    #######################

    elif re.search("Oarfish.byReads.quant$", input_filename) is not None:
        # Oarfish - byReads
        output_filename = os.path.join(output_dir, "Oarfish_reads.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            2,
        )
        return ("Oarfish_reads_quant_file", output_filename)

    ##############
    ## ESPRESSO ##
    ##############

    elif re.search("espresso.counts.tsv$", input_filename) is not None:
        # ESPRESSO
        output_filename = os.path.join(output_dir, "espresso.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            3,
        )
        return ("ESPRESSO_quant_file", output_filename)

    elif re.search("espresso.gtf$", input_filename) is not None:
        return ("ESPRESSO_gtf_file", input_filename)

    ###########
    ## FLAIR ##
    ###########

    elif re.search("flair.counts.tsv$", input_filename) is not None:
        # flair
        output_filename = os.path.join(output_dir, "flair.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            1,
        )
        return ("FLAIR_quant_file", output_filename)

    elif re.search("flair.gtf$", input_filename) is not None:
        return ("FLAIR_gtf_file", input_filename)

    ###############
    ## Isosceles ##
    ###############

    elif re.search("isosceles.*.counts$", input_filename) is not None:
        # Isosceles
        output_filename = os.path.join(output_dir, "isosceles.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            1,
        )
        return ("Isosceles_quant_file", output_filename)

    elif re.search("isosceles.*.gtf$", input_filename) is not None:
        return ("Isosceles_gtf_file", input_filename)

    ###########
    ## Bambu ##
    ###########

    elif re.search("bambu.counts.txt$", input_filename) is not None:
        # bambu
        output_filename = os.path.join(output_dir, "bambu.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            0,
            2,
        )
        return ("Bambu_quant_file", output_filename)

    elif re.search("bambu.gtf$", input_filename) is not None:
        return ("Bambu_gtf_file", input_filename)

    ###############
    ## Stringtie ##
    ###############

    elif re.search("stringtie.quant.tsv$", input_filename) is not None:
        # StringTie
        output_filename = os.path.join(output_dir, "stringtie.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            1,
            4,
            True,
        )
        return ("StringTie_quant_file", output_filename)

    elif re.search("stringtie.gtf$", input_filename) is not None:
        return ("StringTie_gtf_file", input_filename)

    elif re.search("talon_abundance_filtered.tsv$", input_filename) is not None:
        # TALON
        output_filename = os.path.join(output_dir, "talon.tsv")
        make_tsv(
            input_filename,
            output_filename,
            "\t",
            3,
            11,
        )
        return ("TALON_quant_file", output_filename)

    elif re.search("talon.gtf$", input_filename) is not None:
        return ("TALON_gtf_file", input_filename)

    # nothing recognized
    return (None, None)


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
            header = next(fh)
        for line in fh:
            line = line.rstrip()
            vals = line.split(delim)
            transcript_id, tpm = vals[transcript_id_field], vals[tpm_field]
            tpm_val_dict[transcript_id] = float(tpm)

    df = pd.DataFrame(list(tpm_val_dict.items()), columns=["transcript_id", "TPM"])
    df["TPM"] = df["TPM"] / df["TPM"].sum() * 1e6

    df.to_csv(output_filename, sep="\t", index=False)

    print("Done writing {}".format(output_filename))

    return
