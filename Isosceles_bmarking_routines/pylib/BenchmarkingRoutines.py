#!/usr/bin/env python3

# benchmarking methods were ported over from
# https://github.com/Genentech/Isosceles_Paper/blob/devel/reports/simulated_bulk_benchmarks.ipynb
# of
# https://www.nature.com/articles/s41467-024-51584-3
# Kabza et al., Nature Communications volume 15, Article number: 7316 (2024)


import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import csv
import sklearn.metrics
import numpy as np
import pandas as pd
from scipy import stats as stat
from gtfparse import read_gtf
from collections import OrderedDict


__PALETTE = dict()  # stores (color, linetype) for each sampleName
# visit for color options:
# https://matplotlib.org/stable/gallery/color/named_colors.html


def getFiles(path, dataType):
    """
    Identify GTFs and count files in a directory.
    - dataType: a string, either 'sample' or 'reference'.
    """
    gtfList, countList = [], []
    for x in os.listdir(path):
        if x.endswith(".gtf"):
            gtf = path + x
            gtfList.append(gtf)
        if x.endswith(".tsv"):
            count = path + x
            countList.append(count)

    gtfList = sorted(gtfList)
    countList = sorted(countList)
    if dataType == "sample":
        return list(zip(gtfList, countList))
    elif dataType == "reference":
        return list(zip(gtfList, countList * len(gtfList)))
    else:
        print("Must specify dataType argument.")


def processGtf(gtf):
    """
    Read in GTF file, subset for only multi-exon transcripts and
    add a column that shows the next row's transcript ID.
    """
    print("-processGtf( {} )\n".format(gtf))
    df_polar = read_gtf(gtf)
    df = pd.DataFrame(df_polar)
    df.columns = df_polar.columns
    df = df[df["feature"] == "exon"]
    # print(df.head())

    df_exons = df.sort_values(by=["transcript_id", "start"])

    df_exons = df_exons[df_exons.transcript_id.duplicated(keep=False)]
    df_exons["next_transcript"] = df_exons["transcript_id"].shift(-1)
    df_exons = df_exons.astype({"next_transcript": str})

    return df_exons


def intronIds(df, include_strand=True):
    """Derive intron strings from the processed GTF dataframe; df_exons."""
    intronDict = {}  # Dict; {transcript ID: intron string}.
    transcriptGeneDict = {}
    for i, j in df.iterrows():
        txId, nextTx, chrom, strand, start, end = (
            j["transcript_id"],
            j["next_transcript"],
            j["seqname"],
            j["strand"],
            j["start"],
            j["end"],
        )

        if not include_strand:
            strand = "?"

        if "gene_id" in j:
            transcriptGeneDict[txId] = j["gene_id"]

        if txId not in intronDict.keys():  # If the exon is the first of the tx,
            intronDict[txId] = (
                str(chrom) + ":" + str(end) + "-"
            )  # start the intron string.
        elif (
            txId in intronDict.keys() and txId != nextTx
        ):  # If the exon is the last of the tx,
            intronDict[txId] += (
                str(start - 1) + ":" + str(strand)
            )  # finish the intron string.
        else:  # Else continue the intron string.
            intronDict[txId] += (
                str(start - 1)
                + ":"
                + str(strand)
                + ","
                + str(chrom)
                + ":"
                + str(end)
                + "-"
            )
    intronDf = pd.DataFrame.from_dict(
        intronDict, orient="index"
    )  # Convert intron dictionary to dataframe.
    intronDf.columns = ["intronId"]

    if len(transcriptGeneDict) > 0:
        transcriptGeneDf = pd.DataFrame.from_dict(transcriptGeneDict, orient="index")
        transcriptGeneDf.columns = ["gene_id"]
        intronDf = intronDf.join(transcriptGeneDf, how="left")

    intronDf.index.name = "transcript_id"

    return intronDf


def parseGTFtoIntronIDs(gtf_filename, include_strand=True):

    # Produce a dataframe of transcript IDs and their new intron string IDs.
    df_exons = processGtf(gtf_filename)
    intronDf = intronIds(df_exons, include_strand)

    intronDf.reset_index(inplace=True)

    return intronDf


def parseGTFtoIntronIDsandQuants(
    gtf_filename,
    quant_tsv,
    transcript_id_modifier_func=None,
    include_strand_in_intronId=True,
):
    """
    Merge counts with newly assigned intron string IDs.
    - 'gtf': path to GTF.
    - 'tsv': path to TSV of counts.
    """

    intronDf = parseGTFtoIntronIDs(
        gtf_filename, include_strand=include_strand_in_intronId
    )

    countDf = parseQuantsTSV(quant_tsv)

    if transcript_id_modifier_func is not None:
        intronDf["transcript_id"] = intronDf["transcript_id"].apply(
            transcript_id_modifier_func
        )
        countDf["transcript_id"] = countDf["transcript_id"].apply(
            transcript_id_modifier_func
        )

    df = intronDf.merge(countDf, how="inner", on="transcript_id")

    return df


def indexDfByIntronId(df):

    concat = lambda x: ",".join(x)

    agg_operations = {
        "transcript_ids": ("transcript_id", concat),
    }

    if "gene_id" in df.columns:
        agg_operations["gene_ids"] = ("gene_id", concat)

    if "tpm" in df.columns:
        agg_operations["tpm"] = ("tpm", "sum")

    df = df.groupby("intronId").agg(**agg_operations)

    return df


def parseQuantsTSV(quant_tsv):

    # Read in a count file and join it to the "intronDf".
    countDf = pd.read_csv(quant_tsv, sep="\t")
    colnames = list(countDf.columns.values)
    countDf = countDf[
        [colnames[0], colnames[-1]]
    ]  # colnames[0] = txIds, colnames[-1] = counts.
    countDf.columns = ["transcript_id", "tpm"]
    countDf = countDf.groupby(
        "transcript_id"
    ).sum()  # sets the transcript_id as the index automatically.

    countDf.columns = ["tpm"]

    countDf.reset_index(inplace=True)

    return countDf


def absRelativeDiffEach(ground_truth_TPMs, estimated_TPMs):
    """
    Calculate absolute relative difference between ground truth
    and estimated quantifications.  Metric: 'mean' or 'median'.
    """

    x = np.array(ground_truth_TPMs)  # Ground truth TPMs.
    y = np.array(estimated_TPMs)  # Estimated TPMs.

    absRelDiff = abs(x - y) / ((x + y) / 2)

    return absRelDiff


def relativeDiff(ground_truth_TPMs, estimated_TPMs, metric):

    absRelDiff = absRelativeDiffEach(ground_truth_TPMs, estimated_TPMs)

    if metric == "mean":
        return np.mean(absRelDiff)
    elif metric == "median":
        return np.median(absRelDiff)
    else:
        raise RuntimeError("Specify 'mean' or 'median' for argument: 'metric'.")


def assign_binned_expr_quantile(i_ref_df, i_sample_df, num_bins):
    """
    Takes a reference dataframe, outer-joins it to a sample's counts
    and then bins the dataframe based on ground truth TPMs.
    """

    # Outer join the sample counts to the reference dataframe by intron string ID.
    # bigDf's column 4 is ground truth TPMs and column 5 is the sample's estimated TPMs.
    bigDf = (
        i_ref_df[["tpm"]]
        .copy()
        .rename(columns={"tpm": "ref_tpm"})
        .join(i_sample_df, how="outer")
        .fillna(0)
    )

    for counts in [bigDf["ref_tpm"], bigDf["tpm"]]:
        counts = counts.astype(
            float
        )  # Convert ground truth & estimated counts to float type.
        counts = (counts / np.sum(counts)) * 1000000  # Re-normalize to TPMs.

    groundTruth = "ref_tpm"
    estimated = "tpm"
    # Remove any rows in which both ground truth and estimated TPM is 0.
    dropRows = bigDf[(bigDf[groundTruth] == 0) & (bigDf[estimated] == 0)].index
    bigDf.drop(dropRows, inplace=True)
    # Subset bigDf for rows where ground truth is 0; these are all false positives in the sample.
    zeros = bigDf[bigDf[groundTruth] == 0].copy()
    # Subset bigDf for rows where ground truth isn't 0, then order them by ground truth expression.
    nonzeros = bigDf[bigDf[groundTruth] > 0].sort_values(by=[groundTruth]).copy()
    # Bin the data into n bins based on ground truth expression.
    # .rank(method='first') is required because otherwise we have issues with nonunique bin edge values.
    nonzeros["quantile"] = pd.qcut(
        nonzeros[groundTruth].rank(method="first"), q=num_bins, labels=False
    )
    # Assign false positives to bins based on their detected expression
    quantile_expr = nonzeros[groundTruth].groupby(nonzeros["quantile"]).agg("max")
    quantile_expr = quantile_expr + np.arange(len(quantile_expr)) / 1e10
    quantile_expr = [0] + list(quantile_expr)
    zeros["quantile"] = pd.cut(zeros[estimated], quantile_expr, labels=False)
    n_0 = (zeros["quantile"] == 0).sum()
    zeros["quantile"][zeros["quantile"] == 0] = np.floor(np.arange(n_0) / n_0 * 3)
    n_3 = (zeros["quantile"] == 3).sum()
    zeros["quantile"][zeros["quantile"] == 3] = np.floor(np.arange(n_3) / n_3 * 2) + 3
    # Recombine the false and true positives into the full bigDf.  The reason bigDf is first split into
    # true and false positives is to ensure that the first quantile bin isn't all zeros for ground truth.
    # bigDf = zeros.append(nonzeros)  # deprecated in pandas
    bigDf = pd.concat([zeros, nonzeros], ignore_index=False)
    # bigDf.to_csv("ladeda2", sep="\t", index=False, mode="a")
    return bigDf


def measure_rel_diff_by_expr_quantile(
    i_ref_df, i_sample_df, errorType, num_bins, intronIds_use=None
):
    """Calculate and plot error metrics for a given sample."""
    i_bigDf = assign_binned_expr_quantile(i_ref_df, i_sample_df, num_bins)

    # incorporate into the i_bigDf an indicator of a reference intron set.
    i_ref_introns = i_ref_df.copy()
    i_ref_introns["is_ref"] = True
    i_ref_introns = i_ref_introns[["is_ref"]]
    i_bigDf = i_bigDf.join(i_ref_introns)
    i_bigDf["is_ref"] = i_bigDf["is_ref"].fillna(False)

    i_bigDf["use"] = False  # init
    if intronIds_use is None:
        i_bigDf["use"] = i_bigDf["is_ref"]
    else:
        i_bigDf.loc[intronIds_use, "use"] = True

    # 'knownDf': the subset of tx that were included in the downsampled GTF.
    i_use_df = i_bigDf[i_bigDf["use"]].copy()
    knownGroundTruth = np.array(i_use_df["ref_tpm"])
    knownEstimated = np.array(i_use_df["tpm"])

    # Calculate the heights of the horizontal marks in the sidebar which denote either
    # the median or the mean of the big picture of counts in annotated/all transcripts.
    sidebarError = None
    if errorType == "mean":
        meanRD = relativeDiff(knownGroundTruth, knownEstimated, "mean")
        sidebarError = meanRD
    elif errorType == "median":
        medianRD = relativeDiff(knownGroundTruth, knownEstimated, "median")
        sidebarError = medianRD
    else:
        raise RuntimeError("Specify 'mean' or 'median' for argument: 'errorType'.")

    # Calculate error for each quantile, 0 to n - 1.
    # Each error value represents the error for its quantile bin, 'i'
    # and is stored in a list (knownError or allError) to be plotted later.

    i_use_df["relDiffs"] = -1

    relDiffs = []
    for i in range(num_bins):
        i_quantile_df = i_use_df[i_use_df["quantile"] == i]

        quantile_rel_diffs = None

        rel_diffs = absRelativeDiffEach(i_quantile_df["ref_tpm"], i_quantile_df["tpm"])

        i_use_df.loc[i_quantile_df.index, "relDiffs"] = rel_diffs

        if errorType == "mean":
            quantile_rel_diffs = np.mean(rel_diffs)

        else:
            quantile_rel_diffs = np.median(rel_diffs)

        relDiffs.append(quantile_rel_diffs)

    return relDiffs, sidebarError, i_use_df


def sensitivity(df):
    """Calculate TPR from dataframe assuming column 4 is ground truth and column 5 is estimated."""
    # Equation for Sensitivity: TPR = TP / (TP + FN).
    # TP is any instance of a transcript with TPM > 0 in both ground truth & estimated.
    # FN is any instance of a transcript with TPM > 0 in only ground truth but not estimated.

    # TP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] != 0))
    # FN = np.sum((df.iloc[:, 5] == 0) & (df.iloc[:, 4] != 0))

    TP = np.sum((df["class"] == "TP"))
    FN = np.sum((df["class"] == "FN"))

    if TP + FN != 0:
        TPR = TP / (TP + FN)
    else:
        TPR = 0  # In case there are no TP or FN in a bin, avoid 0-division error.

    return TPR


def fdr(df):
    """Calculate FDR from dataframe assuming column 4 is ground truth and column 5 is estimated."""
    # Equation for False Discovery Rate: FDR = FP / (FP + TP).
    # TP is any instance of a transcript with TPM > 0 in both ground truth & estimated.
    # FP is any instance of a transcript with TPM > 0 in only estimated but not ground truth.

    # TP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] != 0))
    # FP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] == 0))

    TP = np.sum((df["class"] == "TP"))
    FP = np.sum((df["class"] == "FP"))

    if FP + TP != 0:
        FDR = FP / (FP + TP)
    else:
        FDR = 0  # In case there are no TP or FP in a bin, avoid 0-division error.

    return FDR


def precision(df):
    """Calculate precision from dataframe assuming column 4 is ground truth and column 5 is estimated."""
    # Equation for Precision: PPV = TP / (FP + TP).
    # TP is any instance of a transcript with TPM > 0 in both ground truth & estimated.
    # FP is any instance of a transcript with TPM > 0 in only estimated but not ground truth.

    # TP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] != 0))
    # FP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] == 0))

    TP = np.sum((df["class"] == "TP"))
    FP = np.sum((df["class"] == "FP"))

    if FP + TP != 0:
        PPV = TP / (FP + TP)
    else:
        PPV = 1  # In case there are no TP or FP in a bin, avoid 0-division error.

    return PPV


def f1_score(df):
    """Calculate F1 score from dataframe assuming column 4 is ground truth and column 5 is estimated."""
    # Equation for F1 score: F1 = 2*TP / (2*TP + FP + FN).
    # TP is any instance of a transcript with TPM > 0 in both ground truth & estimated.
    # FP is any instance of a transcript with TPM > 0 in only estimated but not ground truth.
    # FN is any instance of a transcript with TPM > 0 in only ground truth but not estimated.

    # TP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] != 0))
    # FP = np.sum((df.iloc[:, 5] != 0) & (df.iloc[:, 4] == 0))
    # FN = np.sum((df.iloc[:, 5] == 0) & (df.iloc[:, 4] != 0))

    TP = np.sum((df["class"] == "TP"))
    FP = np.sum((df["class"] == "FP"))
    FN = np.sum((df["class"] == "FN"))

    F1 = 2 * TP / (2 * TP + FP + FN)

    return F1


def calc_TP_FP_FN(i_ref_df, i_sample_df):
    """Calculate TP, FP and FN to get Sensitivty and False Discovery Rate of a sample."""

    i_ref_tmp_df = i_ref_df[["tpm"]].copy().rename(columns={"tpm": "ref_tpm"})
    i_ref_tmp_df["is_ref"] = True

    i_merged_df = i_ref_tmp_df.join(i_sample_df, how="outer")
    i_merged_df["is_ref"] = i_merged_df["is_ref"].fillna(False)
    i_merged_df["ref_tpm"] = i_merged_df["ref_tpm"].fillna(0)
    i_merged_df["tpm"] = i_merged_df["tpm"].fillna(0)

    # set TP, FP, FN
    i_merged_df["class"] = "?"

    # TP
    i_merged_df.loc[
        (
            (i_merged_df["is_ref"] == True)
            & (i_merged_df["ref_tpm"] > 0)
            & (i_merged_df["tpm"] > 0)
        ),
        "class",
    ] = "TP"

    # FP
    i_merged_df.loc[
        ((i_merged_df["is_ref"] == False) & (i_merged_df["tpm"] > 0)), "class"
    ] = "FP"

    # FN
    i_merged_df.loc[
        (i_merged_df["is_ref"] == True)
        & (i_merged_df["ref_tpm"] > 0)
        & (i_merged_df["tpm"] == 0),
        "class",
    ] = "FN"

    return i_merged_df


def calc_TP_FP_FN_and_expr_quantiles(i_ref_df, i_sample_df, num_bins):

    # calc TP, FP, and FN
    i_sample_TP_FP_FN_df = calc_TP_FP_FN(i_ref_df, i_sample_df)

    # incorporate the expr quantiles info
    i_expr_quantiles = assign_binned_expr_quantile(i_ref_df, i_sample_df, num_bins)
    i_sample_TP_FP_FN_df = i_sample_TP_FP_FN_df.join(
        i_expr_quantiles["quantile"], how="left"
    )

    return i_sample_TP_FP_FN_df


def calc_knownTPR_novelTPR_and_FDR_for_quantiles(
    i_sample_TP_FP_FN_df, num_bins, novel_intron_ids=None
):

    novel_introns = set()
    if novel_intron_ids is not None:
        novel_introns.update(novel_intron_ids)

    # Calculate TPR and FDR for each quantile bin 'i'.
    knownTPR = []
    novelTPR = []
    novelFDR = []
    # Subset each dataframe for rows in quantile bin 'i'.
    for i in range(num_bins):
        subBigDf = i_sample_TP_FP_FN_df[i_sample_TP_FP_FN_df["quantile"] == i]

        stats_quantile_i = calc_accuracy_stats_for_TP_FP_FN_df(subBigDf)

        # fdr based on full df
        fdr_val = stats_quantile_i["FDR"]
        novelFDR.append(fdr_val)

        sensitivity_known = stats_quantile_i["Sensitivity"]

        if novel_intron_ids is None:
            # sensitivity known based on all ref transcripts.
            sensitivity_novel = 0
            novelTPR.append(sensitivity_novel)

        else:
            sensitivity_novel = 0  # init
            introns_in_quantile = set(subBigDf.index.values)
            novel_introns_in_quantile = introns_in_quantile & novel_introns
            not_novel_introns_in_quantile = (
                introns_in_quantile - novel_introns_in_quantile
            )
            if len(novel_introns_in_quantile) > 0:
                novel_stats_quantile_i = calc_accuracy_stats_for_TP_FP_FN_df(
                    subBigDf.loc[list(novel_introns_in_quantile), :]
                )
                sensitivity_novel = novel_stats_quantile_i["Sensitivity"]
                novelTPR.append(sensitivity_novel)

                # redo known sensitivity value based on the non-novel introns.
                not_novel_stats_quantile_i = calc_accuracy_stats_for_TP_FP_FN_df(
                    subBigDf.loc[list(not_novel_introns_in_quantile), :]
                )
                sensitivity_known = not_novel_stats_quantile_i["Sensitivity"]
            else:
                novelTPR.append(np.nan)

        knownTPR.append(sensitivity_known)

    statistics = [knownTPR, novelTPR, novelFDR]
    meanStats = [np.mean(knownTPR), np.mean(novelTPR), np.mean(novelFDR)]

    return statistics, meanStats


def calc_TPR_PPV_F1_for_quantiles(i_sample_TP_FP_FN_df, num_bins, novel_intron_ids):
    """Calculate TP, FP and FN to get Sensitivty, F1 Score and Precision of a sample."""

    novel_introns = set()
    if novel_intron_ids is not None:
        novel_introns.update(novel_intron_ids)

    # Calculate TPR and FDR for each quantile bin 'i'.
    allTPR = []
    allPPV = []
    allF1 = []

    # Subset each dataframe for rows in quantile bin 'i'.
    for i in range(num_bins):
        subBigDf = i_sample_TP_FP_FN_df[i_sample_TP_FP_FN_df["quantile"] == i]

        if novel_intron_ids is not None:
            introns_in_quantile = set(subBigDf.index.values)
            novel_introns_in_quantile = introns_in_quantile & novel_introns
            if len(novel_introns_in_quantile) > 0:
                FPs_df = subBigDf[subBigDf["class"] == "FP"]
                subBigDf = subBigDf.loc[list(novel_introns_in_quantile), :]
                # above would only incorporate TP and FNs according to the novel intron ids.
                # must incorporate all the FPs too
                subBigDf = pd.concat([subBigDf, FPs_df])

        sensitivity_val = sensitivity(subBigDf) if len(subBigDf) > 0 else np.nan
        allTPR.append(sensitivity_val)

        # Measure de novo detection metrics and append to lists.
        F1_val = f1_score(subBigDf) if len(subBigDf) > 0 else np.nan
        allF1.append(F1_val)

        PPV_val = precision(subBigDf) if len(subBigDf) > 0 else np.nan
        allPPV.append(PPV_val)

    statistics = [allTPR, allPPV, allF1]
    meanStats = [np.mean(allTPR), np.mean(allPPV), np.mean(allF1)]

    return statistics, meanStats


def set_color_palette(sampleName, color, linetype):
    __PALETTE[sampleName] = (color, linetype)
    return


def colorAndLabel(sampleName):
    """Assign kwargs to a sample's plotting assets."""
    # Assign the sample its legend label (name), color (c), & line style (l).

    if sampleName not in __PALETTE:
        raise RuntimeError(
            "Error, {} not in palette. Use BenchmarkingRoutines.set_color_palette() for it first.".format(
                sampleName
            )
        )

    color, linetype = __PALETTE[sampleName]

    return sampleName, color, linetype

    """
    name = None

    if "flair" in sampleName:
        c = "mediumpurple"
        if "stringtie" not in sampleName:
            name = "FLAIR"
        else:
            name = "FLAIR + StringTie2"
    if "isosceles" in sampleName:
        c = "black"
        if "stringtie" in sampleName:
            name = "Isosceles + StringTie2"
        elif "IsoQuant" in sampleName:
            name = "Isosceles + IsoQuant"
        else:
            name = "Isosceles"
    if "isoquant" in sampleName:
        c = "gold"
        if "stringtie" not in sampleName:
            name = "IsoQuant"
        else:
            name = "IsoQuant + StringTie2"
    if "bambu" in sampleName:
        c = "mediumseagreen"
        if "stringtie" not in sampleName:
            name = "Bambu"
        else:
            name = "Bambu + StringTie2"
    if "nanocount" in sampleName:
        c = "cornflowerblue"
        if "stringtie" not in sampleName:
            name = "NanoCount"
        else:
            name = "NanoCount + StringTie2"
    if "liqa" in sampleName:
        c = "deeppink"
        if "stringtie" not in sampleName:
            name = "LIQA"
        else:
            name = "LIQA + StringTie2"
    if "espresso" in sampleName:
        c = "lightsalmon"
        if "stringtie" not in sampleName:
            name = "ESPRESSO"
        else:
            name = "ESPRESSO + StringTie2"

    if "LRAA" in sampleName:
        c = "purple"
        name = "LRAA"

    if is_best:
        if ("isosceles" in sampleName) and ("stringtie" in sampleName):
            l = (0, (2, 1))
        elif ("isosceles" in sampleName) and ("IsoQuant" in sampleName):
            l = (0, (1, 1))
        else:
            l = "solid"
    else:
        if "stringtie" in sampleName:
            l = (0, (2, 1))
        elif "IsoQuant" in sampleName:
            l = (0, (1, 1))
        else:
            l = "solid"



    if name is None:
        raise RuntimeError(f"Error, {sampleName} not recognized for color assignment")

    return [name, c, l]
    """


def percentileTicks(n):
    """
    Takes the number of bins, 'n', and generates a range for
    plotting across an x-axis of range(0, 99).
    """

    m = int(99 / n)
    ticks = [m * l for l in range(0, n)]
    ticks[-1] += m

    return ticks


def quantileTpmDict(ref, n):
    ref = ref.copy()
    counts = ref[ref.columns[4]]
    counts = counts.astype(float)
    counts = (counts / np.sum(counts)) * 1000000
    ref[ref.columns[4]] = counts
    ref = ref.sort_values(by=[ref.columns[4]])
    ref["quantile"] = pd.qcut(
        ref[ref.columns[4]].rank(method="first"), q=n, labels=False
    )
    quantile_df = ref.groupby("quantile").agg(tpm=(ref.columns[4], "mean"))
    quantile_df["quantile"] = percentileTicks(n)
    quantile_df.index.name = None
    quantile_tpm_dict = {}
    for quantile, tpm in zip(quantile_df["quantile"], quantile_df["tpm"]):
        quantile_tpm_dict[quantile] = tpm
    return quantile_tpm_dict


def scatterplot_adj(i_ref_df, progname_to_df_dict):
    """
    Generates side-by-side scatterplots comparing observed and
    ground truth TPMs from all tools.
    'ref': reference dataframe
    'dfs': list of dataframes made with countsWithIntronIds().
    'n': an integer representing the number of tools being plotted.
    """

    assert (
        i_ref_df.index.name == "intronId"
    ), "Error, i_ref_df input not indexed on intronId"

    ref_quants = i_ref_df.copy().rename(columns={"tpm": "ref_tpm"})

    countDict = {}
    for progname, df in progname_to_df_dict.items():

        assert (
            df.index.name == "intronId"
        ), "Error, df for {} not indexed on intronId".format(progname)

        prog_tpm_colname = progname + "_tpm"

        prog_quants = df[["tpm"]].copy().rename(columns={"tpm": prog_tpm_colname})

        bigDf = prog_quants.join(ref_quants, how="inner").fillna(0)
        colnames = list(bigDf.columns.values)

        name, c, l = colorAndLabel(progname)
        for counts in [bigDf["ref_tpm"], bigDf[prog_tpm_colname]]:
            counts = counts.astype(
                float
            )  # Convert ground truth & estimated counts to float type.
            counts = (counts / np.sum(counts)) * 1000000  # Re-normalize to TPMs.

        groundTruth = np.array(bigDf["ref_tpm"])
        estimated = np.array(bigDf[prog_tpm_colname])
        countDict[name] = [groundTruth, estimated, c]  # Append counts to a dictionary.

        bigDf.copy().reset_index().to_csv(
            progname + ".ref_quant_compare.tsv", sep="\t", index=False
        )

    countDict = OrderedDict(sorted(countDict.items()))  # Sort count dictionary by keys

    # Define figure and panel dimensions
    n = len(progname_to_df_dict)
    fig, ax = plt.subplots(1, n, figsize=(2.7 * n, 2.69), tight_layout=True)

    subplotIndex = 0  # to index the subplots
    for program, tpm in countDict.items():  # by program name, alphabetically.
        print(program)
        print(tpm)
        groundTruth = tpm[0].copy()
        estimated = tpm[1].copy()
        estimated[estimated <= 0.001] = 0.001  # Set all estimated values that are
        # <= 0.001 to 0.001 to show instances of
        # undetected transcripts at x = 0.001.
        color = tpm[2]
        corr = r"$\rho$ = " + str(round(stat.spearmanr(tpm[1], tpm[0]).correlation, 3))
        pcorr = "R = " + str(
            round(stat.pearsonr(np.log(tpm[1] + 1), np.log(tpm[0] + 1)).statistic, 3)
        )

        ax[subplotIndex].plot(groundTruth, groundTruth, color="red", lw=1)
        ax[subplotIndex].scatter(estimated, groundTruth, 0.25, c=color, alpha=0.5)
        ax[subplotIndex].text(0.002, 3000, corr)
        ax[subplotIndex].text(0.002, 1000, pcorr)
        ax[subplotIndex].text(0.002, 100, program)

        ax[subplotIndex].set_xlim(0.0005, 10000)
        ax[subplotIndex].set_ylim(0.1, 10000)
        ax[subplotIndex].set_yscale("log")
        ax[subplotIndex].set_xscale("log")
        ax[subplotIndex].set_xticks([0.001, 0.01, 1, 100, 10000])
        ax[subplotIndex].set_xticklabels(
            ["0", r"$10^{-2}$", r"$10^0$", r"$10^2$", r"$10^4$"]
        )
        ax[subplotIndex].set_xlabel("estimated TPM")
        ax[subplotIndex].set_ylabel("ground truth TPM")
        subplotIndex += 1


def ma_plot_adj(i_ref_df, progname_to_df_dict):
    """
    Generates side-by-side MA plots comparing observed and
    ground truth TPMs from all tools.
    'ref': reference dataframe
    'dfs': list of dataframes made with countsWithIntronIds().
    'n': an integer representing the number of tools being plotted.
    """

    assert (
        i_ref_df.index.name == "intronId"
    ), "Error, i_ref_df input not indexed on intronId"

    ref_quants = i_ref_df.copy().rename(columns={"tpm": "ref_tpm"})

    countDict = {}
    for progname, df in progname_to_df_dict.items():

        assert (
            df.index.name == "intronId"
        ), "Error, df for {} not indexed on intronId".format(progname)

        prog_tpm_colname = progname + "_tpm"

        prog_quants = df[["tpm"]].copy().rename(columns={"tpm": prog_tpm_colname})

        bigDf = prog_quants.join(ref_quants, how="inner").fillna(0)

        colnames = list(bigDf.columns.values)

        name, c, l = colorAndLabel(progname)
        for counts in [bigDf["ref_tpm"], bigDf[prog_tpm_colname]]:
            counts = counts.astype(
                float
            )  # Convert ground truth & estimated counts to float type.
            counts = (counts / np.sum(counts)) * 1000000  # Re-normalize to TPMs.

        groundTruth = np.array(bigDf["ref_tpm"])
        estimated = np.array(bigDf[prog_tpm_colname])
        countDict[name] = [groundTruth, estimated, c]  # Append counts to a dictionary.

    countDict = OrderedDict(sorted(countDict.items()))  # Sort count dictionary by keys

    # Define figure and panel dimensions
    n = len(progname_to_df_dict)
    fig, ax = plt.subplots(1, n, figsize=(2.7 * n, 2.69), tight_layout=True)

    subplotIndex = 0  # to index the subplots
    for program, tpm in countDict.items():  # by program name, alphabetically.
        groundTruth = tpm[0].copy()
        estimated = tpm[1].copy()
        estimated[estimated <= 0.001] = 0.001  # Set all estimated values that are
        # <= 0.001 to 0.001 to show instances of
        # undetected transcripts at x = 0.001.
        color = tpm[2]
        ax[subplotIndex].plot(groundTruth, groundTruth / groundTruth, color="red", lw=1)
        ax[subplotIndex].scatter(
            groundTruth, estimated / groundTruth, 0.25, c=color, alpha=0.5
        )
        ax[subplotIndex].set_xlim(1e-1, 1e4)
        ax[subplotIndex].set_ylim(1e-8, 1e8)
        ax[subplotIndex].set_yscale("log")
        ax[subplotIndex].set_xscale("log")
        ax[subplotIndex].set_yticks([1e-8, 1e-5, 1e-2, 1e2, 1e5, 1e8])
        ax[subplotIndex].set_xticks([1e-1, 1, 1e1, 1e2, 1e3, 1e4])
        ax[subplotIndex].set_xlabel("ground truth TPM")
        ax[subplotIndex].set_ylabel("fold change")
        subplotIndex += 1


def cor_spearman_barplot(i_ref_df, progname_to_df_dict):

    assert (
        i_ref_df.index.name == "intronId"
    ), "Error, i_ref_df input not indexed on intronId"

    ref_quants = i_ref_df.copy().rename(columns={"tpm": "ref_tpm"})

    program_names = []
    program_colors = []
    cor_values = []
    for progname, df in progname_to_df_dict.items():

        assert (
            df.index.name == "intronId"
        ), "Error, df for {} not indexed on intronId".format(progname)

        program_tuple = colorAndLabel(progname)
        program_names.append(program_tuple[0])
        program_colors.append(program_tuple[1])
        program_df = df[["tpm"]].join(ref_quants, how="inner").fillna(0)
        cor_values.append(
            stat.spearmanr(program_df[["ref_tpm"]], program_df[["tpm"]]).correlation
        )
    plot_df = pd.DataFrame(
        {"name": program_names, "color": program_colors, "cor_value": cor_values}
    )
    plot_df = plot_df.sort_values(by="cor_value", ascending=False)
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(4, 3), dpi=100, layout="constrained"
    )
    ax.barh(y=plot_df["name"], width=plot_df["cor_value"], color=plot_df["color"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Spearman correlation")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for bars in ax.containers:
        ax.bar_label(bars, fmt="%.2f", fontsize=12)
    fig.show()

    return plot_df


def cor_pearson_barplot(i_ref_df, progname_to_df_dict):

    assert (
        i_ref_df.index.name == "intronId"
    ), "Error, i_ref_df input not indexed on intronId"

    ref_quants = i_ref_df.copy().rename(columns={"tpm": "ref_tpm"})

    program_names = []
    program_colors = []
    cor_values = []

    for progname, df in progname_to_df_dict.items():

        assert (
            df.index.name == "intronId"
        ), "Error, df for {} not indexed on intronId".format(progname)

        program_tuple = colorAndLabel(progname)
        program_names.append(program_tuple[0])
        program_colors.append(program_tuple[1])
        program_df = df[["tpm"]].join(ref_quants, how="inner").fillna(0)

        # Apply log transformation (adding a small constant to avoid log(0))
        log_col1 = np.log1p(program_df["ref_tpm"])
        log_col2 = np.log1p(program_df["tpm"])

        cor_values.append(stat.pearsonr(log_col1, log_col2).statistic)

    plot_df = pd.DataFrame(
        {"name": program_names, "color": program_colors, "cor_value": cor_values}
    )
    plot_df = plot_df.sort_values(by="cor_value", ascending=False)
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(4, 3), dpi=100, layout="constrained"
    )
    ax.barh(y=plot_df["name"], width=plot_df["cor_value"], color=plot_df["color"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Pearson correlation")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for bars in ax.containers:
        ax.bar_label(bars, fmt="%.2f", fontsize=12)
    fig.show()

    return plot_df


def rel_diff_barplot(i_ref_df, progname_to_df_dict, relDiffType):

    assert (
        i_ref_df.index.name == "intronId"
    ), "Error, i_ref_df input not indexed on intronId"

    ref_quants = i_ref_df.copy().rename(columns={"tpm": "ref_tpm"})

    program_names = []
    program_colors = []
    rel_diffs = []
    median_rel_diffs = []

    rellDiffsEachProg = None

    for progname, df in progname_to_df_dict.items():

        assert (
            df.index.name == "intronId"
        ), "Error, df for {} not indexed on intronId".format(progname)

        program_tuple = colorAndLabel(progname)
        program_names.append(program_tuple[0])
        program_colors.append(program_tuple[1])

        program_df = df[["tpm"]].join(ref_quants[["ref_tpm"]], how="right").fillna(0)

        prog_rel_diffs = relativeDiff(
            program_df["ref_tpm"], program_df["tpm"], relDiffType
        )

        rel_diffs.append(prog_rel_diffs)
        median_rel_diffs.append(
            relativeDiff(program_df["ref_tpm"], program_df["tpm"], "median")
        )

        df2 = program_df.copy()
        df2["relativeDiff"] = absRelativeDiffEach(
            program_df["ref_tpm"], program_df["tpm"]
        )
        df2["progname"] = progname
        rellDiffsEachProg = pd.concat([rellDiffsEachProg, df2])

    plot_df = pd.DataFrame(
        {
            "name": program_names,
            "color": program_colors,
            "rel_diff": rel_diffs,
            "median_rel_diff": median_rel_diffs,
        }
    )
    plot_df = plot_df.sort_values(by="median_rel_diff")
    plot_df = plot_df.drop("median_rel_diff", axis=1)
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(4, 3), dpi=100, layout="constrained"
    )
    ax.barh(y=plot_df["name"], width=plot_df["rel_diff"], color=plot_df["color"])
    ax.set_xlim(0, 2)
    ax.set_xlabel(relDiffType.title() + " rel. diff.")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for bars in ax.containers:
        ax.bar_label(bars, fmt="%.2f", fontsize=12)
    fig.show()

    rellDiffsEachProg.to_csv("absRelativeDiffs.tsv", sep="\t", quoting=csv.QUOTE_NONE)

    return plot_df


def rel_diff_vs_expr_percentile_plot(
    i_ref_df, progname_to_df_dict, num_bins, errorType, plot_title, intron_ids_use=None
):
    """
    Generates an error plot for one specific downsampled percentage.
    'ref': The reference dataframe of intron string IDs and ground truth counts.
    'dfs': a list of sample dataframes from getFiles().
    'n': an integer representing the number of bins to stratify by.
    'errorType': 'median' or 'mean'.
    'downsamplePercentage': a 2-character integer representing the downsample % of interest.
    """

    tickRange = percentileTicks(num_bins)

    figureWidth = 5
    figureHeight = 4.3
    panelWidth = 3.2
    panelHeight = 3
    relativePanelWidth = panelWidth / figureWidth
    relativePanelHeight = panelHeight / figureHeight
    plt.figure(figsize=(figureWidth, figureHeight))

    panel = plt.axes([0.8 / 11, 2.4 / 11, relativePanelWidth, relativePanelHeight])
    panel.set_title(plot_title, fontsize=10.0)
    panel.set_xlabel("Ground truth\nexpression percentile", fontsize=10.0)

    panel.set_xticks([0, 25, 50, 75, 99])
    panel.set_xlim(0, 108)
    panel.set_ylim(0, 2.1)
    panel.vlines(x=99.1, ymin=0, ymax=2.1, color="black", linewidth=1, zorder=5)
    rectangle = mplpatches.Rectangle(
        (99.3, 0.05), 107.5 - 99.3, 1.12, color="white", zorder=3
    )
    panel.add_patch(rectangle)

    if errorType == "median":
        panel.set_ylabel("Median rel. diff.", fontsize=10.0)

        panel.set_xticks([0, 25, 50, 75, 99])
        panel.set_xlim(0, 108)
        panel.set_ylim(0, 2.1)
        panel.vlines(x=99.1, ymin=0, ymax=2.1, color="black", linewidth=1, zorder=5)
        rectangle = mplpatches.Rectangle(
            (99.3, 0.05), 107.5 - 99.3, 1.12, color="white", zorder=3
        )
        panel.add_patch(rectangle)

    elif errorType == "mean":
        panel.set_ylabel("Mean rel. diff.", fontsize=10.0)
        panel.set_xticks([0, 25, 50, 75, 99])
        panel.set_xlim(0, 108)
        panel.set_ylim(0, 2.1)
        panel.vlines(x=99.1, ymin=0, ymax=2.1, color="black", linewidth=1, zorder=5)
        rectangle = mplpatches.Rectangle(
            (99.3, 0.05), 107.5 - 99.3, 1.12, color="white", zorder=3
        )
        panel.add_patch(rectangle)
    else:
        raise RuntimeError("Specify 'mean' or 'median' for argument: 'errorType'.")

    audit_reldiffs_df = None

    errorAnnotations = []
    for progname, i_sample_df in progname_to_df_dict.items():

        ###################################
        ## Mesaure the relative differences
        rel_diffs, sidebarError, i_used_df = measure_rel_diff_by_expr_quantile(
            i_ref_df, i_sample_df, errorType, num_bins, intron_ids_use
        )

        i_used_df2 = i_used_df.copy()
        i_used_df2["progname"] = progname
        audit_reldiffs_df = pd.concat([audit_reldiffs_df, i_used_df2])

        name, c, l = colorAndLabel(progname)
        # Assign the sample to its respective panels based on downsample %.

        # Plot the error line and sidebar error.
        panel.plot(
            tickRange,
            rel_diffs,
            marker="None",
            alpha=0.8,
            markersize=1,
            color=c,
            label=progname,
        )
        panel.hlines(
            y=sidebarError,
            xmin=99.1,
            xmax=107.9,
            color=c,
            linewidth=1.5,
            zorder=4,
        )
        # Add info for sidebar error value annotation to be added later.
        errorAnnotations.append(
            [
                name,
                c,
                str("{0:.2f}".format(round(sidebarError, 2))),
            ]
        )
    # Add the second x axis
    ax2 = panel.twiny()
    ax2.set_xlim(0, 108)
    ax2.set_xticks([0, 30, 66, 93])
    ax2.set_xticklabels(["$10^{-1}$", "$10^0$", "$10^1$", "$10^2$"])
    ax2.tick_params(axis="x", which="major", labelsize=8)

    # Sort the sidebar error value by panel number and sample name.
    errorAnnotations = sorted(errorAnnotations, key=lambda x: x[0])

    # initial height, 'h', and incremental decline, 'j', must be adjusted to fit the number of
    # samples in the comparison and will be dependent on the y-axis range.
    if errorType == "median":
        h = 0.55
        j = 0.10
    elif errorType == "mean":
        h = 0.55
        j = 0.10
    # Plot the sidebar error value as a text annotation in the bottom left-hand corner.
    for i in range(len(errorAnnotations)):
        label, color, err = errorAnnotations[i]
        panel.annotate(err, (2.5, h - j * i), size=6, color=color)

    panel.legend(loc="center left", bbox_to_anchor=(2.5, 0.5), frameon=False)

    return audit_reldiffs_df


def IsoformIdentificationSensitivityPlot(
    i_ref_df,
    progname_to_df_dict,
    num_bins,
    errorType,
    plot_title,
    novel_intron_ids=None,
):
    """
    Generates a sensitivity plot for one specific downsampled percentage.
    'ref': The reference dataframe of intron string IDs and ground truth counts.
    'dfs': a list of sample dataframes from getFiles().
    'n': an integer representing the number of bins to stratify by.
    'downsamplePercentage': a 2-character integer representing the downsample % of interest.
    'is_best': a boolean value affecting the style of plotted lines.
    """

    tickRange = percentileTicks(num_bins)

    figureWidth = 12
    figureHeight = 3
    panelWidth = 2.2
    panelHeight = 2
    relativePanelWidth = panelWidth / figureWidth
    relativePanelHeight = panelHeight / figureHeight
    plt.figure(figsize=(figureWidth, figureHeight))

    panel1 = plt.axes([0.6 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    panel1.set_ylim(0.0, 1.1)
    panel1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.9, color="white", zorder=3)
    panel1.add_patch(rectangle)
    title1 = (
        "All ref trans (annot&unannot)"
        if novel_intron_ids is None
        else "Annotated trans subset"
    )
    panel1.set_title(title1, fontsize=10.0)
    panel1.set_ylabel("Sensitivity (TPR)", fontsize=10.0)

    panel2 = plt.axes([3.4 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    title2 = "empty" if novel_intron_ids is None else "Unannotated trans subset"
    panel2.set_title(title2, fontsize=10.0)
    panel2.set_ylabel("Sensitivity (TPR)", fontsize=10.0)
    panel2.set_xlabel(r"Ground truth expression percentile", fontsize=10.0)
    panel2.set_ylim(-0.025, 1.1)
    panel2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.55, color="white", zorder=3)
    panel2.add_patch(rectangle)

    panel3 = plt.axes([6.2 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    panel3.set_title("All detected transcripts", fontsize=10.0)
    panel3.set_ylabel("False Discovery Rate (FDR)", fontsize=10.0)
    panel3.set_ylim(-0.025, 1.1)
    panel3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.5, color="white", zorder=3)
    panel1.add_patch(rectangle)

    panels = [panel1, panel2, panel3]

    for panel in panels:
        panel.set_xlim(0, 108)
        panel.set_xticks([0, 25, 50, 75, 99])
        panel.vlines(x=99.1, ymin=-1, ymax=1.1, color="black", linewidth=1, zorder=5)

    for progname, i_sample_df in progname_to_df_dict.items():

        # calc TP, FP, and FN and incorporate the expr quantiles info
        i_sample_TP_FP_FN_df = calc_TP_FP_FN_and_expr_quantiles(
            i_ref_df, i_sample_df, num_bins
        )

        statistics, meanStats = calc_knownTPR_novelTPR_and_FDR_for_quantiles(
            i_sample_TP_FP_FN_df, num_bins, novel_intron_ids
        )

        name, c, l = colorAndLabel(progname)

        # Plot the error line and sidebar error.
        for i in range(len(panels)):
            panels[i].plot(
                tickRange,
                statistics[i],
                marker="None",
                alpha=0.8,
                markersize=1,
                color=c,
                linestyle=l,
                label=name,
            )
            panels[i].hlines(
                y=meanStats[i],
                xmin=99.1,
                xmax=107.9,
                color=c,
                linestyle=l,
                linewidth=1.5,
                zorder=4,
            )
            # Add the second x axis
            ax2 = panels[i].twiny()
            ax2.set_xlim(0, 108)
            ax2.set_xticks([0, 30, 66, 93])
            ax2.set_xticklabels(["$10^{-1}$", "$10^0$", "$10^1$", "$10^2$"])
            ax2.tick_params(axis="x", which="major", labelsize=8)

    panel1.legend(loc="center left", bbox_to_anchor=(3.85, 0.5), frameon=False)


def TPR_F1_PPV_plot(
    i_ref_df,
    progname_to_df_dict,
    num_bins=33,
    plot_title="TPR_F1_PPV by expr quintile",
    novel_intron_ids=None,
):

    tickRange = percentileTicks(num_bins)

    figureWidth = 12
    figureHeight = 3
    panelWidth = 2.2
    panelHeight = 2
    relativePanelWidth = panelWidth / figureWidth
    relativePanelHeight = panelHeight / figureHeight
    plt.figure(figsize=(figureWidth, figureHeight))

    panel1 = plt.axes([0.6 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    panel1.set_ylim(0.0, 1.1)
    panel1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.9, color="white", zorder=3)
    panel1.add_patch(rectangle)
    panel1.set_title(plot_title, fontsize=10.0)
    panel1.set_ylabel("Sensitivity (TPR)", fontsize=10.0)

    panel2 = plt.axes([3.4 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    panel2.set_title(plot_title, fontsize=10.0)
    panel2.set_ylabel("Precision (PPV)", fontsize=10.0)
    panel2.set_xlabel(r"Ground truth expression percentile", fontsize=10.0)
    panel2.set_ylim(0.0, 1.1)
    panel2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.55, color="white", zorder=3)
    panel2.add_patch(rectangle)

    panel3 = plt.axes([6.2 / 11, 1.8 / 11, relativePanelWidth, relativePanelHeight])
    panel3.set_title(plot_title, fontsize=10.0)
    panel3.set_ylabel("F1 Score", fontsize=10.0)
    panel3.set_ylim(0.0, 1.1)
    panel3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    rectangle = mplpatches.Rectangle((99.3, 0.05), 8.2, 0.5, color="white", zorder=3)
    panel3.add_patch(rectangle)

    panels = [panel1, panel2, panel3]

    for panel in panels:
        panel.set_xlim(0, 108)
        panel.set_xticks([0, 25, 50, 75, 99])
        panel.vlines(x=99.1, ymin=-1, ymax=1.1, color="black", linewidth=1, zorder=5)

    for progname, i_sample_df in progname_to_df_dict.items():

        # calc TP, FP, and FN and incorporate the expr quantiles info
        i_sample_TP_FP_FN_df = calc_TP_FP_FN_and_expr_quantiles(
            i_ref_df, i_sample_df, num_bins
        )

        statistics, meanStats = calc_TPR_PPV_F1_for_quantiles(
            i_sample_TP_FP_FN_df, num_bins, novel_intron_ids
        )

        name, c, l = colorAndLabel(progname)

        # Plot the error line and sidebar error.
        for i in range(len(panels)):
            panels[i].plot(
                tickRange,
                statistics[i],
                marker="None",
                alpha=0.8,
                markersize=1,
                color=c,
                linestyle=l,
                label=name,
            )
            panels[i].hlines(
                y=meanStats[i],
                xmin=99.1,
                xmax=107.9,
                color=c,
                linestyle=l,
                linewidth=1.5,
                zorder=4,
            )
            # Add the second x axis
            ax2 = panels[i].twiny()
            ax2.set_xlim(0, 108)
            ax2.set_xticks([0, 30, 66, 93])
            ax2.set_xticklabels(["$10^{-1}$", "$10^0$", "$10^1$", "$10^2$"])
            ax2.tick_params(axis="x", which="major", labelsize=8)

    panel1.legend(loc="center left", bbox_to_anchor=(3.85, 0.5), frameon=False)


def calc_accuracy_stats_for_TP_FP_FN_df(i_sample_TP_FP_FN_df, novel_intron_ids=None):

    stats = {"Sensitivity": np.nan, "PPV": np.nan, "FDR": np.nan, "F1": np.nan}

    if novel_intron_ids is not None:
        novel_introns = set(novel_intron_ids)

    if novel_intron_ids is not None:
        introns_in_df = set(i_sample_TP_FP_FN_df.index.values)
        novel_introns = introns_in_df & novel_introns
        if len(novel_introns) > 0:
            FPs_df = i_sample_TP_FP_FN_df[i_sample_TP_FP_FN_df["class"] == "FP"]
            i_sample_TP_FP_FN_df = i_sample_TP_FP_FN_df.loc[list(novel_introns), :]
            # above would only incorporate TP and FNs according to the novel intron ids.
            # must incorporate all the FPs too
            i_sample_TP_FP_FN_df = pd.concat([i_sample_TP_FP_FN_df, FPs_df])
        else:
            # nothing to get stats for.
            return stats

    stats["Sensitivity"] = sensitivity(i_sample_TP_FP_FN_df)
    stats["PPV"] = precision(i_sample_TP_FP_FN_df)
    stats["FDR"] = fdr(i_sample_TP_FP_FN_df)
    stats["F1"] = f1_score(i_sample_TP_FP_FN_df)

    return stats


def overall_knownTPR_novelTPR_and_FDR_barplot(
    i_ref_df, progname_to_df_dict, novel_intron_ids=None
):
    """
    Generates sensitivity statistics bar plots for downsampled data.
    'ref': a pandas dataframe of reference transcript IDs and ground truth counts.
    'dfs': a list of pandas dataframes of counts.
    'n': an integer representing the number of bins out of 99 to stratify by.
    """

    df_data = []

    all_TP_FP_FN_df = None

    for progname, i_sample_df in progname_to_df_dict.items():

        i_sample_TP_FP_FN_df = calc_TP_FP_FN(i_ref_df, i_sample_df)

        all_TP_FP_FN_df = pd.concat([all_TP_FP_FN_df, i_sample_TP_FP_FN_df])

        accuracy_stats = calc_accuracy_stats_for_TP_FP_FN_df(i_sample_TP_FP_FN_df)

        FDR_val = accuracy_stats["FDR"]
        sensitivity_val = accuracy_stats["Sensitivity"]
        novel_trans_sensitivity = np.nan
        F1_val = accuracy_stats["F1"]

        if novel_intron_ids is not None:
            novel_intron_ids = set(novel_intron_ids)
            accuracy_stats_for_novel = calc_accuracy_stats_for_TP_FP_FN_df(
                i_sample_TP_FP_FN_df, novel_intron_ids
            )
            novel_trans_sensitivity = accuracy_stats_for_novel["Sensitivity"]

            # redo sensitivity_val calc based on the non-novel introns.
            not_novel_intron_ids = (
                set(i_sample_TP_FP_FN_df.index.values) - novel_intron_ids
            )

            not_novel_accuracy_stats = calc_accuracy_stats_for_TP_FP_FN_df(
                i_sample_TP_FP_FN_df.loc[list(not_novel_intron_ids), :]
            )
            sensitivity_val = not_novel_accuracy_stats["Sensitivity"]

            F1_val = accuracy_stats_for_novel["F1"]

        name, color, line = colorAndLabel(progname)
        opacity = 1.0

        df_data.append(
            (
                progname,
                color,
                line,
                opacity,
                progname,
                "knownTPR",
                sensitivity_val,
            )
        )
        df_data.append(
            (
                progname,
                color,
                line,
                opacity,
                progname,
                "novelTPR",
                novel_trans_sensitivity,
            )
        )
        df_data.append(
            (
                progname,
                color,
                line,
                opacity,
                progname,
                "novelFDR",
                FDR_val,
            )
        )
        df_data.append(
            (
                progname,
                color,
                line,
                opacity,
                progname,
                "F1",
                F1_val,
            )
        )

    df = pd.DataFrame(
        df_data,
        columns=[
            "Name",
            "Color",
            "Line",
            "Opacity",
            "Program",
            "Metric",
            "Value",
        ],
    )

    fig, axs = plt.subplots(
        ncols=4, nrows=1, figsize=(7, 2.8), dpi=300, layout="constrained", sharey=True
    )
    acc_type = ["knownTPR", "novelTPR", "novelFDR", "F1"]
    for j in range(4):
        plt_df = df[df["Metric"] == acc_type[j]]
        plt_df["Value"] *= 100
        axs[j].barh(y=plt_df["Name"], width=plt_df["Value"], color=plt_df["Color"])
        axs[j].set_xlim(0, 100)
        axs[j].spines["right"].set_visible(False)
        axs[j].spines["top"].set_visible(False)
        var = (
            acc_type[j]
            .replace("known", "Annotated ")
            .replace("novelFDR", "FDR")
            .replace("novel", "Unannotated ")
        )
        var = var.replace("FDR", "False Discovery Rate (FDR)").replace(
            "TPR", "Sensitivity (TPR)"
        )
        axs[j].set_xlabel(var, fontsize=7)
        for bars in axs[j].containers:
            axs[j].bar_label(bars, fmt="%.1f%%")

    return df, all_TP_FP_FN_df


# report writing
def progname_to_i_sample_df_dict_to_tsv(progname_i_sample_df_dict, output_tsv_filename):

    concat_df = None

    for progname, i_sample_df in progname_i_sample_df_dict.items():
        i_sample_df2 = i_sample_df.copy()
        i_sample_df2[["progname"]] = progname
        concat_df = pd.concat([concat_df, i_sample_df2])

    concat_df.to_csv(output_tsv_filename, sep="\t", quoting=csv.QUOTE_NONE)
