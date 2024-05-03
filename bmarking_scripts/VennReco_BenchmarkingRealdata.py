#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


UTILDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../workflows/lr_isoform_custom_docker"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Experimental design: (see IsoQuant paper for details)
#
# Isoforms are reconstructed by various methods
#
# The subset of isoforms in the reference annotation that are identified as expressed: --ref_expressed_gtf
# - Reference transcripts (you decide (if reduced reference, leverage the kept set for comparison of novel)
# - Isoforms reconstructed that are entirely novel (missing from the --ref_expressed_gtf) are compared across different reco methods
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    parser = argparse.ArgumentParser(description="run benchmarking for reconstructions on real data examining known and novel isoforms reconstructed",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--ref_expressed_or_kept_gtf", type=str, required=True, help="reference isoforms to compare - you decide (if reduced reference, leverage the kept set for comparison of novel)")

    parser.add_argument("--reco_gtfs", type=str, required=True, help="transcripts reconstructed by methods", nargs='+')

    parser.add_argument("--dataset_name", type=str, required=True, help="name of dataset")
    

    args = parser.parse_args()


    ref_gtf = os.path.abspath(args.ref_expressed_or_kept_gtf)
    reco_gtfs = args.reco_gtfs
    dataset_name = args.dataset_name


    for i, gtf in enumerate(reco_gtfs):
        reco_gtfs[i] = os.path.abspath(gtf)
    
    workdir = os.getcwd()

    tool_names = list()
    tracking_files = list()
    
    for reco_gtf in reco_gtfs:
        tool_basename = os.path.basename(reco_gtf)
        tool_basename = re.sub("\\.(gtf|gff3|gff)$", "", tool_basename, flags=re.I)

        tool_names.append(tool_basename)

        
    cmd = f"gffcompare -o {tool_basename}.denovo {ref_gtf} " + " ".join(reco_gtfs)
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
            

    # SummarizeAnalysis
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "summarize_denovo_analysis.py"),
                     "--tracking " + f"{tool_basename}.denovo.tracking",
                     "--tool-names " + " ".join(tool_names),
                     "--num-tools " + str(len(tool_names)), 
                     "--dataset-name " + dataset_name ])
    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
    
    

    # PlotDenovoAnalysisSummaryKnown
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "plot_denovo_analysis_summary.py"),
                     "--input " + f"{dataset_name}_denovo_analysis_summary_known.tsv",
                     "--dataset-name " + dataset_name,
                     "--type known",
                     "--save" ])

    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
                     


    # PlotDenovoAnalysisSummaryNovel
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "plot_denovo_analysis_summary.py"),
                     "--input " + f"{dataset_name}_denovo_analysis_summary_novel.tsv",
                     "--dataset-name " + dataset_name,
                     "--type novel",
                     "--save" ])

    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
                     
    
    sys.exit(0)
                        
    



if __name__=='__main__':
    main()
