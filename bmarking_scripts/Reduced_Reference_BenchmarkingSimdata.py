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
# A reduced reference gtf (by 20%) is used during isoform reconstruction by the methods
# In the full reference gtf, the subset of isoforms that are expressed (have read support) are set to --ref_included_gtf
# The subset of the expressed transcrips that were excluded via the reduced refernece gtf are set to --ref_excluded_gtf
# - Those excluded expressed transcripts found reconstructed by the method are considered True Positives.
# - Isoforms reconstructed that are entirely novel (missing from the --ref_included_gtf) are treated as False Positives.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#################################################
## Definitions for categories:
#
#             in_ref_expr  in_ref_expr_kept  reco
# novel_FN:      +                 -          -
#
# known_FN:      +                 +          -
#
# known_TP:      +                 +          +
#
# novel_TP:      +                 -          +
#
# novel_FP:      -                 -          +
#
#################################################



def main():
    parser = argparse.ArgumentParser(description="run benchmarking for ref-based reconstruction where a subset of reference transcripts were withheld",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--ref_expressed_gtf", type=str, required=True, help="all expressed transcripts included during isoform reconstruction")
    parser.add_argument("--ref_expressed_kept_gtf", type=str, required=True, help="the expressed transcripts excluded during isoform reconstruction and examined here for benchmarking accuracy")

    parser.add_argument("--reco_gtfs", type=str, required=True, help="transcripts reconstructed by methods", nargs='+')

    parser.add_argument("--dataset_name", type=str, required=True, help="name of dataset")
    

    args = parser.parse_args()


    ref_expressed_gtf = os.path.abspath(args.ref_expressed_gtf)
    ref_expressed_kept_gtf = os.path.abspath(args.ref_expressed_kept_gtf)
    reco_gtfs = args.reco_gtfs
    dataset_name = args.dataset_name


    for i, gtf in enumerate(reco_gtfs):
        reco_gtfs[i] = os.path.abspath(gtf)
    
    workdir = os.getcwd()

    tool_names = list()
    tracking_files = list()
    
    for reco_gtf in reco_gtfs:
        os.chdir(workdir)
        tool_basename = os.path.basename(reco_gtf)
        tool_basename = re.sub("\\.(gtf|gff3)$", "", tool_basename, flags=re.I)

        tool_names.append(tool_basename)

        if not os.path.exists(tool_basename):
            os.makedirs(tool_basename)
            
        os.chdir(tool_basename)

        tracking_file = os.path.abspath(f"{tool_basename}.tracking")
        tracking_files.append(tracking_file)

        checkpoint = f"{tool_basename}.ok"

        if os.path.exists(tracking_file) and os.path.exists(checkpoint):
            continue
        
        cmd = f"gffcompare -o {tool_basename} {ref_expressed_gtf} {ref_expressed_kept_gtf} {reco_gtf}"
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        cmd = f"touch {checkpoint}"
        subprocess.check_call(cmd, shell=True)
        

    os.chdir(workdir)

    # SummarizeAnalysis
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "summarize_analysis.py"),
                     "--tracking " + " ".join(tracking_files),
                     "--tool-names " + " ".join(tool_names),
                     "--dataset-name " + dataset_name ])
    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
    


    # PlotAnalysisSummary
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "plot_analysis_summary.py"),
                     "--input " + f"{dataset_name}_analysis_summary.tsv",
                     "--dataset-name " + dataset_name,
                     "--type reduced",
                     "--save" ])

    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
                     
    
    
    sys.exit(0)
                        
    



if __name__=='__main__':
    main()
