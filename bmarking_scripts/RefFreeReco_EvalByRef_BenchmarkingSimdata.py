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
# Isoforms are reconstructed just based on the read alignments (bam) and not leveraging any annotation as a guide.
#
# The subset of isoforms in the reference annotation that are expressed and candidates for reconstruction are provided as: --ref_expressed_gtf
# - Refernce transcripts reconstructed by the method are considered True Positives.
# - Isoforms reconstructed that are entirely novel (missing from the --ref_expressed_gtf) are treated as False Positives.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    parser = argparse.ArgumentParser(description="run benchmarking for ref-free reconstruction based on the reference transcript structures",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--ref_expressed_gtf", type=str, required=True, help="all expressed transcripts included during isoform reconstruction")

    parser.add_argument("--reco_gtfs", type=str, required=True, help="transcripts reconstructed by methods", nargs='+')

    parser.add_argument("--dataset_name", type=str, required=True, help="name of dataset")
    

    args = parser.parse_args()


    ref_expressed_gtf = os.path.abspath(args.ref_expressed_gtf)
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
        tool_basename = re.sub("\\.(gtf|gff3|gff)$", "", tool_basename, flags=re.I)

        tool_names.append(tool_basename)

        if not os.path.exists(tool_basename):
            os.makedirs(tool_basename)
            
        os.chdir(tool_basename)

        tracking_file = os.path.abspath(f"{tool_basename}.refFree")
        tracking_files.append(tracking_file)

        checkpoint = f"{tool_basename}.ok"

        if os.path.exists(tracking_file) and os.path.exists(checkpoint):
            continue
        
        cmd = f"gffcompare -r {ref_expressed_gtf} -o {tool_basename}.refFree {reco_gtf}"
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        cmd = f"touch {checkpoint}"
        subprocess.check_call(cmd, shell=True)
        

    os.chdir(workdir)

    # SummarizeAnalysis
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "summarize_reffree_analysis.py"),
                     "--input-list " + " ".join(tracking_files),
                     "--tool-names " + " ".join(tool_names),
                     "--dataset-name " + dataset_name ])
    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
    


    # PlotAnalysisSummary
    cmd = " ".join([ "python3",
                     os.path.join(UTILDIR, "plot_analysis_summary.py"),
                     "--input " + f"{dataset_name}_analysis_summary_reffree.tsv",
                     "--dataset-name " + dataset_name,
                     "--type reffree",
                     "--save" ])

    logger.info(cmd)

    subprocess.check_call(cmd, shell=True)
                     
    
    
    sys.exit(0)
                        
    



if __name__=='__main__':
    main()
