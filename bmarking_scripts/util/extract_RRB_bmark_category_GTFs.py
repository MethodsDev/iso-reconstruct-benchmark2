#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="extract GTF records according to benchmarking categories",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--reco_gtf", type=str, required=True, help="reco gtf file")
    parser.add_argument("--truth_gtf", type=str, required=True, help="reco gtf file")

    parser.add_argument("--tracking", type=str, required=True, help="tracking file from gffcompare")

    args = parser.parse_args()
    
    reco_gtf_file = args.reco_gtf
    truth_gtf_file = args.truth_gtf
    
    input_tracking_file = args.tracking
    

    reco_transcript_to_gtf_txt = parse_transcripts_from_gtf(reco_gtf_file)
    
    truth_transcript_to_gtf_txt = parse_transcripts_from_gtf(truth_gtf_file)
    
    base_gtf_name = os.path.basename(reco_gtf_file)

    types = ('novel_fn', 'known_fn', 'known_tp', 'novel_tp', 'novel_fp')

    fhs = dict()
    for ctype in types:
        fhs[ctype] = open(base_gtf_name + "." + ctype + ".gtf", "wt")


    ofh_audit = open(base_gtf_name + ".refcompare-audit.txt", "wt")
    
        
    # code chunk adapted from isoquant bmarking below:
    for line in open(input_tracking_file, "r"):
        # column[0]: unique internal id for the transfrag
        # column[1]: unique internal id for the super-locus containing these transcripts across all samples and the reference annotation
        # column[2]: gene name and transcript id of the reference record associated to this transcript
        # column[3]: type of overlap or relationship between the reference transcripts and the transcript structure represented by this row
        # columns[4:]: each following column showns the transcript for each sample/tool
        transcript_columns = line.strip().split()

        
        reco_transcript_id = None
        if transcript_columns[6] != '-':
            reco_transcript_id = transcript_columns[6].split("|")[1]


        ref_gene_id = None
        ref_transcript_id = None
        if transcript_columns[4] != '-':
            ref_transcript_id = transcript_columns[4].split("|")[1]
            ref_gene_id = transcript_columns[4].split("|")[0].split(":")[1]
            
        vals = [ref_gene_id, ref_transcript_id, reco_transcript_id]
            
        if transcript_columns[4] != '-' and transcript_columns[5] == '-' and transcript_columns[6] == '-':
            #novel_fn += 1
            print(truth_transcript_to_gtf_txt[ref_transcript_id], file=fhs['novel_fn'])
            vals.insert(0, "novel_fn")
            
        elif transcript_columns[4] != '-' and transcript_columns[5] != '-' and transcript_columns[6] == '-':
            #known_fn += 1
            print(truth_transcript_to_gtf_txt[ref_transcript_id], file=fhs['known_fn'])
            vals.insert(0, "known_fn")
            
        elif transcript_columns[4] != '-' and transcript_columns[5] != '-' and transcript_columns[6] != '-':
            #known_tp += 1
            print(reco_transcript_to_gtf_txt[reco_transcript_id], file=fhs['known_tp'])
            vals.insert(0, "known_tp")
                        
        elif transcript_columns[4] != '-' and transcript_columns[5] == '-' and transcript_columns[6] != '-':
            #novel_tp += 1
            print(reco_transcript_to_gtf_txt[reco_transcript_id], file=fhs['novel_tp'])
            vals.insert(0, "novel_tp")
            
        elif transcript_columns[4] == '-' and transcript_columns[5] == '-' and transcript_columns[6] != '-':
            #novel_fp += 1
            print(reco_transcript_to_gtf_txt[reco_transcript_id], file=fhs['novel_fp'])
            vals.insert(0, "novel_fp")
            
        else:
            print("WARNING: This should not have happened! Current line: " + str(transcript_columns))



        print("\t".join([str(x) for x in vals]), file=ofh_audit)


    sys.exit(0)


def parse_transcripts_from_gtf(gtf_file):

    transcript_to_gtf_txt = defaultdict(str)
    
    with open(gtf_file) as fh:
        for line in fh:
            m = re.search("transcript_id \"([^\"]+)\"", line)
            if m:
                transcript_id = m.group(1)
                transcript_to_gtf_txt[transcript_id] += line
                
    return transcript_to_gtf_txt

    

if __name__=='__main__':
    main()
