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
    

    reco_transcript_to_gtf_txt, reco_transcript_to_gene_id = parse_transcripts_from_gtf(reco_gtf_file)
    
    truth_transcript_to_gtf_txt, truth_transcript_to_gene_id = parse_transcripts_from_gtf(truth_gtf_file)
    
    base_gtf_name = os.path.basename(reco_gtf_file)


    ref_trans_id_seen = set()
    reco_trans_id_seen = set()


    results = []

    types = set()

    
    # code chunk adapted from isoquant bmarking below:
    for line in open(input_tracking_file, "r"):

        """
        0       TCONS_00000001
        1       XLOC_000001
        2       SIRV1|SIRV108
        3       =
        4       q1:g:SIRV1:+:comp-1|t:SIRV1:+:comp-1:iso-2|3|0.000000|0.000000|0.000000|702
        """
        
        transcript_columns = line.strip().split()

        compare_code = transcript_columns[3]
        
        reco_transcript_id = None
        if transcript_columns[4] != '-':
            reco_transcript_id = transcript_columns[4].split("|")[1]


        ref_gene_id = None
        ref_transcript_id = None
        if transcript_columns[2] != '-':
            ref_transcript_id = transcript_columns[2].split("|")[1]
            ref_gene_id = transcript_columns[2].split("|")[0]


            
        vals = ['known_tp', ref_gene_id, ref_transcript_id, reco_transcript_id, compare_code]

        results.append(vals)

        types.add('known_tp-' + compare_code)
        
        ref_trans_id_seen.add(ref_transcript_id)
        reco_trans_id_seen.add(reco_transcript_id)

        
    

    for ref_trans_id in truth_transcript_to_gtf_txt:
        if ref_trans_id not in ref_trans_id_seen:
            results.append( ['known_fn', truth_transcript_to_gene_id[ref_trans_id], ref_trans_id, None, None] )
            types.add('known_fn')
            
    for reco_trans_id in reco_transcript_to_gtf_txt:
        if reco_trans_id not in reco_trans_id_seen:
            results.append(['novel_fp', None, None, reco_trans_id, None]) 
            types.add('novel_fp')

    
    fhs = dict()
    for ctype in types:
        fhs[ctype] = open(base_gtf_name + "." + ctype + ".gtf", "wt")


    ofh_audit = open(base_gtf_name + ".refcompare-audit.txt", "wt")

    for result in results:
        (ctype, ref_gene_id, ref_trans_id, reco_trans_id, class_code) = result

        gtf = None
         
        if ctype in ('known_tp', 'novel_fp'):
            gtf = reco_transcript_to_gtf_txt[reco_trans_id]
        else:
            gtf = truth_transcript_to_gtf_txt[ref_trans_id]
        
        if class_code is not None:
            ctype += "-" + class_code
             
        print(gtf, file=fhs[ctype])
             
        print("\t".join([str(x) for x in result]), file=ofh_audit)

         


         
    sys.exit(0)


def parse_transcripts_from_gtf(gtf_file):

    transcript_to_gtf_txt = defaultdict(str)
    transcript_to_gene_id = dict()

    
    with open(gtf_file) as fh:
        for line in fh:
            m = re.search("transcript_id \"([^\"]+)\"", line)
            if m:
                transcript_id = m.group(1)
                transcript_to_gtf_txt[transcript_id] += line

                m = re.search("gene_id \"([^\"]+)\"", line)
                gene_id = m.group(1)
                transcript_to_gene_id[transcript_id] = gene_id
                
    return transcript_to_gtf_txt, transcript_to_gene_id

    

if __name__=='__main__':
    main()
