#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import random

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="extract GTF records according to benchmarking categories",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--gtf", type=str, required=True, help="gtf input file")
    parser.add_argument("--frac_keep", type=float, required=True, help="fraction of transcripts to keep")
    parser.add_argument("--output_prefix", type=str, required=True, help="prefix for .kept.gtf and .excluded.gtf files")
    
    args = parser.parse_args()

    gtf_input_file = args.gtf
    frac_keep = args.frac_keep
    output_prefix = args.output_prefix
    

    transcript_to_gtf_txt = parse_transcripts_from_gtf(gtf_input_file)

    kept_output_gtf_file = output_prefix + ".kept.gtf"
    excluded_output_gtf_file= output_prefix + ".excluded.gtf"

    ofh_kept = open(kept_output_gtf_file, "wt")
    ofh_excluded = open(excluded_output_gtf_file, "wt")

    for transcript_id, gtf_txt in transcript_to_gtf_txt.items():
        
        ofh = ofh_kept if random.random() < frac_keep else ofh_excluded
        print(gtf_txt, end='', file=ofh)

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
