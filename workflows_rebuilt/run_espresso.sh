#!/bin/bash

# Parameters
input_sam=$1
genome=$2
annotation=$3
out_prefix=$4

# Create output directory
mkdir -p ${out_prefix}/ID/

# Define espresso path
espresso=/home/jupyter/projects/iso_quantification_benchmark/sirvs_from_raw/processed_reads/Espresso_OUT/espresso/src/

# Create espresso samples file
echo -e "${input_sam}\tespresso" > ${out_prefix}/ID/espresso_samples.tsv

# Run ESPRESSO_S.pl
perl ${espresso}ESPRESSO_S.pl --sort_buffer_size 64 -L ${out_prefix}/ID/espresso_samples.tsv -F ${genome} -A ${annotation} -O ${out_prefix}/ID -T 8

# Run ESPRESSO_C.pl
perl ${espresso}ESPRESSO_C.pl --sort_buffer_size 64 -I ${out_prefix}/ID -F ${genome} -X 0 -T 8

# Run ESPRESSO_Q.pl
perl ${espresso}ESPRESSO_Q.pl -L ${out_prefix}/ID/espresso_samples.tsv.updated -A ${annotation} -T 8

# Move output file
mv ${out_prefix}/ID/espresso_samples_N2_R0_updated.gtf ${out_prefix}/ID/ESPRESSO_reduced.gtf