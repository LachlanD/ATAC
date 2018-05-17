#!/bin/bash

# find_motif.sh <homer motif file *.motif> <reference genome *.fa/fasta> <footprint positions *.bed> <number of cpus>

# requires homer and bedtools

all_output=$(basename $1 .motif)"_position.bed"
footprint_output=$(basename $1 .motif)"_footprints.bed" 

scanMotifGenomeWide.pl $1 $2 -bed -p $4 > "$all_output"

bedtools intersect -a $3 -b "$all_output" -wa -wb > "$footprint_output"

./motif_offset.R $footprint_output
