#!/bin/bash

# search_motif.sh <footprints bed file *.bed> <window size fraction or fixed> <reference genome .fa/fasta> <homer motif set eg. insect> <number of cpus> 

# eg. >$ ./search_motif.sh footprint_positions.bed 60 dm6.fa insect 4 

#require homer and python 3+

window=$2"_window.bed"

python binding_window.py $1 $2 > "$window" 
findMotifsGenome.pl "$window" $3 $2"window_homer_output" -mset $4 -size given -p $5
