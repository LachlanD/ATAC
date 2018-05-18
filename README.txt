ATACFactor

Footprint Calling


Find the footprints within called peaks.

Requires peak positions in bed file format. MACS2 calls peaks in a bed file.  

To include the summit positions use xls_to_bed.py to convert the MACS peaks xls file into a bed file. 

Create a reads file with extract_reads.py.  
extract_reads.py -p [peak position BED files] -a [sam/bam alignment files].
Set -s option if summit positions are included in bed file.
Set -o with a filename to write output to file
If you wish to extract reads from a custom window width set using -w option (100).

The reads file and peaks bed file will be the input to load_reads.R. 

load_reads.R -r [reads file] -b [peaks BED file]
Set -s if the BED file includes summit positions.
Use -c [int n] to cut the data to the first n peaks.

load_reads.R creates data structures which are loaded by footprint.R.

footprint.R run an algorithm to find footprints.  
Footprint is defined as the minimum position between 2 summits.  Summits are defined as a local maximum that is above a threshold and that rises far enough above surrounding local minimums (prominence).

Can specify:
	Minimum number of reads for the peak to be considered.
	Summit threshold relative to max.
	Absolute summit threshold.
	Minimum prominence.

A footprints bed file wil be produced with the position of each footprint.  The R data structures are also saved to allow you to explore your results.




Data Analysis


Using Homer to find motifs which are over-represented in footprints.

Note: Requires Homer, bedtools, python 3+, R.

Choose a window from the centre to search.  Can either be a fraction 0.0 to 1.0 or an integer > 1.  Integers will give a fixed window size.  A fraction will search that proportion of the footprint to neighbouring peak.

Run the search_motifs.sh script with the following arguments:
1. The footprints bed file
2. Your selected window size
3. The reference genome in fasta format
4. The set of known motifs to use <vertebrates|insects|worms|plants|yeast|all>
5. Number of cpu cores to utilise

This will produce an output folder where the motifs can be browsed.

To investigate a motif, save the motif matrix *.motif file from homer.

Run the script find_motif.sh with the following arguments:
1. The motif matrix file from homer
2. The reference genome in fasta format
3. The footprints bed file
4. The number of cpu cores to utilise

This script will out put a number of images showing the relationship between motif and the footprints that it occupies. A <motif>_footprints.bed file is also output.

To investigte interactions between motifs run the compare_motifs.R script with the 2 <motif>_footprints.bed files as arguments.
This will again output a number of images showing the relationship between the 2 motifs.

To investigate whether a motif is repeated run the repeated_motif.R with the <motif>_footprints.bed file as the argument.
