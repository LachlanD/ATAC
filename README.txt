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


