Find the footprints within called peaks.

Requires peak positions in bed file format. MACS2 calls peaks in a bed file.  

To include the summitpositions use xls_to_bed.py. 

Create a reads file with extract_reads.py.  
extract_reads.py -p [peak position bed files] -a [sam/bam alignment files].
Set -s option if summit positions are included in bed file.
Set -o with a filename to write output to file
If you wish to extract reads from a custom window width set using -w option (100).

The reads file and peaks bed file will be the input to load_reads.R. 

load_reads.R creates data structures which are loaded by footprint.R

footprint.R run an algorithm to find footprints.  
Can specify:
	Relative max peak threshold
	Absolute peak threshold
	Minimum footprint fall from peak


