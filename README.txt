Find the footprints within called peaks.

Requires peak positions in bed file format. MACS2 calls peaks in a bed file.  

To include the summitpositions use xls_to_bed.py. 

Create a reads file with extract_reads.py.  
extract_reads.py -p [peak position bed files] -a [sam/bam alignment files].
Set -s option if summit positions are included in bed file.
Set -o with a filename to write output to file
If you wish to extract reads from a custom window width set using -w option (100).


