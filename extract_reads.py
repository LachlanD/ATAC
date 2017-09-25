#!/usr/bin/env python
import pysam
import re
import sys
import argparse

parser = argparse.ArgumentParser(description="usage: %prog [options] -p [peaks bed files] -a [alignment bam files] -o [optional output file] -s [if bed file contains a summit position] -w [specify a read window width]")
parser.add_argument("-p", "--peaks", nargs='+', help="specify peaks in bed files")
parser.add_argument("-a", "--alignments", nargs='+', help="specify sam/bam alignment files")
parser.add_argument("-o", "--output", nargs='?', type=argparse.FileType('w'), help="specify output file", default=sys.stdout)
parser.add_argument("-s", "--summit", action="store_true", help="specify if the bed files have a summit position")
parser.add_argument("-w", "--width", help="specify the window width for extracting reads", default=100)

args = parser.parse_args()

out = args.output

out.write("id\tchr\tpos\tstrand\tinsert\n")



for s in args.alignments:
    samfile = pysam.AlignmentFile(s, 'rb')
    refs = samfile.references
    lens = samfile.lengths
    for p in args.peaks:
        f = open(p, 'r')
        for line in f.readlines():
            sp = line.split()

            chrm = sp[0]
            if (args.summit):
                pos = int(sp[3])
                idn = sp[4]
            else:
                pos = int(sp[1])+(int(sp[2])-int(sp[1]))/2
                idn = sp[3]

            for read in samfile.fetch(chrm, max(0, pos-args.width), min(pos+args.width, lens[refs.index(chrm)])):
                if read.is_reverse:
                    strand="-"
                    start = read.reference_end-pos
                else:
                    strand = "+"
                    start = read.reference_end + 1 - pos
                out.write("%s\t%s\t%d\t%s\t%d\n" %(idn,chrm,start,strand,read.template_length)) 
        f.close()

out.close()
