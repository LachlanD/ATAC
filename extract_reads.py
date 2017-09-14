#!/usr/bin/env python
import pysam
import re
import sys

samfile1 = pysam.AlignmentFile("JV-22-PE_adapter_trim_s_120_MD.bam", 'rb')
samfile2 = pysam.AlignmentFile("JV-24-PE_adapt_trim_s_120_MD.bam", 'rb') 
samfile3 = pysam.AlignmentFile("JV-28-PE_cutadapt_s_120_MD.bam", 'rb')

print "id\tchr\tpos\tstrand\tinsert"

peaks = open(sys.argv[1])

refs = samfile1.references
lens = samfile1.lengths



for m in peaks.readlines():
  sp = m.split()

  if (len(sp)<5):
    break

  pos = int(sp[3])
  chr = sp[0]
  
  refs = samfile1.references
  lens = samfile1.lengths

  for read in samfile1.fetch(chr, max(0,pos-100),min(pos+100, lens[refs.index(chr)])):
    if read.is_reverse:
      strand= "-"
      start = read.reference_end -pos

    else:
      strand = "+"
      start = read.reference_start+1 -pos
    print sp[4] + "\t" + chr + " \t%d\t%s\t%d" %(start, strand, read.template_length)
  
  refs = samfile2.references
  lens = samfile2.lengths

  for read in samfile2.fetch(chr, max(0,pos-100),min(pos+100, lens[refs.index(chr)])):
    if read.is_reverse:
      strand= "-"
      start = read.reference_end -pos

    else:
      strand = "+"
      start = read.reference_start+1 -pos
  
    print sp[4] + "\t" + chr + " \t%d\t%s\t%d" %(start, strand, read.template_length)

  refs = samfile3.references
  lens = samfile3.lengths

  for read in samfile3.fetch(chr, max(0,pos-100),min(pos+100, lens[refs.index(chr)])):
    if read.is_reverse:
      strand= "-"
      start = read.reference_end -pos

    else:
      strand = "+"
      start = read.reference_start+1 -pos

    print sp[4] + "\t" + chr + " \t%d\t%s\t%d" %(start, strand, read.template_length)
