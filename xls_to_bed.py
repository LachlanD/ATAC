#!/usr/bin/env python

import xlrd
import sys

xl = open(sys.argv[1], 'r')
bed = open(sys.argv[2], 'w')

for line in xl:
  if line[0:3] == "chr":
    l = line.split()
    if l[1].isdigit():
      bed.write(l[0] + '\t' + l[1] + '\t' + l[2] + '\t' + l[4] + '\t' + l[9] + '\n')
    




