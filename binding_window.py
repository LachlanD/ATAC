import math
import sys

#Usage: python binding_window.py <footprints bed file *.bed> <window size 0-1 fraction of peak-to-peak footprint length, >1 fixed window length> 

f = open(sys.argv[1], 'r') 

window = float(sys.argv[2])

if (window > 1):
    # Fixed window
    window = int((window+1)/2)

    
    for line in f:
        l = line.split()
        print(l[0] + "\t" + str(int(l[3])-window) + "\t" + str(int(l[3])+window) + "\t" + l[3] + "\t" + l[4])
else:
    # footprint length fraction
    for line in f:
        l = line.split()
        left = math.floor(((int(l[1])-int(l[3]))*window)+int(l[3]))
        right = math.ceil(((int(l[2])-int(l[3]))*window)+int(l[3]))
        print(l[0] + "\t" + str(left) + "\t" + str(right) + "\t" + l[3] + "\t" + l[4])

