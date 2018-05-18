#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args= c("M1BP_footprints.bed")

name = unlist(strsplit(args[1], "_"))[1]
footprints = read.delim(args[1], header = FALSE)

offsets = (footprints$V7+footprints$V8+1)/2-footprints$V4
trim = abs(offsets)<=100
offsets = offsets[trim]

png(filename = paste(name, "_uncorrected.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(offsets, breaks =-100:100, main = paste(name, "motif offset from footprint centre"), xlab="offset (bp)")
invisible(dev.off())


strand = (as.numeric(footprints$V11=="+")*2-1)[trim]
correct = (offsets*strand)
pos = offsets[footprints$V11=="+"]

png(filename = paste(name, "_corrected.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(correct, breaks =-100:100, col = "red", main = paste(name, "motif offset from footprint centre \n (strand corrected)"), xlab="offset (bp)")
hist(pos, breaks =-100:100, col = "blue", add = TRUE)
legend("topright", c("Positive strand", "Negative strand"), fill=c("blue", "red"))
invisible(dev.off())

fp_length = footprints$V3-footprints$V2

png(filename = paste(name, "_footprint_length.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(fp_length, breaks= 0:200, main = paste(name, "footprint length"), xlab = "length (bp)")
invisible(dev.off())