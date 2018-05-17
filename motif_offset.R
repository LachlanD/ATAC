#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args= c("M1BP_footprints.bed")

name = unlist(strsplit(args[1], "_"))[1]
footprints = read.delim(args[1], header = FALSE)

offsets = (footprints$V7+footprints$V8+1)/2-footprints$V4

png(filename = paste(name, "_uncorrected.png", sep = ""))
hist(offsets, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre"), xlab="offset (bp)")
invisible(dev.off())

strand = as.numeric(footprints$V11=="+")*2-1
correct = offsets*strand

png(filename = paste(name, "_corrected.png", sep = ""))
hist(correct, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre \n (strand corrected)"), xlab="offset (bp)")
invisible(dev.off())

fp_length = footprints$V3-footprints$V2

png(filename = paste(name, "_footprint_length.png", sep = ""))
hist(fp_length, breaks= 200, xlim = c(0, 200), main = paste(name, "footprint length"), xlab = "length (bp)")
invisible(dev.off())