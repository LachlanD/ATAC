#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args= c("MSC_footprints.bed", "90")


split_length = as.integer(args[2])

name = unlist(strsplit(args[1], "_"))[1]
footprints = read.delim(args[1], header = FALSE)

offsets = (footprints$V7+footprints$V8+1)/2-footprints$V4
trim = abs(offsets)<=100
offsets = offsets[trim]

lengths = footprints$V3-footprints$V2
lenghts = lengths[trim]

strand = (as.numeric(footprints$V11=="+")*2-1)[trim]
correct = (offsets*strand)
split = correct[lengths>=split_length]

png(filename = paste(name, "_split_length.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(correct, breaks =-100:100, col = "red", main = paste(name, "motif offset from footprint centre\n (strand corrected)"), xlab="offset (bp)")
hist(split, breaks =-100:100, col = "blue", add = TRUE)
legend("topright", c(paste("footprint length >=",split_length), paste("footprint length <",split_length)), fill=c("blue", "red"))
invisible(dev.off())
