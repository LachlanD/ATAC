#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args = c("Trlp_footprints.bed")

name = unlist(strsplit(args[1], "_"))[1]

footprint = read.delim(args[1], header = FALSE)

dp = duplicated(footprint$V5)

di = footprint$V5[dp] 
d = footprint$V5 %in% di

df = footprint[d,]
ndf = footprint[!d,]

offsets = (df$V7+df$V8+1)/2-df$V4
trim = abs(offsets)<=100
offsets = offsets[trim]
strand = (as.numeric(df$V11=="+")*2-1)[trim]

correct = offsets*strand 

first = correct[duplicated(df[trim,]$V5)]

png(filename = paste(name, "repeated.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(correct, breaks=-100:100, col= "red", main = paste("Repeated", name,  "offset from footprint centre\n(strand corrected)"), xlab = "offset (bp)")
hist(first, breaks = -100:100, col= "blue", add = TRUE)
legend("topright", c("First", "Second"), fill=c("blue", "red"))
invisible(dev.off())

offsetsn = (ndf$V7+ndf$V8+1)/2-ndf$V4
trimn = abs(offsetsn)<=100
offsetsn = offsetsn[trimn]
strandn = (as.numeric(ndf$V11=="+")*2-1)[trimn]

correctn = offsetsn*strandn

png(filename = paste(name, "unrepeated.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(correctn, breaks=-100:100, main = paste("Unrepeated", name, "offset from footprint centre"), xlab = "offset (bp)")
invisible(dev.off())




id = df$V5
gaps = unlist(lapply(id, function(x){diff(df[df$V5==x,]$V7)}))

pos = unlist(lapply(id, function(x){diff(df[df$V5==x & df$V11=="+",]$V7)}))

png(filename = paste(name, "repeats_gap.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(gaps, breaks = -100:100, col="red", main = paste("distances between repeats of", name), xlab = "distance (bp)")
hist(pos, breaks = -100:100, col="blue", add=TRUE)
legend("topright", c("Positive strand", "Negative strand"), fill=c("blue", "red"))
invisible(dev.off())
