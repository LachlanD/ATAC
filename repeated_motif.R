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

offset = (df$V7+df$V8+1)/2-df$V4
trim = abs(offsets)<=100
offsets = offsets[trim]
strand = (as.numeric(df$V11=="+")*2-1)[trim]

correct = offset*strand 

png(filename = paste(name, "repeated.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(correct, breaks=-100:100, main = paste("Repeated", name,  "offset from footprint centre"), xlab = "offset (bp)")
invisible(dev.off())

offsetn = (ndf$V7+ndf$V8+1)/2-ndf$V4
trimn = abs(offsetsn)<=100
offsetsn = offsetsn[trimn]
strandn = (as.numeric(ndf$V11=="+")*2-1)[trimn]

correctn = offsetn*strandn

png(filename = paste(name, "unrepeated.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(correctn, breaks=-100:100, main = paste("Unrepeated", name, "offset from footprint centre"), xlab = "offset (bp)")
invisible(dev.off())




id = df$V5
gaps = lapply(id, function(x){diff(df[df$V5==x,]$V7)})

png(filename = paste(name, "repeats_gap.png" , sep = "_"), width = 1920, height = 1080, pointsize = 24)
hist(unlist(gaps), breaks = -100:100, main = paste("distances between repeats of", name), xlab = "distance (bp)")
invisible(dev.off())
