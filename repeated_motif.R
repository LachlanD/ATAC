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
strand = as.numeric(df$V11=="+")*2-1

correct = offset*strand 

png(filename = paste(name, "repeated.png" , sep = "_"))
hist(correct, breaks=200, xlim=c(-100,100), main = paste("Repeated", name,  "offset from footprint centre"), xlab = "offset (bp)")
invisible(dev.off())

offsetn = (ndf$V7+ndf$V8+1)/2-ndf$V4
strandn = as.numeric(ndf$V11=="+")*2-1

correctn = offsetn*strandn

png(filename = paste(name, "unrepeated.png" , sep = "_"))
hist(correctn, breaks=200, xlim=c(-100,100), main = paste("Unrepeated", name, "offset from footprint centre"), xlab = "offset (bp)")
invisible(dev.off())




id = df$V5
gaps = lapply(id, function(x){diff(df[df$V5==x,]$V7)})

png(filename = paste(name, "repeats_gap.png" , sep = "_"))
hist(unlist(gaps), breaks = 200, xlim = c(-100,100), main = paste("gap between repeats of", name), xlab = "gap (bp)")
invisible(dev.off())
