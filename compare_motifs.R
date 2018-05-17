#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args = c("tll_footprints.bed", "MSC_footprints.bed")

name1 = unlist(strsplit(args[1], "_"))[1]
name2 = unlist(strsplit(args[2], "_"))[1]

footprint1 = read.delim(args[1], header = FALSE)
footprint2 = read.delim(args[2], header = FALSE)

merged = merge(footprint1, footprint2, by="V5")

offsets1 = (merged$V7.x+merged$V8.x+1)/2-merged$V4.x
#hist(offsets, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre", sep =" "), xlab="offset (bp)")

strand1 = as.numeric(merged$V11.x=="+")*2-1
correct1 = offsets1*strand1

png(filename = paste(name1, "_with_", name2, ".png", sep = ""))
hist(correct1, xlim = c(-100,100), breaks =200, main = paste(name1, "motif offset from footprint centre \n when ", name2, " is present\n (strand corrected)", sep =" "), xlab="offset (bp)")
invisible(dev.off())

offsets2 = (merged$V7.y+merged$V8.y+1)/2-merged$V4.y
#hist(offsets2, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre", sep =" "), xlab="offset (bp)")

strand2 = as.numeric(merged$V11.y=="+")*2-1
correct2 = offsets2*strand2

png(filename = paste(name2, "_with_", name1, ".png", sep = ""))
hist(correct2, xlim = c(-100,100), breaks =200, main = paste(name2, "motif offset from footprint centre \n when ", name1, " is present\n (strand corrected)", sep =" "), xlab="offset (bp)")
invisible(dev.off())

gap = (merged$V7.x+merged$V8.x+1)/2-(merged$V7.y+merged$V8.y+1)/2
diff = as.numeric(merged$V11.x=="+")*2-1

png(paste(name1,"_", name2, "_gap_uncorrected.png", sep = ""))
hist(gap, breaks =200, xlim = c(-100,100), main = paste(name1, " to ", name2, " gap"), xlab="bp gap")
invisible(dev.off())

png(paste(name1,"_", name2, "_gap_corrected.png", sep = ""))
hist(gap*diff, breaks =200, xlim = c(-100,100), main = paste(name1, " to ", name2, " gap\n(strand corrected)"), xlab="bp gap")
invisible(dev.off())

fp_length = merged$V3.x-merged$V2.x

png(paste(name1, "_", name2,"_footprint_length.png", sep=""))
hist(fp_length, breaks= 200, xlim=c(0,200), main = paste(name1, name2, "shared footprint length"), xlab = "lenght (bp)")
invisible(dev.off())

fp1_only = footprint1[!(footprint1$V5 %in% merged$V5),]
fp2_only = footprint2[!(footprint2$V5 %in% merged$V5),]



fp1_length = fp1_only$V3-fp1_only$V2

png(paste(name1, "_without_", name2,"_footprint_length.png", sep=""))
hist(fp1_length, breaks= 200, xlim=c(0,200), main = paste(name1, "without", name2, "footprint length"), xlab = "lenght (bp)")
invisible(dev.off())

fp2_length = fp2_only$V3-fp2_only$V2

png(paste(name2, "_without_", name1,"_footprint_length.png", sep=""))
hist(fp2_length, breaks= 200, xlim=c(0,200), main = paste(name2, "without", name1, "footprint length"), xlab = "lenght (bp)")
invisible(dev.off())