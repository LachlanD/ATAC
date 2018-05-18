#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args = c("tll_footprints.bed", "MSC_footprints.bed")

name1 = unlist(strsplit(args[1], "_"))[1]
name2 = unlist(strsplit(args[2], "_"))[1]

footprint1 = read.delim(args[1], header = FALSE)
footprint2 = read.delim(args[2], header = FALSE)

merged = merge(footprint1, footprint2, by="V5")

offsets1 = (merged$V7.x+merged$V8.x+1)/2-merged$V4.x
trim1 = abs(offsets1)<=100
offsets1 = offsets1[trim1]
#hist(offsets, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre", sep =" "), xlab="offset (bp)")

strand1 = (as.numeric(merged$V11.x=="+")*2-1)[trim1]
correct1 = (offsets1*strand1)
pos1 = correct1[merged$V11.x=="+"]

png(filename = paste(name1, "_with_", name2, ".png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(correct1, col = "red", breaks =-100:100, main = paste(name1, "motif offset from footprint centre\nwhen", name2, "is present\n(strand corrected)"), xlab="offset (bp)")
hist(pos1, col = "blue", breaks = -100:100, add= TRUE)
legend("topright", c("Positive strand", "Negative strand"), fill=c("blue", "red"))
invisible(dev.off())

offsets2 = (merged$V7.y+merged$V8.y+1)/2-merged$V4.y
trim2 = abs(offsets2)<=100
offsets2 = offsets2[trim2]
#hist(offsets2, xlim = c(-100,100), breaks =200, main = paste(name, "motif offset from footprint centre", sep =" "), xlab="offset (bp)")

strand2 = (as.numeric(merged$V11.y=="+")*2-1)[trim2]
correct2 = (offsets2*strand2)
pos2 = correct2[merged$V11.y=="+"]

png(filename = paste(name2, "_with_", name1, ".png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(correct2, breaks =-100:100, col = "red", main = paste(name2, "motif offset from footprint centre\n when", name1, "is present\n (strand corrected)"), xlab="offset (bp)")
hist(pos2, col = "blue", breaks = -100:100, add= TRUE)
legend("topright", c("Positive strand", "Negative strand"), fill=c("blue", "red"))
invisible(dev.off())

gap = (merged$V7.x+merged$V8.x+1)/2-(merged$V7.y+merged$V8.y+1)/2
gap = gap[abs(gap)<=100]
diff = as.numeric(merged$V11.x=="+")*2-1
pos_g = gap[merged$V11.x=="+"]

png(paste(name1,"_", name2, "_gap_uncorrected.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(gap, breaks =-100:100, main = paste(name1, "to", name2, "motif centres distance"), xlab="distance (bp)")
invisible(dev.off())

png(paste(name1,"_", name2, "_gap_corrected.png", sep = ""), width = 1920, height = 1080, pointsize = 24)
hist(gap*diff, breaks =-100:100, col= "red", main = paste(name1, "to", name2, "motif centres distance\n(strand corrected)"), xlab="distance (bp)")
hist(pos_g, breaks = -100:100, col="blue", add = TRUE)
legend("topright", c("Positive strand", "Negative strand"), fill=c("blue", "red"))
invisible(dev.off())

fp_length = merged$V3.x-merged$V2.x

png(paste(name1, "_", name2,"_footprint_length.png", sep=""), width = 1920, height = 1080, pointsize = 24)
hist(fp_length, breaks= 0:200, main = paste(name1, name2, "shared footprint lengths"), xlab = "length (bp)")
invisible(dev.off())

fp1_only = footprint1[!(footprint1$V5 %in% merged$V5),]
fp2_only = footprint2[!(footprint2$V5 %in% merged$V5),]



fp1_length = fp1_only$V3-fp1_only$V2

png(paste(name1, "_without_", name2,"_footprint_length.png", sep=""), width = 1920, height = 1080, pointsize = 24)
hist(fp1_length, breaks= 0:200, main = paste(name1, "without", name2, "footprint lengths"), xlab = "length (bp)")
invisible(dev.off())

fp2_length = fp2_only$V3-fp2_only$V2

png(paste(name2, "_without_", name1,"_footprint_length.png", sep=""), width = 1920, height = 1080, pointsize = 24)
hist(fp2_length, breaks= 0:200, main = paste(name2, "without", name1, "footprint lengths"), xlab = "length (bp)")
invisible(dev.off())