#!/usr/bin/env Rscript

############################
# Load
############################
require(optparse)
require(data.table)

option_list = list(
  make_option(c("-r","--reads"), type = "character", default="test.counts", help="reads file from extract_reads.py"),
  make_option(c("-b","--bed"), type = "character", default="test.bed", help="peaks BED files"),
  make_option(c("-c", "--cut"), type = "integer", default=NA, help="cut the dataset to this number of peaks"),
  make_option(c("-s", "--summit"), action = "store_true", default = FALSE, help="Use this option if summit is specified in BED file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



######################################################
# Load reads in data.frame, change to data.table? 
######################################################
#reads = read.delim('Data/MODEL_120_peaks.reads', header = TRUE, sep = '\t')
reads = fread(opt$reads, header = TRUE)
peaks = fread(opt$bed, header = FALSE)

if(!opt$summit)
{
  colnames(peaks) = c("chr", "start", "end",  "id")
  peaks[,summit := start + (end-start)%/%2]
} else 
{
  colnames(peaks) = c("chr", "start", "end", "summit", "id") 
}


###################################################
#Create a limited dataset if wanted
###################################################
if(!is.na(opt$cut))
{
  peaks = peaks[1:opt$cut]
  reads=reads[reads$id %in% peaks$id,]
} 

peaks[,centre:=0]
peaks[,nfootprints:=0]
peaks[,nreads:=0]

###################
# Create a file
###################
save(reads, file="Data/reads_data.Rba")
save(peaks, file="Data/peaks_data.Rba")

