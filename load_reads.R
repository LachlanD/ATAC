############################
# Load
############################
library(data.table)


#############################
# Parameters
#############################
# Use a small subset of peaks
small.set = FALSE
first.n.peaks = 1000


######################################################
# Load reads in data.frame, change to data.table? 
######################################################
#reads = read.delim('Data/MODEL_120_peaks.reads', header = TRUE, sep = '\t')
reads = fread('test.counts', header = TRUE)
peaks = fread('test.bed', header = FALSE)
colnames(peaks) = c("chr", "start", "end", "summit", "id")


###################################################
#Create a limited dataset if wanted
###################################################
if(small.set)
{
  peaks = peaks[1:first.n.peaks]
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

