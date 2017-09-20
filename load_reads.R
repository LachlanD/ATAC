############################
# Load
############################
library(data.table)


#############################
# Parameters
#############################
# Use a small subset of peaks
small.set = TRUE
first.n.peaks = 100

######################################################
# Load reads in data.frame, change to data.table? 
######################################################
#reads = read.delim('Data/MODEL_120_peaks.reads', header = TRUE, sep = '\t')
reads = fread('Data/MODEL_120_peaks.reads', header = TRUE, sep = '\t')
ids = unique(reads$id)

###################################################
#Create a limited dataset if wanted
###################################################
if(small.set)
{
  ids = ids[1:first.n.peaks]
  reads=reads[reads$id %in% ids,]
}

###################
# Create a file
###################
save(reads, file="Data/reads_data.Rba")
