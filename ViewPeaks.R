###########################################
#Visualise loop for debugging and tuning
##########################################

###################
# Tuning Parameters
###################

# Minimum read coverage for a peak to be considered  
min.reads = 200

# Threshold relative to max
rel.threshold = 0.5

# Absolute threshold
abs.threshold = 0.0

# Footprint must fall this far below the smaller of the two peaks
footprint.fraction = 0.9


#############################
#Setup
#############################
library(data.table)
source("ATAC_functions.R")
load("Data/reads_data.Rba")
ids = unique(reads$id)

for (i in 1:length(ids))
{
  r = reads[reads$id==ids[i],]
  if(nrow(r) > 100)
  {
    d = density(r$pos)
    pos = density(r[r$strand=='+',]$pos)
    neg = density(r[r$strand=='-',]$pos)
    
    plot(d,col='blue', ylim=c(0,max(pos$y,pos$y,neg$y)*1.1), xlim=c(-200,200))
    lines(pos,col='green')
    lines(neg,col='red')
  }
  
  threshold = max(c(max(d$y)*rel.threshold , abs.threshold))
  
  peaks = sort(find_peaks(d$y, threshold, footprint.fraction))
  for(p in peaks)
  {
    abline(v= d$x[p], col='yellow')
  }
  
  footprints = find_footprints(d$y, peaks)
  
  for(f in footprints)
  {
    abline(v=d$x[f], col='orange')
  }
  
  centre = find_centre(peaks, footprints)
  
  abline(v=d$x[centre], col='purple')
  
  invisible(readline(prompt="Press [enter] to continue"))
}
