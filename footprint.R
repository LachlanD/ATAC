#########################
# Load
#########################
library(parallel)
library(data.table)
source("ATAC_functions.R")
load("Data/reads_data.Rba")
load("Data/peaks_data.Rba")


###################
# Parameters
###################

# Minimum read coverage for a peak to be considered  
min.reads = 200

# Threshold relative to max
rel.threshold = 0.5

# Absolute threshold
 abs.threshold = 0.0

# Footprint must fall this far below the smaller of the two peaks
footprint.fraction = 0.9

# Distributed execution (Needed for Windows)
dist_mem = TRUE

# number of cores to use
no_cores = 0
if (no_cores==0)
{
  no_cores = detectCores(logical = FALSE, all.tests = T)
}


################
# Run
################
ids = peaks$id

if(dist_mem)
{
  #####################################################
  #Distributed memory !!!Use for multicore Windows!!!
  #####################################################
  cl = makeCluster(no_cores)
  
  # Find the footprints and centres
  clusterExport(cl=cl, list("process_peak", "find_peaks", "find_footprints", "find_centre", "reads", "peaks", "min.reads", "rel.threshold", "abs.threshold", "footprint.fraction"))
  clusterEvalQ(cl, library(data.table))
  fp.dist = parLapply(cl, ids, function(x){process_peak(x, reads, peaks, min.reads = min.reads, rel.threshold = rel.threshold, abs.threshold = abs.threshold, footprint.frac = footprint.fraction )})
  fp = do.call("rbind", fp.dist)
  fp = as.data.table(fp[complete.cases(fp),])
  
  parLapply(cl, ids, function(x){reads[id==x, cent.pos := pos-peaks[id==x,centre]]})
  
  stopCluster(cl=cl)

} else 
{
  ################################################################
  # Undistribute memory Versions  !!! Windows can only use 1 !!!
  ################################################################
  if (.Platform$OS.type=="windows")
  {
    no_cores = 1
  } else
  {
    no_cores = detectCores(logical = FALSE, all.tests = T)
  }
  
  # Find centres and footprints
  fp = mclapply(ids, function(x){process_peak(x, reads, peaks, min.reads = min.reads, rel.threshold = rel.threshold, abs.threshold = abs.threshold, footprint.frac = footprint.fraction )}, mc.cores=no_cores, mc.preschedule = TRUE)
  fp = do.call("rbind", fp)
  fp = as.data.table(fp[complete.cases(fp),])
  
  reads = mclapply(ids, function(x){reads[id==x, cent.pos := pos-peaks[id==x,centre]]}, mc.cores=no_cores, mc.preschedule = FALSE)
}



###########################
# Save Data structures
###########################
save(fp, file="Data/footprints_data.Rba")
save(peaks, file="Data/peak_data.Rba")
save(reads, file="Data/reads_data.Rba")
