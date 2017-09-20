#########################
# Load
#########################
library(parallel)
library(data.table)
source("ATAC_functions.R")
load("Data/reads_data.Rba")


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
dist_mem = FALSE

# number of cores to use
no_cores = 0
if (no_cores==0)
{
  no_cores = detectCores(logical = FALSE, all.tests = T)
}


################
# Run
################
ids = unique(reads$id)

if(dist_mem)
{
  #####################################################
  #Distributed memory !!!Use for multicore Windows!!!
  #####################################################
  cl = makeCluster(no_cores)
  
  # Find the footprints and centres
  clusterExport(cl=cl, list("process_peak", "find_peaks", "find_footprints", "find_centre", "reads", "min.reads", "min.peak.fraction", "footprint.fraction"))
  fp.dist = parLapply(cl, ids, function(x){process_peak(id = x, reads = reads[reads$id==x,], min.reads = min.reads, rel.threshold = rel.threshold, abs.threshold = abs.threshold, footprint.frac = footprint.fraction )})
  fp = do.call("rbind", fp.dist)
  
  # Create a data frame of peaks with the centre of the peak
  clusterExport(cl=cl, list("fp"))
  centres = parLapply(cl, ids, function(x){data.frame(x,fp[fp$id==x,][1,]$centre, fp[fp$id==x,][1,]$nfootprints, fp[fp$id==x,][1,]$nreads)})
  centres = do.call("rbind", centres)
  colnames(centres) =  c("id", "centre", "no.footprints", "no.reads")
  
  # Add the centre paramater to the reads
  clusterExport(cl=cl, list("centres"))
  centred_reads = parLapply(cl, ids, function(x){reads[reads$id==x,]$pos-centres[centres$id==x,]$centre})
  reads$centred.pos=unlist(centred_reads)
  
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
  fp = mclapply(ids, function(x){process_peak(id = x, reads = reads[reads$id==x,], min.reads = min.reads, rel.threshold = rel.threshold, abs.threshold = abs.threshold, footprint.frac = footprint.fraction )}, mc.cores=no_cores, mc.preschedule = TRUE)
  fp = do.call("rbind", fp)
  
  # Centre of peaks
  centres = mclapply(ids, function(x){data.frame(x,fp[fp$id==x,][1,2], fp[fp$id==x,][1,3], fp[fp$id==x,][1,5])}, mc.cores = no_cores, mc.preschedule = TRUE)
  centres = do.call("rbind", centres)
  colnames(centres) =  c("id", "centre", "no.footprints", "no.reads")
  
  # Centre the reads
  centred_reads = mclapply(ids, function(x){reads[reads$id==x,]$pos-centres[centres$id==x,]$centre},mc.cores = no_cores, mc.preschedule = TRUE)
  reads$centred.pos=unlist(centred_reads)  
}



###########################
# Save Data structures
###########################
save(fp, file="Data/footprints_data.Rba")
save(centres, file="Data/peak_centre_data.Rba")
save(reads, file="Data/centred_reads_data.Rba")
