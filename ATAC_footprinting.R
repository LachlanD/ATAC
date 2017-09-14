#################################
#Functions
#################################


# #############################
# #Old footprints
# #############################
# localmins <- function(y)
# {
#   tfs  = c()
#   dy=0
#   maxy = y[1]
#   miny = y[1]
#   for (i in 2:length(y))
#   {
#     maxy = max(c(maxy, y[i]))
#     miny = min(c(miny, y[i]))
#     ndy = y[i] - y[i-1]
#     if (dy<0 && ndy>=0 && maxy>1.1*miny)
#     {
#       tfs = c(tfs, i)
#       maxy=y[i]
#       miny=y[i]
#     }
#     dy = ndy
#   }
#   return(tfs)
# }


#################################
#Current peak finding
#################################
find_peaks <- function(y, peak.threshold, footprint.fraction)
{
  local.max = max(y)
  if (local.max<peak.threshold)
  {
    return(c())
  }
  max.index = which.max(y)
  peaks = c(max.index)
  
  w = length(y)
  
  
  # Descend forward until reach footprint threshold
  i=max.index
  while(i < w  && y[i]>local.max*footprint.fraction)
  {
    i = i + 1
  }
  # start Tracking local min
  local.min = y[i]
  while(i < w)
  {
    i = i + 1
    local.min = min(c(local.min,y[i]))
    
    # Continue until back above the footprint threshold
    if(y[i]>local.min/footprint.fraction)
    {
      # recursively call on the remaining curve
      peaks = c(peaks, i+find_peaks(y[i:w], peak.threshold, footprint.fraction))
      break
    }
  }
  
  
  # Descend backward until reach footprint threshold
  i=max.index
  while(i > 1  && y[i]>local.max*footprint.fraction)
  {
    i = i - 1
  }
  # start tracking local min
  local.min = y[i]
  while(i > 1)
  {
    i = i - 1
    local.min = min(c(local.min, y[i]))
    
    # Continue auntil back above the footprint threshold
    if(y[i]>local.min/footprint.fraction)
    {
      # recursively call on the remaining curve
      peaks = c(peaks, find_peaks(y[0:i],peak.threshold, footprint.fraction))
      break
    }
  }
  return(peaks)
}


#############################
#Current footprint finding
#############################
find_footprints <- function(y, peaks)
{
  npeaks = length(peaks)
  footprints = rep(0,npeaks-1)
  if(npeaks>1)
  {
    for(i in 1:(npeaks-1))
    {
      footprints[i] =  peaks[i]+which.min(y[peaks[i]:peaks[i+1]])
    }
  }
  return(footprints)
}

###########################
#Current centres
###########################
find_centre <- function(peaks, footprints)
{
  # Centre on centre peak if the number of peaks is odd 
  if(length(footprints)%%2==0)
  {
    return(peaks[length(footprints)/2+1])
  }
  # Centre on centre footprint if the number of footprints is odd
  return(footprints[(length(footprints)+1)/2])
}


# ###########################
# #Old peaks
# ###########################
# peaks <- function(y, frac_of_max)
# {
#   miny = y[1]
#   maxy = y[1]
#   
#   m = max(y)
#   
#   p = c()
#   
#   dy =  0
#   for (i in 2:length(y))
#   {
#     maxy = max(c(maxy, y[i]))
#     miny = min(c(miny, y[i]))
#     ndy = y[i] - y[i-1]
#     
#     if(y[i] > frac_of_max*m && dy>0 && ndy<=0)  
#     {
#       maxy = y[i]
#       if(maxy >= 1.1*miny)
#       {
#         p = c(p, i)
#         miny = y[i]
#       }
#       else if(maxy > p[length(p)])
#       {
#         p[length(p)] = i
#         miny = y[i]
#       }
#     }
#     dy = ndy
#   }
#   return(p)
# }

####################
# Wrap functions
####################
process_peak <- function(id, reads, min.reads, min.peak.frac, footprint.frac)
{
  r = reads[reads$id==id,]
  centre = 0
  nfootprints=NA
  footprint.pos = NA
  nreads = nrow(r)
  
  d = density(r$pos)
  if(nreads >= min.reads)
  {
    peaks = sort(find_peaks(d$y, max(d$y)*min.peak.frac, footprint.frac))
    footprints = find_footprints(d$y, peaks)
    nfootprints = length(footprints)
    if (nfootprints>0)
    {
      footprint.pos = round(d$x[footprints])
    }
    centre = round(d$x[find_centre(peaks, footprints)])
  }
  
  return(data.frame(id, centre, nfootprints, footprint.pos, nreads))
}


# ########################
# #Old find centre
# ########################
# findcentres <- function(id, reads, peaks, localmins)
# {
#   r = reads[reads$id==id,]
#   centre = 0
#   ntfs=0
#   positions=0
#   nreads=nrow(r)
#   if(nreads > 100)
#   {
#     d = density(r$pos)
#     pos = density(r[r$strand=='+',]$pos)
#     neg = density(r[r$strand=='-',]$pos)
#     
#     pos_peaks = peaks(pos$y)
#     neg_peaks = peaks(rev(neg$y))
#     
#     left = pos$x[pos_peaks[1]]
#     right = neg$x[length(neg$x)-neg_peaks[1]]
#     
#     if(left<right)
#     {
#       c = d$x>left&d$x<right
#       tfs = localmins(d$y[c])
#       ntfs = length(tfs)
#       
# 
#       if(ntfs>0)
#       {
#         positions = d$x[c][tfs]
#         centre = mean(c(d$x[c][tfs[ceiling(ntfs/2)]], d$x[c][tfs[floor(ntfs/2+1)]]))
#       }
#     }
#   }
#   return(data.frame(id, centre, ntfs, positions, nreads))
# }



###################
# Parameters
###################
# Use a small subset of peaks
small.set = TRUE
first.n.peaks = 1000

# Minimum read coverage for a peak to be considered  
min.reads = 200

# peaks must be with this fraction of the absolute maximum
min.peak.fraction = 0.5

# Footprint must fall this far below the smaller of the two peaks
footprint.fraction = 0.9

#############################
#Setup
#############################
if(.Platform$OS.type=="windows")
{
  #setwd("D:/") # Windows exHDD
  setwd("C:/FastFiles/ATAC") # Windows Fast ssd
} else
{
  setwd("/media/lachlan/ldryburghEx") # Linux exHDD
}
reads = read.delim('MODEL_120_peaks.reads', header = TRUE, sep = '\t')
ids = unique(reads$id)

#####################################
#Create a limited dataset if wanted
#####################################
if(small.set)
{
  ids = ids[1:first.n.peaks]
  reads=reads[reads$id %in% ids,]
}

################################################################
# Undistribute memory Versions  !!! Windows can only use 1 !!!
################################################################
library(parallel)

if (.Platform$OS.type=="windows")
{
  no_cores = 1
} else
{
  no_cores = detectCores(logical = FALSE, all.tests = T)
}
  
# Find centres and footprints
fp = mclapply(ids, function(x){process_peak(id = x, reads = reads[reads$id==x,], min.reads = min.reads, min.peak.frac = min.peak.fraction, footprint.frac = footprint.fraction )}, mc.cores=no_cores, mc.preschedule = TRUE)
fp = do.call("rbind", fp)

# Centre of peaks
centres = mclapply(ids, function(x){data.frame(x,fp[fp$id==x,][1,2], fp[fp$id==x,][1,3], fp[fp$id==x,][1,5])}, mc.cores = no_cores, mc.preschedule = TRUE)
centres = do.call("rbind", centres)
colnames(centres) =  c("id", "centre", "no.footprints", "no.reads")

# Centre the reads
centred_reads = mclapply(ids, function(x){reads[reads$id==x,]$pos-centres[centres$id==x,]$centre},mc.cores = no_cores, mc.preschedule = TRUE)
reads$centred.pos=unlist(centred_reads)



#####################################################
#Distributed memory !!!Use for multicore Windows!!!
#####################################################
library(parallel)
no_cores = detectCores(logical = FALSE, all.tests = T)


cl = makeCluster(no_cores)

# Find the footprints and centres
clusterExport(cl=cl, list("process_peak", "find_peaks", "find_footprints", "find_centre", "reads", "min.reads", "min.peak.fraction", "footprint.fraction"))
fp.dist = parLapply(cl, ids, function(x){process_peak(id = x, reads = reads[reads$id==x,], min.reads = min.reads, min.peak.frac = min.peak.fraction, footprint.frac = footprint.fraction )})
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


###########################
# Save Data structures
###########################
save(fp, file="footprints_data.Rba")
save(centres, file="peak_centre_data.Rba")
save(reads, file="centred_reads_data.Rba")


#############################
# Find absolute positions
#############################
summits = read.delim2("MODEL_120_peaks.bed", header = FALSE)
colnames(summits)=c("chr", "start", "end", "abs.summit", "id")
load("footprints_data.Rba")


footprints_into_bed <- function(pk, fp)
{
  pos = pk$abs.summit+fp[fp$id==pk$id,]$footprint.pos
  
  return(data.frame(pk$chr,pk$start, pk$end, pos, pk$id))
}

fp = fp[complete.cases(fp),]
ids= unique(fp$id)

fp.bed = mclapply(ids, function(x){footprints_into_bed(summits[summits$id==x,], fp)}, mc.cores = no_cores, mc.preschedule = FALSE)
fp.bed = do.call("rbind",fp.bed)
colnames(fp.bed)=c("chr", "start", "end", "footprint", "id")

write.table(fp.bed, file="MODEL_120_footprints.bed", sep= "\t", quote= FALSE, row.names = FALSE, col.names = FALSE)

###########################################
#Visualise loop for debugging and tuning
##########################################
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
  
  peaks = sort(find_peaks(d$y, max(d$y)*min.peak.fraction, footprint.fraction))
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

# #############################
# #Initial looped version
# #############################
# reads$offset=0
# reads$ntf=0
# for (i in 1:length(ids))
# {
#   r = reads[reads$id==ids[i],]
#   if(nrow(r) > 100)
#   {
#     d = density(r$pos)
#     pos = density(r[r$strand=='+',]$pos)
#     neg = density(r[r$strand=='-',]$pos)
#     
#     #plot(pos,col='blue', ylim=c(0,max(pos$y,pos$y,neg$y)*1.1), xlim=c(-200,200))
#     
#     
#     #plot(d,col='blue', ylim=c(0,max(pos$y,pos$y,neg$y)*1.1), xlim=c(-200,200))
#     #lines(pos,col='green')
#     #lines(neg,col='red')
#     
#     #dy = diff(d$y)
#     #p = which.min(abs(dy))
#     #abline(v= d$x[p])
#     
#     #dy = dy[c(-p, -p+1, -p-1)]
#     
#     #p = which.min(abs(dy))
#     #abline(v= d$x[p])
#     
#     #dy = dy[c(-p, -p+1, -p-1)]
#     
#     #p = which.min(abs(dy))
#     #abline(v= d$x[p])
#     
#     #mp = which.max(pos$y)
#     #left = pos$x[mp]
#     #abline(v=left, col='green')
#     #mn = which.max(neg$y)
#     #right = neg$x[mn]
#     #abline(v=right, col='red')
#     
#     pos_peaks = peaks(pos$y)
#     neg_peaks = peaks(rev(neg$y))
# 
#     left = pos$x[pos_peaks[1]]
#     right = neg$x[length(neg$x)-neg_peaks[1]]
#         
#     if(left<right)
#     {
#       c = d$x>left&d$x<right
#       
#       tfs = localmins(d$y[c])
#       
#       ntfs = length(tfs)
#       
#       #for(j in 1:ntfs)
#       #{
#       #  abline(v=d$x[c][tfs[j]])
#       #}
#       
#       reads[reads$id==ids[i],]$ntf = ntfs
#       if(ntfs>0)
#       {
#         centre = mean(c(d$x[c][tfs[ceiling(ntfs/2)]], d$x[c][tfs[floor(ntfs/2+1)]]))
#         reads[reads$id==ids[i],]$offset = centre
#         #abline(v=centre, col='red')
#       }
#         
#       
#       #centre = which.min(d$y[c])
#       #dy = diff(d$y[c])
#       #if(sum(c) > 0)
#       #{
#       #  if(length(dy)>=centre)
#       #  {
#       #    if(abs(dy[centre])<0.00001)
#       #    {
#       #      #abline(v=d$x[c][centre], col='blue')
#       #      print(d$x[c][centre])
#       #      reads[reads$id==ids[i],]$offset=d$x[c][centre]
#       #      reads[reads$id==ids[i],]$ntf=1
#       #    }
#       #  }
#       #}
#     }
#   }
#   #invisible(readline(prompt="Press [enter] to continue"))
# }

