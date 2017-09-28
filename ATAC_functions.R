#################################
#Functions
#################################

#################################
#Recursive peak finding, sorted
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
    
    # Continue until back above the footprint threshold
    if(y[i]>local.min/footprint.fraction)
    {
      # recursively call on the remaining curve
      peaks = c(find_peaks(y[0:i],peak.threshold, footprint.fraction), peaks)
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




####################
# Wrap functions
####################
process_peak <- function(i, r, peaks, min.reads, rel.threshold, abs.threshold, footprint.frac)
{
  no.reads = nrow(r)
  
  peaks[id==i, nreads := no.reads]
  
  if(no.reads >= min.reads)
  {
    d = density(r[id==i,pos])
    threshold = max(c(abs.threshold, max(d$y)*rel.threshold))
    p = find_peaks(d$y, threshold, footprint.frac)
    if(length(p)>0)
    {
      footprints = find_footprints(d$y, p)
      no.footprints = length(footprints)
      if (no.footprints>0)
      {
        footprint.pos = round(d$x[footprints])
        cent = round(d$x[find_centre(p, footprints)])
        peaks[id==i, centre := cent]
        peaks[id==i, nfootprints := no.footprints]
        return(data.frame(i, footprint.pos, round(d$x[p[1:length(p)-1]]),round(d$x[p[2:length(p)]])))
      } else
      {
        peaks[id==i, centre := round(d$x[p[1]])]
        peaks[id==i, nfootprints := 0]
      }
    } 
  }else
  {
    peaks[id==i, centre := NA]
    peaks[id==i, nfootprints := NA]
  }
  
  return(data.frame())
}
