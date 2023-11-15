getPeaks <- function(frame,chans,tinypeak.removal=tinypeak.removal){
  
  ##=====================================================================================================================
  ## Finds the peaks in the given density
  ## Args:
  ##   dens: density of the channel whose peaks are going to be found. It is a variable of class 'density'.
  ##   w: the length of the window where the function searches for the peaks. If all peaks required, use the default w=1.
  ## Value:
  ##   peaks in the density of the provided channel
  ## Written by Mehrnoush Malek
  # Revised on October 2014
  # Modified by Albina Rahim July 2019
  ##---------------------------------------------------------------------------------------------------------------------
  data <- exprs(frame)[,chans]
  dens <- density(data[which(!is.na(data))])
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  d <- dens$y
  w <- 1
  peaks <- c()
  peaks.ind <- c()
  for(i in 1:(length(d)-w)){
    if(d[i+w] > d[(i+w+1):(i+2*w)] && d[i+w] > d[i:(i+w-1)] && d[i+w] > tinypeak.removal*max(d)){ # also removes tiny artificial peaks less than ~%4 of the max peak
      peaks <- c(peaks, dens$x[i+w])
      peaks.ind <- c(peaks.ind, i+w)
    }
  }
  return(list(Peaks=peaks, Dens=dens,Ind=peaks.ind))
}