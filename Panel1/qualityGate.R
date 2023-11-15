
# Functions written by Albina Rahim

qualityGate <- function(f, scat.chans, channels.ind, centre, panel.no){
  # Applies quality gates on the flow frame f. Returns the live singlets.
  #
  # Args:
  #   f: flow frame
  #   scat.chans: named vector with scatter channel indices
  #   channels.ind: named vector with the other channel indices, channels.ind['Live'] is required
  #   panel.no: integer, 1 for panel 1 and 2 for panel 2
  #
  # Returns: A list containing
  #       scat.chans: updated scatter channel indices (some may be switched)
  #       live: a flowDensity object containing live cells
  #       FSCsinglets: a flowDensity object containing the live/FSC singlets
  #       singlets: a flowDensity object containing live/FSC singlets/SSC singlets 
  
  # Gating Live ---------------------------------------------------------------------------------------
  
  temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(NA, F), gates = c(NA, 1.5))
  
  if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
    #fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5)
    fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F)
  }else{
    fsca.live.gate <- c(deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), all.cuts = T, tinypeak.removal = 0.02, use.upper = T, upper = F, percentile = NA),
                        deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5))
    
    # Find FSC-A peak and pick gate that is just below it 
    maxDens <- density(getflowFrame(temp.flowD)@exprs[,c(scat.chans['FSC-A'])])
    fscapeak.lcn <- findpeaks(maxDens$y)
    fscapeak.lcn.minthreshold <- 35000
    # The next bit is to correct the fscapeak.lcn.minthreshold for some BCM files
    maj.peaks <- fscapeak.lcn[which(fscapeak.lcn[,1] > 0.2*max(fscapeak.lcn[,1])),]
    if(length(maj.peaks) > 4){
      if((maxDens$x[maj.peaks[2,2]] < 35000) & (maxDens$x[maj.peaks[2,2]] > 14000)){
        fscapeak.lcn.minthreshold <- maxDens$x[maj.peaks[2,2]] - 2000
      }
    }
    fscapeak.lcn <- fscapeak.lcn[which(maxDens$x[fscapeak.lcn[,2]] > fscapeak.lcn.minthreshold),]
    if(length(fscapeak.lcn) > 4){
      fscapeak.lcn <- maxDens$x[fscapeak.lcn[which(fscapeak.lcn[,1] > 0.2*max(fscapeak.lcn[,1])), 2]]
      fscapeak.lcn <- fscapeak.lcn[which.min(abs(fscapeak.lcn - 105000))]
    }else{
      fscapeak.lcn <- maxDens$x[fscapeak.lcn[2]]
    }
    
    # subtraction of 7500 is needed for files such as TCP PANEL_B_ACOI_65_F7F07032.fcs
    fsca.live.gate.temp <- max(fsca.live.gate[which(fsca.live.gate < (fscapeak.lcn - 7500))])
    ## Adding this part for UCD Panel 1
    if(fsca.live.gate.temp-50000 > 10000){
      if(min(fsca.live.gate[which(fsca.live.gate < (fscapeak.lcn - 7500))]) < 50000){
        fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']))
        if(fsca.live.gate > 60000){
          fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), tinypeak.removal = 0.01)
          if(fsca.live.gate < 50000){
            fsca.live.gate <- 55000
          }
        }else if(fsca.live.gate <= 55000){
          fsca.live.gate<- mean(c(deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), all.cuts = T, tinypeak.removal = 0.02, upper = F, percentile = NA),
                                  deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5)))
        }
      }else{
        fsca.live.gate <- min(fsca.live.gate[which(fsca.live.gate < (fscapeak.lcn - 7500))])
        if(fsca.live.gate <= 55000){
          fsca.live.gate <- mean(c(deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), all.cuts = T, tinypeak.removal = 0.02, upper = F, percentile = NA),
                                   deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5)))
          if(fsca.live.gate > 70000){## changing this from >60000 to >70000
            fsca.live.gate <- 50000 ## Changing this from 55000 to 50000
          }
        }else if(fsca.live.gate > 60000){
          fsca.live.gate <- 50000 ## Changing this from 55000 to 50000
        }
      }
    }else{
      fsca.live.gate <- c(deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), all.cuts = T, tinypeak.removal = 0.02, upper = F, percentile = NA),
                          deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5))
      
      if(length(fsca.live.gate) > 1){
        fsca.live.gate <- c(min(fsca.live.gate))
      }
    }
    
  }
  
    
  # Find the peak corresponding to the dead cells
  maxDens <- density(f@exprs[,c(channels.ind['Live'])])
  maxDens2 <- density(getflowFrame(temp.flowD)@exprs[,c(channels.ind['Live'])])
  peak.lcns <- findpeaks(maxDens$y)
  n <- nrow(peak.lcns)
  
  if(n > 1){ # if the dead peak is detectable
    if(panel.no == 2){
      if(channels.ind["Live"] == channels.ind['CD19']){ # expect 3 peaks, the dead peak should be the one with the highest Live/Dead value
        largest.peaks <- sort(peak.lcns[,1], index.return = TRUE)
        largest.peaks <- peak.lcns[largest.peaks$ix[(n-2):n], 2]
        dead.peak <- maxDens$x[max(largest.peaks)]
      }else{
        live.peak <- maxDens2$x[which.max(maxDens2$y)]
        mindeadpeak.lcn <- live.peak + 0.6*(quantile(f@exprs[,c(channels.ind['Live'])], c(0.999)) -live.peak)
        peak.lcns <- peak.lcns[which(maxDens$x[peak.lcns[,2]] > mindeadpeak.lcn), ]
        if(length(peak.lcns) > 4){
          dead.peak <- maxDens$x[peak.lcns[which.max(peak.lcns[,1]), 2]]
        }else if(length(peak.lcns) == 4){
          dead.peak <- maxDens$x[peak.lcns[2]]
        }else{
          dead.peak <- 2.4
        }
      }
    }else{ # Panel 1
      live.peak <- maxDens2$x[which.max(maxDens2$y)]
      peak.lcns <- peak.lcns[which(maxDens$x[peak.lcns[,2]] > (live.peak + 1.1)), ]
      if(length(peak.lcns) > 4){
        dead.peak <- maxDens$x[peak.lcns[which.max(peak.lcns[,1]), 2]]
      }else if(length(peak.lcns) == 4){
        dead.peak <- maxDens$x[peak.lcns[2]]
      }else{
        dead.peak <- 2.4
      }
      
    }
  }else{ # if cannot detect dead peak
    dead.peak <- 2.4
  }
  
  
  is.DAPI <- ((length(grep('DAPI', pData(parameters(f))$desc)) > 0) | (length(grep('DAPI', pData(parameters(f))$name)) > 0))
  
  if(!is.DAPI){
    
    if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
      live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.upper = T, upper = T)
      if(live.gate <= 0.5 | round(live.gate,1) > 2.5){
        live.gate <- deGate(f, channel = c(channels.ind["Live"]), tinypeak.removal = 0.01)

      }
    }else{
      live.gate <- deGate(f, channel = c(channels.ind["Live"]), all.cuts = T, tinypeak.removal = 0.001)
      live.gate <- live.gate[which(live.gate > maxDens2$x[which.max(maxDens2$y)])]
      
      if(length(live.gate) > 0){ # if there are any gates above the live peak
        live.gate <- live.gate[which(live.gate < dead.peak)]
        if(length(live.gate) > 0){
          live.gate <- live.gate[which.min(abs(dead.peak - 0.5 - live.gate))]
          if(live.gate > 3){
            live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.percentile = T, percentile = 0.9)
            if(live.gate < 2.5 | live.gate > 3){
              live.gate <- 2.75
            }
          }else if(live.gate < 2.25){
            live.gate <- min(2.25, deGate(f, channel = c(channels.ind["Live"]), use.percentile = T, percentile = 0.9))
            if(live.gate < 2.25){
              live.gate <- 2.5
              
            }
          }
          
          # else if(live.gate < 1){
          #   live.gate <- 2
          # }
        }else{
          live.gate <- dead.peak - 0.5
          if(live.gate > 3){
            live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.percentile = T, percentile = 0.9)
            if(live.gate < 2.5  | live.gate > 3){
              live.gate <- 2.75
            }
          }
          # else if(live.gate < 1){
          #   live.gate <- 2
          # }
        }
      }else{ # For CIPHE 16-Oct-04 there are no dead cells, so there are no gates found above the live peak
        live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.upper = T, upper = T, alpha = 0.005, tinypeak.removal = 0.2)
        if(live.gate > 3){
          live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.percentile = T, percentile = 0.9)
          if(live.gate < 2.5  | live.gate > 3){
            live.gate <- 2.5
          }
        }
      }
      
    }
       
   
  }else{
    
    if(centre == 'BCM'){
      live.gate <- max(dead.peak - 0.5, 2.75)
    }else{
      temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, NA), gates = c(fsca.live.gate + 10000, NA))
      maxDens <- density(getflowFrame(temp.flowD)@exprs[,c(channels.ind['Live'])])
      start.idx <- which.max(maxDens$y)
      stop.idx <- which.min(abs(maxDens$x - dead.peak + 0.2))
      valley.yval <- min(maxDens$y[start.idx:stop.idx])
      live.gate <- min(maxDens$x[which(maxDens$y[start.idx:stop.idx] < 2*valley.yval) + (start.idx -1)]) 
    }
  }
  
  live.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, F), gates = c(fsca.live.gate, live.gate))
  
  
  if(live.flowD@proportion < 80){
    live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.percentile = T, percentile = 0.97)
    if(live.gate > 2.85){
      live.gate <- 2.85
      fsca.live.gate <- mean(c(deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), tinypeak.removal = 0.01),
                                    deGate(temp.flowD@flow.frame, channel = c(scat.chans['FSC-A']), use.upper = T, upper = F, tinypeak.removal = 0.5)))
      if(fsca.live.gate > 70000){
        live.gate <- 2.6
        fsca.live.gate <- 50000 ## Changing this from 55000 to 50000
      }
    }
      
    live.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, F), gates = c(fsca.live.gate, live.gate))
  }
  
  # if(live.gate < 1 & live.flowD@proportion < 80){
  #   live.gate <- live.gate[which.max(abs(dead.peak - 0.5 - live.gate))]
  #   live.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, F), gates = c(fsca.live.gate, live.gate))
  # }
  
  plotDens(f, c(scat.chans['FSC-A'], channels.ind['Live']), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2)
  abline(v=fsca.live.gate)
  lines(live.flowD@filter)
  
  # Gating size ----------------------------------------------------------------------------------------
  ## Temporarily commeting this line to gate the size.flowD
  #size.flowD <- live.flowD
  ## Temporary code
  ssca.gate <- deGate(live.flowD@flow.frame, channel = c(scat.chans["SSC-A"]), use.percentile = T, percentile = 0.99)
  
  # plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'],scat.chans['SSC-A']), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(exprs(live[singlets.flowD@index, c(scat.chans['FSC-A'], scat.chans['SSC-A'])]), col=2, pch=".")
  # abline(h=ssca.gate)
  
  size.flowD <- flowDensity(live.flowD@flow.frame, channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), position = c(NA, F), gates = c(NA, ssca.gate))
  # Gating FSC singlets --------------------------------------------------------------------------------
  
  # Check if FSC-H and FSC-W channel names are switched
  #temp <- getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
  temp <- size.flowD@flow.frame@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
  try({
    FSC.angle <- atan(summary(lm(temp[, c(scat.chans['FSC-H'])] ~ temp[, c(scat.chans['FSC-A'])]))$coefficients[2,1])
    if(FSC.angle < pi/9){ # assume that the FSC-H and FSC-W channels are interchanged
      scat.chans['FSC-H'] <- grep('FS-W|FSC-W', colnames(f))
      scat.chans['FSC-W'] <- grep('FS-H|FSC-H', colnames(f))
      temp <- getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
    }
  })
  
  if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){## This condition is for GMC dataset   
    theta0 <- -0.15
  }else{
    # cannot use the FSC.angle to rotate the data directly b/c it will be rotated towards doublets
    # So I calculate the angle to rotate the data as follows:
    # Pick out the FSC-A versus FSC-H distribution around FSC-A of x0, and find the peak in the FSC-H distribution
    # And set the angle to be equal to arctan(peak in FSC-H / x0)
    
    # This line is for the CIPHE data to get rid of a pile-up close to the FSC-H axis, shouldn't affect other centers much
    temp.flowD <- flowDensity(size.flowD@flow.frame, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA,F), gates = c(NA, 0.98*max(getflowFrame(size.flowD)@exprs[,c(scat.chans['FSC-H'])])))
    #x0 <- 125000
    x0 <- quantile(getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-H'])], c(0.05, 0.99))
    x0 <- x0[1] + 0.75*(x0[2]-x0[1])
    temp <- temp[which(abs(temp[ , scat.chans['FSC-A']] - x0) < 5000),]
    maxDens <- density(temp[, 2])
    peak.lcn <- findpeaks(maxDens$y) # get all peaks
    peak.lcn <- peak.lcn[which(peak.lcn[,1] > 0.1*max(peak.lcn[,1])),] # choose only peaks with significant y values
    if(length(peak.lcn) > 4){
      peak.lcn <- maxDens$x[peak.lcn[,2]] # get x-values
      peak.lcn <- peak.lcn[which.min(abs(peak.lcn - x0))] # assume we are interested in peak closed to diagonal
    }else{
      peak.lcn <- maxDens$x[peak.lcn[2]]
    }
    theta0 <- -atan(peak.lcn/x0)
    
  }
  
   
  # rotate data
  temp.flowD <- flowDensity(size.flowD@flow.frame, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA,F), gates = c(NA, 0.975*max(getflowFrame(size.flowD)@exprs[,c(scat.chans['FSC-H'])])))
  rot <- rotate.data(getflowFrame(temp.flowD), c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = -theta0)$data
  
  # plotDens(rot, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
  # abline(h=gate)
  
  maxDens <- density(rot@exprs[,c(scat.chans['FSC-H'])])
  maxDens <- smooth.spline(maxDens$x, maxDens$y, spar = 0.25)
  
  gate <- deGate(rot, channel = c(scat.chans['FSC-H']), all.cuts = T, tinypeak.removal = 0.001)
  #gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.percentile = T, percentile = 0.025)
  gate0 <- gate
  singlets.peak <- maxDens$x[which.max(maxDens$y)]
  if(is.null(names(gate))){ 
    idx <- which(gate < (singlets.peak - 10000))
    if(length(idx) > 0){
      gate <- gate[idx]
      p05percentile <- min(deGate(rot, channel = c(scat.chans['FSC-H']), use.percentile = T, percentile = 0.0005), 50000) 
      gate <- gate[which.min(abs(singlets.peak - 0.45*(singlets.peak - p05percentile) - gate))]
    }else{
      gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.01)
    }
  }
  
  if((singlets.peak - gate) > 20000){ ##  (which was previously changed from 50000 to 20000)
    gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.5) - 12000 ## Changed this to 12000 from 15000
  }
  
  
  rot <- rotate.data(temp.flowD@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = -theta0)$data
  temp <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, gate))
  #   , ellip.gate = T, scale = .9999999);
  temp@filter <- rotate.data(temp@filter, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
  temp@flow.frame <- rotate.data(temp@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  #temp@flow.frame <- rotate.data(getflowFrame(temp), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  # Gate really closesly - need this to tell if SSC-H/SSC-W channels are interchanged 
  temp2 <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, singlets.peak - 2000))
  temp2@flow.frame <- rotate.data(temp2@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  #temp2@flow.frame <- rotate.data(getflowFrame(temp2), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  FSCsinglets.flowD <- temp
  
  # plotDens(live.flowD, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=singlets.peak-2000)
  # lines(FSCsinglets.flowD@filter)
  
  
  # Gating SSC singlets -------------------------------------------------------------------------
  
  #if(!is.na(scat.chans['SSC-W'])){ # skip second singlet gate if no SSC-W   ## This if condition not needed with the files we are running on right now
  if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){  
    temp <- temp2@flow.frame@exprs[, c(scat.chans['SSC-W'], scat.chans['FSC-H'])]
  }else{
    temp <- temp2@flow.frame@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])]
    
  }
  
 
  
  if(centre == 'BCM'){
    
    # Check if SSC-H and FSC-W channel names are switched
    # For BCM In approx. 2/3 of the FCS files the FSC-H and FSC-W channel names are interchanged. 
    # This occurs randomly so we have a check in the code when gating the FSC singlets to see if the FSC-A vs FSC-H distribution looks as expected 
    SSC.angle <- atan(summary(lm(temp[, c('SSC-H')] ~ temp[, c('SSC-W')]))$coefficients[2,1])
    if(SSC.angle < pi/16){ 
      scat.chans['SSC-H'] <- grep('SSC-W', colnames(f))
      scat.chans['SSC-W'] <- grep('SSC-H', colnames(f))
    }
    
    
    # Better than below for BCM. Doesn't work that well on CIPHE. Haven't tested on other centres
    ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.001, all.cuts = T)
    maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[,c(scat.chans['SSC-W'])])
    temp2.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(NA, T), gates = c(NA, (maxDens$x[which.max(maxDens$y)] + 20000)))
    ss.high <- c(ss.high, deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.001, all.cuts = T))
    
    ss.high <- ss.high[which(ss.high > (maxDens$x[which.max(maxDens$y)] + 5000))]
    if(length(ss.high) > 1){
      ss.high <- min(ss.high)
    }else{
      ss.high <- min(deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T,tinypeak.removal = 0.4),
                     deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T, tinypeak.removal = 0.4))
    }
    if((ss.high - maxDens$x[which.max(maxDens$y)]) > 30000){
      ss.high <- min(deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T, alpha = 0.9, tinypeak.removal = 0.9),
                     deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T, alpha = 0.9, tinypeak.removal = 0.9))
    }
    
  }else{
    
    #This works for CIPHE/TCP/JAX/GMC
    if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
      ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.4, all.cuts = T, use.percentile = T, percentile = 0.95)
    }else{
      ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.4, all.cuts = T, use.percentile = T, percentile = 0.95)
    }
       
    #ss.high2 <- NA
    if(!is.null(names(ss.high)) & !(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3)){ # if there is a 95% gate
      # if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
      #   ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.9, use.upper = T, upper = T, alpha = 0.9)
      #   maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[,c(scat.chans['SSC-A'])])
      #   
      # }else{
      #   ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.9, use.upper = T, upper = T, alpha = 0.9)
      #   maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[,c(scat.chans['SSC-W'])])
      #   
      # }
      
      ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.9, use.upper = T, upper = T, alpha = 0.9)
      maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[,c(scat.chans['SSC-W'])])
      
      if((ss.high - maxDens$x[which.max(maxDens$y)]) < 50000){
        # if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        #   temp.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-A'], scat.chans['FSC-H']), position=c(T, NA), gates = c(ss.high - 12000, NA))
        #   ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.02)
        #   if(!is.null(names(ss.high))){
        #     ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-A']), tinypeak.removal = 0.02, use.upper = T, upper = T) + 3000
        #   }
        #   maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['FSC-H'])])
        #   temp2.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-A'], scat.chans['FSC-H']), position=c(NA, T), gates = c(NA, (maxDens$x[which.max(maxDens$y)] + 20000)))
        #   ss.high2 <- deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-A']), upper = T, tinypeak.removal = 0.4, percentile = NA)
        #   
        # }else{
        #   temp.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(T, NA), gates = c(ss.high - 12000, NA))
        #   ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02)
        #   if(!is.null(names(ss.high))){
        #     ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02, use.upper = T, upper = T) + 3000
        #   }
        #   maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-H'])])
        #   temp2.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(NA, T), gates = c(NA, (maxDens$x[which.max(maxDens$y)] + 20000)))
        #   ss.high2 <- deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), upper = T, tinypeak.removal = 0.4, percentile = NA)
        #   
        # }
        
        temp.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(T, NA), gates = c(ss.high - 12000, NA))
        ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02)
        if(!is.null(names(ss.high))){
            ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02, use.upper = T, upper = T) + 3000
        }
        maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-H'])])
        temp2.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(NA, T), gates = c(NA, (maxDens$x[which.max(maxDens$y)] + 20000)))
        if(temp2.flowD@proportion < 0.1){
          if(ss.high > 80000){## Changing this from > 75000 to > 80000
            ss.high <- 75000
          }else{
            ss.high <- ss.high
          }
          
        }else{
          ss.high2 <- deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T, tinypeak.removal = 0.4, percentile = NA)
          
          ss.high <- c(ss.high, ss.high2)
          #maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'])])
          ss.high <- min(ss.high)
        }
       
        
       
        # ## This part is only for GMC where there is only SSC-A is present.
        # if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        #     if((ss.high) > 50000 | ss.high < 49000){
        #       #ss.high <- 50000
        #       if(maxDens$x[which.max(maxDens$y)] < 60000){
        #         ss.high <- maxDens$x[which.max(maxDens$y)]
        #       }else{
        #         ss.high <- 55000
        #       }
        #     }
        # }  
      }else{  # This is for cases such as BCM "Panel2_M_174601_015.labelled.fcs"
        ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.005)
      }
      
    }
  }
  
  if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
    singlets.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-A'], scat.chans['FSC-H']), position=c(F, NA), gates = c(ss.high, NA))
    #plotDens(FSCsinglets.flowD, c(scat.chans['SSC-A'], scat.chans['FSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=singlets.peak-2000)
    # lines(singlets.flowD@filter)
  }else{
    singlets.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(F, NA), gates = c(ss.high, NA))
    # plotDens(FSCsinglets.flowD, c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=singlets.peak-2000)
    # lines(singlets.flowD@filter)
  }
  
  
  
  # }else{
  #   singlets.flowD <- FSCsinglets.flowD
  # }
  
  # Return results ---------------------------------------------------------------------------------------------
  results <- list()
  results$scat.chans <- scat.chans
  ## Returing the flowDensity objects for event counts + proportions
  results$live <- live.flowD
  results$size <- size.flowD
  results$FSCsinglets <- FSCsinglets.flowD
  results$singlets <- singlets.flowD
  results$dead.peak <- dead.peak
  
  ## Returing the gating thresholds
  results$fsca.live.gate <- fsca.live.gate
  results$live.gate <- live.gate
  results$ss.high <- ss.high
  
  remove(temp.flowD)
  gc()
  return(results)
}


## End of the qualityGate function