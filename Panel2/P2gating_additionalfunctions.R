
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
  fsca.live.gate <- c(deGate(temp.flowD@flow.frame, channel = c('FSC-A'), all.cuts = T, tinypeak.removal = 0.02, upper = F, percentile = NA),
                      deGate(temp.flowD@flow.frame, channel = c('FSC-A'), use.upper = T, upper = F, tinypeak.removal = 0.5))
  
  # Find FSC-A peak and pick gate that is just below it 
  maxDens <- density(getflowFrame(temp.flowD)@exprs[,c('FSC-A')])
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
  fsca.live.gate <- max(fsca.live.gate[which(fsca.live.gate < (fscapeak.lcn - 7500))])
  
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
  
    live.gate <- deGate(f, channel = c(channels.ind["Live"]), all.cuts = T, tinypeak.removal = 0.001)
    live.gate <- live.gate[which(live.gate > maxDens2$x[which.max(maxDens2$y)])]
    
    if(length(live.gate) > 0){ # if there are any gates above the live peak
      live.gate <- live.gate[which(live.gate < dead.peak)]
      if(length(live.gate) > 0){
        live.gate <- live.gate[which.min(abs(dead.peak - 0.5 - live.gate))]
      }else{
        live.gate <- dead.peak - 0.5
      }
    }else{ # For CIPHE 16-Oct-04 there are no dead cells, so there are no gates found above the live peak
      live.gate <- deGate(f, channel = c(channels.ind["Live"]), use.upper = T, upper = T, alpha = 0.005, tinypeak.removal = 0.2)
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
  
  # Gating size ----------------------------------------------------------------------------------------
  size.flowD <- live.flowD
  
  # Gating FSC singlets --------------------------------------------------------------------------------
  
  # Check if FSC-H and FSC-W channel names are switched
  temp <- getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
  try({
    FSC.angle <- atan(summary(lm(temp[, c('FSC-H')] ~ temp[, c('FSC-A')]))$coefficients[2,1])
    if(FSC.angle < pi/9){ # assume that the FSC-H and FSC-W channels are interchanged
      scat.chans['FSC-H'] <- grep('FSC-W', colnames(f))
      scat.chans['FSC-W'] <- grep('FSC-H', colnames(f))
      temp <- getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
    }
  })
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
  
  # rotate data
  temp.flowD <- flowDensity(size.flowD@flow.frame, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA,F), gates = c(NA, 0.975*max(getflowFrame(size.flowD)@exprs[,c(scat.chans['FSC-H'])])))
  rot <- rotate.data(getflowFrame(temp.flowD), c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = -theta0)$data
  
  maxDens <- density(rot@exprs[,c(scat.chans['FSC-H'])])
  maxDens <- smooth.spline(maxDens$x, maxDens$y, spar = 0.25)
  
  gate <- deGate(rot, channel = c(scat.chans['FSC-H']), all.cuts = T, tinypeak.removal = 0.001, percentile = 0)
  gate0 <- gate
  singlets.peak <- maxDens$x[which.max(maxDens$y)]
  if(is.null(names(gate))){ 
    idx <- which(gate < (singlets.peak - 10000))
    if(length(idx > 0)){
      gate <- gate[idx]
      p05percentile <- min(deGate(rot, channel = c(scat.chans['FSC-H']), use.percentile = T, percentile = 0.0005), 50000) 
      gate <- gate[which.min(abs(singlets.peak - 0.45*(singlets.peak - p05percentile) - gate))]
    }else{
      gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.01)
    }
  }
  
  if((singlets.peak - gate) > 50000){
    gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.5) - 15000
  }
  
  temp <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, gate))
  #   , ellip.gate = T, scale = .9999999);
  temp@filter <- rotate.data(temp@filter, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
  temp@flow.frame <- rotate.data(getflowFrame(temp), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  # Gate really closesly - need this to tell if SSC-H/SSC-W channels are interchanged 
  temp2 <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, singlets.peak - 2000))
  temp2@flow.frame <- rotate.data(getflowFrame(temp2), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  FSCsinglets.flowD <- temp
  
  
  # Gating SSC singlets -------------------------------------------------------------------------
  
  #if(!is.na(scat.chans['SSC-W'])){ # skip second singlet gate if no SSC-W   ## This if condition not needed with the files we are running on right now
    
  temp <- getflowFrame(temp2)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])]
  
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
    
    #This works for CIPHE/TCP/JAX
    ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.4, all.cuts = T)
    
    #ss.high2 <- NA
    if(!is.null(names(ss.high))){ # if there is a 95% gate
      ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.9, use.upper = T, upper = T, alpha = 0.9)
      maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[,c(scat.chans['SSC-W'])])
      
      if((ss.high - maxDens$x[which.max(maxDens$y)]) < 50000){
        temp.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(T, NA), gates = c(ss.high - 12000, NA))
        ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02)
        if(!is.null(names(ss.high))){
          ss.high <- deGate(temp.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.02, use.upper = T, upper = T) + 3000
        }
        maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-H'])])
        temp2.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(NA, T), gates = c(NA, (maxDens$x[which.max(maxDens$y)] + 20000)))
        ss.high2 <- deGate(temp2.flowD@flow.frame, channel = c(scat.chans['SSC-W']), upper = T, tinypeak.removal = 0.4, percentile = NA)
        ss.high <- c(ss.high, ss.high2)
        #maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'])])
        ss.high <- min(ss.high)
      }else{  # This is for cases such as BCM "Panel2_M_174601_015.labelled.fcs"
        ss.high <- deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), tinypeak.removal = 0.005)
      }
      
    }
  }
  
  singlets.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(F, NA), gates = c(ss.high, NA))
  
  # }else{
  #   singlets.flowD <- FSCsinglets.flowD
  # }
  
  # Return results ---------------------------------------------------------------------------------------------
  results <- list()
  results$live <- live.flowD
  results$FSCsinglets <- FSCsinglets.flowD
  results$singlets <- singlets.flowD
  results$scat.chans <- scat.chans
  
  # # temporary
  # results$rot <- rot
  # results$gate <- gate0
  # results$ss.high <- ss.high
  # results$ss.high2 <- ss.high2
  
  return(results)
}

compensateIMPC <- function(ff, fname.full, centre, panel.no){
  # Compensates the flow frame data
  #
  # Args:
  #   ff: flow frame
  #   fname.full: FCS file full path 
  #   centre: TCP, CIPHE, BCM or Jax
  #   panel.no: integer (1 for panel 1, 2 for panel 2)
  #
  # Returns: a compensated flow frame
  
  fname <- basename(fname.full)
  compMatrix.dir <- dirname(fname.full)
  
  if(centre == 'BCM'){
    
    compMatrix <- dir(compMatrix.dir, full.names = T, recursive = F, pattern = "Comp.*csv")
    compMatrix <- read.csv(compMatrix, stringsAsFactors = F)
    
    chanNames <- str_match(compMatrix[,1], pattern = '(?<=:: ).*|DAPI-A')
    idx <- which(is.na(chanNames))
    
    if(length(idx) < 7){
      chanNames[idx] <- compMatrix[idx, 1]
      chanNames <- gsub('-A','', chanNames)
      chanNames <- gsub('.*-','', chanNames)
      chanNames <- gsub('_', '/', chanNames)
      chanNames <- gsub('BV786', 'CD23', chanNames)  # There are actually 2 different BV786 labels in pData(parameters(f))$desc. One of which is CD23, I hope this is correct.
      channel.ind <- sapply(chanNames, function(x){
        if(x == 'CD4'){ 
          ind<- grep('CD4(?!4)', pData(parameters(ff))$desc, perl = T)
        }else{
          ind <- grep(x, pData(parameters(ff))$desc)
          if(length(ind) < 1){
            if(x == 'IA/E'){
              ind <- grep("MHCII", pData(parameters(ff))$desc)
            }
          }
        }
        return(ind)
      })
      channel.ind <- unlist(channel.ind)
    }else{ 
      # For P1, there are four days (150123, 150505, 150512, 160920)
      # where the compensation matrix row/column names have fluorophore names instead of Ab names
      # The following channel.ind values work for those days
      channel.ind <- 7:14
    }
    
    compMatrix <- compMatrix[,-c(1)]
    colnames(compMatrix) <- colnames(ff)[channel.ind]
    rownames(compMatrix) <- colnames(ff)[channel.ind]
    ff <- compensate(ff, compMatrix)
    
  }else if (centre == "TCP"){
    
    file.names.0909A <- c("PANEL_A_ABIX_223_C7C07009.fcs", "PANEL_A_ABOV_126_C6C06008.fcs",
                          "PANEL_A_ABQX_108_C4C04006.fcs", "PANEL_A_ABRB_82_C2C02004.fcs",
                          "PANEL_A_ABRB_85_C3C03005.fcs",  "PANEL_A_ABRL_111_C5C05007.fcs",
                          "PANEL_A_ABUN_78_C1C01003.fcs",  "PANEL_A_B6NC_C8C08010.fcs")
    
    file.names.0909B <- c("PANEL_B_ABIX_223_F7F07022.fcs", "PANEL_B_ABOV_126_F6F06021.fcs",
                          "PANEL_B_ABQX_108_F4F04019.fcs", "PANEL_B_ABRB_82_F2F02017.fcs",
                          "PANEL_B_ABRB_85_F3F03018.fcs",  "PANEL_B_ABRL_111_F5F05020.fcs",
                          "PANEL_B_ABUN_78_F1F01016.fcs",  "PANEL_B_B6NC_F8F08023.fcs")
    
    file.names.0915 <- c("PANEL_A_ABMX_153_C7C07009.fcs", "PANEL_A_ABNL_105_C3C03005.fcs",
                         "PANEL_A_ABNL_106_C4C04006.fcs", "PANEL_A_ABRB_94_C1C01003.fcs",
                         "PANEL_A_ABSU_98_C2C02004.fcs",  "PANEL_A_ABTV_105_C5C05007.fcs",
                         "PANEL_A_ABTV_108_C6C06008.fcs", "PANEL_A_B6NC_736_C8C08010.fcs",
                         "PANEL_B_ABMX_153_F7F07021.fcs", "PANEL_B_ABNL_105_F3F03017.fcs",
                         "PANEL_B_ABNL_106_F4F04018.fcs", "PANEL_B_ABRB_94_F1F01015.fcs",
                         "PANEL_B_ABSU_98_F2F02016.fcs",  "PANEL_B_ABTV_105_F5F05019.fcs",
                         "PANEL_B_ABTV_108_F6F06020.fcs", "PANEL_B_B6NC_736_F8F08022.labelled.fcs")
    
    file.names.20160526A <- c("PANEL_A_ACWH_243_C1C01003.fcs", "PANEL_A_ACWH_244_C2C02004.fcs",
                              "PANEL_A_ACWH_245_C3C03011.fcs", "PANEL_A_ACWH_247_C4C04007.fcs",
                              "PANEL_A_ACWH_255_C5C05008.fcs", "PANEL_A_B6NC_950_C6C06009.fcs",
                              "PANEL_A_B6NC_955_C7C07010.fcs")
    
    file.names.20160526B <- c("PANEL_B_ACWH_243_F1F01002.fcs", "PANEL_B_ACWH_244_F2F02003.fcs",
                              "PANEL_B_ACWH_245_F3F03004.fcs", "PANEL_B_ACWH_247_F4F04005.fcs",
                              "PANEL_B_ACWH_255_F5F05006.fcs", "PANEL_B_B6NC_950_F6F06007.fcs",
                              "PANEL_B_B6NC_955_F7F07008.fcs")
    
    file.names.20160623A <- c("PANEL_A_ACRC_43_C1C01003.fcs", "PANEL_A_ACRK_62_C2C02004.fcs", "PANEL_A_ACRK_65_C3C03005.fcs")
    
    file.names.20160623B <- c("PANEL_B_ACRC_43_F1F01012.fcs", "PANEL_B_ACRK_62_F2F02013.fcs", "PANEL_B_ACRK_65_F3F03014.fcs")
    
    
    if(fname %in% file.names.0909A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/COMPENSATION_CONTROL_2015-09-09_PanelA.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.0909B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Controls_2015-09-09_PanelB.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.0915){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/COMPENSATION_2015-09-15_PANELA&B.csv", check.names = F)[,-1]
      if(panel.no == 1) # one channel name is different from Panel 1 to 2.
        colnames(SPILL.matrix)[8] <- "BV786-A"
    }else if(fname %in% file.names.20160526A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/MATRIX_PANELA_2016-05-26.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.20160526B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/MATRIX_PANELB_2016-05-26.csv", check.names = F)[,-1]
      colnames(SPILL.matrix) <- gsub("(.*) ::.*", '\\1', colnames(SPILL.matrix))
      rownames(SPILL.matrix) <- colnames(SPILL.matrix)
    }else if(fname %in% file.names.20160623A){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Matrix_2016-06-23_PANEL_A.csv", check.names = F)[,-1]
    }else if(fname %in% file.names.20160623B){
      SPILL.matrix <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/Tcp_extra/Compensations/Originals/Compensation_Matrix_2016-06-23_PANEL_B.csv", check.names = F)[,-1]
      colnames(SPILL.matrix) <- gsub("(.*) ::.*", '\\1', colnames(SPILL.matrix))
      rownames(SPILL.matrix) <- colnames(SPILL.matrix)
    }else{
      if(det(ff@description$SPILL) == 1){
        cat(paste0("Check the spillover matrix for ", fname, "it's probably an identity matrix!", "\n"))
        SPILL.matrix <- ff@description$SPILL
      }else{
        SPILL.matrix <- ff@description$SPILL
      }
    }
    ff <- compensate(ff, SPILL.matrix)
    
  }else if(centre == 'Jax'){
    
    if(panel.no == 1){
      compMatrix <- dir(compMatrix.dir, full.names = T, recursive = F, pattern = "CompMatrix_Tmem")
    }else{
      compMatrix <- dir(compMatrix.dir, full.names = T, recursive = F, pattern = "CompMatrix_APC")
    }
    
    # Reading the file depends on whether the compensation matrix data is saved as a csv file or not
    if(length(grep(".csv", compMatrix)) > 0){ 
      compMatrix <- read.csv(compMatrix, stringsAsFactors = F)
      chanNames <- compMatrix[2, 2:ncol(compMatrix)]
      if(compMatrix.dir == '/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/Jax_2017-03-29/PKG - JAXKOMP data files/16-0602  KOMP'){
        compMatrix <- compMatrix[3:(2 + length(chanNames)), 2:ncol(compMatrix)]
      }else{
        compMatrix <- compMatrix[5:(4 + length(chanNames)), 2:ncol(compMatrix)]
      }
      compMatrix <- sapply(compMatrix, as.numeric)
      rownames(compMatrix) <- chanNames
      colnames(compMatrix) <- chanNames
    }else{ 
      chanNames <- read.table(compMatrix, skip = 2, nrows = 1, sep ='\t', stringsAsFactors =F)
      compMatrix <- read.table(compMatrix, sep = '\t', skip = 3, nrows = length(chanNames), stringsAsFactors = F)
      rownames(compMatrix) <- chanNames
      colnames(compMatrix) <- chanNames
    }
    
    index.Remove.row <- which(is.na(compMatrix[,1]))
    index.Remove.col <- which(is.na(compMatrix[1,]))
    if(length(index.Remove.row) != 0 & length(index.Remove.col) != 0){
      compMatrix <- compMatrix[-index.Remove.row,]
      compMatrix <- compMatrix[,-index.Remove.col]
    }

    ff <- compensate(ff, compMatrix)
    
  }else{  # CIPHE
    
    if(det(ff@description$SPILL) == 1){
      cat(paste0("Check the spillover matrix, it's probably an identity matrix!", "\n"))
      ff <- compensate(ff, ff@description$SPILL)
    }else{
      ff <-compensate(ff, ff@description$SPILL)
    }
  }
  
  return(ff)
}
 

##Finds markers in the FCS file
Find.markers <- function(frame,marker.list){
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}

removeMargins<- function(f,chans,sens=1, debris=FALSE,return.ind=F,neg=500, verbose = T){
  neg <-cbind(chans,neg)[,2]
  #Size is a vector of size 2, to be passed to mfrow in case of plotting
  data <- exprs(f)
  margins <- c()
  marg.list <-list()
  if(!debris)
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.max <-max(data[,chan])
      margins <- which ( data[, chan] >= stain.max*sens)
      marg.list <- c(marg.list, list(margins))
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
    
  }else
  {
    for(chan in chans)
    {
      stain.min <-min(data[,chan])
      margins <- which ( data[, chan] <= stain.min*sens)
      marg.list <- c(marg.list, list(margins))
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
  }
  for(i in 1:length(chans))
  {
    if (neg[i]<500)
    {
      negs <- which ( data[, chans[i]] < neg[i])
      margins <- negs
    }
    marg.list <- c(marg.list, list(margins))
    if (length(margins)!=0){
      data <- data[ -margins, ]
    }
    if(verbose == T){print(paste(length(margins), "negative events in",colnames(f)[chans[i]], "will be removed.",sep =" "))}
  }
  exprs(f) <- data
  if (!return.ind)
    return(f)
  else
    return(list(frame=f,ind=marg.list))
  
}

rotate.data <- function(data, chans=NULL, theta=NULL){
  if (class(data)== "flowFrame" & !is.null(chans))
  {
    data.new <- exprs(data)[,chans]
    if (is.null(theta))
    {
      reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
      theta <- pi/2 - reg.slope
    }
    data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    exprs(data)[,chans] <- data.new
  }else{
    data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
  }
  return(list(data=data,theta=theta))
}

# Written by Quentin

logiclTransformCiphe <- function(flow.frame, markers.transform){
  
  #######################################################################################################
  #
  #
  #
  #
  #
  #######################################################################################################
  
  #no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
 # markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x) 
      as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
    )	
  } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x) 
      as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
    )	
  } else 
  {
    r.values <- rep(90, length(list.index))
  }
  
  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }
  
  return(flow.frame)
}