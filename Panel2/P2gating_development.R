rm(list = ls())

library('colorRamps')
library('plyr')
library('doMC')
library('e1071') # for flowCut
library('Cairo') # for flowCut
library('flowDensity')
library('pracma') # for findpeaks
library('stringr')
library('flowPeaks')
library('MASS')

setwd('/data/IMPC')

no_cores <- detectCores() - 2
registerDoMC(no_cores)

source("/data/IMPC_Pipeline/helperfunc.R")
source("/data/Universal_Codes/flowCut.R")

# load IMPReSS ID info ------------------------------------------------------------
load('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/immpressIDconversion.Rdata')
colnames(immpressIDconversion) <- c("IMPReSS.id", "Genotype")

# TCP --------------------------------------------------------------------------------
TCPdata.dir <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/TCP_170209/UPLOAD_DUMP"
allFilesP1 <- dir(TCPdata.dir, recursive = T, pattern="PANEL_A", full.names = T)
allFilesP1 <- gsub("__", "_", allFilesP1)
allFilesP1 <- allFilesP1[-grep("FMO",allFilesP1)]

allFilesP2 <- dir(TCPdata.dir, recursive = T, pattern="PANEL_B", full.names = T)
allFilesP2 <- gsub("__", "_", allFilesP2)
allFilesP2 <- allFilesP2[-grep("FMO",allFilesP2)]

corrupted.files <- grep('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/TCP_170209/UPLOAD_DUMP/2016-01-07/PANEL_B_ACFK_63_F1F01014.fcs', allFilesP2)
allFilesP2 <- allFilesP2[-corrupted.files]

allFilesP1_TCP <- allFilesP1
allFilesP2_TCP <- allFilesP2

# BCM --------------------------------------------------------------------------------
allFiles <- dir("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/BCM/BCM", full.names = T, recursive = T, pattern = ".fcs")

allFiles <- allFiles[-grep("FMO",allFiles)]
allFiles <- allFiles[-grep("Neg",allFiles)]
allFiles <- allFiles[-grep("Beads",allFiles)]
allFiles <- allFiles[-grep("DAPI",allFiles)]
allFiles <- allFiles[-grep("SingleStaincells",allFiles)]
allFiles <- allFiles[-grep("_Nikolai",allFiles)]
allFiles <- allFiles[-grep("WT1",allFiles)]
allFiles <- allFiles[-grep("WT2",allFiles)]
allFiles <- allFiles[-grep("WT3",allFiles)]
allFiles <- allFiles[-grep("WT4",allFiles)]
allFiles <- allFiles[-grep("KO1",allFiles)]
allFiles <- allFiles[-grep("KO2",allFiles)]
allFiles <- allFiles[-grep("KO3",allFiles)]
allFiles <- allFiles[-grep("KO4",allFiles)]

allFiles2 <- sapply( allFiles, function(x) {strsplit(x, split = "/")[[1]][length(strsplit(x, split = "/")[[1]])]}); names(allFiles2) <- NULL

allFilesPanel1 <- allFiles[grep("Panel1",allFiles2)]
allFilesPanel2 <- allFiles[grep("Panel2",allFiles2)]

allFilesP1_BCM <- allFilesPanel1
allFilesP2_BCM <- allFilesPanel2

# BCM.datadir <- '/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/BCM/BCM Test/'
# allFiles <- dir(BCM.datadir, full.names = T, recursive = T, pattern = ".fcs")
# allFilesP2_BCM <- allFiles[grep("Panel2", allFiles)]
# allFilesP2_BCM <- allFilesP2_BCM[-grep("FMO", allFilesP2_BCM)]
# allFilesP2_BCM <- allFilesP2_BCM[-grep("Neg", allFilesP2_BCM)]
# allFilesP2_BCM <- allFilesP2_BCM[-grep("DAPI", allFilesP2_BCM)]
# allFilesP2_BCM <- allFilesP2_BCM[-grep("Nikolai", allFilesP2_BCM)]

# CIPHE --------------------------------------------------------------------------------
allFiles <- dir("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/", full.names = T, recursive = T, pattern = ".fcs")
allFilesP1_CIPHE <- allFiles[grep("IMPC1",allFiles)]
allFilesP2_CIPHE <-allFiles[grep("IMPC2",allFiles)]

# Jax--------------------------------------------------------------------------------
# allFilesP1_Jax <- dir('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/PKG-JAX_KOMP_P1', full.names = T, recursive = T, pattern = ".fcs")
# allFilesP2_Jax <- dir('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/PKG-JAX_APC_P2', full.names = T, recursive = T, pattern = ".fcs")

allFiles_Jax <- dir('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/Jax_2017-03-29/PKG - JAXKOMP data files', full.names = T, recursive = T, pattern = ".fcs")
allFilesP1_Jax <- allFiles_Jax[grep("/Tmem",allFiles_Jax)]
allFilesP2_Jax <-allFiles_Jax[grep("/APC",allFiles_Jax)]

# Now I have all FCS data; singlet gate ---------------------------------------------------------------

qualityGate <- function(f, scat.chans, channels.ind){
  # Applies quality gates on the flow frame f. Returns the live singlets.
  #
  # Args:
  #   f: flow frame
  #   scat.chans: named vector with scatter channel indices
  #   channels.ind: namesd vector with the other channel indices, channels.ind['Live'] is required
  #
  # Returns: A list containing
  #       scat.chans: updated scatter channel indices (some may be switched)
  #       live: a flowDensity object containing live cells
  #       FSCsinglets: a flowDensity object containing the live/FSC singlets
  #       singlets: a flowDensity object containing live/FSC singlets/SSC singlets 
  
  # Gating Live ---------------------------------------------------------------------------------------
  
  temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(NA, F), gates = c(NA, 1.5))
  fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c('FSC-A'), all.cuts = T, tinypeak.removal = 0.02, upper = F, percentile = NA)
  
  # Find FSC-A peak and pick gate that is just below it 
  maxDens <- density(getflowFrame(temp.flowD)@exprs[,c('FSC-A')])
  fscapeak.lcn <- findpeaks(maxDens$y)
  fscapeak.lcn <- fscapeak.lcn[which(maxDens$x[fscapeak.lcn[,2]] > 35000),]
  if(length(fscapeak.lcn) > 4){
    #fscapeak.lcn <- maxDens$x[max(fscapeak.lcn[which(fscapeak.lcn[,1] > 0.2*max(fscapeak.lcn[,1])), 2])]
    fscapeak.lcn <- maxDens$x[fscapeak.lcn[which(fscapeak.lcn[,1] > 0.2*max(fscapeak.lcn[,1])), 2]]
    fscapeak.lcn <- fscapeak.lcn[which.min(abs(fscapeak.lcn - 100000))]
  }else{
    fscapeak.lcn <- maxDens$x[fscapeak.lcn[2]]
  }
  
  if(min(fsca.live.gate) > fscapeak.lcn){
    fsca.live.gate <- deGate(temp.flowD@flow.frame, channel = c('FSC-A'), use.upper = T, upper = F)
  }
  
  fsca.live.gate <- max(fsca.live.gate[which(fsca.live.gate < fscapeak.lcn)])
  
  # Find the peak corresponding to the dead cells
  maxDens <- density(f@exprs[,c(channels.ind['Live'])])
  peak.lcns <- findpeaks(maxDens$y)
  # dead.peak <- peak.lcns[which(maxDens$x[peak.lcns[, 2]] > 2.7),]
  dead.threshold <- max(deGate(f, channel = c(channels.ind['Live']), use.percentile = T, percentile = 0.85), 2.4)
  dead.threshold <- min(dead.threshold, 2.8) # needed if there are lots of dead cells
  dead.peak <- peak.lcns[which(maxDens$x[peak.lcns[, 2]] > dead.threshold),]
  if(length(dead.peak) > 4){
    dead.peak <- maxDens$x[dead.peak[which.max(dead.peak[,1]), 2]]
  }else{
    dead.peak <- maxDens$x[dead.peak[2]]
  }
  
  is.DAPI <- ((length(grep('DAPI', pData(parameters(f))$desc)) > 0) | (length(grep('DAPI', pData(parameters(f))$name)) > 0))
  
  if(!is.DAPI){
    live.gate <- deGate(f, channel = c(channels.ind["Live"]), all.cuts = T, tinypeak.removal = 0.001)
    #live.gate <- max(live.gate[which(live.gate < (dead.peak - 0.2))])
    live.gate <- live.gate[which(live.gate < dead.peak)]
    live.gate <- live.gate[which.min(abs(dead.peak - 0.5 - live.gate))]
    
  }else{
    temp.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, NA), gates = c(fsca.live.gate + 10000, NA))
    #live.gate <- deGate(temp.flowD@flow.frame, channel = c(channels.ind['Live']), use.upper = T, upper = T, alpha = 0.01, tinypeak.removal = 0.1) 
    maxDens <- density(getflowFrame(temp.flowD)@exprs[,c(channels.ind['Live'])])
    start.idx <- which.max(maxDens$y)
    stop.idx <- which.min(abs(maxDens$x - dead.peak))
    valley.yval <- min(maxDens$y[start.idx:stop.idx])
    live.gate <- min(maxDens$x[which(maxDens$y[start.idx:stop.idx] < 2*valley.yval) + (start.idx -1)])
  }
  
  live.flowD <- flowDensity(f, channels = c(scat.chans['FSC-A'], channels.ind['Live']), position = c(T, F), gates = c(fsca.live.gate, live.gate))
  
  # Gating size ----------------------------------------------------------------------------------------
  size.flowD <- live.flowD
  
  # Gating FSC singlets --------------------------------------------------------------------------------
  
  # Check if FSC-H and FSC-W channel names are switched
  temp <- getflowFrame(size.flowD)@exprs[, c(scat.chans['FSC-A'], scat.chans['FSC-H'])]
  try({
    FSC.angle <- atan(summary(lm(temp[, c('FSC-H')] ~ temp[, c('FSC-A')]))$coefficients[2,1])
    if(FSC.angle < pi/8){ # assume that the FSC-H and FSC-W channels are interchanged
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
  #x0 <- deGate(temp.flowD@flow.frame, channel = c("FSC-A"), use.percentile = T, percentile = 0.925)
  x0 <- 125000
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
  if(is.null(names(gate))){ 
    idx <- which(gate < maxDens$x[which.max(maxDens$y)])
    if(length(idx > 0)){
     gate <- max(gate[idx])
    }else{
      gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.01)
    }
  }
  
  if((maxDens$x[which.max(maxDens$y)] - gate) > 40000){
    gate <- deGate(rot, channel = c(scat.chans['FSC-H']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.5)
  }
  
  
  temp <- flowDensity(rot, channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), position = c(NA, T), gates = c(NA, gate))
  #   , ellip.gate = T, scale = .9999999);
  temp@filter <- rotate.data(temp@filter, c(scat.chans['FSC-A'], scat.chans['FSC-H']), theta = theta0)$data
  temp@flow.frame <- rotate.data(getflowFrame(temp), c(scat.chans['FSC-A'], scat.chans['FSC-H']),theta = theta0)$data
  
  FSCsinglets.flowD <- temp
  
  
  # Gating SSC singlets -------------------------------------------------------------------------
  
  if(!is.na(scat.chans['SSC-W'])){ # skip second singlet gate if no SSC-W
    
    # Check if SSC-H and FSC-W channel names are switched
    temp <- getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])]
    SSC.angle <- atan(summary(lm(temp[, c('SSC-H')] ~ temp[, c('SSC-W')]))$coefficients[2,1])
    if(SSC.angle < pi/16){ # assume that the SSC-H and SSC-W channels are interchanged
      scat.chans['SSC-H'] <- grep('SSC-W', colnames(f))
      scat.chans['SSC-W'] <- grep('SSC-H', colnames(f))
    }
    
    ss.high <- c(deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W']), use.upper = T, upper = T, tinypeak.removal = 0.4),
                 deGate(FSCsinglets.flowD@flow.frame, channel = c(scat.chans['SSC-W'])))
    maxDens <- density(getflowFrame(FSCsinglets.flowD)@exprs[ , c(scat.chans['SSC-W'])])
    ss.high <- min(ss.high[which(ss.high > maxDens$x[which.max(maxDens$y)])])
    singlets.flowD <- flowDensity(FSCsinglets.flowD, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), position=c(F, NA), gates = c(ss.high, NA))
  }else{
    singlets.flowD <- FSCsinglets.flowD
  }
  
  # Return results ---------------------------------------------------------------------------------------------
  results <- list()
  results$live <- live.flowD
  results$FSCsinglets <- FSCsinglets.flowD
  results$singlets <- singlets.flowD
  results$scat.chans <- scat.chans
  
  return(results)
}

compensateIMPC <- function(ff, fname.full, centre, panel.no){
  
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
    
  }else{  # all other centers (Jax, CIPHE)
    
    if(det(ff@description$SPILL) == 1){
      cat(paste0("Check the spillover matrix, it's probably an identity matrix!", "\n"))
      ff <- compensate(ff, ff@description$SPILL)
    }else{
      ff <-compensate(ff, ff@description$SPILL)
    }
  }
  
  return(ff)
}

switchCenterData <- function(type){
  switch(type,
         Jax = { 
           #allFCS <<- c(allFilesP2_Jax[c(685)])
           allFCS <<- allFilesP2_Jax
           scatterplot.dir <<- "/data/IMPC/ScatterPlots/Jax"
         },
         TCP = { 
           allFCS <<- c(allFilesP2_TCP[c(50:55, 200:203, 455:460, 500:502, 683:689, 700:710, 1000:1005)])
           scatterplot.dir <<- "/data/IMPC/ScatterPlots/TCP"
         },
         BCM = { 
           allFCS <<- c(allFilesP2_BCM[c(168:172, 300:305, 472:474, 500:503, 670:673, 801:810, 1200:1210, 1250:1253, 1261, 1400:1422)])
           scatterplot.dir <<- "/data/IMPC/ScatterPlots/BCM"
         },
         CIPHE = {
           #allFCS <<- c(195, 800, 819:832) # problematic files '16-Apr-12';'15-Sep-15_IMPC2_08'; '14-Jul-22_IMPC2_11'
           allFCS <<- c(allFilesP2_CIPHE[c(168:172, 300:305, 472:474, 500:503, 670:673, 801:810, 1200:1210, 1250:1253, 1261)])
           #allFCS <<- c(allFilesP2_CIPHE[c(819:832, 1200:1205)])
           scatterplot.dir <<- "/data/IMPC/ScatterPlots/CIPHE"
         }
  )
}


for(j in c('TCP', 'CIPHE', 'Jax', 'BCM')){
#for(j in c('Jax', 'BCM')){
  
  switchCenterData(j)
  print(scatterplot.dir)
  
  foreach(i = 1:length(allFCS)) %dopar% {
  #for(i in 28:length(allFCS)){ 
    try({
      
      f <- read.FCS(allFCS[i])
      fname <- basename(allFCS[i])
      fname.full <- allFCS[i]
      print(fname)
      
      # Get channel and scatter indices -------------------------------------------------------------------------------
     
      # Panel 2
      markers <- c("Live|I515*|Syto*|DAPI", "CD5|Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                   "CD11c", "CD21", "CD23")
      channels.ind <-Find.markers(f, markers)
      names(channels.ind)[grep(names(channels.ind), pattern = "Live*")] <- "Live"
      names(channels.ind)[grep(names(channels.ind), pattern = "*II")] <- "MHCII"
      names(channels.ind)[grep(names(channels.ind), pattern = "*Ly6G")] <- "CD5/Ly6G"
      names(channels.ind)[grep(names(channels.ind), pattern = "*6C")] <- "Ly6C"
      names(channels.ind)[grep(names(channels.ind), pattern = "F4")] <- "F4/80"
      names(channels.ind)[grep(names(channels.ind), pattern = "CD21")] <- "CD21/CD35"
      if(is.na(channels.ind["Live"])){  
        # BCM (new data, DAPI is under the $name descriptor 
        if(length(grep('DAPI', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('DAPI', pData(parameters(f))$name)
        }else{
          # For TCP, Sytox Blue is read in the BV-510 channel
          channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
        }
      }
      if(is.na(channels.ind["Time"])){  
        channels.ind["Time"] <- grep('Time', pData(parameters(f))$name)
      }
      
      
      scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
      names(scat.chans) <- colnames(f)[scat.chans]
      
      # Remove scatter margins and compensate ---------------------------------------------
      # Removing margin events in Scatter channels
      f <- removeMargins(f, chans = scat.chans, verbose = F)
      #Removing negative values in scatter channels
      f <- removeMargins(f, chans = scat.chans, debris = T, neg = T, verbose = F)
      

      if(j == 'Jax'){
        compMatrix <- read.csv('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/Jax/PKG - Jax comp info/CompMatrix_APC 15-0129.csv', stringsAsFactors = F, header = T)
        markers.compMatrix <- compMatrix[3,2:ncol(compMatrix)]
        markers.compMatrix <- gsub('\\+.*', '', markers.compMatrix)
        markers.compMatrix<- gsub('CD21/35', 'CD21_35', markers.compMatrix)
        markers.compMatrix<- gsub('F4/80', 'F4-80', markers.compMatrix)
        markers.compMatrix<- gsub('KLRG-1', 'KLRG1', markers.compMatrix)
        names(markers.compMatrix) <- markers.compMatrix
        compMatrix <- compMatrix[5:16,2:ncol(compMatrix)]
        compMatrix <- apply(compMatrix, 2, as.numeric)
        channel.ind <- sapply(markers.compMatrix, function(x){
          ind <- grep(x, pData(parameters(f))$desc)
          if(length(ind) < 1){
            if(x == 'DAPI'){
              ind <- grep(x, pData(parameters(f))$name)
            }
          }
          return(ind)
        })
        channel.ind <- unlist(channel.ind)
        colnames(compMatrix) <- colnames(f)[channel.ind]
        rownames(compMatrix) <- colnames(f)[channel.ind]
        f <- compensate(f, compMatrix)

        #f <- compensate(f, f@description$SPILL)
      }else{ # CIPHE, TCP and BCM
        f <- compensateIMPC(f, fname.full, centre = j, panel.no = 2)
      }
      
      # Transformation
      if(j == 'Jax'){ # just haven't set up a global transformation matrix for BCM yet
        no.transform <- union(scat.chans, channels.ind['Time'])
        lgl <- estimateLogicle(f, channels = colnames(f)[-no.transform])
      }else{
        load(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/', j, "/Panel2/Results/lgl.Rdata"))
      }
      
      fT <- transform(f, lgl)
      #fT <- logiclTransformCiphe(f)
      f.Clean <- flowCut(fT, segment = floor(nrow(f)*5/3000), FileID = fname, directory = '/data/IMPC_PanelB/flowCut/', Plot = 'Failed Only')
      
      f <- transform(f, lgl)
      #f <- logiclTransformCiphe(f)
      if(length(f.Clean$ind) > 0){
        f@exprs <- f@exprs[-f.Clean$ind, ]
      }
      
      
      # Quality Gate (live/dead, singlets gating) ---------------------------------------------------
      
      results <- qualityGate(f, scat.chans, channels.ind)
      live.flowD <- results$live
      FSCsinglets.flowD <- results$FSCsinglets
      singlets.flowD <- results$singlets
      scat.chans <- results$scat.chans
      
      # Panel 2 Gating -------------------------------------------------------------------------------
      
      singlets <- getflowFrame(singlets.flowD)
      
      # Gating singlets to get Granulocytes ------------------------------------------------------------
      
      CD5.gate <- deGate(singlets, channel = channels.ind['CD5/Ly6G']) - 0.3
      CD5pos.flowD <- flowDensity(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, NA), gates = c(CD5.gate, NA))
      
      temp <- getflowFrame(CD5pos.flowD)@exprs[, c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"])]
      flowPeaks.Res <- flowPeaks(temp)
      cluster.size <- sapply(flowPeaks.Res$peaks$cid, function(x){ length(which(flowPeaks.Res$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res$peaks$cid[which(cluster.size > 0.001)]
      
      mu.short <- flowPeaks.Res$peaks$mu[sig.clusters,]
      granulo.clusterid <- sig.clusters[which.max(5*mu.short[,2]^2 + mu.short[,1]^2)]
      granulo.idx <- which(flowPeaks.Res$peaks.cluster ==  granulo.clusterid)
      granulocytes <- getflowFrame(CD5pos.flowD)
      granulocytes@exprs <- granulocytes@exprs[granulo.idx, ]
      
      CD11b.granulo.gate <- min(granulocytes@exprs[, channels.ind["CD11b"]])
      granulocytes.flowD <- flowDensity(granulocytes, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, T), gates = c(min(granulocytes@exprs[, channels.ind["CD5/Ly6G"]]), CD11b.granulo.gate) , ellip.gate = T, alpha = 0.85)
      CD11b.granulo.gate <- deGate(granulocytes.flowD, channel = c(channels.ind['CD11b']), tinypeak.removal = 0.01)
      maxDens <- density(getflowFrame(granulocytes.flowD)@exprs[, c(channels.ind['CD11b'])])
      if(CD11b.granulo.gate < (maxDens$x[which.max(maxDens$y)] - 0.1)){
        granulocytes.flowD <- flowDensity(granulocytes, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, T), gates = c(min(granulocytes@exprs[, channels.ind["CD5/Ly6G"]]), CD11b.granulo.gate) , ellip.gate = T, alpha = 0.85)
      }
      NOT.granulocytes <- notSubFrame(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), filter = granulocytes.flowD@filter )
      
      # Gating NOT(Granulocytes) to get Monocytes -----------------------------------------------------------------------------------------------------
     
      temp <- getflowFrame(NOT.granulocytes)@exprs[,  c(channels.ind["Ly6C"], channels.ind["CD11b"])]
      flowPeaks.Res2 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res2$peaks$cid, function(x){ length(which(flowPeaks.Res2$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res2$peaks$cid[which(cluster.size > 0.001)]
      
      mu.short <- flowPeaks.Res2$peaks$mu[sig.clusters,]
      monocytes.clusterid <- sig.clusters[which.max(1.5*mu.short[,2]^2 + mu.short[,1]^2)]
      monocytes.idx <- which(flowPeaks.Res2$peaks.cluster ==  monocytes.clusterid)
      monocytes <- getflowFrame(NOT.granulocytes)
      monocytes@exprs <- monocytes@exprs[monocytes.idx, ]
      
      Ly6C.monocyte.gate <- deGate(monocytes, channel = channels.ind["Ly6C"], upper = F, tinypeak.removal = 0.5, percentile = NA)
      CD11b.monocyte.gate <- deGate(monocytes, channel = channels.ind["CD11b"], upper = F, tinypeak.removal = 0.5, percentile = NA)
      monocytes.flowD <- flowDensity(monocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), position = c(T, T), gates = c(Ly6C.monocyte.gate, CD11b.monocyte.gate), ellip.gate = T, alpha = 0.9)
      
      NOT.monocytes <- notSubFrame(NOT.granulocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), filter = monocytes.flowD@filter)
      
      # Gating NOT(Monocytes) to get Eosinophils -------------------------------------------------------
      
      # # Jax does not always have SSC-H, so use SSC-A instead then.
      # if(is.na(scat.chans['SSC-H'])){
      #   print('SSC-H does not exist')
      #   scat.chans['SSC-H'] <- scat.chans['SSC-A']
      # }
      
      # cut off main peak otherwise might only get one cluster
      maxDens <- density(getflowFrame(NOT.monocytes)@exprs[,  c(channels.ind["CD11b"])])
      temp.flowD <- flowDensity(NOT.monocytes, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, NA), gates = c(maxDens$x[which.max(maxDens$y)] + 0.1, NA))
      
      temp <- getflowFrame(temp.flowD)@exprs[,  c(channels.ind["CD11b"], scat.chans['SSC-H'])]
      # temp <- getflowFrame(NOT.monocytes)@exprs[,  c(channels.ind["CD11b"], scat.chans['SSC-H'])]
      temp[,2] <- temp[,2]/100000 # scale SSC-H so flowPeaks works properly
      flowPeaks.Res3 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res3$peaks$cid, function(x){ length(which(flowPeaks.Res3$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res3$peaks$cid[which(cluster.size > 0.0001)]
      
      mu.short <- flowPeaks.Res3$peaks$mu[sig.clusters, ]
      potential.eosin.clusters <-sig.clusters[which((mu.short[,2] > (min(mu.short[,2]) + 0.5)) & (mu.short[,1] > (min(mu.short[,1]) + 0.5))  )]
      
      if(length(potential.eosin.clusters) > 1){
        #mu.short <- flowPeaks.Res3$peaks$mu[potential.eosin.clusters, ]
        #eosinophils.clusterid <- potential.eosin.clusters[which.max(3*mu.short[,2]^2 + mu.short[,1]^2)]
        potential.eosin.clusters <- potential.eosin.clusters[which(cluster.size[potential.eosin.clusters] > 0.3*max(cluster.size[potential.eosin.clusters]))]
        if(length(potential.eosin.clusters) > 1){
          mu.short <- flowPeaks.Res3$peaks$mu[potential.eosin.clusters, ]
          eosinophils.clusterid <- potential.eosin.clusters[which.max(2.5*mu.short[,2]^2 + mu.short[,1]^2)]
        }else{
          eosinophils.clusterid <- potential.eosin.clusters
        }
      }else{
        eosinophils.clusterid <- potential.eosin.clusters
      }
      eosinophils.idx <- which(flowPeaks.Res3$peaks.cluster ==  eosinophils.clusterid)
      #eosinophils <- getflowFrame(NOT.monocytes)
      eosinophils <- getflowFrame(temp.flowD)
      eosinophils@exprs <- eosinophils@exprs[eosinophils.idx, ]
      
      CD11b.eosino.gate <- min(eosinophils@exprs[, channels.ind["CD11b"]])
      eosinophils.flowD <- flowDensity(eosinophils, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, T), gates = c(CD11b.eosino.gate, min(eosinophils@exprs[, scat.chans['SSC-H']])), ellip.gate = T, alpha = 0.85)
      
      CD11b.eosino.gate <- deGate(eosinophils.flowD, channel = c(channels.ind["CD11b"]), twin.factor = 0.9)
      if(is.null(names(CD11b.eosino.gate))){
        eosinophils.flowD <- flowDensity(eosinophils, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, T), gates = c(CD11b.eosino.gate, min(eosinophils@exprs[, scat.chans['SSC-H']])), ellip.gate = T, alpha = 0.85)
      }
      
      NOT.eosinophils <- notSubFrame(NOT.monocytes, channels =  c(channels.ind["CD11b"], scat.chans['SSC-H']), filter = eosinophils.flowD@filter )
      
      # clear some memory up - I don't need these if I am not plotting them
      rm(flowPeaks.Res)
      rm(flowPeaks.Res2)
      rm(flowPeaks.Res3)
      
      # Gating NOT(Eosinophils) to get CD161+ and CD161- ------------------------------------------------
      
      # TCP, CIPHE, BCM all have edge events at the CD161+ side. Remove
      #maxCD161.value <- max(getflowFrame(NOT.eosinophils)@exprs[, c(channels.ind['CD161'])])
      #NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(0.99*maxCD161.value, NA))
      
      CD19lo.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD19']), use.percentile = T, percentile = 0.005) - 0.5
      CD161hi.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD161']), use.percentile = T, percentile = 0.995) + 0.5
      NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, T), gates = c(CD161hi.gate, CD19lo.gate))
      
      # For Jax the CD161 and CD317 channels are shared. I need to remove CD317+ before gating CD161+
      if(j == 'Jax'){
        CD161CD317.gate <- max(deGate(NOT.eosinophils, channel = c(channels.ind['CD161']), all.cuts = T, tinypeak.removal = 0.001))
        NOT.eosinophils.preCD317gate <- NOT.eosinophils
        NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(CD161CD317.gate, NA))
      }
      
      
      # maxCD161.value <- deGate(NOT.eosinophils, channel = c(channels.ind['CD161']), use.upper = T, upper = T) + 0.4
      # NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(maxCD161.value, NA))
      # 
      CD19.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD19']))
      
      temp <- getflowFrame(NOT.eosinophils)@exprs[,  c(channels.ind["CD161"], channels.ind["CD19"])]
      flowPeaks.Res4 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res4$peaks$cid, function(x){ length(which(flowPeaks.Res4$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res4$peaks$cid[which(cluster.size > 0.001)]
      
      mu.short <- flowPeaks.Res4$peaks$mu[sig.clusters,]
      CD19pos.clusters <- flowPeaks.Res4$peaks$cid[which(mu.short[,2] > (max(mu.short[,2]) - 0.2))]
      CD19neg.idx <- which(!(flowPeaks.Res4$peaks.cluster %in%  CD19pos.clusters))
      CD19neg <- getflowFrame(NOT.eosinophils)
      CD19neg@exprs <- CD19neg@exprs[CD19neg.idx, ]
      
      # use the CD19- population to get a CD161 gate
      # CD19neg.flowD <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(NA, F), gates = c(NA, CD19.gate))
      #deGate(CD19neg, channel = c(channels.ind["CD161"]))
      
      CD161.gate <- deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.05) 
      # FOR BCM in one case (Panel2_M_232049_015.labelled.fcs) the above does not work, but the following does. I'd need to test if I can use this in general 
      # max(deGate(CD19neg, channel = c(channels.ind["CD161"]), all.cuts = T) 
      CD161.gate2 <- deGate(CD19neg, channel = c(channels.ind["CD161"]), after.peak = T) 
      if(is.null(names(CD161.gate2))){
        CD161.gate <- min(CD161.gate, CD161.gate2)
      }
      
      # This is Jax-specific
      if(j == 'Jax'){
        CD161.gate <- deGate(CD19neg, channel = c(channels.ind["CD161"]), after.peak = T) 
      }
      
      
      
      # maxDens <- density(getflowFrame(CD19neg.flowD)@exprs[, c(channels.ind['CD161'])])
      # peak.idx <- which.max(maxDens$y)
      # yval <- maxDens$y[which.min(abs(maxDens$x - CD161.gate))]
      # CD161.gate <- maxDens$x[min(which(maxDens$y[peak.idx:length(maxDens$y)] < (yval + 0.01))) + peak.idx - 1]
      # 
      # # Try to improve on the CD19 gate
      # temp.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, NA), gates = c(CD161.gate, NA))
      # CD19.gate <- deGate(temp.flowD, channel = c(channels.ind['CD19']), use.upper = T, upper = T, alpha = 0.05, tinypeak.removal = 0.01)
      
      # CD161pos.flowD <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      # NOT.CD161pos <- notSubFrame(NOT.eosinophils, channels =  c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      
      CD161pos.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      NOT.CD161pos <- notSubFrame(NOT.eosinophils, channels =  c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      
      rm(flowPeaks.Res4)
      
      # Gating CD161+ to get NK/NKT -----------------------------------------------------------------------------------
      
      CD5Ly6G.gate <- deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), tinypeak.removal = 0.01)
      maxDens <- density(getflowFrame(CD161pos.flowD)@exprs[, c(channels.ind["CD5/Ly6G"])])
      #NK.CD5.peak <- maxDens$x[which.max(maxDens$y[1:maxDens$x[which.min(abs(maxDens$x - 2))]])]
      NK.CD5.peak <- findpeaks(maxDens$y)
      NK.CD5.peak <- NK.CD5.peak[which(NK.CD5.peak[,1] > 0.2*max(NK.CD5.peak[,1])),]
      if(length(NK.CD5.peak) > 4){
        NK.CD5.peak <- maxDens$x[min(NK.CD5.peak[,2])] # pick left-most peak
      }else{
        NK.CD5.peak <- maxDens$x[NK.CD5.peak[2]]
      }
      if(CD5Ly6G.gate < NK.CD5.peak){
        CD5Ly6G.gate <- deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), all.cuts = T, tinypeak.removal = 0.01)
        CD5Ly6G.gate <- CD5Ly6G.gate[which(CD5Ly6G.gate > maxDens$x[which.max(maxDens$y)])]
        if(length(CD5Ly6G.gate) < 1){
          CD5Ly6G.gate <- deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), use.upper = T, upper = T)
        }else{
          CD5Ly6G.gate <- max(CD5Ly6G.gate)
        }
      } 
      if(!is.null(names(CD5Ly6G.gate))){
        CD5Ly6G.gate <-  deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), use.upper = T, upper = T)
      } 
      
      NKcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(F, NA), gates = c(CD5Ly6G.gate, NA))
      CD161.NK.gate <- deGate(NKcells.flowD, channel = c(channels.ind["CD161"]))
      maxDens <- density(getflowFrame(CD161pos.flowD)@exprs[, c(channels.ind['CD161'])])
      NK.CD161.peak <- maxDens$x[which.max(maxDens$y)]
      if(is.null(names(CD161.NK.gate))){
        if(CD161.NK.gate < NK.CD161.peak){
          NKcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(F, T), gates = c(CD5Ly6G.gate, CD161.NK.gate))
        }
      }
      NKTcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(T, NA), gates = c(CD5Ly6G.gate, NA))
      
      # Try to use NK population to get the CD11b and Ly6C gates (tried using whole CD161+ population and this seems better)
      CD11b.gate <- deGate(getflowFrame(NKcells.flowD), channel = c(channels.ind["CD11b"]))
      Ly6C.gate <- min(deGate(getflowFrame(NKcells.flowD), channel = c(channels.ind["Ly6C"]), all.cuts = T))
      
      # Gating NOT(CD161+) to get T cells ect. --------------------------------------- ------------------------------
      
      temp <- getflowFrame(NOT.CD161pos)@exprs[,  c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"])]
      flowPeaks.Res5 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res5$peaks$cid, function(x){ length(which(flowPeaks.Res5$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res5$peaks$cid[which(cluster.size > 0.01)]
      
      mu.short <- flowPeaks.Res5$peaks$mu[sig.clusters,]
      
      Tcell.clusterid <- sig.clusters[which.max(2*mu.short[, 2]^2 + (4.5 - mu.short[, 1])^2)]
      Tcell.clusterid <- union(Tcell.clusterid, sig.clusters[which((mu.short[, 2] > (min(mu.short[, 2] + 0.2)) & (mu.short[, 1] < mu.short[Tcell.clusterid, 1] + 0.2)))])
      
      # remove low density regions of the clusters
      fpc <- assign.flowPeaks(flowPeaks.Res5, flowPeaks.Res5$x, tol = 0.015, fc = 0)
      
      #Tcell.idx <- which(flowPeaks.Res5$peaks.cluster %in%  Tcell.clusterid)
      Tcell.idx <- which(fpc %in% Tcell.clusterid)
      Tcells <- getflowFrame(NOT.CD161pos)
      Tcells@exprs <- Tcells@exprs[Tcell.idx, ]
      
      # temp <- getflowFrame(NOT.CD161pos)
      # temp.idx <- which((temp@exprs[, c(channels.ind['MHCII'])] < min(mu.short[Tcell.clusterid, 1])) & (temp@exprs[, c(channels.ind['CD5/Ly6G'])] < min(mu.short[, 2])))
      # 
      maxDens <- density(Tcells@exprs[, c(channels.ind['CD5/Ly6G'])])
      Tcell.peak <- maxDens$x[which.max(maxDens$y)]
      
      # # If The CD5-MHCII- population is included in the T cell cluster, then it is easier to find if I don't remove the low density regions.
      Tcell2.idx <- which(flowPeaks.Res5$peaks.cluster %in%  Tcell.clusterid)
      Tcells2 <- getflowFrame(NOT.CD161pos)
      Tcells2@exprs <- Tcells2@exprs[Tcell2.idx, ]
      
      CD5.gate <- deGate(Tcells2, channel = c(channels.ind["CD5/Ly6G"]), tinypeak.removal = 0.05)
      
      if(CD5.gate < (Tcell.peak - 0.5)){
        Tcells <- getflowFrame(NOT.CD161pos)
        Tcell.idx <- which((flowPeaks.Res5$peaks.cluster %in%  Tcell.clusterid) & (Tcells@exprs[, c(channels.ind["CD5/Ly6G"])] > CD5.gate))
        Tcell.idx <- which((fpc %in%  Tcell.clusterid) & (Tcells@exprs[, c(channels.ind["CD5/Ly6G"])] > CD5.gate))
        Tcells@exprs <- Tcells@exprs[Tcell.idx, ]
      }
      
      #MHCII.Tcell.gate <- deGate(Tcells, channel = c(channels.ind['MHCII']), use.upper = T, upper = T, alpha = 0.05)
      #Tcells.flowD <- flowDensity(Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(F, NA), gates = c(MHCII.Tcell.gate, NA))
    
      
      Ly6Cdensity.Tcells <- density(Tcells@exprs[, c(channels.ind['Ly6C'])])

      NOT.Tcells <- getflowFrame(NOT.CD161pos)
      temp2 <- Tcells@exprs[which.max(Tcells@exprs[, channels.ind['MHCII']]), c(channels.ind['MHCII'], channels.ind['CD5/Ly6G'])]
      temp <- Tcells@exprs[which.min(Tcells@exprs[, channels.ind['CD5/Ly6G']]), c(channels.ind['MHCII'], channels.ind['CD5/Ly6G'])]
      
      idx.to.remove <- union(which((NOT.Tcells@exprs[, c(channels.ind['CD5/Ly6G'])] > temp2[2]) & (NOT.Tcells@exprs[, c(channels.ind['MHCII'])] < temp2[1])),
                             which((NOT.Tcells@exprs[, c(channels.ind['CD5/Ly6G'])] > temp[2]) & (NOT.Tcells@exprs[, c(channels.ind['MHCII'])] < temp[1])))
      idx.to.remove <- union(idx.to.remove, Tcell.idx)
      #idx.to.remove <- union(idx.to.remove, Tcell.idx[which(Tcells@exprs[, c(channels.ind['MHCII'])] < MHCII.Tcell.gate)])
      NOT.Tcells@exprs <- NOT.Tcells@exprs[-idx.to.remove, ]
      #Tcells <- getflowFrame(Tcells.flowD)
      
      MHCII.gate <- max(deGate(NOT.Tcells, channel = c(channels.ind['MHCII']), all.cuts = T, tinypeak.removal = 0.005, upper = F, percentile = NA))
      MHCII.gate.alt <- deGate(NOT.Tcells, channel = c(channels.ind['MHCII']), use.upper = T, upper = F, tinypeak.removal = 0.5, alpha = 0.05)
      if(!is.null(names(MHCII.gate))){
        MHCII.gate <- MHCII.gate.alt
      }else if(abs(MHCII.gate - MHCII.gate.alt) > 0.15){
        MHCI.gate <- MHCII.gate.alt
      }

      MHCIIpos.flowD <- flowDensity(NOT.Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(T, NA), gates = c(MHCII.gate, NA))

      # CD5.MHCIIpos.gate <- deGate(MHCIIpos.flowD, channel = channels.ind['CD5/Ly6G'], use.upper = T, upper = T)
      # MHCIIpos.flowD <- flowDensity(NOT.Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(T, F), gates = c(MHCII.gate, CD5.MHCIIpos.gate ))
      #

      # if(length(intersect(temp.idx, Tcell.idx))/length(temp.idx) > 0.1){
      #     # use flowDensity to remove CD5-MHCII- cells
      #     CD5.gate <- deGate(Tcells, channel = c( channels.ind["CD5/Ly6G"]), upper = F, percentile = NA)
      #     Tcells <- getflowFrame(flowDensity(Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(NA, T), gates = c(NA, CD5.gate)))
      # }


      # Gating MHCII+ to get B cells, RPmacs and cDcs --------------------------------------------------

      # get RP macs, if F4/80 exists
      if(!is.na(channels.ind['F4/80'])){

        F4.gate <- deGate(MHCIIpos.flowD@flow.frame, channel = c(channels.ind['F4/80']), use.upper = T, upper = T, tinypeak.removal = 0.9) - 0.1  # I don't really like this. Change later
        #ff <- flowDensity(MHCIIpos.flowD, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(F, NA), gates = c(F4.gate, NA))
        RPmac.flowD <- flowDensity(MHCIIpos.flowD, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(T, NA), gates = c(F4.gate, NA))
        RPmac.idx <- RPmac.flowD@index

        RPmac <- getflowFrame(RPmac.flowD)

        # May need to remove some additional events
        temp <- RPmac@exprs[, c(channels.ind["F4/80"], channels.ind["MHCII"])]
        flowPeaks.Res6 <- flowPeaks(temp)

        cluster.size <- sapply(flowPeaks.Res6$peaks$cid, function(x){ length(which(flowPeaks.Res6$peaks.cluster == x))})
        cluster.size <- cluster.size/nrow(temp)
        sig.clusters <- flowPeaks.Res6$peaks$cid[which(cluster.size > 0.01)]

        ff <- MHCIIpos.flowD@flow.frame

        mu.short <- flowPeaks.Res6$peaks$mu[sig.clusters,]
        if(length(sig.clusters) > 1){
          clusterid.toremove <- sig.clusters[which.max(mu.short[, 2]^2 + (4.5 - mu.short[, 1])^2)]
          toremove.idx <- which(flowPeaks.Res6$peaks.cluster %in%  clusterid.toremove)
          RPmac@exprs <- RPmac@exprs[-toremove.idx, ]
          RPmac.idx <- RPmac.idx[-toremove.idx]
          ff@exprs <- ff@exprs[-RPmac.idx,]
        }

      }else{
        ff <- MHCIIpos.flowD@flow.frame
      }

      temp <- na.omit(ff@exprs[, c(channels.ind["CD19"], channels.ind["CD11c"])])
      flowPeaks.Res7 <- flowPeaks(temp)

      cluster.size <- sapply(flowPeaks.Res7$peaks$cid, function(x){ length(which(flowPeaks.Res7$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res7$peaks$cid[which(cluster.size > 0.001)]

      mu.short <- flowPeaks.Res7$peaks$mu[sig.clusters, ]

      # cDC.clusterid <- sig.clusters[which.min(abs((mu.short[,1] - min(mu.short[,1]))^2 + (mu.short[,2] - max(mu.short[,2]))^2))]
      # cDCs.idx <- which(flowPeaks.Res7$peaks.cluster %in%  cDC.clusterid)
      #
      # cDCs <- ff
      # cDCs@exprs <- na.omit(cDCs@exprs)
      # cDCs@exprs <- cDCs@exprs[cDCs.idx, ]
      #
      # CD11c.cDC.gate <- deGate(cDCs, channel = c(channels.ind["CD11c"]), use.upper = T, upper = F)
      # cDCs.flowD <- flowDensity(cDCs, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, T), gates = c(max(cDCs@exprs[,c(channels.ind['CD19'])]), CD11c.cDC.gate), ellip.gate = T, scale = 0.95)
      #
      Bcells.clusterid <- sig.clusters[which(mu.short[, 1] > (max(mu.short[, 1]) - 0.2))]
      Bcells.idx <- which(flowPeaks.Res7$peaks.cluster %in%  Bcells.clusterid)

      Bcells <- ff
      Bcells@exprs <- na.omit(Bcells@exprs)
      Bcells@exprs <- Bcells@exprs[Bcells.idx, ]

      CD19.Bcell.gate <- deGate(Bcells, channel = c(channels.ind['CD19']), tinypeak.removal = 0.01, upper = F, percentile = NA)
      Bcells.flowD <- flowDensity(Bcells, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(T, NA), gates = c(CD19.Bcell.gate, NA))

      # cannot use the CD5Ly6G.gate found for NK/NKT gating for BCM (tho it works really well for TCP./CIPHE)
      CD5Ly6G.Bcell.gate <- deGate(getflowFrame(Bcells.flowD), channel = c(channels.ind['CD5/Ly6G']), use.upper = T, upper = T, tinypeak.removal = 0.5)
      B2Bcells.flowD <- flowDensity(Bcells.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD19"]), position = c(T, NA), gates = c(CD5Ly6G.Bcell.gate, NA))


      temp <- getflowFrame(B2Bcells.flowD)@exprs[, c(channels.ind["CD21/CD35"], channels.ind["CD23"])]
      flowPeaks.Res8 <- flowPeaks(temp)

      
      # Plots -----------------------------------------------------------------------------------------
      
      png(file = paste0(scatterplot.dir, '/', fname, '.png'), width = 2000*6/5, height = 2000)
      par(mfrow = c(5,6), mar = (c(5, 5, 4, 2) + 0.1))
      
      # Plot flowCut data
      plotDens(fT, channels = c(channels.ind["Time"], channels.ind["CD19"]), main = " All Events", cex.lab = 2, cex.axis = 2, cex.main=2)
      points(fT@exprs[f.Clean$ind, c(channels.ind["Time"], channels.ind["CD19"])], pch ='.')
      
      
      plotDens(f, c(scat.chans['FSC-A'], channels.ind["Live"]), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(live.flowD@filter)
      
      plotDens(live.flowD@flow.frame, c('FSC-A','SSC-A'), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
      
      # Plot FSC singlet gating
      plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "Size", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(FSCsinglets.flowD@filter)
      
      # Plot SSC singlet gating
      if(!is.na(scat.chans['SSC-W'])){ 
        
        col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = primary.colors, nbin = 1000)
        plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
        lines(singlets.flowD@filter, col = 'blue')
        
        col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = matlab.like2, nbin = 1000)
        plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
        lines(singlets.flowD@filter)
        
      }else{
        plot(1, type="n", axes=F, xlab="", ylab="")
        plot(1, type="n", axes=F, xlab="", ylab="")
      }
      
      # Panel 2 gating
      
      # # Gating singlets to get Granulocytes
      plotDens(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), main = "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(granulocytes.flowD@filter)
      text(1, max(granulocytes.flowD@filter[, 2]), labels =  strcat('Granulocytes\n', toString(signif(granulocytes.flowD@cell.count/nrow(singlets)*100, 2))), cex = 2)
      
  
      # Gating NOT(Granulocytes) to get Monocytes
      plotDens(NOT.granulocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), main = "NOT(granulocytes)", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(monocytes.flowD@filter)
      text(2, max(monocytes.flowD@filter[, 2]), labels =  strcat('Monocytes\n', toString(signif(monocytes.flowD@cell.count/NOT.granulocytes@cell.count*100, 3))), cex = 2)
      

      # Gating NOT(Monocytes) to get Eosinophils
      plotDens(NOT.monocytes, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), main = "NOT(monocytes)", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(eosinophils.flowD@filter)
      text(min(eosinophils.flowD@filter[, 1]) - 1, max(eosinophils.flowD@filter[, 2]) - 5000, labels =  strcat('Eosinophils\n', toString(signif(eosinophils.flowD@cell.count/NOT.monocytes@cell.count*100, 3))), cex = 2)
      
      #plot(flowPeaks.Res3)
      
      plotDens(getflowFrame(eosinophils.flowD), channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), main = "Eosinophils",
               xlim = c(0,4.5), ylim = c(0,4.5), cex.lab = 2, cex.axis = 2, cex.main=2)
      
      # For Jax the CD161 and CD317 channels are shared. I need to remove CD317+ before gating CD161+
      if(j == 'Jax'){
        plotDens(NOT.eosinophils.preCD317gate, channels = c(channels.ind["CD161"], channels.ind["CD19"]), main = "NOT(eosinophils)", cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(CD161pos.flowD@filter)
        text(min(CD161pos.flowD@filter[, 1]) + 1, max(CD161pos.flowD@filter[, 2]) + 0.2, labels =  strcat('CD161+\n', toString(signif(CD161pos.flowD@cell.count/NOT.eosinophils.preCD317gate@cell.count*100, 3))), cex = 2)
      }else{
        plotDens(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), main = "NOT(eosinophils)", cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(CD161pos.flowD@filter)
        text(min(CD161pos.flowD@filter[, 1]) + 1, max(CD161pos.flowD@filter[, 2]) + 0.2, labels =  strcat('CD161+\n', toString(signif(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count*100, 3))), cex = 2)
      }
      
      # Plot NK and NKT cells
      plotDens(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), main = "CD161+", cex.lab = 2, cex.axis = 2, cex.main=2)
      abline(v = CD5Ly6G.gate)
      lines(NKcells.flowD@filter)
      text(CD5Ly6G.gate - 0.5, NK.CD161.peak + 0.5, labels =  strcat('NK cells\n', toString(signif(NKcells.flowD@cell.count/CD161pos.flowD@cell.count*100, 3))), cex = 2)
      text(CD5Ly6G.gate + 1, NK.CD161.peak + 0.5, labels =  strcat('NKT cells\n', toString(signif(NKTcells.flowD@cell.count/CD161pos.flowD@cell.count*100, 3))), cex = 2)
      
      
      plotDens(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), main = "NK cells", cex.lab = 2, cex.axis = 2, cex.main=2)
      abline(v = CD11b.gate)
      abline(h = Ly6C.gate)
      plotDens(NKTcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), main = "NKT cells", cex.lab = 2, cex.axis = 2, cex.main=2)
      abline(v = CD11b.gate)
      abline(h = Ly6C.gate)
      
      plotDens(getflowFrame(NOT.CD161pos), channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), main = "NOT(CD161+)", cex.lab = 2, cex.axis = 2, cex.main=2)
      plot(flowPeaks.Res5)
      
      plotDens(Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), main = "T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
      
      plot(Ly6Cdensity.Tcells, main = "T cells",  xlab = paste0("<", pData(parameters(f))$desc[channels.ind['Ly6C']], ">:", pData(parameters(f))$name[channels.ind['Ly6C']]),cex.lab = 2, cex.axis = 2, cex.main=2)
      polygon(Ly6Cdensity.Tcells, col = 'gray', border = 'black')
      abline(v = Ly6C.gate, lty = 2)
      
      plotDens(NOT.Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), main = "NOT(T cells)", cex.lab = 2, cex.axis = 2, cex.main=2)
      abline(v = MHCII.gate)
      lines(MHCIIpos.flowD@filter)
      #abline(h = CD5.MHCIIpos.gate, lty = 2)
      
      if(!is.na(channels.ind['F4/80'])){
        plotDens(getflowFrame(MHCIIpos.flowD), channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), main = "MHCII+", cex.lab = 2, cex.axis = 2, cex.main=2)
        abline(v = F4.gate)
        
        #plot(flowPeaks.Res6)
        plotDens(RPmac, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), main = "RP mac", cex.lab = 2, cex.axis = 2, cex.main=2)

        plotDens(RPmac, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "RP mac", cex.lab = 2, cex.axis = 2, cex.main=2)
        abline(v = F4.gate)
        
        plotDens(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+ - RP mac", cex.lab = 2, cex.axis = 2, cex.main=2)
        data.new <- exprs(getflowFrame(MHCIIpos.flowD))[, c(channels.ind["CD19"], channels.ind["CD11c"])]
        z <- kde2d(data.new[, 1], data.new[, 2])
        contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
        
        # lines(cDCs.flowD@filter)
        # text(max(cDCs.flowD@filter[,1]) -0.5, max(cDCs.flowD@filter[,2]) + 0.25 , labels =  strcat('cDCs \n', toString(signif(cDCs.flowD@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)
        # 
        # 
      }else{
        
        plotDens(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+", cex.lab = 2, cex.axis = 2, cex.main=2)
        data.new <- exprs(getflowFrame(MHCIIpos.flowD))[, c(channels.ind["CD19"], channels.ind["CD11c"])]
        z <- kde2d(data.new[, 1], data.new[, 2])
        contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
        
        # lines(cDCs.flowD@filter)
        # text(max(cDCs.flowD@filter[,1]) -0.5, max(cDCs.flowD@filter[,2]) + 0.25 , labels =  strcat('cDCs \n', toString(signif(cDCs.flowD@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)
        # 
        
      }
      
      
      plot(flowPeaks.Res7)
      
      plotDens(Bcells, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "B cell cluster", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(Bcells.flowD@filter)
      
      # plotDens(cDCs, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "cDC cluster", cex.lab = 2, cex.axis = 2, cex.main=2)
      # lines(cDCs.flowD@filter)
      
      plotDens(getflowFrame(Bcells.flowD), channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD19"]), main = "Bcells", cex.lab = 2, cex.axis = 2, cex.main=2)
      abline(v = CD5Ly6G.Bcell.gate, lty = 2)
      
      plotDens(getflowFrame(B2Bcells.flowD), channels = c(channels.ind["CD21/CD35"], channels.ind["CD23"]), main = "B2Bcells", cex.lab = 2, cex.axis = 2, cex.main=2)
      data.new <- exprs(getflowFrame(B2Bcells.flowD))[, c(channels.ind["CD21/CD35"], channels.ind["CD23"])]
      z <- kde2d(data.new[, 1], data.new[, 2])
      contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
      
      plot(flowPeaks.Res8)
      
      dev.off()
      
    })
  }
  
}
