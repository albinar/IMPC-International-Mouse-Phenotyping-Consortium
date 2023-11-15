# Testing FMO controls for gating cDC and RPMacs
# This script gets the flow frame just prior to Bcell/RPmac and cDC subdivision
# Saves the resultant flow frame so I can use it later to gate the RPmacs and cDCs

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

source("/data/IMPC/scripts/P2gating_additionalfunctions.R")
source("/data/IMPC/scripts/flowCut_20170331.R")


TCPdata.dir <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/TCP/TCP_170209/UPLOAD_DUMP"
allFilesP2 <- dir(TCPdata.dir, recursive = T, pattern="PANEL_B", full.names = T)
CD11cFMOs <- allFilesP2[grep("CD11C",allFilesP2)]
centre <- 'TCP'

BCMdata.dir <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/BCM/BCM"
allFilesP2 <- dir(BCMdata.dir, recursive = T, pattern="Panel2", full.names = T)
CD11cFMOs <- allFilesP2[grep("CD11c",allFilesP2)]
centre <- 'BCM'

foreach(i = 1:length(CD11cFMOs)) %dopar% {
  
  try({
    
    f <- read.FCS(CD11cFMOs[i])
    fname <- basename(CD11cFMOs[i])
    fname.date <-  basename(dirname(CD11cFMOs[i]))
    fname.full <- CD11cFMOs[i]
    
    
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
    
    # there is atleast one file where all the pData(parameters(f))$desc were <NA>, so channels.ind contains only Time and Live/Dead (BV-510A)
    # "FMO_PANEL_B_CD11C_039.fcs"
    
    # For now I'll skip this file
    if(length(channels.ind) > 2){
      
      scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
      names(scat.chans) <- colnames(f)[scat.chans]
      
      # Remove scatter margins and compensate ---------------------------------------------
      # Removing margin events in Scatter channels
      f <- removeMargins(f, chans = scat.chans, verbose = F)
      #Removing negative values in scatter channels
      f <- removeMargins(f, chans = scat.chans, debris = T, neg = T, verbose = F)
      
      f <- compensateIMPC(f, fname.full, centre = centre, panel.no = 2)
      
      # Transformation
      # if(centre == 'Jax'){ # just haven't set up a global transformation matrix for BCM yet
      #   no.transform <- union(scat.chans, channels.ind['Time'])
      #   lgl <- estimateLogicle(f, channels = colnames(f)[-no.transform])
      # }else{
      load(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/', centre, "/Panel2/Results/lgl.Rdata"))
      # }
      
      fT <- transform(f, lgl)
      channels.to.clean <- setdiff(channels.ind, channels.ind['Time'])
      f.Clean <- flowCut(fT, segment = floor(nrow(f)*5/3000), CleanChan.loc = channels.to.clean, FileID = fname, directory = paste0('/data/IMPC/ScatterPlots/', centre, '/flowCut/'), Plot = 'Failed Only')
      passed.flowCut <- f.Clean$data['Has the file passed',1]
      
      # For CIPHE The F4/80 channel looks best if I use the parameters from the FCS file for transformation so I'll do that
      if(centre == 'CIPHE'){
        f.temp <- logiclTransformCiphe(f, markers.transform = colnames(f)[channels.ind['F4/80']])
      }
      f <- transform(f, lgl)
      if(centre == 'CIPHE'){
        f@exprs[,c(channels.ind['F4/80'])] <- f.temp@exprs[,c(channels.ind['F4/80'])]
        rm(f.temp)
      }
      if(length(f.Clean$ind) > 0){
        f@exprs <- f@exprs[-f.Clean$ind, ]
      }
      
      # Quality Gate (live/dead, singlets gating) ---------------------------------------------------
      
      results <- qualityGate(f, scat.chans, channels.ind, centre, panel.no = 2)
      live.flowD <- results$live
      FSCsinglets.flowD <- results$FSCsinglets
      singlets.flowD <- results$singlets
      scat.chans <- results$scat.chans
      
      # Panel 2 Gating -------------------------------------------------------------------------------
      
      singlets <- getflowFrame(singlets.flowD)
      
      # Gating singlets to get Granulocytes ------------------------------------------------------------
      
      singlets0 <- singlets # Just to check the CIPHE/BCM filtering of edge events
      
      if((centre == 'CIPHE') | (centre == 'BCM')){ # For the CIPHE there are some days where there are edge events along the top edge of the CD11b channel
        #max.CD11b.value <- max(singlets@exprs[, c(channels.ind["CD11b"])])
        max.CD11b.value <- quantile(singlets@exprs[, c(channels.ind["CD11b"])], c(.999)) + 0.25
        singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["CD11b"])] < max.CD11b.value), ]
      }
      
      # This helps remove noisy events below the CD5peak - which is problematic on some days for CIPHE
      CD5.gate <- deGate(singlets, channel = channels.ind["CD5/Ly6G"], use.upper = T, upper = F)
      CD5pos.flowD <- flowDensity(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, NA), gates = c(CD5.gate, NA))
      
      temp <- getflowFrame(CD5pos.flowD)@exprs[, c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"])]
      
      flowPeaks.Res <- flowPeaks(temp)
      cluster.size <- sapply(flowPeaks.Res$peaks$cid, function(x){ length(which(flowPeaks.Res$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res$peaks$cid[which(cluster.size > 0.001)]
      
      mu.short <- flowPeaks.Res$peaks$mu[sig.clusters,]
      granulo.clusterid <- sig.clusters[which((mu.short[, 2] >  (min(mu.short[,2]) + 0.5*(max(mu.short[,2]) - min(mu.short[,2])))) & (mu.short[, 1] >  (min(mu.short[, 1]) + 0.5)))]
      if(length(granulo.clusterid) > 1){
        granulo.clusterid <- granulo.clusterid[which(cluster.size[granulo.clusterid] > 0.2*max(cluster.size[granulo.clusterid]))]
        granulo.clusterid <- granulo.clusterid[which.max(4.5*flowPeaks.Res$peaks$mu[granulo.clusterid,2]^2 + flowPeaks.Res$peaks$mu[granulo.clusterid,1]^2)]
      }
      
      granulo.idx <- which(flowPeaks.Res$peaks.cluster ==  granulo.clusterid)
      granulocytes <- getflowFrame(CD5pos.flowD)
      #granulocytes <- singlets
      granulocytes@exprs <- granulocytes@exprs[granulo.idx, ]
      
      CD11b.granulo.gate <- min(granulocytes@exprs[, channels.ind["CD11b"]])
      granulocytes.flowD <- flowDensity(granulocytes, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, T), gates = c(min(granulocytes@exprs[, channels.ind["CD5/Ly6G"]]), CD11b.granulo.gate) , ellip.gate = T, alpha = 0.85)
      CD11b.granulo.gate <- deGate(granulocytes.flowD, channel = c(channels.ind['CD11b']), tinypeak.removal = 0.01)
      maxDens <- density(getflowFrame(granulocytes.flowD)@exprs[, c(channels.ind['CD11b'])])
      if(CD11b.granulo.gate < (maxDens$x[which.max(maxDens$y)] - 0.1)){
        granulocytes.flowD <- flowDensity(granulocytes, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), position = c(T, T), gates = c(min(granulocytes@exprs[, channels.ind["CD5/Ly6G"]]), CD11b.granulo.gate) , ellip.gate = T, alpha = 0.85)
      }
      NOT.granulocytes <- notSubFrame(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), filter = granulocytes.flowD@filter )
      
      # clear some memory up - I don't need these if I am not plotting them
      rm(flowPeaks.Res)
      
      # Gating NOT(Granulocytes) to get Monocytes -----------------------------------------------------------------------------------------------------
      
      maxDens <- density(getflowFrame(NOT.granulocytes)@exprs[,  c(channels.ind["Ly6C"])])
      temp.flowD <- flowDensity(getflowFrame(NOT.granulocytes), channels = c(channels.ind["Ly6C"], channels.ind['CD11b']), position = c(T, NA), gates = c(maxDens$x[which.max(maxDens$y)] + 0.1, NA))
      # 
      #temp <- getflowFrame(NOT.granulocytes)@exprs[,  c(channels.ind["Ly6C"], channels.ind["CD11b"])]
      temp <- getflowFrame(temp.flowD)@exprs[,  c(channels.ind["Ly6C"], channels.ind["CD11b"])]
      flowPeaks.Res2 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res2$peaks$cid, function(x){ length(which(flowPeaks.Res2$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res2$peaks$cid[which(cluster.size > 0.001)]
      
      mu.short <- flowPeaks.Res2$peaks$mu[sig.clusters,]
      monocyte.clusterid <-sig.clusters[which((mu.short[, 2] > (min(mu.short[,2]) + 0.5*(max(mu.short[, 2]) -  min(mu.short[,2])))) & (mu.short[, 1] > (max(mu.short[, 1]) - 1))  )]
      if(length(monocyte.clusterid) > 1){
        monocyte.clusterid <- monocyte.clusterid[which(cluster.size[monocyte.clusterid] > 0.3*max(cluster.size[monocyte.clusterid]))]
        monocyte.clusterid <- monocyte.clusterid[which.max(1.1*flowPeaks.Res2$peaks$mu[monocyte.clusterid, 2]^2 + flowPeaks.Res2$peaks$mu[monocyte.clusterid, 1]^2)]
      }
      
      monocytes.idx <- which(flowPeaks.Res2$peaks.cluster ==  monocyte.clusterid)
      monocytes <- getflowFrame(temp.flowD)
      monocytes@exprs <- monocytes@exprs[monocytes.idx, ]
      
      Ly6C.monocyte.gate <- deGate(monocytes, channel = channels.ind["Ly6C"], upper = F, tinypeak.removal = 0.5, percentile = NA)
      Ly6C.monocyte.alt.gate <- deGate(monocytes, channel = channels.ind["Ly6C"], tinypeak.removal = 0.05)
      if(is.null(names(Ly6C.monocyte.alt.gate))){
        Ly6C.monocyte.gate <- max(Ly6C.monocyte.alt.gate, Ly6C.monocyte.gate)
      }
      maxDens <- density(monocytes@exprs[, c(channels.ind['Ly6C'])])
      Ly6C.monocyte.alt.gate <- 2*maxDens$x[which.max(maxDens$y)] - deGate(monocytes, channel = channels.ind["Ly6C"], use.upper = T, upper = T) 
      if((maxDens$x[which.max(maxDens$y)] - Ly6C.monocyte.gate) > 2*(maxDens$x[which.max(maxDens$y)] - Ly6C.monocyte.alt.gate)){
        Ly6C.monocyte.gate <- Ly6C.monocyte.alt.gate
      }
      #CD11b.monocyte.gate <- deGate(monocytes, channel = channels.ind["CD11b"], upper = F, tinypeak.removal = 0.5, percentile = NA)
      CD11b.monocyte.gate <- deGate(monocytes, channel = channels.ind["CD11b"], use.upper = T, upper = F, percentile = NA)
      monocytes.flowD <- flowDensity(monocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), position = c(T, T), gates = c(Ly6C.monocyte.gate, CD11b.monocyte.gate), ellip.gate = T, alpha = 0.9)
      
      NOT.monocytes <- notSubFrame(NOT.granulocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), filter = monocytes.flowD@filter)
      
      rm(flowPeaks.Res2)
      
      # Gating NOT(Monocytes) to get Eosinophils -------------------------------------------------------
      
      # # Jax does not always have SSC-H, so use SSC-A instead then.
      # if(is.na(scat.chans['SSC-H'])){
      #   print('SSC-H does not exist')
      #   scat.chans['SSC-H'] <- scat.chans['SSC-A']
      # }
      
      # cut off main peak otherwise might only get one cluster
      maxDens <- density(getflowFrame(NOT.monocytes)@exprs[,  c(channels.ind["CD11b"])])
      temp.flowD <- flowDensity(NOT.monocytes, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, NA), gates = c(maxDens$x[which.max(maxDens$y)] + 0.1, NA))
      
      # This is an approximate location for where I want the eosinophil cluster to be centered (actually 0.2 below this value)
      CD11b.temp.gate <- deGate(NOT.monocytes, channel = c(channels.ind["CD11b"]), upper = T, tinypeak.removal = 0.9)
      temp2.flowD <- flowDensity(NOT.monocytes, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, NA), gates = c(CD11b.temp.gate, NA))
      CD11b.upper.gate <- deGate(temp2.flowD, channel = c(channels.ind["CD11b"]), use.upper = T, upper = T) 
      
      temp <- getflowFrame(temp.flowD)@exprs[,  c(channels.ind["CD11b"], scat.chans['SSC-H'])]
      # temp <- getflowFrame(NOT.monocytes)@exprs[,  c(channels.ind["CD11b"], scat.chans['SSC-H'])]
      temp[,2] <- temp[,2]/100000 # scale SSC-H so flowPeaks works properly
      flowPeaks.Res3 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res3$peaks$cid, function(x){ length(which(flowPeaks.Res3$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res3$peaks$cid[which(cluster.size > 0.0002)]
      
      mu.short <- flowPeaks.Res3$peaks$mu[sig.clusters, ]
      eosinophils.clusterid <-sig.clusters[which((mu.short[, 2] > (min(mu.short[,2]) + 0.5*(max(mu.short[, 2]) -  min(mu.short[,2])))) & (mu.short[, 1] > (max(mu.short[, 1]) - 1.5))  )]
      
      if(length(eosinophils.clusterid) > 1){
        eosinophils.clusterid <- eosinophils.clusterid[which(cluster.size[eosinophils.clusterid] > 0.15*max(cluster.size[eosinophils.clusterid]))]
        if(length(eosinophils.clusterid) > 1){
          mu.short <- flowPeaks.Res3$peaks$mu[eosinophils.clusterid, ]
          eosinophils.clusterid <- eosinophils.clusterid[which.min((mu.short[,1] - CD11b.upper.gate + 0.3)^2 + (mu.short[, 2] - 1.8)^2)]
        }
      }
      
      eosinophils.idx <- which(flowPeaks.Res3$peaks.cluster ==  eosinophils.clusterid)
      #eosinophils <- getflowFrame(NOT.monocytes)
      eosinophils <- getflowFrame(temp.flowD)
      eosinophils@exprs <- eosinophils@exprs[eosinophils.idx, ]
      
      CD11b.eosino.gate <- min(eosinophils@exprs[, channels.ind["CD11b"]])
      eosinophils.flowD <- flowDensity(eosinophils, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, T), gates = c(CD11b.eosino.gate, min(eosinophils@exprs[, scat.chans['SSC-H']])), ellip.gate = T, alpha = 0.85)
      
      CD11b.eosino.gate <- deGate(eosinophils.flowD, channel = c(channels.ind["CD11b"]), twin.factor = 0.9)
      if(is.null(names(CD11b.eosino.gate)) & ((max(getflowFrame(eosinophils.flowD)@exprs[, channels.ind["CD11b"]]) - CD11b.eosino.gate) > 0.8)){
        eosinophils.flowD <- flowDensity(eosinophils, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, T), gates = c(CD11b.eosino.gate, min(eosinophils@exprs[, scat.chans['SSC-H']])), ellip.gate = T, alpha = 0.85)
      }
      
      NOT.eosinophils <- notSubFrame(NOT.monocytes, channels =  c(channels.ind["CD11b"], scat.chans['SSC-H']), filter = eosinophils.flowD@filter )
      
      rm(flowPeaks.Res3)
      
      # Gating NOT(Eosinophils) to get CD161+ and CD161- ------------------------------------------------
      
      # TCP, CIPHE, BCM all have edge events at the CD161+ side. Remove
      #maxCD161.value <- max(getflowFrame(NOT.eosinophils)@exprs[, c(channels.ind['CD161'])])
      #NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(0.99*maxCD161.value, NA))
      
      CD19lo.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD19']), use.percentile = T, percentile = 0.005) - 0.5
      CD161hi.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD161']), use.percentile = T, percentile = 0.995) + 0.5
      NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, T), gates = c(CD161hi.gate, CD19lo.gate))
      
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
      
      
      CD161.gate <- deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.05) 
      
      # Additional code mainly for BCM where the compensation is such that there are a significant amount of events below the main peak
      maxDens <- density(CD19neg@exprs[, c(channels.ind["CD161"])])
      CD161.gate2 <- deGate(CD19neg, channel = c(channels.ind["CD161"]), after.peak = T, twin.factor = 0.5) 
      if(is.null(names(CD161.gate2))){
        CD161.gate <- min(CD161.gate, CD161.gate2)
      } 
      # because after.peak = T does not alway give the correct result if the CD161+ population is not recognized as a peak
      if(CD161.gate <  maxDens$x[which.max(maxDens$y)]){ 
        CD161.gate <- deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9) + 0.05
      }else if((CD161.gate - maxDens$x[which.max(maxDens$y)]) > 1.5){ # This suggests that the large CD161 peak is slowly decreasing and thus we need to increase alpha
        CD161.gate <- deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9) 
      }
      
      
      CD161pos.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      if(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count > 0.085){
        CD161.gate2 <- c(deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.05),
                         deGate(CD19neg, channel = c(channels.ind["CD161"]), all.cuts = T, upper = T, percentile = NA, tinypeak.removal = 0.05))
        if(max(CD161.gate2) > CD161.gate){
          CD161.gate <- min(CD161.gate2[which(CD161.gate2 > (CD161.gate + 0.1))])
          CD161pos.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
        }
      }
      NOT.CD161pos <- notSubFrame(NOT.eosinophils, channels =  c(channels.ind["CD161"], channels.ind["CD19"]), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      
      # For Jax the CD161 and CD317 channels are shared. I need to remove CD317+ from CD161+
      if(centre == 'Jax'){
        maxDens <- density(getflowFrame(CD161pos.flowD)@exprs[, c(channels.ind['CD161'])])
        CD161.peak <- findpeaks(maxDens$y) # maxDens$x[which.max(maxDens$y)]
        if(length(CD161.peak) > 4){
          CD161.peak <- max(maxDens$x[CD161.peak[which(CD161.peak[,1] > 0.6*max(CD161.peak[,1])), 2]])
        }else{
          CD161.peak <- maxDens$x[CD161.peak[2]]
        }
        CD161CD317.gate <- c(deGate(CD19neg, channel = c(channels.ind['CD161']), all.cuts = T, tinypeak.removal = 0.0005),
                             deGate(CD19neg, channel = c(channels.ind['CD161']), use.upper = T, upper = T))
        CD161CD317.gate <- min(CD161CD317.gate[which(CD161CD317.gate > CD161.peak)])
        CD161pos.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(CD161CD317.gate, NA))
        NOT.CD161pos <- flowDensity(NOT.CD161pos, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(CD161CD317.gate, NA))
      }
      
      
      rm(flowPeaks.Res4)
      
      # Gating CD161+ to get NK/NKT -----------------------------------------------------------------------------------
      
      CD5Ly6G.gate <- deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), tinypeak.removal = 0.01, upper = T, percentile = NA)
      maxDens <- density(getflowFrame(CD161pos.flowD)@exprs[, c(channels.ind["CD5/Ly6G"])])
      maxDens <- smooth.spline(maxDens, spar = 0.7) # smoothing pretty aggressively b/c for BCm sometimes I get a double peak in the NK population
      NK.CD5.peak <- findpeaks(maxDens$y)
      NK.CD5.peak <- NK.CD5.peak[which(NK.CD5.peak[,1] > 0.2*max(NK.CD5.peak[,1])),]
      if(length(NK.CD5.peak) > 4){
        NK.CD5.peak <- maxDens$x[min(NK.CD5.peak[,2])] # pick left-most peak
      }else{
        NK.CD5.peak <- maxDens$x[NK.CD5.peak[2]]
      }
      CD5Ly6G.alt.gate <- 2*NK.CD5.peak - deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), use.upper = T, upper = F) + 0.05
      if((CD5Ly6G.gate < NK.CD5.peak)|((CD5Ly6G.gate - CD5Ly6G.alt.gate) > 0.3)){
        CD5Ly6G.gate <- CD5Ly6G.alt.gate
      } 
      
      
      additionalNOTCD161pos <- NULL
      NKcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(F, NA), gates = c(CD5Ly6G.gate, NA))
      CD161.NK.gate <- deGate(NKcells.flowD, channel = c(channels.ind["CD161"]))
      maxDens <- density(getflowFrame(CD161pos.flowD)@exprs[, c(channels.ind['CD161'])])
      NK.CD161.peak <- findpeaks(maxDens$y)   #maxDens$x[which.max(maxDens$y)]
      if(length(NK.CD161.peak) > 4){ # needed in case there is a CD161+ peak below the NK peak that is larger in magnitude
        NK.CD161.peak <- maxDens$x[max(NK.CD161.peak[which(NK.CD161.peak[,1] > 0.4*max(NK.CD161.peak[,1])), 2])]
      }else{
        NK.CD161.peak <- maxDens$x[NK.CD161.peak[2]]
      }
      if(is.null(names(CD161.NK.gate))){
        if(CD161.NK.gate < (NK.CD161.peak - 0.1)) {
          NKcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(F, T), gates = c(CD5Ly6G.gate, CD161.NK.gate))
          additionalNOTCD161pos <-flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(F, F), gates = c(CD5Ly6G.gate, CD161.NK.gate))
        }
      }
      NKTcells.flowD <- flowDensity(CD161pos.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD161"]), position = c(T, NA), gates = c(CD5Ly6G.gate, NA))
      
      # Try to use NK population to get the CD11b and Ly6C gates (tried using whole CD161+ population and using only NK cells seems better)
      Ly6C.gate <- min(deGate(getflowFrame(NKcells.flowD), channel = c(channels.ind["Ly6C"]), all.cuts = T))
      
      temp.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(NA, F), gates = c(NA, Ly6C.gate))
      CD11b.gate <- deGate(getflowFrame(temp.flowD), channel = c(channels.ind["CD11b"]))
      if(!is.null(names(CD11b.gate))){ # sometimes I need to find this man
        maxDens <- density(getflowFrame(temp.flowD)@exprs[, c(channels.ind['CD11b'])])
        peak.lcns <- findpeaks(maxDens$y)
        peak.lcns <- peak.lcns[which(peak.lcns[,1] > 0.1*max(peak.lcns[,1])), ]
        if(length(peak.lcns) > 4){
          n <- nrow(peak.lcns)
          largest.peaks <- sort(peak.lcns[,1], index.return = TRUE)
          largest.peaks <- peak.lcns[largest.peaks$ix[(n-1):n], 2]
          CD11b.gate <- maxDens$x[which.min(maxDens$y[min(largest.peaks):max(largest.peaks)]) + min(largest.peaks) - 1]
        }else{
          CD11b.gate <- 2*maxDens$x[peak.lcns[2]] - deGate(getflowFrame(temp.flowD), channel = c(channels.ind["CD11b"]), use.upper = T, upper = T, alpha = 0.5)
        }
      }
      
      # Testing... NOT(CD161) = NOT(NK or NKT cells) Adding cells back helps me to get a better/cleaner RP mac population
      if(!is.null( additionalNOTCD161pos)){
        NOT.CD161pos@flow.frame@exprs <- rbind(getflowFrame(NOT.CD161pos)@exprs, getflowFrame(additionalNOTCD161pos)@exprs)
        NOT.CD161pos@cell.count <- additionalNOTCD161pos@cell.count + NOT.CD161pos@cell.count
      }
      
      # Gating NOT(CD161+) to get T cells ---------------------------------------------------------------------
      
      temp <- getflowFrame(NOT.CD161pos)@exprs[,  c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"])]
      flowPeaks.Res5 <- flowPeaks(temp)
      
      cluster.size <- sapply(flowPeaks.Res5$peaks$cid, function(x){ length(which(flowPeaks.Res5$peaks.cluster == x))})
      cluster.size <- cluster.size/nrow(temp)
      sig.clusters <- flowPeaks.Res5$peaks$cid[which(cluster.size > 0.01)]
      
      # sometimes in the CIPHE data there are some MHCII-- events (not sure why or what they are), but it screws up the T cell
      # gating, so I am removing any clusters corresponding to those events, if they exist
      MHCIIlo.gate <- deGate(NOT.CD161pos, channel = c(channels.ind['MHCII']), use.upper = T, upper = F)
      sig.clusters <- sig.clusters[which(flowPeaks.Res5$peaks$mu[sig.clusters, 1] > MHCIIlo.gate)]
      
      mu.short <- flowPeaks.Res5$peaks$mu[sig.clusters,]
      
      MHCIImid.gate <- deGate(NOT.CD161pos, channel = c(channels.ind['MHCII']))
      if(is.null(names(MHCIImid.gate)) & (MHCIImid.gate > MHCIIlo.gate + 0.2)){
        Tcell.clusterid <- sig.clusters[which((mu.short[, 2] > (min(mu.short[, 2] + 0.2)) & (mu.short[, 1] < MHCIImid.gate)))]
      }else{
        Tcell.clusterid <- sig.clusters[which.max(2*mu.short[, 2]^2 + (4.5 - mu.short[, 1])^2)]
        Tcell.clusterid <- union(Tcell.clusterid, sig.clusters[which((mu.short[, 2] > (min(mu.short[, 2] + 0.2)) & (mu.short[, 1] < mu.short[Tcell.clusterid, 1] + 0.2)))])
      }
      
      # remove low density regions of the clusters
      fpc <- assign.flowPeaks(flowPeaks.Res5, flowPeaks.Res5$x, tol = 0.015, fc = 0)
      
      Tcell.idx <- which(fpc %in% Tcell.clusterid)
      Tcells <- getflowFrame(NOT.CD161pos)
      Tcells@exprs <- Tcells@exprs[Tcell.idx, ]
      
      CD5lo.threshold <- min(mu.short[, 2]) + 0.4*(max(mu.short[,2] - min(mu.short[,2]))) 
      maxDens <- density(Tcells@exprs[, c(channels.ind['CD5/Ly6G'])])
      CD5lo.threshold.idx <- which.min(abs(maxDens$x - CD5lo.threshold))
      Tcell.peak <- maxDens$x[which.max(maxDens$y[CD5lo.threshold.idx:length(maxDens$y)]) + CD5lo.threshold.idx - 1]
      
      # # If The CD5-MHCII- population is included in the T cell cluster, then it is easier to find if I don't remove the low density regions.
      Tcell2.idx <- which(flowPeaks.Res5$peaks.cluster %in%  Tcell.clusterid)
      Tcells2 <- getflowFrame(NOT.CD161pos)
      Tcells2@exprs <- Tcells2@exprs[Tcell2.idx, ]
      
      CD5.gate <- deGate(Tcells2, channel = c(channels.ind["CD5/Ly6G"]), tinypeak.removal = 0.01)
      rm(Tcells2)
      
      if(CD5.gate < (Tcell.peak - 0.5)){
        Tcells <- getflowFrame(NOT.CD161pos)
        Tcell.idx <- which((flowPeaks.Res5$peaks.cluster %in%  Tcell.clusterid) & (Tcells@exprs[, c(channels.ind["CD5/Ly6G"])] > CD5.gate))
        Tcell.idx <- which((fpc %in%  Tcell.clusterid) & (Tcells@exprs[, c(channels.ind["CD5/Ly6G"])] > CD5.gate))
        Tcells@exprs <- Tcells@exprs[Tcell.idx, ]
      }
      
      # Don't use T cells as is, use the convex hull to get the T cells (this avoids some strange effects of how flowPeaks combines clusters)
      data.new <- exprs(Tcells)[, c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"])]
      hpts <- chull(x = data.new[, 1], y = data.new[, 2])
      hpts <- c(hpts, hpts[1])
      Tcells.filter <- data.new[hpts, ]
      
      Tcells.flowD <- SubFrame(NOT.CD161pos, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), filter = Tcells.filter)
      Tcells <- getflowFrame(Tcells.flowD)
      
      Ly6Cdensity.Tcells <- density(Tcells@exprs[, c(channels.ind['Ly6C'])])
      # Using the Ly6C gate from the NK/NKT cells for this!!!
      Ly6CposTcells.flowD <- flowDensity(Tcells, channels = c(channels.ind['Ly6C'], scat.chans['SSC-A']), position = c(T, NA), gates = c(Ly6C.gate, NA))
      
      NOT.Tcells.flowD <- notSubFrame(NOT.CD161pos, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), filter = Tcells.filter)
      NOT.Tcells <- getflowFrame(NOT.Tcells.flowD)
      
      # remove some debris from NOT.Tcells
      temp2 <- Tcells@exprs[which.max(Tcells@exprs[, channels.ind['MHCII']]), c(channels.ind['MHCII'], channels.ind['CD5/Ly6G'])]
      temp <- Tcells@exprs[which.min(Tcells@exprs[, channels.ind['CD5/Ly6G']]), c(channels.ind['MHCII'], channels.ind['CD5/Ly6G'])]
      idx.to.remove <- union(which((NOT.Tcells@exprs[, c(channels.ind['CD5/Ly6G'])] > temp2[2]) & (NOT.Tcells@exprs[, c(channels.ind['MHCII'])] < temp2[1])),
                             which((NOT.Tcells@exprs[, c(channels.ind['CD5/Ly6G'])] > temp[2]) & (NOT.Tcells@exprs[, c(channels.ind['MHCII'])] < temp[1])))
      NOT.Tcells@exprs <- NOT.Tcells@exprs[-idx.to.remove, ]
      
      rm(flowPeaks.Res5)
      
      # Gate MHCII+ events -------------------------------------------------------------------------------------------
      
      MHCII.gate <- max(deGate(NOT.Tcells, channel = c(channels.ind['MHCII']), all.cuts = T, tinypeak.removal = 0.005, upper = F, percentile = NA))
      MHCII.gate.alt <- deGate(NOT.Tcells, channel = c(channels.ind['MHCII']), use.upper = T, upper = F, tinypeak.removal = 0.5, alpha = 0.05)
      if(!is.null(names(MHCII.gate))){
        MHCII.gate <- MHCII.gate.alt
      }else if(abs(MHCII.gate - MHCII.gate.alt) > 0.15){
        MHCI.gate <- MHCII.gate.alt
      }
      
      MHCIIpos.flowD <- flowDensity(NOT.Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(T, NA), gates = c(MHCII.gate, NA))
      
      
      
      # Gating MHCII+ to get B cells, RPmacs and cDcs --------------------------------------------------
      
      # get RP macs, if F4/80 exists
      if(!is.na(channels.ind['F4/80'])){
        
        F4.gate <- deGate(MHCIIpos.flowD@flow.frame, channel = c(channels.ind['F4/80']), use.upper = T, upper = T, tinypeak.removal = 0.9) - 0.12 
        RPmac.flowD <- flowDensity(getflowFrame(MHCIIpos.flowD), channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(T, NA), gates = c(F4.gate, NA))
        
        RPmac <- getflowFrame(RPmac.flowD)
        RPmac.ff.idx <- RPmac.flowD@index
        ff <- getflowFrame(MHCIIpos.flowD)
        
        temp <- RPmac@exprs[, c(channels.ind["F4/80"], channels.ind["MHCII"])]
        flowPeaks.Res6 <- flowPeaks(temp)
        
        cluster.size <- sapply(flowPeaks.Res6$peaks$cid, function(x){ length(which(flowPeaks.Res6$peaks.cluster == x))})
        cluster.size <- cluster.size/nrow(temp)
        sig.clusters <- flowPeaks.Res6$peaks$cid[which(cluster.size > 0.01)]
        
        if(length(sig.clusters) == 1){
          # Try to deGate again to remove some of the NOT(RPmac) cells
          F4.gate2 <-  deGate(RPmac.flowD, channel = c(channels.ind['F4/80']), use.upper = T, upper = T, tinypeak.removal = 0.2)
          if((F4.gate2 - F4.gate) < 0.61){
            RPmac.flowD <- flowDensity(MHCIIpos.flowD, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(T, NA), gates = c(F4.gate2, NA))
            
            RPmac <- getflowFrame(RPmac.flowD)
            RPmac.ff.idx <- RPmac.flowD@index
            
            temp <- RPmac@exprs[, c(channels.ind["F4/80"], channels.ind["MHCII"])]
            flowPeaks.Res6 <- flowPeaks(temp)
            
            cluster.size <- sapply(flowPeaks.Res6$peaks$cid, function(x){ length(which(flowPeaks.Res6$peaks.cluster == x))})
            cluster.size <- cluster.size/nrow(temp)
            sig.clusters <- flowPeaks.Res6$peaks$cid[which(cluster.size > 0.01)]
          }
        }
        
        mu.short <- flowPeaks.Res6$peaks$mu[sig.clusters,]
        maxDens <- density(getflowFrame(MHCIIpos.flowD)@exprs[, c(channels.ind['MHCII'])])
        if(length(sig.clusters) > 1){
          
          RPmac.cluster <- sig.clusters[which.min((max(mu.short[, 1]) - mu.short[, 1])^2 + (mu.short[, 2] - maxDens$x[which.max(maxDens$y)] + 0.4)^2)]
          additional.clusters <- which((mu.short[, 2] < flowPeaks.Res6$peaks$mu[RPmac.cluster, 2]) & (mu.short[, 1] > (F4.gate + 0.25*(flowPeaks.Res6$peaks$mu[RPmac.cluster, 1] - F4.gate))))
          if(length(additional.clusters) > 0){
            RPmac.cluster <- union(RPmac.cluster, sig.clusters[additional.clusters])
          }
          RPmac.idx <- which(flowPeaks.Res6$peaks.cluster %in%  RPmac.cluster)
          RPmac@exprs <- RPmac@exprs[RPmac.idx, ]
          ff@exprs <- ff@exprs[-RPmac.ff.idx[RPmac.idx],]
          
        }else{
          # Either the 1 cluster is the RPmac population, OR it is NOT(RPmac) and there are barely any RP mac events
          
          if(quantile(getflowFrame(RPmac.flowD)@exprs[,c(channels.ind['MHCII'])], c(0.5)) > (maxDens$x[which.max(maxDens$y)] - 0.1)){
            # Assume the cluster is NOT(RPmac)
            rot <- rotate.data(getflowFrame(MHCIIpos.flowD), c(channels.ind["F4/80"], channels.ind["MHCII"]), theta = -pi/6)$data
            rot.gate <- deGate(rot, channel = c(channels.ind["F4/80"]), use.upper = T, upper = T, alpha = 0.05)
            
            RPmac.flowD <- flowDensity(rot, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(T, NA), gates = c(rot.gate, NA))
            
            RPmac.flowD@filter <- rotate.data(RPmac.flowD@filter, c(channels.ind["F4/80"], channels.ind["MHCII"]), theta = pi/6)$data
            RPmac.flowD@flow.frame <- rotate.data(getflowFrame(RPmac.flowD), c(channels.ind["F4/80"], channels.ind["MHCII"]),theta = pi/6)$data
            RPmac <- getflowFrame(RPmac.flowD)
            ff@exprs <- ff@exprs[-RPmac.flowD@index,]
            
          }else{
            # Assume whole cluster is RPmac
            rot <- rotate.data(getflowFrame(RPmac.flowD), c(channels.ind["F4/80"], channels.ind["MHCII"]), theta = -pi/4)$data
            rot.gate <- deGate(rot, channel = c(channels.ind["F4/80"]), use.upper = T, upper = F, alpha = 0.5, tinypeak.removal = 0.5)
            
            RPmac.flowD <- flowDensity(rot, channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), position = c(T, NA), gates = c(rot.gate, NA))
            
            RPmac.flowD@filter <- rotate.data(RPmac.flowD@filter, c(channels.ind["F4/80"], channels.ind["MHCII"]), theta = pi/4)$data
            RPmac.flowD@flow.frame <- rotate.data(getflowFrame(RPmac.flowD), c(channels.ind["F4/80"], channels.ind["MHCII"]),theta = pi/4)$data
            RPmac <- getflowFrame(RPmac.flowD)
            ff@exprs <- ff@exprs[-RPmac.flowD@index,]
            
            #ff@exprs <- ff@exprs[-RPmac.ff.idx, ]
          }
        }
        
        rm(flowPeaks.Res6)
        
      }else{
        ff <- MHCIIpos.flowD@flow.frame
      }
      
      ff_CD11CFMO <- ff
      save(ff_CD11CFMO, file = paste0('/data/IMPC/FMOs/', centre, '/', fname.date, '.Rdata'))
      
    }
  })
  
}

    