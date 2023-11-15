rm(list = ls())

library('colorRamps')
library('flowCore')
library('plyr')
library('doMC')
library('e1071') # for flowCut
library('Cairo') # for flowCut
library('flowDensity')
library('pracma') # for findpeaks
library('stringr')
library('flowPeaks')
library('MASS')
library('tclust') # to get FoB cluster

setwd('/data/IMPC')

#no_cores <- detectCores() - 2
no_cores <- 12
registerDoMC(no_cores)

source("/data/IMPC/scripts/P2gating_additionalfunctions.R")
source("/data/IMPC/scripts/flowCut_20170331.R")


# load IMPReSS ID info ------------------------------------------------------------
load('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/immpressIDconversion.Rdata')
colnames(immpressIDconversion) <- c("IMPReSS.id", "Genotype")
immpressIDconversion <- immpressIDconversion[47:75, ] # limit to only Panel 2 

centre <- 'TCP'
load('/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel2/Results/store.allFCS.Rdata')
scatterplot.dir <- '/data/IMPC/ScatterPlots/TCP'
file.names <- data.frame(store.allFCS, stringsAsFactors = F)

centre <- 'CIPHE'
load('/mnt/f/FCS data/IMPC/IMPC-Results/CIPHE/Panel2/Results/store.allFCS.Rdata')
scatterplot.dir <- '/data/IMPC/ScatterPlots/CIPHE'
file.names <- data.frame(store.allFCS, stringsAsFactors = F)

centre <- 'BCM'
load('/mnt/f/FCS data/IMPC/IMPC-Results/BCM/Panel2/Results/store.allFCS.Rdata')
scatterplot.dir <- '/data/IMPC/ScatterPlots/BCM'
file.names <- data.frame(store.allFCS, stringsAsFactors = F)

centre <- 'Jax'
load('/mnt/f/FCS data/IMPC/IMPC-Results/Jax/Panel2/Results/store.allFCS.Rdata')
scatterplot.dir <- '/data/IMPC/ScatterPlots/Jax'
file.names <- data.frame(store.allFCS, stringsAsFactors = F)

Sys.time()

# # Get wildtypes for BCM
wildtypef.idx <- which((store.allFCS[, 'Genotype'] == '+_+') & (store.allFCS[, 'Gender'] == 'Female'))
wildtypem.idx <- which(((store.allFCS[, 'Genotype'] == '+_+') & (store.allFCS[, 'Gender'] == 'Male'))|(store.allFCS[, 'Genotype'] == 'Y_+'))
sample.idx <- union(union(853:nrow(file.names), wildtypef.idx), wildtypem.idx)

# Get wildtypes for TCP
#wildtype.idx <- which((store.allFCS[, 'Genotype'] == 'WT') & (store.allFCS[, 'Gender'] == 'F'))
#wildtype.idx <- which((store.allFCS[, 'Genotype'] == 'WT'))


immpressData <- ldply(sample.idx, function(i){ 
  
  x <- file.names[i,]  
  immpressData.FCSfile <- matrix(data = NA, nrow = 1, ncol = 29)
  colnames(immpressData.FCSfile) <- immpressIDconversion[, 'IMPReSS.id']
  passed.flowCut <- F
  
  try({
    
    fname <- x$FCS.files
    if((centre == 'TCP')|(centre == 'BCM')|(centre == 'Jax')){
      fname.full <- paste(x$Path, x$Panel.Organ.Folder, x$FCS.files, sep = '/')
    }else{ # CIPHE
      fname.full <- paste(x$Path, x$FCS.files, sep = '/')
    }
    f <- read.FCS(fname.full)
 
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
    
    original.no.cells <- nrow(f)
    
    # Remove scatter margins and compensate ---------------------------------------------
    # Removing margin events in Scatter channels
    f <- removeMargins(f, chans = scat.chans, verbose = F)
    #Removing negative values in scatter channels
    f <- removeMargins(f, chans = scat.chans, debris = T, neg = T, verbose = F)
    
    f <- compensateIMPC(f, fname.full, centre = centre, panel.no = 2)
    
    # Transformation
    load(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/', centre, "/Panel2/Results/lgl.Rdata"))

    fT <- transform(f, lgl)
    channels.to.clean <- setdiff(channels.ind, channels.ind['Time'])
    fclean <- flowCut(fT, segment = floor(nrow(f)*5/3000), CleanChan.loc = channels.to.clean, FileID = fname, directory = paste0('/data/IMPC/ScatterPlots/', centre, '/flowCut/'), Plot = 'Failed Only')
    fclean.ind <- fclean$ind
    passed.flowCut <- fclean$data['Has the file passed',1]
    rm(fclean)
    
    # For CIPHE The F4/80 channel looks best if I use the parameters from the FCS file for transformation so I'll do that
    if(centre == 'CIPHE'){
      f.temp <- logiclTransformCiphe(f, markers.transform = colnames(f)[channels.ind['F4/80']])
    }
    f <- transform(f, lgl)
    if(centre == 'CIPHE'){
      f@exprs[,c(channels.ind['F4/80'])] <- f.temp@exprs[,c(channels.ind['F4/80'])]
      rm(f.temp)
    }
    if(length(fclean.ind) > 0){
      f@exprs <- f@exprs[-fclean.ind, ]
    }
    
    # Quality Gate (live/dead, singlets gating) ---------------------------------------------------
    
    results <- qualityGate(f, scat.chans, channels.ind, centre, panel.no = 2)
    live.flowD <- results$live
    FSCsinglets.flowD <- results$FSCsinglets
    singlets.flowD <- results$singlets
    scat.chans <- results$scat.chans
    rm(results)
    
    # Panel 2 Gating -------------------------------------------------------------------------------
    
    singlets <- getflowFrame(singlets.flowD)
    
    # Gating singlets to get Granulocytes ------------------------------------------------------------

    singlets0 <- singlets # Just to check the CIPHE/BCM filtering of edge events
    
    if((centre == 'CIPHE')){ # For the CIPHE there are some days where there are edge events along the top edge of the CD11b channel
      #max.CD11b.value <- max(singlets@exprs[, c(channels.ind["CD11b"])])
      max.CD11b.value <- quantile(singlets@exprs[, c(channels.ind["CD11b"])], c(.999)) + 0.25
      singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["CD11b"])] < max.CD11b.value), ]
    }else if(centre == 'TCP'){
      max.CD21.value <- max(singlets@exprs[, c(channels.ind["CD21/CD35"])])
      singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["CD21/CD35"])] < 0.99*(max.CD21.value)), ]
    }else if(centre == 'BCM'){
      max.CD11c.value <- max(singlets@exprs[, c(channels.ind["CD11c"])])
      singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["CD11c"])] < 0.99*(max.CD11c.value)), ]
      max.CD11b.value <- max(singlets@exprs[, c(channels.ind["CD11b"])])
      singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["CD11b"])] < 0.99*(max.CD11b.value)), ]
      max.Ly6C.value <- min(quantile(singlets@exprs[, c(channels.ind["Ly6C"])], c(.999)) + 0.25, max(singlets@exprs[, c(channels.ind["Ly6C"])]))
      singlets@exprs <- singlets@exprs[which(singlets@exprs[, c(channels.ind["Ly6C"])] < (max.Ly6C.value)), ]
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
      granulo.clusterid <- granulo.clusterid[which(cluster.size[granulo.clusterid] < 0.15)]
      granulo.clusterid <- granulo.clusterid[which(cluster.size[granulo.clusterid] > 0.174*max(cluster.size[granulo.clusterid]))]
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
    monocyte.Ly6C.target <- deGate(NOT.granulocytes, channel = c(channels.ind['Ly6C']), use.upper = T, upper = T) - 0.25
    if(monocyte.Ly6C.target < 2){
      maxDens2 <- density(getflowFrame(NOT.granulocytes)@exprs[, channels.ind['CD11b']])
      temp2.flowD <- flowDensity(getflowFrame(NOT.granulocytes), channels = c(channels.ind["Ly6C"], channels.ind['CD11b']), position = c(T, F), gates = c(monocyte.Ly6C.target, (maxDens2$x[which.max(maxDens2$y)] + 0.2)))
      monocyte.Ly6C.target <- deGate(temp2.flowD, channel = c(channels.ind['Ly6C']), use.upper = T, upper = T, tinypeak.removal = 0.001)
    }
    #monocyte.CD11b.target <- deGate(NOT.granulocytes, channel = c(channels.ind['CD11b']), use.upper = T, upper = T) - 0.25
    
    #temp <- getflowFrame(NOT.granulocytes)@exprs[,  c(channels.ind["Ly6C"], channels.ind["CD11b"])]
    temp <- getflowFrame(temp.flowD)@exprs[,  c(channels.ind["Ly6C"], channels.ind["CD11b"])]
    flowPeaks.Res2 <- flowPeaks(temp)
    
    cluster.size <- sapply(flowPeaks.Res2$peaks$cid, function(x){ length(which(flowPeaks.Res2$peaks.cluster == x))})
    cluster.size <- cluster.size/nrow(temp)
    sig.clusters <- flowPeaks.Res2$peaks$cid[which(cluster.size > 0.001)]
    
    mu.short <- flowPeaks.Res2$peaks$mu[sig.clusters,]
    monocyte.clusterid <-sig.clusters[which((mu.short[, 2] > (min(mu.short[,2]) + 0.5*(max(mu.short[, 2]) -  min(mu.short[,2])))) & (mu.short[, 1] > (max(mu.short[, 1]) - 1))  )]
    if(length(monocyte.clusterid) > 1){
      monocyte.clusterid <- monocyte.clusterid[which(cluster.size[monocyte.clusterid] > 0.15*max(cluster.size[monocyte.clusterid]))]
      monocyte.clusterid <- monocyte.clusterid[which.min((flowPeaks.Res2$peaks$mu[monocyte.clusterid, 2] - max(mu.short[, 2]))^2 + (flowPeaks.Res2$peaks$mu[monocyte.clusterid, 1] - monocyte.Ly6C.target)^2)]
    }
    
    monocytes.idx <- which(flowPeaks.Res2$peaks.cluster ==  monocyte.clusterid)
    monocytes <- getflowFrame(temp.flowD)
    monocytes@exprs <- monocytes@exprs[monocytes.idx, ]
    
    Ly6C.monocyte.gate <- deGate(monocytes, channel = channels.ind["Ly6C"], upper = F, tinypeak.removal = 0.5, percentile = NA)
    Ly6C.monocyte.alt.gate <- deGate(monocytes, channel = channels.ind["Ly6C"], tinypeak.removal = 0.05)
    if(is.null(names(Ly6C.monocyte.alt.gate))){ # So if the monocytes currently have 2 Ly6C peaks give option to split them up
      if(Ly6C.monocyte.alt.gate <  monocyte.Ly6C.target){ # If the Ly6C gate is not above the approximate monocyte location (this avoids gating to get only noise)
        Ly6C.monocyte.gate <- max(Ly6C.monocyte.alt.gate, Ly6C.monocyte.gate)
      }
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
    CD11b.temp.gate <- deGate(NOT.monocytes, channel = c(channels.ind["CD11b"]), upper = T, percentile = NA, tinypeak.removal = 0.9)
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
    
    if(eosinophils.flowD@cell.count > 50){
      CD11b.eosino.gate <- deGate(eosinophils.flowD, channel = c(channels.ind["CD11b"]), twin.factor = 0.9)
      if(is.null(names(CD11b.eosino.gate)) & ((max(getflowFrame(eosinophils.flowD)@exprs[, channels.ind["CD11b"]]) - CD11b.eosino.gate) > 0.8)){
        eosinophils.flowD <- flowDensity(eosinophils, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), position = c(T, T), gates = c(CD11b.eosino.gate, min(eosinophils@exprs[, scat.chans['SSC-H']])), ellip.gate = T, alpha = 0.85)
      }
    }
    NOT.eosinophils <- notSubFrame(NOT.monocytes, channels =  c(channels.ind["CD11b"], scat.chans['SSC-H']), filter = eosinophils.flowD@filter )
   
    rm(flowPeaks.Res3)
    
    # Gating NOT(Eosinophils) to get CD161+ and CD161- ------------------------------------------------
    
    # TCP, CIPHE, BCM all have edge events at the CD161+ side. Remove
    #maxCD161.value <- max(getflowFrame(NOT.eosinophils)@exprs[, c(channels.ind['CD161'])])
    #NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), position = c(F, NA), gates = c(0.99*maxCD161.value, NA))
    
    if(centre == 'TCP'){ # Strategy is the same as the other centers except this usese the MHCII channel, while other centres use the CD19 channel
      
      CD19lo.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['MHCII']), use.percentile = T, percentile = 0.005) - 0.5
      CD161hi.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['CD161']), use.percentile = T, percentile = 0.995) + 0.5
      NOT.eosinophils <- flowDensity(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["MHCII"]), position = c(F, T), gates = c(CD161hi.gate, CD19lo.gate))
      
      CD19.gate <- deGate(NOT.eosinophils, channel = c(channels.ind['MHCII']))
      
      temp <- getflowFrame(NOT.eosinophils)@exprs[,  c(channels.ind["CD161"], channels.ind['MHCII'])]
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
      
      
      CD161pos.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind['MHCII']), position = c(T, F), gates = c(CD161.gate, CD19.gate))
      if(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count > 0.086){
        CD161.gate2 <- c(deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.05),
                         deGate(CD19neg, channel = c(channels.ind["CD161"]), all.cuts = T, upper = T, percentile = NA, tinypeak.removal = 0.05))
        if(max(CD161.gate2) > (CD161.gate + 0.1)){
          CD161.gate <- min(CD161.gate2[which(CD161.gate2 > (CD161.gate + 0.1))])
          CD161pos.flowD <- flowDensity(CD19neg, channels = c(channels.ind["CD161"], channels.ind['MHCII']), position = c(T, F), gates = c(CD161.gate, CD19.gate))
        }
      }
      NOT.CD161pos <- notSubFrame(NOT.eosinophils, channels =  c(channels.ind["CD161"], channels.ind['MHCII']), position = c(T, F), gates = c(CD161.gate, CD19.gate))

      
    }else{
    
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
      if(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count > 0.086){
        CD161.gate2 <- c(deGate(CD19neg, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.05),
                         deGate(CD19neg, channel = c(channels.ind["CD161"]), all.cuts = T, upper = T, percentile = NA, tinypeak.removal = 0.05))
        if(max(CD161.gate2) > (CD161.gate + 0.1)){
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
    CD5Ly6G.alt.gate <- 2*NK.CD5.peak - deGate(CD161pos.flowD, channel = c(channels.ind["CD5/Ly6G"]), use.upper = T, upper = F, tinypeak.removal = 0.15) + 0.05
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
    Ly6C.gate <- deGate(getflowFrame(NKcells.flowD), channel = c(channels.ind["Ly6C"]), all.cuts = T, tinypeak.removal = 0.15, upper = T, percentile = NA)
    maxDens <- density(getflowFrame(NKcells.flowD)@exprs[, c(channels.ind['Ly6C'])])
    NKcell.peaks <- findpeaks(maxDens$y)
    NKcell.peaks <- NKcell.peaks[which(NKcell.peaks[,1] > 0.2*max(NKcell.peaks[,1])),]
    if(length(NKcell.peaks) > 4){
      NKcell.peaks <- maxDens$x[NKcell.peaks[1,2]] # pick peak with lowest x value
    }else{
      NKcell.peaks <- maxDens$x[NKcell.peaks[2]]
    }
    Ly6C.gate <- Ly6C.gate[which.min(abs(Ly6C.gate - NKcell.peaks - 0.4))]
    
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
        CD11b.gate <- deGate(getflowFrame(NKTcells.flowD), channel = c(channels.ind["CD11b"]), use.upper = T, upper = T, tinypeak.removal = 0.5, alpha = 0.5)
        #CD11b.gate <- 2*maxDens$x[peak.lcns[2]] - deGate(getflowFrame(temp.flowD), channel = c(channels.ind["CD11b"]), use.upper = T, upper = T, alpha = 0.5))
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
    
    CD5.gate <- deGate(Tcells2, channel = c(channels.ind["CD5/Ly6G"]), all.cuts = T, tinypeak.removal = 0.01)
    rm(Tcells2)
    
    if(min(CD5.gate) < (Tcell.peak - 0.25)){
      CD5.gate <- CD5.gate[which(CD5.gate < Tcell.peak)]
      CD5.gate <- CD5.gate[which.min(abs(CD5.gate - Tcell.peak + 0.5))]
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
    
    # get pDCs 
    # Commented out for now because this gating sucks as written.
    # if(centre == 'CIPHE'){
    #   CD317.gate <- deGate(NOT.Tcells, channel = c(channels.ind['CD317']), all.cuts = T, tinypeak.removal = 0.0001)
    #   maxDens <- density(NOT.Tcells@exprs[ ,c(channels.ind['CD317'])])
    #   peak.lcn <- maxDens$x[which.max(maxDens$y)] 
    #   CD317.extent<- deGate(NOT.Tcells, channel = c(channels.ind['CD317']), use.percentile = T, percentile = 0.9999)
    #   target.gate.lcn <- peak.lcn + 0.5*(CD317.extent - peak.lcn)
    #   CD317.gate <- CD317.gate[which.min(abs(CD317.gate - target.gate.lcn))]
    #   
    #   pDC.flowD <- flowDensity(NOT.Tcells, channels = c(channels.ind['CD317'], channels.ind['Ly6C']), position = c(T, NA), gates = c(CD317.gate, NA))
    #   temp <- NOT.Tcells
    #   temp@exprs <- temp@exprs[-pDC.flowD@index, ] 
    #   pDC.cellcount <- pDC.flowD@cell.count
    #   rm(pDC.flowD)
    # }else{
    #   temp <- NOT.Tcells
    # }
    
    # Gate MHCII+ events -------------------------------------------------------------------------------------------
    
    temp <- NOT.Tcells
    MHCII.gate <- max(deGate(temp, channel = c(channels.ind['MHCII']), all.cuts = T, tinypeak.removal = 0.005, upper = F, percentile = NA))
    MHCII.gate.alt <- deGate(temp, channel = c(channels.ind['MHCII']), use.upper = T, upper = F, tinypeak.removal = 0.5, alpha = 0.05)
    if(!is.null(names(MHCII.gate))){
      MHCII.gate <- MHCII.gate.alt
    }else if(abs(MHCII.gate - MHCII.gate.alt) > 0.15){
      MHCI.gate <- MHCII.gate.alt
    }

    MHCIIpos.flowD <- flowDensity(temp, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), position = c(T, NA), gates = c(MHCII.gate, NA))

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

    # Get cDCs
    
    # Try to get cDCs from flowPeaks clusters
    temp <- na.omit(ff@exprs[, c(channels.ind["CD19"], channels.ind["CD11c"])])
    flowPeaks.Res7 <- flowPeaks(temp)

    cluster.size <- sapply(flowPeaks.Res7$peaks$cid, function(x){ length(which(flowPeaks.Res7$peaks.cluster == x))})
    cluster.size <- cluster.size/nrow(temp)
    sig.clusters <- flowPeaks.Res7$peaks$cid[which(cluster.size > 0.001)]

    if((length(sig.clusters) == 1) | ((length(sig.clusters) == 2) & (min(cluster.size[sig.clusters]) < 0.008))){ # ONLY HAPPENS RARELY: In this case rely on solution based on deGate() only
      
      # first remove noise (any clusters  not in sig.clusters)
      cDC.idx <- which(flowPeaks.Res7$peaks.cluster %in%  sig.clusters) 
      cDC <- ff
      cDC@exprs <- na.omit(cDC@exprs)
      cDC@exprs <- cDC@exprs[cDC.idx, ]
      
      CD19.cDC.gate  <- deGate(cDC, channel = c(channels.ind['CD19']), upper = F, percentile = NA, tinypeak.removal = 0.9)
      
      cDC@exprs <- cDC@exprs[which(cDC@exprs[, c(channels.ind['CD19'])] < CD19.cDC.gate), ]
      
      temp.flowD <- flowDensity(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(T, NA), gates = c(CD19.cDC.gate, NA))
      RPmac.cDC.CD11c.target <- deGate(temp.flowD, channel = c(channels.ind['CD11c']), upper = T, percentile = NA, tinypeak.removal = 0.9, alpha = 0.3)
      if(centre != 'CIPHE'){
        RPmac.cDC.CD11c.target <- RPmac.cDC.CD11c.target + 0.18
      }
      
      # For TCP use FMOs to get the target cDC location
      # BCM also has daily FMOs, but there is too much within day variation to use the FMO control data as is. I would need to normalize
      # I cannot use the FMO data directly to set my cD11c gate because sometimes it sets the gate at a position that is obviously wrong (cuts cDC cluster in two)
      CD19.Bcell.gate <- NA
      ff_CD11CFMO <- NA
      if(centre == 'TCP'){
        try({
          load(file = paste0('/data/IMPC/FMOs/', centre, '/', x$Panel.Organ.Folder, '.Rdata'))
          CD19FMO.channel <- grep('CD19', pData(parameters(ff))$desc)
          CD11cFMO.channel <- grep('cd11c', tolower(pData(parameters(ff))$desc))
          CD19.CD11CFMO.gate <- deGate(ff_CD11CFMO, channel = CD19FMO.channel, tinypeak.removal = 0.005, upper = F, percentile = NA, after.peak = F)
          temp.flowD <- flowDensity(ff_CD11CFMO, channels = c(CD19FMO.channel, CD11cFMO.channel), position = c(F, NA), gates = c(CD19.CD11CFMO.gate, NA))
          temp.flowD <- flowDensity(temp.flowD, channels = c(CD19FMO.channel, CD11cFMO.channel), position = c(F, NA), gates = c(CD19.CD11CFMO.gate, NA), ellip.gate = T, scale = 0.97)
          RPmac.cDC.CD11c.target <- max(temp.flowD@filter[, 2])
        })
      }
        
      CD11c.cDC.gate <- c(deGate(cDC, channel = c(channels.ind['CD11c']), all.cuts = T, upper = F, percentile = NA),
                          deGate(cDC, channel = c(channels.ind['CD11c']), use.upper = T, upper = F , tinypeak.removal = 0.9))
      CD11c.cDC.gate <- CD11c.cDC.gate[which.min(abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target))]
      if(centre == 'TCP'){
        if((abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target) > 0.6) & !is.na(ff_CD11CFMO)){
          
          temp.flowD <- flowDensity(cDC, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(NA, T), gates = c(NA, CD11c.cDC.gate))
          temp.gate <- deGate(temp.flowD, channel = c(channels.ind['CD11c']), upper = F, percentile = NA)
          if((temp.gate > CD11c.cDC.gate) & (temp.gate < RPmac.cDC.CD11c.target)){
            CD11c.cDC.gate <- temp.gate
          }
          if(abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target) > 0.9){
            CD11c.cDC.gate <- RPmac.cDC.CD11c.target - 0.1
          }
        }
      }
      
      cDCs.flowD <- flowDensity(cDC, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, T), gates = c(CD19.cDC.gate + 0.05, CD11c.cDC.gate))
      cDCs.flowD <- flowDensity(cDCs.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, T), gates = c(CD19.cDC.gate + 0.05, CD11c.cDC.gate), ellip.gate = T, scale = 0.95)
        
    }else{ # NORMAL CASE: more than one significantly cluster found, can separate cDCs and B cells using flowPeaks 
      
      mu.short <- flowPeaks.Res7$peaks$mu[sig.clusters, ]
      
      Bcells.clusterid <- sig.clusters[which(cluster.size[sig.clusters] > 0.1)]
      Bcells.clusterid <- Bcells.clusterid[which(mu.short[Bcells.clusterid, 1] > (max(mu.short[Bcells.clusterid, 1]) - 0.4))]
      largestBcell.clusterid <- Bcells.clusterid[which.max(cluster.size[Bcells.clusterid])]
      
      # Won't use these B cells because sometimes things that are part of the cDC cluster should really be in the B cell cluster
      # But I need to get these cells to get a target gate location between RP mac and cDCs
      Bcells.temp <- ff
      Bcells.temp@exprs <- na.omit(Bcells.temp@exprs)
      Bcells.temp.idx <- which(flowPeaks.Res7$peaks.cluster %in% Bcells.clusterid)
      Bcells.temp@exprs <- Bcells.temp@exprs[Bcells.temp.idx, ]
      RPmac.cDC.CD11c.target <- deGate(Bcells.temp, channel = c(channels.ind['CD11c']), upper = T, percentile = NA, tinypeak.removal = 0.9, alpha = 0.3)
      if(centre != 'CIPHE'){
        RPmac.cDC.CD11c.target <- RPmac.cDC.CD11c.target + 0.18
      }
      
      # Start with NOT(B cells)
      cDC.clusterids <- setdiff(sig.clusters, Bcells.clusterid)
      cDC.clusterids <- cDC.clusterids[which(mu.short[cDC.clusterids, 1] < (mu.short[largestBcell.clusterid, 1] - 0.1))]
      
      cDC.idx <- which(flowPeaks.Res7$peaks.cluster %in%  cDC.clusterids)
      
      cDC <- ff
      cDC@exprs <- na.omit(cDC@exprs)
      if(length(cDC.clusterids) > 1){
        CD11cneg.CD19neg.cluster <- cDC.clusterids[which.min(mu.short[cDC.clusterids, 2])]
        CD11c.cDC.gate <- max(cDC@exprs[which(flowPeaks.Res7$peaks.cluster %in% CD11cneg.CD19neg.cluster), channels.ind['CD11c']])
      }else{
        CD11c.cDC.gate <- NA
      }
      cDC@exprs <- cDC@exprs[cDC.idx, ]
      
      # For TCP use FMOs to get the target cDC location
      # BCM also has daily FMOs, but there is too much within day variation to use the FMO control data as is. I would need to normalize
      # I cannot use the FMO data directly to set my cD11c gate because sometimes it sets the gate at a position that is obviously wrong (cuts cDC cluster in two)
      CD19.Bcell.gate <- NA
      ff_CD11CFMO <- NA
      if(centre == 'TCP'){
        try({
          load(file = paste0('/data/IMPC/FMOs/', centre, '/', x$Panel.Organ.Folder, '.Rdata'))
          CD19FMO.channel <- grep('CD19', pData(parameters(ff))$desc)
          CD11cFMO.channel <- grep('cd11c', tolower(pData(parameters(ff))$desc))
          CD19.CD11CFMO.gate <- deGate(ff_CD11CFMO, channel = CD19FMO.channel, tinypeak.removal = 0.005, upper = F, percentile = NA, after.peak = F)
          temp.flowD <- flowDensity(ff_CD11CFMO, channels = c(CD19FMO.channel, CD11cFMO.channel), position = c(F, NA), gates = c(CD19.CD11CFMO.gate, NA))
          temp.flowD <- flowDensity(temp.flowD, channels = c(CD19FMO.channel, CD11cFMO.channel), position = c(F, NA), gates = c(CD19.CD11CFMO.gate, NA), ellip.gate = T, scale = 0.97)
          RPmac.cDC.CD11c.target <- max(temp.flowD@filter[, 2])
        })
      }
      
      CD11c.cDC.gate <- c(CD11c.cDC.gate,
                          deGate(cDC, channel = c(channels.ind['CD11c']), all.cuts = T, upper = F, percentile = NA),
                          deGate(cDC, channel = c(channels.ind['CD11c']), use.upper = T, upper = F , tinypeak.removal = 0.9))
      CD11c.cDC.gate <- CD11c.cDC.gate[which.min(abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target - 0.05))]
      if(centre == 'TCP'){
        if((abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target) > 0.6) & !is.na(ff_CD11CFMO)){

          temp.flowD <- flowDensity(cDC, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(NA, T), gates = c(NA, CD11c.cDC.gate))
          temp.gate <- deGate(temp.flowD, channel = c(channels.ind['CD11c']), upper = F, percentile = NA)
          if((temp.gate > CD11c.cDC.gate) & (temp.gate < RPmac.cDC.CD11c.target)){
            CD11c.cDC.gate <- temp.gate
          }
          if(abs(CD11c.cDC.gate - RPmac.cDC.CD11c.target) > 0.9){
            CD11c.cDC.gate <- RPmac.cDC.CD11c.target - 0.1
          }
        }
      }else if(centre == 'BCM'){
        maxDens <- density(na.omit(ff@exprs[,c(channels.ind['CD11c'])]))
        if(CD11c.cDC.gate < (maxDens$x[which.max(maxDens$y)])){
          CD11c.cDC.gate <- RPmac.cDC.CD11c.target - 0.1
        }
      }
      
      CD19.cDC.gate <- deGate(cDC, channel = c(channels.ind['CD19']), upper = T, percentile = NA)
      maxDens <- density(cDC@exprs[,c(channels.ind['CD19'])])
      if(CD19.cDC.gate < (maxDens$x[which.max(maxDens$y)] + 0.2)){
        CD19.cDC.gate <- max(cDC@exprs[, c(channels.ind['CD19'])])
      }
      
      cDCs.flowD <- flowDensity(cDC, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(NA, T), gates = c(NA, CD11c.cDC.gate))
      cDCs.flowD <- flowDensity(cDCs.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, T), gates = c(CD19.cDC.gate , CD11c.cDC.gate), ellip.gate = T, scale = 0.95)
    
    } # end of getting cDCs
    
    # Get CD11b+ cDCs
    cDCs.CD11bdensity <- density(getflowFrame(cDCs.flowD)@exprs[, c(channels.ind['CD11b'])])
    if(centre == 'TCP'){
      CD11b.cDC.gate <- CD11b.gate
    }else{ # needed for BCM, but maybe not for other centres
      CD11b.cDC.gate <- deGate(cDCs.flowD, channel = c(channels.ind['CD11b']), all.cuts = T, tinypeak.removal = 0.01)
      if(!is.null(names(CD11b.cDC.gate))){
        CD11b.cDC.gate <- CD11b.gate 
      }else{
        CD11b.cDC.gate <- CD11b.cDC.gate[which.min(abs(CD11b.cDC.gate - CD11b.gate))]
      }
    }
    
    cDCs.CD11bpos.flowD <- flowDensity(cDCs.flowD, channels = c(channels.ind["CD11b"], channels.ind["CD19"]), position = c(T, NA), gates = c(CD11b.cDC.gate, NA))
    
    
    # Get B cells
    Bcells0.flowD <- notSubFrame(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), filter = cDCs.flowD@filter)

    if(is.na(CD19.Bcell.gate)){
      # Make it easier to detect a CD19-CD11c- peak, if it exists
      CD19.temp.gate <- deGate(Bcells0.flowD, channel = c(channels.ind['CD19']), use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.5)
      temp.flowD <- flowDensity(Bcells0.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, NA), gate = c(CD19.temp.gate - 0.05, NA))
      maxDens <- density(getflowFrame(temp.flowD)@exprs[,c(channels.ind['CD19'])])
      CD19.Bcell.gate <- max(deGate(temp.flowD, channel = c(channels.ind['CD19']), all.cuts = T, upper = F, percentile = NA, adjust.dens = 0.5))
      if(CD19.Bcell.gate < quantile(cDC@exprs[, c(channels.ind['CD19'])], c(0.5))){
        CD19.Bcell.gate <- deGate(Bcells0.flowD, channel = c(channels.ind['CD19']), use.upper = T, upper = F, tinypeak.removal = 0.9)
        #CD19.Bcell.gate <- deGate(temp.flowD, channel = c(channels.ind['CD19']), use.upper = T, upper = F, tinypeak.removal = 0.9, adjust.dens = 0.5, alpha = 0.5)
      }
    }
    
    RPmac.CD19negCD11cneg <- flowDensity(Bcells0.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(F, F), gates = c(CD19.Bcell.gate, min(cDCs.flowD@filter[, 2])))
    Bcells.flowD <- flowDensity(Bcells0.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), position = c(T, F), gates = c(CD19.Bcell.gate, max(cDCs.flowD@filter[, 2])))
    
    # Gating B cells to get B1B/B2B cells and the B2B cell subsetss FoB, MZB and Pre-B ----------------------------------------------------------
    
    # cannot use the CD5Ly6G.gate found for NK/NKT gating for BCM (tho it works really well for TCP./CIPHE)
    CD5Ly6G.Bcell.gate <- deGate(getflowFrame(Bcells.flowD), channel = c(channels.ind['CD5/Ly6G']), use.upper = T, upper = T, tinypeak.removal = 0.5)
    B2Bcells.flowD <- flowDensity(Bcells.flowD, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD19"]), position = c(F, NA), gates = c(CD5Ly6G.Bcell.gate, NA))
    

    if(!(centre == 'BCM')){
      FoB.gating.iO <- T
    }else{ # For BCM some of the B2B have strange CD21/CD35 vs CD23 distributions - I am finding the files with strange distributions and ignoring them
      # Next double check that the FoB/MZB/pre0b distribution looks as expected.
      temp.gate <- deGate(B2Bcells.flowD, channel = c(channels.ind['CD21/CD35']), use.upper = T, upper = F, alpha = 0.9, tinypeak.removal = 0.9)
      CD21neg.idx <- which(getflowFrame(B2Bcells.flowD)@exprs[, c(channels.ind['CD21/CD35'])] < (temp.gate - 0.2))
      if(length(CD21neg.idx) > 10){
        CD21pos.CD23median <- median(getflowFrame(B2Bcells.flowD)@exprs[-CD21neg.idx, c(channels.ind['CD23'])])
        CD21neg.CD23median <- median(getflowFrame(B2Bcells.flowD)@exprs[CD21neg.idx, c(channels.ind['CD23'])])
        
        # Below fails for CIPHE files such as "16-Jun-14_IMPC2_03_labelled.fcs" because the deGate above is not ok
        FoB.gating.iO <- CD21neg.CD23median < (CD21pos.CD23median - 0.35)
      }else{
        FoB.gating.iO <- F
      }
    }
    
    theta0 <- pi/6
    rot <- rotate.data(getflowFrame(B2Bcells.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = -theta0)$data
    rot.gate <- deGate(rot, channel = c(channels.ind['CD21/CD35']), use.upper = T, upper = T, tinypeak.removal = 0.5)
    rot.alt.gate <- deGate(rot, channel = c(channels.ind['CD21/CD35']), use.upper = T, upper = T, tinypeak.removal = 0.5, alpha = 0.9)
    if((rot.gate - rot.alt.gate) > 0.4){ # this suggests that the MZB peak is hidden under the FoB peak
      rot.gate <- rot.alt.gate
    }
    rot.gate0 <- deGate(rot, channel = c(channels.ind['CD21/CD35']))
    maxDens <- density(rot@exprs[, c(channels.ind['CD21/CD35'])])
    if(is.null(names(rot.gate0)) & (rot.gate0 > (maxDens$x[which.max(maxDens$y)] + 0.3))){
      rot.gate <- min(rot.gate0, rot.gate)
    }
    
    rot2 <- rotate.data(getflowFrame(B2Bcells.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = theta0)$data
    rot2.gate <- deGate(rot2, channel = c(channels.ind['CD21/CD35']), use.upper = T, upper = F, tinypeak.removal = 0.5)
    rot2.alt.gate <- deGate(rot2, channel = c(channels.ind['CD21/CD35']), use.upper = T, upper = F, tinypeak.removal = 0.5, alpha = 0.5)
    if((rot2.alt.gate - rot2.gate) > 1){ # this suggests that the MZB peak is hidden under the FoB peak
      rot2.gate <- rot2.alt.gate
    }
    # example of file the above does not work for:  "16-Mar-08_IMPC2_08_labelled.fcs"
    # The below script makes the pre-B peak stand out, if it exists.
    temp.gate <- deGate(rot2, channel = c(channels.ind['CD21/CD35']), use.percentile = T, percentile = 0.1)
    temp.flowD <- flowDensity(rot2, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(F, NA), gates = c(temp.gate, NA))
    maxDens <- density(getflowFrame(temp.flowD)@exprs[,channels.ind["CD23"]])
    temp.flowD <- flowDensity(rot2, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(NA, F), gates = c(NA, maxDens$x[which.max(maxDens$y)] + 0.2))
    temp.flowD <- flowDensity(temp.flowD, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(NA, T), gates = c(NA, maxDens$x[which.max(maxDens$y)] - 0.2))
    rot2.alt.gate <- deGate(temp.flowD, channel = c(channels.ind['CD21/CD35']))
    if(is.null(names(rot2.alt.gate))){
      maxDens <- density(rot2@exprs[, c(channels.ind['CD21/CD35'])])
      if(rot2.alt.gate < (maxDens$x[which.max(maxDens$y)] - 0.3)){
        rot2.gate <- max(rot2.gate, rot2.alt.gate)
      }
    }
    
    MZB.flowD <- flowDensity(rot, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(T, NA), gates = c(rot.gate, NA))
    rot.gate2 <- deGate(MZB.flowD, channel = c(channels.ind['CD23']), upper = T)
    MZB.flowD <- flowDensity(rot, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(T, F), gates = c(rot.gate, rot.gate2))
    MZB.flowD@filter <- rotate.data(MZB.flowD@filter, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = theta0)$data
    MZB.flowD@flow.frame <- rotate.data(getflowFrame(MZB.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]),theta = theta0)$data
    
    preB.flowD <- flowDensity(rot2, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(F, NA), gates = c(rot2.gate, NA))
    preB.flowD@filter <- rotate.data(preB.flowD@filter, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = -theta0)$data
    preB.flowD@flow.frame <- rotate.data(getflowFrame(preB.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]),theta = -theta0)$data
    
    if(centre == 'TCP'){
      theta0 <- pi/48
      rot <- rotate.data(getflowFrame(preB.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = theta0)$data
      CD23.gate <- deGate(rot, channel = c(channels.ind['CD23']), sd.threshold = T, percentile = NA, upper = NA)
      preB.flowD <- flowDensity(rot, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), position = c(NA, F), gates = c(NA, CD23.gate))
      preB.flowD@filter <- rotate.data(preB.flowD@filter, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), theta = -theta0)$data
      preB.flowD@flow.frame <- rotate.data(getflowFrame(preB.flowD), c(channels.ind["CD21/CD35"], channels.ind["CD23"]),theta = -theta0)$data
    }
    
    maxDens <- density(getflowFrame(B2Bcells.flowD)@exprs[, c(channels.ind['CD21/CD35'])])
    CD21.peak <- maxDens$x[which.max(maxDens$y)]
    maxDens <- density(getflowFrame(B2Bcells.flowD)@exprs[, c(channels.ind['CD23'])])
    CD23.peak <- maxDens$x[which.max(maxDens$y)]
    
    
    # Get NK and NKT Subsets ----------------------------------------------------------------------
    
    NKcells.SubsetQ1.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(F, T), gates = c(CD11b.gate, Ly6C.gate))
    NKcells.SubsetQ2.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(T, T), gates = c(CD11b.gate, Ly6C.gate))
    NKcells.SubsetQ3.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(F, F), gates = c(CD11b.gate, Ly6C.gate))
    NKcells.SubsetQ4.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(T, F), gates = c(CD11b.gate, Ly6C.gate))
    
    NKTcells.SubsetQ1.flowD <- flowDensity(NKTcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(F, T), gates = c(CD11b.gate, Ly6C.gate))
    NKTcells.SubsetQ3.flowD <- flowDensity(NKTcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), position = c(F, F), gates = c(CD11b.gate, Ly6C.gate))
    
    # Plots -----------------------------------------------------------------------------------------
    
    png(file = paste0(scatterplot.dir, '/', fname, '.png'), width = 2000*7/4, height = 2000)
    par(mfrow = c(4,7), mar = (c(5, 5, 4, 2) + 0.1))
    
    # Plot flowCut data
    plotDens(fT, channels = c(channels.ind["Time"], channels.ind["CD19"]), main = " All Events", cex.lab = 2, cex.axis = 2, cex.main=2)
    points(fT@exprs[fclean.ind, c(channels.ind["Time"], channels.ind["CD19"])], pch ='.')
    
    plotDens(f, c(scat.chans['FSC-A'], channels.ind["Live"]), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(live.flowD@filter)
    
    plotDens(live.flowD@flow.frame, c('FSC-A','SSC-A'), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    
    # Plot FSC singlet gating
    plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "Size", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(FSCsinglets.flowD@filter)
    
    # Plot SSC singlet gating
    #if(!is.na(scat.chans['SSC-W'])){ 
      
      # col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = primary.colors, nbin = 1000)
      # plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
      # lines(singlets.flowD@filter, col = 'blue')
      # 
    col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = matlab.like2, nbin = 1000)
    plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    lines(singlets.flowD@filter)
      
    # }else{
    #   plot(1, type="n", axes=F, xlab="", ylab="")
    #   plot(1, type="n", axes=F, xlab="", ylab="")
    # }
    
    # Panel 2 gating
    
    # if((centre == 'BCM')){
    #   plotDens(singlets0, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), main = "Singlets (prior to removing CD11b edge events)", cex.lab = 2, cex.axis = 2, cex.main=2)
    #   abline(v = max.Ly6C.value, lty = 2)
    # }
    
    # # Gating singlets to get Granulocytes
    plotDens(singlets, channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD11b"]), main = "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(granulocytes.flowD@filter)
    text(1, max(granulocytes.flowD@filter[, 2]), labels =  strcat('Granulocytes\n', toString(signif(granulocytes.flowD@cell.count/nrow(singlets)*100, 2))), cex = 2)
    
    
    # Gating NOT(Granulocytes) to get Monocytes
    plotDens(NOT.granulocytes, channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), main = "NOT(granulocytes)", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(monocytes.flowD@filter)
    text(2, max(monocytes.flowD@filter[, 2]), labels =  strcat('Monocytes\n', toString(signif(monocytes.flowD@cell.count/NOT.granulocytes@cell.count*100, 3))), cex = 2)
    abline(v = monocyte.Ly6C.target, lty = 2) 
    
    # Gating NOT(Monocytes) to get Eosinophils
    plotDens(NOT.monocytes, channels = c(channels.ind["CD11b"], scat.chans['SSC-H']), main = "NOT(monocytes)", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(eosinophils.flowD@filter)
    text(min(eosinophils.flowD@filter[, 1]) - 1, max(eosinophils.flowD@filter[, 2]) - 5000, labels =  strcat('Eosinophils\n', toString(signif(eosinophils.flowD@cell.count/NOT.monocytes@cell.count*100, 3))), cex = 2)
    abline(v = CD11b.upper.gate, lty = 2)
    
    plotDens(getflowFrame(eosinophils.flowD), channels = c(channels.ind["Ly6C"], channels.ind["CD11b"]), main = "Eosinophils",
             xlim = c(0,4.5), ylim = c(0,4.5), cex.lab = 2, cex.axis = 2, cex.main=2)
    
    if(centre == 'TCP'){
      plotDens(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind['MHCII']), main = "NOT(eosinophils)", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(CD161pos.flowD@filter)
      text(min(CD161pos.flowD@filter[, 1]) + 1, max(CD161pos.flowD@filter[, 2]) + 0.2, labels =  strcat('CD161+\n', toString(signif(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count*100, 3))), cex = 2)
    }else{
      plotDens(NOT.eosinophils, channels = c(channels.ind["CD161"], channels.ind["CD19"]), main = "NOT(eosinophils)", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(CD161pos.flowD@filter)
      text(min(CD161pos.flowD@filter[, 1]) + 1, max(CD161pos.flowD@filter[, 2]) + 0.2, labels =  strcat('CD161+\n', toString(signif(CD161pos.flowD@cell.count/NOT.eosinophils@cell.count*100, 3))), cex = 2)
      if(centre == "Jax"){ # CD161 and CD317 channels are shared
        abline(v = CD161CD317.gate)
      }
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
    text(CD11b.gate - 0.5, Ly6C.gate + 0.5, labels = 'Q1', cex = 2)
    text(CD11b.gate + 0.5, Ly6C.gate + 0.5, labels = 'Q2', cex = 2)
    text(CD11b.gate - 0.5, Ly6C.gate - 0.5, labels = 'Q3', cex = 2)
    text(CD11b.gate + 0.5, Ly6C.gate - 0.5, labels = 'Q4', cex = 2)
    plotDens(NKTcells.flowD, channels = c(channels.ind["CD11b"], channels.ind["Ly6C"]), main = "NKT cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11b.gate)
    abline(h = Ly6C.gate)
    text(CD11b.gate - 0.5, Ly6C.gate + 0.5, labels = 'Q1', cex = 2)
    text(CD11b.gate - 0.5, Ly6C.gate - 0.5, labels = 'Q3', cex = 2)
    
    plotDens(getflowFrame(NOT.CD161pos), channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), main = "NOT(CD161+)", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(Tcells.filter)
    text(1, 4, labels = 'T cells', cex = 2)
    
    plot(Ly6Cdensity.Tcells, main = "T cells",  xlab = paste0("<", pData(parameters(f))$desc[channels.ind['Ly6C']], ">:", pData(parameters(f))$name[channels.ind['Ly6C']]),cex.lab = 2, cex.axis = 2, cex.main=2)
    polygon(Ly6Cdensity.Tcells, col = 'gray', border = 'black')
    abline(v = Ly6C.gate, lty = 2)
    text(Ly6C.gate + 1, 0.1, labels = 'T Subset', cex = 2)
    
    # if(centre == 'CIPHE'){ #pDC gates
    #   plotDens(NOT.Tcells, channels = c(channels.ind['CD317'], channels.ind['Ly6C']), main = "NOT(T cells)", cex.lab = 2, cex.axis = 2, cex.main=2)
    #   abline(v = CD317.gate)
    # }
    
    plotDens(NOT.Tcells, channels = c(channels.ind["MHCII"], channels.ind["CD5/Ly6G"]), main = "NOT(T cells)", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = MHCII.gate)
    lines(MHCIIpos.flowD@filter)
    

    if(!is.na(channels.ind['F4/80'])){
      
      plotDens(getflowFrame(MHCIIpos.flowD), channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), main = "MHCII+", cex.lab = 2, cex.axis = 2, cex.main=2)
      
      plotDens(getflowFrame(MHCIIpos.flowD), channels = c(channels.ind["F4/80"], channels.ind["MHCII"]), main = "MHCII+", cex.lab = 2, cex.axis = 2, cex.main=2)
      data.new <- exprs(RPmac)[, c(channels.ind["F4/80"], channels.ind["MHCII"])]
      points(data.new, pch = '.')
      text(max(data.new[, 1]) - 0.5, max(data.new[,2]) + 0.2, labels = 'RP mac (F4/80+)', cex = 2)
      
      plotDens(RPmac, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "RP mac", cex.lab = 2, cex.axis = 2, cex.main=2)

      plotDens(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+ - RP mac", cex.lab = 2, cex.axis = 2, cex.main=2)
      data.new <- exprs(getflowFrame(MHCIIpos.flowD))[, c(channels.ind["CD19"], channels.ind["CD11c"])]
      z <- kde2d(data.new[, 1], data.new[, 2])
      contour(z, drawlabels = FALSE, add = TRUE, lty = 2)

      lines(cDCs.flowD@filter)
      text(max(cDCs.flowD@filter[,1]) -0.5, max(cDCs.flowD@filter[,2]) + 0.25 , labels =  strcat('cDCs \n', toString(signif(cDCs.flowD@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)

    }else{

      # if(!is.na(ff_CD11CFMO)){
      #     plotDens(ff_CD11CFMO, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+ for CD11c FMO control", 
      #              ylim = c(-0.5, 4.5), cex.lab = 2, cex.axis = 2, cex.main=2)
      #     abline(h = CD11c.cDC.gate)
      # }
    
      plotDens(ff, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+", ylim = c(-0.5, 4.5), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      data.new <- exprs(getflowFrame(MHCIIpos.flowD))[, c(channels.ind["CD19"], channels.ind["CD11c"])]
      z <- kde2d(data.new[, 1], data.new[, 2])
      contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
      
      lines(cDCs.flowD@filter)
      text(max(cDCs.flowD@filter[,1]) -0.5, max(cDCs.flowD@filter[,2]) + 0.25 , labels =  strcat('cDCs \n', toString(signif(cDCs.flowD@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)
      #lines(Bcells.flowD@filter)
      lines(RPmac.CD19negCD11cneg@filter)
      text(max(RPmac.CD19negCD11cneg@filter[,1]) - 0.25, min(RPmac.CD19negCD11cneg@filter[,2]), labels =  strcat('RP mac (CD11c-CD19-)\n', toString(signif(RPmac.CD19negCD11cneg@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)
    }

    plotDens(Bcells0.flowD, channels = c(channels.ind["CD19"], channels.ind["CD11c"]), main = "MHCII+ - cDCs", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(Bcells.flowD@filter)
    text(max(Bcells.flowD@filter[,1]) -0.5, max(Bcells.flowD@filter[,2]) - 0.5 , labels =  strcat('B cells \n', toString(signif(Bcells.flowD@cell.count/nrow(na.omit(ff@exprs))*100, 3))), cex = 2)
    
    plotDens(cDCs.flowD, channels = c(channels.ind['CD11b'], channels.ind['MHCII']), main = "cDCs", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11b.cDC.gate)
    data.new <- exprs(getflowFrame(cDCs.flowD))[, c(channels.ind["CD11b"], channels.ind["MHCII"])]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
   
    plot(cDCs.CD11bdensity, main = "cDCs",  xlab = paste0("<", pData(parameters(f))$desc[channels.ind['CD11b']], ">:", pData(parameters(f))$name[channels.ind['CD11b']]),cex.lab = 2, cex.axis = 2, cex.main=2)
    polygon(cDCs.CD11bdensity, col = 'gray', border = 'black')
    #abline(v = CD11b.gate, lty = 2)
    abline(v = CD11b.cDC.gate)
    text(CD11b.cDC.gate + 1, 0.1, labels = 'cDCs CD11b Type', cex = 2)
  
    plotDens(getflowFrame(Bcells.flowD), channels = c(channels.ind["CD5/Ly6G"], channels.ind["CD19"]), main = "Bcells", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD5Ly6G.Bcell.gate)
    text(CD5Ly6G.Bcell.gate - 0.5, 3, labels = 'B2B cells', cex = 2)
    text(CD5Ly6G.Bcell.gate + 0.75, 3, labels = 'B1B cells', cex = 2)
    
    plotDens(getflowFrame(B2Bcells.flowD), channels = c(channels.ind["CD21/CD35"], channels.ind["CD23"]), main = "B2Bcells", cex.lab = 2, cex.axis = 2, cex.main=2)
    data.new <- exprs(getflowFrame(B2Bcells.flowD))[, c(channels.ind["CD21/CD35"], channels.ind["CD23"])]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
    
    
    if(FoB.gating.iO){
      
      plotDens(getflowFrame(B2Bcells.flowD), channels = c(channels.ind["CD21/CD35"], channels.ind["CD23"]), main = "B2Bcells",
               cex.lab = 2, cex.axis = 2, cex.main=2)
      points(getflowFrame(MZB.flowD)@exprs[, c(channels.ind["CD21/CD35"], channels.ind["CD23"])], pch = '.')
      text(CD21.peak + 1, CD23.peak, labels = 'MZB', cex = 2)
      points(getflowFrame(preB.flowD)@exprs[, c(channels.ind["CD21/CD35"], channels.ind["CD23"])], pch = '.', col = 'red')
      text(CD21.peak - 1.5, CD23.peak + 0.3, labels = 'pre-B', cex = 2)

      text(CD21.peak, CD23.peak + 1, labels = 'Follicular B', cex = 2)
      
      # if(centre == 'TCP'){
      #   plotDens(preB.flowD, c(channels.ind["CD21/CD35"], channels.ind["CD23"]), main = "pre-B cells",
      #            cex.lab = 2, cex.axis = 2, cex.main=2)
      #   maxDens <- density(getflowFrame(preB.flowD)@exprs[,c(channels.ind['CD21/CD35'])])
      #   plot(maxDens, main = "pre-B cells",  xlab = paste0("<", pData(parameters(f))$desc[channels.ind['CD21/CD35']], ">:", pData(parameters(f))$name[channels.ind['CD21/CD35']]),cex.lab = 2, cex.axis = 2, cex.main=2)
      #   polygon(maxDens, col = 'gray', border = 'black')
      #   
      # }

    }
    
    dev.off()

    # Save IMPReSS populations----------------------------------------------------
    
    immpressData.FCSfile[1:5] <- c(original.no.cells, nrow(singlets)/original.no.cells*100, granulocytes.flowD@cell.count, monocytes.flowD@cell.count, 
                                   eosinophils.flowD@cell.count)
    immpressData.FCSfile[6:10] <- c(NKcells.flowD@cell.count, NKcells.SubsetQ1.flowD@cell.count, NKcells.SubsetQ2.flowD@cell.count, 
                                    NKcells.SubsetQ3.flowD@cell.count, NKcells.SubsetQ4.flowD@cell.count)
    immpressData.FCSfile[11:13] <- c(NKTcells.flowD@cell.count, NKTcells.SubsetQ1.flowD@cell.count, NKcells.SubsetQ3.flowD@cell.count)
    immpressData.FCSfile[14:15] <- c(nrow(Tcells), Ly6CposTcells.flowD@cell.count) 
    immpressData.FCSfile[16:18] <- c(Bcells.flowD@cell.count, Bcells.flowD@cell.count - B2Bcells.flowD@cell.count, B2Bcells.flowD@cell.count)
    
    immpressData.FCSfile[25:26] <- c(cDCs.flowD@cell.count, cDCs.CD11bpos.flowD@cell.count)
    
    if(FoB.gating.iO){
      immpressData.FCSfile[21] <- c(preB.flowD@cell.count)
      immpressData.FCSfile[23] <- c(MZB.flowD@cell.count)
      immpressData.FCSfile[19] <- c(B2Bcells.flowD@cell.count - preB.flowD@cell.count - MZB.flowD@cell.count)
    }
    
    # if(centre == 'CIPHE'){
    #   immpressData.FCSfile[27] <- c(pDC.cellcount)
    # }
    
    if(!is.na(channels.ind['F4/80'])){
      immpressData.FCSfile[28] <- nrow(RPmac)
    }else{
      immpressData.FCSfile[29] <- RPmac.CD19negCD11cneg@cell.count 
    }
  })
  
  data.frame(x$FCS.files, passed.flowCut = passed.flowCut, immpressData.FCSfile)
    
}, .parallel = TRUE) # end ldply

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv") 
save(immpressData, file = paste0('immpressData', centre, '.RData'))

colnames(immpressData)[1] <- "FCS.files"
IMMPReSS.data <- join(file.names, immpressData)
cols2remove <- c('passed.flowCut', 'Path', 'Panel.Organ.Folder')
IMMPReSS.data <- IMMPReSS.data[, -c(which(colnames(IMMPReSS.data) %in% cols2remove))]
write.csv(IMMPReSS.data, file =  paste('/data/IMPC/DCCResults_', centre, '_Panel2', date.time, sep=""), row.names = F)

