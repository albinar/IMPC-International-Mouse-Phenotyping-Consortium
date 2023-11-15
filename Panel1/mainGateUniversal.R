## Developed by Albina Rahim

remove(list=ls())

## Setting the path based on docker
setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")

## This version of the script follows the new gating strategy given to us on September 2018.
## Datasets analysed will be from TCP, UCD, BCM, and Harwell.
## I am using the qualityGate function written by Sibyl for obtaining the Singlets population from Ungated population.

###############################################################################################
## NOTE: Run this section functions either separately before you run the rest of the script or you can run the entire script in just one go
## source("Codes/preProcessingMain.R"). The following two function prompts user for input: one for the centre name and the other panel number

source("helperFunc.R")
source("qualityGate.R") 
if (interactive() ){
  centre <- readCentreFunc()
  if(centre != "sanger"){
    if (interactive() ){
      panel <- readPanelFunc()
    }
  }
}



library('colorRamps')
library('plyr')
library('doMC')
library('e1071') # for flowCut
library('Cairo') # for flowCut
library('flowCore')
library('flowDensity')
library('pracma') # for findpeaks
library('tools')
library('MASS')
library('stringr')## for str_match used in compensateIMPC function
library('flowViz')
library('knitr')
library("flowCut")



start <- Sys.time()



results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")




load(paste0(results.dir,"/lgl.Rdata"))
load(paste0(results.dir,"/Genotype.Rdata"))
load(paste0(results.dir,"/uniqueGT.Rdata"))
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/store.allFMO.Rdata"))
load(paste0(results.dir,"/Props.Events.Gates.Rdata"))
load(paste0(results.dir,"/Gates.Filter.list.Rdata"))
load(paste0(results.dir,"/all.gthres.store.Rdata"))
## Temporary line for working on the failed files
failedFCS <- read.csv(paste0(results.dir,"/Failed_Gating_FilesTCP_Panel1_20220112_0156.csv"))


# ## Loading th Impress ID Table for naming the populations accordingly
load('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/immpressIDconversion.Rdata')
colnames(immpressIDconversion) <- c('Impress IDs', 'Population Name')

# Create directories
suppressWarnings(dir.create(paste0(results.dir,"/Live-Singlets-Figures")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlots/")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/FailedFiles/")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/flowCut/")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/flowCut-FMO/")))
suppressWarnings(dir.create(paste0(results.dir,"/SingletsRdata")))
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create(paste0(results.dir,"/SingletsRdata/", uniqueGT[x])))
}))



rownames(store.allFCS) <- 1:nrow(store.allFCS)

store.allFCS.Original <- store.allFCS
store.allFCS <- store.allFCS[1266:nrow(store.allFCS),]
rownames(store.allFCS) <- 1:nrow(store.allFCS)

## temporary 

#file.names <- data.frame(failedFCS, stringsAsFactors = F)

#file.names <- data.frame(store.allFCS[2767:2775,], stringsAsFactors = F)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)


FMO.cd5.flag <- 0
FMO.cd161.flag <- 0 
FMO.cd44.flag <- 0
FMO.cd25.flag <- 0
FMO.cd62l.flag <- 0
FMO.tcrd.flag <- 0
FMO.cd45.flag <- 0
FMO.cd4.flag <- 0
FMO.cd8a.flag <- 0

indexAssayDate <- 0
indexGenotype <- 0


no_cores <- detectCores() - 4
registerDoMC(no_cores)

print("Starting Gating & Plotting")
props.events.gates <- llply(1:nrow(file.names), function(i){ 
  #props.events.gates <- llply(1:length(index), function(i){ 
  #props.events.gates <- llply(1:20, function(i){ 
  
  x<- file.names[i,]
  #x<- file.names[index[i],]
  
  all.events <- matrix(nrow = 1, ncol = 35, data = NA)# matrix for saving the event counts
  all.props <- matrix(nrow = 1, ncol = 35, data = NA)# matrix for saving the proportions
  all.gthres <- matrix(nrow = 1, ncol = 23, data = NA) # matrix for saving the gating thresholds
  all.fClean <-  matrix(nrow = 1, ncol = 1, data = NA) # matrix for saving the outputs from flowCut
  flaggedFile <- matrix(nrow = 1, ncol = 1, data = NA)
  ## Saving the flowDensity filter for flowType
  Filters.list <- list("NKcells.filter" = 0, "NK.Effector.filter" =0 , "NK.Resting.filter" =0, "cd5.filter" =0, "cd161pos.cd8neg.filter" = 0, "cd4pos.NKT.filter" = 0, "cd4neg.NKT.filter" = 0,
                       "cd4pos.NKT.Effector.filter" = 0, "cd4pos.NKT.Resting.filter" = 0, "cd4neg.NKT.Effector.filter" =0, "cd4neg.NKT.Resting.filter" =0, "NOT.Tcells.filter" =0, "cd4.Tcells.filter" =0,
                       "cd8.Tcells.filter" =0, "Tregs.filter" =0, "Effector.Tregs.filter" =0, "Resting.Tregs.filter" =0, "Effector.T.helper.filter" =0, "Resting.T.helper.filter" =0,
                       "Effector.cd8.filter" =0, "Resting.cd8.filter" =0, "Naive.cd8.filter" =0)
  
  
  tryCatch({
    if(centre == "sanger" | centre == "ciphe"){
      f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
      fpath <- x$Path
    }else if(centre == "tcp" | centre == "bcm"| centre == "jax" | centre == "gmc" | centre == "ucd" | centre == "ccp"){
      f <- try(read.FCS(filename = paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)), silent = TRUE)
      fpath <- paste0(x$Path,"/", x$Panel.Organ.Folder)
      if(centre == "gmc"){
        store.allFMO <- cbind(store.allFMO, NA)
        colnames(store.allFMO) <- c("Path", "Panel/Organ/Folder", "FMO", "FMO-Names")
        store.allFMO[1,c('FMO-Names')] <- "Unstain"
        store.allFMO[2,c('FMO-Names')] <- "TCRD"
        store.allFMO[3,c('FMO-Names')] <- "CD161"
        store.allFMO[4,c('FMO-Names')] <- "CD4"
        store.allFMO[5,c('FMO-Names')] <- "CD62L"
        store.allFMO[6,c('FMO-Names')] <- "CD25"
        store.allFMO[7,c('FMO-Names')] <- "CD45"
        store.allFMO[8,c('FMO-Names')] <- "CD8A"
        store.allFMO[9,c('FMO-Names')] <- "CD5"
        store.allFMO[10,c('FMO-Names')] <- "CD44"
        
      }else if(centre != "ccp"){
        if(centre == "ucd"){
          FMO.index <- which(x$Panel.Organ.Folder == store.allFMO[,c("Panel/Organ/Folder")])
        }else{
          FMO.index <- which(x$Assay.Date == store.allFMO[,c("Assay Date")])
        }
      }
    }
    
    ## For CCP we are removing the SPILLOVER Matrix and replacing it with the Compensation matrix inside each sub-folder
    if(centre == "ccp"){
       compMatrix <- dir(fpath, full.names = T, recursive = F, pattern = ".*csv")
       SPILL.matrix <- read.csv(compMatrix, stringsAsFactors = F, check.names = F)[,-1]
       # #colnames(SPILL.matrix) <- f.colnames
       # f@description$`$SPILLOVER` <- NULL
       # f@description$`$SPILLOVER` <- SPILL.matrix
     
       if(length(f@description$SPILL) > 0){
         colnames(SPILL.matrix) <- colnames(f@description$SPILL)
         f@description$SPILL <- NULL
         f@description$SPILL <- SPILL.matrix
       }else{
         colnames(SPILL.matrix) <- colnames(f@description$`$SPILLOVER`)
         f@description$`$SPILLOVER` <- NULL
         f@description$`$SPILLOVER` <- SPILL.matrix
       }
     }

    
    # totalEvents <- nrow(f@exprs)
    # impressID.populations[1] <- totalEvents ## Total number of acquired events in Panel A
    
    ## UCDavis has some files in panel 2 in which the marker labels are missing.
    if(centre == "ucd" & panel == 2){
      if(is.na(f@parameters@data$desc[7])){
        f@parameters@data$desc <- markers.temp
      }
      markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5", "Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                   "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "gdTCR|TCR*|Tcr-d", "CD45", "CD43")
    }else{
      markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5|Ly6G|CD5V450+Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                   "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD4", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "CD3e|CD3", "gdTCR|TCR*|Tcr-d", "CD45", "CD43", "via")
    }
    
    channels.ind <-Find.markers(f, markers)
    names(channels.ind)[grep(names(channels.ind), pattern = "Live*")] <- "Live"
    names(channels.ind)[grep(names(channels.ind), pattern = "*II")] <- "MHCII"
    if(centre != "ucd"){
      names(channels.ind)[grep(names(channels.ind), pattern = "*Ly6G")] <- "CD5/Ly6G"
    }
    
    names(channels.ind)[grep(names(channels.ind), pattern = "*6C")] <- "Ly6C"
    names(channels.ind)[grep(names(channels.ind), pattern = "F4")] <- "F4/80"
    names(channels.ind)[grep(names(channels.ind), pattern = "CD21")] <- "CD21/CD35"
    names(channels.ind)[grep(names(channels.ind), pattern = "CD62*")] <- "CD62L"
    names(channels.ind)[grep(names(channels.ind), pattern = "KLRG*|Klrg")] <- "KLRG1"
    names(channels.ind)[grep(names(channels.ind), pattern = "gdTCR|TCR*|Tcr*")] <- "TCRD"
    names(channels.ind)[grep(names(channels.ind), pattern = "CD11c|CD11C")] <- "CD11c"
    names(channels.ind)[grep(names(channels.ind), pattern = "via")] <- "Indo-1 Violet" ## Marker used by CCP instead of Live
    
    if(panel == 1){
      names(channels.ind)[grep(names(channels.ind), pattern = "CD5|Ly6G")] <- "CD5"
    }
    
    if(is.na(channels.ind["Live"])){  
      # BCM (new data, DAPI is under the $name descriptor 
      if(length(grep('DAPI', pData(parameters(f))$name)) > 0){
        channels.ind["Live"] <- grep('DAPI', pData(parameters(f))$name)
      }else if(length(grep('TCRD', pData(parameters(f))$desc)) > 0){
        # For TCP, Sytox Blue is read in the BV-510 channel. TCRD is under $desc for TCP
        channels.ind["TCRD"] <- NA
        channels.ind <- channels.ind[!is.na(channels.ind)]
        channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
      }else if(length(grep('Sytox Blue', pData(parameters(f))$desc)) > 0){
        # For TCP, Sytox Blue is read in the BV-510, V525-A, V480-A channels.
        channels.ind["Sytox Blue"] <- NA
        channels.ind <- channels.ind[!is.na(channels.ind)]
        if(length(grep('BV510', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
        }
        if(length(grep('V525', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
        }
        if(length(grep('V480', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
        }
      }else if(length(grep('SYTOX BLUE', pData(parameters(f))$desc)) > 0){
        # For TCP, SYTOX BLUE is read in the BV-510, V525-A, V480-A channels.
        channels.ind["SYTOX BLUE"] <- NA
        channels.ind <- channels.ind[!is.na(channels.ind)]
        if(length(grep('BV510', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
        }
        if(length(grep('V525', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
        }
        if(length(grep('V480', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
        }
        
      }else if(length(grep('via', pData(parameters(f))$desc)) > 0){ 
        ## CCP doesn't have Live marker. They want to use Indo-1 (Violet)-A marker in place of Live 
        channels.ind["Indo-1 Violet"] <- NA
        channels.ind <- channels.ind[!is.na(channels.ind)]
        channels.ind["Live"] <- grep('Indo*', pData(parameters(f))$name)
      }else{
        # For TCP, Sytox Blue is read in the BV510 & V525-A channels
        if(length(grep('V525', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
        }
        if(length(grep('BV510', pData(parameters(f))$name)) > 0){
          channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
        }
        
      }
    }
    
    if(is.na(channels.ind["Time"])){  
      channels.ind["Time"] <- grep('TIME|Time', pData(parameters(f))$name)
    }
    
    channels.ind <- sort(channels.ind)
    
    
    scat.chans <- c(grep(colnames(f),pattern = "FS-*|FSC*"), grep(colnames(f),pattern = "SS-*|SSC*"))
    names(scat.chans) <- colnames(f)[scat.chans]
    
    ## This part is added because of GMC
    names(scat.chans)[grep(names(scat.chans), pattern = "FS-H")] <- "FSC-H"
    names(scat.chans)[grep(names(scat.chans), pattern = "FS-A")] <- "FSC-A"
    names(scat.chans)[grep(names(scat.chans), pattern = "FS-W")] <- "FSC-W"
    names(scat.chans)[grep(names(scat.chans), pattern = "SS-H")] <- "SSC-H"
    names(scat.chans)[grep(names(scat.chans), pattern = "SS-A")] <- "SSC-A"
    names(scat.chans)[grep(names(scat.chans), pattern = "SS-W")] <- "SSC-W"
    
    # Remove scatter margins and compensate ---------------------------------------------
    # Removing margin events in Scatter channels
    f <- removeMargins(f, chans = scat.chans, verbose = F)
    #Removing negative values in scatter channels
    f <- removeMargins(f, chans = scat.chans, debris = T, neg = T, verbose = F)
    
    #f <- compensateIMPC(f, basename(x$FCS.files), fpath, centre = centre, panel.no = panel)
    # # Adding the Assay Date to the compensateIMPC function
    # f <- compensateIMPC(f, basename(x$FCS.files), fpath, centre = centre, panel.no = panel, assayDate = x$Assay.Date)
    
    
    load(paste0(results.dir,"/lgl.Rdata"))
    lgl.name <- names(lgl@transforms)
    
    if(centre == "ucd"){
      indexGenotype<-which(store.allFCS[,c('Genotype')]==store.allFCS[i,c('Genotype')])
      
      ## Creating gFrame based on Genotype
      store.allFCS.Genotype <- store.allFCS[indexGenotype,]
      file.names.Genotype <- data.frame(store.allFCS.Genotype, stringsAsFactors = F)
      
      
      gFrameGenotype <- ddply(file.names.Genotype, "FCS.files", function(y){
        
        ##Uncomment this part when testing individual file
        #y<- file.names.Genotype[j,]
        
        if(centre == "sanger" | centre == "ciphe"){
          f <- try(read.FCS(filename = paste0(y$Path, "/", y$FCS.files)), silent = TRUE)
          fpath <- y$Path
        }else if(centre == "tcp" | centre == "bcm" | centre == "jax" | centre == "gmc" | centre == "ucd"){
          f <- try(read.FCS(filename = paste0(y$Path, "/", y$Panel.Organ.Folder, "/", y$FCS.files)), silent = TRUE)
          fpath <- paste0(y$Path,"/", y$Panel.Organ.Folder)
        }
        
        ## Adding the Assay Date to the compensateIMPC function
        #f <- compensateIMPC(f, basename(x$FCS.files), fpath, centre = centre, panel.no = panel, assayDate = y$Assay.Date)
        
        ## UCDavis has some files in panel 2 in which the marker labels are missing.
        if(centre == "ucd" & panel == 2){
          if(is.na(f@parameters@data$desc[7])){
            f@parameters@data$desc <- markers.temp
          }
          markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5", "Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                       "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "gdTCR|TCR*|Tcr-d", "CD45", "CD43")
        }else{
          markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5|Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                       "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD4", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "CD3e|CD3", "gdTCR|TCR*|Tcr-d", "CD45", "CD43")
        }
        
        channels.ind <-Find.markers(f, markers)
        names(channels.ind)[grep(names(channels.ind), pattern = "Live*")] <- "Live"
        names(channels.ind)[grep(names(channels.ind), pattern = "*II")] <- "MHCII"
        if(centre != "ucd"){
          names(channels.ind)[grep(names(channels.ind), pattern = "*Ly6G")] <- "CD5/Ly6G"
        }
        
        names(channels.ind)[grep(names(channels.ind), pattern = "*6C")] <- "Ly6C"
        names(channels.ind)[grep(names(channels.ind), pattern = "F4")] <- "F4/80"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD21")] <- "CD21/CD35"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD62*")] <- "CD62L"
        names(channels.ind)[grep(names(channels.ind), pattern = "KLRG*")] <- "KLRG1"
        names(channels.ind)[grep(names(channels.ind), pattern = "gdTCR|TCR*|Tcr*")] <- "TCRD"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD11c|CD11C")] <- "CD11c"
        if(panel == 1){
          names(channels.ind)[grep(names(channels.ind), pattern = "CD5|Ly6G")] <- "CD5"
        }
        
        if(is.na(channels.ind["Live"])){  
          # BCM (new data, DAPI is under the $name descriptor 
          if(length(grep('DAPI', pData(parameters(f))$name)) > 0){
            channels.ind["Live"] <- grep('DAPI', pData(parameters(f))$name)
          }else if(length(grep('TCRD', pData(parameters(f))$desc)) > 0){
            # For TCP, Sytox Blue is read in the BV-510 channel. TCRD is under $desc for TCP
            channels.ind["TCRD"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
          }else if(length(grep('Sytox Blue', pData(parameters(f))$desc)) > 0){
            # For TCP, Sytox Blue is read in the BV-510, V525-A, V480-A channels.
            channels.ind["Sytox Blue"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('V480', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
            }
          }else if(length(grep('SYTOX BLUE', pData(parameters(f))$desc)) > 0){
            # For TCP, SYTOX BLUE is read in the BV-510, V525-A, V480-A channels.
            channels.ind["SYTOX BLUE"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('V480', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
            }
            
          }else{
            # For TCP, Sytox Blue is read in the BV510 & V525-A channels
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            
          }
        }
        
        if(is.na(channels.ind["Time"])){  
          channels.ind["Time"] <- grep('TIME|Time', pData(parameters(f))$name)
        }
        
        channels.ind <- sort(channels.ind)
        
        temp <- f@exprs[sample(1:length(f@exprs[,1]), 1000), as.numeric(channels.ind)]
        
        return(data.frame(temp, check.names = F))
        
      }, .parallel = TRUE) # end ddply
      
      ## Remove rows in gFrame with NAs and also saving files which has NA in their flurophore/channels
      fileswNAs <- NULL
      colswNAs <- which(is.na(gFrameGenotype), arr.ind = TRUE)
      if(length(colswNAs) > 0){
        rowswNAs <- unique(colswNAs[,'row'])
        colswNAs <- unique(colswNAs[,'col'])
        fileswNAs <- unique(gFrameGenotype[rowswNAs,1])
        gFrameGenotype <- gFrameGenotype[-rowswNAs,]
      }
      
      ## Reading the first FCS file in the storage matrix as a template for creating the global frame
      if(centre == "sanger" | centre == "ciphe" ){
        g <- read.FCS(filename = paste0(store.allFCS.Genotype[1,c('Path')], "/", store.allFCS.Genotype[1,c('FCS files')]))
      }else if(centre == "tcp" | centre == "gmc" | centre == "ucd"){ ## TCP, GMC, and UCD have inconsistency in channel number. So we are taking the one which has the maximum channel number
        channelNum <- which(store.allFCS.Genotype[,c('Number of Channels')] == max(store.allFCS.Genotype[,c('Number of Channels')]))[1]
        g <- read.FCS(filename = paste0(store.allFCS.Genotype[channelNum,c('Path')], "/", store.allFCS.Genotype[channelNum,c('Panel/Organ/Folder')],"/", store.allFCS.Genotype[channelNum,c('FCS files')]))
      }else if(centre == "bcm" | centre == "jax"){
        g <- read.FCS(filename = paste0(store.allFCS.Genotype[1,c('Path')], "/", store.allFCS.Genotype[1,c('Panel/Organ/Folder')],"/", store.allFCS.Genotype[1,c('FCS files')]))
      }
      
      Transform.idx <- unlist(sapply(colnames(gFrameGenotype), function(x) {grep(x, colnames(g))})) 
      gexprs.temp <- matrix(0, nrow = nrow(gFrameGenotype), ncol = ncol(g@exprs))
      gexprs.temp[, Transform.idx] <- as.matrix(gFrameGenotype[, 2:ncol(gFrameGenotype)])
      g@exprs <- gexprs.temp
      colnames(g@exprs) <- colnames(g) 
      
      lgl <- estimateLogicle(g, channels = colnames(g)[Transform.idx])
      fT <- transform(f, lgl)
    }else{
      fT <- transform(f, lgl)
    }
    
    
    if (!all(lgl.name %in% colnames(f)) && centre != "ucd"){
      #if(centre == "jax"){
      indexAssayDate<-which(store.allFCS[,c('Panel/Organ/Folder')]==store.allFCS[i,c('Panel/Organ/Folder')])
      
      ## Creating gFrame based on Assay Date
      store.allFCS.AssayDate <- store.allFCS[indexAssayDate,]
      file.names.AssayDate <- data.frame(store.allFCS.AssayDate, stringsAsFactors = F)
      
      gFrameAssayDate <- ddply(file.names.AssayDate, "FCS.files", function(x){
        x<- file.names.AssayDate[i,]
        
        if(centre == "sanger" | centre == "ciphe"){
          f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
          fpath <- x$Path
        }else if(centre == "tcp" | centre == "bcm" | centre == "jax" | centre == "gmc" | centre == "ucd"){
          f <- try(read.FCS(filename = paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)), silent = TRUE)
          fpath <- paste0(x$Path,"/", x$Panel.Organ.Folder)
        }
        
        ## Adding the Assay Date to the compensateIMPC function
        f <- compensateIMPC(f, basename(x$FCS.files), fpath, centre = centre, panel.no = panel, assayDate = x$Assay.Date)
        
        ## UCDavis has some files in panel 2 in which the marker labels are missing.
        if(centre == "ucd" & panel == 2){
          if(is.na(f@parameters@data$desc[7])){
            f@parameters@data$desc <- markers.temp
          }
          markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5", "Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                       "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "gdTCR|TCR*|Tcr-d", "CD45", "CD43")
        }else{
          markers <- c("Live|I515*|Syto*|SYTO*|DAPI|L/D", "CD5|Ly6G", "CD11b", "Ly6C|Lyc6C*", "CD161", "CD19", "*II|IA/E", "CD317", "F4|F4/80",
                       "CD11c|CD11C", "CD21|CD21/CD35", "CD23", "CD4", "CD62*", "CD25", "CD44", "CD8", "GITR", "KLRG*", "CD24", "CD3e|CD3", "gdTCR|TCR*|Tcr-d", "CD45", "CD43")
        }
        
        channels.ind <-Find.markers(f, markers)
        names(channels.ind)[grep(names(channels.ind), pattern = "Live*")] <- "Live"
        names(channels.ind)[grep(names(channels.ind), pattern = "*II")] <- "MHCII"
        if(centre != "ucd"){
          names(channels.ind)[grep(names(channels.ind), pattern = "*Ly6G")] <- "CD5/Ly6G"
        }
        
        names(channels.ind)[grep(names(channels.ind), pattern = "*6C")] <- "Ly6C"
        names(channels.ind)[grep(names(channels.ind), pattern = "F4")] <- "F4/80"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD21")] <- "CD21/CD35"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD62*")] <- "CD62L"
        names(channels.ind)[grep(names(channels.ind), pattern = "KLRG*")] <- "KLRG1"
        names(channels.ind)[grep(names(channels.ind), pattern = "gdTCR|TCR*|Tcr*")] <- "TCRD"
        names(channels.ind)[grep(names(channels.ind), pattern = "CD11c|CD11C")] <- "CD11c"
        if(panel == 1){
          names(channels.ind)[grep(names(channels.ind), pattern = "CD5|Ly6G")] <- "CD5"
        }
        
        if(is.na(channels.ind["Live"])){  
          # BCM (new data, DAPI is under the $name descriptor 
          if(length(grep('DAPI', pData(parameters(f))$name)) > 0){
            channels.ind["Live"] <- grep('DAPI', pData(parameters(f))$name)
          }else if(length(grep('TCRD', pData(parameters(f))$desc)) > 0){
            # For TCP, Sytox Blue is read in the BV-510 channel. TCRD is under $desc for TCP
            channels.ind["TCRD"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
          }else if(length(grep('Sytox Blue', pData(parameters(f))$desc)) > 0){
            # For TCP, Sytox Blue is read in the BV-510, V525-A, V480-A channels.
            channels.ind["Sytox Blue"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('V480', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
            }
          }else if(length(grep('SYTOX BLUE', pData(parameters(f))$desc)) > 0){
            # For TCP, SYTOX BLUE is read in the BV-510, V525-A, V480-A channels.
            channels.ind["SYTOX BLUE"] <- NA
            channels.ind <- channels.ind[!is.na(channels.ind)]
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('V480', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V480', pData(parameters(f))$name)
            }
            
          }else{
            # For TCP, Sytox Blue is read in the BV510 & V525-A channels
            if(length(grep('V525', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('V525', pData(parameters(f))$name)
            }
            if(length(grep('BV510', pData(parameters(f))$name)) > 0){
              channels.ind["Live"] <- grep('BV510', pData(parameters(f))$name)
            }
            
          }
        }
        
        if(is.na(channels.ind["Time"])){  
          channels.ind["Time"] <- grep('TIME|Time', pData(parameters(f))$name)
        }
        
        channels.ind <- sort(channels.ind)
        
        temp <- f@exprs[sample(1:length(f@exprs[,1]), 1000), as.numeric(channels.ind)]
        
        return(data.frame(temp, check.names = F))
        
      }, .parallel = TRUE) # end ddply
      
      
      ## Remove rows in gFrame with NAs and also saving files which has NA in their flurophore/channels
      fileswNAs <- NULL
      colswNAs <- which(is.na(gFrameAssayDate), arr.ind = TRUE)
      if(length(colswNAs) > 0){
        rowswNAs <- unique(colswNAs[,'row'])
        colswNAs <- unique(colswNAs[,'col'])
        fileswNAs <- unique(gFrameAssayDate[rowswNAs,1])
        gFrameAssayDate <- gFrameAssayDate[-rowswNAs,]
      }
      
      ## Reading the first FCS file in the storage matrix as a template for creating the global frame
      if(centre == "sanger" | centre == "ciphe" ){
        g <- read.FCS(filename = paste0(store.allFCS.AssayDate[1,c('Path')], "/", store.allFCS.AssayDate[1,c('FCS files')]))
      }else if(centre == "tcp" | centre == "gmc" | centre == "ucd"){ ## TCP, GMC, and UCD have inconsistency in channel number. So we are taking the one which has the maximum channel number
        channelNum <- which(store.allFCS.AssayDate[,c('Number of Channels')] == max(store.allFCS.AssayDate[,c('Number of Channels')]))[1]
        g <- read.FCS(filename = paste0(store.allFCS.AssayDate[channelNum,c('Path')], "/", store.allFCS.AssayDate[channelNum,c('Panel/Organ/Folder')],"/", store.allFCS.AssayDate[channelNum,c('FCS files')]))
      }else if(centre == "bcm" | centre == "jax"){
        g <- read.FCS(filename = paste0(store.allFCS.AssayDate[1,c('Path')], "/", store.allFCS.AssayDate[1,c('Panel/Organ/Folder')],"/", store.allFCS.AssayDate[1,c('FCS files')]))
      }
      
      Transform.idx <- unlist(sapply(colnames(gFrameAssayDate), function(x) {grep(x, colnames(g))})) 
      gexprs.temp <- matrix(0, nrow = nrow(gFrameAssayDate), ncol = ncol(g@exprs))
      gexprs.temp[, Transform.idx] <- as.matrix(gFrameAssayDate[, 2:ncol(gFrameAssayDate)])
      g@exprs <- gexprs.temp
      colnames(g@exprs) <- colnames(g) 
      
      lgl <- estimateLogicle(g, channels = colnames(g)[Transform.idx])
      fT <- transform(f, lgl)
    }else{
      fT <- transform(f, lgl)
    }
    
    
    #channels.to.clean <- setdiff(channels.ind, channels.ind['Time'])
    channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
    
    f.Clean <- flowCut(fT, Segment = floor(nrow(f)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut/"), Plot = 'Flagged Only')
    #f.Clean <- flowCut(fT, Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut/"), Plot = 'Flagged Only') #Failed Only')
    all.fClean[1] <- f.Clean$data['Has the file passed',1]
    
    f <- transform(f, lgl)
    
    if(length(f.Clean$ind) > 0){
      f@exprs <- f@exprs[-f.Clean$ind, ]
    }
    
    ## Cell counts/proportions after cleaning and margin event removal
    all.events[1] <- nrow(f)
    all.props[1] <- (nrow(f)/nrow(f))*100
    
    #######################################################################################################
    
    # Quality Gate (live/dead, singlets gating) ---------------------------------------------------
    results <- qualityGate(f, scat.chans, channels.ind, centre, panel.no = panel)
    scat.chans <- results$scat.chans
    
    live.flowD <- results$live
    live <- live.flowD@flow.frame
    all.events[2] <- live.flowD@cell.count
    all.props[2] <- live.flowD@proportion 
    
    size.flowD <- results$size
    FSCsinglets.flowD <- results$FSCsinglets
    all.events[3] <- FSCsinglets.flowD@cell.count
    all.props[3] <- FSCsinglets.flowD@proportion 
    
    singlets.flowD <- results$singlets
    all.events[4] <- singlets.flowD@cell.count
    all.props[4] <- singlets.flowD@proportion 
    singlets <- getflowFrame(singlets.flowD)
    singlets.temp <- singlets.flowD@flow.frame
    
    
    #gateAssignments <- 1:nrow(singlets)%in%live.flowD@index
    
    ## Saving Singlets flowFrame for running flowType later
    save(singlets, file = paste0(results.dir, "/SingletsRdata/",x$Genotype, "/", x$FCS.files, ".Rdata"))
    print("pre-gating finished")
    #impressID.populations[2] <-  (singlets.flowD@cell.count/totalEvents)*100 ## Percentage of live gated events in Panel A
    
    ## Saving the gating thresholds from the qualityGate.R
    all.gthres[1] <- results$fsca.live.gate
    all.gthres[2] <- results$live.gate
    all.gthres[3] <- results$ss.high
    
    plotDens(f, c(scat.chans['FSC-A'], channels.ind["Live"]), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h=results$dead.peak, lwd=2)
    lines(live.flowD@filter)
    text(max(live.flowD@filter[, 1])-100000, max(live.flowD@filter[, 2]), labels =  strcat('Live \n', toString(signif(live.flowD@proportion, 3))), cex = 1.5)
    
    plotDens(f, c(scat.chans['FSC-A'], channels.ind["Live"]), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h=results$dead.peak, lwd=2)
    points(exprs(f[singlets.flowD@index, c(scat.chans['FSC-A'], channels.ind["Live"])]), col=2, pch=".")
    
    
    plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'],scat.chans['SSC-A']), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    points(exprs(live[singlets.flowD@index, c(scat.chans['FSC-A'], scat.chans['SSC-A'])]), col=2, pch=".")
    
    # Plot size gate
    plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'],scat.chans['SSC-A']), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    text(max(f@exprs[, c(scat.chans['FSC-A'])])-100000, max(f@exprs[, c(scat.chans['SSC-A'])])-100000, labels =  strcat('Size \n 100.0'), cex = 1.5)
    
    
    #Plot FSC singlet gate
    plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(FSCsinglets.flowD@filter)
    text(max(live.flowD@filter[, 1])-100000, max(FSCsinglets.flowD@filter[, 2])-10000, labels =  strcat('FSC Singlets\n', toString(signif(FSCsinglets.flowD@proportion, 3))), cex = 1.5)
    #
    #
    # # Plot FSC singlet gate
    # #if(!is.na(scat.chans['SSC-W'])){
    #
    # # col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = primary.colors, nbin = 1000)
    # # plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    # # lines(singlets.flowD@filter, col = 'blue')
    if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
      col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-A'], scat.chans['FSC-H'])], colramp = matlab.like2, nbin = 1000)
      plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-A'], scat.chans['FSC-H']), main = "SSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    }else{
      col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = matlab.like2, nbin = 1000)
      plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "SSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    }
    lines(singlets.flowD@filter)
    text(max(FSCsinglets.flowD@filter[, 1])-200000, max(FSCsinglets.flowD@filter[, 2])-20000, labels =  strcat('SSC Singlets\n', toString(signif(singlets.flowD@proportion, 3))), cex = 1.5)
    
    # #plotDens(FSC.singlets@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    # #lines(temp2@filter)
    
    #}
    
    
    #################################################################################################################
    #################################################################################################################
    
    ## Pre-processing the CD5 FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        #FMO.index <- which(store.allFMO[,3] == "CD5")
        FMO.cd5.index <- grep("CD5", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd5.index <- grep("CD5", store.allFMO[FMO.index,3], value = FALSE)
      }
      
      
      if(length(FMO.cd5.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          #FMO.index <- which(store.allFMO[,3] == "CD5")
          FMO.cd5 <- read.FCS(filename = paste0(store.allFMO[FMO.cd5.index,"Path"], "/", store.allFMO[FMO.cd5.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd5.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd5.index <- FMO.index[FMO.cd5.index]
          }
          FMO.cd5 <- read.FCS(filename = paste0(store.allFMO[FMO.cd5.index,"Path"], "/", store.allFMO[FMO.cd5.index,"Assay Date"], "/", store.allFMO[FMO.cd5.index,"FMO"]))
          
        }
        
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd5 <- removeMargins(FMO.cd5, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd5 <- removeMargins(FMO.cd5, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd5 <- compensateIMPC(FMO.cd5, basename(store.allFMO[FMO.cd5.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd5.fT <- transform(FMO.cd5, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd5@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd5.index) <- 0
        }else{
          FMO.cd5.Clean <- flowCut(FMO.cd5.fT, Segment = floor(nrow(FMO.cd5)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          #FMO.cd5.Clean <- flowCut(FMO.cd5.fT, Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only') #Failed Only')
          # FMO.cd5.Clean$data['Has the file passed',1]
          FMO.cd5 <- transform(FMO.cd5, lgl)
          
          if(length(FMO.cd5.Clean$ind) > 0){
            FMO.cd5@exprs <- FMO.cd5@exprs[-FMO.cd5.Clean$ind, ]
          }
          remove(FMO.cd5.Clean)
          gc()
        }
      }
      
    }
    
    #############################################################################################
    ## Pre-processing the CD161 FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd161.index <- grep("CD161", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd161.index <- grep("CD161", store.allFMO[FMO.index,3], value = FALSE)
      }
      
      
      if(length(FMO.cd161.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          #FMO.index <- which(store.allFMO[,3] == "CD5")
          FMO.cd161 <- read.FCS(filename = paste0(store.allFMO[FMO.cd161.index,"Path"], "/", store.allFMO[FMO.cd161.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd161.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd161.index <- FMO.index[FMO.cd161.index]
          }
          FMO.cd161 <- read.FCS(filename = paste0(store.allFMO[FMO.cd161.index,"Path"], "/", store.allFMO[FMO.cd161.index,"Assay Date"], "/", store.allFMO[FMO.cd161.index,"FMO"]))
          
        }
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd161 <- removeMargins(FMO.cd161, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd161 <- removeMargins(FMO.cd161, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd161 <- compensateIMPC(FMO.cd161, basename(store.allFMO[FMO.cd161.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd161.fT <- transform(FMO.cd161, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd161@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd161.index) <- 0
        }else{
          FMO.cd161.Clean <- flowCut(FMO.cd161.fT, Segment = floor(nrow(FMO.cd161)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          
          FMO.cd161 <- transform(FMO.cd161, lgl)
          
          if(length(FMO.cd161.Clean$ind) > 0){
            FMO.cd161@exprs <- FMO.cd161@exprs[-FMO.cd161.Clean$ind, ]
          }
          
          remove(FMO.cd161.Clean)
          gc()
        }
        
      }
      
    }
    
    
    #############################################################################################
    ## Pre-processing the CD44 FMO
    
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd44.index <- grep("CD44", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd44.index <- grep("CD44", store.allFMO[FMO.index,3], value = FALSE)
        
      }
      
      if(length(FMO.cd44.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd44 <- read.FCS(filename = paste0(store.allFMO[FMO.cd44.index,"Path"], "/", store.allFMO[FMO.cd44.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd44.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd44.index <- FMO.index[FMO.cd44.index]
          }
          FMO.cd44 <- read.FCS(filename = paste0(store.allFMO[FMO.cd44.index,"Path"], "/", store.allFMO[FMO.cd44.index,"Assay Date"], "/", store.allFMO[FMO.cd44.index,"FMO"]))
          
        }
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd44 <- removeMargins(FMO.cd44, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd44 <- removeMargins(FMO.cd44, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd44 <- compensateIMPC(FMO.cd44, basename(store.allFMO[FMO.cd44.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd44.fT <- transform(FMO.cd44, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd44@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd44.index) <- 0
        }else{
          FMO.cd44.Clean <- flowCut(FMO.cd44.fT, Segment = floor(nrow(FMO.cd44)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          
          FMO.cd44 <- transform(FMO.cd44, lgl)
          
          if(length(FMO.cd44.Clean$ind) > 0){
            FMO.cd44@exprs <- FMO.cd44@exprs[-FMO.cd44.Clean$ind, ]
          }
          remove(FMO.cd44.Clean)
          gc()
        }
      }
      
      
    }
    
    #############################################################################################
    ## Pre-processing the CD25 FMO
    
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd25.index <- grep("CD25", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd25.index <- grep("CD25", store.allFMO[FMO.index,3], value = FALSE)
        
      }
      
      if(length(FMO.cd25.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd25 <- read.FCS(filename = paste0(store.allFMO[FMO.cd25.index,"Path"], "/", store.allFMO[FMO.cd25.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd25.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd25.index <- FMO.index[FMO.cd25.index]
          }
          FMO.cd25 <- read.FCS(filename = paste0(store.allFMO[FMO.cd25.index,"Path"], "/", store.allFMO[FMO.cd25.index,"Assay Date"], "/", store.allFMO[FMO.cd25.index,"FMO"]))
          
        }
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd25 <- removeMargins(FMO.cd25, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd25 <- removeMargins(FMO.cd25, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd25 <- compensateIMPC(FMO.cd25, basename(store.allFMO[FMO.cd25.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd25.fT <- transform(FMO.cd25, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd25@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd25.index) <- 0
        }else{
          FMO.cd25.Clean <- flowCut(FMO.cd25.fT, Segment = floor(nrow(FMO.cd25)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          
          FMO.cd25 <- transform(FMO.cd25, lgl)
          
          if(length(FMO.cd25.Clean$ind) > 0){
            FMO.cd25@exprs <- FMO.cd25@exprs[-FMO.cd25.Clean$ind, ]
          }
          remove(FMO.cd25.Clean)
          gc()
        }
      }
      
    }
    
    #############################################################################################
    ## Pre-processing the CD8A FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd8a.index <- grep("CD8A", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd8a.index <- grep("CD8A", store.allFMO[FMO.index,3], value = FALSE)
        
      }
      
      if(length(FMO.cd8a.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd8a <- read.FCS(filename = paste0(store.allFMO[FMO.cd8a.index,"Path"], "/", store.allFMO[FMO.cd8a.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd8a.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd8a.index <- FMO.index[FMO.cd8a.index]
          }
          FMO.cd8a <- read.FCS(filename = paste0(store.allFMO[FMO.cd8a.index,"Path"], "/", store.allFMO[FMO.cd8a.index,"Assay Date"], "/", store.allFMO[FMO.cd8a.index,"FMO"]))
          
        }
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd8a <- removeMargins(FMO.cd8a, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd8a <- removeMargins(FMO.cd8a, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd8a <- compensateIMPC(FMO.cd8a, basename(store.allFMO[FMO.cd8a.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd8a.fT <- transform(FMO.cd8a, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd8a@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd8a.index) <- 0
        }else{
          FMO.cd8a.Clean <- flowCut(FMO.cd8a.fT, Segment = floor(nrow(FMO.cd8a)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          
          FMO.cd8a <- transform(FMO.cd8a, lgl)
          
          if(length(FMO.cd8a.Clean$ind) > 0){
            FMO.cd8a@exprs <- FMO.cd8a@exprs[-FMO.cd8a.Clean$ind, ]
          }
          remove(FMO.cd8a.Clean)
          gc()
        }
        
      }
      
    }
    
    #############################################################################################
    ## Pre-processing the CD62L FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd62l.index <- grep("CD62L", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd62l.index <- grep("CD62L", store.allFMO[FMO.index,3], value = FALSE)
        
      }
      
      if(length(FMO.cd62l.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd62l <- read.FCS(filename = paste0(store.allFMO[FMO.cd62l.index,"Path"], "/", store.allFMO[FMO.cd62l.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd62l.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd62l.index <- FMO.index[FMO.cd62l.index]
          }
          FMO.cd62l <- read.FCS(filename = paste0(store.allFMO[FMO.cd62l.index,"Path"], "/", store.allFMO[FMO.cd62l.index,"Assay Date"], "/", store.allFMO[FMO.cd62l.index,"FMO"]))
        } 
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd62l <- removeMargins(FMO.cd62l, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd62l <- removeMargins(FMO.cd62l, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd62l <- compensateIMPC(FMO.cd62l, basename(store.allFMO[FMO.cd62l.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd62l.fT <- transform(FMO.cd62l, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd62l@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd62l.index) <- 0
        }else{
          FMO.cd62l.Clean <- flowCut(FMO.cd62l.fT, Segment = floor(nrow(FMO.cd62l)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          
          FMO.cd62l <- transform(FMO.cd62l, lgl)
          
          if(length(FMO.cd62l.Clean$ind) > 0){
            FMO.cd62l@exprs <- FMO.cd62l@exprs[-FMO.cd62l.Clean$ind, ]
          }
          remove(FMO.cd62l.Clean)
          gc()
        }
      }
      
    }  
    
    ###################################################################
    
    ## Pre-processing the TCRD FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        #FMO.index <- which(store.allFMO[,3] == "CD5")
        FMO.tcrd.index <- grep("TCRD", store.allFMO[,4], value = FALSE)
      }else{
        FMO.tcrd.index <- grep("TCRD", store.allFMO[FMO.index,3], value = FALSE)
      }
      
      if(length(FMO.tcrd.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          #FMO.index <- which(store.allFMO[,3] == "CD5")
          FMO.tcrd <- read.FCS(filename = paste0(store.allFMO[FMO.tcrd.index,"Path"], "/", store.allFMO[FMO.tcrd.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.tcrd.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.tcrd.index <- FMO.index[FMO.tcrd.index]
          }
          FMO.tcrd <- read.FCS(filename = paste0(store.allFMO[FMO.tcrd.index,"Path"], "/", store.allFMO[FMO.tcrd.index,"Assay Date"], "/", store.allFMO[FMO.tcrd.index,"FMO"]))
          
        }
        
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.tcrd <- removeMargins(FMO.tcrd, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.tcrd <- removeMargins(FMO.tcrd, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.tcrd <- compensateIMPC(FMO.tcrd, basename(store.allFMO[FMO.tcrd.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.tcrd.fT <- transform(FMO.tcrd, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd5@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.tcrd.index) <- 0
        }else{
          FMO.tcrd.Clean <- flowCut(FMO.tcrd.fT, Segment = floor(nrow(FMO.tcrd)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          FMO.tcrd <- transform(FMO.tcrd, lgl)
          
          if(length(FMO.tcrd.Clean$ind) > 0){
            FMO.tcrd@exprs <- FMO.tcrd@exprs[-FMO.tcrd.Clean$ind, ]
          }
          remove(FMO.tcrd.Clean)
          gc()
        }
      }
      
      
    }  
    
    ###################################################################
    
    ## Pre-processing the CD45 FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd45.index <- grep("CD45", store.allFMO[,4], value = FALSE)
      }else{
        FMO.cd45.index <- grep("CD45", store.allFMO[FMO.index,3], value = FALSE)
      }
      
      if(length(FMO.cd45.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd45 <- read.FCS(filename = paste0(store.allFMO[FMO.cd45.index,"Path"], "/", store.allFMO[FMO.cd45.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd45.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd45.index <- FMO.index[FMO.cd45.index]
          }
          FMO.cd45 <- read.FCS(filename = paste0(store.allFMO[FMO.cd45.index,"Path"], "/", store.allFMO[FMO.cd45.index,"Assay Date"], "/", store.allFMO[FMO.cd45.index,"FMO"]))
          
        }
        
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd45 <- removeMargins(FMO.cd45, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd45 <- removeMargins(FMO.cd45, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd45 <- compensateIMPC(FMO.cd45, basename(store.allFMO[FMO.cd45.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd45.fT <- transform(FMO.cd45, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd45@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd45.index) <- 0
        }else{
          FMO.cd45.Clean <- flowCut(FMO.cd45.fT, Segment = floor(nrow(FMO.cd45)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          FMO.cd45 <- transform(FMO.cd45, lgl)
          
          if(length(FMO.cd45.Clean$ind) > 0){
            FMO.cd45@exprs <- FMO.cd45@exprs[-FMO.cd45.Clean$ind, ]
          }
          remove(FMO.cd45.Clean)
          gc()
        }
      }
      
      
      
    }
    ###################################################################
    
    ## Pre-processing the CD4 FMO
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
        FMO.cd4.index <- grep("^CD4$", store.allFMO[,4], value = FALSE)  ## Using ^ to assert that it is the start of the pattern and $ to assert that it is the end.
      }else{
        FMO.cd4.index <- grep("^CD4$", store.allFMO[FMO.index,3], value = FALSE)
      }
      
      if(length(FMO.cd4.index) == 1){
        if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
          FMO.cd4 <- read.FCS(filename = paste0(store.allFMO[FMO.cd4.index,"Path"], "/", store.allFMO[FMO.cd4.index,"Panel/Organ/Folder"], "/", store.allFMO[FMO.cd4.index,"FMO"]))
          
        }else{
          if(centre == "tcp"){
            FMO.cd4.index <- FMO.index[FMO.cd4.index]
          }
          FMO.cd4 <- read.FCS(filename = paste0(store.allFMO[FMO.cd4.index,"Path"], "/", store.allFMO[FMO.cd4.index,"Assay Date"], "/", store.allFMO[FMO.cd4.index,"FMO"]))
          
        }
        
        
        # Remove scatter margins and compensate of the FMO---------------------------------------------
        # Removing margin events in Scatter channels
        FMO.cd4 <- removeMargins(FMO.cd4, chans = scat.chans, verbose = F)
        #Removing negative values in scatter channels
        FMO.cd4 <- removeMargins(FMO.cd4, chans = scat.chans, debris = T, neg = T, verbose = F)
        
        FMO.cd4 <- compensateIMPC(FMO.cd4, basename(store.allFMO[FMO.cd4.index,"FMO"]), fpath, centre = centre, panel.no = panel, assayDate = store.allFMO[FMO.cd5.index,"Assay Date"])
        
        FMO.cd4.fT <- transform(FMO.cd4, lgl)
        
        channels.to.clean <- which(complete.cases(FMO.cd45@parameters@data$desc) == TRUE)
        
        if(length(channels.to.clean) == 0){
          length(FMO.cd4.index) <- 0
        }else{
          FMO.cd4.Clean <- flowCut(FMO.cd4.fT, Segment = floor(nrow(FMO.cd4)*5/3000), Channels = channels.to.clean, Directory = paste0(results.dir,"/Figures/flowCut-FMO/"), Plot = 'Flagged Only')
          FMO.cd4 <- transform(FMO.cd4, lgl)
          
          if(length(FMO.cd4.Clean$ind) > 0){
            FMO.cd4@exprs <- FMO.cd4@exprs[-FMO.cd4.Clean$ind, ]
          }
          remove(FMO.cd4.Clean)
          gc()
        }
      }
      
      
    }
    #################################################################################################################
    #################################################################################################################
    
    
    ## FMO Controls: Gating All Events
    ## CD5 FMO
    
    if(centre != "ccp"){ ## Centre CCP has no FMOs
      if(length(FMO.cd5.index) == 1){
        # Quality Gate of CD5 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd5 <- qualityGate(FMO.cd5, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd5.flag <- 1
        if(results.FMO.cd5$live@proportion < 80){
          FMO.cd5.flag <- 0
        }
        else{
          scat.chans.FMO.cd5 <- results.FMO.cd5$scat.chans
          
          live.flowD.FMO.cd5 <- results.FMO.cd5$live
          
          FSCsinglets.flowD.FMO.cd5 <- results.FMO.cd5$FSCsinglets
          
          singlets.flowD.FMO.cd5 <- results.FMO.cd5$singlets
        }
        
      }
      
      ## CD161 FMO
      if(length(FMO.cd161.index) == 1){
        
        # Quality Gate of CD161 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd161 <- qualityGate(FMO.cd161, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd161.flag <- 1
        if(results.FMO.cd161$live@proportion < 80){
          FMO.cd161.flag <- 0
        }else{
          scat.chans.FMO.cd161 <- results.FMO.cd161$scat.chans
          
          live.flowD.FMO.cd161 <- results.FMO.cd161$live
          
          FSCsinglets.flowD.FMO.cd161 <- results.FMO.cd161$FSCsinglets
          
          singlets.flowD.FMO.cd161 <- results.FMO.cd161$singlets
        }
        
        
      }
      
      ## CD44 FMO
      if(length(FMO.cd44.index) == 1){
        # Quality Gate of CD44 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd44 <- qualityGate(FMO.cd44, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd44.flag <- 1
        if(results.FMO.cd44$live@proportion < 80){
          FMO.cd44.flag <- 0
        }else{
          scat.chans.FMO.cd44 <- results.FMO.cd44$scat.chans
          
          live.flowD.FMO.cd44 <- results.FMO.cd44$live
          
          FSCsinglets.flowD.FMO.cd44 <- results.FMO.cd44$FSCsinglets
          
          singlets.flowD.FMO.cd44 <- results.FMO.cd44$singlets
        }
        
      }
      
      ## CD25 FMO
      if(length(FMO.cd25.index) == 1){
        # Quality Gate of CD25 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd25 <- qualityGate(FMO.cd25, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd25.flag <- 1
        if(results.FMO.cd25$live@proportion < 80){
          FMO.cd25.flag <- 0
        }else{
          scat.chans.FMO.cd25 <- results.FMO.cd25$scat.chans
          
          live.flowD.FMO.cd25 <- results.FMO.cd25$live
          
          FSCsinglets.flowD.FMO.cd25 <- results.FMO.cd25$FSCsinglets
          
          singlets.flowD.FMO.cd25 <- results.FMO.cd25$singlets
        }
      }
      
      ## CD8A FMO
      if(length(FMO.cd8a.index) == 1){
        # Quality Gate of CD25 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd8a <- qualityGate(FMO.cd8a, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd8a.flag <- 1
        if(results.FMO.cd8a$live@proportion < 80){
          FMO.cd8a.flag <- 0
        }else{
          scat.chans.FMO.cd8a <- results.FMO.cd8a$scat.chans
          
          live.flowD.FMO.cd8a <- results.FMO.cd8a$live
          
          FSCsinglets.flowD.FMO.cd8a <- results.FMO.cd8a$FSCsinglets
          
          singlets.flowD.FMO.cd8a <- results.FMO.cd8a$singlets
        }
        
      }
      
      ## CD62L FMO
      if(length(FMO.cd62l.index) == 1){
        # Quality Gate of CD62l FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd62l <- qualityGate(FMO.cd62l, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd62l.flag <- 1
        if(results.FMO.cd62l$live@proportion < 80){
          FMO.cd62l.flag <- 0
        }else{
          scat.chans.FMO.cd62l <- results.FMO.cd62l$scat.chans
          
          live.flowD.FMO.cd62l <- results.FMO.cd62l$live
          
          FSCsinglets.flowD.FMO.cd62l <- results.FMO.cd62l$FSCsinglets
          
          singlets.flowD.FMO.cd62l <- results.FMO.cd62l$singlets
        }
      }
      
      ## TCRD FMO
      if(length(FMO.tcrd.index) == 1){
        # Quality Gate of TCRD FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.tcrd <- qualityGate(FMO.tcrd, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.tcrd.flag <- 1
        if(results.FMO.tcrd$live@proportion < 80){
          FMO.tcrd.flag <- 0
        }else{
          scat.chans.FMO.tcrd <- results.FMO.tcrd$scat.chans
          
          live.flowD.FMO.tcrd <- results.FMO.tcrd$live
          
          FSCsinglets.flowD.FMO.tcrd <- results.FMO.tcrd$FSCsinglets
          
          singlets.flowD.FMO.tcrd <- results.FMO.tcrd$singlets
        }
      }
      
      ## CD45 FMO
      if(length(FMO.cd45.index) == 1){
        # Quality Gate of CD45 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd45 <- qualityGate(FMO.cd45, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd45.flag <- 1
        if(results.FMO.cd45$live@proportion < 80){
          FMO.cd45.flag <- 0
        }else{
          scat.chans.FMO.cd45 <- results.FMO.cd45$scat.chans
          
          live.flowD.FMO.cd45 <- results.FMO.cd45$live
          
          FSCsinglets.flowD.FMO.cd45 <- results.FMO.cd45$FSCsinglets
          
          singlets.flowD.FMO.cd45 <- results.FMO.cd45$singlets
        }
      }
      
      
      ## CD4 FMO
      if(length(FMO.cd4.index) == 1){
        # Quality Gate of CD4 FMO (live/dead, singlets gating) ---------------------------------------------------
        results.FMO.cd4 <- qualityGate(FMO.cd4, scat.chans, channels.ind, centre, panel.no = panel)
        FMO.cd4.flag <- 1
        if(results.FMO.cd4$live@proportion < 80){
          FMO.cd4.flag <- 0
        }else{
          scat.chans.FMO.cd4 <- results.FMO.cd4$scat.chans
          
          live.flowD.FMO.cd4 <- results.FMO.cd4$live
          
          FSCsinglets.flowD.FMO.cd4 <- results.FMO.cd4$FSCsinglets
        }
        singlets.flowD.FMO.cd4 <- results.FMO.cd4$singlets
        
      }
    }
    ## COMMENT FROM HERE
    ###########################################################################################################
    ###########################################################################################################
    #
    # ## Gating SSC Singlets to obtain Lymphocytes.
    # ## This step is carried out only when CD45 marker is there.
    # ## Plotting CD45_CD161 to obtain Lymphocytes
    #
    # cd45marker.index <- which(names(channels.ind) == "CD45")
    # if(length(cd45marker.index) == 1){
    #   cd45.gate <-  deGate(singlets.flowD, channel = channels.ind["CD45"])
    #   lymphocytes.flowD <- flowDensity(singlets.flowD, channels = c(channels.ind["CD45"], channels.ind["CD161"]), position = c(T,NA), gates = c(cd45.gate, NA))
    #
    # }
    #
    # singlets <-rotate.data(singlets.flowD@flow.frame,c(channels.ind['CD45'], channels.ind['CD161']),theta = pi/4)$data
    # ## Plotting CD161_CD5 to obtain NK-cells
    # plotDens(singlets.flowD, channels = c(channels.ind['CD161'], channels.ind['CD5']), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    # lines(lymphocytes.flowD@filter, lwd=2)
    #
    # plotDens(singlets.flowD, channels = c(channels.ind['TCRD'], channels.ind['CD4']), main = "TCRD+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    #
    
    #########################################################################################
    # ## If TCRD marker is present, then gate SSC Singlets to obtain TCRD+ and TCRD-. Plotting TCRD_CD4 to obtain TCRD+ and
    # 
    # if(length(which(names(channels.ind)=="TCRD")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
    #     
    # }
    # plotDens(singlets.flowD, channels = c(channels.ind['TCRD'], channels.ind['CD4']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    
    #########################################################################################
    
    ## Gating SSC Singlets to obtain NK and NOT NK-cells. Plotting CD161_CD5 to obtain NOT NK-cells and NK-cells
    
    # NK-cells-----
    
    if(centre == "ucd"){
      
      theta0 <- 0.60
      rot <- rotate.data(getflowFrame(singlets.flowD), c(channels.ind["CD161"], channels.ind["CD5"]), theta = -theta0)$data
      
      cd5.gate.temp <- deGate(rot, channel = channels.ind["CD5"])
      rot.flowD.temp <- flowDensity(rot, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(NA, F), gates = c(NA, cd5.gate.temp))
      
      cd5.gate <- deGate(rot.flowD.temp, channel = channels.ind["CD5"])
      if(cd5.gate > 2){
        cd5.gate <- 2
      }
      cd161.gate <- deGate(rot.flowD.temp, channel = channels.ind["CD161"])
      if(cd161.gate < 1){
        cd161.gate <- deGate(rot.flowD.temp, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.97)
      }
      
      cd161.gate.NK <- cd161.gate
      
      
      # plotDens(rot.flowD.temp, c(channels.ind['CD161'], channels.ind['CD5']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
      # abline(h=cd5.gate); abline(v=cd161.gate)
      
      
      
      NKcells.flowD.temp <- flowDensity(rot.flowD.temp, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,T), gates = c(cd161.gate, cd5.gate), ellip.gate = T)
      NKcells.flowD.temp@filter <- rotate.data(NKcells.flowD.temp@filter, c(channels.ind["CD161"], channels.ind["CD5"]), theta = theta0)$data
      NKcells.flowD.temp@flow.frame <- rotate.data(NKcells.flowD.temp@flow.frame, c(channels.ind["CD161"], channels.ind["CD5"]),theta = theta0)$data
      
      NKcells.flowD <- NKcells.flowD.temp
      
      plotDens(singlets, c(channels.ind['CD161'], channels.ind['CD5']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(NKcells.flowD@filter)
      # abline(h=cd5.gate)
      # abline(v=cd161.gate)
      NKcells <- getflowFrame(NKcells.flowD)
      
    }else{
      
      cd5.gate <- deGate(singlets.flowD, channel = channels.ind["CD5"])
      if(cd5.gate >= 3){
        cd5.gate <- deGate(singlets.flowD, channel = channels.ind["CD5"], tinypeak.removal = 0.1)
      }
      
      cd161.gate.NK <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.999)
      if(cd161.gate.NK > 3.25){
        cd161.gate.NK <- max(c(deGate(singlets.flowD, channel = channels.ind["CD161"], use.upper = T, upper = T), deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.99)))+0.1
        
      }
      
      numPeaks.cd161.gate <- density(singlets@exprs[,c(channels.ind["CD161"])])
      
      cd161.gate.peak.lcn <- findpeaks(numPeaks.cd161.gate$y)
      cd161.gate.peak.lcn.index <- order(cd161.gate.peak.lcn[,1], decreasing = TRUE)
      cd161.gate.secondPeak.lcn <- numPeaks.cd161.gate$x[cd161.gate.peak.lcn[cd161.gate.peak.lcn.index[2],3]] ## Using findpeaks() to find the location where the second max peak starts
      
      if(cd161.gate.secondPeak.lcn < 1.2){
        ## I am changing this line of deGate. I need to check if it works for TCP
        ##cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"])
        cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.upper = T, upper = T)
        if(cd161.gate > cd161.gate.NK | abs(cd161.gate-cd161.gate.NK) < 0.5){
          if( abs(cd161.gate-cd161.gate.NK) < 0.5){
            cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile =0.75)
          }else{
            cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"])
          }
          
          
        }
        if(round(cd161.gate,1) < 1.4){ # round(cd161.gate,1) <= 1.5
          cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.96)
          if(cd161.gate < 1.6){
            cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.975)
          }
        }
      }else{
        cd161.gate <- mean(c(cd161.gate.secondPeak.lcn, deGate(singlets.flowD, channel = channels.ind["CD161"], use.upper = T, upper = T, tinypeak.removal = 0.1)))
        if(cd161.gate > 2.5){
          cd161.gate <- cd161.gate.secondPeak.lcn
        }else if(round(cd161.gate) < 2){
          cd161.gate <- cd161.gate.secondPeak.lcn
        }
      }
      
      
      # ## CD5 FMO & CD161 FMO
      # if(length(FMO.cd5.index) == 1 & length(FMO.cd161.index) == 1){
      #
      #
      #   NKcells.flowD.temp <- flowDensity(singlets, channels = c('APC-A', 'BV421-A'), position = c(T,NA), use.control = c(T,F), control = c(singlets.flowD.FMO.cd161, NA),  gates = c(cd161.gate, NA))
      #   NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c('APC-A', 'BV421-A'), position = c(NA,F), use.control = c(F,T), control = c(NA, singlets.flowD.FMO.cd5), use.upper = c(NA,T), upper = c(NA,T),  gates = c(NA, cd5.gate), ellip.gate = T)
      #   NKcells.flowD@proportion <- (NKcells.flowD@cell.count/nrow(singlets))*100
      #
      #
      # }else if(length(FMO.cd5.index) == 1 & length(FMO.cd161.index) == 0){
      #
      #   NKcells.flowD.temp <- flowDensity(singlets, channels = c('APC-A', 'BV421-A'), position = c(NA,F), use.control = c(F,T), control = c(NA, singlets.flowD.FMO.cd5), use.upper = c(NA,T), upper = c(NA,T), gates = c(NA, cd5.gate-0.2))
      #   NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c('APC-A', 'BV421-A'), position = c(T,NA),  gates = c(cd161.gate+0.2, NA), ellip.gate = T)
      #   NKcells.flowD@proportion <- (NKcells.flowD@cell.count/nrow(singlets))*100
      #
      #
      #
      # }else if(length(FMO.cd5.index) == 0 & length(FMO.cd161.index) == 1){
      #
      #   NKcells.flowD.temp <- flowDensity(singlets, channels = c('APC-A', 'BV421-A'), position = c(T,NA), use.control = c(T,F), control = c(singlets.flowD.FMO.cd161, NA),  gates = c(cd161.gate, NA),use.percentile = c(T,F), percentile = 0.9998)
      #   NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c('APC-A', 'BV421-A'), position = c(NA,F), gates = c(NA, cd5.gate-0.2), ellip.gate = T)
      #   NKcells.flowD@proportion <- (NKcells.flowD@cell.count/nrow(singlets))*100
      #
      # }else{
      #   NKcells.flowD <- flowDensity(singlets.flowD, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,F), gates = c(cd161.gate+0.2, cd5.gate-0.2), ellip.gate = T)
      #
      # }
      if(cd5.gate > 1.8){
        NKcells.flowD.temp <- flowDensity(singlets.flowD, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,F), gates = c(cd161.gate, cd5.gate-0.7), ellip.gate = F)
        NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(F,NA), gates = c(cd161.gate.NK, NA), ellip.gate = T)
      }else{
        NKcells.flowD.temp <- flowDensity(singlets.flowD, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,F), gates = c(cd161.gate, cd5.gate-0.2), ellip.gate = F)
        NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(F,NA), gates = c(cd161.gate.NK, NA), ellip.gate = T)
        
      }
      
      
      NKcells <- getflowFrame(NKcells.flowD)
      if((nrow(NKcells)/nrow(singlets))*100 < 1){
        cd161.gate.temp <- 0
        cd161.gate.temp <- cd161.gate
        cd161.gate <- deGate(singlets.flowD, channel = channels.ind["CD161"], use.percentile = T, percentile = 0.98)
        if(cd161.gate < 1.25 & cd161.gate.temp > 1.5 | round(cd5.gate-0.2,2) < 1.1){
          cd161.gate <- cd161.gate.temp
        }
        
        if( round(cd5.gate-0.2,2) < 1.1){
          NKcells.flowD.temp <- flowDensity(singlets.flowD, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,F), gates = c(cd161.gate, cd5.gate-0.2), ellip.gate = F)
        }else{
          NKcells.flowD.temp <- flowDensity(singlets.flowD, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(T,F), gates = c(cd161.gate, cd5.gate-0.5), ellip.gate = F)
          
        }
        NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c(channels.ind["CD161"], channels.ind["CD5"]), position = c(F,NA), gates = c(cd161.gate.NK, NA), ellip.gate = T)
        NKcells <- getflowFrame(NKcells.flowD)
      }
      
      ## Plotting CD161_CD5 to obtain NK-cells
      
      # #plotDens(singlets.flowD, channels = c(channels.ind['CD161'], channels.ind['CD5']), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      min.x <-  min(exprs(singlets)[, c(channels.ind["CD161"])])-1
      max.x <- max(exprs(singlets)[, c(channels.ind["CD161"])])+2
      plotDens(singlets, channels = c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x, max.x))
      lines(NKcells.flowD@filter, lwd=2); abline(v=cd161.gate)
      
    }
    
    
    
    
    all.gthres[4] <- cd5.gate
    all.gthres[5] <- cd161.gate
    all.gthres[6] <- cd161.gate.NK
    
    all.events[5] <- nrow(NKcells)
    all.props[5] <- (nrow(NKcells)/nrow(singlets))*100
    Filters.list$NKcells.filter <- NKcells.flowD@filter
    
    
    
    ################################
    # NOT NK-cells-----
    not.NKcells <- singlets
    not.NKcells@exprs <- not.NKcells@exprs[-NKcells.flowD@index,]
    
    all.events[6] <- nrow(not.NKcells)
    all.props[6] <- (nrow(not.NKcells)/nrow(singlets))*100
    
    if(centre != "ccp"){ ## CCP Centre has no FMO
      ## CD5 FMO
      if(length(FMO.cd5.index) == 1 & FMO.cd5.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd5, channel = grep('CD5', FMO.cd5@parameters@data$desc), tinypeak.removal = 0.1)
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd5, channel = grep('CD161', FMO.cd5@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.5
        
        NKcells.flowD.FMO.cd5 <- flowDensity(singlets.flowD.FMO.cd5, channels = c(grep('CD161', FMO.cd5@parameters@data$desc),grep('CD5', FMO.cd5@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        # NOT NK-cells-----
        not.NKcells.FMO.cd5 <- getflowFrame(singlets.flowD.FMO.cd5)
        not.NKcells.FMO.cd5@exprs <- not.NKcells.FMO.cd5@exprs[-NKcells.flowD.FMO.cd5@index,]
        
        #cd5.gate.FMO <- deGate(not.NKcells.FMO.cd5, channel = c(grep('CD5', FMO.cd5@parameters@data$desc)))
        
        # cd5.flowD.FMO.cd5 <- flowDensity(not.NKcells.FMO.cd5, channels = c(grep('CD25', FMO.cd5@parameters@data$desc),grep('CD5', FMO.cd5@parameters@data$desc)), position = c(NA, T), gates = c(NA, cd5.gate.FMO))
        
        # # plot(not.NKcells.FMO.cd5, cd5.flowD.FMO.cd5)
        # min.x <-  min(exprs(singlets.flowD.FMO.cd5@flow.frame)[, c(channels.ind["CD161"])])-1
        # max.x <- max(exprs(singlets.flowD.FMO.cd5@flow.frame)[, c(channels.ind["CD161"])])+2
        # plotDens(singlets.flowD.FMO.cd5, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NKcells.flowD.FMO.cd5@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd5)
        gc()
      }
      
      
      ## CD161 FMO
      if(length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd161, channel = grep('CD5', FMO.cd161@parameters@data$desc), tinypeak.removal = 0.1)
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd161, channel = grep('CD161', FMO.cd161@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.5
        
        NKcells.flowD.FMO.cd161 <- flowDensity(singlets.flowD.FMO.cd161, channels = c(grep('CD161', FMO.cd161@parameters@data$desc),grep('CD5', FMO.cd161@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        # NOT NK-cells-----
        not.NKcells.FMO.cd161 <- getflowFrame(singlets.flowD.FMO.cd161)
        not.NKcells.FMO.cd161@exprs <- not.NKcells.FMO.cd161@exprs[-NKcells.flowD.FMO.cd161@index,]
        
        if(NKcells.flowD.FMO.cd161@proportion < 0.05){
          length(FMO.cd161.index) <- 0
        }
        # cd5.flowD.FMO.cd161 <- flowDensity(not.NKcells.FMO.cd161, channels = c(grep('CD25', FMO.cd161@parameters@data$desc),grep('CD5', FMO.cd161@parameters@data$desc)), position = c(NA, T), gates = c(NA, cd5.gate.FMO))
        
        # min.x <-  min(exprs(singlets.flowD.FMO.cd161@flow.frame)[, c(channels.ind["CD161"])])-1
        # max.x <- max(exprs(singlets.flowD.FMO.cd161@flow.frame)[, c(channels.ind["CD161"])])+2
        # plotDens(singlets.flowD.FMO.cd161, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NKcells.flowD.FMO.cd161@filter, lwd=2)
        remove(singlets.flowD.FMO.cd161)
        gc()
      }
      
      
      ## CD44 FMO
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd44, channel = grep('CD5', FMO.cd44@parameters@data$desc), tinypeak.removal = 0.1)-0.3
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd44, channel = grep('CD161', FMO.cd44@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.5
        
        
        if(cd161.gate.FMO > 2.5){
          cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd44, channel = grep('CD161', FMO.cd44@parameters@data$desc), tinypeak.removal = 0.1)
          
        }
        NKcells.flowD.FMO.cd44 <- flowDensity(singlets.flowD.FMO.cd44, channels = c(grep('CD161', FMO.cd44@parameters@data$desc),grep('CD5', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        
        if(NKcells.flowD.FMO.cd44@proportion > 2.25 & cd161.gate.FMO < 1){ ## Changing the cd161.gate.FMO < 2 to cd161.gate.FMO < 1
          cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd44, channel = grep('CD161', FMO.cd44@parameters@data$desc), use.percentile = T, percentile = 0.98)
          
          NKcells.flowD.FMO.cd44 <- flowDensity(singlets.flowD.FMO.cd44, channels = c(grep('CD161', FMO.cd44@parameters@data$desc),grep('CD5', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
          
        }
        # NOT NK-cells-----
        not.NKcells.FMO.cd44 <- getflowFrame(singlets.flowD.FMO.cd44)
        not.NKcells.FMO.cd44@exprs <- not.NKcells.FMO.cd44@exprs[-NKcells.flowD.FMO.cd44@index,]
        
        plotDens(singlets.flowD.FMO.cd44, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(NKcells.flowD.FMO.cd44@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd44)
        gc()
      }
      
      ## CD25 FMO
      if(length(FMO.cd25.index) == 1 & FMO.cd25.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd25, channel = grep('CD5', FMO.cd25@parameters@data$desc), tinypeak.removal = 0.1)-0.3
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd25, channel = grep('CD161', FMO.cd25@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.5
        
        NKcells.flowD.FMO.cd25 <- flowDensity(singlets.flowD.FMO.cd25, channels = c(grep('CD161', FMO.cd25@parameters@data$desc),grep('CD5', FMO.cd25@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        if(NKcells.flowD.FMO.cd25@proportion < 2){
          cd161.gate.FMO <- 2
          NKcells.flowD.FMO.cd25 <- flowDensity(singlets.flowD.FMO.cd25, channels = c(grep('CD161', FMO.cd25@parameters@data$desc),grep('CD5', FMO.cd25@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
          
        }
        # NOT NK-cells-----
        not.NKcells.FMO.cd25 <- getflowFrame(singlets.flowD.FMO.cd25)
        not.NKcells.FMO.cd25@exprs <- not.NKcells.FMO.cd25@exprs[-NKcells.flowD.FMO.cd25@index,]
        
        plotDens(singlets.flowD.FMO.cd25, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(NKcells.flowD.FMO.cd25@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd25)
        gc()
      }
      
      
      ## CD8A FMO
      if(length(FMO.cd8a.index) == 1 & FMO.cd8a.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd8a, channel = grep('CD5', FMO.cd8a@parameters@data$desc), tinypeak.removal = 0.1)
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd8a, channel = grep('CD161', FMO.cd8a@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)
        
        NKcells.flowD.FMO.cd8a <- flowDensity(singlets.flowD.FMO.cd8a, channels = c(grep('CD161', FMO.cd8a@parameters@data$desc),grep('CD5', FMO.cd8a@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        if(NKcells.flowD.FMO.cd8a@proportion > 3){
          cd161.gate.FMO <- 2
          
          NKcells.flowD.FMO.cd8a <- flowDensity(singlets.flowD.FMO.cd8a, channels = c(grep('CD161', FMO.cd8a@parameters@data$desc),grep('CD5', FMO.cd8a@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
          
        }
        # NOT NK-cells-----
        not.NKcells.FMO.cd8a <- getflowFrame(singlets.flowD.FMO.cd8a)
        not.NKcells.FMO.cd8a@exprs <- not.NKcells.FMO.cd8a@exprs[-NKcells.flowD.FMO.cd8a@index,]
        
        # plotDens(singlets.flowD.FMO.cd8a, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NKcells.flowD.FMO.cd8a@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd8a)
        gc()
      }
      
      ## CD62L FMO
      if(length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd62l, channel = grep('CD5', FMO.cd62l@parameters@data$desc), tinypeak.removal = 0.1)-0.3
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd62l, channel = grep('CD161', FMO.cd62l@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.5
        
        if(cd161.gate.FMO > 2.5){
          cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd62l, channel = grep('CD161', FMO.cd62l@parameters@data$desc), tinypeak.removal = 0.1)
          
        }
        
        NKcells.flowD.FMO.cd62l <- flowDensity(singlets.flowD.FMO.cd62l, channels = c(grep('CD161', FMO.cd62l@parameters@data$desc), grep('CD5', FMO.cd62l@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        
        # NOT NK-cells-----
        not.NKcells.FMO.cd62l <- getflowFrame(singlets.flowD.FMO.cd62l)
        not.NKcells.FMO.cd62l@exprs <- not.NKcells.FMO.cd62l@exprs[-NKcells.flowD.FMO.cd62l@index,]
        
        plotDens(singlets.flowD.FMO.cd62l, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(NKcells.flowD.FMO.cd62l@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd62l)
        gc()
        
      }
      
      
      ## CD4 FMO
      if(length(FMO.cd4.index) == 1 & FMO.cd4.flag !=0){
        cd5.gate.FMO <- deGate(singlets.flowD.FMO.cd4, channel = grep('CD5', FMO.cd4@parameters@data$desc), tinypeak.removal = 0.1)-0.3
        cd161.gate.FMO <- deGate(singlets.flowD.FMO.cd4, channel = grep('CD161', FMO.cd4@parameters@data$desc), use.upper = T, upper = T, tinypeak.removal = 0.1)+0.2
        
        NKcells.flowD.FMO.cd4 <- flowDensity(singlets.flowD.FMO.cd4, channels = c(grep('CD161', FMO.cd4@parameters@data$desc), grep('CD5', FMO.cd4@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd5.gate.FMO), ellip.gate = T)
        
        
        # NOT NK-cells-----
        not.NKcells.FMO.cd4 <- getflowFrame(singlets.flowD.FMO.cd4)
        not.NKcells.FMO.cd4@exprs <- not.NKcells.FMO.cd4@exprs[-NKcells.flowD.FMO.cd4@index,]
        
        
        # plotDens(singlets.flowD.FMO.cd4, channels =  c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NKcells.flowD.FMO.cd4@filter, lwd=2)
        
        remove(singlets.flowD.FMO.cd4)
        gc()
        
      }
    }
    #########################################################################################
    ## Gating NK cells to obtain Effector NK cells and Resting NK cells.
    ## Plotting CD62L_CD44 to obtain Effector NK cells and Resting NK cells
    
    cd62.gate <- 0
    
    numPeaks.cd62.gate <- density(NKcells@exprs[,c(channels.ind["CD62L"])])
    
    cd62.gate.peak.lcn <- findpeaks(numPeaks.cd62.gate$y)
    if(nrow(cd62.gate.peak.lcn) > 1 & which.max(cd62.gate.peak.lcn[,1])==1){## If there are two major peaks and the max peak is the first peak, then we take the location where the first peak ends
      cd62.gate.maxPeak.lcn <- numPeaks.cd62.gate$x[cd62.gate.peak.lcn[which.max(cd62.gate.peak.lcn[,1]),4]] ## Using findpeaks() to find the location where the max peak ends
      if(cd62.gate.maxPeak.lcn > 3){
        cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.4)
        #cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.upper = T, upper = F)
        if(cd62.gate < 0){
          cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.15)
        }
      }
    }else{
      cd62.gate.maxPeak.lcn <- numPeaks.cd62.gate$x[cd62.gate.peak.lcn[which.max(cd62.gate.peak.lcn[,1]),3]] ## Using findpeaks() to find the location where the max peak starts
    }
    
    if(round(cd62.gate.maxPeak.lcn,2) < 1 & cd62.gate == 0){
      cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.1)
      if((1-cd62.gate) > 0.1){
        cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.upper = T, upper = F)
        if(cd62.gate < 1){
          cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.1)
          if(cd62.gate < 1){
            cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.25)
            if(round(cd62.gate,1) == 1.5 | round(cd62.gate,1) == 1.4){
              cd62.gate <- 1
            }
          }
        }
      }else if(cd62.gate > 1.5){
        cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.05)
      }
    }else{
      if(cd62.gate == 0){
        cd62.gate <- cd62.gate.maxPeak.lcn
        if(cd62.gate > 1.9 & cd62.gate < 2){
          cd62.gate <- round(cd62.gate.maxPeak.lcn)
        }
        if(cd62.gate > 3){
          cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.4)
        }
      }
      
    }
    
    
    cd44.gate.Peaks <- getPeaks(NKcells.flowD, channel = channels.ind["CD44"])
    
    if(length(cd44.gate.Peaks$Peaks) == 2){
      cd44.gate <- deGate(NKcells.flowD, channel = channels.ind["CD44"], use.upper = T, upper = F)
      
      if(round(cd44.gate,1) > 1.2){
        #if(cd44.gate > 2.4){
        cd44.gate <- deGate(NKcells.flowD, channel = channels.ind["CD44"], use.percentile = T, percentile = 0.001)
        
        NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
        NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate, cd44.gate))
        if(NK.Effector.flowD@proportion >= 45){
          cd44.gate <- deGate(NKcells.flowD, channel = channels.ind["CD44"])
          
          NK.temp <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(NA,T), gates = c(NA, cd44.gate))
          #plot(NKcells, NK.temp)
          cd62.gate <-  deGate(NK.temp, channel = channels.ind["CD62L"])
          
          NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
          
          NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate, cd44.gate))
        }
        
      }else{
        NK.temp <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(NA,T), gates = c(NA, cd44.gate))
        #plot(NKcells, NK.temp)
        cd62.gate <-  deGate(NK.temp, channel = channels.ind["CD62L"])
        
        NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
        
        NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate, cd44.gate))
        
        if(NK.Effector.flowD@proportion < 5 & NK.Resting.flowD@proportion < 30){
          cd44.gate <- deGate(NKcells.flowD, channel = channels.ind["CD44"], use.upper = T, upper = F)
          
          NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
          
          NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate, cd44.gate))
          
        }
      }
      
    }else{
      cd44.gate <- deGate(NKcells.flowD, channel = channels.ind["CD44"], use.percentile = T, percentile = 0.001)
      
      NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
      
      NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate, cd44.gate))
      
      
    }
    
    
    if(NK.Effector.flowD@proportion < 2){
      cd62.gate <- max(deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.1), deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.upper = T, upper = F))
      NK.Effector.flowD <- flowDensity(NKcells.flowD, channels =c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
      
      NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels =c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate, cd44.gate))
      
    }else if(NK.Effector.flowD@proportion > 20){
      cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.upper = T, upper = F, tinypeak.removal = 0.9)
      if(cd62.gate < 0.5){
        #cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.4)
        #cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"])
        cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], all.cuts=T)
        if(length(cd62.gate) == 2){
          cd62.gate <- max(deGate(NKcells.flowD, channel = channels.ind["CD62L"], all.cuts=T))
        }else{
          cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"])
          if(cd62.gate > 3){
            cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.5)
          }
        }
        
        if(cd62.gate < 0.5){
          cd62.gate <- deGate(NKcells.flowD, channel = channels.ind["CD62L"], use.percentile = T, percentile = 0.3)
        }
      }
      NK.Effector.flowD <- flowDensity(NKcells.flowD, channels =c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate))
      
      NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels =c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate, cd44.gate))
      
    }
    
    
    
    all.gthres[7] <- cd62.gate
    all.gthres[8] <- cd44.gate
    
    ## CD44 FMO
    # if(length(FMO.cd44.index) == 1 & length(FMO.cd62l.index) == 0){
    #
    #   NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(F,T), use.control = c(F,T), control = c(NA, NKcells.flowD.FMO.cd44), gates = c(cd62.gate, NA), use.percentile = c(F,T), percentile = c(NA,0.0001))
    #
    #   NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(T,NA), use.control = c(F,T), control = c(NA, NKcells.flowD.FMO.cd44), gates = c(cd62.gate, NA), use.percentile = c(F,T), percentile = c(NA,0.0001))
    #
    # }else if(length(FMO.cd44.index) == 0 & length(FMO.cd62l.index) == 1){ # FMO CD62L
    #
    #   NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(F,T), use.control = c(T,F), control = c(NKcells.flowD.FMO.cd62l, NA), gates = c(NA, cd44.gate))
    #
    #   NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(T,NA), use.control = c(T,F), control = c(NKcells.flowD.FMO.cd62l, NA))
    #
    #
    # }else if(length(FMO.cd44.index) == 1 & length(FMO.cd62l.index) == 1){ # FMO CD44 & FMO CD62L
    #
    #   NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(F,T), use.control = c(T,T), control = c(NKcells.flowD.FMO.cd62l, NKcells.flowD.FMO.cd44))
    #
    #   NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(T,NA), use.control = c(T,F), control = c(NKcells.flowD.FMO.cd62l, NA))
    #
    # }else{ ## Neither FMO CD44 or CD62L are present
    #
    #   NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(F,T), gates = c(cd62.gate, cd44.gate))
    #
    #   NK.Resting.flowD <-  flowDensity(NKcells.flowD, channels = c('APC-Cy7-A', 'PE-A'), position = c(T,NA), gates = c(cd62.gate, cd44.gate))
    #
    # }
    
    all.events[7] <- NK.Effector.flowD@cell.count
    all.props[7] <- NK.Effector.flowD@proportion
    all.events[8] <- NK.Resting.flowD@cell.count
    all.props[8] <- NK.Resting.flowD@proportion
    
    Filters.list$NK.Effector.filter <- NK.Effector.flowD@filter
    Filters.list$NK.Resting.filter <- NK.Resting.flowD@filter
    
    ## Plotting CD62L_CD44 to obtain Effector NK-cells & Resting NK-cells
    plotDens(NKcells.flowD, channels =  c(channels.ind['CD62L'], channels.ind['CD44']), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(v=cd62.gate, lwd= 2); #abline(h=NK.Effector.flowD@gates[2], lwd = 2)
    lines(NK.Effector.flowD@filter, lwd=2)
    lines(NK.Resting.flowD@filter, lwd=2)
    
    ## Removing the flowD objects of FMO NK cells and freeing up space
    remove(NKcells.flowD.FMO.cd5, NKcells.flowD.FMO.cd161, NKcells.flowD.FMO.cd44, NKcells.flowD.FMO.cd25,
           NKcells.flowD.FMO.cd8a, NKcells.flowD.FMO.cd62l, NKcells.flowD.FMO.cd4)
    gc()
    
    #########################################################################################
    ## Gating NK cells to obtain KLRG1+ NK cells, if the marker KLRG1 is present
    
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        klrg1.gate <- deGate(NKcells.flowD, channel = c(channels.ind['KLRG1']))
        
        klrg1.NK.flowD <- flowDensity(NKcells.flowD, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate, NA))
        klrg1.NK <- getflowFrame(klrg1.NK.flowD)
        
        all.gthres[21] <- klrg1.gate
        all.events[30] <- klrg1.NK.flowD@cell.count
        all.props[30] <- klrg1.NK.flowD@proportion
        
        
        plotDens(NKcells.flowD, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "NK cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=klrg1.gate, lwd=2)
        lines(klrg1.NK.flowD@filter, lwd=2)
        
      }else{
        all.gthres[21] <- NA
        all.events[30] <- NA
        all.props[30] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        klrg1.gate <- deGate(NKcells.flowD, channel = c(channels.ind['KLRG1']), use.upper = T, upper = 0.9)
        
        klrg1.NK.flowD <- flowDensity(NKcells.flowD, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate, NA))
        klrg1.NK <- getflowFrame(klrg1.NK.flowD)
        
        all.gthres[21] <- klrg1.gate
        all.events[30] <- klrg1.NK.flowD@cell.count
        all.props[30] <- klrg1.NK.flowD@proportion
        
        
        plotDens(NKcells.flowD, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "NK cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=klrg1.gate, lwd=2)
        lines(klrg1.NK.flowD@filter, lwd=2)
        
      }else{
        all.gthres[21] <- NA
        all.events[30] <- NA
        all.props[30] <- NA
      }
      
    }
    
    
    ##########################################################################################
    ## Gating NOT NK-cells to obtain CD5+ cells
    
    cd5.gate.cd5pos <- deGate(not.NKcells, channel = c(channels.ind['CD5']), tinypeak.removal = 0.9)
    if(cd5.gate.cd5pos < 1.97 & cd5.gate.cd5pos > 1.75){
      cd5.gate.cd5pos <- cd5.gate.cd5pos + 0.25
    }else if(cd5.gate.cd5pos < 1.75){ ## Changing this from < 1.5 to < 1.75
      cd5.gate.cd5pos <- deGate(not.NKcells, channel = c(channels.ind['CD5']))
    }
    all.gthres[9] <- cd5.gate.cd5pos
    
    ## CD5 FMO
    if(centre != "ccp"){
      if(length(FMO.cd5.index) == 1 & FMO.cd5.flag !=0){
        if(nrow(not.NKcells.FMO.cd5) < 20){
          length(FMO.cd5.index) <- 0
        }else{
          cd5.gate.FMO <- deGate(not.NKcells.FMO.cd5, channel = c(grep('CD5', FMO.cd5@parameters@data$desc)))
          if(cd5.gate.FMO < 1.5){
            cd5.gate.FMO <- 2
          }
          cd5.flowD.FMO.cd5 <- flowDensity(not.NKcells.FMO.cd5, channels = c(grep('CD25', FMO.cd5@parameters@data$desc),grep('CD5', FMO.cd5@parameters@data$desc)), position = c(NA, T), gates = c(NA, cd5.gate.FMO))
          if(cd5.flowD.FMO.cd5@proportion < 0.5){
            length(FMO.cd5.index) <- 0
          }
        }
        
        #plotDens(not.NKcells.FMO.cd5, channels =c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        #lines(cd5.flowD.FMO.cd5@filter, lwd =2)
      }
    }
    
    # NOT NK-cells-----
    if(centre == "ccp"){
      cd5.flowD <- flowDensity(not.NKcells, channels =c(channels.ind["CD25"], channels.ind["CD5"]), position = c(NA,T), gates = c(NA, cd5.gate.cd5pos))
      cd5 <- getflowFrame(cd5.flowD)
      
    }else{ ##CCP doesn't have any FMO 
      ## FMO CD5
      if(length(FMO.cd5.index) == 1 & FMO.cd5.flag !=0){
        
        #cd5.flowD <- flowDensity(not.NKcells, channels = c('PE-Cy7-A', 'BV421-A'), position = c(NA,T), use.control = c(F,T), control = c(NA, cd5.flowD.FMO.cd5), use.upper = c(F,T), upper = c(NA,T), gates = c(NA, cd5.gate))
        cd5.flowD <- flowDensity(not.NKcells, channels = c(channels.ind["CD25"], channels.ind["CD5"]), position = c(NA,T), use.control = c(F,T), control = c(NA, cd5.flowD.FMO.cd5), use.percentile = c(NA, T), percentile = c(NA, 0.4))
        if(cd5.flowD@gates[2] > 2.0 | cd5.flowD@gates[2] < 1){
          cd5.flowD <- flowDensity(not.NKcells, channels = c(channels.ind["CD25"], channels.ind["CD5"]), position = c(NA,T), gates = c(NA, cd5.gate.cd5pos))
          
        }
        cd5 <- getflowFrame(cd5.flowD)
        
      }else{
        
        cd5.flowD <- flowDensity(not.NKcells, channels =c(channels.ind["CD25"], channels.ind["CD5"]), position = c(NA,T), gates = c(NA, cd5.gate.cd5pos))
        cd5 <- getflowFrame(cd5.flowD)
        
      }
    }
    
    all.events[9] <- nrow(cd5)
    all.props[9] <- (nrow(cd5)/nrow(not.NKcells))*100
    
    Filters.list$cd5.filter <- cd5.flowD@filter
    
    ## Plotting CD25_CD5 to obtain CD5+ cells
    plotDens(not.NKcells, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
    lines(cd5.flowD@filter, lwd =2)
    
    
    if(centre != "ccp"){
      ## CD161 FMO
      if(length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0){
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd161, channel = c(grep('CD5', FMO.cd161@parameters@data$desc)))
        
        cd5.flowD.FMO.cd161 <- flowDensity(not.NKcells.FMO.cd161, channels = c(grep('CD25', FMO.cd161@parameters@data$desc),grep('CD5', FMO.cd161@parameters@data$desc)), position = c(NA, T), gates = c(NA, cd5.gate.FMO))
        # plotDens(not.NKcells.FMO.cd161, channels =c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        # lines(cd5.flowD.FMO.cd161@filter, lwd =2)
      }
      
      ## CD44 FMO
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd44, channel = grep('CD5', FMO.cd44@parameters@data$desc), tinypeak.removal = 0.1)-0.1
        
        cd5.flowD.FMO.cd44 <- flowDensity(not.NKcells.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc),grep('CD5', FMO.cd44@parameters@data$desc)), position = c(NA,T), gates = c(NA, cd5.gate.FMO))
        plotDens(not.NKcells.FMO.cd44, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        lines(cd5.flowD.FMO.cd44@filter, lwd =2)
      }
      
      ## CD25 FMO
      if(length(FMO.cd25.index) == 1 & FMO.cd25.flag !=0){
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd25, channel = grep('CD5', FMO.cd25@parameters@data$desc))
        
        cd5.flowD.FMO.cd25 <- flowDensity(not.NKcells.FMO.cd25, channels = c(grep('CD25', FMO.cd25@parameters@data$desc),grep('CD5', FMO.cd25@parameters@data$desc)), position = c(NA,T), gates = c(NA, cd5.gate.FMO))
        
        plotDens(not.NKcells.FMO.cd25, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        lines(cd5.flowD.FMO.cd25@filter, lwd =2)
      }
      
      ## CD8A FMO
      if(length(FMO.cd8a.index) == 1 & FMO.cd8a.flag !=0){
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd8a, channel = grep('CD5', FMO.cd8a@parameters@data$desc), tinypeak.removal = 0.1)-0.1
        
        cd5.flowD.FMO.cd8a <- flowDensity(not.NKcells.FMO.cd8a, channels = c(grep('CD25', FMO.cd8a@parameters@data$desc),grep('CD5', FMO.cd8a@parameters@data$desc)), position = c(NA,T), gates = c(NA, cd5.gate.FMO))
        # plotDens(not.NKcells.FMO.cd8a, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        # lines(cd5.flowD.FMO.cd8a@filter, lwd =2)
      }
      
      
      ## CD62L FMO
      if(length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0){
        
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd62l, channel = grep('CD5', FMO.cd62l@parameters@data$desc))
        
        cd5.flowD.FMO.cd62l <- flowDensity(not.NKcells.FMO.cd62l, channels = c(grep('CD25', FMO.cd62l@parameters@data$desc),grep('CD5', FMO.cd62l@parameters@data$desc)), position = c(NA,T), gates = c(NA, cd5.gate.FMO))
        # plotDens(not.NKcells.FMO.cd62l, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        # lines(cd5.flowD.FMO.cd62l@filter, lwd =2)
        
      }
      
      ## CD4 FMO
      if(length(FMO.cd4.index) == 1 & FMO.cd4.flag !=0){
        
        cd5.gate.FMO <- deGate(not.NKcells.FMO.cd4, channel = grep('CD5', FMO.cd4@parameters@data$desc))
        
        cd5.flowD.FMO.cd4 <- flowDensity(not.NKcells.FMO.cd4, channels = c(grep('CD25', FMO.cd4@parameters@data$desc),grep('CD5', FMO.cd4@parameters@data$desc)), position = c(NA,T), gates = c(NA, cd5.gate.FMO))
        # plotDens(not.NKcells.FMO.cd4, channels = c(channels.ind["CD25"], channels.ind["CD5"]), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
        # lines(cd5.flowD.FMO.cd4@filter, lwd =2)
        
      }
      
    }
    #########################################################################################
    ## Gating CD5+ cells to obtain CD161+ CD8- cells and NOT CD161+ CD8- cells
    
    # CD161+ CD8- cells-----
    
    # #numPeaks.cd5pos.cd8 <- getPeaks(cd5, channel = channels.ind['CD8'], tinypeak.removal = 0.9)
    # numPeaks.cd5pos.cd8 <- getPeaks(cd5, channel = channels.ind['CD8'])
    #
    # if(length(numPeaks.cd5pos.cd8$Peaks) != 2){
    #   flaggedFile[1] <- x$FCS.files
    #   #next
    # }
    numPeaks.cd5cells <- density(cd5@exprs[,c(channels.ind["CD8"])])
    
    cd5cells.peak.lcn <- findpeaks(numPeaks.cd5cells$y)
    cd5cells.firstPeak.lcn <- numPeaks.cd5cells$x[cd5cells.peak.lcn[1,4]] ## Using findpeaks() to find the location where the first peak ends
    
    if(cd5cells.firstPeak.lcn < 1.5){ ## Changing this from < 0.5 to < 1.5
      cd8a.gate <- max(c(cd5cells.firstPeak.lcn, deGate(cd5, channel = c(channels.ind["CD8"]))))
    }else{
      cd8a.gate <- cd5cells.firstPeak.lcn
    }
    
    
    
    all.gthres[10] <- cd8a.gate
    
    
    # ## I had to add this extra step because for some files in GMC there is a staining issue in CD8 marker which is part of NOT CD161+CD8+ population.
    # cd5pos.cd8.flowD.temp <- flowDensity(cd5, channels = c(channels.ind['CD161'], channels.ind['CD8']), position = c(NA, T))
    # if(cd5pos.cd8.flowD.temp@proportion < 10){
    #   flaggedFile[1] <- x$FCS.files
    #   #next()
    # }else{
    #   flaggedFile[1] <- NA
    # }
    
    cd161.gate.high <- deGate(cd5, channel = c(channels.ind["CD161"]), use.percentile = T, percentile = 0.9995)
    if(cd161.gate.high > 4){
      cd161.gate.high <- deGate(cd5, channel = c(channels.ind["CD161"]), use.upper = T, upper = T, tinypeak.removal = 0.01)
    }
    all.gthres[11] <- cd161.gate.high
    
    numPeaks.cd161.gate.low <- density(cd5@exprs[,c(channels.ind["CD161"])])
    
    cd161.gate.low.peak.lcn <- findpeaks(numPeaks.cd161.gate.low$y)
    cd161.gate.low <- numPeaks.cd161.gate.low$x[cd161.gate.low.peak.lcn[which.max(cd161.gate.low.peak.lcn[,1]),4]] ## Using findpeaks() to find the location where the max peak ends
    
    cd161.gate.low.temp <- deGate(cd5, channel = c(channels.ind["CD161"]), use.upper = T, upper=T)
    
    if(cd161.gate.low > 2.2 | (cd161.gate.low != cd161.gate.low.temp)){
      cd161.gate.low <- min(c(cd161.gate.low,cd161.gate.low.temp))
      if(cd161.gate.low < 1.25){
        cd161.gate.low <- deGate(cd5, channel = c(channels.ind["CD161"]))
        if(cd161.gate.low < 1.1){
          cd161.gate.low <- deGate(cd5, channel = c(channels.ind["CD161"]), use.upper = T, upper=T)
        }
      }
    }
    
    
    all.gthres[12] <- cd161.gate.low
    
    if(centre == "ccp"){
      cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
      if(cd161pos.cd8neg.flowD.temp@proportion < 0.6){ ## changing from <0.5 to 0.6
        cd161.gate.low <- deGate(cd5, channel = c(channels.ind["CD161"]), use.percentile = T, percentile = 0.975)
        cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
        
      }
      cd161pos.cd8neg.flowD <- flowDensity(cd161pos.cd8neg.flowD.temp, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(F,NA), gates = c(cd161.gate.high, NA))
      cd161pos.cd8neg.flowD@proportion <- (cd161pos.cd8neg.flowD@cell.count/nrow(cd5))*100
      
    }else{
      ## FMO CD161 & FMO CD8A
      if(length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0 | (length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0) & (length(FMO.cd8a.index) == 1 & FMO.cd8a.flag !=0)){
        
        if(cd5.flowD.FMO.cd161@gates[1] < 1){
          cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
          if(cd161pos.cd8neg.flowD.temp@proportion < 0.6){ ## changing from <0.5 to 0.6
            cd161.gate.low <- deGate(cd5, channel = c(channels.ind["CD161"]), use.percentile = T, percentile = 0.975)
            cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
            
          }
          cd161pos.cd8neg.flowD <- flowDensity(cd161pos.cd8neg.flowD.temp, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(F,NA), gates = c(cd161.gate.high, NA))
          cd161pos.cd8neg.flowD@proportion <- (cd161pos.cd8neg.flowD@cell.count/nrow(cd5))*100
          
        }else{
          cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels = c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), use.control = c(T,F), control = c(cd5.flowD.FMO.cd161, NA))
          
          cd161pos.cd8neg.flowD <- flowDensity(cd161pos.cd8neg.flowD.temp, channels = c(channels.ind["CD161"], channels.ind["CD8"]), position = c(F,NA), gates = c(cd161.gate.high, cd8a.gate))
          cd161pos.cd8neg.flowD@proportion <- (cd161pos.cd8neg.flowD@cell.count/nrow(cd5))*100
          
        }
        # cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels = c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), use.control = c(T,F), control = c(cd5.flowD.FMO.cd161, NA), gates = c(cd161.gate.low, NA))
        
        
      }else if(length(FMO.cd161.index) == 0 & length(FMO.cd8a.index) == 1){
        
        # cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), use.control = c(F,T), control = c(NA, cd5.flowD.FMO.cd8a), gates = c(cd161.gate.low, cd8a.gate))
        
        cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), use.control = c(F,T), control = c(NA, cd5.flowD.FMO.cd8a), gates = c(cd161.gate.low, NA))
        if(cd161pos.cd8neg.flowD.temp@gates[2] < 2){## Changing this to < 2 from <0
          cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
          
        }
        cd161pos.cd8neg.flowD <- flowDensity(cd161pos.cd8neg.flowD.temp, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(F,NA), gates = c(cd161.gate.high, NA))
        cd161pos.cd8neg.flowD@proportion <- (cd161pos.cd8neg.flowD@cell.count/nrow(cd5))*100
        
        
      }else{
        cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
        if(cd161pos.cd8neg.flowD.temp@proportion < 0.6){ ## changing from <0.5 to 0.6
          cd161.gate.low <- deGate(cd5, channel = c(channels.ind["CD161"]), use.percentile = T, percentile = 0.975)
          cd161pos.cd8neg.flowD.temp <- flowDensity(cd5, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(T,F), gates = c(cd161.gate.low, cd8a.gate))
          
        }
        cd161pos.cd8neg.flowD <- flowDensity(cd161pos.cd8neg.flowD.temp, channels =  c(channels.ind["CD161"], channels.ind["CD8"]), position = c(F,NA), gates = c(cd161.gate.high, NA))
        cd161pos.cd8neg.flowD@proportion <- (cd161pos.cd8neg.flowD@cell.count/nrow(cd5))*100
        
        
      }
      
    }
    
    
    cd161pos.cd8neg <- getflowFrame(cd161pos.cd8neg.flowD)
    
    
    all.events[10] <- nrow(cd161pos.cd8neg)
    all.props[10] <- (nrow(cd161pos.cd8neg)/nrow(cd5))*100
    Filters.list$cd161pos.cd8neg.filter <- cd161pos.cd8neg.flowD@filter
    
    # NOT CD161+ CD8- cells------
    not.cd161pos.cd8neg <- cd5
    not.cd161pos.cd8neg@exprs <- not.cd161pos.cd8neg@exprs[-cd161pos.cd8neg.flowD@index,]
    
    all.events[11] <- nrow(not.cd161pos.cd8neg)
    all.props[11] <- (nrow(not.cd161pos.cd8neg)/nrow(cd5))*100
    
    Filters.list$cd161pos.cd8neg.filter <- cd161pos.cd8neg.flowD@filter
    
    
    ## Plotting CD161_CD8a to obtain CD161+ CD8- cells and NOT CD161+ CD8- cells
    plotDens(cd5.flowD, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd161pos.cd8neg.flowD@filter, lwd=2); abline(v=cd161.gate.low)
    
    #################################
    
    ## CD5 FMO ####Temporary Comment: For GMC issue with FMO CD5. So not using any further. No change in the code for FMO CD5
    if(centre != "ccp"){
      if(length(FMO.cd5.index) == 1 & FMO.cd5.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd5, channel = c(grep('CD8a', FMO.cd5@parameters@data$desc)))
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd5, channel = c(grep('CD161', FMO.cd5@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd5, channel = c(grep('CD161', FMO.cd5@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd5, channel = c(grep('CD161', FMO.cd5@parameters@data$desc)), use.percentile = T, percentile = 0.9995)
        
        
        cd161pos.cd8neg.flowD.FMO.cd5.temp <- flowDensity(cd5.flowD.FMO.cd5, channels = c(grep('CD161', FMO.cd5@parameters@data$desc),grep('CD8a', FMO.cd5@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd5 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd5.temp, channels = c(grep('CD161', FMO.cd5@parameters@data$desc),grep('CD8a', FMO.cd5@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd5@proportion <- (cd161pos.cd8neg.flowD.FMO.cd5@cell.count/cd5.flowD.FMO.cd5@cell.count)*100
        
        # NOT CD161+ CD8- cells------
        not.cd161pos.cd8neg.FMO.cd5 <- getflowFrame(cd5.flowD.FMO.cd5)
        not.cd161pos.cd8neg.FMO.cd5@exprs <- not.cd161pos.cd8neg.FMO.cd5@exprs[-cd161pos.cd8neg.flowD.FMO.cd5@index,]
        
        # plotDens(cd5.flowD.FMO.cd5, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(cd161pos.cd8neg.flowD.FMO.cd5@filter, lwd=2)
        
      }
      
      
      
      ## CD161 FMO.  ## CD161 is not used in subsequent gating for TCP and GMC since there were issues for FMOs in both centres.
      ## This part of the code needs to be worked on when we find useful FMO CD161 for other centres.
      if(length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0){
        
        ## Noticed problems in CD8a marker in the FMO CD161 for some files, hence added tryCatch here to catch the error
        ## Example of files with problem: PANEL_A_ABDP_111_C3C03005.fcs
        cd8a.gate.FMO <- tryCatch({deGate(cd5.flowD.FMO.cd161, channel = c(grep('CD8a', FMO.cd161@parameters@data$desc)))},error = function(err) {
          
          return(0)
        }
        ) # end of tryCatch
        if(cd8a.gate.FMO == 0){
          FMO.cd161.index <- 0
        }else{
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd161, channel = c(grep('CD161', FMO.cd161@parameters@data$desc)))
          cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd161, channel = c(grep('CD161', FMO.cd161@parameters@data$desc)), use.percentile = T, percentile = 0.9995)
          
          
          cd161pos.cd8neg.flowD.FMO.cd161.temp <- flowDensity(cd5.flowD.FMO.cd161, channels = c(grep('CD161', FMO.cd161@parameters@data$desc),grep('CD8a', FMO.cd161@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
          cd161pos.cd8neg.flowD.FMO.cd161 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd161.temp, channels = c(grep('CD161', FMO.cd161@parameters@data$desc),grep('CD8a', FMO.cd161@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
          cd161pos.cd8neg.flowD.FMO.cd161@proportion <- (cd161pos.cd8neg.flowD.FMO.cd161@cell.count/cd5.flowD.FMO.cd161@cell.count)*100
          
          
          # NOT CD161+ CD8- cells------
          not.cd161pos.cd8neg.FMO.cd161 <- getflowFrame(cd5.flowD.FMO.cd161)
          not.cd161pos.cd8neg.FMO.cd161@exprs <- not.cd161pos.cd8neg.FMO.cd161@exprs[-cd161pos.cd8neg.flowD.FMO.cd161@index,]
          
          # plotDens(cd5.flowD.FMO.cd161, channels = c('APC-A', 'PE-CF594-A'), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(cd161pos.cd8neg.flowD.FMO.cd161@filter, lwd=2)
        }
        
      }
      
      
      ## CD44 FMO
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd44, channel = c(grep('CD8a', FMO.cd44@parameters@data$desc)))
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd44, channel = c(grep('CD161', FMO.cd44@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd44, channel = c(grep('CD161', FMO.cd44@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd44, channel = c(grep('CD161', FMO.cd44@parameters@data$desc)), use.percentile = T, percentile = 0.9995)
        
        
        cd161pos.cd8neg.flowD.FMO.cd44.temp <- flowDensity(cd5.flowD.FMO.cd44, channels = c(grep('CD161', FMO.cd44@parameters@data$desc),grep('CD8a', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd44 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd44.temp, channels = c(grep('CD161', FMO.cd44@parameters@data$desc),grep('CD8a', FMO.cd44@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd44@proportion <- (cd161pos.cd8neg.flowD.FMO.cd44@cell.count/cd5.flowD.FMO.cd44@cell.count)*100
        
        # NOT CD161+ CD8- cells------
        not.cd161pos.cd8neg.FMO.cd44 <- getflowFrame(cd5.flowD.FMO.cd44)
        not.cd161pos.cd8neg.FMO.cd44@exprs <- not.cd161pos.cd8neg.FMO.cd44@exprs[-cd161pos.cd8neg.flowD.FMO.cd44@index,]
        
        # plotDens(cd5.flowD.FMO.cd44, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(cd161pos.cd8neg.flowD.FMO.cd44@filter, lwd=2)
        
      }
      
      
      ## FMO CD25
      if(length(FMO.cd25.index) == 1 & FMO.cd25.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd25, channel = c(grep('CD8a', FMO.cd25@parameters@data$desc)))-0.1
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd25, channel = c(grep('CD161', FMO.cd25@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd25, channel = c(grep('CD161', FMO.cd25@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd25, channel = c(grep('CD161', FMO.cd25@parameters@data$desc)), use.percentile = T, percentile = 0.9999)
        
        
        cd161pos.cd8neg.flowD.FMO.cd25.temp <- flowDensity(cd5.flowD.FMO.cd25, channels = c(grep('CD161', FMO.cd25@parameters@data$desc),grep('CD8a', FMO.cd25@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd25 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd25.temp, channels = c(grep('CD161', FMO.cd25@parameters@data$desc),grep('CD8a', FMO.cd25@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd25@proportion <- (cd161pos.cd8neg.flowD.FMO.cd25@cell.count/cd5.flowD.FMO.cd25@cell.count)*100
        
        #FMO.cd25.flag <- 1
        if(cd161pos.cd8neg.flowD.FMO.cd25@proportion <= 0.5){
          #FMO.cd25.index <- 0
          FMO.cd25.flag <- 0
        }else{
          # NOT CD161+ CD8- cells------
          not.cd161pos.cd8neg.FMO.cd25 <- getflowFrame(cd5.flowD.FMO.cd25)
          not.cd161pos.cd8neg.FMO.cd25@exprs <- not.cd161pos.cd8neg.FMO.cd25@exprs[-cd161pos.cd8neg.flowD.FMO.cd25@index,]
          
          # plotDens(cd5.flowD.FMO.cd25, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(cd161pos.cd8neg.flowD.FMO.cd25@filter, lwd=2)
        }
        
        
      }
      
      
      ## CD8a FMO seemed to have the same staining problem that was seen in CD62l FMO for TCP
      ## CD62L FMO-----Review by Mehrnoush
      ## Mehrnoush mentioned that CD62L may not have been stained properly with CD8 marker; hence the problem with NOT CD161+Cd8- population
      ## So we will only use CD161+Cd8- population from  FMO CD8a. Also we will not be using not.cd161pos.cd8neg.FMO.cd8a as control for gating
      ## NOT CD161+CD8- population.
      
      ## CD8a FMO
      if(length(FMO.cd8a.index) == 1 & FMO.cd8a.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd8a, channel = c(grep('CD8a', FMO.cd8a@parameters@data$desc)))-0.1
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd8a, channel = c(grep('CD161', FMO.cd8a@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd8a, channel = c(grep('CD161', FMO.cd8a@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd8a, channel = c(grep('CD161', FMO.cd8a@parameters@data$desc)), use.percentile = T, percentile = 0.9999)
        
        
        cd161pos.cd8neg.flowD.FMO.cd8a.temp <- flowDensity(cd5.flowD.FMO.cd8a, channels = c(grep('CD161', FMO.cd8a@parameters@data$desc),grep('CD8a', FMO.cd8a@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd8a <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd8a.temp, channels = c(grep('CD161', FMO.cd8a@parameters@data$desc),grep('CD8a', FMO.cd8a@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd8a@proportion <- (cd161pos.cd8neg.flowD.FMO.cd8a@cell.count/cd5.flowD.FMO.cd8a@cell.count)*100
        
        cd5.FMO.cd8a <- getflowFrame(cd5.flowD.FMO.cd8a)
        
        ## Finding the number of peaks across the CD161 channel of the FMO CD8. Setting tinypeak.removal to 0.9 in order to remove tiny peaks.
        numPeaks.cd5.FMO.cd8a.cd161 <- getPeaks(cd5.FMO.cd8a, channel = c(channels.ind['CD161']), tinypeak.removal = 0.9)
        
        
        if(cd161pos.cd8neg.flowD.FMO.cd8a@proportion <= 0.5 | length(numPeaks.cd5.FMO.cd8a.cd161$Peaks) == 1){
          length(FMO.cd8a.index) <- 0
        }else{
          # NOT CD161+ CD8- cells------
          not.cd161pos.cd8neg.FMO.cd8a <- getflowFrame(cd5.flowD.FMO.cd8a)
          not.cd161pos.cd8neg.FMO.cd8a@exprs <- not.cd161pos.cd8neg.FMO.cd8a@exprs[-cd161pos.cd8neg.flowD.FMO.cd8a@index,]
          
          # plotDens(cd5.flowD.FMO.cd8a, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(cd161pos.cd8neg.flowD.FMO.cd8a@filter, lwd=2)
        }
      }
      
      
      ## CD62L FMO-----Review by Mehrnoush. The following comment is for TCP FMO CD62L
      ## Mehrnoush mentioned that CD62L may not have been stained properly with CD8 marker; hence the problem with NOT CD161+Cd8- population
      ## So we will only use CD161+Cd8- population from gating FMO CD62L
      if(length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd62l, channel = c(grep('CD8a', FMO.cd62l@parameters@data$desc)))-0.1
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd62l, channel = c(grep('CD161', FMO.cd62l@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd62l, channel = c(grep('CD161', FMO.cd62l@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd62l, channel = c(grep('CD161', FMO.cd62l@parameters@data$desc)), use.percentile = T, percentile = 0.9999)
        
        
        cd161pos.cd8neg.flowD.FMO.cd62l.temp <- flowDensity(cd5.flowD.FMO.cd62l, channels = c(grep('CD161', FMO.cd62l@parameters@data$desc),grep('CD8a', FMO.cd62l@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd62l <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd62l.temp, channels = c(grep('CD161', FMO.cd62l@parameters@data$desc),grep('CD8a', FMO.cd62l@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd62l@proportion <- (cd161pos.cd8neg.flowD.FMO.cd62l@cell.count/cd5.flowD.FMO.cd62l@cell.count)*100
        
        #FMO.cd62l.flag <- 1
        if(cd161pos.cd8neg.flowD.FMO.cd62l@proportion <= 0.5){
          #FMO.cd25.index <- 0
          FMO.cd62l.flag <- 0
        }else{
          # NOT CD161+ CD8- cells------
          not.cd161pos.cd8neg.FMO.cd62l <- getflowFrame(cd5.flowD.FMO.cd62l)
          not.cd161pos.cd8neg.FMO.cd62l@exprs <- not.cd161pos.cd8neg.FMO.cd62l@exprs[-cd161pos.cd8neg.flowD.FMO.cd62l@index,]
          
          # plotDens(cd5.flowD.FMO.cd62l, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(cd161pos.cd8neg.flowD.FMO.cd62l@filter, lwd=2)
        }
        
        
      }
      
      ## CD4 FMO
      if(length(FMO.cd4.index) == 1 & FMO.cd4.flag !=0){
        cd8a.gate.FMO <- deGate(cd5.flowD.FMO.cd4, channel = c(grep('CD8a', FMO.cd4@parameters@data$desc)))-0.1
        cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd4, channel = c(grep('CD161', FMO.cd4@parameters@data$desc)))
        
        if(cd161.gate.FMO < 1){
          cd161.gate.FMO <- deGate(cd5.flowD.FMO.cd4, channel = c(grep('CD161', FMO.cd4@parameters@data$desc)), use.upper = T, upper = T)
          
        }
        cd161.gate.high.FMO <- deGate(cd5.flowD.FMO.cd4, channel = c(grep('CD161', FMO.cd4@parameters@data$desc)), use.percentile = T, percentile = 0.9999)
        
        
        cd161pos.cd8neg.flowD.FMO.cd4.temp <- flowDensity(cd5.flowD.FMO.cd4, channels = c(grep('CD161', FMO.cd4@parameters@data$desc),grep('CD8a', FMO.cd4@parameters@data$desc)), position = c(T,F), gates = c(cd161.gate.FMO, cd8a.gate.FMO))
        cd161pos.cd8neg.flowD.FMO.cd4 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd4.temp, channels = c(grep('CD161', FMO.cd4@parameters@data$desc),grep('CD8a', FMO.cd4@parameters@data$desc)), position = c(F,NA), gates = c(cd161.gate.high.FMO, NA))
        cd161pos.cd8neg.flowD.FMO.cd4@proportion <- (cd161pos.cd8neg.flowD.FMO.cd4@cell.count/cd5.flowD.FMO.cd4@cell.count)*100
        
        FMO.cd4.flag <- 1
        if(cd161pos.cd8neg.flowD.FMO.cd4@proportion <= 0.5){
          #FMO.cd25.index <- 0
          FMO.cd4.flag <- 0
        }else{
          # NOT CD161+ CD8- cells------
          not.cd161pos.cd8neg.FMO.cd4 <- getflowFrame(cd5.flowD.FMO.cd4)
          not.cd161pos.cd8neg.FMO.cd4@exprs <- not.cd161pos.cd8neg.FMO.cd4@exprs[-cd161pos.cd8neg.flowD.FMO.cd4@index,]
          
          # plotDens(cd5.flowD.FMO.cd4, channels = c(channels.ind["CD161"], channels.ind["CD8"]), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(cd161pos.cd8neg.flowD.FMO.cd4@filter, lwd=2)
        }
        
        
      }
    }
    
    
    #########################################################################################
    ## Gating CD161+ CD8- cells to obtain CD4+ NKT cells and CD4- NKT cells
    ## Using CD25 marker instead of GITR for TCP and GMC (since they don't use GITR marker)
    ## Using GITR marker for UCD
    
    cd4.gate <- deGate(not.cd161pos.cd8neg, channel = c(channels.ind["CD4"]))
    all.gthres[13] <- cd4.gate
    
    if(length(which(names(channels.ind)=="GITR")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
      
      gitr.gate <- deGate(not.cd161pos.cd8neg, channel = c(channels.ind['GITR']), use.upper = T, upper=F)
      
      cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), position = c(T,T), gates = c(gitr.gate, cd4.gate), ellip.gate = T)
      
      cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), position = c(T,F), gates = c(gitr.gate, cd4.gate-0.55), ellip.gate = T)
      
      if(round(cd4neg.NKT.flowD@proportion) <= 25 | sum(cd4pos.NKT.flowD@proportion, cd4neg.NKT.flowD@proportion) > 98){
        cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), position = c(T,T), gates = c(gitr.gate, cd4.gate))
        
        cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), position = c(T,F), gates = c(gitr.gate, cd4.gate))
        
      }
      
      all.gthres[14] <- NA ## cd25
      all.gthres[22] <- gitr.gate
      
      
      ## Plotting GITR_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
      plotDens(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd4pos.NKT.flowD@filter, lwd = 2)
      lines(cd4neg.NKT.flowD@filter, lwd = 2)
      
      
    }else{
      
      cd25.gate.low <- deGate(cd161pos.cd8neg, channel = c(channels.ind["CD25"]), use.upper = T, upper = F)
      if(is.na(cd25.gate.low)){
        cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD161'], channels.ind['CD4']), position = c(T,T), gates = c(cd161.gate.low, cd4.gate), ellip.gate = T)
        
        cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD161'], channels.ind['CD4']), position = c(T,F), gates = c(cd161.gate.low, cd4.gate-0.55), ellip.gate = T)
        
        all.gthres[14] <- NA ## CD25 marker issue
        all.gthres[22] <- NA ## gitr
        
        
        ## Plotting CD25_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
        plotDens(cd161pos.cd8neg, channels = c(channels.ind['CD161'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd4pos.NKT.flowD@filter, lwd = 2)
        lines(cd4neg.NKT.flowD@filter, lwd = 2)
      }else{
        ## FMO 25
        if(length(FMO.cd25.index) == 1 &  FMO.cd25.flag!= 0){
          cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,T), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate), ellip.gate = T)
          
          cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,F), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate-0.52), ellip.gate = T)
          
          if(round(cd4neg.NKT.flowD@proportion) <= 25 | sum(cd4pos.NKT.flowD@proportion, cd4neg.NKT.flowD@proportion) > 95){
            cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,T), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate))
            
            cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,F), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate))
            
            if(cd4pos.NKT.flowD@proportion < 30){
              cd4.gate <- deGate(not.cd161pos.cd8neg, channel = c(channels.ind["CD4"]), use.percentile = T, percentile = 0.55)
              cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,T), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate))
              
              cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,F), use.control = c(T,F), control = c(cd161pos.cd8neg.flowD.FMO.cd25, NA), use.upper = c(T, F), upper = F, gates = c(NA, cd4.gate))
              
            }
          }
          
          
        }else{
          cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,T), gates = c(cd25.gate.low, cd4.gate), ellip.gate = T)
          
          cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,F), gates = c(cd25.gate.low, cd4.gate-0.55), ellip.gate = T)
          
          if(round(cd4neg.NKT.flowD@proportion) <= 25 | sum(cd4pos.NKT.flowD@proportion, cd4neg.NKT.flowD@proportion) > 98){
            cd4pos.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,T), gates = c(cd25.gate.low, cd4.gate))
            
            cd4neg.NKT.flowD <- flowDensity(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), position = c(T,F), gates = c(cd25.gate.low, cd4.gate))
            
          }
          
        }
        
        
        all.gthres[14] <- cd25.gate.low
        all.gthres[22] <- NA ## gitr
        
        
        ## Plotting CD25_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
        plotDens(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd4pos.NKT.flowD@filter, lwd = 2)
        lines(cd4neg.NKT.flowD@filter, lwd = 2)
        
        
        
      }
      
    }
    
    
    cd4pos.NKT <- getflowFrame(cd4pos.NKT.flowD)
    
    
    all.events[12] <- nrow(cd4pos.NKT)
    all.props[12] <- (nrow(cd4pos.NKT)/nrow(cd161pos.cd8neg))*100
    Filters.list$cd4pos.NKT.filter <- cd4pos.NKT.flowD@filter
    
    cd4neg.NKT <- getflowFrame(cd4neg.NKT.flowD)
    
    all.events[13] <- nrow(cd4neg.NKT)
    all.props[13] <- (nrow(cd4neg.NKT)/nrow(cd161pos.cd8neg))*100
    
    Filters.list$cd4neg.NKT.filter <- cd4neg.NKT.flowD@filter
    
    ## Flagging files
    numPeaks.cd4 <- getPeaks(cd161pos.cd8neg, channel = channels.ind['CD4'])
    if (length(numPeaks.cd4$Peaks) == 1 | cd4pos.NKT.flowD@proportion < 25){
      flaggedFile[1] <- x$FCS.files
      #next()
    }else{
      flaggedFile[1] <- NA
    }
    
    if(centre != "ccp"){
      ## CD161 FMO ####Temporary Comment: For GMC issue with FMO CD161. So not using any further
      if(length(FMO.cd161.index) == 1 & FMO.cd161.flag !=0){
        ## Add code when working with dataset from other centres (not GMC)
      }
      
      
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd4.gate.high.FMO <- deGate(cd161pos.cd8neg.flowD.FMO.cd44, channel = c(grep('CD4', FMO.cd44@parameters@data$desc)))
        cd4.gate.low.FMO <- cd4.gate.high.FMO - 0.25
        numPeaks.cd25.FMO <- getPeaks(cd161pos.cd8neg.flowD.FMO.cd44, channel = channels.ind["CD25"])
        if(length(numPeaks.cd25.FMO$Peaks) < 2){
          #FMO.cd44.flag <- 1
          FMO.cd44.flag <- 0
        }else{
          cd25.gate.FMO <- deGate(cd161pos.cd8neg.flowD.FMO.cd44, channel = c(grep('CD25', FMO.cd44@parameters@data$desc)), use.upper = T, upper = F)
          cd4pos.NKT.flowD.FMO.cd44 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,T), gates = c(cd25.gate.FMO, cd4.gate.high.FMO), ellip.gate = T)
          
          cd4neg.NKT.flowD.FMO.cd44 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd25.gate.FMO, cd4.gate.low.FMO), ellip.gate = T)
          
          if((cd4pos.NKT.flowD.FMO.cd44@proportion + cd4neg.NKT.flowD.FMO.cd44@proportion) > 100){
            cd4pos.NKT.flowD.FMO.cd44 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,T), gates = c(cd25.gate.FMO, cd4.gate.high.FMO), ellip.gate = F)
            
            cd4neg.NKT.flowD.FMO.cd44 <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd25.gate.FMO, cd4.gate.high.FMO), ellip.gate = F)
            
          }
          ## Plotting CD25_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
          plotDens(cd161pos.cd8neg.flowD.FMO.cd44, channels = c(channels.ind['CD25'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          lines(cd4pos.NKT.flowD.FMO.cd44@filter, lwd = 2)
          lines(cd4neg.NKT.flowD.FMO.cd44@filter, lwd = 2)
        }
        
      }
      
      
      
      # ## CD8A FMO ####Temporary Comment: For GMC issue with FMO CD8A. So not using any further
      # if(length(FMO.cd8a.index) == 1){
      #   ## Add code when working with dataset from other centres (not GMC)
      # }
      
      
      ## FMO CD62L
      if(length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0){
        cd4.gate.high.FMO <- deGate(cd161pos.cd8neg.flowD.FMO.cd62l, channel = c(grep('CD4', FMO.cd62l@parameters@data$desc)))
        cd4.gate.low.FMO <- cd4.gate.high.FMO - 0.5
        cd25.gate.FMO <- deGate(cd161pos.cd8neg.flowD.FMO.cd62l, channel = c(grep('CD25', FMO.cd62l@parameters@data$desc)), use.upper = T, upper = F)
        if(is.na(cd25.gate.FMO)){
          FMO.cd62l.flag <- 0
        }else{
          cd4pos.NKT.flowD.FMO.cd62l <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd62l, channels = c(grep('CD25', FMO.cd62l@parameters@data$desc), grep('CD4', FMO.cd62l@parameters@data$desc)), position = c(T,T), gates = c(cd25.gate.FMO, cd4.gate.high.FMO), ellip.gate = T)
          
          cd4neg.NKT.flowD.FMO.cd62l <- flowDensity(cd161pos.cd8neg.flowD.FMO.cd62l, channels = c(grep('CD25', FMO.cd62l@parameters@data$desc), grep('CD4', FMO.cd62l@parameters@data$desc)), position = c(T,F), gates = c(cd25.gate.FMO, cd4.gate.low.FMO), ellip.gate = T)
          
          ## Plotting CD25_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
          plotDens(cd161pos.cd8neg.flowD.FMO.cd62l, channels = c(channels.ind['CD25'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          lines(cd4pos.NKT.flowD.FMO.cd62l@filter, lwd = 2)
          lines(cd4neg.NKT.flowD.FMO.cd62l@filter, lwd = 2)
        }
        
      }
    }   
    
    #########################################################################################
    ## Gating CD4+ NKT cells to obtain Effector CD4+ NKT and Resting CD4+ NKT cells
    cd62.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind["CD62L"]), use.upper = T, upper = T)
    if(cd62.gate.NKT > 1.9){ ## Changing this from >2 to >1.9
      cd62.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind["CD62L"]))
      if(cd62.gate.NKT < 1 | cd62.gate.NKT > 2){
        cd62.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind["CD62L"]), tinypeak.removal = 0.9)
      }
    }
    
    if(centre == "ccp"){
      cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
      
      cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
      
      
      if(cd4pos.NKT.Resting.flowD@proportion > 15){
        cd62.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind["CD62L"]))
        
        cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
        
        cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
        
        # cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c('APC-Cy7-A', 'PE-A'), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
        #
        #
        # cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c('APC-Cy7-A', 'PE-A'), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
      }
    }else{
      if((length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0) & (length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0)){
        cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(T,T), control = c(cd4pos.NKT.flowD.FMO.cd62l, cd4pos.NKT.flowD.FMO.cd44))
        
        cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), use.control = c(T,NA), control = c(cd4pos.NKT.flowD.FMO.cd62l, NA))
        
        
        if(cd4pos.NKT.Resting.flowD@proportion > 15){
          cd62.gate.NKT <- deGate(cd4neg.NKT, channel = c(channels.ind["CD62L"]))
          
          cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
          
          if(cd4pos.NKT.Effector.flowD@proportion < 5 & cd4pos.NKT.Effector.flowD@gates[2] > 1){
            cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
            
          }
          cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
        }
        
      }else if(length(FMO.cd62l.index) == 1 & length(FMO.cd44.index) == 0){
        cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(T,NA), control = c(cd4pos.NKT.flowD.FMO.cd62l, NA), gates = c(NA, cd44.gate))
        
        cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), use.control = c(T,NA), control = c(cd4pos.NKT.flowD.FMO.cd62l, NA))
        
        
        if(cd4pos.NKT.Resting.flowD@proportion > 15){
          cd62.gate.NKT <- deGate(cd4neg.NKT, channel = c(channels.ind["CD62L"]))
          
          cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
          
          
          cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
        }
        
        
        
      }else if(length(FMO.cd62l.index) == 0 & length(FMO.cd44.index) == 1){
        cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
        
        cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
        
        if(cd4pos.NKT.Resting.flowD@proportion > 15){
          cd62.gate.NKT <- deGate(cd4neg.NKT, channel = c(channels.ind["CD62L"]))
          
          cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
          
          cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
        }
        
        
      }else{
        cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
        
        cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
        
        if(cd4pos.NKT.Resting.flowD@proportion > 15){
          cd62.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind["CD62L"]))
          
          cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
          
          cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
          
          # cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKT, channels = c('APC-Cy7-A', 'PE-A'), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4pos.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
          #
          #
          # cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKT, channels = c('APC-Cy7-A', 'PE-A'), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
        }
        
        
      }
      
      
    }
    
    all.gthres[15] <- cd62.gate.NKT
    
    cd4pos.NKT.Effector <- getflowFrame(cd4pos.NKT.Effector.flowD)
    
    
    all.events[14] <- nrow(cd4pos.NKT.Effector)
    all.props[14] <- (nrow(cd4pos.NKT.Effector)/nrow(cd4pos.NKT))*100
    Filters.list$cd4pos.NKT.Effector.filter <- cd4pos.NKT.Effector.flowD@filter
    
    cd4pos.NKT.Resting <- getflowFrame(cd4pos.NKT.Resting.flowD)
    
    
    all.events[15] <- nrow(cd4pos.NKT.Resting)
    all.props[15] <- (nrow(cd4pos.NKT.Resting)/nrow(cd4pos.NKT))*100
    Filters.list$cd4pos.NKT.Resting.filter <- cd4pos.NKT.Resting.flowD@filter
    
    ## Plotting CD62L_CD44 to obtain Effector CD4+ NKT cells and Resting CD4+ NKT cells
    plotDens(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4pos.NKT.Effector.flowD@filter, lwd = 2)
    lines(cd4pos.NKT.Resting.flowD@filter, lwd = 2)
    
    #########################################################################################
    ## Gating CD4+NKT cells to obtain KLRG1+ CD4+NKT cells, if the marker KLRG1 is present
    ##  For UCD I am using KLRG1 vs CD44
    
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        if(klrg1.gate < 2){
          klrg1.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind['KLRG1']))
          if(klrg1.gate.NKT > 3.5){
            klrg1.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind['KLRG1']), use.percentile = T, percentile = 0.95)
          }
          klrg1.cd4pos.NKT.flowD <- flowDensity(cd4pos.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate.NKT, NA))
          
          
        }else{
          klrg1.cd4pos.NKT.flowD <- flowDensity(cd4pos.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate, NA))
          
        }
        
        
        klrg1.cd4pos.NKT <- getflowFrame(klrg1.cd4pos.NKT.flowD)
        
        
        all.events[31] <- klrg1.cd4pos.NKT.flowD@cell.count
        all.props[31] <- klrg1.cd4pos.NKT.flowD@proportion
        
        
        plotDens(cd4pos.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=klrg1.gate)
        lines(klrg1.cd4pos.NKT.flowD@filter, lwd=2)
        
      }else{
        
        all.events[31] <- NA
        all.props[31] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        klrg1.cd4pos.NKT.flowD <- flowDensity(cd4pos.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate, NA))
        klrg1.cd4pos.NKT <- getflowFrame(klrg1.cd4pos.NKT.flowD)
        
        if(klrg1.cd4pos.NKT.flowD@proportion < 1){
          klrg1.gate.NKT <- deGate(cd4pos.NKT, channel = c(channels.ind['KLRG1']), use.percentile = T, percentile = 0.99)
          klrg1.cd4pos.NKT.flowD <- flowDensity(cd4pos.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate.NKT, NA))
          klrg1.cd4pos.NKT <- getflowFrame(klrg1.cd4pos.NKT.flowD)
          
        }
        all.events[31] <- klrg1.cd4pos.NKT.flowD@cell.count
        all.props[31] <- klrg1.cd4pos.NKT.flowD@proportion
        
        
        plotDens(cd4pos.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);#(abline(v=klrg1.gate))
        lines(klrg1.cd4pos.NKT.flowD@filter, lwd=2)
        
      }else{
        
        all.events[31] <- NA
        all.props[31] <- NA
      }
      
    }
    
    
    #########################################################################################
    ## Gating CD4- NKT cells to obtain Effector CD4- NKT and Resting CD4- NKT cells
    if(centre == "ccp"){
      cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
      
      cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
      
    }else{
      if((length(FMO.cd62l.index) == 1 & FMO.cd62l.flag !=0) & (length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0)){
        cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(T,T), control = c(cd4neg.NKT.flowD.FMO.cd62l, cd4neg.NKT.flowD.FMO.cd44))
        
        cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), use.control = c(T,NA), control = c(cd4neg.NKT.flowD.FMO.cd62l, NA))
        
        if(cd4neg.NKT.Resting.flowD@proportion > 55){
          cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
          
          cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
          
        }
        
      }else if(length(FMO.cd62l.index) == 1 & length(FMO.cd44.index) == 0){
        cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(T,NA), control = c(cd4neg.NKT.flowD.FMO.cd62l, NA), gates = c(NA, cd44.gate))
        
        cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), use.control = c(T,NA), control = c(cd4neg.NKT.flowD.FMO.cd62l, NA))
        
        if(cd4neg.NKT.Resting.flowD@proportion > 55){
          cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
          
          cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
          
          
        }
        
      }else if(length(FMO.cd62l.index) == 0 & length(FMO.cd44.index) == 1){
        cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), use.control = c(NA,T), control = c(NA, cd4neg.NKT.flowD.FMO.cd44), gates = c(cd62.gate.NKT, NA))
        
        cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
        
      }else{
        cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F, T), gates = c(cd62.gate.NKT, cd44.gate))
        
        cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T, NA), gates = c(cd62.gate.NKT, NA))
        
        
      }
      
      
    }
    
    cd4neg.NKT.Effector <- getflowFrame(cd4neg.NKT.Effector.flowD)
    
    
    all.events[16] <- nrow(cd4neg.NKT.Effector)
    all.props[16] <- (nrow(cd4neg.NKT.Effector)/nrow(cd4neg.NKT))*100
    Filters.list$cd4neg.NKT.Effector.filter <- cd4neg.NKT.Effector.flowD@filter
    
    cd4neg.NKT.Resting <- getflowFrame(cd4neg.NKT.Resting.flowD)
    
    all.events[17] <- nrow(cd4neg.NKT.Resting)
    all.props[17] <- (nrow(cd4neg.NKT.Resting)/nrow(cd4neg.NKT))*100
    Filters.list$cd4neg.NKT.Resting.filter <- cd4neg.NKT.Resting.flowD@filter
    
    ## Plotting CD62L_CD44 to obtain Effector CD4- NKT cells and Resting CD4- NKT cells
    plotDens(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD4- NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd4neg.NKT.Effector.flowD@filter, lwd = 2)
    lines(cd4neg.NKT.Resting.flowD@filter, lwd = 2)
    
    
    #########################################################################################
    ## Gating CD4-NKT cells to obtain KLRG1+ CD4-NKT cells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        if(klrg1.gate < 2){
          klrg1.gate.NKT.neg <- deGate(cd4neg.NKT, channel = c(channels.ind['KLRG1']))
          klrg1.cd4neg.NKT.flowD <- flowDensity(cd4neg.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate.NKT.neg, NA))
          
          
        }else{
          klrg1.cd4neg.NKT.flowD <- flowDensity(cd4neg.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate, NA))
          
        }
        
        
        klrg1.cd4neg.NKT <- getflowFrame(klrg1.cd4neg.NKT.flowD)
        
        
        all.events[32] <- klrg1.cd4neg.NKT.flowD@cell.count
        all.props[32] <- klrg1.cd4neg.NKT.flowD@proportion
        
        
        plotDens(cd4neg.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd4neg.NKT.flowD@filter, lwd=2)
        
      }else{
        
        all.events[32] <- NA
        all.props[32] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        #klrg1.gate <- deGate(cd4neg.NKT, channel = c(channels.ind['KLRG1']))
        
        klrg1.cd4neg.NKT.flowD <- flowDensity(cd4neg.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate, NA))
        klrg1.cd4neg.NKT <- getflowFrame(klrg1.cd4neg.NKT.flowD)
        
        if(klrg1.cd4neg.NKT.flowD@proportion < 1){
          klrg1.gate.NKT <- deGate(cd4neg.NKT, channel = c(channels.ind['KLRG1']), use.percentile = T, percentile = 0.99)
          klrg1.cd4neg.NKT.flowD <- flowDensity(cd4neg.NKT, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate.NKT, NA))
          klrg1.cd4neg.NKT <- getflowFrame(klrg1.cd4neg.NKT.flowD)
          
        }
        
        
        all.events[32] <- klrg1.cd4neg.NKT.flowD@cell.count
        all.props[32] <- klrg1.cd4neg.NKT.flowD@proportion
        
        
        plotDens(cd4neg.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4- NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(v=klrg1.gate)
        lines(klrg1.cd4neg.NKT.flowD@filter, lwd=2)
        
      }else{
        
        all.events[32] <- NA
        all.props[32] <- NA
      }
      
    }
    
    
    
    
    ########################################################################################
    ## Gating NOT CD161+ CD8- cells to obtain T-cells
    ## TCP & GMC: Because of the staining problem of FMO CD8a we will not be using FMOs for gating NOT CD161+ CD8- and its subsequent populations
    
    NOT.Tcells.flowD <- flowDensity(not.cd161pos.cd8neg, channels = c(channels.ind['CD8'], channels.ind['CD4']), position = c(F,F), gates = c(cd8a.gate, cd4.gate))
    
    Tcells <- not.cd161pos.cd8neg
    Tcells@exprs <- Tcells@exprs[-NOT.Tcells.flowD@index,]
    
    all.events[18] <- nrow(Tcells)
    all.props[18] <- (nrow(Tcells)/nrow(not.cd161pos.cd8neg))*100
    
    Filters.list$NOT.Tcells.filter <- NOT.Tcells.flowD@filter
    
    ## Plotting CD8_CD4 to obtain T cells
    plotDens(not.cd161pos.cd8neg, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(NOT.Tcells.flowD@filter, lwd = 2)
    
    if(centre != "ccp"){
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd8a.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd44, channel = c(grep('CD8', FMO.cd44@parameters@data$desc)))
        cd4.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd44, channel = c(grep('CD4', FMO.cd44@parameters@data$desc)))
        
        NOT.Tcells.flowD.FMO.cd44 <- flowDensity(not.cd161pos.cd8neg.FMO.cd44, channels = c(grep('CD8', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(F,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        
        Tcells.FMO.cd44 <- not.cd161pos.cd8neg.FMO.cd44
        Tcells.FMO.cd44@exprs <- Tcells.FMO.cd44@exprs[-NOT.Tcells.flowD.FMO.cd44@index,]
        
        # ## Plotting CD8_CD4 to obtain T cells
        # plotDens(not.cd161pos.cd8neg.FMO.cd44, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NOT.Tcells.flowD.FMO.cd44@filter, lwd = 2)
        
      }
      
      
      ## FMO CD25
      if(length(FMO.cd25.index) == 1 & FMO.cd25.flag != 0){
        cd8a.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd25, channel = c(grep('CD8', FMO.cd25@parameters@data$desc)))
        cd4.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd25, channel = c(grep('CD4', FMO.cd25@parameters@data$desc)))
        
        NOT.Tcells.flowD.FMO.cd25 <- flowDensity(not.cd161pos.cd8neg.FMO.cd25, channels = c(grep('CD8', FMO.cd25@parameters@data$desc),grep('CD4', FMO.cd25@parameters@data$desc)), position = c(F,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        
        Tcells.FMO.cd25 <- not.cd161pos.cd8neg.FMO.cd25
        Tcells.FMO.cd25@exprs <- Tcells.FMO.cd25@exprs[-NOT.Tcells.flowD.FMO.cd25@index,]
        
        
        # ## Plotting CD8_CD4 to obtain T cells
        # plotDens(not.cd161pos.cd8neg.FMO.cd25, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NOT.Tcells.flowD.FMO.cd25@filter, lwd = 2)
        
      }
      
      ## FMO CD62L
      if(length(FMO.cd62l.index) == 1 & FMO.cd62l.flag != 0){
        cd8a.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd62l, channel = c(grep('CD8', FMO.cd62l@parameters@data$desc)))
        cd4.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd62l, channel = c(grep('CD4', FMO.cd62l@parameters@data$desc)))
        
        NOT.Tcells.flowD.FMO.cd62l <- flowDensity(not.cd161pos.cd8neg.FMO.cd62l, channels = c(grep('CD8', FMO.cd62l@parameters@data$desc),grep('CD4', FMO.cd62l@parameters@data$desc)), position = c(F,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        
        Tcells.FMO.cd62l <- not.cd161pos.cd8neg.FMO.cd62l
        Tcells.FMO.cd62l@exprs <- Tcells.FMO.cd62l@exprs[-NOT.Tcells.flowD.FMO.cd62l@index,]
        
        
        # ## Plotting CD8_CD4 to obtain T cells
        # plotDens(not.cd161pos.cd8neg.FMO.cd62l, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NOT.Tcells.flowD.FMO.cd62l@filter, lwd = 2)
        
      }
      
      
      ## GMC: Staining issues in FMO CD4. So we will not be using FMO CD4 any further
      ## FMO CD4
      if(length(FMO.cd4.index) == 1 & FMO.cd4.flag != 0){
        cd8a.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd4, channel = c(grep('CD8', FMO.cd4@parameters@data$desc)))
        cd4.gate.FMO <- deGate(not.cd161pos.cd8neg.FMO.cd4, channel = c(grep('CD4', FMO.cd4@parameters@data$desc)))
        
        ## Finding the number of peaks across the CD4 channel of the FMO CD4.
        numPeaks.not.cd161pos.cd8neg.FMO.cd4.cd4 <- getPeaks(not.cd161pos.cd8neg.FMO.cd4, channel = c(channels.ind['CD4']))
        
        if(length(numPeaks.not.cd161pos.cd8neg.FMO.cd4.cd4$Peaks) == 1){
          length(FMO.cd4.index) <- 0
        }
        NOT.Tcells.flowD.FMO.cd4 <- flowDensity(not.cd161pos.cd8neg.FMO.cd4, channels = c(grep('CD8', FMO.cd4@parameters@data$desc),grep('CD4', FMO.cd4@parameters@data$desc)), position = c(F,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        
        Tcells.FMO.cd4 <- not.cd161pos.cd8neg.FMO.cd4
        Tcells.FMO.cd4@exprs <- Tcells.FMO.cd4@exprs[-NOT.Tcells.flowD.FMO.cd4@index,]
        
        
        # ## Plotting CD8_CD4 to obtain T cells
        # plotDens(not.cd161pos.cd8neg.FMO.cd4, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(NOT.Tcells.flowD.FMO.cd4@filter, lwd = 2)
        
      }
      
    }
    
    ###########################################################################################################
    
    ## Gating T-cells to obtain CD4 T-cells and CD8 T-cells
    ## For both TCP and GMC we are not using FMOs CD8 and CD4 for gating this population due to issues with the FMOs
    cd4.gate.high <- deGate(Tcells, channel = c(channels.ind["CD4"]), use.upper = T, upper = T, tinypeak.removal = 0.9)
    if(cd4.gate.high < cd4.gate){
      cd4.gate.high <- deGate(Tcells, channel = c(channels.ind["CD4"]), use.upper = T, upper = T)
      
    }
    all.gthres[16] <- cd4.gate.high
    
    cd4.Tcells.flowD.temp <- flowDensity(Tcells, channels = c(channels.ind['CD8'], channels.ind['CD4']), position = c(F,T), gates = c(cd8a.gate, cd4.gate-0.25))
    cd4.Tcells.flowD <- flowDensity(cd4.Tcells.flowD.temp, channels = c(channels.ind['CD8'], channels.ind['CD4']), position = c(NA,F), gates = c(NA, cd4.gate.high))
    
    cd8.Tcells.flowD <- flowDensity(Tcells, channels = c(channels.ind['CD8'], channels.ind['CD4']), position = c(T,F), gates = c(cd8a.gate, cd4.gate+0.5))
    
    
    cd4.Tcells <- getflowFrame(cd4.Tcells.flowD)
    cd8.Tcells <- getflowFrame(cd8.Tcells.flowD)
    
    all.events[19] <- nrow(cd4.Tcells)
    all.props[19] <- (nrow(cd4.Tcells)/nrow(Tcells))*100
    Filters.list$cd4.Tcells.filter <- cd4.Tcells.flowD@filter
    
    all.events[20] <- nrow(cd8.Tcells)
    all.props[20] <- (nrow(cd8.Tcells)/nrow(Tcells))*100
    Filters.list$cd8.Tcells.filter <- cd8.Tcells.flowD@filter
    
    ## Plotting CD8_CD4 to obtain CD4 T cells and CD8 T cells
    plotDens(Tcells, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "CD4 Tcells & CD8 Tcells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(h=cd4.gate.high, lwd=2)
    lines(cd4.Tcells.flowD.temp@filter, lwd = 2)
    lines(cd8.Tcells.flowD@filter, lwd = 2)
    
    if(centre != "ccp"){
      ## FMO CD25
      if(length(FMO.cd25.index) == 1 & FMO.cd25.flag != 0){
        cd8a.gate.FMO <- deGate(Tcells.FMO.cd25, channel = c(grep('CD8', FMO.cd25@parameters@data$desc)))
        cd4.gate.FMO <- deGate(Tcells.FMO.cd25, channel = c(grep('CD4', FMO.cd25@parameters@data$desc)))
        
        cd4.Tcells.flowD.FMO.cd25 <- flowDensity(Tcells.FMO.cd25, channels = c(grep('CD8', FMO.cd25@parameters@data$desc),grep('CD4', FMO.cd25@parameters@data$desc)), position = c(F,T), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        cd8.Tcells.flowD.FMO.cd25 <- flowDensity(Tcells.FMO.cd25, channels = c(grep('CD8', FMO.cd25@parameters@data$desc),grep('CD4', FMO.cd25@parameters@data$desc)), position = c(T,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO+0.5))
        
        # ## Plotting CD8_CD4 to obtain CD4 T cells and CD8 T cells
        # plotDens(Tcells.FMO.cd25, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "CD4 Tcells & CD8 Tcells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(cd4.Tcells.flowD.FMO.cd25@filter, lwd = 2)
        # lines(cd8.Tcells.flowD.FMO.cd25@filter, lwd = 2)
        
      }
      
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        cd8a.gate.FMO <- deGate(Tcells.FMO.cd44, channel = c(grep('CD8', FMO.cd44@parameters@data$desc)))
        cd4.gate.FMO <- deGate(Tcells.FMO.cd44, channel = c(grep('CD4', FMO.cd44@parameters@data$desc)))
        
        cd4.Tcells.flowD.FMO.cd44 <- flowDensity(Tcells.FMO.cd44, channels = c(grep('CD8', FMO.cd44@parameters@data$desc),grep('CD4', FMO.cd44@parameters@data$desc)), position = c(F,T), gates = c(cd8a.gate.FMO, cd4.gate.FMO))
        cd8.Tcells.flowD.FMO.cd44 <- flowDensity(Tcells.FMO.cd44, channels = c(grep('CD8', FMO.cd44@parameters@data$desc),grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,F), gates = c(cd8a.gate.FMO, cd4.gate.FMO+0.5))
        
        # ## Plotting CD8_CD4 to obtain CD4 T cells and CD8 T cells
        # plotDens(Tcells.FMO.cd44, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "CD4 Tcells & CD8 Tcells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(cd4.Tcells.flowD.FMO.cd44@filter, lwd = 2)
        # lines(cd8.Tcells.flowD.FMO.cd44@filter, lwd = 2)
        
      }
    }
    
    
    
    
    #########################################################################################
    ## Gating CD4 T-cells to obtain Tregs and T helper cells
    ## For UCD I am using CD25 vs CD44 instead of CD25 vs CD4
    
    if(centre == "ucd"){
      #cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.upper = T, upper = T, tinypeak.removal = 0.9)
      cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.percentile = T, percentile = 0.99)
      # if(cd25.gate.high < 1.75){
      #   cd25.gate.high <- max(c(deGate(cd4.Tcells, channel = c(channels.ind["CD25"])), deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), tinypeak.removal = 0.9)))
      #   
      # }else if(cd25.gate.high > 2.5){
      #   cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.percentile = T, percentile = 0.9)
      # }
      
      all.gthres[17] <- cd25.gate.high
      
      Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd25.gate.high, NA))
      
      
      T.helper <- cd4.Tcells
      T.helper@exprs <- T.helper@exprs[-Tregs.flowD@index,]
      
      Tregs <- getflowFrame(Tregs.flowD)
      
      all.events[21] <- nrow(Tregs)
      all.props[21] <- (nrow(Tregs)/nrow(cd4.Tcells))*100
      Filters.list$Tregs.filter <- Tregs.flowD@filter
      
      all.events[22] <- nrow(T.helper)
      all.props[22] <- (nrow(T.helper)/nrow(cd4.Tcells))*100
      
      
      ## Plotting CD25_CD4 to obtain Tregs and T helper cells
      # min.y <-  min(exprs(cd4.Tcells)[, c(channels.ind["CD44"])])-1
      # max.y <- max(exprs(cd4.Tcells)[, c(channels.ind["CD44"])])
      # min.x <- min(exprs(cd4.Tcells)[, c(channels.ind["CD25"])])
      # max.x <- max(exprs(cd4.Tcells)[, c(channels.ind["CD25"])])
      plotDens(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD44"]), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(Tregs.flowD@filter, lwd=2)
      
      
    }else {
      
      if(centre == "ccp"){
        cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.upper = T, upper = T, tinypeak.removal = 0.9)
        all.gthres[17] <- cd25.gate.high
        
        Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), gates = c(cd25.gate.high, NA))
        
      }else{
        cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.upper = T, upper = T, tinypeak.removal = 0.9)
        if(cd25.gate.high < 1.75){
          cd25.gate.high <- max(c(deGate(cd4.Tcells, channel = c(channels.ind["CD25"])), deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), tinypeak.removal = 0.9)))
          
        }else if(cd25.gate.high > 2.5){
          cd25.gate.high <- deGate(cd4.Tcells, channel = c(channels.ind["CD25"]), use.percentile = T, percentile = 0.9)
        }
        all.gthres[17] <- cd25.gate.high
        
        ## FMO CD25
        if(length(FMO.cd25.index) == 1 & FMO.cd25.flag != 0){
          ## Checking the CD25 marker gating threshold of FMO CD25
          if(cd4.Tcells.flowD.FMO.cd25@gates[1] > 2){
            Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), gates = c(cd25.gate.high, NA))
            
          }else{
            Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), use.control = c(T,F), control = c(cd4.Tcells.flowD.FMO.cd25, NA), use.percentile = c(T,F), percentile = 0.9955)
            
            if(round(Tregs.flowD@proportion) < 11){
              Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), use.control = c(T,F), control = c(cd4.Tcells.flowD.FMO.cd25, NA), use.percentile = c(T,F), percentile = 0.991)
            }else if(round(Tregs.flowD@proportion) >= 14){
              Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), use.control = c(T,F), control = c(cd4.Tcells.flowD.FMO.cd25, NA), use.percentile = c(T,F), percentile = 0.997)
              if(round(Tregs.flowD@proportion) >= 14){
                ## Not using FMOs anymore
                Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), gates = c(cd25.gate.high, NA))
                
                
              }
            }
            
          }
          
          
        }else{
          Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), position = c(T,NA), gates = c(cd25.gate.high, NA))
          
          
        }
        
      }
      
      T.helper <- cd4.Tcells
      T.helper@exprs <- T.helper@exprs[-Tregs.flowD@index,]
      
      Tregs <- getflowFrame(Tregs.flowD)
      
      all.events[21] <- nrow(Tregs)
      all.props[21] <- (nrow(Tregs)/nrow(cd4.Tcells))*100
      
      # Filters.list$Tregs.filter <- Tregs.flowD@filter
      
      all.events[22] <- nrow(T.helper)
      all.props[22] <- (nrow(T.helper)/nrow(cd4.Tcells))*100
      
      
      ## Plotting CD25_CD4 to obtain Tregs and T helper cells
      # min.y <-  min(exprs(cd4.Tcells)[, c(channels.ind["CD4"])])-1
      # max.y <- max(exprs(cd4.Tcells)[, c(channels.ind["CD4"])])
      # min.x <- min(exprs(cd4.Tcells)[, c(channels.ind["CD25"])])
      # max.x <- max(exprs(cd4.Tcells)[, c(channels.ind["CD25"])])
      plotDens(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(Tregs.flowD@filter, lwd=2)
      
      
    }
    
    if(centre != "ccp"){
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        numPeaks.cd25.FMO <- getPeaks(cd4.Tcells.flowD.FMO.cd44, channel = channels.ind["CD25"])
        if(length(numPeaks.cd25.FMO$Peaks) < 2){
          #FMO.cd44.flag <- 1
          FMO.cd44.flag <- 0
        }else{
          cd25.gate.FMO <- deGate(cd4.Tcells.flowD.FMO.cd44, channel = c(grep('CD25', FMO.cd44@parameters@data$desc)), use.percentile = T, percentile = 0.9)
          
          Tregs.flowD.FMO.cd44 <- flowDensity(cd4.Tcells.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(T,NA), gates = c(cd25.gate.FMO, NA))
          T.helper.flowD.FMO.cd44 <- flowDensity(cd4.Tcells.flowD.FMO.cd44, channels = c(grep('CD25', FMO.cd44@parameters@data$desc), grep('CD4', FMO.cd44@parameters@data$desc)), position = c(F,NA), gates = c(cd25.gate.FMO, NA))
          
          plotDens(cd4.Tcells.flowD.FMO.cd44, channels = c(channels.ind["CD25"], channels.ind["CD4"]), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
          lines(Tregs.flowD.FMO.cd44@filter, lwd=2)
          
        }      
        
      }
      
      
    }
    
    #########################################################################################
    ## Gating Tregs to obtain Effector Treg cells and Resting Treg cells
    
    cd62.gate.Tregs <- deGate(Tregs, channel = c(channels.ind["CD62L"]))
    if(centre != "ccp"){
      if(cd62.gate.Tregs > 1.75){
        cd62.gate.Tregs <- deGate(Tregs, channel = c(channels.ind["CD62L"]), use.percentile = T, percentile = 0.95)
        if(cd62.gate.Tregs > 1.75){
          cd62.gate.Tregs <- deGate(Tregs, channel = c(channels.ind["CD62L"]), all.cuts = T)
          if(length(cd62.gate.Tregs) > 1){
            #cd62.gate.Tregs <- min(deGate(Tregs, channel = c(channels.ind["CD62L"]), all.cuts = T))
            cd62.gate.Tregs <- mean(deGate(Tregs, channel = c(channels.ind["CD62L"]), all.cuts = T))
          }
          
        }
      }else if(cd62.gate.Tregs < 1){## changing this from <0.5 to <1
        cd62.gate.Tregs <- deGate(Tregs, channel = c(channels.ind["CD62L"]), use.upper=T, upper=T)
      }
    }
    
    all.gthres[18] <- cd62.gate.Tregs
    
    numPeaks.Tregs <- density(Tregs@exprs[,c(channels.ind["CD44"])])
    
    Tregs.peak.lcn <- findpeaks(numPeaks.Tregs$y)
    if(nrow(Tregs.peak.lcn) > 1){
      Tregs.secondPeak.lcn <- numPeaks.Tregs$x[Tregs.peak.lcn[ which.max(Tregs.peak.lcn[,1]),3]] ## Using findpeaks() to find the location where the max peak starts
    }else if(nrow(Tregs.peak.lcn) == 1){
      Tregs.secondPeak.lcn <- numPeaks.Tregs$x[Tregs.peak.lcn[1,3]] ## Using findpeaks() to find the location where the first peak starts
    }
    
    if(round(Tregs.secondPeak.lcn,1) > 1.5){
      Tregs.secondPeak.lcn <- deGate(Tregs, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.1)
      if(Tregs.secondPeak.lcn > 1.5){
        #Tregs.secondPeak.lcn <- numPeaks.Tregs$x[Tregs.peak.lcn[ which.max(Tregs.peak.lcn[,1]),3]] ## Using findpeaks() to find the location where the max peak starts
        Tregs.secondPeak.lcn <- 1.5
      }
    }
    
    if(centre == "ccp"){
      Effector.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F,T), gates = c(cd62.gate.Tregs, Tregs.secondPeak.lcn))
      
      Resting.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T,NA), gates = c(cd62.gate.Tregs,NA))
      
    }else{
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        #Effector.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F,T), use.control = c(F,T), control = c(NA, Tregs.flowD.FMO.cd44), gates = c(NA, Tregs.secondPeak.lcn))
        Effector.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F,T), use.control = c(F,T), control = c(NA, Tregs.flowD.FMO.cd44), use.upper = c(F,T), upper=c(NA,F))
        
        Resting.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T,NA),  use.control = c(F,T), control = c(NA, Tregs.flowD.FMO.cd44))
        
      }else{
        Effector.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(F,T), gates = c(cd62.gate.Tregs, Tregs.secondPeak.lcn))
        
        Resting.Tregs.flowD <- flowDensity(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), position = c(T,NA), gates = c(cd62.gate.Tregs,NA))
        
      }
    }
    
    
    Effector.Tregs <- getflowFrame(Effector.Tregs.flowD)
    
    
    Resting.Tregs <- getflowFrame(Resting.Tregs.flowD)
    
    
    all.events[23] <- nrow(Effector.Tregs)
    all.props[23] <- (nrow(Effector.Tregs)/nrow(Tregs))*100
    Filters.list$Effector.Tregs.filter <- Effector.Tregs.flowD@filter
    
    all.events[24] <- nrow(Resting.Tregs)
    all.props[24] <- (nrow(Resting.Tregs)/nrow(Tregs))*100
    Filters.list$Resting.Tregs.filter <- Resting.Tregs.flowD@filter
    
    ## Plotting CD62L_CD44 to obtain Effector Treg and Resting Treg cells
    plotDens(Tregs, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), main = "Effector & Resting Treg cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(Effector.Tregs.flowD@filter, lwd=2)
    lines(Resting.Tregs.flowD@filter, lwd=2)
    
    #########################################################################################
    ## Gating Tregs to obtain KLRG1+ Tregs cells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44 otherwise for other centres like CCP, I used KLRG1 vs CD4
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        klrg1.gate.Tregs <- deGate(Tregs, channel = c(channels.ind['KLRG1']))
        klrg1.Tregs.flowD <- flowDensity(Tregs, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate.Tregs, NA))
        
        klrg1.Tregs <- getflowFrame(klrg1.Tregs.flowD)
        
        
        all.events[33] <- klrg1.Tregs.flowD@cell.count
        all.props[33] <- klrg1.Tregs.flowD@proportion
        
        plotDens(Tregs, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ Tregs cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Tregs.flowD@filter, lwd=2)
        
      }else{
        
        all.events[33] <- NA
        all.props[33] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        klrg1.gate.Tregs <- deGate(Tregs, channel = c(channels.ind['KLRG1']))
        
        
        klrg1.Tregs.flowD <- flowDensity(Tregs, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate.Tregs, NA))
        
        klrg1.Tregs <- getflowFrame(klrg1.Tregs.flowD)
        
        
        all.events[33] <- klrg1.Tregs.flowD@cell.count
        all.props[33] <- klrg1.Tregs.flowD@proportion
        
        
        
        plotDens(Tregs, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ Tregs cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Tregs.flowD@filter, lwd=2)
        
      }else{
        
        all.events[33] <- NA
        all.props[33] <- NA
      }
      
    }
    
    
    
    
    #########################################################################################
    ## Gating T-helper cells to obtain Effector T-helper and Resting T-helper
    
    Effector.T.helper.flowD.temp <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,NA), gates = c(cd62.gate.Tregs, NA))
    Effector.T.helper.temp <- getflowFrame(Effector.T.helper.flowD.temp)
    
    numPeaks.cd44.gate.T.helper <- density(Effector.T.helper.temp@exprs[,c(channels.ind["CD44"])])
    
    cd44.gate.T.helper.peak.lcn <- findpeaks(numPeaks.cd44.gate.T.helper$y)
    if(nrow(Tregs.peak.lcn) == 1){
      cd44.gate.T.helper <- numPeaks.cd44.gate.T.helper$x[cd44.gate.T.helper.peak.lcn[1,3]] ## Using findpeaks() to find the location where the first peak starts
    }else if(nrow(Tregs.peak.lcn) > 1){
      cd44.gate.T.helper <- numPeaks.cd44.gate.T.helper$x[cd44.gate.T.helper.peak.lcn[1,4]] ## Using findpeaks() to find the location where the first peak ends
      
    }
    
    if(cd44.gate.T.helper > 2){
      cd44.gate.T.helper <- deGate(Effector.T.helper.temp, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.1)
    }else if(cd44.gate.T.helper < 2){
      cd44.gate.T.helper <- deGate(Effector.T.helper.temp, channel = c(channels.ind["CD44"]), use.upper = T, upper = F, tinypeak.removal = 0.9)
      if(cd44.gate.T.helper < 2){
        cd44.gate.T.helper <- deGate(Effector.T.helper.temp, channel = c(channels.ind["CD44"]))
        
      }
    }
    
    all.gthres[19] <- cd44.gate.T.helper
    
    if(centre == "ccp"){
      Effector.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate.T.helper))
      
      Resting.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate,NA))
      
    }else{
      ## FMO CD44
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        # Effector.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), use.control = c(F,T), control = c(NA, T.helper.flowD.FMO.cd44), gates = c(NA, cd44.gate.T.helper))
        
        Effector.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), use.control = c(F,T), control = c(NA, T.helper.flowD.FMO.cd44))
        Resting.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA),  use.control = c(F,T), control = c(NA, T.helper.flowD.FMO.cd44))
        
        if(Effector.T.helper.flowD@proportion < 10){
          Effector.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate.T.helper))
          Resting.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA),  use.control = c(F,T), control = c(NA, T.helper.flowD.FMO.cd44), gates = c(cd62.gate,NA))
          
        }
        
        
      }else{
        Effector.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate, cd44.gate.T.helper))
        
        Resting.T.helper.flowD <- flowDensity(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA), gates = c(cd62.gate,NA))
        
      }
    }
    
    
    Effector.T.helper <- getflowFrame(Effector.T.helper.flowD)
    
    Resting.T.helper <- getflowFrame(Resting.T.helper.flowD)
    
    
    all.events[25] <- nrow(Effector.T.helper)
    all.props[25] <- (nrow(Effector.T.helper)/nrow(T.helper))*100
    Filters.list$Effector.T.helper.filter <- Effector.T.helper.flowD@filter
    
    all.events[26] <- nrow(Resting.T.helper)
    all.props[26] <- (nrow(Resting.T.helper)/nrow(T.helper))*100
    Filters.list$Resting.T.helper.filter <- Resting.T.helper.flowD@filter
    
    ## Plotting CD62L_CD44 to obtain Effector T-helper and Resting T-helper cells
    plotDens(T.helper, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), main = "Effector & Resting T-helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(Effector.T.helper.flowD@filter, lwd=2)
    lines(Resting.T.helper.flowD@filter, lwd=2)
    
    
    
    #########################################################################################
    ## Gating T helper to obtain KLRG1+ T helper cells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        klrg1.gate.Thelper <- deGate(T.helper, channel = c(channels.ind['KLRG1']), use.upper = T, upper = T)
        klrg1.Thelper.flowD <- flowDensity(T.helper, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate.Thelper, NA))
        
        klrg1.Thelper <- getflowFrame(klrg1.Thelper.flowD)
        
        
        all.events[34] <- klrg1.Thelper.flowD@cell.count
        all.props[34] <- klrg1.Thelper.flowD@proportion
        
        plotDens(T.helper, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ T helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Thelper.flowD@filter, lwd=2)
        
      }else{
        
        all.events[34] <- NA
        all.props[34] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        klrg1.gate.Thelper <- deGate(T.helper, channel = c(channels.ind['KLRG1']), use.upper = T, upper = T)
        klrg1.Thelper.flowD <- flowDensity(T.helper, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate.Thelper, NA))
        
        klrg1.Thelper <- getflowFrame(klrg1.Thelper.flowD)
        
        
        all.events[34] <- klrg1.Thelper.flowD@cell.count
        all.props[34] <- klrg1.Thelper.flowD@proportion
        
        plotDens(T.helper, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ T helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Thelper.flowD@filter, lwd=2)     
        
      }else{
        
        all.events[34] <- NA
        all.props[34] <- NA
      }
      
    }
    
    
    
    #########################################################################################
    ## Gating CD8 T-cells to obtain CD8 Effector cells and CD8 Resting/Naive cells
    
    ## FMO CD44
    staticGate <- 0
    
    cd8.Tcells.flowD.temp <- flowDensity(cd8.Tcells, channels =  c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,NA),  gates = c(cd62.gate.Tregs, NA))
    if(centre == "ccp"){
      cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), tinypeak.removal = 0.1)
      
    }else{
      cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), tinypeak.removal = 0.1)-0.2
      if((cd44.gate.cd8.Tcells-3) > 0.1){
        cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.7)
        if((cd44.gate.cd8.Tcells-3) > 0.1){
          cd44.gate.cd8.Tcells <- 2.4 ## Changing this to 2.4 from 2.75
          staticGate <- 1
        }
      }
      
    }
    
    if(centre == "ccp"){
      Effector.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate.Tregs, cd44.gate.T.helper))
      
    }else{
      if(length(FMO.cd44.index) == 1 & FMO.cd44.flag !=0){
        Effector.cd8.flowD <- flowDensity(cd8.Tcells, channels =  c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), use.control = c(F,T), control = c(NA, cd8.Tcells.flowD.FMO.cd44),  gates = c(cd62.gate.Tregs, cd44.gate.T.helper))
        
      }else{
        Effector.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate.Tregs, cd44.gate.T.helper))
        
      }
    }
    
    
    
    ## For Resting CD8 and Naive CD8 T cells I am not using any FMO controls.
    
    Resting.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
    if(centre != "ccp"){
      if(Resting.cd8.flowD@gates[2] < cd44.gate.T.helper | Resting.cd8.flowD@proportion > 35 ){
        if(staticGate == 0){
          cd44.gate.cd8.Tcells <- max(deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), all.cuts = T))
          if(round(cd44.gate.cd8.Tcells,1) <= round(cd44.gate.T.helper,1)){
            cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.5)-0.1## Changing this from -0.1 to -0.2
            if(round(cd44.gate.cd8.Tcells,1) > 2.7){
              cd44.gate.cd8.Tcells <- min(2.9, deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.565))
              if(round(cd44.gate.cd8.Tcells,1) > 2.7){
                cd44.gate.cd8.Tcells <- 2.5
              }
            }
          }
          Resting.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
          
          
        }
      }else if((Resting.cd8.flowD@gates[2]-cd44.gate.T.helper) <= 0.1){
        cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), use.percentile = T, percentile = 0.565)
        if(round(cd44.gate.cd8.Tcells,1) <= round(cd44.gate.T.helper,1)){
          cd44.gate.cd8.Tcells <- deGate(cd8.Tcells.flowD.temp, channel = c(channels.ind["CD44"]), tinypeak.removal = 0.1)
          
        }
        Resting.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
        
      }
      
    }
    
    Naive.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,F), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
    
    # Plotting CD62L_CD44 to obtain CD8 Effector cells and CD8 Resting/Naive cells
    min.y <-  min(exprs(cd8.Tcells)[, c(channels.ind["CD44"])])-1
    max.y <- max(exprs(cd8.Tcells)[, c(channels.ind["CD44"])])
    if(min.y < 0 | max.y > 100){
      plotDens(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), main = "CD8 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v= cd62.gate.Tregs, lwd=2)
      
    }else{
      plotDens(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), main = "CD8 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y)); abline(v= cd62.gate.Tregs, lwd=2)
      
    }
    lines(Effector.cd8.flowD@filter, lwd = 2)
    lines(Resting.cd8.flowD@filter, lwd=2)
    lines(Naive.cd8.flowD@filter, lwd=2)
    # text(min(Effector.cd8.flowD@filter[, 1])+0.5, min(Effector.cd8.flowD@filter[, 2])-0.5, labels =  strcat('CD8 Effector cells \n', toString(signif(Effector.cd8.flowD@proportion, 3))), cex = 1.5)
    # text(max(Resting.cd8.flowD@filter[, 1])-0.5, max(Resting.cd8.flowD@filter[, 2])-0.5, labels =  strcat('CD8 Resting cells \n', toString(signif(Resting.cd8.flowD@proportion, 3))), cex = 1.5)
    # text(min(Naive.cd8.flowD@filter[, 1])+0.5, min(Naive.cd8.flowD@filter[, 2])-0.5, labels =  strcat('CD8 Naive cells \n', toString(signif(Naive.cd8.flowD@proportion, 3))), cex = 1.5)
    
    
    all.gthres[20] <- cd44.gate.cd8.Tcells
    
    Effector.cd8 <- getflowFrame(Effector.cd8.flowD)
    
    
    Resting.cd8 <- getflowFrame(Resting.cd8.flowD)
    
    
    Naive.cd8 <- getflowFrame(Naive.cd8.flowD)
    
    
    if (Resting.cd8.flowD@proportion < 5){
      flaggedFile[1] <- x$FCS.files
      #next()
    }else{
      flaggedFile[1] <- NA
    }
    
    
    all.events[27] <- nrow(Effector.cd8)
    all.props[27] <- (nrow(Effector.cd8)/nrow(cd8.Tcells))*100
    Filters.list$Effector.cd8.filter <- Effector.cd8.flowD@filter
    
    all.events[28] <- nrow(Resting.cd8)
    all.props[28] <- (nrow(Resting.cd8)/nrow(cd8.Tcells))*100
    Filters.list$Resting.cd8.filter <- Resting.cd8.flowD@filter
    
    all.events[29] <- nrow(Naive.cd8)
    all.props[29] <- (nrow(Naive.cd8)/nrow(cd8.Tcells))*100
    Filters.list$Naive.cd8.filter <- Naive.cd8.flowD@filter
    
    #########################################################################################
    ## Gating CD8 T cells to obtain KLRG1+ CD8 Tcells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        klrg1.gate.cd8.Tcells <- deGate(cd8.Tcells, channel = c(channels.ind['KLRG1']), use.upper = T, upper = T)
        klrg1.cd8.Tcells.flowD <- flowDensity(cd8.Tcells, channels =c(channels.ind["KLRG1"], channels.ind["CD44"]), position = c(T,NA), gates = c(klrg1.gate.cd8.Tcells, NA))
        
        klrg1.cd8.Tcells <- getflowFrame(klrg1.cd8.Tcells.flowD)
        
        
        all.events[35] <- klrg1.cd8.Tcells.flowD@cell.count
        all.props[35] <- klrg1.cd8.Tcells.flowD@proportion
        
        plotDens(cd8.Tcells, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD8 T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd8.Tcells.flowD@filter, lwd=2)
        
      }else{
        
        all.events[35] <- NA
        all.props[35] <- NA
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        klrg1.gate.cd8.Tcells <- deGate(cd8.Tcells, channel = c(channels.ind['KLRG1']), use.upper = T, upper = T)
        klrg1.cd8.Tcells.flowD <- flowDensity(cd8.Tcells, channels =c(channels.ind["KLRG1"], channels.ind["CD4"]), position = c(T,NA), gates = c(klrg1.gate.cd8.Tcells, NA))
        
        klrg1.cd8.Tcells <- getflowFrame(klrg1.cd8.Tcells.flowD)
        
        
        all.events[35] <- klrg1.cd8.Tcells.flowD@cell.count
        all.props[35] <- klrg1.cd8.Tcells.flowD@proportion
        
        plotDens(cd8.Tcells, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD8 T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd8.Tcells.flowD@filter, lwd=2)
        
      }else{
        
        all.events[35] <- NA
        all.props[35] <- NA
      } 
    }
    
    
    print("Gating finished")
    ##########################################################################################################################
    ##########################################################################################################################
    
    ## Start Plots-----
    
    if(is.na(flaggedFile[1])){
      png ( file = paste0(results.dir,"/Figures/ScatterPlots/", x$FCS.files, "_Automated.png"), width=2100, height=2100*6/4)
      #png ( file = paste0(results.dir,"/Figures/FailedFiles/", x$FCS.files, ".png"), width=2100, height=2100*6/4)
      
    }else{
      png ( file = paste0(results.dir,"/Figures/ScatterPlots/", x$FCS.files, "_Automated.png"), width=2100, height=2100*6/4)
      
    }
    #png ( file = paste0(results.dir,"/Live-Singlets-Figures/", x$FCS.files, ".png"), width=2100, height=2100*4/4)
    if(centre == "tcp"){
      par(mfrow=c(5,4),mar=(c(5, 5, 4, 2) + 0.1))
    }else{
      par(mfrow=c(6,4),mar=(c(5, 5, 4, 2) + 0.1))
    }
    
    plotDens(f, c(channels.ind["Time"], scat.chans['SSC-A']), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2)
    text(max(f@exprs[, c(channels.ind["Time"])])-25000, max(f@exprs[, c(scat.chans['SSC-A'])])-20000, labels =  strcat('All Events'), cex = 1.5)
    
    plotDens(f, c(scat.chans['FSC-A'], channels.ind["Live"]), main = "All Events, post flowCut", cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h=results$dead.peak, lwd=2)
    lines(live.flowD@filter)
    text(max(live.flowD@filter[, 1])-100000, max(live.flowD@filter[, 2]), labels =  strcat('Live \n', toString(signif(live.flowD@proportion, 3))), cex = 1.5)
    
    #text(max(f@exprs[, c(scat.chans['FSC-A'])])-18000, max(f@exprs[, c(channels.ind['Live'])])-1.5, labels =  strcat('Live \n', toString(signif(live.flowD@proportion, 3))), cex = 1.5)
    
    
    # Plot size gate
    plotDens(live.flowD@flow.frame, c(scat.chans['FSC-A'],scat.chans['SSC-A']), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(size.flowD@filter, lwd=2)
    text(max(f@exprs[, c(scat.chans['FSC-A'])])-100000, max(f@exprs[, c(scat.chans['SSC-A'])])-100000, labels =  strcat('Size \n 100.0'), cex = 1.5)
    
    
    #Plot FSC singlet gate
    plotDens(size.flowD@flow.frame, c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(FSCsinglets.flowD@filter)
    text(max(live.flowD@filter[, 1])-100000, max(FSCsinglets.flowD@filter[, 2])-10000, labels =  strcat('FSC Singlets\n', toString(signif(FSCsinglets.flowD@proportion, 3))), cex = 1.5)
    
    
    # Plot FSC singlet gate
    #if(!is.na(scat.chans['SSC-W'])){ 
    
    # col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = primary.colors, nbin = 1000)
    # plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    # lines(singlets.flowD@filter, col = 'blue')
    if(length(grep(colnames(f),pattern = "SSC|SS-*")) < 3){
      col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-A'], scat.chans['FSC-H'])], colramp = matlab.like2, nbin = 1000)
      plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-A'], scat.chans['FSC-H']), main = "SSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    }else{
      col <- densCols(getflowFrame(FSCsinglets.flowD)@exprs[, c(scat.chans['SSC-W'], scat.chans['SSC-H'])], colramp = matlab.like2, nbin = 1000)
      plotDens(FSCsinglets.flowD@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "SSC Singlets", cex.lab = 2, cex.axis = 2, cex.main=2, col = col)
    }
    lines(singlets.flowD@filter)
    text(max(FSCsinglets.flowD@filter[, 1])-200000, max(FSCsinglets.flowD@filter[, 2])-20000, labels =  strcat('SSC Singlets\n', toString(signif(singlets.flowD@proportion, 3))), cex = 1.5)
    
    #plotDens(FSC.singlets@flow.frame, channels = c(scat.chans['SSC-W'], scat.chans['SSC-H']), main = "FSC singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(temp2@filter)
    
    #}
    
    # ## COMMENT FROM HERE
    # ## Plotting TCRD_CD4 
    # plotDens(singlets.flowD, channels = c(channels.ind['TCRD'], channels.ind['CD4']), main = "TCRD", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    
    ## Plotting CD161_CD5 to obtain NK-cells
    plotDens(singlets.flowD, channels = c(channels.ind['CD161'], channels.ind['CD5']), main = "SSC Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(NKcells.flowD@filter, lwd=2)
    # text(max(NKcells.flowD@filter[, 1])+0.25, max(NKcells.flowD@filter[, 2])+0.25, labels =  strcat('NK cells\n', toString(signif((nrow(NKcells)/nrow(singlets))*100, 3))), cex = 1.5)
    
    
    ## Plotting CD62L_CD44 to obtain Effector NK-cells & Resting NK-cells
    plotDens(NKcells.flowD, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(v=cd62.gate, lwd= 2); #abline(h=NK.Effector.flowD@gates[2], lwd = 2)
    lines(NK.Effector.flowD@filter, lwd=2)
    lines(NK.Resting.flowD@filter, lwd=2)
    text(min(NK.Effector.flowD@filter[, 1])+1, min(NK.Effector.flowD@filter[, 2]), labels =  strcat('Effector NK \n', toString(signif(NK.Effector.flowD@proportion, 3))), cex = 1.5)
    text(min(NK.Resting.flowD@filter[, 1])+1, min(NK.Resting.flowD@filter[, 2])-1, labels =  strcat('Resting NK \n', toString(signif(NK.Resting.flowD@proportion, 3))), cex = 1.5)
    
    ## Plotting KLRG1_CD4 to obtain KLRG1+ CD4+ NK
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        plotDens(NKcells.flowD, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "NK cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=klrg1.gate, lwd=2)
        lines(klrg1.NK.flowD@filter, lwd=2)
        
      }
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        plotDens(NKcells.flowD, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "NK cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=klrg1.gate, lwd=2)
        lines(klrg1.NK.flowD@filter, lwd=2)
      }
      
    }
    
    
    ## Plotting CD25_CD5 to obtain CD5+ cells
    plotDens(not.NKcells, channels = c(channels.ind['CD25'], channels.ind['CD5']), main = "NOT NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); #abline(h= cd5.flowD@gates[2], lwd=2)
    lines(cd5.flowD@filter, lwd=2)
    text(max(cd5.flowD@filter[, 1])-0.5, max(cd5.flowD@filter[, 2])-2, labels =  strcat('CD5+ \n', toString(signif(cd5.flowD@proportion, 3))), cex = 1.5)
    
    ## Plotting CD161_CD8a to obtain CD161+ CD8- cells and NOT CD161+ CD8- cells
    plotDens(cd5.flowD, channels = c(channels.ind['CD161'], channels.ind['CD8']), main = "CD5+ cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd161pos.cd8neg.flowD@filter, lwd=2)
    text(max(cd161pos.cd8neg.flowD@filter[, 1]), max(cd161pos.cd8neg.flowD@filter[, 2])+0.5, labels =  strcat('CD161+ CD8- \n', toString(signif(cd161pos.cd8neg.flowD@proportion, 3))), cex = 1.5)
    
    
    ## Plotting CD25_CD4 to obtain CD4+ NKT cells and CD4- NKT cells
    if(length(which(names(channels.ind)=="GITR")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
      plotDens(cd161pos.cd8neg, channels = c(channels.ind['GITR'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd4pos.NKT.flowD@filter, lwd = 2)
      lines(cd4neg.NKT.flowD@filter, lwd = 2)
      text(max(cd4pos.NKT.flowD@filter[, 1])+0.25, max(cd4pos.NKT.flowD@filter[, 2]), labels =  strcat('CD4+ NKT \n', toString(signif(cd4pos.NKT.flowD@proportion, 3))), cex = 1.5)
      text(max(cd4neg.NKT.flowD@filter[, 1])+0.25, max(cd4neg.NKT.flowD@filter[, 2])-1, labels =  strcat('CD4- NKT \n', toString(signif(cd4neg.NKT.flowD@proportion, 3))), cex = 1.5)
      
    }else{
      if(is.na(cd25.gate.low)){
        plotDens(cd161pos.cd8neg, channels = c(channels.ind['CD161'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd4pos.NKT.flowD@filter, lwd = 2)
        lines(cd4neg.NKT.flowD@filter, lwd = 2)
        text(max(cd4pos.NKT.flowD@filter[, 1])+0.25, max(cd4pos.NKT.flowD@filter[, 2]), labels =  strcat('CD4+ NKT \n', toString(signif(cd4pos.NKT.flowD@proportion, 3))), cex = 1.5)
        text(max(cd4neg.NKT.flowD@filter[, 1])+0.25, max(cd4neg.NKT.flowD@filter[, 2])-1, labels =  strcat('CD4- NKT \n', toString(signif(cd4neg.NKT.flowD@proportion, 3))), cex = 1.5)
      }else{
        plotDens(cd161pos.cd8neg, channels = c(channels.ind['CD25'], channels.ind['CD4']), main = "NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd4pos.NKT.flowD@filter, lwd = 2)
        lines(cd4neg.NKT.flowD@filter, lwd = 2)
        text(max(cd4pos.NKT.flowD@filter[, 1])+0.25, max(cd4pos.NKT.flowD@filter[, 2]), labels =  strcat('CD4+ NKT \n', toString(signif(cd4pos.NKT.flowD@proportion, 3))), cex = 1.5)
        text(max(cd4neg.NKT.flowD@filter[, 1])+0.25, max(cd4neg.NKT.flowD@filter[, 2])-1, labels =  strcat('CD4- NKT \n', toString(signif(cd4neg.NKT.flowD@proportion, 3))), cex = 1.5)
        
      }
      
    }
    
    ## Plotting CD62L_CD44 to obtain Effector CD4+ NKT and Resting CD4+ NKT cells
    plotDens(cd4pos.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    lines(cd4pos.NKT.Effector.flowD@filter, lwd = 2)
    lines(cd4pos.NKT.Resting.flowD@filter, lwd = 2)
    text(min(cd4pos.NKT.Effector.flowD@filter[, 1])+0.75, min(cd4pos.NKT.Effector.flowD@filter[, 2])-0.25, labels =  strcat('Effector CD4+ NKT \n', toString(signif(cd4pos.NKT.Effector.flowD@proportion, 3))), cex = 1.5)
    text(max(cd4pos.NKT.Resting.flowD@filter[, 1])-0.5, max(cd4pos.NKT.Resting.flowD@filter[, 2]), labels =  strcat('Resting CD4+ NKT \n', toString(signif(cd4pos.NKT.Resting.flowD@proportion, 3))), cex = 1.5)
    
    
    
    ##########################################################################################
    
    # TEMPORARY COMMENT STARTS
    ## Plotting KLRG1_CD4 to obtain KLRG1+ CD4+ NKT
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        plotDens(cd4pos.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd4pos.NKT.flowD@filter, lwd=2)
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        plotDens(cd4pos.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd4pos.NKT.flowD@filter, lwd=2)
        #plotDens(cd4pos.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);(abline(v=klrg1.gate))
        
      }
      
    }
    
    ## Plotting CD62L_CD44 to obtain Effector CD4- NKT and Resting CD4- NKT cells
    plotDens(cd4neg.NKT, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD4- NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    lines(cd4neg.NKT.Effector.flowD@filter, lwd = 2)
    lines(cd4neg.NKT.Resting.flowD@filter, lwd = 2)
    text(min(cd4neg.NKT.Effector.flowD@filter[, 1])+0.75, min(cd4neg.NKT.Effector.flowD@filter[, 2])+0.75, labels =  strcat('Effector CD4- NKT \n', toString(signif(cd4neg.NKT.Effector.flowD@proportion, 3))), cex = 1.5)
    text(max(cd4neg.NKT.Resting.flowD@filter[, 1])-0.5, max(cd4neg.NKT.Resting.flowD@filter[, 2]), labels =  strcat('Resting CD4- NKT \n', toString(signif(cd4neg.NKT.Resting.flowD@proportion, 3))), cex = 1.5)
    
    ## Plotting KLRG1_CD4 to obtain KLRG1+ CD4- NKT
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        plotDens(cd4neg.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd4neg.NKT.flowD@filter, lwd=2)
      }
      
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        plotDens(cd4neg.NKT, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ NKT cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);#abline(v=klrg1.gate)
        lines(klrg1.cd4neg.NKT.flowD@filter, lwd=2)
      }
      
    }
    
    ## Plotting CD8a_CD4 to obtain T-cells
    plotDens(Tcells, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "NOT CD161+CD8-", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    text(max(Tcells@exprs[, c(channels.ind['CD8'])])-0.5, max(Tcells@exprs[, c(channels.ind['CD4'])])-0.5, labels =  strcat('Tcells \n', toString(signif((nrow(Tcells)/nrow(not.cd161pos.cd8neg))*100, 3))), cex = 1.5)
    
    ## Plotting CD8a_CD4 to obtain CD4 T-cells and CD8 T-cells
    plotDens(Tcells, channels = c(channels.ind['CD8'], channels.ind['CD4']), main = "T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    lines(cd4.Tcells.flowD@filter, lwd = 2)
    lines(cd8.Tcells.flowD@filter, lwd = 2)
    text(max(cd4.Tcells.flowD@filter[, 1]), max(cd4.Tcells.flowD@filter[, 2])+0.5, labels =  strcat('CD4 T-cells \n', toString(signif((nrow(cd4.Tcells)/nrow(Tcells))*100, 3))), cex = 1.5)
    text(max(cd8.Tcells.flowD@filter[, 1])-0.3, min(cd8.Tcells.flowD@filter[, 2])+0.5, labels =  strcat('CD8 T-cells \n', toString(signif(cd8.Tcells.flowD@proportion, 3))), cex = 1.5)
    
    
    ## Plotting CD25_CD4 to obtain Tregs and T helper cells
    ## For UCD I am using CD25 vs CD44 instead of CD25 vs CD4
    if(centre == "ucd"){
      
      plotDens(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD44"]), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(Tregs.flowD@filter, lwd=2)
      
    }else if(centre == "ccp"){
      plotDens(cd4.Tcells, channels = c(channels.ind["CD25"], channels.ind["CD4"]), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(Tregs.flowD@filter, lwd=2)
      
    }else if(length(which(names(channels.ind)=="GITR")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
      
      min.y <-  min(exprs(cd4.Tcells)[, c(channels.ind['CD4'])])-1
      max.y <- max(exprs(cd4.Tcells)[, c(channels.ind['CD4'])])
      plotDens(cd4.Tcells, channels = c(channels.ind['CD25'], channels.ind['GITR']), main = "CD4 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y)); #abline(v=cd4.Tcells.flowD.FMO.cd25@gates[2], lwd=2);
      lines(Tregs.flowD@filter, lwd=2)
      text(max(Tregs.flowD@filter[, 1])-1, min(Tregs.flowD@filter[, 2]), labels =  strcat('Tregs \n', toString(signif(Tregs.flowD@proportion, 3))), cex = 1.5)
      
    }
    
    
    ## Plotting CD62L_CD44 to obtain Effector Treg cells and Resting Treg cells
    min.y <-  min(exprs(Tregs)[, c(channels.ind['CD44'])])-1
    max.y <- max(exprs(Tregs)[, c(channels.ind['CD44'])])
    plotDens(Tregs, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y))
    lines(Effector.Tregs.flowD@filter)
    lines(Resting.Tregs.flowD@filter)
    text(min(Effector.Tregs.flowD@filter[, 1])+0.5, min(Effector.Tregs.flowD@filter[, 2])-0.3, labels =  strcat('Effector Treg cells \n', toString(signif(Effector.Tregs.flowD@proportion, 3))), cex = 1.5)
    text(min(Resting.Tregs.flowD@filter[, 1])+0.5, min(Resting.Tregs.flowD@filter[, 2])-0.3, labels =  strcat('Resting Treg cells \n', toString(signif(Resting.Tregs.flowD@proportion, 3))), cex = 1.5)
    
    
    ## Gating Tregs to obtain KLRG1+ Tregs cells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        plotDens(Tregs, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ Tregs cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Tregs.flowD@filter, lwd=2)
        
      }
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        plotDens(Tregs, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ Tregs cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Tregs.flowD@filter, lwd=2)
        
      }
    }
    
    
    ## Plotting CD62L_CD44 to obtain Effector T-helper and Resting T-helper
    min.y <-  min(exprs(T.helper)[, c(channels.ind['CD44'])])-1
    max.y <- max(exprs(T.helper)[, c(channels.ind['CD44'])])
    plotDens(T.helper, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "T-helper", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y))
    lines(Effector.T.helper.flowD@filter)
    lines(Resting.T.helper.flowD@filter)
    text(min(Effector.T.helper.flowD@filter[, 1])+0.5, min(Effector.T.helper.flowD@filter[, 2])-0.5, labels =  strcat('Effector T-helper \n', toString(signif(Effector.T.helper.flowD@proportion, 3))), cex = 1.5)
    text(min(Resting.T.helper.flowD@filter[, 1])+0.5, min(Resting.T.helper.flowD@filter[, 2])-0.5, labels =  strcat('Resting T-helper \n', toString(signif(Resting.T.helper.flowD@proportion, 3))), cex = 1.5)
    
    ## Gating T helper to obtain KLRG1+ T helper cells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        plotDens(T.helper, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD4+ T helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Thelper.flowD@filter, lwd=2)
        
      }
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        plotDens(T.helper, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD4+ T helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.Thelper.flowD@filter, lwd=2)
        
      }
    }
    
    
    # Plotting CD62L_CD44 to obtain CD8 Effector cells and CD8 Resting/Naive cells
    min.y <-  min(exprs(cd8.Tcells)[, c(channels.ind["CD44"])])-1
    max.y <- max(exprs(cd8.Tcells)[, c(channels.ind["CD44"])])
    if(min.y < 0 | max.y > 100){
      plotDens(cd8.Tcells, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD8 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v= cd62.gate.Tregs, lwd=2)
      
    }else{
      plotDens(cd8.Tcells, channels = c(channels.ind['CD62L'], channels.ind['CD44']), main = "CD8 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y)); abline(v= cd62.gate.Tregs, lwd=2)
      
    }
    lines(Effector.cd8.flowD@filter, lwd = 2)
    lines(Resting.cd8.flowD@filter, lwd=2)
    lines(Naive.cd8.flowD@filter, lwd=2)
    text(min(Effector.cd8.flowD@filter[, 1])+0.5, min(Effector.cd8.flowD@filter[, 2])-0.5, labels =  strcat('CD8 Effector cells \n', toString(signif(Effector.cd8.flowD@proportion, 3))), cex = 1.5)
    text(max(Resting.cd8.flowD@filter[, 1])-0.6, max(Resting.cd8.flowD@filter[, 2])-0.5, labels =  strcat('CD8 Resting cells \n', toString(signif(Resting.cd8.flowD@proportion, 3))), cex = 1.5)
    text(min(Naive.cd8.flowD@filter[, 1])+0.8, min(Naive.cd8.flowD@filter[, 2])+0.25, labels =  strcat('CD8 Naive cells \n', toString(signif(Naive.cd8.flowD@proportion, 3))), cex = 1.5)
    
    
    ## Gating CD8 T cells to obtain KLRG1+ CD8 Tcells, if the marker KLRG1 is present
    ## For UCD I am using KLRG1 vs CD44
    if(centre == "ucd"){
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD44")) == 1){
        
        plotDens(cd8.Tcells, channels = c(channels.ind['KLRG1'], channels.ind['CD44']), main = "CD8 T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd8.Tcells.flowD@filter, lwd=2)
        
      }
    }else{
      if(length(which(names(channels.ind)=="KLRG1")) == 1 && length(which(names(channels.ind)=="CD4")) == 1){
        
        plotDens(cd8.Tcells, channels = c(channels.ind['KLRG1'], channels.ind['CD4']), main = "CD8 T cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(klrg1.cd8.Tcells.flowD@filter, lwd=2)
        
      }
    }
    
    
    ##TEMPORARY COMMENT ENDS
    
    dev.off()
    par(mfrow=c(1,1))
    
    ## End Plots-----
    
  },error = function(err) {
    
    return(0)
  }
  ) # end of tryCatch
  
  list(x$FCS.files, all.props, all.events, all.gthres, all.fClean, Filters.list)
  #list(x$FCS.files, all.props, all.events, all.gthres, all.fClean, Filters.list, flaggedFile)
  
}, .parallel = TRUE) # end llply


print("Finished Gating & Plotting")

# all.props.store[i,] <- all.props
# all.events.store[i,] <- all.events

# all.props.store[i,12:ncol(all.props.store)] <- all.props
# all.events.store[i,12:ncol(all.events.store)] <- all.events
# all.gthres.store[i,12:ncol(all.gthres.store)] <- all.gthres
# all.filters.store[[i]] <- Filters.list


# 

# Saving the big dataframe of Proportions, Events, and Gating thresholds 
save(props.events.gates, file = paste0(results.dir,"/Props.Events.Gates.Rdata"))  

# Extracting the flowCut result from the large list 
f.Clean.store <- as.matrix(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][5])}))
#f.Clean.store <- props.events.gates[,ncol(props.events.gates)]
f.Clean.store <- as.matrix(f.Clean.store)
f.Clean.store <- cbind(store.allFCS[,c('Genotype', 'FCS files', 'Strain Code/Ear Tag/Mouse Number', 'Mouse ID/Animal ID', 'Colony ID/Line/Run Number','Assay Date', 'Gender', 'Number of Channels')], f.Clean.store)


# Extracting the FCS Files names from the large list 
all.FCS.store <- as.matrix(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][1])}))

# Extracting the Proportions from the large list & finding the files which failed the gating
all.props.store.temp <- as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][2])})))
all.props.store.temp <- cbind(all.FCS.store, all.props.store.temp)
if(centre == "tcp"){
  all.props.store.temp <- all.props.store.temp[,-c(31:36)]
}



failedGating.files.index.temp <- sapply(1:nrow(all.props.store.temp), function(x){which(is.na(all.props.store.temp[x,2:ncol(all.props.store.temp)]))})
failedGating.files.index <- which(sapply(1:length(failedGating.files.index.temp), function(x){length(failedGating.files.index.temp[[x]])})!=0)
failedGating.files <- store.allFCS[failedGating.files.index,]


save(failedGating.files.index, file = paste0(results.dir,"/failedGating.files.index.Rdata"))

all.props.store.temp <- all.props.store.temp[-failedGating.files.index,]
na.index <- which(is.na(all.props.store.temp[,c(2:ncol(all.props.store.temp))]))



all.props.store <- all.props.store.temp[,2:ncol(all.props.store.temp)]
rownames(all.props.store) <- 1:nrow(all.props.store)

if(centre == "tcp"){
  colnames(all.props.store) <- c("All Events", "Live", "FSC Singlets", "SSC Singlets", "NK cells", "NOT NK cells",
                                 "Effector NK cells", "Resting NK cells", "CD5+", "CD161+CD8-", "NOT CD161+CD8-",
                                 "CD4+ NKT cells", "CD4- NKT cells", "Effector CD4+ NKT cells", "Resting CD4+ NKT cells",
                                 "Effector CD4- NKT cells", "Resting CD4- NKT cells", "T-cells", "CD4 T-cells", "CD8 T-cells",
                                 "Tregs", "T-helper cells", "Effector Treg cells", "Resting Treg cells", "Effector T-helper cells",
                                 "Resting T-helper cells", "Effector CD8 T-cells", "Resting CD8 T-cells", "Naive CD8 T-cells")
  
}else{
  colnames(all.props.store) <- c("All Events", "Live", "FSC Singlets", "SSC Singlets", "NK cells", "NOT NK cells",
                                 "Effector NK cells", "Resting NK cells", "CD5+", "CD161+CD8-", "NOT CD161+CD8-",
                                 "CD4+ NKT cells", "CD4- NKT cells", "Effector CD4+ NKT cells", "Resting CD4+ NKT cells",
                                 "Effector CD4- NKT cells", "Resting CD4- NKT cells", "T-cells", "CD4 T-cells", "CD8 T-cells",
                                 "Tregs", "T-helper cells", "Effector Treg cells", "Resting Treg cells", "Effector T-helper cells",
                                 "Resting T-helper cells", "Effector CD8 T-cells", "Resting CD8 T-cells", "Naive CD8 T-cells",
                                 "KLRG1+ NK", "KLRG1+ CD4+ NKT", "KLRG1+ CD4- NKT", "KLRG1+ Tregs", "KLRG1+ T helper", "KLRG1+ CD8 T cells")
  
}

store.allFCS.temp <- store.allFCS[-failedGating.files.index,]
rownames(store.allFCS.temp) <- 1:nrow(store.allFCS.temp)
all.props.store <- cbind(store.allFCS.temp, all.props.store)


# # Extracting the Proportions from the large list & finding the files which failed the gating
props.events.gates.temp <- props.events.gates
# all.props.store.temp <- as.matrix(t(sapply(1:length(props.events.gates.temp), function(x){unlist(props.events.gates.temp[[x]][2])})))
# 
# colnames(all.props.store.temp) <- c("All Events", "Live", "FSC Singlets", "SSC Singlets", "NK cells", "NOT NK cells",
#                                "Effector NK cells", "Resting NK cells", "CD5+", "CD161+CD8-", "NOT CD161+CD8-",
#                                "CD4+ NKT cells", "CD4- NKT cells", "Effector CD4+ NKT cells", "Resting CD4+ NKT cells",
#                                "Effector CD4- NKT cells", "Resting CD4- NKT cells", "T-cells", "CD4 T-cells", "CD8 T-cells",
#                                "Tregs", "T-helper cells", "Effector Treg cells", "Resting Treg cells", "Effector T-helper cells",
#                                "Resting T-helper cells", "Effector CD8 T-cells", "Resting CD8 T-cells", "Naive CD8 T-cells",
#                                "KLRG1+ NK", "KLRG1+ CD4+ NKT", "KLRG1+ CD4- NKT", "KLRG1+ Tregs", "KLRG1+ T helper", "KLRG1+ CD8 T cells")
# 
# all.props.store.temp <- cbind(all.FCS.store, all.props.store.temp)
# all.props.store.temp <- cbind(failedGating.files[34:37,], all.props.store.temp)
# 
# for(i in 1:nrow(all.props.store.temp)){
#   index <- which(all.props.store.temp[i,c('FCS files')] == all.props.store[,c('FCS files')])
#   all.props.store[index, ] <- all.props.store.temp[i,]
# }



# Extracting the Event counts from the large list & finding the files which failed the gating
all.events.store.temp <-  as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][3])})))
all.events.store.temp <- cbind(all.FCS.store, all.events.store.temp)
if(centre == "tcp"){
  all.events.store.temp <- all.events.store.temp[,-c(31:36)]
}


all.events.store.temp <- all.events.store.temp[-failedGating.files.index,]
na.index <- which(is.na(all.events.store.temp[,c(2:ncol(all.events.store.temp))]))



all.events.store <- all.events.store.temp[,2:ncol(all.events.store.temp)]
rownames(all.events.store) <- 1:nrow(all.events.store)
if(centre == "tcp"){
  colnames(all.events.store) <- c("All Events", "Live", "FSC Singlets", "SSC Singlets", "NK cells", "NOT NK cells",
                                  "Effector NK cells", "Resting NK cells", "CD5+", "CD161+CD8-", "NOT CD161+CD8-",
                                  "CD4+ NKT cells", "CD4- NKT cells", "Effector CD4+ NKT cells", "Resting CD4+ NKT cells",
                                  "Effector CD4- NKT cells", "Resting CD4- NKT cells", "T-cells", "CD4 T-cells", "CD8 T-cells",
                                  "Tregs", "T-helper cells", "Effector Treg cells", "Resting Treg cells", "Effector T-helper cells",
                                  "Resting T-helper cells", "Effector CD8 T-cells", "Resting CD8 T-cells", "Naive CD8 T-cells")
  
}else{
  colnames(all.events.store) <- c("All Events", "Live", "FSC Singlets", "SSC Singlets", "NK cells", "NOT NK cells",
                                  "Effector NK cells", "Resting NK cells", "CD5+", "CD161+CD8-", "NOT CD161+CD8-",
                                  "CD4+ NKT cells", "CD4- NKT cells", "Effector CD4+ NKT cells", "Resting CD4+ NKT cells",
                                  "Effector CD4- NKT cells", "Resting CD4- NKT cells", "T-cells", "CD4 T-cells", "CD8 T-cells",
                                  "Tregs", "T-helper cells", "Effector Treg cells", "Resting Treg cells", "Effector T-helper cells",
                                  "Resting T-helper cells", "Effector CD8 T-cells", "Resting CD8 T-cells", "Naive CD8 T-cells",
                                  "KLRG1+ NK", "KLRG1+ CD4+ NKT", "KLRG1+ CD4- NKT", "KLRG1+ Tregs", "KLRG1+ T helper", "KLRG1+ CD8 T cells")
  
}
all.events.store <- cbind(store.allFCS.temp, all.events.store)

# # Extracting the Event counts from the large list & finding the files which failed the gating
# all.events.store.temp <- as.matrix(t(sapply(1:length(props.events.gates.temp), function(x){unlist(props.events.gates.temp[[x]][3])})))
# 
# all.events.store.temp <- cbind(all.props.store.temp[,1:11], all.events.store.temp)
# 
# for(i in 1:nrow(all.events.store.temp)){
#   index <- which(all.events.store.temp[i,c('FCS files')] == all.events.store[,c('FCS files')])
#   all.events.store[index, ] <- all.events.store.temp[i,]
# }
# 

#all.events.store[!is.na(all.events.store[,2]),]
# Extracting the Gating Thresholds from the large list 
all.gthres.store.temp <-  as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][4])})))
if(centre == "tcp"){
  all.gthres.store.temp <- all.gthres.store.temp[,-c(9, 21:23)]
}

colnames(all.gthres.store.temp) <- c("fsca.live.gate", "live.gate", "ss.high.gate", "cd5.gate", "cd161.gate", "cd161.gate.NK", "cd62.gate", "cd44.gate", 
                                     "cd8a.gate", "cd161.gate.high", "cd161.gate.low", "cd4.gate", "cd25.gate.low", "cd62.gate.NKT", "cd4.gate.high", "cd25.gate.high",
                                     "cd62.gate.Tregs", "cd44.gate.Thelper", "cd44.gate.cd8Tcells")
all.gthres.store.temp <- cbind(store.allFCS, all.gthres.store.temp)


all.gthres.store.temp <- all.gthres.store.temp[-failedGating.files.index,]
na.index <- which(is.na(all.gthres.store.temp[,c(2:ncol(all.gthres.store.temp))]))



all.gthres.store <- all.gthres.store.temp[,2:ncol(all.gthres.store.temp)]
rownames(all.gthres.store) <- 1:nrow(all.gthres.store)

# # Extracting the GThres from the large list & finding the files which failed the gating
# all.gthres.store.temp <- as.matrix(t(sapply(1:length(props.events.gates), function(x){unlist(props.events.gates[[x]][4])})))
# 
# all.gthres.store.temp <- cbind(all.props.store[,1:11], all.gthres.store.temp)
# 
# for(i in 1:nrow(all.gthres.store.temp)){
#   index <- which(all.gthres.store.temp[i,c('FCS files')] == all.gthres.store[,c('FCS files')])
#   all.gthres.store[index, ] <- all.gthres.store.temp[i,]
# }


# Extracting the Filters from the large list 
all.filters.store.temp <-  sapply(1:length(props.events.gates), function(x){props.events.gates[[x]][6]})
names(all.filters.store.temp) <- all.FCS.store


all.filters.store.temp <- all.filters.store.temp[-failedGating.files.index]
all.filters.store <- all.filters.store.temp

all.filters.store.temp <-  sapply(1:length(props.events.gates.temp), function(x){props.events.gates.temp[[x]][6]})
names(all.filters.store.temp) <- all.FCS.store

for(i in 1:nrow(all.gthres.store)){
  index <- which(all.gthres.store[i,c('FCS files')] == store.allFCS[,c('FCS files')])
  all.filters.store[[index]] <- all.filters.store.temp[[i]]
}



## Saving thresholds and filters in a list (will need them later for flowType)
Gates.Filter.list.new <- list()
for (i in 1:nrow(all.gthres.store)) {
  Gates.Filter.list.new[[i]] <- list(all.filters.store[i], all.gthres.store[i,c(12:ncol(all.gthres.store))])  
}

Gates.Filter.list.combined <- append(Gates.Filter.list, Gates.Filter.list.new)
save(all.props.store, file = paste0(results.dir,"/all.props.store.Rdata"))
save(all.events.store, file = paste0(results.dir,"/all.events.store.Rdata"))
save(all.gthres.store, file = paste0(results.dir,"/all.gthres.store.Rdata"))
save(all.filters.store, file = paste0(results.dir,"/all.filters.store.Rdata"))
save(Gates.Filter.list.combined, file = paste0(results.dir,"/Gates.Filter.list.combined.Rdata"))

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.props.store, file =  paste0(results.dir, "/DCCResults_Proportions_",toupper(centre),"_Panel1", date.time), row.names = FALSE)

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.events.store, file =  paste0(results.dir, "/DCCResults_EventCounts_",toupper(centre),"_Panel1", date.time), row.names = FALSE)

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(failedGating.files, file =  paste0(results.dir, "/Failed_Gating_Files",toupper(centre),"_Panel1", date.time), row.names = FALSE)

cat("Total time is: ",TimeOutput(start),sep="")

results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results/Figures/Flagged-ScatterPlots/")


file.rename(list.files(paste0(results.dir, "/Figures/FLagged-Scatterplots/", pattern=".fcs.png")), paste0(".fcs_Automated.png", 1:13))


filez <- list.files(results.dir)
tempCol <- sapply(1:length(filez), function(x){unlist(strsplit(filez[x], split = ".fcs.png"))})

file.rename(from=filez, to=sub(pattern=".fcs.png", replacement="_fcs.png", filez))
test <-"test"


### Processing the failed files
outliers.index <- sapply(1:nrow(failedFCS), function(j){which(failedFCS[j,c('FCS.files')] == store.allFCS[,c('FCS files')])})
outliers.store <- store.allFCS[outliers.index,]
# write.csv(outliers.store, file =  paste0(results.dir, "/outliers_Panel2"), row.names = FALSE)
failedFCS <- outliers.store[2:nrow(outliers.store),]
rownames(failedFCS) <- 1:nrow(failedFCS)
