# Written by Albina 
# Date: 2018-03-20

# NB: this script is keeping too many flowSets, which is a 'waste' of memory. Needs to work on this and re-write using gatingSets, but right now that is too much work.


remove(list=ls())

setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel2")

###############################################################################################
## NOTE: Run this section functions either separately before you run the rest of the script or you can run the entire script in just one go
## source("Codes/preProcessingMain.R"). The following two function prompts user for input: one for the centre name and the other panel number

source("helperFunc.R")
if (interactive() ){
  centre <- readCentreFunc()
  if(centre != "sanger"){
    if (interactive() ){
      panel <- readPanelFunc()
    }
  }
}

#############################################################################################

library('colorRamps')
library('plyr')
library('doMC')
library('e1071') # for flowCut
library('flowCore')
library('flowDensity')
library('pracma') # for findpeaks
library('tools')
library('MASS')
library('stringr')## for str_match used in compensateIMPC function
#library('flowViz')
library('flowCut')
library('Cairo')



#source("/data/Albina_IMPC-Universal/flowCut_20170331.R")
# source("/data/Albina_IMPC-Universal/getPeaks.R")
source("qualityGate.R")
# source('flowPrep.R')
# source('test.R')
start <- Sys.time()


results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")



load(paste0(results.dir,"/lgl.Rdata"))
load(paste0(results.dir,"/Genotype.Rdata"))
load(paste0(results.dir,"/uniqueGT.Rdata"))
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/store.allFMO.Rdata"))

if(panel == 2){
  store.allFMO <- store.allFMO[grep('PANEL_B', store.allFMO[,'FMO']), ]
}

colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
rownames(store.allFCS) <- 1:nrow(store.allFCS)
# scatterplot.dir <- '/data/projects/Codes/IMPC/scatterplots'
# ## Loading th Impress ID Table for naming the populations accordingly
# load('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/immpressIDconversion.Rdata')
# colnames(immpressIDconversion) <- c('Impress IDs', 'Population Name')

# Create directories
suppressWarnings(dir.create(paste0(results.dir,"/Cell_Counts")))
suppressWarnings(dir.create(paste0(results.dir,"/Cell_Proportions")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlots/")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/flowCut/")))
suppressWarnings(dir.create(paste0(results.dir,"/Figures/flowCut-FMO/")))

scatterplot.dir <- paste0(results.dir,"/Figures/ScatterPlots/")
#scatterplot.dir <- '/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel2/Results/Figures/ScatterPlots/'
# invisible(sapply(1:length(uniqueGT), function(x){
#   suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlots/", uniqueGT[x])))
# }))


file.names <- data.frame(store.allFCS, stringsAsFactors = F)
unique.run.inds <- as.matrix(unique(file.names[,'Panel.Organ.Folder']))

print("Starting Gating & Plotting")

# Panel 2 Issues coming up: On "2017-03-02"  "FMO_PANEL_B_REPEAT_CD11C_G09034.fcs" does not have the 
# channels labelled properly

startloop <- Sys.time()
# 187 <- no fmos

# run.inds <- c(9:186)
#run.inds <- c(2,18,46,70,100,120,140,160,180)
# Note: Data from 2016-12-13 looks bad. 
# run.inds <- c(1:187)

# # #use this code for checking for missing dates
# missingDates <- CheckmissingDate(allFileDates = unique(file.names[,'Panel.Organ.Folder']), OutDir = paste0(getwd(),"/results"))
# run.inds <- sapply(missingDates,function(x){ which(unique(file.names[,'Panel.Organ.Folder']) == x)})

#use this code to run outlier files
# plotFiles <- list.files(path="/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel2/Results/Figures/outlierPlots", pattern = ".png")
# run.inds <- unique(substr(plotFiles, 1,10))
# run.inds <- sapply(run.inds, function(x){which(unique(file.names[,'Panel.Organ.Folder'])== x)})

run.inds <- c(123, 186) #2018-01-11 files
run.inds <- c(4,5)
run.inds <- c(1)
# tempCol <- NA
# for(i in 1:nrow(store.allFCS)){
#   tempCol[i] <- gsub("/home/rstudio/results/IMPC/", "/home/rstudio/data/IMPC/", store.allFCS[i,1])
# }
# store.allFCS <- cbind(tempCol, store.allFCS)
# store.allFCS <-store.allFCS[,-2]
# colnames(store.allFCS)[1] <- "Path"
# 
# save(store.allFCS, file = paste0(results.dir,'/store.allFCS.Rdata'))

no_cores <- detectCores() - 10
registerDoMC(no_cores)

impressID.props.events.gates <- ldply(unique(file.names[,'Panel.Organ.Folder'])[run.inds], function(i1){ 
#impressID.props.events.gates <- ldply(unique(file.names[,'Panel.Organ.Folder']), function(i1){ 
  #impressID.props.events.gates <- ldply(1:nrow(file.names), function(i){ 
  
  # Changing the code to read in one day's worth of data at a time. This is because for TCP there are FMOs on a daily
  # basis so it makes sense to run the FMOs and all files on 1 day in 1 flowSet
  
  #x <- file.names[which(file.names[,'Panel.Organ.Folder'] == i1), ]
  x <- file.names[which(file.names[,'Panel.Organ.Folder'] == unique(file.names[,'Panel.Organ.Folder'])[run.inds]), ]
  
  print(i1)
  
  #all.gthres <- matrix(nrow = 1, ncol = 17, data = NA) # matrix for saving the gating thresholds
  
  
  try({
    
    if(centre == "sanger" | centre == "ciphe"){
      fs <- read.flowSet(paste0(x$Path, "/", x$FCS.files))
      fpath <- x$Path
    }else if(centre == "tcp" | centre == "bcm"| centre == "jax"){
      fs <- read.flowSet(paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)) # This assumes that all FCS files from that day can be put in a flowSet (same # channels ect)
      #fs <- read.ncdfFlowSet(paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)) 
      fpath <- paste0(x$Path,"/", x$Panel.Organ.Folder)
      FMO.index <- which(x$Assay.Date[1] == store.allFMO[,c("Assay Date")]) # This assumes all FCS files on 1 day have 1 set of FMOs
        
      
     }
    FCS.cellcount <- fsApply(fs, function(x) return(nrow(x)))
    gco <- gc(reset=T)
    channels.ind <- lapply(as(fs, Class = "list"), function(f){get.P2.channels.ind(f)})
    
    scat.chans <- c(grep(colnames(fs[[1]]),pattern = "FSC*"), grep(colnames(fs[[1]]),pattern = "SSC*"))
    names(scat.chans) <- colnames(fs[[1]])[scat.chans]
    
    # Remove scatter margins and compensate ---------------------------------------------
    
    # Removing margin events in Scatter channels
    fs <- fsApply(fs, function(f) removeMargins(f, chans = scat.chans, verbose = F))
    #Removing negative values in scatter channels
    fs <- fsApply(fs, function(f) removeMargins(f, chans = scat.chans, debris = T, neg = T, verbose = F))
    
    fs <- fsApply(fs, function(f) compensateIMPC(f, basename(x$FCS.files[which(x$FCS.files == basename(description(f)$FIL))]), fpath, centre = centre, panel.no = panel))
    
    fs <- fsApply(fs, function(f) transform(f, lgl))
    fs.preflowCut <- fs
    
    fcut <- fsApply(fs, function(f){
      channels.to.clean <- which(complete.cases(f@parameters@data$desc) == TRUE)
      f.Clean <- flowCut(f, Channels = channels.to.clean, MaxPercCut=0.5, FileID = identifier(f),
                         Directory = paste0(results.dir,"/Figures/flowCut/"), Plot = 'None')
      return(list(ind = f.Clean$ind, passed.flowCut = f.Clean$data['Has the file passed',]))
    })
    fs <- fsApply(fs.preflowCut, function(f){
      if(length(fcut[[which(x$FCS.files == basename(description(f)$FIL))]]$ind) > 0){
        f@exprs <- f@exprs[-fcut[[which(x$FCS.files == basename(description(f)$FIL))]]$ind, ]
      }
      return(f)
    })
    passed.flowCut <- lapply(fcut, function(x) return(x$passed.flowCut))
    gco <- gc(reset=T)
    
    #############################################################################################
    # Load and pre-process (compensate, transform, flowCut) all of the FMOs
    
    # Interesting: For TCP for the last FCS file, the channel names are all the same - ie same fluorophores. But different
    # order. So that wouldn't let me set up a flowFrame immediately.
    
    fs.temp <- list()
    for(iter in 1:length(fs)){
      fs.temp[[iter]] <- fs[[iter]]
    }
    channels.ind.list <- channels.ind
    counter <- length(fs) + 1
    # For the FMOS to make sense the markers need to be the same. This is just a check
    chans.name <- pData(parameters(fs[[1]]))$name[channels.ind[[1]]]
    #names(chans.name) <- names(channels.ind)
    names(chans.name) <- names(channels.ind[[1]])
  
    fmo.inds <- sapply(c("CD5","CD11B","CD11C","CD21","CD23","CD161","LY6C","MHCII"), function(fmo){ #marker names except CD19
        temp <- grep(x = store.allFMO[FMO.index, 'FMO'], pattern = fmo)
        return(ifelse(length(temp) == 0, yes = -1, no = temp))
      })
      
      
      if (length(which(fmo.inds==-1))<length(fmo.inds)){
        for(idx in fmo.inds){
          if (idx!=-1) { #if this is fmo file
            fs.temp[[counter]] <- preprocess.FMO.FCSfile(FMO.path = store.allFMO[FMO.index[idx],"Path"], FMO.assay.date = store.allFMO[FMO.index[idx],"Assay Date"], 
                                                         FMO.FCS.filename = store.allFMO[FMO.index[idx],"FMO"],fpath =  fpath, centre = centre, panel.no = panel.no, 
                                                         results.dir = results.dir, scat.chans = scat.chans, f.FMO.index = idx)
            fmo.channels.ind <- get.P2.channels.ind(fs.temp[[counter]])
            
            # ## Temporary Commented part
            # if(i1 =="2014-11-04"){ #fmo on this data does not have marker desc, force labelling channels
            #   if(length(fmo.channels.ind) < length(channels.ind.list[[1]])){
            #     if(all(fs.temp[[counter]]@parameters@data$name == fs.temp[[1]]@parameters@data$name)){
            #       fs.temp[[counter]]@parameters@data$desc <- fs.temp[[1]]@parameters@data$desc
            #       fmo.channels.ind <- get.P2.channels.ind(fs.temp[[counter]])
            #     }
            #   }
            # }
            
            chans.name.fmo <- pData(parameters(fs.temp[[counter]]))$name[fmo.channels.ind]
            names(chans.name.fmo) <- names(fmo.channels.ind)
            channels.ind.list[[counter]] <- fmo.channels.ind
            counter <- counter+1
          }
        }
        fs <- fs.temp
        names(fs) <- c(x$FCS.files, names(fmo.inds[which(fmo.inds > 0)]))
        names(channels.ind.list)[(length(fs.preflowCut)+1):length(channels.ind.list)] <- names(fmo.inds[which(fmo.inds > 0)])
      }
      
      if(is.null(names(fs))){
        fs <- fs.temp
        names(fs) <- c(x$FCS.files)
      }
      
      
      
      rm(fs.temp)
      # Panel 2 Issues coming up: On "2017-03-02"  "FMO_PANEL_B_REPEAT_CD11C_G09034.fcs" does not have the 
      # channels labelled properly. I'm just going to force rename the channels 
      if(i1 == "2017-03-02"){
        channels.ind.list['CD11C'] <- channels.ind.list[1]
      }
      
      # Albina added this part especially for the date 2014-11-04 (TCP data) since fmo on this date does not have marker desc
      if (length(which(fmo.inds==-1))<length(fmo.inds)){
        if(length(fmo.channels.ind) < length(channels.ind.list[[1]])){
          ## Because of FMO files problems, will only analyse the FCS files and skip the FMOs for the dates with problem
          fs.temp <- fs[-(length(fmo.inds)+1:length(fs))]
          fs <- fs.temp
          
          for(j in 1:length(fmo.inds)){
            fmo.inds[[j]] <- -1
          }
          rm(fs.temp)       
        }
      }
    
  
     #########################################################################################
    ## Panel2 Main Gating Starts
    
    #Quality Gate(Live/dead, singlets gating)
    
    results_QG <- lapply(1:length(fs), function(idx){
      f <- fs[[idx]]
      result <- qualityGate(f,scat.chans,channels.ind.list[[idx]],centre, panel.no = panel)
      return(result)
      
    })
    
    
    #plotting live
    #par(mfrow = c(3,3))
    
    # a<-lapply(1:length(fs),function(x1){
    #   plotDens(fs[[x1]], channels = c("FSC-A", "BV510-A"))
    #   lines((results_QG[[x1]]$live)@filter)
    # })
    
    live.filter <- lapply(results_QG,function(x1){
      return((x1$live)@filter)
    })
    
    fs.live <- lapply(results_QG,function(x1){
      return((x1$live)@flow.frame)
    })
    
    fs.live.proportion <- lapply(results_QG,function(x1){
      return((x1$live)@proportion)
    })
    
    #plotting FSC singlets 
    #par(mfrow = c(3,3))
    # a<-lapply(results_QG, function(x1){
    #   plotDens(x1$live,channels = c('FSC-A','FSC-H'))
    #   lines((x1$FSCsinglets)@filter)
    # })
    
    fs.fcssinglets <- lapply(results_QG,function(x1){
      return((x1$FSCsinglets)@flow.frame)
    })
    
    fs.fcssinglets.proportion <- lapply(results_QG,function(x1){
      return((x1$FSCsinglets)@proportion)
    })
    
    fcssinglets.filter <- lapply(results_QG,function(x1){
      return((x1$FSCsinglets)@filter)
    })
    
    #plotting SSC singlets 
    #par(mfrow = c(3,3))
    # b<-lapply(results_QG, function(x1){
    #   plotDens(x1$FSCsinglets,channels = c('SSC-W','SSC-H'))
    #   lines((x1$singlets)@filter)
    # })
    
    fs.singlets <- lapply(results_QG, function(x1){
      x1$singlets@flow.frame
    })
    
    fs.singlets.proportion <- lapply(results_QG, function(x1){
      x1$singlets@proportion
    })
    
    fs.singlets2 <- fs.singlets
    
    singlet.filter <-lapply(results_QG,function(x1){
      return((x1$singlets)@filter)
    }) 
    
    #########################################################################################
    ## If the marker F4/80 is present then use it to identify macrophages
    ## F4/80- and Macrophages
    
    F4marker.index <- which(names(channels.ind.list[[1]]) == "F4/80")
    if(length(F4marker.index) == 1){
      #f <- fs.singlets2[[idx]]
      
      F4.gate.list <- lapply(1:length(fs.singlets2), function(idx){
       
        f <- fs.singlets2[[idx]]
        
        tryCatch({
        
          
          f4.gate <- deGate(f, c(channels.ind.list[[idx]]["F4/80"]), tinypeak.removal = 0.001)
          
          if(f4.gate > 2.05 | f4.gate < 1){
            f4.gate <- 2
          }
          
          # print(f4.gate)
          # plotDens(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), main = paste0("Lymph: ", basename(description(f)$FIL)),
          #          cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=f4.gate, lwd=2)
          
          return(f4.gate)
        }, error=function(e){cat("Error in if :",conditionMessage(e), "\n")})
        
      })
      
      
      ## Macrophages
      macrophages <- lapply(1:length(fs.singlets2), function(idx){
        f <- fs.singlets2[[idx]]
        
        macrophages.flowD <- flowDensity(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), position = c(T, NA),
                                         gates = c(F4.gate.list[[idx]], NA))
        
        
        # plotDens(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), main = paste0("F4/80 Neg & Mac: ", basename(description(f)$FIL)),
        #          cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(macrophages.flowD@filter, type="l",lwd=2)
        
        return(macrophages.flowD)
      })
      fs.macrophages <- lapply(macrophages, function(x) return(x@flow.frame))
      names(fs.macrophages)<-names(fs)
      
      fs.macrophages.proportion <- lapply(macrophages, function(x) return(x@proportion))
      
      
      ## F4/80-
      f4.80neg <- lapply(1:length(fs.singlets2), function(idx){
        f <- fs.singlets2[[idx]]
        
        f4.80neg.flowD <- flowDensity(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), position = c(F, NA),
                                      gates = c(F4.gate.list[[idx]], NA))
        
        
        # plotDens(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), main = paste0("F4/80 Neg & Mac: ", basename(description(f)$FIL)),
        #          cex.lab = 2, cex.axis = 2, cex.main=2)
        # lines(f4.80neg.flowD@filter, type="l",lwd=2)
        
        return(f4.80neg.flowD)
      })
      fs.f4.80neg <- lapply(f4.80neg, function(x) return(x@flow.frame))
      names(fs.f4.80neg)<-names(fs)
      
      fs.f4.80neg.proportion <- lapply(f4.80neg, function(x) return(x@proportion))
      
    }
    
    
    
    #########################################################################################
    ## CD19-/+
    #par(mfrow = c(3,3))
    CD19.gate.list <- lapply(1:length(fs.singlets2), function(idx){
      print(idx)
      if(length(F4marker.index) == 1){
        f <- fs.f4.80neg[[idx]]
      }else{
        f <- fs.singlets2[[idx]]
      }
     
      
      tryCatch({
        
        
        cd19.gate <- deGate(f, c(channels.ind.list[[idx]]["CD19"]))
        # plotDens(f, c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), main = paste0("Lymph: ", basename(description(f)$FIL)),
        #          cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=cd19.gate, lwd=2)
        
        return(cd19.gate)
      }, error=function(e){cat("Error in if :",conditionMessage(e), "\n")})
      
    })
    
    
    ## CD19+
    #par(mfrow = c(3,3))
    cd19pos <- lapply(1:length(fs.singlets2), function(idx){
      
      if(length(F4marker.index) == 1){
        f <- fs.f4.80neg[[idx]]
      }else{
        f <- fs.singlets2[[idx]]
      }
      
      cd19 <- flowDensity(f, c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), position = c(T, NA),
                          gates = c(CD19.gate.list[[idx]], NA))
      
      
      plotDens(f, c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), main = paste0("CD19+: ", basename(description(f)$FIL)),
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd19@filter, type="l",lwd=2)
      
      return(cd19)
    })
    fs.19p <- lapply(cd19pos, function(x) return(x@flow.frame))
    names(fs.19p)<-names(fs)
    
    fs.19p.proportion <- lapply(cd19pos, function(x) return(x@proportion))
    
    
    ## CD19-
    cd19neg<- lapply(1:length(fs.singlets2), function(idx){
      
      if(length(F4marker.index) == 1){
        f <- fs.f4.80neg[[idx]]
      }else{
        f <- fs.singlets2[[idx]]
      }
      
      
      cd19 <- flowDensity(f, c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), position = c(F, NA),
                          gates = c(CD19.gate.list[[idx]], NA))
      
      plotDens(f, c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), main = paste0("CD19-: ", basename(description(f)$FIL)),
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd19@filter, type="l",lwd=2)
      return(cd19)
    })
    fs.19n <- lapply(cd19neg, function(x) return(x@flow.frame))
    names(fs.19n)<-names(fs)
    
    fs.19n.proportion <- lapply(cd19neg, function(x) return(x@proportion))
    
    #########################################################################################
    ## Gating CD19- to obtain T and NKT
    
    # Very simple deGate on the non-FMO population
    tnkt <- lapply(1:length(fs.19n), function(idx){
      f<-fs.19n[[idx]]
      if(fmo.inds["CD11B"]!=-1){
        cd11b.gate <- deGate(fs.19n[["CD11B"]], channels.ind.list[['CD11B']]["CD11b"], use.percentile = T, percentile = .999)+0.6
        if(cd11b.gate < 2){
          cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
        }
        # For TCP 2015-08-27, there is a small cluster of cells at CD11b ~3. So I need to chop off the
        # cluster before using the above gate
        if(cd11b.gate > 2.5){
          cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
          if(cd11b.gate > 3){
            cd11b.gate <- deGate(fs.19n[["CD11B"]], channels.ind.list[['CD11B']]["CD11b"], use.percentile = T, percentile = .999)+0.6
          }
        }
        
      }else{
        #cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.05)
        cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
      }
      
      if(fmo.inds["CD5"]!=-1){
        temp <- flowDensity(fs.19n[["CD5"]], c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                            position = c(NA, F), gates = c(NA, cd11b.gate))
        
        cd5.gate <- deGate(temp, channels.ind.list[['CD5']]["CD5"], tinypeak.removal = 0.0001)
        if(round(cd5.gate,2) >= 1.5 | cd5.gate < 1){
          temp <- flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                              position = c(NA, F), gates = c(NA, cd11b.gate))
          cd5.gate <- deGate(temp,channels.ind.list[[idx]]["CD5"], tinypeak.removal = 0.001)
          if(cd5.gate > 2){
            cd5.gate <- deGate(temp,channels.ind.list[[idx]]["CD5"], use.upper = T, upper = F)
            if(cd5.gate < 0.65){
              cd5.gate <- deGate(temp,channels.ind.list[[idx]]["CD5"])
              if(cd5.gate > 2){
                ## NOT using FMOs for this
                cd5.gate <- deGate(f,channels.ind.list[[idx]]["CD5"])
              }
            }
          }else if(cd5.gate < 0.55){
            cd5.gate <-1.75
          }
        }
        # maxDens.cd5.high <- density(getflowFrame(fs.19n[["CD5"]])@exprs[,c(channels.ind.list[[idx]]["CD5"])])
        # #plot(maxDens.cd5.high)
        # cd5.high.peak.lcn <- findpeaks(maxDens.cd5.high$y, npeaks = 2)
        # cd5.high.gate.temp1 <- maxDens.cd5.high$x[cd5.high.peak.lcn[which.max(cd5.high.peak.lcn[,1]),3]] ## Using findpeaks() to find the location where the highest peak ends
        
        
        
        #cd5.gate <- deGate(temp, channels.ind.list[['CD5']]["CD5"], use.upper = T, upper = T, alpha = 0.05)
        #cd5.gate <- deGate(temp, channels.ind.list[['CD5']]["CD5"], use.upper = T, upper = T, alpha = 0.005)+0.2
        
        temp1 <- flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                             position = c(T, F), gates = c(cd5.gate, cd11b.gate))
        
        # if(cd5.gate < 1.6){
        #   cd5.gate <- deGate(temp, channels.ind.list[['CD5']]["CD5"], use.upper = T, upper = T, alpha = 0.001)+0.2
        #  
        # }
        #check cell proportions
        if(temp1@proportion < 2){
          cd5.gate <-deGate(temp,channels.ind.list[[idx]]["CD5"], use.upper = T, upper=F, alpha=.08)+0.2
          if(cd5.gate < 1.5){
            cd5.gate <- deGate(temp, channels.ind.list[['CD5']]["CD5"], use.upper = T, upper = T, alpha = 0.001)
            
          }
        }
        rm(temp1)
        
      }else{
        temp <- flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                            position = c(NA, F), gates = c(NA, cd11b.gate))
        cd5.gate <-deGate(temp,channels.ind.list[[idx]]["CD5"], tinypeak.removal = 0.001)
        # cd5.gate <-deGate(temp,channels.ind.list[[idx]]["CD5"], use.upper= T, upper=F, alpha=.08)+0.2
        if(cd5.gate < 1.6){
          cd5.gate <- deGate(temp,channels.ind.list[[idx]]["CD5"])+0.2
        }
      }
      T.NKT<-flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                         position = c(T, F), gates = c(cd5.gate, cd11b.gate))
      plotDens(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]),
               main = paste0("CD19- : ", basename(description(f)$FIL)),
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(T.NKT@filter)
      if(T.NKT@proportion > 0){
        N.TNKT<-notSubFrame(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                            filter = T.NKT@filter)
        f <- getflowFrame(N.TNKT)
      }
     
      return(list(pop=T.NKT,frame=f))
    })
    
    fs.tnkt <- lapply(tnkt, function(x) return(x$pop@flow.frame))
    fs.tnkt.proportion <- lapply(tnkt, function(x) return(x$pop@proportion))
    fs.tnkt.flowD <- lapply(tnkt, function(x) return(x$pop))
    
    fs.not.tnkt <- lapply(tnkt, function(x) return(x$frame))
    tnkt.filter <- lapply(tnkt, function(x) return(x$pop@filter))
    rm(tnkt)
    
    names(fs.tnkt)<-names(fs.not.tnkt) <- names(fs)
    
    #####################################################################
    
    ## Gating NOT (T and NKT) cells to obtain NK cells.
    #par(mfrow = c(3,3))
    cd161 <- lapply(1:length(fs.not.tnkt), function(idx){
      f <- fs.not.tnkt[[idx]]
      
      if (fmo.inds["MHCII"]!=-1){
        #mhc.gate <- deGate(fs.not.tnkt[["MHCII"]], channels.ind.list[["MHCII"]]["MHCII"], tinypeak.removal = .9,upper=T,alpha=0.01)
        mhc.gate <- deGate(fs.not.tnkt[["MHCII"]], channels.ind.list[["MHCII"]]["MHCII"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.001)
      }else{
        mhc.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], use.upper = T, upper=T, tinypeak.removal = .9, alpha=.05)
      }
      
      if(fmo.inds["CD161"]!=-1){
        temp <- flowDensity(fs.not.tnkt[["CD161"]], c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),
                            position = c(NA, F), gates = c(NA, mhc.gate))
        cd161.gate <- deGate(temp, channels.ind.list[['CD161']]["CD161"], tinypeak.removal = .9, use.upper = T, upper=T, alpha=0.05)
      }else{
        temp <- flowDensity(f, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),
                            position = c(NA, F), gates = c(NA, mhc.gate))
        # cd161.gate <- deGate(temp, channels.ind.list[[idx]]["CD11b"], upper=T, tinypeak.removal = .9, alpha=.05)
        cd161.gate <- deGate(temp, channels.ind.list[[idx]]["CD161"])
        if(cd161.gate > 3 | cd161.gate < 2){
          numPeaks.cd161.gate <- density(f@exprs[,c(channels.ind.list[[idx]]["CD161"])])
          
          cd161.gate.peak.lcn <- findpeaks(numPeaks.cd161.gate$y)
          if(nrow(cd161.gate.peak.lcn) >= 3){
            cd161.gate.temp <- numPeaks.cd161.gate$x[cd161.gate.peak.lcn[,4]]
            cd161.gate.index <- which(cd161.gate.temp < 3 & cd161.gate.temp > 2)
            cd161.gate <- numPeaks.cd161.gate$x[cd161.gate.peak.lcn[cd161.gate.index,4]]
          }
    
        }
      }
      
      # For TCP the mhc gate sometimes looks like it is cutting off the CD161pos population so I'm going to try
      # To do something more complicated. Also, it depends on how "slanted" the mhcii+ populaiton gets. Is this
      # a compensation issue?
      
      #cd161.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),position = c(T,F),gates=c(cd161.gate,mhc.gate))
      #cd161.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),position = c(F,NA),gates=c(cd161.gate,mhc.gate))
      
      
      ## Commenting this part out temporarily.
      # if(!(names(fs.not.tnkt)[idx] %in% c('CD161'))){
      #   cd161.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),position = c(T,NA),gates=c(cd161.gate, NA))
      #   mhc.gate <- deGate(cd161.pos, channels.ind.list[[idx]]["MHCII"], use.upper = T, upper = T,tinypeak.removal = 0.9)
      #   if(mhc.gate > 3){
      #     mhc.gate <- deGate(cd161.pos, channels.ind.list[[idx]]["MHCII"])
      #     
      #   }
      # }
      cd161.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),position = c(T,F),gates=c(cd161.gate,mhc.gate))
      cd161.neg <- notSubFrame(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]), filter = cd161.pos@filter)
      
      plotDens(f, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),
               main = paste0("NOT(T and NKT): ", basename(description(f)$FIL)),
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd161.pos@filter)
      
      return(list(pop1=cd161.pos,pop2=cd161.neg))
    })
    
    fs.161p <- lapply(cd161, function(x) x$pop1@flow.frame)
    fs.161p.proportion <- lapply(cd161, function(x) x$pop1@proportion)
    fs.161p.flowD <- lapply(cd161, function(x) x$pop1)
    
    fs.161n <- lapply(cd161, function(x) x$pop2@flow.frame)
    fs.161n.flowD <- lapply(cd161, function(x) x$pop2)
    
    names(fs.161n)<- names(fs.161p.flowD) <- names(fs.161p)<- names(fs.161n.flowD) <- names(fs)
    
    
    ##################################################################################################
    # -------------------------------------------------------------------------------------
    # Myeloid cells ------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    
    # In certain cases, eg some of the  "2017-05-09" files, there is not a clear granulocyte cluster
    # In that case I'm putting the granulocyte CD5/Ly6G gate at the mean location
    #par(mfrow = c(3,3))
    
    cd5.granulo.gate <- lapply(1:length(fs.161n), function(idx){
      
      if(!(names(fs.161p)[idx] %in% c('CD11B'))){
        f <- fs.161n[[idx]]
        
        if(fmo.inds["CD5"]!=-1){
          cd5.gate <- deGate(fs.161n[["CD5"]], channels.ind.list[["CD5"]]["CD5"], tinypeak.removal = .9,  use.upper = T, upper=T, alpha=0.01)
        }else{
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"], tinypeak.removal = .9, alpha=.9)
        }
        
        temp <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                            position = c(NA, T), gates = c(NA, cd5.gate))
        if(temp@proportion < 5){
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"], alpha=0.9)
        }
        
        cd11b.gate <- deGate(temp, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F)
        granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                  position = c(T, T), gates = c(cd11b.gate, cd5.gate))
        #cd5.gate <- c(cd5.gate, deGate(granulo.fD, channels.ind.list[[idx]]["CD5"], all.cuts = T, upper = F))
        granulo.cd5.dens <- density(getflowFrame(granulo.fD)@exprs[, channels.ind.list[[idx]]["CD5"]])
        granulo.peak <- granulo.cd5.dens$x[which.max(granulo.cd5.dens$y)]
        cd5.gate2 <- 2*granulo.peak - deGate(granulo.fD,  channels.ind.list[[idx]]["CD5"], use.upper = T, upper = T, tinypeak.removal = 0.3)
        if(cd5.gate2 < (cd5.gate - 0.22)){
          cd5.gate <- NA
        }else{
          cd5.gate <- cd5.gate2
        }
        
      }else{
        cd5.gate <- NA
      }
      return(cd5.gate)
    })
    
    mean.cd5.granulo.gate <- mean(unlist(cd5.granulo.gate), na.rm = T)
    
    ## Adding this part for Jax but need to check if this works for TCP and other centres
    if(mean.cd5.granulo.gate < 2){
      index <- which(unlist(cd5.granulo.gate) > 2)
      mean.cd5.granulo.gate <- mean(unlist(cd5.granulo.gate)[index])
    }
    
    cd5.granulo.gate <- unlist(cd5.granulo.gate)
    cd5.granulo.gate[which(is.na(cd5.granulo.gate))] <- mean.cd5.granulo.gate
    
    
    cd11b.granulo.gate <- lapply(1:length(fs.161n), function(idx){
      
      if(!(names(fs.161p)[idx] %in% c('CD11B'))){
        f <- fs.161n[[idx]]
        
        if(fmo.inds["CD5"]!=-1){
          cd5.gate <- deGate(fs.161n[["CD5"]], channels.ind.list[["CD5"]]["CD5"], tinypeak.removal = .9,  use.upper = T, upper=T, alpha=0.01)
        }else{
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"],  use.upper = T, upper=F, tinypeak.removal = .9, alpha=.9)
        }       
        temp <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                            position = c(NA, T), gates = c(NA, cd5.gate))

         
       
          cd11b.gate2 <- deGate(temp, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F)
          
          granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate2, cd5.gate))
          
          cd5.gate <- cd5.granulo.gate[idx]
          
          granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(NA, T), gates = c(NA, cd5.gate))
          
          cd11b.gate <- deGate(granulo.fD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F, tinypeak.removal = 0.4, magnitude = 0.4)
          
          if(cd11b.gate == -Inf){
            cd11b.gate <- cd11b.gate2
          }else{
            peaks <- getPeaks(granulo.fD, channels.ind.list[[idx]]['CD11b'], tinypeak.removal = 0.9)
            
            if(cd11b.gate > min(peaks$Peaks)){
              cd11b.gate <- deGate(granulo.fD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F, twin.factor = 0.8)
            }
            
          }
          
        
               
       
      }else{
        cd11b.gate <- NA
      }
      
      return(cd11b.gate)
    })
    
    cd11b.granulo.gate <- unlist(cd11b.granulo.gate)
    inf.index <- which(is.infinite(cd11b.granulo.gate))
    gate.index <- which(cd11b.granulo.gate < 0.75)
    
    if(length(inf.index)!= 0 | length(gate.index) != 0){
      fs.temp <- fs[-(nrow(x)+1:length(fs))]
      fs <- fs.temp
      
      for(j in 1:length(fmo.inds)){
        fmo.inds[[j]] <- -1
      }
      
      fs.temp <- fs.19n[-(nrow(x)+1:length(fs.19n))]
      fs.19n <- fs.temp
      
      fs.temp <- fs.19p[-(nrow(x)+1:length(fs.19p))]
      fs.19p <- fs.temp
      
      fs.temp <- fs.tnkt[-(nrow(x)+1:length(fs.tnkt))]
      fs.tnkt <- fs.temp
      
      fs.temp <- fs.161n[-(nrow(x)+1:length(fs.161n))]
      fs.161n <- fs.temp
      
      fs.temp <- fs.161p[-(nrow(x)+1:length(fs.161p))]
      fs.161p <- fs.temp
      rm(fs.temp) 
      
      cd11b.granulo.gate <- cd11b.granulo.gate[-(nrow(x)+1:length(cd11b.granulo.gate))]
    }
    
    mean.cd11b.granulo.gate <- mean(unlist(cd11b.granulo.gate), na.rm = T)
  
    ## Adding this part for Jax but need to check if this works for TCP and other centres
    if(mean.cd11b.granulo.gate < 2.9){
      index <- which(unlist(cd11b.granulo.gate) > 2.9)
      mean.cd11b.granulo.gate <- mean(unlist(cd11b.granulo.gate)[index])
    }
    cd11b.granulo.gate[which(is.na(cd11b.granulo.gate))] <- mean.cd11b.granulo.gate
    cd11b.granulo.gate[which(cd11b.granulo.gate < 2)] <- mean.cd11b.granulo.gate
    
    counter <- 0
    ## Back tracking from T and NKT population gating. I will later find an efficient way of doing this part but for now will rewrite previous code with modification to cd5 gate 
    cd11b.granulo.gate.backtrack <- lapply(1:length(fs.19n), function(idx){
      
      cd11b.gate.backtrack <- cd11b.granulo.gate[idx]
      if(cd11b.gate.backtrack < 2){
        print("Start Backtracking")
        counter <- counter+1
        ## Gating CD19- to obtain T and NKT
        f<-fs.19n[[idx]]
        if(fmo.inds["CD11B"]!=-1){
          cd11b.gate <- deGate(fs.19n[["CD11B"]], channels.ind.list[['CD11B']]["CD11b"], use.percentile = T, percentile = .999)+0.6
          if(cd11b.gate < 2){
            cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
            if(cd11b.gate < 2){
              cd11b.gate <- 2
            }
          }
          # For TCP 2015-08-27, there is a small cluster of cells at CD11b ~3. So I need to chop off the
          # cluster before using the above gate
          if(cd11b.gate > 2.5){
            cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
            if(cd11b.gate > 3){
              cd11b.gate <- deGate(fs.19n[["CD11B"]], channels.ind.list[['CD11B']]["CD11b"], use.percentile = T, percentile = .999)+0.6
            }
          }
          
        }else{
          #cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.05)
          cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.01)
        }
        
        ## For CD5 gating, I am not using FMOs while Back Tracking
        
        temp <- flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                            position = c(NA, F), gates = c(NA, cd11b.gate))
        cd5.gate <-deGate(temp,channels.ind.list[[idx]]["CD5"], tinypeak.removal = 0.001)
        if(cd5.gate > 2.5){
          cd5.gate <-deGate(temp,channels.ind.list[[idx]]["CD5"], use.upper = T, upper = F)
        }
        T.NKT<-flowDensity(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                           position = c(T, F), gates = c(cd5.gate, cd11b.gate))
        plotDens(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]),
                 main = paste0("CD19- : ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(T.NKT@filter)
        N.TNKT<-notSubFrame(f, c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
                            filter = T.NKT@filter)
        fs.N.TNKT <- getflowFrame(N.TNKT)
        
        ## Gating NOT (T and NKT) cells to obtain NK cells. Not using FMOs
        
        mhc.gate <- deGate(fs.N.TNKT, channels.ind.list[[idx]]["MHCII"], use.upper = T, upper=T, tinypeak.removal = .9, alpha=.05)
        temp <- flowDensity(fs.N.TNKT, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),
                            position = c(NA, F), gates = c(NA, mhc.gate))
        cd161.gate <- deGate(temp, channels.ind.list[[idx]]["CD161"])
        
        cd161.POS<- flowDensity(fs.N.TNKT,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),position = c(T,F),gates=c(cd161.gate,mhc.gate))
        cd161.NEG <- notSubFrame(fs.N.TNKT,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]), filter = cd161.POS@filter)
        
        plotDens(fs.N.TNKT, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]),
                 main = paste0("NOT(T and NKT): ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd161.POS@filter)
        
        cd11b.gran.gate <- deGate(cd161.NEG, channels.ind.list[[idx]]['CD11b'], use.upper = T, upper = F, tinypeak.removal = 0.9, alpha = 0.9)
        cd5.gran.gate <- deGate(cd161.NEG, channels.ind.list[[idx]]["CD5"], alpha=0.9)
        # plotDens(cd161.neg, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']),
        #          main = paste0("NOT(NK): ", basename(description(f)$FIL)),
        #          cex.lab = 2, cex.axis = 2, cex.main=2)
        
      }else{
        T.NKT <- fs.tnkt.flowD[[idx]]
        fs.N.TNKT <- fs.not.tnkt[[idx]]
        cd161.POS <- fs.161p.flowD[[idx]]
        cd161.NEG <- fs.161n.flowD[[idx]]
        cd11b.gran.gate <- cd11b.granulo.gate[idx]
        cd5.gran.gate <- cd5.granulo.gate[idx]
        
      }
      
      
      return(list(pop1=T.NKT,frame=fs.N.TNKT, pop2=cd161.POS, pop3=cd161.NEG, gate1 = cd11b.gran.gate, gate2 = cd5.gran.gate, counter = counter))
    }) ## end of Back Tracking
    
    counter <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$counter))
    counter <- sum(unlist(counter))
    
    if(counter > 0){
      fs.tnkt <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$pop1@flow.frame))
      fs.tnkt.proportion <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$pop1@proportion))
      
      fs.not.tnkt <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$frame))
      tnkt.filter <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$pop1@filter))
      
      
      fs.161p <- lapply(cd11b.granulo.gate.backtrack, function(x) x$pop2@flow.frame)
      fs.161p.proportion <- lapply(cd11b.granulo.gate.backtrack, function(x) x$pop2@proportion)
      
      fs.161n <- lapply(cd11b.granulo.gate.backtrack, function(x) x$pop3@flow.frame)
      
      names(fs.tnkt)<-names(fs.not.tnkt) <-  names(fs.161n) <- names(fs.161p) <- names(fs)
      
      rm(cd11b.granulo.gate)
      cd11b.granulo.gate <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$gate1))
      
      mean.cd11b.granulo.gate <- mean(unlist(cd11b.granulo.gate), na.rm = T)
      if(mean.cd11b.granulo.gate < 3){
        mean.cd11b.granulo.gate <- mean(unlist(cd11b.granulo.gate[which(cd11b.granulo.gate > 3)]), na.rm = T)
      }
      cd11b.granulo.gate <- unlist(cd11b.granulo.gate)
      cd11b.granulo.gate[which(cd11b.granulo.gate < 3)] <- mean.cd11b.granulo.gate
      
      rm(cd5.granulo.gate)
      cd5.granulo.gate <- lapply(cd11b.granulo.gate.backtrack, function(x) return(x$gate2))
      
      mean.cd5.granulo.gate <- mean(unlist(cd5.granulo.gate), na.rm = T)
      cd5.granulo.gate <- unlist(cd5.granulo.gate)
      #cd5.granulo.gate[which(cd5.granulo.gate < 3)] <- mean.cd5.granulo.gate
      
    }
    
  
    granulo <- lapply(1:length(fs.161n), function(idx){
      
      if(!(names(fs.161p)[idx] %in% c('CD11B'))){
        f <- fs.161n[[idx]]
        
        cd5.gate <- cd5.granulo.gate[idx]
        if(is.na(cd5.gate) | cd5.gate < 1.5){
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"], alpha=0.9)
          if(cd5.gate < 2.25){
            cd5.gate <- 2.5
          }
        }
        
        cd11b.gate <- cd11b.granulo.gate[idx]
        granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                  position = c(T, T), gates = c(cd11b.gate, cd5.gate))
        if(granulo.fD@proportion < 0.5){
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"], alpha=0.9)
          granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate, cd5.gate))
          
        }
        granulo.fD <- flowDensity(granulo.fD, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                  position = c(T, T), gates = c(cd11b.gate, cd5.gate), ellip.gate = T, scale = 0.99)
              
        #if(min(granulo.fD@filter[,1]) < (mean.cd11b.granulo.gate - 0.4)){
        
        if(round(min(granulo.fD@filter[,1]),3) < round(mean.cd11b.granulo.gate - 0.6,3)){
          print('changing cd11b.gate')
          if(round(cd11b.gate,3) != 3.193){
            cd11b.gate <- mean.cd11b.granulo.gate
          }
          
          granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate, cd5.gate))
          granulo.fD <- flowDensity(granulo.fD, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate, cd5.gate), ellip.gate = T, scale = 0.99)
        }else if(granulo.fD@proportion > 90 & cd5.gate > 2){
          cd5.gate<- deGate(f, channel = channels.ind.list[[idx]]['CD5'])
          if(cd5.gate < 1){
            cd5.gate <- 2
          }
          granulo.fD <- flowDensity(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate, cd5.gate))
          granulo.fD <- flowDensity(granulo.fD, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                    position = c(T, T), gates = c(cd11b.gate, cd5.gate), ellip.gate = T, scale = 0.99)
                  
        }
        
        notgranulo <- getflowFrame(notSubFrame(f, channels =  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
                                               filter = granulo.fD@filter))
        
        
        plotDens(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']),
                 main = paste0("NOT(NK): ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(granulo.fD@filter)
      }else{
        granulo.fD <- NA
        notgranulo <- NA
      }
      
      return(list(granulo = granulo.fD, notgranulo = notgranulo))
    })
    
    fs.notgranulo <- lapply(granulo, function(x){ 
      if(!is.na(x)){
        return(x$notgranulo)
      }else{
        return(NA)
      }
    })
    
    granulo.filter <-lapply(granulo, function(x){
      if(!is.na(x)){
        return(x$granulo@filter)
      }else{
        return(NA)
      }
    })
    
    granulo.cellcount <-lapply(granulo, function(x){
      if(!is.na(x)){
        return(x$granulo@cell.count)
      }else{
        return(NA)
      }
    })
    
    granulo.proportion <- lapply(1:length(granulo.cellcount), function(x){(granulo.cellcount[[x]]/nrow(fs.161n[[x]])*100)})
    
    rm(granulo)
    gco <- gc(reset=T)
    names(fs.notgranulo) <-names(fs)
    
    ############################################################################################
    
    ## Gating Non Granulocytes to obtain Eosinophils
    #par(mfrow = c(3,3))
    eosino <- lapply(1:length(fs.notgranulo), function(idx){
      #print(idx)
      f <- fs.notgranulo[[idx]]
      
      if(!(names(fs.notgranulo)[idx] %in% c('CD11B'))){
        
        cd11b.gate <- cd11b.granulo.gate[idx]
        
        # if (fmo.inds["CD11B"]!=-1){
        #   cd11b.gate <- deGate(fs.161n[["CD11B"]], channels.ind.list[["CD11B"]]["CD11b"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.01)
        # }else{
        #   cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"],  use.upper = T, upper=F,tinypeak.removal = .9, alpha=.9)
        # }
        
        temp <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
                            position = c(T, NA), gates = c(cd11b.gate, NA))
        # ssch.gate <- c(deGate(temp, 'SSC-H',  use.upper = T, upper = T), deGate(temp, 'SSC-H',  use.upper = T, upper = T, all.cuts = T))
        # ssch.gate <- ssch.gate[which(ssch.gate > 50000)]
        # ssch.gate <- ssch.gate[which.min(abs(ssch.gate - 75000))]   
        ssch.gate <- deGate(f, 'SSC-H')
        # if(ssch.gate < 70000){
        #   ssch.gate <- deGate(f, 'SSC-H', tinypeak.removal = 0.001)
        # }
        if(ssch.gate < 50000){
          ssch.gate <- deGate(f, 'SSC-H', tinypeak.removal = 0.01)
          if(ssch.gate < 50000){
            ssch.gate <- deGate(f, 'SSC-H', use.percentile = T, percentile = 0.99)
          }
        }
        
        eosino.fD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
                                 position = c(T, T), gates = c(cd11b.gate-0.2, ssch.gate), ellip.gate = T)
        
        if(eosino.fD@proportion > 15){
          ssch.gate <- deGate(f, 'SSC-H', tinypeak.removal = 0.01)
          eosino.fD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
                                   position = c(T, T), gates = c(cd11b.gate-0.2, ssch.gate), ellip.gate = T)
          
        }else if(eosino.fD@proportion < 1.5 & ssch.gate > 130000){
          ssch.gate <- deGate(f, 'SSC-H', use.percentile = T, percentile = 0.97)
          eosino.fD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
                                   position = c(T, T), gates = c(cd11b.gate-0.2, ssch.gate), ellip.gate = T)
          
        }
        # cd11b.gate <- deGate(eosino.fD, channels.ind.list[[idx]]['CD11b'], use.upper = T, upper = F, tinypeak.removal = 0.5)
        # eosino.fD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
        #                          position = c(T, T), gates = c(cd11b.gate, ssch.gate))
        
        plotDens(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')),
                 main = paste0("NOT(granulocytes): ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(eosino.fD@filter, lwd =2)
        
        noteosino <- getflowFrame(notSubFrame(f, c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
                                              filter = eosino.fD@filter))
        
        
      }else{
        noteosino <- NA
        eosino.fD <- NA
      }
      
      return(list(eosino.fD = eosino.fD, noteosino = noteosino))
    })
    
    fs.eosino.cellcount <- lapply(eosino, function(x){
      if(!is.na(x)){
        return(x$eosino.fD@cell.count)
      }else{
        return(NA)
      }
    })
    
    fs.eosino.proportion <- lapply(eosino, function(x){
      if(!is.na(x)){
        return(x$eosino.fD@proportion)
      }else{
        return(NA)
      }
    })
    
    fs.noteosino <- lapply(eosino, function(x){ 
      if(!is.na(x)){
        return(x$noteosino)
      }else{
        return(NA)
      }
    })
    
    
    names(fs.noteosino) <-names(fs)
    
    
    
    ##################################################################################
    
    ##Gating Not Eosinophils to obtain Monocytes
    
    #par(mfrow = c(3,3))
    monocytes <- lapply(1:length(fs.noteosino), function(idx){
      f <- fs.noteosino[[idx]]
      
      if(!(names(fs.noteosino)[idx] %in% c('CD11B', 'LY6C'))){
        
        
        if (fmo.inds["LY6C"]!=-1){
          ly6c.gate <- deGate(fs.noteosino[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.01)
          if(ly6c.gate < 1.75){
            ly6c.gate <- deGate(fs.noteosino[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.00001)
            if(ly6c.gate < 1.9){
              ly6c.gate <- 2
            }
          }
        }else{
          ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"],tinypeak.removal = .8, alpha=.9)
          if(ly6c.gate < 1){
            ly6c.gate <- 2
          }else if(ly6c.gate > 2.5){
            ly6c.gate <- 2
          }
        }
        
        if (fmo.inds["CD11B"]!=-1){
          cd11b.gate <- deGate(fs.161p[["CD11B"]], channels.ind.list[["CD11B"]]["CD11b"], tinypeak.removal = .9,  use.upper = T, upper=T,alpha=0.01)
        }else{
          temp1 <- flowDensity(f,c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                               position = c(NA,T), gates=c(NA,ly6c.gate))
          cd11b.gate <- deGate(temp1, channels.ind.list[[idx]]["CD11b"], use.percentile = T, percentile = 0.25)
          if(cd11b.gate < 2.9) {
            cd11b.gate <- deGate(temp1, channels.ind.list[[idx]]["CD11b"], use.percentile = T, percentile = 0.35)
            
          }
        }
        
        temp1 <- flowDensity(f,c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                             position = c(NA,F), gates=c(NA,ly6c.gate))
        if(temp1@proportion < 12){
          temp1 <- flowDensity(f,c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                               position = c(NA,T), gates=c(NA,ly6c.gate))
          ly6c.gate <- deGate(temp1, channels.ind.list[[idx]]["Ly6C"],  use.upper = T, upper=F, alpha=.9)
          # temp1
        }
        rm(temp1)
        temp.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                                  position = c(T, T), gates = c(cd11b.gate, ly6c.gate))
        
        # plotDens(temp.flowD, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
        #          main = paste0("NOT(eosinophils): ", basename(description(f)$FIL)), 
        #          cex.lab = 2, cex.axis = 2, cex.main=2)
        
        cd11b.gate <- c(deGate(temp.flowD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F, tinypeak.removal = 0.9),
                        deGate(temp.flowD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F))
        cd11b.gate <- cd11b.gate[which.min(abs(cd11b.gate - 3))]
        
        if(cd11b.gate < 2.75){
          cd11b.gate <- deGate(temp.flowD, channels.ind.list[[idx]]['CD11b'])
        }
        
        monocytes.flowD <- flowDensity(f,  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                                       position = c(T, T), gates =  c(cd11b.gate, ly6c.gate))
        if(monocytes.flowD@proportion < 12){
          cd11b.gate <- c(deGate(temp.flowD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F, tinypeak.removal = 0.9),
                          deGate(temp.flowD, channels.ind.list[[idx]]['CD11b'],  use.upper = T, upper = F))
          cd11b.gate <- cd11b.gate[which.min(abs(cd11b.gate - 3))]
          if(cd11b.gate > 2){
            cd11b.gate <- 2
          }
          monocytes.flowD <- flowDensity(f,  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                                         position = c(T, T), gates =  c(cd11b.gate, ly6c.gate))
          
          
        }
        notmonocytes <- getflowFrame(notSubFrame(f,  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                                                 position = c(T, T), gates =  c(cd11b.gate, ly6c.gate)))
        
        plotDens(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
                 main = paste0("NOT(eosinophils): ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(monocytes.flowD@filter, lwd =2)
        abline(h = ly6c.gate)
        abline(v = cd11b.gate)
        
      }else{
        monocytes.flowD <- NA
        notmonocytes <- NA
      }
      
      return(list(monocytes.flowD = monocytes.flowD, notmonocytes = notmonocytes))
    })
    
    fs.notmono <- lapply(monocytes, function(x){ 
      if(!is.na(x)){
        return(x$notmonocytes)
      }else{
        return(NA)
      }
    })
    
    fs.monocytes.cellcount <- lapply(monocytes, function(x){ 
      if(!is.na(x)){
        return(x$monocytes.flowD@cell.count)
      }else{
        return(NA)
      }
    })
    
    fs.monocytes.proportions <- lapply(monocytes, function(x){ 
      if(!is.na(x)){
        return(x$monocytes.flowD@proportion)
      }else{
        return(NA)
      }
    })
    
    names(fs.notmono) <-names(fs)
    
    ##################################################################################################
    # ## Using CD317 to identify pDC if the marker CD317 is present.
    # 
    # CD317marker.index <- which(names(channels.ind.list[[idx]]) == "CD317")
    # if(length(CD317marker.index) == 1){
    #   pDCs <- lapply(1:length(fs.notmono), function(idx){
    #     f <- fs.notmono[[idx]]
    #     
    #     if(!(names(fs.noteosino)[idx] %in% c('CD11B', 'LY6C', 'MHCII','CD11C', 'Ly6c', 'CD317'))){
    #       
    #       # if (fmo.inds["MHCII"]!=-1 & !is.null(fs.notmono[["MHCII"]])){
    #       #   mhcii.gate <- deGate(fs.notmono[["MHCII"]], channels.ind.list[["MHCII"]]['MHCII'], tinypeak.removal = .9,  use.upper = T, upper=T,alpha=0.01)
    #       # }else{
    #       #   mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], tinypeak.removal = 0.9)
    #       #   if(mhcii.gate > 3){
    #       #     mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"])
    #       #   }
    #       # }
    #       
    #       ## Not using FMOs for this population
    #      
    #       cd317.gate <- deGate(f, channels.ind.list[[idx]]["CD317"], use.upper = T, upper = T)
    #       ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"],tinypeak.removal = 0.001)
    #       pDCs.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['Ly6C'], channels.ind.list[[idx]]['CD317']),
    #                             position = c(T, T), gates = c(ly6c.gate, cd317.gate), ellip.gate = F)
    #      
    #       plotDens(f, c(channels.ind.list[[idx]]['Ly6C'], channels.ind.list[[idx]]['CD317']),
    #                main = paste0("pDc: ", basename(description(f)$FIL)),
    #                cex.lab = 2, cex.axis = 2, cex.main=2);abline(v=ly6c.gate, lwd=2); abline(h=cd317.gate, lwd=2)
    #       lines(pDCs.flowD@filter, lwd=2)
    #     }else{
    #       pDCs.flowD <- NA
    #     }
    #     
    #     return(pDCs.flowD)
    #   })
    # }    
    # 
    # fs.pDCs <- lapply(pDCs, function(x){ 
    #   if(!is.na(x)){
    #     return(x@flow.frame)
    #   }else{
    #     return(NA)
    #   }
    # })
    # 
    # fs.pDCs.proportions <- lapply(pDCs, function(x){ 
    #   if(!is.na(x)){
    #     return(x@proportion)
    #   }else{
    #     return(NA)
    #   }
    # })
    # 
    # 
    # names(fs.pDCs) <-names(fs)
    # 
    ########################################################################################################
    
    ## Gating Non Monocytes to obtain Conventional Dendritic cells (cDC)
    #par(mfrow = c(3,3))
    cDCs <- lapply(1:length(fs.notmono), function(idx){
      f <- fs.notmono[[idx]]
      
      if(!(names(fs.noteosino)[idx] %in% c('CD11B', 'LY6C', 'MHCII','CD11C'))){
        
        # if (fmo.inds["MHCII"]!=-1 & !is.null(fs.notmono[["MHCII"]])){
        #   mhcii.gate <- deGate(fs.notmono[["MHCII"]], channels.ind.list[["MHCII"]]['MHCII'], tinypeak.removal = .9,  use.upper = T, upper=T,alpha=0.01)
        # }else{
        #   mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], tinypeak.removal = 0.9)
        #   if(mhcii.gate > 3){
        #     mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"])
        #   }
        # }
        
        ## Not using FMOs for this population
        flag <- 0
        mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], tinypeak.removal = 0.9)
        if(mhcii.gate > 2.85 | mhcii.gate < 1.75){
          mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"])
          if(mhcii.gate < 1.75){
            mhcii.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], all.cuts = T)
            if(length(mhcii.gate > 1)){
              flag <- 1
              mhcii.gate <- max(deGate(f, channels.ind.list[[idx]]["MHCII"], all.cuts = T))
            }
          }
        }
        
        if(flag == 1){
          cd11c.gate <- deGate(f, channels.ind.list[[idx]]["CD11c"], use.percentile = T, percentile = 0.85)
          if(cd11c.gate >= 2.75){
            cd11c.gate <- deGate(f, channels.ind.list[[idx]]["CD11c"], use.percentile = T, percentile = 0.65)
            if(cd11c.gate < 2.5){
              cd11c.gate <- 2.5
            }
          }
          cdcs <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                              position = c(T, T), gates = c(cd11c.gate, mhcii.gate), ellip.gate = T)
          flag <- 0
          
          if(cdcs@proportion < 22){
            cd11c.gate <- 1.9
            cdcs <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                                position = c(T, T), gates = c(cd11c.gate, mhcii.gate), ellip.gate = T)
          }
        }else{
          temp.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                                    position = c(NA, T), gates = c(NA, mhcii.gate))
          cd11c.gate <- deGate(temp.flowD, channels.ind.list[[idx]]["CD11c"])
          temp.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                                    position = c(T, T), gates = c(cd11c.gate, mhcii.gate))
          cd11c.gate <- deGate(temp.flowD, channels.ind.list[[idx]]["CD11c"],  use.upper = T, upper = F)
          
          if(is.infinite(cd11c.gate)){
            cd11c.gate <- deGate(f, channels.ind.list[[idx]]["CD11c"], use.percentile = T, percentile = 0.65)
            if(cd11c.gate < 2.15){
              cd11c.gate <- 2.5
            }
          }else if(cd11c.gate < 2.2){
            cd11c.gate <- deGate(temp.flowD, channels.ind.list[[idx]]["CD11c"],  use.percentile = T, percentile = 0.15)
            if(is.infinite(cd11c.gate)){
              cd11c.gate <- deGate(f, channels.ind.list[[idx]]["CD11c"])
              
            }else if(round(cd11c.gate,1) < 2.15){
              cd11c.gate <- 2.5
            }
          }else if(cd11c.gate > 3.5){
            # temp.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
            #                           position = c(NA, T), gates = c(NA, mhcii.gate))
            # cd11c.gate <- deGate(temp.flowD, channels.ind.list[[idx]]["CD11c"], tinypeak.removal = 0.1)
            # temp.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
            #                           position = c(T, T), gates = c(cd11c.gate, mhcii.gate))
            # cd11c.gate <- deGate(temp.flowD, channels.ind.list[[idx]]["CD11c"],  use.upper = T, upper = F, tinypeak.removal = 0.9)+0.2
            cd11c.gate <- 2.5
          }
           
          
          cdcs <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                              position = c(T, T), gates = c(cd11c.gate, mhcii.gate), ellip.gate = T)
          ## Changing this condition for JAX. Need to check if it works for TCP and other centres 
          #if(cdcs@proportion < 22 & cd11c.gate > 2.85){
          if(cdcs@proportion < 22){
            cd11c.gate <- 1.9
            cdcs <- flowDensity(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                                position = c(T, T), gates = c(cd11c.gate, mhcii.gate), ellip.gate = T)
          }
          
        }
        
        plotDens(f, c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
                 main = paste0("NOT(monocytes): ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cdcs@filter)
      }else{
        cdcs <- NA
      }
      
      return(cdcs)
    })
    
    
    fs.cDCs <- lapply(cDCs, function(x){ 
      if(!is.na(x)){
        return(x@flow.frame)
      }else{
        return(NA)
      }
    })
    
    fs.cDCs.proportions <- lapply(cDCs, function(x){ 
      if(!is.na(x)){
        return(x@proportion)
      }else{
        return(NA)
      }
    })
    
    
    names(fs.cDCs) <-names(fs)
    
    
    ##################################################################################################
    
    ## Gating cDCs to obtain CD11b+ cDC and CD11b- cDC
    
    #par(mfrow = c(3,3))
    #cd11b.cDCs.gate <- lapply(1:length(fs.cDCs), function(idx){
    cd11b.cDCs.gate <- lapply(1:nrow(x), function(idx){
      f <- fs.cDCs[[idx]]
      
      if(!(names(fs.cDCs)[idx] %in% c('CD11B', 'LY6C', 'MHCII','CD11C'))){
        
        # FMO doesn't really help here. The gate appears to be at a higher location in the universal gating strategy
        #cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], upper=F,tinypeak.removal =.9, alpha=.9)
        cd11b.peak <- flowDensity::getPeaks(f, channels.ind.list[[idx]]["CD11b"],tinypeak.removal = 0.2)$Peaks
        cd11b.peak <- cd11b.peak[which.min(abs(cd11b.peak - 3.25))]
        
        cd11b.kdedens <- kde2d.density(f,  c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']))
        #cd11b.gate <- 2*cd11b.peak - deGate(f, channels.ind.list[[idx]]["CD11b"], use.upper = T, upper = T, alpha = 0.2, tinypeak.removal = 0.2)
        cd11b.gate <- 2*cd11b.peak - deGate(cd11b.kdedens, use.upper = T, upper = T, alpha = 0.2, tinypeak.removal = 0.2)
        
        ## Changing this for Jax. Need to check if it works for TCp and other centres.
        #if(cd11b.gate < 2){
        if(cd11b.gate < 1.9){
          cd11b.gate <- NA
        }
      }else{
        cd11b.gate <- NA
      }
      return(cd11b.gate)
    })
    # Move to average location if needed
    mean.cd11b <- mean(unlist(cd11b.cDCs.gate), na.rm = T)
    
    #cd11b.cDCs <- lapply(1:length(fs.cDCs), function(idx){
    cd11b.cDCs <- lapply(1:nrow(x), function(idx){  
      f <- fs.cDCs[[idx]]
      cd11b.gate <- cd11b.cDCs.gate[[idx]]
      
      if(abs(cd11b.cDCs.gate[[idx]] - mean.cd11b) > 0.4 | is.na(cd11b.cDCs.gate[[idx]])){
        cd11b.gate <- mean.cd11b
      }
      if(!is.na(f)){
        cd11b.cDCs.pos.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']),
                                            position = c(T, NA), gates = c(cd11b.gate, NA))
        if(cd11b.cDCs.pos.flowD@proportion < 35){
          cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], all.cuts = T)
          if(length(cd11b.gate) > 1){
            cd11b.gate <- cd11b.gate[which(cd11b.gate < 2.5)]
          }
          
          cd11b.cDCs.pos.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']),
                                              position = c(T, NA), gates = c(cd11b.gate, NA))
          
        }
        cd11b.cDCs.neg.flowD <- flowDensity(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']),
                                            position = c(F, NA), gates = c(cd11b.gate, NA))
        plotDens(f, c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']),
                 main = paste0("cDCs: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2);abline(v = cd11b.gate, lwd=2)
        lines(cd11b.cDCs.pos.flowD@filter, lwd = 2)
        lines(cd11b.cDCs.neg.flowD@filter, lwd = 2)
      }else{
        cd11b.cDCs.pos.flowD <- NA
        cd11b.cDCs.neg.flowD <- NA
      }
      
      
      
      return(list(cd11b.cDCs.pos.flowD = cd11b.cDCs.pos.flowD, cd11b.cDCs.neg.flowD = cd11b.cDCs.neg.flowD))
    })
    
    fs.cd11b.cDCs.pos.cellcount <- lapply(cd11b.cDCs, function(x){ 
      if(!is.na(x)){
        return(x$cd11b.cDCs.pos.flowD@cell.count)
      }else{
        return(NA)
      }
    })
    
    fs.cd11b.cDCs.pos.proportions <- lapply(cd11b.cDCs, function(x){ 
      if(!is.na(x)){
        return(x$cd11b.cDCs.pos.flowD@proportion)
      }else{
        return(NA)
      }
    })
    
    fs.cd11b.cDCs.neg.cellcount <- lapply(cd11b.cDCs, function(x){ 
      if(!is.na(x)){
        return(x$cd11b.cDCs.neg.flowD@cell.count)
      }else{
        return(NA)
      }
    })
    
    fs.cd11b.cDCs.neg.proportions <- lapply(cd11b.cDCs, function(x){ 
      if(!is.na(x)){
        return(x$cd11b.cDCs.neg.flowD@proportion)
      }else{
        return(NA)
      }
    })
    
    
    ######################################################################################
    ######################################################################################
    
    ## Gating CD19+: B cells core panel with optional CD23
    # This will assume all FCS files on a given day have CD23 or not, otherwise will fail
    #par(mfrow = c(3,3))
    if(length(grep(pattern = "CD23",x = names(channels.ind[[1]])))>0){
      cd19.161<- lapply(1:length(fs.19p), function(idx){
        f <- fs.19p[[idx]]
        
        if(fmo.inds["CD161"]!=-1){
          cd161.gate <- deGate(fs.19p[["CD161"]], channels.ind.list[['CD161']]["CD161"], tinypeak.removal = .9,use.upper = T, upper=T, alpha=0.05)
        }else{
          cd161.gate <- deGate(f, channels.ind.list[[idx]]["CD161"], use.upper = T, upper=T, tinypeak.removal = .9)
        }
        
        cd161.pos <-flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,T),gates=c(NA,cd161.gate)) 
        cd161.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,F),gates=c(NA,cd161.gate))
        
        if(cd161.pos@proportion >= 1.45){
          if(cd161.gate < 1.45){
            cd161.pos <-flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,T),gates=c(NA,cd161.gate+0.5)) 
            
            cd161.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,F),gates=c(NA,cd161.gate+0.5))
            
          }else{
            cd161.pos <-flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,T),gates=c(NA,cd161.gate+0.3)) 
            
            cd161.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD161"]),position = c(NA,F),gates=c(NA,cd161.gate+0.3))
            
          }
          
        }
        
        plotDens(f, c( channels.ind.list[[idx]]["CD19"],channels.ind.list[[idx]]["CD161"]),
                 main = paste0("CD19+: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd161.neg@filter)
        
        return(list(pop1=cd161.neg))
      })
    }
    fs.19.n161 <- lapply(cd19.161,function(x) x$pop1@flow.frame)
    fs.19.n161.proportion <- lapply(cd19.161,function(x) x$pop1@proportion)
    
    names(fs.19.n161)<-names(fs)
    
    
    ##########################################################################
    
    
    
    # For these days I need to remove debris for the B cell panel to gate properly
    # If I had time, I'd check if this works on all files and then apply to all, but I don't.
    # Also doing this here is not good because it messes up the CD161+CD19+ to CD19+ ratio (ie the 
    # former doesn't contain noise but the later does)
    if((i1 ==  "2016-07-28")|(i1 == '2016-07-14')){
      fs.19.n161 <- lapply(fs.19.n161,function(x){
        x@exprs <- x@exprs[which(x@exprs[, channels.ind.list[["CD21"]]["CD21/CD35"]] < 4.45), ]
        return(x)
      })   
    }
    names(fs.19.n161)<-names(fs)
    
    # outlier_csv <- read.csv("Outlier_cellprops.csv")
    # 
    # outlier_nonmz <- outlier_csv[which(outlier_csv[,'marker']=="Non-MZ or FO - fraction of Live"),]
    # 
    # outlier_nonmz2018 <- outlier_nonmz[which(outlier_nonmz[,'Panel.Organ.Folder']=="2018-01-11"),]
    # rm(outlier_csv,outlier_nonmz)
    
    
    
    # Sometimes the location of the 2nd CD23 gate is not obvious from the 1D density so I'll move the ones that are 
    # far from the mean to the mean location.
    # Also added in a check for the 1st CD23 gate location
    
    
    ####################################################################
    
    ## Gating CD19+CD161- to obtain MZ B, FO B, and Non-FO/ MZ B
    
    fs.rot <- lapply(1:length(fs.19.n161),function(idx){
      f <- rotate.data(fs.19.n161[[idx]],c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),theta = -pi/5)$data
      # if(names(fs.19.n161[idx]) %in% outlier_nonmz2018[,1]){
      #   f <- fs.19.n161[[idx]]
      # }
      return(f)
    }) 
    
    names(fs.rot)<-names(fs.19.n161)
    
    
    # plotDens(fs.rot[[idx]], c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
    #          main = paste0("CD19+CD161-: ", basename(description(f)$FIL)),
    #          cex.lab = 2, cex.axis = 2, cex.main=2)
    # abline(v=cd23.gate1.mean)
    # abline(v=cd23.gate2.mean[[idx]])
    
    cd23.gate1 <-  lapply(1:length(fs.19.n161), function(idx){
      if(!(names(fs.rot)[idx] %in% c('CD23', 'CD21'))){
        f <- fs.rot[[idx]]
        cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], alpha=.05, tinypeak.removal = 1/40)
        cd23.gate<-ifelse(names(fs.19.n161)[idx]=="CD23",yes =deGate(f, channels.ind.list[["CD23"]]["CD23"], use.upper = T, upper=T, tinypeak.removal = .99) ,no = cd23.gate)
      }else{
        cd23.gate <- NA
      }
      return(cd23.gate = cd23.gate)
    })
    
    cd23.gate1.mean <- mean(unlist(cd23.gate1), na.rm = T)
    
    cd23.gate2.mean <-  lapply(1:length(fs.19.n161), function(idx){
      if(!(names(fs.rot)[idx] %in% c('CD23', 'CD21'))){
        f <- fs.rot[[idx]]
        cd23.gate <- cd23.gate1[[idx]]
        if(abs(cd23.gate - cd23.gate1.mean) > 0.5){
          cd23.gate <- min(deGate(f, channels.ind.list[[idx]]["CD23"], all.cuts = T, alpha=.05, tinypeak.removal = 1/40))
          cd23.gate <- ifelse(names(fs.19.n161)[idx]=="CD23", yes = min(deGate(f, channels.ind.list[["CD23"]]["CD23"], use.upper = T, upper=T, all.cuts = T, tinypeak.removal = .99)),no = cd23.gate)
          if(abs(cd23.gate - cd23.gate1.mean) > 0.5){
            cd23.gate <- cd23.gate1.mean
          }
          cd23.gate1[[idx]] <- cd23.gate
        }
        mzb <- rotate.fd(flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                     position = c(F,F),gates=c(cd23.gate, NA),use.upper=c(F,T),upper=c(NA,T)),angle = -pi/5)
        nonmzb <- rotate.fd(flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                        position = c(T,NA),gates=c(cd23.gate, NA)),angle = -pi/5) 
        f.rot <-rotate.data(getflowFrame(nonmzb), c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),theta = pi/7)$data
        cd23.gate <-deGate(f.rot, channels.ind.list[[idx]]["CD23"], use.upper = T, upper = F, tinypeak.removal = 0.8, alpha = 0.2)
        
        
      }else{
        cd23.gate <- NA
      }
      return(cd23.gate = cd23.gate)
    })
    
    cd23.gate2.mean <- mean(unlist(cd23.gate2.mean), na.rm = T)
    
    MZB <- lapply(1:length(fs.19.n161), function(idx){
      
      if(!(names(fs.rot)[idx] %in% c('CD23', 'CD21'))){
        f <- fs.rot[[idx]]
        cd23.gate <- cd23.gate1[[idx]]
        
        
        mzb <- rotate.fd(flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                     position = c(F,F),gates=c(cd23.gate, NA),use.upper=c(F,T),upper=c(NA,T)),angle = -pi/5)
        
        # plotDens(mzb, c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
        #          main = paste0("CD19+CD161-: ", basename(description(f)$FIL)),
        #          cex.lab = 2, cex.axis = 2, cex.main=2)
        # abline(v=cd23.gate1.mean)
        # abline(v=cd23.gate2.mean)
        
        if(mzb@proportion < 1.35){ #1.4 is an optimal proportion, do not change this number. 
          #write code for not rotating if proportions of mzb is less than 1.4. 
          f <- fs.19.n161[[idx]]
          #cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], all.cuts = T, upper = F)
          cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], all.cuts = T, use.upper = T, upper = F)
          cd23.gate <- cd23.gate[which.min(abs(cd23.gate2.mean - cd23.gate))]
          if(cd23.gate < (cd23.gate2.mean - 0.35)){
            cd23.gate <- cd23.gate2.mean - 0.5
          }
          
          #cd21.gate <- deGate(f, channels.ind.list[[idx]]["CD21/CD35"], upper=F, alpha=.1, tinypeak.removal = 1/40)
          cd21.gate <- deGate(f, channels.ind.list[[idx]]["CD21/CD35"],use.upper = T, upper=F, alpha=.1, tinypeak.removal = 1/40)
          
          mzb <- flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                             position = c(F,T),gates=c(cd23.gate, cd21.gate))
          
          nonmzb <- notSubFrame(f,c(channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                position = c(F,F), filter=mzb@filter)
          
          FOB<-  flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                             position = c(T,T),gates=c(cd23.gate,cd21.gate))
          
          nonMZB.FOB <- notSubFrame(nonmzb,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                    position = c(NA,F),filter=FOB@filter)
          
          # plotDens(fs.19.n161[[idx]], c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
          #          main = paste0("CD19+CD161-: ", basename(description(f)$FIL)),
          #          cex.lab = 2, cex.axis = 2, cex.main=2)
          # lines(mzb@filter)
          # lines(FOB@filter)
          # lines(nonMZB.FOB@filter)
          
          
        }else{
          if(mzb@proportion > 20 | cd23.gate > -0.5){
            
            #cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], all.cuts = T, upper = F)
            cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"])
            if(cd23.gate > 0.1){
              cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], use.upper = T, upper = F)
              if(cd23.gate < -0.85){
                cd23.gate <- deGate(f, channels.ind.list[[idx]]["CD23"], use.percentile = T, percentile = 0.1)
              }
              
            }
            
            #cd21.gate <- deGate(f, channels.ind.list[[idx]]["CD21/CD35"],use.upper = T, upper=F, alpha=.1, tinypeak.removal = 1/40)
            
            
            
            mzb <- rotate.fd(flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                         position = c(F,F),gates=c(cd23.gate, NA),use.upper=c(F,T),upper=c(NA,T)),angle = -pi/5)
            
          }
          
          nonmzb <- rotate.fd(flowDensity(f,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                          position = c(T,NA),gates=c(cd23.gate, NA)),angle = -pi/5) 
          f.rot <-rotate.data(getflowFrame(nonmzb), c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),theta = pi/7)$data
          
          #cd23.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD23"], all.cuts = T, upper = F, tinypeak.removal = 0.8, alpha = 0.2)
          cd23.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD23"], all.cuts = T, tinypeak.removal = 0.8, alpha = 0.2)
          # plotDens(f.rot, c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
          #          main = paste0("CD19+CD161-: ", basename(description(f)$FIL)),
          #          cex.lab = 2, cex.axis = 2, cex.main=2)
          # abline(v=cd23.gate)
          
          if(length(cd23.gate)==1){
            #cd23.gate2 <- deGate(f.rot, channels.ind.list[[idx]]["CD23"], all.cuts = T, upper = F)
            cd23.gate2 <- deGate(f.rot, channels.ind.list[[idx]]["CD23"], all.cuts = T, use.upper = T, upper = F)
            cd23.gate2 <- cd23.gate2[which.min(abs(cd23.gate2.mean - cd23.gate2))]
            if(abs(cd23.gate2.mean-cd23.gate) < abs(cd23.gate2.mean-cd23.gate2)){
              cd23.gate <- cd23.gate
            }else{
              cd23.gate <- cd23.gate2
            }
            
          }else{
            cd23.gate <- cd23.gate[which.min(abs(cd23.gate2.mean - cd23.gate))]
          }
          
          if(cd23.gate < (cd23.gate2.mean - 0.35) | cd23.gate > 3.25){
            cd23.gate <- cd23.gate2.mean
            if(cd23.gate < 2.5 ){
              cd23.gate <- 2.5
            }
          }
          
          cd21.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD21/CD35"], use.upper = T, upper=F, alpha=.1, tinypeak.removal = 1/40)
          if(cd21.gate < 0){
            cd21.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD21/CD35"], alpha=.1, tinypeak.removal = 0.025)+0.25
            
          }
          
          FOB<-  rotate.fd(flowDensity(f.rot,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                       position = c(T,T),gates=c(cd23.gate,cd21.gate)),angle = pi/7) 
          FOB.prop <- nrow(getflowFrame((FOB)))/nrow(f)*100
          
          if(FOB.prop < 0.5 & cd23.gate > 2.5 & cd21.gate > 2.5){
            cd23.gate <- 2.5
            cd21.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD21/CD35"], use.upper = T, upper=F, alpha=.1, tinypeak.removal = 1/40)
            FOB<-  rotate.fd(flowDensity(f.rot,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                         position = c(T,T),gates=c(cd23.gate,cd21.gate)),angle = pi/7) 
            FOB.prop <- nrow(getflowFrame((FOB)))/nrow(f)*100
            
          }else if (FOB.prop < 1 & cd23.gate < 1){
            cd23.gate <- 1.75
            cd21.gate <- deGate(f.rot, channels.ind.list[[idx]]["CD21/CD35"], use.upper = T, upper=F, alpha=.1, tinypeak.removal = 1/40)
            
            FOB<-  rotate.fd(flowDensity(f.rot,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                         position = c(T,T),gates=c(cd23.gate,cd21.gate)),angle = pi/7) 
            FOB.prop <- nrow(getflowFrame((FOB)))/nrow(f)*100
            
            
          }else if(FOB.prop < 20){
            cd23.gate <- 1.9
            FOB<-  rotate.fd(flowDensity(f.rot,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                         position = c(T,T),gates=c(cd23.gate,cd21.gate)),angle = pi/7)
            FOB.prop <- nrow(getflowFrame((FOB)))/nrow(f)*100
            
          }
          nonMZB.FOB<- notSubFrame(nonmzb,c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                                   position = c(F,NA),filter=FOB@filter)
          
          plotDens(fs.19.n161[[idx]], c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]),
                   main = paste0("CD19+CD161-: ", basename(description(f)$FIL)),
                   cex.lab = 2, cex.axis = 2, cex.main=2)
          lines(mzb@filter)
          lines(FOB@filter)
          lines(nonMZB.FOB@filter)
        }
        
        
      }else{
        mzb <- NA
        FOB <- NA
        nonMZB.FOB <- NA
      }
      return(list(mzb=mzb, follicular=FOB,nonfomzb= nonMZB.FOB))
    })
    fs.nmzb<-lapply(MZB, function(x){
      if(!is.na(x)){
        return(x$nonfomzb@flow.frame)
      }else{
        return(NA) 
      }
    })
    names(fs.nmzb)<- names(fs)
    
    
    #######################################################
    
    ## Gating Non-FO/MZB to obtain B1a using MHCII marker
    B1a<- lapply(1:length(fs.nmzb), function(idx){
      
      if(!(names(fs.nmzb)[[idx]]  %in% c('CD23', 'CD21'))){
        f<- fs.nmzb[[idx]]
        if(fmo.inds["MHCII"]!=-1) {
          
          mhc.peaks <- getPeaks(f,c(channels.ind.list[[idx]]["MHCII"]), tinypeak.removal = 0.9)
          mhc.gate.temp <- deGate(fs.nmzb[["MHCII"]], channels.ind.list[["MHCII"]]["MHCII"], tinypeak.removal = .9, use.upper = T, upper=T, alpha=0.01)
          mhc.gate <- mean(c(mhc.gate.temp, mhc.peaks$Peaks))
          # temp <- flowDensity(fs.nmzb[["MHCII"]], c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]), 
          #                     position = c(T,NA), gates=c(mhc.gate,NA))
          
          # if(temp@proportion < 5){
          #   mhc.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], use.upper = T, upper=F, tinypeak.removal = .9, alpha=.9)
          # }
          
          # if(mhc.gate < 0.75){
          #   mhc.gate <-  deGate(fs.nmzb[["MHCII"]], channels.ind.list[["MHCII"]]["MHCII"], tinypeak.removal = .9, use.upper = T, upper=T, alpha=0.01)
          # }
          
          if(mhc.gate < 2 | round(mhc.gate,2) > 2.99){
            mhc.gate <- 2
          }
          
        }else{
          mhc.gate <- deGate(f, channels.ind.list[[idx]]["MHCII"], use.percentile = T, percentile = 0.4)
          # if(mhc.gate < 0.75){
          #   mhc.gate <- deGate(f, channels.ind.list[["MHCII"]]["MHCII"], tinypeak.removal = .9,upper=T,alp)
          # }
        }
        
        if(fmo.inds["CD5"]!=-1){
          cd5.gate <- deGate(fs.nmzb[["CD5"]], channels.ind.list[[idx]]["CD5"], use.upper = T, upper=T, tinypeak.removal = .9)
          temp <- flowDensity(fs.nmzb[["CD5"]], c(channels.ind.list[[idx]]["MHCII"], channels.ind.list[[idx]]["CD5"]), 
                              position = c(NA, F), gates = c(NA, mhc.gate))
          cd5.gate <- deGate(temp, channels.ind.list[["CD5"]]["CD5"], tinypeak.removal = .9, use.upper = T, upper=T)+0.2
          
          #add a condition here to check if fmo is reliable.
          #if the fmo gate falls to the left of the first peak of the actual file, then ditch fmo gate and use deGate.
          #However, this is to assume that peak/population is always there. Which is true for all files analyzed here, but might not true for 
          #certain cases. 
          
          cd5.peaks <- getPeaks(f,c(channels.ind.list[[idx]]["CD5"]), tinypeak.removal = 0.9)
          
          if(length(grep("CD5",f@parameters@data$desc)) ==0 ){
            value <- 0.5* sd(na.omit(f@exprs[,"BV421-A"]))
          }else{
            value <- 0.5* sd(na.omit(f@exprs[,grep("CD5",f@parameters@data$desc)]))
          }
          
          if(cd5.gate <= cd5.peaks$Peaks[1] + value){
            cd5.gate <-  deGate(f, channels.ind.list[[idx]]["CD5"], use.upper = T, upper=T,tinypeak.removal = .9, alpha=.9)+0.2
          }
          if(cd5.gate < 1.5){
            ## Not using FMO
            cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"],  use.upper = T, upper=T,tinypeak.removal = .9, alpha=.95)+0.2
            if(cd5.gate > 2.75){
            #if(cd5.gate > 2.2){
              cd5.gate <- 1.5
              
            }
          }
          
        }else{
          cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"],  use.upper = T, upper=T,tinypeak.removal = .9, alpha=.95)+0.2
        }
        
        B1A <- flowDensity(f,  c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]),position = c(T,T),gates=c(mhc.gate,cd5.gate))
        ## Changed this condition. Need to check if this works for TCP and other centres.
        if(B1A@proportion < 10 & cd5.gate > 1.5){
          cd5.gate <- 1.5
          B1A <- flowDensity(f,  c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]),position = c(T,T),gates=c(mhc.gate,cd5.gate))
          
          if(B1A@proportion > 15){
            cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"],  use.upper = T, upper=T,tinypeak.removal = .9, alpha=.95)+0.2
            B1A <- flowDensity(f,  c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]),position = c(T,T),gates=c(mhc.gate,cd5.gate))
            
            
          }
        }
        plotDens(f, c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]),
                 main = paste0("Non-MZ/FO: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(B1A@filter)
        #abline(h=cd5.gate)
        
        return(B1A)
      }
      
    })
    
    # #######################################################
    # ## For now I am commenting this part because I think the gating will be MHCII vs CD5 and not CD19 vs CD5.
    # ## Gating Non-FO/MZB to obtain B1a using CD19 marker instead of MHCII
    # B1a<- lapply(1:length(fs.nmzb), function(idx){
    #   
    #   if(!(names(fs.nmzb)[[idx]]  %in% c('CD23'))){
    #     f<- fs.nmzb[[idx]]
    #     
    #     cd19.gate <- deGate(f, channels.ind.list[[idx]]["CD19"], upper=F,tinypeak.removal = .9, alpha=.9)
    #     
    #     
    #     if (fmo.inds["CD5"]!=-1){
    #       cdf.gate <- deGate(fs.nmzb[["CD5"]], channels.ind.list[[idx]]["CD11b"], upper=T,tinypeak.removal = .9)
    #       temp <- flowDensity(fs.nmzb[["CD5"]], c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
    #                           position = c(NA, F), gates = c(NA, cd19.gate))
    #       cd5.gate <- deGate(temp, channels.ind.list[["CD5"]]["CD5"], tinypeak.removal = .9, upper=T)
    #       
    #       #add a condition here to check if fmo is reliable.
    #       #if the fmo gate falls to the left of the first peak of the actual file, then ditch fmo gate and use deGate.
    #       #However, this is to assume that peak/population is always there. Which is true for all files analyzed here, but might not true for 
    #       #certain cases. 
    #       
    #       cd5.peaks <- getPeaks(f,c(channels.ind.list[[idx]]["CD5"]), tinypeak.removal = 0.9)
    #       
    #       if(length(grep("CD5",f@parameters@data$desc)) ==0 ){
    #         value <- 0.5* sd(na.omit(f@exprs[,"BV421-A"]))
    #       }else{
    #         value <- 0.5* sd(na.omit(f@exprs[,grep("CD5",f@parameters@data$desc)]))
    #       }
    #       
    #       if(cd5.gate <= cd5.peaks$Peaks[1] + value){
    #         cd5.gate <-  deGate(f, channels.ind.list[[idx]]["CD5"], upper=T,tinypeak.removal = .9, alpha=.9)
    #       }
    #       
    #     }else{
    #       cd5.gate <- deGate(f, channels.ind.list[[idx]]["CD5"], upper=T,tinypeak.removal = .9, alpha=.95)
    #     }
    #     
    #     B1A <- flowDensity(f,  c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]),position = c(T,T),gates=c(mhc.gate,cd5.gate))
    #     plotDens(f, c( channels.ind.list[[idx]]["CD19"],channels.ind.list[[idx]]["CD5"]), 
    #              main = paste0("Non-MZ/FO: ", basename(description(f)$FIL)), 
    #              cex.lab = 2, cex.axis = 2, cex.main=2)
    #     lines(B1A@filter)
    #     
    #     return(B1A)
    #   }
    #   else if( "CD43" %in% names(fmo.inds)){
    #     stop("Need to write codes based on 2C in gating strategy.")
    #   }else if(!("CD43" %in% names(fmo.inds)) & !("CD23" %in% names(fmo.inds))){
    #     stop("Need to write codes based on 2A in gating strategy.")
    #   }
    # })
    # 
    # 
    # 
    # 
    # 
    #######################################################################################
    
    
    ## Gating T and NKT cells to obtain NKT cells
    #par(mfrow = c(3,3))
    cd161.t <- lapply(1:length(fs.tnkt), function(idx){
      f <- fs.tnkt[[idx]]
      
      if (fmo.inds["CD161"]!=-1){
        cd161.gate <- deGate(fs.tnkt[["CD161"]], channels.ind.list[["CD161"]]["CD161"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.05)
        if(cd161.gate > 1.7){
          cd161.gate <- deGate(fs.tnkt[["CD161"]], channels.ind.list[["CD161"]]["CD161"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.05) - 0.25
        }
      }else{
        cd161.gate <- deGate(f, channels.ind.list[[idx]]["CD161"],  use.upper = T, upper=T,tinypeak.removal = .9)
        # if(cd161.gate > 2){
        #   cd161.gate <- deGate(f, channels.ind.list[[idx]]["CD161"],  use.upper = T, upper=T,tinypeak.removal = .9) - 0.25
        # }
          
      }
      cd161.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(cd161.gate,NA))
      
      plotDens(f, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),
               main = paste0("T and NKT: ", basename(description(f)$FIL)),
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd161.pos@filter)
      return(cd161.pos)
    })
    
    fs.nkt161p <- lapply(cd161.t, function(x) x@flow.frame)
    fs.nkt161p.proportion <- lapply(cd161.t, function(x) x@proportion)
    names(fs.nkt161p)<-names(fs)
    
    
    #######################################################################################
    ## Gating NKT cells to obtain Ly6C+ NKT cells
    
    #par(mfrow = c(3,3))
    ly6.nkt <- lapply(1:length(fs.nkt161p ), function(idx){
      
      if(!(names(fs.161p)[idx] %in% c('CD161', 'CD5', 'LY6C'))){
        f <- fs.nkt161p[[idx]]
        
        # if(fmo.inds["LY6C"]!=-1){
        #   ly6.gate <- deGate(fs.nkt161p[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9,upper=T,alpha=0.01)
        # }else{
        #   ly6.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"],upper=T,alpha=0.05)
        # }
        
        # The gate in the gating strategy is obviously at a higher point than the FMO
        ly6.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"],  alpha=0.05)
        
        ly6.pos<- flowDensity(f,c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(ly6.gate,NA))
        
        plotDens(f, c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),
                 main = paste0("T and NKT: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(ly6.pos@filter, lwd = 2)
      }else{
        ly6.pos <- NA
      }
      
      return(ly6.pos)
    })
    
    fs.ly6.nkt <- lapply(ly6.nkt, function(x) x@flow.frame)
    fs.ly6.nkt.proportion <- lapply(ly6.nkt, function(x) x@proportion)
    names(fs.ly6.nkt)<-names(fs)
    
    #######################################################################################
    ## Gating NK Cells to obtain CD11b+ NK cells
    #par(mfrow = c(3,3))
    cd11b.pos<- lapply(1:length(fs.161p), function(idx){
      
      if(names(fs.161p)[idx] != "CD161"){
        f <- fs.161p[[idx]]
        
        #We can add an if condition to ignore fmo  gates when there are 2 distinct peaks. 
        ## Need to work on the FMO part more
        if (fmo.inds["LY6C"]!=-1){
          cd11b.gate <- deGate(fs.161p[["CD11B"]], channels.ind.list[["CD11B"]]["CD11b"], tinypeak.removal = .9, use.upper = T, upper=T, alpha=0.005)
          if(cd11b.gate > 1.5){
            cd11b.gate <- deGate(fs.161p[["CD11B"]], channels.ind.list[["CD11B"]]["CD11b"], alpha = 0.0001)
            
          }
          cd11b.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),position = c(NA,T),gates=c( NA, cd11b.gate))
          cd11b.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),position = c(NA,F),gates=c( NA, cd11b.gate))
          
      }else{
          cd161.gate <- deGate(f, channels.ind.list[[idx]]["CD161"])
      
          cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"])
          if(cd11b.gate > 3){
            cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], tinypeak.removal = 0.1)
            if(cd11b.gate > 3){
              cd11b.gate <- deGate(f, channels.ind.list[[idx]]["CD11b"], use.percentile = T, percentile = 0.15)
            }
          }
          
          cd11b.pos<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),position = c(F,T),gates=c(cd161.gate, cd11b.gate))
          cd11b.neg<- flowDensity(f,c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),position = c(F,F),gates=c(cd161.gate, cd11b.gate))
          
        }
        
               
        plotDens(f, c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]),
                 main = paste0("NK: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd11b.pos@filter)
      }else{
        cd11b.pos <- NA
        cd11b.neg <- NA
      }
      
      return(list(pop1=cd11b.pos,pop2=cd11b.neg))
    })
    
    
    fs.nk11b.pos <- lapply(cd11b.pos, function(x) { 
      if(!is.na(x)){
        return(x$pop1@flow.frame)
      }else{
        return(NA)
      }
    })
    
    fs.nk11b.pos.proportion <- lapply(cd11b.pos, function(x) { 
      if(!is.na(x)){
        return(x$pop1@proportion)
      }else{
        return(NA)
      }
    })
    
    fs.nk11b.neg <- lapply(cd11b.pos, function(x) {
      if(!is.na(x)){
        return( x$pop2@flow.frame)
      }else{
        return(NA)
      }
    })
    
    fs.nk11b.neg.proportion <- lapply(cd11b.pos, function(x) {
      if(!is.na(x)){
        return( x$pop2@proportion)
      }else{
        return(NA)
      }
    })
    
    names(fs.nk11b.pos)<-names(fs.nk11b.neg)<-names(fs)
    
    #######################################################################
    
    #Gating CD11b+NK cells to obtain Ly6C+CD11b+NK cells
    ly6c.nk11.pos <- lapply(1:length(fs.nk11b.pos), function(idx){
      
      if(!(names(fs.161p)[idx] %in% c('CD11B', 'CD161'))){#if this is not a fmo file
        f <- fs.nk11b.pos[[idx]]
        
        if (fmo.inds["LY6C"]!=-1){
          # ly6.gate <- deGate(fs.nk11b.pos[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9,upper=T,alpha=0.01)
          # ly6.gate <- deGate(fs.nk11b.pos[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .99,upper=T,alpha=0.05)
          # ly6.gate <- deGate(fs.nk11b.pos[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9,upper=T,alpha=0.01, 
          #                    twin.factor = 0.5)
          ly6c.gate <-  deGate(fs.nk11b.pos[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9,use.upper = T, upper=T,alpha=0.001)
          if(ly6c.gate < 1.65){
            ## Not using FMO for this condition
            ly6c.gate.peaks <- getPeaks(f, channel =  channels.ind.list[[idx]]["Ly6C"])
            if(length(ly6c.gate.peaks$Peaks) > 1){
              ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"])
              if(ly6c.gate > 2){
                ly6c.gate <- 1.5
              }
            }else{
              ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"], use.upper = T, upper = F)
            }
              
            
            
          }else if(round(ly6c.gate,2) >= 1.70){
            ly6c.gate <- deGate(fs.nk11b.pos[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], use.percentile = T, percentile = 0.998)
            if(ly6c.gate > 1.75){
              ## Not using FMO for this condition
              ly6c.gate.peaks <- getPeaks(f, channel =  channels.ind.list[[idx]]["Ly6C"])
              if(length(ly6c.gate.peaks$Peaks) > 1){
                ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"])
              }else{
                ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"], use.upper = T, upper = F)
              }
            }
          }
          
        }else{
          ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"])
          if(ly6c.gate > 3){
            ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"], tinypeak.removal = 0.9)
          } 
        }
        
        temp <-flowDensity(f,c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(ly6c.gate,NA))
        
        if(temp@proportion < 10){
          ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"], use.upper = T, upper=T,alpha=0.05)
          if(ly6c.gate > 3){
            ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"],alpha=0.05)
            
          }
        }
        
        rm(temp)
        
        ly6c.pos<- flowDensity(f,c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(ly6c.gate,NA))
        
        plotDens(f, c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),
                 main = paste0("NK CD11b+: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(ly6c.pos@filter)
      }else{
        ly6c.pos <- NA
      }
      
      return(ly6c.pos)
    })
    
   
    
    fs.ly6c.nk11.pos <- lapply(ly6c.nk11.pos, function(x) x@flow.frame)
    fs.ly6c.nk11.pos.proportion <- lapply(ly6c.nk11.pos, function(x) x@proportion)
    names(fs.ly6c.nk11.pos)<-names(fs)
    
    #######################################################################
    
    #######################################################################
    
    #Gating CD11b-NK cells to obtain Ly6C+CD11b-NK cells
    ## Excluding the FMOs from this step since this is the last step of gating strategy
    ly6c.nk11.neg <- lapply(1:nrow(x), function(idx){
      f <- fs.nk11b.neg[[idx]]
      
      if(!(names(fs.161p)[idx] %in% c('CD161'))){
        
        if (fmo.inds["LY6C"]!=-1){
          # ly6.gate <- deGate(fs.nk11b.neg[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9,
          #                    upper=T, alpha=0.005, twin.factor = 0.2)
          
          ly6c.gate <- deGate(fs.nk11b.neg[["LY6C"]], channels.ind.list[["LY6C"]]["Ly6C"], tinypeak.removal = .9, use.upper = T, upper=T,alpha=0.0001)
        }else{
          ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"])-0.5
        }
        
        ly6c.pos<- flowDensity(f,c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(ly6c.gate,NA))
        
        if(ly6c.pos@proportion > 50){
          ly6c.gate <- deGate(f, channels.ind.list[[idx]]["Ly6C"])
          ly6c.pos<- flowDensity(f,c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),position = c(T,NA),gates=c(ly6c.gate,NA))
          
        }
        if(nrow(na.omit(f@exprs)) <= 2){
          print("Zero neg population, skipping plotting...")
        }else{
          plotDens(f, c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]),
                   main = paste0("NK CD11b-: ", basename(description(f)$FIL)),
                   cex.lab = 2, cex.axis = 2, cex.main=2)
          lines(ly6c.pos@filter)
        }
        
        
      }else{
        ly6c.pos <- NA
      }
      return(ly6c.pos)
    })
    
    
    fs.ly6c.nk11.neg <- lapply(ly6c.nk11.neg, function(x) x@flow.frame)
    fs.ly6c.nk11.neg.proportion <- lapply(ly6c.nk11.neg, function(x) x@proportion)
    names(fs.ly6c.nk11.neg)<-names(fs)
    
    graphics.off()
    # -------------------------------------------------------------------------------------
    # Plots ------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------
    
    for(idx in 1:length(fs.preflowCut)){
      CairoPNG(file = paste0(scatterplot.dir, i1, '_', gsub('.fcs', '.png', names(fs)[[idx]])), width = 2000*4/4, height = 2000*6/4)
      par(mfrow = c(6,4), mar = (c(5, 5, 4, 2) + 0.1))
      
      f <- fs[[idx]]
      
      plotDens(fs.preflowCut[[idx]], c(channels.ind.list[[idx]]['Time'], channels.ind.list[[idx]]['CD19']), main = paste0("Post scatter margin removal: ", basename(description(f)$FIL)), cex.lab = 2, cex.axis = 2, cex.main = 2)
      points(fs.preflowCut[[idx]]@exprs[fcut[[idx]]$ind, c(channels.ind.list[[idx]]['Time'], channels.ind.list[[idx]]['CD19'])], pch = '.')
      
      plotDens(fs[[idx]], c(scat.chans["FSC-A"], channels.ind.list[[idx]]["Live"]), main = paste0("CLOG- ", basename(description(f)$FIL)), cex.lab = 2, cex.axis = 2, cex.main = 2)
      lines(live.filter[[idx]], type="l",lwd=2)
      
      plotDens(fs.live[[idx]], channels = c("FSC-A","SSC-A"), main = paste0("Live: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      # lines(lymph.filter[[idx]])
      
      plotDens(fs.live[[idx]], channels = c("FSC-A","FSC-H"), main = paste0("FCS SINGLETS: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(fcssinglets.filter[[idx]], type="l",lwd=2)
      
      
      plotDens(fs.fcssinglets[[idx]], channels = c("SSC-W","SSC-H"), main = paste0("SSC SINGLETS: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(singlet.filter[[idx]], type="l",lwd=2)
      

      if(length(F4marker.index) == 1){
        plotDens(f, c(channels.ind.list[[idx]]["F4/80"], channels.ind.list[[idx]]["CD11b"]), main = paste0("F4/80 Neg & Mac: ", basename(description(f)$FIL)),
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(macrophages[[idx]]@filter, type="l",lwd=2)
        
        plotDens(fs.f4.80neg[[idx]], c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), main = paste0("Singlets: ", basename(description(f)$FIL)), 
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd19pos[[idx]]@filter, type="l",lwd=2)
        lines(cd19neg[[idx]]@filter, type="l",lwd=2)
      }else{
        plotDens(fs.singlets2[[idx]], c(channels.ind.list[[idx]]["CD19"], channels.ind.list[[idx]]["CD5"]), main = paste0("Singlets: ", basename(description(f)$FIL)), 
                 cex.lab = 2, cex.axis = 2, cex.main=2)
        lines(cd19pos[[idx]]@filter, type="l",lwd=2)
        lines(cd19neg[[idx]]@filter, type="l",lwd=2)
        
      }
      
    
      plotDens(fs.19n[[idx]], c(channels.ind.list[[idx]]["CD5"], channels.ind.list[[idx]]["CD11b"]), 
               main = paste0("CD19- : ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(tnkt.filter[[idx]], type="l",lwd=2) 
      
      plotDens(fs.not.tnkt[[idx]], c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["MHCII"]), 
               main = paste0("NOT(T and NKT): ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd161[[idx]]$pop1@filter, type="l",lwd=2)
      
      plotDens(fs.19p[[idx]], c( channels.ind.list[[idx]]["CD19"],channels.ind.list[[idx]]["CD161"]), 
               main = paste0("CD19+: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd19.161[[idx]]$pop1@filter, type="l",lwd=2)
      
      plotDens(fs.19.n161[[idx]], c( channels.ind.list[[idx]]["CD23"],channels.ind.list[[idx]]["CD21/CD35"]), 
               main = paste0("CD19+CD161-: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(MZB[[idx]]$mzb@filter, type="l",lwd=2)
      lines(MZB[[idx]]$follicular@filter, type="l",lwd=2)
      lines(MZB[[idx]]$nonfomzb@filter, type="l",lwd=2)
      
      plotDens(fs.nmzb[[idx]], c( channels.ind.list[[idx]]["MHCII"],channels.ind.list[[idx]]["CD5"]), 
               main = paste0("Non-MZ/FO: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(B1a[[idx]]@filter, type="l",lwd=2)
      
      plotDens(fs.tnkt[[idx]], c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]), 
               main = paste0("T and NKT: ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd161.t[[idx]]@filter, type="l",lwd=2)
      
      plotDens2(fs.nkt161p[[idx]], c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]), 
                main = paste0("T and NKT/CD161+: ", basename(description(f)$FIL)), 
                cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(ly6.nkt[[idx]]@filter, type="l",lwd=2)
      
      plotDens2(fs.161p[[idx]], c(channels.ind.list[[idx]]["CD161"], channels.ind.list[[idx]]["CD11b"]), 
                main = paste0("NK: ", basename(description(f)$FIL)), 
                cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cd11b.pos[[idx]]$pop1@filter, type="l",lwd=2)
      
      plotDens2(fs.nk11b.pos[[idx]], c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]), 
                main = paste0("NK CD11b+: ", basename(description(f)$FIL)), 
                cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(ly6c.nk11.pos[[idx]]@filter, type="l",lwd=2)
      
      plotDens2(fs.nk11b.neg[[idx]],  c(channels.ind.list[[idx]]["Ly6C"], channels.ind.list[[idx]]["CD11b"]), 
                main = paste0("NK CD11b-: ", basename(description(f)$FIL)), 
                cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(ly6c.nk11.neg[[idx]]@filter, type="l",lwd=2)
      
      plotDens(fs.161n[[idx]], c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD5']), 
               main = paste0("NOT(NK): ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(granulo.filter[[idx]], type="l",lwd=2)
      
      plotDens(fs.notgranulo[[idx]], c(channels.ind.list[[idx]]['CD11b'], which(colnames(f) == 'SSC-H')), 
               main = paste0("NOT(granulocytes): ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(eosino[[idx]]$eosino.fD@filter, type="l",lwd=2)
      
      plotDens(fs.noteosino[[idx]], c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['Ly6C']),
               main = paste0("NOT(eosinophils): ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(monocytes[[idx]]$monocytes.flowD@filter, type="l",lwd=2)
      
      plotDens(fs.notmono[[idx]], c(channels.ind.list[[idx]]['CD11c'], channels.ind.list[[idx]]['MHCII']),
               main = paste0("NOT(monocytes): ", basename(description(f)$FIL)), 
               cex.lab = 2, cex.axis = 2, cex.main=2)
      lines(cDCs[[idx]]@filter,lwd=2)
      
      plotDens2(fs.cDCs[[idx]], c(channels.ind.list[[idx]]['CD11b'], channels.ind.list[[idx]]['CD11c']),
                main = paste0("cDCs: ", basename(description(f)$FIL)), 
                cex.lab = 2, cex.axis = 2, cex.main=2); abline(v = cd11b.cDCs.gate[[idx]],lwd=2)
      lines(cd11b.cDCs[[idx]]$cd11b.cDCs.pos.flowD@filter, lwd = 2)
      lines(cd11b.cDCs[[idx]]$cd11b.cDCs.neg.flowD@filter, lwd = 2)
      
      
      dev.off()
    }
    
    # # -------------------------------------------------------------------------------------
    # # Cell counts ------------------------------------------------------------------------------
    # # -------------------------------------------------------------------------------------
    # 
    # if(length(F4marker.index) == 1){
    #   
    # }else{
    #   
    # }
    # names.cell.counts <- c('FCS.files','Passed flowCut',
    #                        'All events', 'Live','FSC Singlets', 'SSC Singlets', 'B cells', 
    #                        'CD19-', 'T and NKT', 'NK cells', 'CD19 CD161-', 'MZ B', 'FO B', 'Non-FO/MZ B',
    #                        'B1a', 'NKT cells', 'Ly6C+ NKT', 'CD11b+ NK', 'CD11b- NK', 'Ly6C+ CD11b+ NK',
    #                        'Ly6C+ CD11b- NK', 'Granulocytes', 'Eosinophils', 'Monocytes', 'cDCs',
    #                        'CD11b+ cDC', 'CD11b- cDC')
    # 
    # cell.counts <- ldply(1:length(fs.preflowCut), function(idx){
    #   cc.vector <- c(names(fs)[[idx]], passed.flowCut[[idx]],
    #                  FCS.cellcount[idx], 
    #                  nrow(na.omit(fs.live[[idx]]@exprs)),
    #                  nrow(na.omit(fs.fcssinglets[[idx]]@exprs)),
    #                  nrow(na.omit(fs.singlets2[[idx]]@exprs)), nrow(na.omit(fs.19p[[idx]]@exprs)),
    #                  nrow(na.omit(fs.19n[[idx]]@exprs)), nrow(na.omit(fs.tnkt[[idx]]@exprs)), 
    #                  nrow(na.omit(fs.161p[[idx]]@exprs)), nrow(na.omit(fs.19.n161[[idx]]@exprs)),
    #                  MZB[[idx]]$mzb@cell.count, MZB[[idx]]$follicular@cell.count, MZB[[idx]]$nonfomzb@cell.count,
    #                  B1a[[idx]]@cell.count, nrow(na.omit(fs.nkt161p[[idx]]@exprs)), ly6.nkt[[idx]]@cell.count,
    #                  nrow(na.omit(fs.nk11b.pos[[idx]]@exprs)), nrow(na.omit(fs.nk11b.neg[[idx]]@exprs)), 
    #                  ly6c.nk11.pos[[idx]]@cell.count, ly6c.nk11.neg[[idx]]@cell.count, granulo.cellcount[[idx]],
    #                  fs.eosino.cellcount[[idx]], fs.monocytes.cellcount[[idx]], nrow(na.omit(fs.cDCs[[idx]]@exprs)),
    #                  fs.cd11b.cDCs.pos.cellcount[[idx]], fs.cd11b.cDCs.neg.cellcount[[idx]])
    #   
    #   names(cc.vector) <- names.cell.counts
    #   return(cc.vector)
    # })
    # 
    # cell.counts <- merge(x, cell.counts, by = c('FCS.files'))
    # 
    # cell.proportions <- ldply(1:length(fs.preflowCut), function(idx){
    #   cc.vector <- c(names(fs)[[idx]], passed.flowCut[[idx]],
    #                  FCS.cellcount[idx], fs.live.proportion[[idx]], fs.fcssinglets.proportion[[idx]],
    #                  fs.singlets.proportion[[idx]], fs.19p.proportion[[idx]], fs.19n.proportion[[idx]],
    #                  fs.tnkt.proportion[[idx]], fs.161p.proportion[[idx]], fs.19.n161.proportion[[idx]],
    #                  MZB[[idx]]$mzb@proportion, MZB[[idx]]$follicular@proportion, MZB[[idx]]$nonfomzb@proportion,
    #                  B1a[[idx]]@proportion, fs.nkt161p.proportion[[idx]], ly6.nkt[[idx]]@proportion, fs.nk11b.pos.proportion[[idx]],
    #                  fs.nk11b.neg.proportion[[idx]], ly6c.nk11.pos[[idx]]@proportion, ly6c.nk11.neg[[idx]]@proportion, 
    #                  granulo.proportion[[idx]], fs.eosino.proportion[[idx]], fs.monocytes.proportions[[idx]], fs.cDCs.proportions[[idx]],
    #                  fs.cd11b.cDCs.pos.proportions[[idx]], fs.cd11b.cDCs.neg.proportions[[idx]])
    #   
    #   
    #   
    #   names(cc.vector) <- names.cell.counts
    #   return(cc.vector)
    # })
    # 
    # cell.proportions <- merge(x, cell.proportions, by = c('FCS.files'))
    # 
    # 
    # #save(channels.ind.list, file = paste0(results.dir, "/test.RData"))
    # save(cell.counts, file = paste0(results.dir, "/Cell_Counts/", i1, "_cellcounts.RData"))
    # write.csv(cell.counts, file = paste0(results.dir, "/Cell_Counts/", i1, "_cellcounts.csv"), row.names = F)
    # 
    # save(cell.proportions, file = paste0(results.dir, "/Cell_Proportions/", i1, "_cellproportions.RData"))
    # write.csv(cell.proportions, file = paste0(results.dir, "/Cell_Proportions/", i1, "_cellproportions.csv"), row.names = F)
    # 
  }) # end of tryCatch
  
  return(data.frame(x$FCS.files))
  
}, .parallel = TRUE) # end ldply

###########################################################################

## Saving the event counts in one spreadsheet for the centre
store.allCSV.temp <- NULL
pathCSV <- paste0(results.dir,"/Cell_Counts/") 
# Reads all folders and files in current path folder and makes a list of all of their paths
allCSV <- dir(pathCSV, full.names=T, recursive=T, pattern = "*.csv") 
store.allCSV.temp <- sapply(1:length(allCSV), function(x){pathCSV})
#store.allCSV.temp <- cbind(store.allCSV.temp, sapply(1:length(allCSV), function(x){unlist(strsplit(allCSV[x], split = "/"))[length(unlist(strsplit(allCSV[x], split = "/")))-1]}))
store.allCSV.temp <- cbind(store.allCSV.temp, sapply(1:length(allCSV), function(x){unlist(strsplit(allCSV[x], split = "/"))[length(unlist(strsplit(allCSV[1], split = "/")))]}))

duplicate.index <- which(duplicated(store.allCSV.temp[,2])==TRUE)
all.events.store <- NULL
if(length(duplicate.index) == 0){
  for(i in 1:nrow(store.allCSV.temp)){
    all.events.store.temp <- read.csv(paste0(results.dir,"/Cell_Counts/", store.allCSV.temp[i,2]) )
    all.events.store <- rbind(all.events.store, all.events.store.temp)
  }
}

## Saving the proportions in one spreadsheet for the centre
store.allCSV.temp <- NULL
pathCSV <- paste0(results.dir,"/Cell_Proportions/") 
# Reads all folders and files in current path folder and makes a list of all of their paths
allCSV <- dir(pathCSV, full.names=T, recursive=T, pattern = "*.csv") 
store.allCSV.temp <- sapply(1:length(allCSV), function(x){pathCSV})
#store.allCSV.temp <- cbind(store.allCSV.temp, sapply(1:length(allCSV), function(x){unlist(strsplit(allCSV[x], split = "/"))[length(unlist(strsplit(allCSV[x], split = "/")))-1]}))
store.allCSV.temp <- cbind(store.allCSV.temp, sapply(1:length(allCSV), function(x){unlist(strsplit(allCSV[x], split = "/"))[length(unlist(strsplit(allCSV[1], split = "/")))]}))

duplicate.index <- which(duplicated(store.allCSV.temp[,2])==TRUE)
all.props.store <- NULL
if(length(duplicate.index) == 0){
  for(i in 1:nrow(store.allCSV.temp)){
    all.props.store.temp <- read.csv(paste0(results.dir,"/Cell_Proportions/", store.allCSV.temp[i,2]) )
    all.props.store <- rbind(all.props.store, all.props.store.temp)
  }
}


names.columns <- c('FCS files', 'Path', 'Panel/Organ/Folder','Genotype','Strain Code/Ear Tag','Mouse ID/Animal ID',
                   'Colony ID/Line','Assay Date','Gender','Number of Channels','Number of Cells','Passed flowCut',
                   'All Events', 'Live','FSC Singlets', 'SSC Singlets', 'B cells', 
                   'CD19-', 'T and NKT', 'NK cells', 'CD19 CD161-', 'MZ B', 'FO B', 'Non-FO/MZ B',
                   'B1a', 'NKT cells', 'Ly6C+ NKT', 'CD11b+ NK', 'CD11b- NK', 'Ly6C+ CD11b+ NK',
                   'Ly6C+ CD11b- NK', 'Granulocytes', 'Eosinophils', 'Monocytes', 'cDCs',
                   'CD11b+ cDC', 'CD11b- cDC')

colnames(all.events.store) <- names.columns
colnames(all.props.store) <- names.columns

dockerImage <-c("0b1d76dacd88")
all.events.store <- cbind(all.events.store, dockerImage)
all.props.store <- cbind(all.props.store, dockerImage)

save(all.props.store, file = paste0(results.dir,"/all.props.store.Rdata"))
save(all.events.store, file = paste0(results.dir,"/all.events.store.Rdata"))

## Saving the spreadsheets for DCC and Centres
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.events.store, file =  paste0(results.dir, "/DCCResults_EventCounts_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.props.store, file =  paste0(results.dir, "/DCCResults_Proportions_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

## Saving the spreadsheets for us
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.events.store, file =  paste0(results.dir, "/EventCounts_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(all.props.store, file =  paste0(results.dir, "/Proportions_",toupper(centre),"_Panel2", date.time), row.names = FALSE)



cat("One acquistion date completed in: ", TimePrint(startloop), "\n", sep="")
print("Finished Gating & Plotting")
