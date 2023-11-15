## Developed by Albina Rahim
## Date: December 13, 2016
## This function creates GlobalFrame by storing 1000 random cells from each FCS file. 
## As input it requires the store.allFCS matrix and the number of CPUs to use while parallelization
## The store.allFCS matrix is created as an output from the preProcessingFunc.R function, which contains
## information on the Paths, Name of the Files, Genotypes, Barcodes, Assay Dates, Gender, Number of Channels, Number of Cells for each FCS file


globalFrameFunc <- function(store.allFCS, centre, panel, outputPath){
    library("flowCore")
    library("flowBin")
    
    centre <- tolower(centre)
    
    print("Start creating the Global Frame")
    rownames(store.allFCS) <- 1:nrow(store.allFCS) 
    file.names <- data.frame(store.allFCS, stringsAsFactors = F)
    #file.names <- data.frame(store.allFCS[1:2500,], stringsAsFactors = F)
    
  
    if(centre == "ucd" & panel == 2){
        f <- try(read.FCS(filename = paste0(store.allFCS[1,c('Path')], "/", store.allFCS[1,c('Panel/Organ/Folder')], "/", store.allFCS[1,c('FCS files')])), silent = TRUE)
         markers.temp <- f@parameters@data$desc
    }
    # index <- 0
    # for(i in 1:nrow(store.allFCS)){
    #   f <- read.FCS(paste0(store.allFCS[i,1], "/", store.allFCS[i,2], "/", store.allFCS[i,4]))
    #   if(is.na(f@parameters@data$desc[7])){
    #     index <- c(index,i)
    #   }
    # }
    # index <-index[index !=0]
    # 
    # store.allFCS.NA <- store.allFCS[index,]
    # write.csv(store.allFCS.NA, file =  paste0(outputPath, "store.allFCS.NA.csv"), row.names = FALSE)
    
    gFrame <- ddply(file.names, "FCS.files", function(x){
      x<- file.names[i,]
     
      if(centre == "sanger" | centre == "ciphe"){
        f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
        fpath <- x$Path
      }else if(centre == "tcp" | centre == "bcm" | centre == "jax" | centre == "gmc" | centre == "ucd"){
        f <- try(read.FCS(filename = paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)), silent = TRUE)
        fpath <- paste0(x$Path,"/", x$Panel.Organ.Folder)
      }
      
      
      f <- compensateIMPC(f, basename(x$FCS.files), fpath, centre = centre, panel.no = panel)
      
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
   
      channels.ind <- sort(channels.ind)
      
      temp <- f@exprs[sample(1:length(f@exprs[,1]), 1000), as.numeric(channels.ind)]
     
      return(data.frame(temp, check.names = F))
      
    }, .parallel = TRUE) # end ddply
    
    # ## Arranging the output from ddply() in order of the file.names
    # gFrame <- join(file.names, gFrame)
    library(tidyr)
    
   
    if(centre == "tcp" & panel == 2){
     temp <- gFrame
     temp[is.na(temp)] = ''
     temp1 <- unite(temp, "FITC-A", c('B525-A', 'FITC-A'), sep='')
     temp2 <- unite(temp1, "APC-A", c('R670-A','R660-A', 'APC-A'), sep='')
     temp3 <- unite(temp2, "APC-Cy7-A", c('R780-A', 'APC-Cy7-A'), sep='')
     temp4 <- unite(temp3, "BV421-A", c('V450-A', 'V431-A', 'BV421-A'), sep='')
     temp5 <- unite(temp4, "BV510-A", c('V525-A', 'BV510-A'), sep='') ## Sytox Blue/CD19
    
     temp5 <- temp5 %>% 
       mutate("BV510-A2" = temp5$`BV510-A`)
     temp6 <- unite(temp5, "BV510-A", c('V480-A', 'BV510-A'), sep='') ## SYTOX BLUE/CD19
     temp7 <- unite(temp6, "V660-A", c('V660-A', 'BV510-A2'), sep='') ## CD19
     temp8 <- unite(temp7, "PE-A", c('Y586-A', 'PE-A'), sep='')
     temp9 <- unite(temp8, "PE-CF594-A", c('Y610-A', 'PE-CF594-A'), sep='')
     temp10 <- unite(temp9, "BV786-A", c('V780-A', 'BV786-A'), sep='')
     temp11 <- unite(temp10, "PE-Cy7-A", c('Y780-A', 'PE-Cy7-A'), sep='')

     temp12 <- matrix(nrow = nrow(temp11), ncol = ncol(temp11)-1, data = NA)
     
     for(i in 1:ncol(temp12)){
       temp12[,i] <- temp11[,i+1]
     }
     mode(temp12) = "numeric"
     colnames(temp12) <- c("B525-A", "R660-A", "R780-A", 
                           "V431-A", "V480-A",  "V660-A", 
                           "V780-A", "Y586-A", "Y610-A", "Y780-A")
     gFrame <- temp12
     
     temp13 <-temp12
     temp12 <- cbind(temp11[,1], temp12)
     
     temp12 <- temp11[, c("FCS.files", "FITC-A", "APC-A", "APC-Cy7-A", 
                          "BV421-A", "BV510-A",  "V660-A", 
                          "BV786-A", "PE-A", "PE-CF594-A", "PE-Cy7-A")]
     # colnames(temp12) <- c("FCS.files", "FITC-A/B525-A", "APC-A/R670-A/R660-A", "APC-Cy7-A/R780-A", 
     #                       "BV421-A/V431-A/V450-A", "BV510-A/V525-A/V480-A",  "BV510-A/V525-A/V660-A", 
     #                       "BV786-A/V780-A", "PE-A/Y586-A", "PE-CF594-A/Y610-A", "PE-Cy7-A/Y780-A")
     
    
    }
    ## Remove rows in gFrame with NAs and also saving files which has NA in their flurophore/channels
    fileswNAs <- NULL
    colswNAs <- which(is.na(gFrame), arr.ind = TRUE)
    if(length(colswNAs) > 0){
      rowswNAs <- unique(colswNAs[,'row'])
      colswNAs <- unique(colswNAs[,'col'])
      fileswNAs <- unique(gFrame[rowswNAs,1])
      gFrame <- gFrame[-rowswNAs,]
    }
    
    ## Reading the first FCS file in the storage matrix as a template for creating the global frame
    if(centre == "sanger" | centre == "ciphe" ){
      g <- read.FCS(filename = paste0(store.allFCS[1,c('Path')], "/", store.allFCS[1,c('FCS files')]))
    }else if(centre == "tcp" | centre == "gmc" | centre == "ucd"){ ## TCP, GMC, and UCD have inconsistency in channel number. So we are taking the one which has the maximum channel number
     channelNum <- which(store.allFCS[,c('Number of Channels')] == max(store.allFCS[,c('Number of Channels')]))[1]
     g <- read.FCS(filename = paste0(store.allFCS[channelNum,c('Path')], "/", store.allFCS[channelNum,c('Panel/Organ/Folder')],"/", store.allFCS[channelNum,c('FCS files')]))
    }else if(centre == "bcm" | centre == "jax"){
      g <- read.FCS(filename = paste0(store.allFCS[1,c('Path')], "/", store.allFCS[1,c('Panel/Organ/Folder')],"/", store.allFCS[1,c('FCS files')]))
    }
    
    Transform.idx <- unlist(sapply(colnames(gFrame), function(x) {grep(x, colnames(g))})) 
    gexprs.temp <- matrix(0, nrow = nrow(gFrame), ncol = ncol(g@exprs))
    if(centre == "tcp"){
      gexprs.temp[, Transform.idx] <- as.matrix(gFrame[, 1:ncol(gFrame)])
    }else {
      gexprs.temp[, Transform.idx] <- as.matrix(gFrame[, 2:ncol(gFrame)])
    }
   
    g@exprs <- gexprs.temp
    colnames(g@exprs) <- colnames(g) 
    
    print("End of creating the Global Frame")
    
    ## Since I am already removing rows with NAs from gFrame the following section is not needed
    # # The following is needed since for centres TCP and for CIPHE (CD23 has 2 possible channels/fluorophores)
    # if(centre == "tcp"){
    #   colswNAs <- which(is.na(g@exprs[,Transform.idx]), arr.ind = TRUE)
    #   fileswNAs <- gFrame[colswNAs[,1],]
    #   fileswNAs <- unique(fileswNAs[,1])
    #   colswNAs <- unique(colswNAs[,'col'])
    #   
    # }else{
    #   colswNAs <- which(is.na(g@exprs), arr.ind = TRUE)
    #   fileswNAs <- gFrame[colswNAs[,1],]
    #   fileswNAs <- unique(fileswNAs[,1])
    #   colswNAs <- unique(colswNAs[,'col'])
    # }
    
    ## Tranformation
    print("Start computing the transform using the estimateLogicle()")
    lgl <- estimateLogicle(g, channels = colnames(g)[Transform.idx])
    print("End of computing the transform using the estimateLogicle()")
    
    return(list(globalFrame = g, globalFrame.Matrix = gFrame, lgl = lgl, files.w.NAs = fileswNAs))
}



