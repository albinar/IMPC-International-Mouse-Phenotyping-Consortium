## Created by Albina Rahim
## Date: August 07, 2019
## This function will relabell the files based on the naming convention that was decided upon by IMPC.
## In addition, if the file is in .LMD format it will convert it into FCS files.
## If the description field of the files are missing (i.e. NA), it will map the channels to the correct marker names.
## The marker names will have to be provided by the centre. 


labelFCSFunc <- function(inputPath, outputPath, markerNames, panelPrefix){
 
  ## This tool requires the flowCore library
  if(!require("flowCore")) {
    cat("flowCore library not present, trying to install flowCore\n")
    source("http://bioconductor.org/biocLite.R")
    biocLite("flowCore")
    if(require("flowCore")) {
      cat("flowCore library installed successfully\n")
    } else {
      stop("Failed to install the flowCore library, please try installing flowCore manually and try again.", call. = FALSE)
    }
  }

  ## This tool requires the stringr package
  install.packages("stringr")

  #################################################################################
  ## Main code starts----

  fileExt <- c("LMD", "lmd")
  
  allLMD <- dir(inputPath, full.names=T, recursive=T)
  
  tempCol <- sapply(1:length(allLMD), function(x){unlist(strsplit(allLMD[x], split = "/"))[length(unlist(strsplit(allLMD[x], split = "/")))]})
  
  for (i in 1:length(tempCol)){
    matches <- grep(paste(fileExt,collapse="|"), tempCol[i], value=FALSE)
    if(length(matches) != 0){
      f <- read.FCS(allLMD[i], transformation = FALSE, dataset = 2)
    }else{
      f <- read.FCS(allLMD[i])
    }
    
    if(is.na(f@parameters@data$desc[1])){
      for(j in 1:ncol(f)){
        f@parameters@data$desc[j] <- markerNames[j]
      }
    }
    
    if(panelPrefix == "FMO_"){
      write.FCS(f, file.path(outputPath, paste0(strsplit(identifier(f), split = ".LMD")[[1]],".fcs")))
    }else{
      write.FCS(f, file.path(outputPath, paste0(panelPrefix, strsplit(identifier(f), split = " ")[[1]][1],".fcs")))
      
    }
 
  }
}