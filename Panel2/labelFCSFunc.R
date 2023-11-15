## Created by Albina Rahim
## Date: August 07, 2019
## This function will relabell the files based on the naming convention that was decided upon by IMPC.
## In addition, if the file is in .LMD format it will convert it into FCS files.
## If the description field of the files are missing (i.e. NA), it will map the channels to the correct marker names.
## The marker names will have to be provided by the centre. 


labelFCSFunc <- function(inputPath, outputPath, markerNames, panelSuffix){
  library("flowCore")
  library("stringr")
  
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
    write.FCS(f, file.path(outputPath, paste0(panelSuffix, strsplit(identifier(f), split = " ")[[1]][1],".fcs")))

  }
}