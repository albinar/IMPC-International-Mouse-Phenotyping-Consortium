## Written by Albina Rahim
## Last Updated: Feb 22, 2016

## Note: source("Codes/Panel1-CoreGating.R")
## listKnownIlarCodes(): "BCM"  "Gmc"  "H"    "Ics"  "J"    "Ning" "Rbrc" "Rbri" "Tcp"  "Ucd"  "WTSI"

## Accesses the FlowRepositoryR database, retrieves the datasets for different centres. Checks for duplicates and corrupted flowFrames and gets rid of them. 
## Also checks if the flowFrame has already been downloaded or not. If yes, then skips it, else downloads it.


## Setting the working directory-----
setwd("/home/CRC/arahim/Documents/IMPC-Pipeline/")

source('Codes/changeChannels.R')
source('Codes/swapChannels.R')
source('Codes/helperfunc.R')

# Loading the FlowRepositoryR & flowCore packages-----
library(FlowRepositoryR)
library(flowCore)
library(flowClean)

# Make sure you have the latest version of FlowRepositoryR, which at this point is 1.3.38-----
sessionInfo()

# At this point we can test the data against a local FlowRepository that is running in Josef's computer.
# Therefore, we set the FlowRepository URL as follows -----
setFlowRepositoryURL("http://bioinfo17l.bccrc.ca")

ptm <- proc.time()

# Matrix containing information of all the corrupted FCS files
allCorruptedFCS <- matrix(ncol = 3, byrow = TRUE)
colnames(allCorruptedFCS) <- c("Center Code", "ID", "FCS File Name")

# Variables for tracking how many Knockout FCS files for Panel1 we are working with
countFCS.KO.P1 <- 0

# Center Name
centre <- "Tcp"

# Panel
Panel <- "Panel1"

# Creating a Data directory and then separate for each center
data.path <- paste("Data",paste0("Center-", centre),sep = "/")
suppressWarnings ( dir.create ( paste("Data",paste0("Center-", centre),sep = "/") ,recursive = T))

## Two matrices for Wild Types:
## 1. allWildTypes contains information for all the Wild types of a center
## 2. oneWildType contains information for one Wild type at a time
allWildTypes <- matrix(ncol = 6, byrow = TRUE)
colnames(allWildTypes) <- c("Center Code", "ID", "Name", "FCS File Index", "FCS File Name", "md5sum")


## Downloading the Wild Types-----
# This function connects to FlowRepository (flowrepository.org) via an XML-based API and retrieves a vector of identifiers of available datasets.
# It is assumed that by baseline it refers to WildType (not completely sure at the moment)-----
wildTypeList <- flowRep.ls(impc.only = TRUE, impc.centre = centre, impc.specimen.baseline = TRUE, impc.unanalyzed.only = TRUE)

# The flowRep.get() function connects to FlowRepository (flowrepository.org) via an XML-based API and
# retrieves metadata about a specified dataset in the form of a flowRepData object
# We are retrieving the Wild Type dataset-----
wildTypeData <- NULL
countCorruptedWT <- 0
for ( k1 in 1:length(wildTypeList) ){
  wildTypeData <- flowRep.get(wildTypeList[k1], use.credentials=TRUE,impc.details=TRUE)
  summary(wildTypeData) # Gives a summary of the dataset, including the number of FCS files and Experiment name

  ## We want to download one Wild type Line at a time
  ## NOTE: sapply() adds a vector to generate the matrix oneWildType. The vectors are added by column. Hence, we had to use the transpose(t()) function after sapply to obtain the final oneWildType matrix.
  oneWildType <- t(sapply(1:length(wildTypeData@fcs.files), function(l1){c(centre, wildTypeData@id, unlist(strsplit(wildTypeData@name," "))[3], l1, wildTypeData@fcs.files[[l1]]@name, wildTypeData@fcs.files[[l1]]@md5sum)}))
  colnames(oneWildType) <- c("Center Code", "ID", "Name", "FCS File Index", "FCS File Name", "md5sum")


  ## Removing duplicate FCS files, since, FlowRepository seemed to have duplicate FCS files.
  ## The duplicate files have been removed based on their "md5sum" values.
  oneWildType <- oneWildType[!duplicated(oneWildType[,6]),]
  # Storing all the Wild types in the allWildTypes matrix
  allWildTypes <- rbind(allWildTypes, oneWildType)

  ## Checks if the FCS files already exists. If not, download the FCS files for Panel1 of Wild Types
  fcsPath <- paste0(data.path,"/WildType/",oneWildType[1,2],"/",Panel,"/")
  dir.create(file.path(fcsPath), recursive = TRUE)
  setwd(fcsPath)
  wildTypeFCS <- sapply(1:nrow(oneWildType), function(i) {
    if(file.exists(oneWildType[i,5])){
      print("FCS file already exists")
    } else{
      download(wildTypeData@fcs.files[[as.integer(oneWildType[i,4])]], only.files="pA.*fcs")
    }
  })


  ##Locating the corrupted FCS files among the Wild Types for both Panels
  print(paste("Finding the Wild Type Corrupted FCS Files in Panel1"))
  setwd("/home/CRC/arahim/Documents/IMPC-Pipeline/")
  fs.path <- paste0(data.path,"/WildType/",oneWildType[1,2],"/",Panel)
  allFCS <- dir(fs.path, full.names = T)
  f <- sapply(1:length(allFCS), function(i){read.FCS(filename = allFCS[i])})
  corruptedFCS <- grep("Corrupted", sapply(1:length(allFCS), function(i){
    f <- try(read.FCS(filename = allFCS[i]), silent = TRUE)
    if(class(f)=="try-error"){
      print(paste0("Corrupted FCS File:", allFCS[i]))
      countCorruptedWT <- countCorruptedWT+1
    }
  }))
  print(paste("Number of Wild Type Corrupted FCS files in Panel1:", length(corruptedFCS)))
  if(length(corruptedFCS) != 0){
    allCorruptedFCS <- t(sapply(1:length(corruptedFCS), function(i){
      allCorruptedFCS <- rbind(allCorruptedFCS, c(unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[2], unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[4], unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[6]))
    }))
    ## removeCorrupted <- sapply(1:length(corruptedFCS), function(i){file.remove(allFCS[corruptedFCS][i])})
  }
  print(paste("Total number of Corrupted FCS files in Panel1 of WT:", countCorruptedWT))


}# end of outer for-loop

allWildTypes <- na.omit(allWildTypes)



##****************************************************************************************************************************************************************************

setwd("/home/CRC/arahim/Documents/IMPC-Pipeline/")


countCorruptedKO <- 0
## Downloading the Knockout Types-----
## Two matrices for Knockout Types:
## 1. allKnockoutTypes contains information for all the Knockouts of a center
## 2. oneKnockoutType contains information for one Knockout at a time
allKnockoutTypes <- matrix(ncol = 6, byrow = TRUE)
colnames(allKnockoutTypes) <- c("Center Code", "ID", "Name", "FCS File Index", "FCS File Name", "md5sum")

knockoutTypeList <- flowRep.ls(impc.only = TRUE, impc.centre = centre, impc.specimen.baseline = FALSE, impc.unanalyzed.only = TRUE)
knockoutTypeData <- NULL

# We are retrieving the Knockout Type dataset
for ( k2 in 1:length(knockoutTypeList) ){
  setwd("/home/CRC/arahim/Documents/IMPC-Pipeline/")
  knockoutTypeData <- flowRep.get(knockoutTypeList[k2], use.credentials=TRUE,impc.details=TRUE)
  if(length(knockoutTypeData@fcs.files) == 0){
    print("The Knockout folder is empty; hence, move to the next Knockout")
    break
  }
  summary(knockoutTypeData)

  ## We want to download one Knockout Line at a time
  ## NOTE: spply() adds a vector to generate the matrix oneKnockoutType. The vectors are added by column. Hence, we had to use the transpose(t()) function after sapply to obtain the final oneKnockoutType matrix.
  oneKnockoutType <- t(sapply(1:length(knockoutTypeData@fcs.files), function(l2){c(centre, knockoutTypeData@id, unlist(strsplit(knockoutTypeData@name," "))[4], l2, knockoutTypeData@fcs.files[[l2]]@name, knockoutTypeData@fcs.files[[l2]]@md5sum)}))
  colnames(oneKnockoutType) <- c("Center Code", "ID", "Name", "FCS File Index", "FCS File Name", "md5sum")

  # Removing duplicate FCS files based on their unique md5sum values
  oneKnockoutType <- oneKnockoutType[!duplicated(oneKnockoutType[,6]),]
  # Storing all the Knockouts in the allKnockoutTypes matrix
  allKnockoutTypes <- rbind(allKnockoutTypes, oneKnockoutType)

  ## Downloading the FCS files for Knockouts of Panel1 (if the number of FCS files is >= 6)
  if(length(grep("pA", oneKnockoutType[,5])) >= 6){
    ## Checks if the FCS files already exists. If not, download the FCS files for Panel1 of Wild Types
    fcsPath <- paste0(data.path,"/KnockoutType/",oneKnockoutType[1,2],"-",oneKnockoutType[1,3],"/",Panel,"/")
    dir.create(file.path(fcsPath), recursive = TRUE)
    setwd(fcsPath) # The working directory is changed so that file.exists() works
    knockoutTypeFCS <- sapply(1:nrow(oneKnockoutType), function(i) {
      if(file.exists(oneKnockoutType[i,5])){
        print("FCS file already exists")
      } else{
        download(knockoutTypeData@fcs.files[[as.integer(oneKnockoutType[i,4])]], only.files="pA.*fcs")
      }
    })
    countFCS.KO.P1 <- countFCS.KO.P1 + length(grep("pA", oneKnockoutType[,5]))
    print(paste("Number of KO files so far", countFCS.KO.P1))
    ## Locating the corrupted FCS files among the Knockout Types for Panel1
      print(paste("Finding the Knockout Type Corrupted FCS Files in Panel1"))
      setwd("/home/CRC/arahim/Documents/IMPC-Pipeline/")
      fs.path <- paste0(data.path,"/KnockoutType/",oneKnockoutType[1,2],"-",oneKnockoutType[1,3],"/",Panel)
      allFCS <- dir(fs.path, full.names = T)
      f <- sapply(1:length(allFCS), function(i){read.FCS(filename = allFCS[i])})
        corruptedFCS <- grep("Corrupted", sapply(1:length(allFCS), function(i){
          f <- try(read.FCS(filename = allFCS[i]), silent = TRUE)
          if(class(f)=="try-error"){
            print(paste0("Corrupted FCS File:", allFCS[i]))
            countCorruptedKO <- countCorruptedKO + 1
          }
        }))
        print(paste("Number of Knockout Type Corrupted FCS files in Panel1:", length(corruptedFCS)))
        if(length(corruptedFCS) != 0){
          allCorruptedFCS <- t(sapply(1:length(corruptedFCS), function(i){
            allCorruptedFCS <- rbind(allCorruptedFCS, c(unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[2], unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[4], unlist(strsplit(as.character(allFCS[corruptedFCS][i]), "/"))[6]))
          }))
          # removeCorrupted <- sapply(1:length(corruptedFCS), function(i){file.remove(allFCS[corruptedFCS][i])})
        }

  } # end of outer if condition


  ## Removing the folder containg the FCS files for each Knockout Type for both Panels
     unlink(paste0("./Center-",centre,"/KnockoutType"), recursive = T)

} # end of outer for-loop

allKnockoutTypes <- na.omit(allKnockoutTypes)
allCorruptedFCS <- na.omit(allCorruptedFCS)

print(paste("Worked with a total of Knockout FCS files:", countFCS.KO.P1))
print(paste("Total number of Corrupted FCS files in Panel1 of KO:", countCorruptedKO))



####*********************************************************************

