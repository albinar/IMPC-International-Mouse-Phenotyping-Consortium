## Originally written by Albina Rahim 
## Date: December 09, 2016
## This is the Main Script which calls the two functions:
## preProcessingFunc.R : for pre-Processing of the datasets
## globalFrameFunc.R: for creating the Global Frame

remove(list=ls())

#setwd("/code/Projects/IMPC-Universal/Panel1/Codes/")
## Setting the path based on docker
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
##############################################################################################

library("plyr")
library("doMC")


## Function for the pre-Processing of the datasets
source("preProcessingFunc.R")
## Function for creating a Global Frame
source("globalFrameFunc.R")


centreInfo <- NULL

# centreInfo <- rbind(centreInfo, c("Sanger", "L[0-9]+", "Label.Barcode","NA" ,"Assay.Date", "Gender", "Genotype"), 
#                     c("TCP", "[[:alnum:]]{4}[_][0-9]+", "Strain_Code", "Ear_Tag", "Experiment_Date","Gender", "Genotype"), c("CIPHE", "NA", "Strain_Code", "Ear_Tag", "Experiment_Date", "Sex", "Genotype"), c("BCM", "[:digit:]{5,}", "ID", "NA", "Date.of.Sample.Prep", "Sex", "Group"), c("Jax","[:digit:]{4,}", "TestCode", "NA", "NA", "Sex", "NA"))
# colnames(centreInfo) <- c("Centre", "Barcode_Expression", "Barcode/Strain_Code","Barcode/Ear_Tag","Assay_Date", "Gender", "Genotype")

#centreInfo <- rbind(centreInfo, c("Sanger", "L[0-9]+"), c("TCP", "[[:alnum:]]{4}[_][0-9]+"), c("CIPHE", "NA"), c("BCM", "[:digit:]{5,}"), c("Jax","[:digit:]{4,}"),  c("GMC","[:digit:]{8,}"), c("UCD", "[[:alnum:]]{6,}[_][:digit:]{2,}"))
centreInfo <- rbind(centreInfo, c("Sanger", "L[0-9]+"), c("TCP", "[[:alnum:]]{4}[_][0-9]+"), c("CIPHE", "NA"), c("BCM", "[:digit:]{5,}"), c("Jax","[:digit:]{4,}"),  c("GMC","[:digit:]{8,}"), c("UCD", "[^.fcs]+"))## UCDavis has inconsistency in the way they have named the files. So for Barcodes, I am grabbing everything before the .fcs and later in the preProcessingFunc.R I will be removing the _PANEL_A or _PANEL_B suffix for Barcodes
colnames(centreInfo) <- c("Centre", "Barcode_Expression")

barcodeExpression <- centreInfo[which(tolower(centreInfo[,1]) == centre),c('Barcode_Expression')]


## This part of the script was taken from Sibyl for the purpose of parallelizing the execution of this script
## no_cores to determine how many CPUs to use while implenting the parallelization
no_cores <- detectCores() - 2
registerDoMC(no_cores)

start <- Sys.time()

## Paths to the FCS files, metadata spreadsheets, and output folder in the Bioinformatics drive
if(centre == "tcp"){
  inPath1 <- "/home/rstudio/data/IMPC/TCP/TCP_170209/UPLOAD_DUMP"
  inCSVpath1 <- "/home/rstudio/data/IMPC/TCP/TCP_170209/Joined_Flow_2017_04_28.csv"
  inPath2 <- "/home/rstudio/data/IMPC/TCP/TCP_180215/UPLOAD_DUMP"
  inCSVpath2 <- "/home/rstudio/data/IMPC/TCP/TCP_180215/FlowData_2018_02_20.csv"
  
  inputPath <- list(inPath1, inPath2)
  inputCSV <- list(inCSVpath1, inCSVpath2)
  if(panel == 1){
    outputPath <- "/home/rstudio/results/IMPC/TCP/Panel1/"
  }else if(panel == 2){
    outputPath <- "/home/rstudio/results/IMPC/TCP/Panel2/"
  }
}else if(centre == "ciphe"){
  inCSVpath1 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/IMPC_2014_Raw_codeICS-2.csv"
  inCSVpath2 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/IMPC_2015_Raw_codeICSV3.csv"
  inCSVpath3 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/IMPC_2016_Raw.csv"
  inputCSV <- list(inCSVpath1, inCSVpath2, inCSVpath3)
  if(panel == 1){
    inPath1 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2014/IMPC1"
    inPath2 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2015/IMPC1"
    inPath3 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2016/IMPC1"
    inputPath <- list(inPath1, inPath2, inPath3)
    outputPath <- "/mnt/f/FCS data/IMPC/IMPC-Results/CIPHE/Panel1/"
  }else if(panel == 2){
    inPath1 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2014/IMPC2"
    inPath2 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2015/IMPC2"
    inPath3 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/CIPHE/CIPHE/IMPC_CIPHE_2016/IMPC2"
    inputPath <- list(inPath1, inPath2, inPath3)
    outputPath <- "/mnt/f/FCS data/IMPC/IMPC-Results/CIPHE/Panel2/"
  }
}else if(centre == "bcm"){
  inPath1 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/BCM/BCM"
  inCSVpath1 <- "/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/BCM/immunophenotyping_experiments_short_3.csv"
  inputPath <- list(inPath1)
  inputCSV <- list(inCSVpath1)
  if(panel == 1){
    outputPath <- "/mnt/f/FCS data/IMPC/IMPC-Results/BCM/Panel1/"
  }else if(panel == 2){
    outputPath <- "/mnt/f/FCS data/IMPC/IMPC-Results/BCM/Panel2/"
  }
}else if(centre == "jax"){
  inPath1 <- "/home/rstudio/data/IMPC/JAX/Jax_2017-03-29/PKG - JAXKOMP data files"
  inputPath <- list(inPath1)
  inCSVpath1 <- "/home/rstudio/data/IMPC/JAX/Jax_2017-03-29/Jax_metadata_20170331.csv"
  inputCSV <- list(inCSVpath1)
  if(panel == 1){
    outputPath <- "/home/rstudio/results/IMPC/JAX/Panel1/"
  }else if(panel == 2){
    outputPath <- "/home/rstudio/results/IMPC/JAX/Panel2/"
  }
}else if(centre == "gmc"){
  
  if(panel == 1){
    inPath1 <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_A"
    outputPath <- "/home/rstudio/results/IMPC/GMC/Panel1/"
    inCSVpath1 <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_A/Metadata_GMC_2013-2018_PanelA.csv"
    
  }else if(panel == 2){
    inPath1 <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_B"
    outputPath <- "/home/rstudio/results/IMPC/GMC/Panel2/"
    inCSVpath1 <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_B/Metadata_GMC_2013-2018_PanelB.csv"
    
  }
  
  
  inputPath <- list(inPath1)
  inputCSV <- list(inCSVpath1)
  
}else if(centre == "ucd"){
  if(panel == 1){
    inPath1 <- "/home/rstudio/data/IMPC/UCD/FCS_Files/Panel_A"
    outputPath <- "/home/rstudio/results/IMPC/UCD/Panel1/"
  }else if(panel == 2){
    inPath1 <- "/home/rstudio/data/IMPC/UCD/FCS_Files/Panel_B"
    outputPath <- "/home/rstudio/results/IMPC/UCD/Panel2/"
    #inCSVpath1 <- "/home/rstudio/data/IMPC/UCD/Metadata/Combined_Metadata_2020.csv"
    
  }
  
  inCSVpath1 <- "/home/rstudio/data/IMPC/UCD/Metadata/Combined_Metadata_2020.csv"
  
  inputPath <- list(inPath1)
  inputCSV <- list(inCSVpath1)
}


## Calling the preProcessing function
preProcessing.Output <- preProcessingFunc(inputPath, inputCSV, centre, barcodeExpression)
store.allFCS <- preProcessing.Output$store.allFCS
store.allFMO <- preProcessing.Output$store.allFMO
NObarcodes.FCS <- preProcessing.Output$NObarcodes.FCS
Genotype <- preProcessing.Output$Genotype
uniqueGT <- preProcessing.Output$uniqueGT
notListed.FCS <- preProcessing.Output$notListed.FCS
Barcodes.NoFCS <- preProcessing.Output$Barcodes.NoFCS
corrupted.FCS <- preProcessing.Output$corrupted.FCS
lessCells.FCS <- preProcessing.Output$lessCells.FCS
duplicate.FCS <- preProcessing.Output$duplicate.FCS
numChannels <- preProcessing.Output$numChannels

panel1.Names <- c("Panel1", "IMPC1", "PANEL_A", "Tmem panel")
panel2.Names <- c("Panel2", "IMPC2", "PANEL_B", "APC panel")

if(centre == "tcp" | centre == "jax"){
  panel1.Matches <- grep(paste(panel1.Names,collapse="|"), store.allFCS[,c('FCS files')], value=FALSE)
  panel2.Matches <- grep(paste(panel2.Names,collapse="|"), store.allFCS[,c('FCS files')], value=FALSE)
  store.allFCS.panel1 <- store.allFCS[panel1.Matches,]
  store.allFCS.panel2 <- store.allFCS[panel2.Matches,]
  
  panel1.FMO.Names <- c("FMO_PANEL_A")
  panel2.FMO.Names <- c("FMO_PANEL_B")
  panel1.FMO.Matches <- grep(paste(panel1.FMO.Names,collapse="|"), store.allFMO[,c('FMO')], value=FALSE)
  panel2.FMO.Matches <- grep(paste(panel2.FMO.Names,collapse="|"), store.allFMO[,c('FMO')], value=FALSE)
  store.allFMO.panel1 <- store.allFMO[panel1.FMO.Matches,]
  store.allFMO.panel2 <- store.allFMO[panel2.FMO.Matches,]
  
  # ## NOTE: Temporary action since the compensation matrix for the folder "16-0602 KOMP" is missing for Panel 1
  # index.remove.temp <- which(store.allFCS.panel1[,c('Panel/Organ/Folder')] == "16-0602  KOMP")
  # store.allFCS.panel1 <- store.allFCS.panel1[-index.remove.temp,]
  
  ## Writing the Summary of the Pre-Processing output in a text file
  print("Printing the Summary of the Pre-Processing Output: ")
  suppressWarnings(dir.create (paste0(outputPath,"Results/")))
  
  sink(file = paste0(outputPath,"Results/preProcessing-Summary.txt"), split = TRUE)
  print(paste0("There are in total ", nrow(store.allFCS), " Panels 1 & 2 files for analysis."))
  print(paste0("There are in total ", nrow(store.allFCS.panel1), " Panels 1 files for analysis."))
  print(paste0("There are in total ", nrow(store.allFCS.panel2), " Panels 2 files for analysis."))
  print(paste0("There are in total ", nrow(store.allFMO), " FMOs for Panels 1 & 2."))
  print(paste0("There are in total ", nrow(store.allFMO.panel1), " FMOs for Panels 1."))
  print(paste0("There are in total ", nrow(store.allFMO.panel2), " FMOs for Panels 2."))
  print(paste0("Number of FCS files with NO Barcodes: ", NObarcodes.FCS))
  print(paste0("Number of FCS files which were sent to us but whose information was missing in the metadata spreadsheet: ", nrow(notListed.FCS)))
  print(paste0("Number of Barcodes which were there in the metadata spreadsheet but for which we received no FCS files: ", length(Barcodes.NoFCS)))
  print(paste0("Number of Corrupted files: ", nrow(corrupted.FCS)))
  print(paste0("Number of files with < 20,000 cells: ", nrow(lessCells.FCS)))
  print(paste0("Number of Duplicate FCS files: ", nrow(duplicate.FCS)))
  if(length(numChannels) == 1){
    print(paste0("All files have the same number of channels: ", numChannels))
  } else{
    #paste(c("The first three notes are: ", notes), collapse=" ")
    print(paste0(c("Files have different number of channels:", numChannels), collapse = " "))
  }
  sink()
  
  if(panel == 1){
    store.allFCS <- store.allFCS.panel1
    store.allFMO <- store.allFMO.panel1
  }else if(panel == 2) {
    store.allFCS <- store.allFCS.panel2
    store.allFMO <- store.allFMO.panel2
  }
  
}else if (centre == "sanger" | centre == "ciphe" | centre == "bcm" | centre == "gmc" | centre == "ucd"){
  ## Writing the Summary of the Pre-Processing output in a text file
  print("Printing the Summary of the Pre-Processing Output: ")
  suppressWarnings(dir.create (paste0(outputPath,"Results/")))
  sink(file = paste0(outputPath,"Results/preProcessing-Summary.txt"), split = TRUE)
  print(paste0("There are in total ", nrow(store.allFCS), " files for analysis."))
  print(paste0("Number of FCS files with NO Barcodes: ", NObarcodes.FCS))
  print(paste0("Number of FCS files which were sent to us but whose information was missing in the metadata spreadsheet: ", nrow(notListed.FCS)))
  print(paste0("Number of Barcodes which were there in the metadata spreadsheet but for which we received no FCS files: ", length(Barcodes.NoFCS)))
  print(paste0("Number of Corrupted files: ", nrow(corrupted.FCS)))
  print(paste0("Number of files with < 20,000 cells: ", nrow(lessCells.FCS)))
  print(paste0("Number of Duplicate FCS files: ", nrow(duplicate.FCS)))
  if(length(numChannels) == 1){
    print(paste0("All files have the same number of channels: ", numChannels))
  } else{
    #paste(c("The first three notes are: ", notes), collapse=" ")
    print(paste0(c("Files have different number of channels:", numChannels), collapse = " "))
  }
  sink()
}

## Calling the Function for creating the Global Frame
globalFrame.Output <- globalFrameFunc(store.allFCS, centre, panel, outputPath)
## Saving the Global Frame
globalFrame <- globalFrame.Output$globalFrame
## Saving the Expression Matrix of the Global Frame with the information of each FCS file
## This matrix can later be used if we need to remove any FCS files and its corresponding expression matrix values from the global frame
globalFrame.Matrix <- globalFrame.Output$globalFrame.Matrix
##Saving the computed transform from estimateLogicle()
lgl <- globalFrame.Output$lgl
## Saving the list of files which had NAs in some of their channels. These files were already removed from the global frame
files.w.NAs <- globalFrame.Output$files.w.NAs


save (store.allFCS, file =  paste0(outputPath,"Results/store.allFCS.Rdata") )
save (store.allFMO, file =  paste0(outputPath,"Results/store.allFMO.Rdata") )
save(Genotype, file = paste0(outputPath,"Results/Genotype.Rdata"))
save(uniqueGT, file = paste0(outputPath,"Results/uniqueGT.Rdata"))
save(notListed.FCS, file = paste0(outputPath, "Results/notListed.FCS.Rdata"))
save(Barcodes.NoFCS, file = paste0(outputPath, "Results/Barcodes.NoFCS.Rdata"))
save(corrupted.FCS , file = paste0(outputPath, "Results/corrupted.FCS.Rdata"))
save(lessCells.FCS, file = paste0(outputPath,"Results/lessCells.FCS.Rdata"))
save(duplicate.FCS, file = paste0(outputPath,"Results/duplicate.FCS.Rdata"))
save(globalFrame, file = paste0(outputPath,"Results/globalFrame.Rdata"))
save(globalFrame.Matrix, file = paste0(outputPath, "Results/globalFrame.Matrix.Rdata"))
save(lgl, file = paste0(outputPath, "Results/lgl.Rdata"))
save(files.w.NAs, file = paste0(outputPath, "Results/files.w.NAs.Rdata"))

cat("Total time is: ",TimeOutput(start),sep="")
