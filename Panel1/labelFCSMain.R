## Created by Albina Rahim
## Date: August 07, 2019
## This is the Main Script which calls the labelFCSFunc function

remove(list=ls())

##This the working directory where the R scripts (both the Main and Function) have been saved
## Change this to your working directory.
setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")

## Loading the function that will process the raw files
source("labelFCSFunc.R")

## Specifying the panel number.
## Change this to panel<-"2" when working on Panel 2
panel <- "1" 

if(panel == 1){
  ##Change the inputPath to the directory where the raw files have been kept for Panel1
  inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel A"
  
  ## Input path for the FMO files
  inputPath.FMO <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/FMO_20180319_KW12_2019/Panel_A_FMO"
  
  ##Change the outputPath to the directory where the processed FCS files will be saved for Panel1
  outputPath <- "/home/rstudio/results/IMPC/GMC/Panel1"
  
  ##Change the outputPath to the directory where the processed FMO files will be saved for Panel1
  outputPath.FMO <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_A/FMO_20180319_KW12_2019"
  
  ## Marker names for Panel1. No need to change anything
  markerNames <- c("FSC-H", "FSC-A", "FSC-W", "SSC-A", "FITC gdTCR", "PE CD161", "PI Live/Death", 
                   "PerCP-Cy5.5 CD4", "PE-Cy7 CD62L", "APC CD25", "A700 CD45", "APC-A750 CD8a", "eF450 CD5", "BV570 CD44", "Time")
  
  panelPrefix <- "Panel_A_"
  
  fmoPrefix <- "FMO_"
  
}else if(panel == 2){
  ##Change the inputPath to the directory where the raw files have been kept for Panel2
  inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel B"
  
  ## Input path for the FMO files
  inputPath.FMO <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/FMO_20180319_KW12_2019/Panel_B_FMO"
  
  ##Change the outputPath to the directory where the processed FCS files will be saved for Panel2
  outputPath <- "/home/rstudio/results/IMPC/GMC/Panel2"
  
  ##Change the outputPath to the directory where the processed FMO files will be saved for Panel1
  outputPath.FMO <- "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_B/FMO_20180319_KW12_2019"
  
  
  ## Marker names for Panel1. No need to change anything
  markerNames <- c("FSC-H", "FSC-A", "FSC-W", "SSC-A", "FITC CD161", "PE CD21/35", "PI Live/Death", "PerCP-Cy5.5 MHCII",
                   "PE-Cy7 CD19", "APC Ly6G and CD5", "A700 CD45", "APC-eF780 CD11c", "PB CD11b", "BV570 Ly6C", "Time")
  
  panelPrefix <- "Panel_B_"
  
  fmoPrefix <- "FMO_"
}

## Calling the labelFCSFunc function for analyzing the raw FCS files
labelFCSFunc.Output <- labelFCSFunc(inputPath, outputPath, markerNames, panelPrefix)

## Calling the labelFCSFunc function for analyzing the raw FMO files
labelFCSFunc.Output.FMO <- labelFCSFunc(inputPath = inputPath.FMO, outputPath = outputPath.FMO, markerNames, panelPrefix = fmoPrefix)
