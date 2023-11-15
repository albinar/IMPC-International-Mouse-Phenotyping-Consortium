## Created by Albina Rahim
## Date: August 07, 2019
## This is the Main Script which calls the labelFCSFunc function

remove(list=ls())

setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")

source("labelFCSFunc.R")

panel <- "1" ## Change this to panel <- 2 when working on Panel 2

if(panel == 1){
  inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel A"
  outputPath <- "/home/rstudio/results/IMPC/GMC/Panel1"
  markerNames <- c("FSC-H", "FSC-A", "FSC-W", "SSC-A", "FITC gdTCR", "PE CD161", "PI Live/Death", 
                   "PerCP-Cy5.5 CD4", "PE-Cy7 CD62L", "APC CD25", "A700 CD45", "APC-A750 CD8a", "eF450 CD5", "BV570 CD44", "Time")
  panelSuffix <- "Panel_A_"
}else if(panel == 2){
  inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel B"
  outputPath <- "/home/rstudio/results/IMPC/GMC/Panel2"
  markerNames <- c("FSC-H", "FSC-A", "FSC-W", "SSC-A", "FITC CD161", "PE CD21/35", "PI Live/Death", "PerCP-Cy5.5 MHCII",
                   "PE-Cy7 CD19", "APC Ly6G and CD5", "A700 CD45", "APC-eF780 CD11c", "PB CD11b", "BV570 Ly6C", "Time")
  
  panelSuffix <- "Panel_B_"
  
}

labelFCSFunc.Output <- labelFCSFunc(inputPath, outputPath, markerNames, panelSuffix)
