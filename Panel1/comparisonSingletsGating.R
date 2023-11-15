## Albina Rahim
## August 19, 2019

## This is a temporary code written for comparing the population proportions if we use FSC-H versus SSC-A for 
## SSC Singlets gating. It was observed that for centre GMC, the SSC-H channel is missing. Therefore we are testing 
## the alternate approach on TCP. The comparison is done on 100 FCS files


remove(list=ls())

#setwd("/code/Projects/IMPC-Universal/Panel1/")
## Setting the path based on docker
setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")
## This version of the script identifies the Worst 30 outliers and random 30  between Manual and Automated results

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

#source("cvFunc.R")

library('colorRamps')
library('plyr')
library('doMC')
library('tools')
library('MASS')
#library('flowViz')

no_cores <- detectCores() - 2
registerDoMC(no_cores)


start <- Sys.time()

# if(centre == "tcp"| centre == "ciphe" | centre == "bcm"){
#   results.dir <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/", toupper(centre), "/Panel", panel, "/Results")
# }else if(centre == "jax"){
#   results.dir <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/", toTitleCase(centre), "/Panel", panel, "/Results")
# }

if(centre == "tcp"| centre == "ciphe" | centre == "bcm"){
  results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")
}else if(centre == "jax"){
  results.dir <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/", toTitleCase(centre), "/Panel", panel, "/Results")
}



load(paste0(results.dir,"/all.events.store.Live-Singlets.Rdata"))
load(paste0(results.dir,"/all.props.store.Live-Singlets.Rdata"))

all.props.store.Live <- all.props.store
all.events.store.Live <- all.events.store

load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/all.events.store.Rdata"))
load(paste0(results.dir,"/all.props.store.Rdata"))

all.props.store.temp <- all.props.store[1:200,]
all.events.store.temp <- all.events.store[1:200,]


for(i in 15:ncol(all.events.store.Live)){
  result <- matrix(nrow = nrow(all.props.store.Live), ncol = 3, data = NA)# matrix for saving the  proportions for comparison
  
  result[,1] <- all.events.store.Live[,c('FCS files')]
  result[,2] <- all.events.store.temp[,i]
  result[,3] <- all.events.store.Live[,i]
  
  colnames(result) <- c('FCS Files', paste0("Original-",colnames(all.events.store.temp)[i]),  paste0("Modified-",colnames(all.events.store.Live)[i]))
  resultDF <- as.data.frame(result)
  
  ## Calculating the line of regression for the scatter plot
  reg1 <- lm(as.numeric(resultDF[[2]])~as.numeric(resultDF[[3]]),data=resultDF) 
  #summary(reg1)
  
  ## Identifying the worst 30 outliers
  resultDF$residuals <- reg1$residuals 
  outliers <- order(resultDF$residuals^2, decreasing=T)
  worst30.outliers <- outliers[1:30]
  resultDF$outliers <- 0
  for(j in 1: length(worst30.outliers)){
    resultDF$outliers[worst30.outliers[j]] <- 1
  }
  
  
  png ( file = paste0(results.dir,"/comparisonSinglets-plots/",colnames(all.events.store.temp)[i],".png"))
  plot(as.numeric(resultDF[[3]]), as.numeric(resultDF[[2]]), main=paste0(colnames(all.events.store.temp)[i],": Original vs Modified"),
       xlab="Modified Event Counts", ylab="Original Event Counts", pch=19, col = ifelse(resultDF$outliers == 1, "red", "blue")) 
  abline(reg1, col="black", lwd=4) 
  
  dev.off()
}










