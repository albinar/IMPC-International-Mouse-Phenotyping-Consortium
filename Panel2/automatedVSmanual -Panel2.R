## Written by Albina Rahim
remove(list=ls())

setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel2")

## This version of the script is for Panel 2
## It identifies the Worst 30 outliers and random 30 between Manual and Automated results

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
library('colorRamps')
library('plyr')
library('doMC')
library('tools')
library('MASS')

no_cores <- detectCores() - 2
registerDoMC(no_cores)


start <- Sys.time()

if(centre == "tcp"| centre == "ciphe" | centre == "bcm"){
  #results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results-Original-March2019")
  results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel)
  
  #results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")
}else if(centre == "jax"){
  results.dir <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/", toTitleCase(centre), "/Panel", panel, "/Results")
}



all.events.store <- read.csv(paste0(results.dir,"/EventCounts_TCP_Panel2_20190710_1858.csv"))

NoFCS.temp <- NULL
#WTindex <- which(store.allFCS[,c('Genotype')] == "WT")

manualResults <- read.csv(paste0(results.dir,"/manualResults.csv"))
allFiles.metadata <- manualResults[,c('FCS_Files')]  
# duplicate.index <- which(duplicated(allFiles.metadata)==TRUE)
# if(length(duplicate.index) != 0){
#   allBarcodes.metadata <- allBarcodes.metadata[-duplicate.index]
# }

NoFCS.temp <- c(NoFCS.temp, setdiff(allFiles.metadata, all.events.store[,c('FCS.files')]))

NoFCS.index <- sapply(1:length(NoFCS.temp), function(x){which(NoFCS.temp[x] == manualResults[,c('FCS_Files')])})
manualResults <- manualResults[-NoFCS.index,]


manualResults <- as.matrix(manualResults)
rownames(manualResults) <- 1:nrow(manualResults)


autoVSmanual.store <- NULL
indexMatch <- 0
for(i in 1:nrow(all.events.store)){
  indexMatch <- which(manualResults[i,c('FCS_Files')] == all.events.store[,c('FCS.files')])
  autoVSmanual.store <- c(autoVSmanual.store, manualResults[indexMatch,c('FCS_Files')])
  
}


manualIndex <- 0
autoIndex <- 0
result <- matrix(nrow = length(autoVSmanual.store), ncol = 3, data = NA)# matrix for saving the event counts
result[,1] <- autoVSmanual.store
for(i in 1:length(autoVSmanual.store)){
  manualIndex[i] <- which(result[i,1] == manualResults[,c('FCS_Files')])
  #result[manualIndex,1] <- manualResults[manualIndex,c('FCS_Files')]
  result[i,2] <- (as.numeric(manualResults[manualIndex[i], c('Monocytes')])/as.numeric(manualResults[manualIndex[i], c('All_Events')]))*100
  
  autoIndex[i] <-   which(result[i,1] == all.events.store[,c('FCS.files')])
  result[i,3] <- (as.numeric(all.events.store[autoIndex[i], c('Monocytes')])/as.numeric(all.events.store[autoIndex[i], c('All.Events')]))*100
  
}

colnames(result) <- c('FCS Files', 'Manual', 'Automated')
resultDF <- as.data.frame(result)
#write.csv(result, paste0(results.dir,"/results.csv"))


## Calculating the line of regression for the scatter plot
reg1 <- lm(as.numeric(resultDF$Manual)~as.numeric(resultDF$Automated),data=resultDF) 
summary(reg1)

## Identifying the worst 30 outliers
resultDF$residuals <- reg1$residuals 
outliers <- order(resultDF$residuals^2, decreasing=T)
worst30.outliers <- outliers[1:30]
resultDF$outliers <- 0
for(i in 1: length(worst30.outliers)){
  resultDF$outliers[worst30.outliers[i]] <- 1
}


plot(as.numeric(resultDF$Automated), as.numeric(resultDF$Manual), main="Monocytes: Automated vs Manual (worst 30 outliers)",
     xlab="Automated Proportion (% of All Events)", ylab="Manual Proportion (% of All Events)", pch=19, col = ifelse(resultDF$outliers == 1, "red", "blue")) 
abline(reg1, col="black", lwd=4) 
# ## Labeling the points
# with(resultDF[worst30.outliers,], text(Manual~Automated, labels = row.names(resultDF[worst30.outliers,]), pos = 4))

## Saving the Worst 30 Outliers in a spreadsheet (with the first one being the worst among the 30)
store.worst30.outliers<- NULL
store.worst30.outliers <- manualResults[worst30.outliers,c('Experiment_Date', 'Genotype', 'Sex', 'Mouse_ID', 'FCS_Files')]
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(store.worst30.outliers, file =  paste0(results.dir, "/Monocytes.Worst.30.Outliers_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

## Copy pasting the gating plots for the Worst 30 outliers
for(i in 1: nrow(store.worst30.outliers)){
  fcsName <-sub(".fcs","", store.worst30.outliers[i,c('FCS_Files')])
  print(paste0("Copy pasting", fcsName))
  file.copy(paste0(results.dir,"/Figures/Scatterplots/", store.worst30.outliers[i,c('Experiment_Date')],"_",fcsName, ".png"), paste0(results.dir,"/autoVSmanual-plots/Monocytes/Worst30/"))
}



## Random 30 files
random30 <- sample(1:nrow(resultDF), 30, replace = F)
#random30.outliers <- outliers[random30.outliers.index]
resultDF$Random <- 0
for(i in 1: length(random30)){
  resultDF$Random[outliers[random30[i]]] <- 1
}


plot(as.numeric(resultDF$Automated), as.numeric(resultDF$Manual), main="Monocytes: Automated vs Manual (random 30 files)",
     xlab="Automated Proportion (% of All Events)", ylab="Manual Proportion (% of All Events)", pch=19, col = ifelse(resultDF$Random == 1, "red", "blue")) 
abline(reg1, col="black", lwd=4) 

## Saving the Random 30 in a spreadsheet 
store.random30 <- NULL
store.random30 <- manualResults[outliers[random30],c('Experiment_Date', 'Genotype', 'Sex', 'Mouse_ID', 'FCS_Files')]
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(store.random30, file =  paste0(results.dir, "/Monocytes.Random.30_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

## Copy pasting the gating plots for the Random 30 files
for(i in 1: nrow(store.random30)){
  fcsName <-sub(".fcs","", store.random30[i,c('FCS_Files')])
  file.copy(paste0(results.dir,"/Figures/Scatterplots/", store.random30[i,c('Experiment_Date')],"_",fcsName, ".png"), paste0(results.dir,"/autoVSmanual-plots/Monocytes/Random30/"))
}



# Saving the Result dataframe of Proportions (as % of All Events, Worst 30 outliers and Random 30 Outliers) 
save(resultDF, file = paste0(results.dir,"/Monocytes.automatedVSmanual.resultsDF.Rdata"))  

###################################################################################

## Adding the Event Counts to the worst30 and random30 spreadsheets

worst30CSV <- as.matrix(read.csv(paste0(results.dir, "/Monocytes.Worst.30.Outliers_TCP_Panel2_20190710_2204.csv")))
worst30CSV <- cbind(worst30CSV, NA)
worst30CSV <- cbind(worst30CSV, NA)

manualIndex <-0
autoIndex <- 0
for(i in 1:nrow(worst30CSV)){
  manualIndex[i] <- which(worst30CSV[i,5] == manualResults[,c('FCS_Files')])
  #result[manualIndex,1] <- manualResults[manualIndex,c('FCS_Files')]
  worst30CSV[i,6] <- (as.numeric(manualResults[manualIndex[i], c('Monocytes')]))
  
  autoIndex[i] <-   which(worst30CSV[i,5] == all.events.store[,c('FCS.files')])
  worst30CSV[i,7] <- (as.numeric(all.events.store[autoIndex[i], c('Monocytes')]))
  
}

colnames(worst30CSV)[6] <- "Event counts-Manual"
colnames(worst30CSV)[7] <- "Event counts-Auto"
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(worst30CSV, file =  paste0(results.dir, "/Monocytes.Worst.30.Outliers_",toupper(centre),"_Panel2", date.time), row.names = FALSE)


random30CSV <- as.matrix(read.csv(paste0(results.dir, "/Monocytes.Random.30_TCP_Panel2_20190710_2204.csv")))
random30CSV <- cbind(random30CSV, NA)
random30CSV <- cbind(random30CSV, NA)

for(i in 1:nrow(random30CSV)){
  manualIndex[i] <- which(random30CSV[i,5] == manualResults[,c('FCS_Files')])
  #result[manualIndex,1] <- manualResults[manualIndex,c('FCS_Files')]
  random30CSV[i,6] <- (as.numeric(manualResults[manualIndex[i], c('Monocytes')]))
  
  autoIndex[i] <-   which(random30CSV[i,5] == all.events.store[,c('FCS.files')])
  random30CSV[i,7] <- (as.numeric(all.events.store[autoIndex[i], c('Monocytes')]))
  
}

colnames(random30CSV)[6] <- "Event counts-Manual"
colnames(random30CSV)[7] <- "Event counts-Auto"
date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(random30CSV, file =  paste0(results.dir, "/Monocytes.Random.30_",toupper(centre),"_Panel2", date.time), row.names = FALSE)


########################################################

## Date: Nov 17. 2020
## Comparing Worst 30 and Random 30 of TCP files based on the gating strategy Version V 2.0

## WORST30
autoResults.worst30 <-  as.matrix(read.csv(paste0(results.dir,"/Automated_Results/Worst30/Eosinophils.Worst.30.Outliers_TCP_Panel2_20190710_2148.csv")))
manualResults.worst30 <- as.matrix(read.csv(paste0(results.dir,"/Manual_Results/Worst30/PanelB_Event_Counts_Eosinophils_Worst.csv")))

#colnames(autoResults.worst30)
#colnames(manualResults.worst30)

for(i in 1:23){
  result <- matrix(nrow = nrow(manualResults.worst30), ncol = 3, data = NA)# matrix for saving the  event counts for comparison
  
  result[,1] <- manualResults.worst30[,c('FCS_Files')]
  result[,2] <- autoResults.worst30[,i+2] ## Starting from Live population
  result[,3] <- manualResults.worst30[,i+2]
  
  colnames(result) <- c('FCS Files', paste0("Automated-",colnames(manualResults.worst30)[i+2]),  paste0("Manual-",colnames(manualResults.worst30)[i+2]))
  resultDF <- as.data.frame(result)
  
  ## Calculating the line of regression for the scatter plot
  reg1 <- lm(as.numeric(resultDF[[2]])~as.numeric(resultDF[[3]]),data=resultDF)
  #summary(reg1)
  
  ## Identifying the worst 30 outliers
  resultDF$residuals <- reg1$residuals
  outliers <- order(resultDF$residuals^2, decreasing=T)
  worst30.outliers <- outliers[1:29] ## replacing outliers[1:30]
  resultDF$outliers <- 0
  for(j in 1: length(worst30.outliers)){
    resultDF$outliers[worst30.outliers[j]] <- 1
  }
  
  
  png ( file = paste0(results.dir,"/Comparison/Worst30/",colnames(manualResults.worst30)[i+2],".png"))
  plot(as.numeric(resultDF[[3]]), as.numeric(resultDF[[2]]), main=paste0(colnames(manualResults.worst30)[i],": Automated vs Manual"),
       xlab="Automated Event Counts", ylab="Manual Event Counts", pch=19, col = ifelse(resultDF$outliers == 1, "red", "blue"))
  abline(reg1, col="black", lwd=4)
  
  dev.off()
}

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(autoResults.worst30, file =  paste0(results.dir, "/Comparison/Worst30/Eosinophils.Worst.30_",toupper(centre),"_Panel2", date.time), row.names = FALSE)

################
## Accessing the FlowJo files for the WORST30
library("flowCore")
library("CytoML")
library("flowWorkspace")
library("flowDensity")
ws <- "/home/rstudio/results/IMPC/TCP/Panel2/Manual_Results/Worst30/PanelB_analysis_Eosinophils_worst.wsp"
ws.file <- open_flowjo_xml(ws)
gs <- flowjo_to_gatingset(ws.file, path = "/home/rstudio/results/IMPC/TCP/Panel2/Manual_Results/Worst30")


######################
######################

## RANDOM30
autoResults.random30 <-  as.matrix(read.csv(paste0(results.dir,"/Automated_Results/Random30/Eosinophils.Random.30_TCP_Panel2_20190710_2149.csv")))
manualResults.random30 <- as.matrix(read.csv(paste0(results.dir,"/Manual_Results/Random30/PanelB_Event_Counts_Eosinophils_Random.csv")))

#colnames(autoResults.random30)
#colnames(manualResults.random30)

for(i in 1:23){
  result <- matrix(nrow = nrow(manualResults.random30), ncol = 3, data = NA)# matrix for saving the  event counts for comparison
  
  result[,1] <- manualResults.random30[,c('FCS_Files')]
  result[,2] <- autoResults.random30[,i+2] ## Starting from Live population
  result[,3] <- manualResults.random30[,i+2]
  
  colnames(result) <- c('FCS Files', paste0("Automated-",colnames(manualResults.random30)[i+2]),  paste0("Manual-",colnames(manualResults.random30)[i+2]))
  resultDF <- as.data.frame(result)
  
  ## Calculating the line of regression for the scatter plot
  reg1 <- lm(as.numeric(resultDF[[2]])~as.numeric(resultDF[[3]]),data=resultDF)
  #summary(reg1)
  
  ## Identifying the random 30 outliers
  resultDF$residuals <- reg1$residuals
  outliers <- order(resultDF$residuals^2, decreasing=T)
  random30.outliers <- outliers[1:27] ## replacing outliers[1:30]
  resultDF$outliers <- 0
  for(j in 1: length(random30.outliers)){
    resultDF$outliers[random30.outliers[j]] <- 1
  }
  
  
  png ( file = paste0(results.dir,"/Comparison/Random30/",colnames(manualResults.random30)[i+2],".png"))
  plot(as.numeric(resultDF[[3]]), as.numeric(resultDF[[2]]), main=paste0(colnames(manualResults.random30)[i],": Automated vs Manual"),
       xlab="Automated Event Counts", ylab="Manual Event Counts", pch=19, col = ifelse(resultDF$outliers == 1, "red", "blue"))
  abline(reg1, col="black", lwd=4)
  
  dev.off()
}

date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
write.csv(autoResults.random30, file =  paste0(results.dir, "/Comparison/Random30/Eosinophils.Random.30_",toupper(centre),"_Panel2", date.time), row.names = FALSE)


#############################################
## Accessing the FlowJo files for the RANDOM30
library("flowWorkspace")
library("CytoML")
ws <- "/home/rstudio/results/IMPC/TCP/Panel2/Manual_Results/Random30/PanelB_analysis_Eosinophils_random.wsp"
ws.file <- open_flowjo_xml(ws)
gs <- flowjo_to_gatingset(ws.file, path = "/home/rstudio/results/IMPC/TCP/Panel2/Manual_Results/Random30")
