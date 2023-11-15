##Testing F1 score for COVID data
##Written by : Albina & Mehrnoush
##Modified: April, 2020

##I looked at UCD panel1, and I think CD8 tcells has 3 gates, 
##so we can use that one.IF you think another popu;ation is better, then do that one.
#Check the code for one file first, to make sure there's no type, as I didn't run anything :-D

remove(list=ls())

## Setting the path based on docker
setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")

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
library('e1071') # for flowCut
library('Cairo') # for flowCut
library('flowCore')
library('flowDensity')
library('pracma') # for findpeaks
library('tools')
library('MASS')
library('stringr')## for str_match used in compensateIMPC function
library('flowViz')
library('knitr')

no_cores <- detectCores() - 3
registerDoMC(no_cores)

results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")

suppressWarnings(dir.create(paste0(results.dir,"/CD8TcellsCLRs-10K")))

## Loading the gating thresholds, store matrix, and channels.ind
load(paste0(results.dir,"/all.gthres.store.Rdata"))
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir, "/channels.ind.Rdata"))

## I am manually removing "PANEL_A_ABDQ_67_C4004.fcs" file 
manual.remove.Ind <- which(store.allFCS[,c('FCS files')] == "PANEL_A_ABDQ_67_C4004.fcs")
store.allFCS <- store.allFCS[-manual.remove.Ind,]
rownames(store.allFCS) <- 1:nrow(store.allFCS)



N <- 10000 # 200000

store.allFCS.cd8Tcells <- store.allFCS
store.allFCS.cd8Tcells[,c('Path')]<- "/home/rstudio/results/IMPC/TCP/Panel1/Results/CD8Tcells"
file.names <- data.frame(store.allFCS.cd8Tcells, stringsAsFactors = F)

get.data <- llply(1:200, function(i){ 
#get.data <- llply(1:nrow(file.names), function(i){ 
  #props.events.gates <- llply(1:length(index), function(i){ 
  #props.events.gates <- llply(1:200, function(i){ 
  
  x<- file.names[i,]
  
  load(paste0(x$Path, "/", x$FCS.files, ".Rdata"))
  #frame <- cd8.Tcells
  
  
  # Gating CD8 T-cells to obtain CD8 Effector cells and CD8 Resting/Naive cells
  cd62.gate.Tregs <- as.numeric(all.gthres.store[i, c('cd62.gate.Tregs')])
  cd44.gate.T.helper <- as.numeric(all.gthres.store[i,c('cd44.gate.Thelper')])
  cd44.gate.cd8.Tcells <- as.numeric(all.gthres.store[i,c('cd44.gate.cd8Tcells')])
  
  Effector.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(F,T), gates = c(cd62.gate.Tregs, cd44.gate.T.helper))
  
  Resting.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,T), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
  
  Naive.cd8.flowD <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), position = c(T,F), gates = c(cd62.gate.Tregs, cd44.gate.cd8.Tcells))
  
  # Effector.cd8 <- getflowFrame(Effector.cd8.flowD)
  # 
  # Resting.cd8 <- getflowFrame(Resting.cd8.flowD)
  # 
  # Naive.cd8 <- getflowFrame(Naive.cd8.flowD)
  
  
  # filter1 <- Effector.cd8.flowD@filter
  # filter2 <- Resting.cd8.flowD@filter
  # filter3 <- Naive.cd8.flowD@filter
  
  
  
  if(nrow(cd8.Tcells)>N){
    frame <- cd8.Tcells[sample (1:nrow(cd8.Tcells), N), c(channels.ind["CD62L"], channels.ind["CD44"])]
  }else{
    frame <- cd8.Tcells[sample (1:nrow(cd8.Tcells), N,replace=T), c(channels.ind["CD62L"], channels.ind["CD44"])]
  }
  
  clr1 <- c(1:nrow(frame)) %in% Effector.cd8.flowD@index 
  save(clr1, file=paste0(results.dir,"/CD8TcellsCLRs-10K/", x$FCS.files, "_Pop1.RData"))
  clr2 <- c(1:nrow(frame)) %in% Resting.cd8.flowD@index
  save(clr2, file=paste0(results.dir,"/CD8TcellsCLRs-10K/", x$FCS.files, "_Pop2.RData"))
  clr3 <- c(1:nrow(frame)) %in% Naive.cd8.flowD@index
  save(clr3, file=paste0(results.dir,"/CD8TcellsCLRs-10K/", x$FCS.files, "_Pop3.RData"))
  #cd8.Tcells <- f
  
  list(x$FCS.files)
  
}, .parallel = TRUE) # end llply


get.data <- function(frame=cd8.Tcells,
                     channels.ind,
                     name=x$FCS.files,
                     filter1=Effector.cd8.flowD@filter,
                     filter2=Resting.cd8.flowD@filter,
                     filter3=Naive.cd8.flowD@filter)
#Arguments: flowFrame of the parent population, and three gates
  
  {
  if (nrow(frame)>N)
  {

    f <- frame[sample (1:nrow(frame), N), c(channels.ind["CD62L"], channels.ind["CD44"])]
  }else{
    f <- frame[sample (1:nrow(frame), N,replace=T), c(channels.ind["CD62L"], channels.ind["CD44"])]
  }

  
  
  
  
  plotDens(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]), main = "CD8 T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(Effector.cd8.flowD@filter)
  lines(Resting.cd8.flowD@filter)
  lines(Naive.cd8.flowD@filter)
 
  
  
#    pop1 <- flowDensity(cd8.Tcells, channels = c(channels.ind["CD62L"], channels.ind["CD44"]),filter=Effector.cd8.flowD@filter) 
#    pop2 <- flowDensity(f, channels = c(channels.ind["CD62L"], channels.ind["CD44"]),filter=pop2@filter) 
#    pop3 <- flowDensity(f, channels = c(channels.ind["CD62L"], channels.ind["CD44"]),filter=pop3@filter) 
   

  #save(f, file=paste0(results.dir,"/CD8TcellsData/",name,".RData"))
 
  return(1)
}



f1score <- function(clr.1,clr.2)
  # Calculates precision, recall and F1 score.
  #
  # Args:
  #  clr of automated method
  # clr of manual method from flowJo
  # Returns:
  #    f1
{
  precision <- sum(clr.1 &clr.2) / sum(clr.1)
  recall <-  sum(clr.1 &clr.2) / sum(clr.2)
  f1 <- 2 * precision * recall / (precision + recall)
  if (f1== "NaN") 
    f1 <- 0
  return(f1)
 
}


###Setting up the code to do the pairwise F1 score

clr.files <- list.files(paste0(results.dir,"/CD8TcellsCLRs400K/"),pattern=".RData")
clr.files <- clr.files[1:90]

pairs <- combn(clr.files,m=2)
t1 <- Sys.time()
#pairs.temp<- pairs[]

F1.all <- apply(pairs, 2, function(p1)
{
  tmp <- load(paste0(results.dir,"/CD8TcellsCLRs400K/",p1[1]))
  clr.1 <- get(tmp)
  tmp <- load(paste0(results.dir,"/CD8TcellsCLRs400K/",p1[2]))
  clr.2 <- get(tmp)
  F1 <- f1score(clr.1,clr.2)
 
}
#

)
t2 <- Sys.time()-t1

t2.2 <- t2
t2.30 <- t2
t2.35 <- t2
t2.50 <- t2
t2.75 <- t2

