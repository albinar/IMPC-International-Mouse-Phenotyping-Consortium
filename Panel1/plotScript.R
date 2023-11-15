## Date: October 21, 2020
## Created by Albina Rahim
## Script for creating boxplots to highlight the distribution of Proportions between the Wild types and Knockouts

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

library("flowCore")
library("reshape2")
library("stringr")
library("ggplot2")

## This part is for comparing the Proportions distributions between WT versus KO
results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")

CSVfile <- read.csv(paste0(results.dir,"/DCCResults_Proportions_UCD_Panel1_20191105_1819.csv"))
CSVfile <- as.matrix(CSVfile)
tempB <- str_replace_all(CSVfile[,c("FCS.files")],fixed("PANEL_A_"), "")
tempC <- str_replace_all(tempB,fixed(".fcs"), "_PANEL_A.fcs")
CSVfile[,c("FCS.files")] <- tempB
tempB <- str_replace_all(CSVfile[1:739,c("FCS.files")],fixed(".fcs"), "_PANEL_A.fcs")
CSVfile[1:739,c("FCS.files")] <- tempB
CSVfile <- CSVfile[,-4:-6]


CSVfile[,c('Assay.Date')] <- CSVfile[,c('Panel.Organ.Folder')]
temp <- CSVfile[complete.cases(CSVfile), ]


write.csv(temp, file =  paste0(results.dir, "/Proportions_",toupper(centre),"_Panel1.csv"), row.names = FALSE)

Live.matrix <- NA
Live.matrix <- cbind("Live", CSVfile[,c(2,9)])
colnames(Live.matrix) <- c("Populations","Genotype","Proportions")

FSC.matrix <- NA
FSC.matrix <- cbind("FSC.Singlets", CSVfile[,c(2,10)])
colnames(FSC.matrix) <- c("Populations","Genotype","Proportions")


SSC.matrix <- NA
SSC.matrix <- cbind("SSC.Singlets", CSVfile[,c(2,11)])
colnames(SSC.matrix) <- c("Populations","Genotype","Proportions")


NK.matrix <- NA
NK.matrix <- cbind("NK", CSVfile[,c(2,12)])
colnames(NK.matrix) <- c("Populations","Genotype","Proportions")

NK.matrix <- NA
NK.matrix <- cbind("NK", CSVfile[,c(2,12)])
colnames(NK.matrix) <- c("Populations","Genotype","Proportions")

EffectorNK.matrix <- NA
EffectorNK.matrix <- cbind("Effector NK", CSVfile[,c(2,14)])
colnames(EffectorNK.matrix) <- c("Populations","Genotype","Proportions")

RestingNK.matrix <- NA
RestingNK.matrix <- cbind("Resting NK", CSVfile[,c(2,15)])
colnames(RestingNK.matrix) <- c("Populations","Genotype","Proportions")

CD5.matrix <- NA
CD5.matrix <- cbind("CD5+", CSVfile[,c(2,16)])
colnames(CD5.matrix) <- c("Populations","Genotype","Proportions")

CD161CD8.matrix <- NA
CD161CD8.matrix <- cbind("CD161+CD8-", CSVfile[,c(2,17)])
colnames(CD161CD8.matrix) <- c("Populations","Genotype","Proportions")

CD4posNKT.matrix <- NA
CD4posNKT.matrix <- cbind("CD4+NKT", CSVfile[,c(2,19)])
colnames(CD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

CD4negNKT.matrix <- NA
CD4negNKT.matrix <- cbind("CD4-NKT", CSVfile[,c(2,20)])
colnames(CD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD4posNKT.matrix <- NA
EffectorCD4posNKT.matrix <- cbind("Effector CD4+NKT", CSVfile[,c(2,21)])
colnames(EffectorCD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

RestingCD4posNKT.matrix <- NA
RestingCD4posNKT.matrix <- cbind("Resting CD4+NKT", CSVfile[,c(2,22)])
colnames(RestingCD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD4negNKT.matrix <- NA
EffectorCD4negNKT.matrix <- cbind("Effector CD4-NKT", CSVfile[,c(2,23)])
colnames(EffectorCD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

RestingCD4negNKT.matrix <- NA
RestingCD4negNKT.matrix <- cbind("Resting CD4-NKT", CSVfile[,c(2,24)])
colnames(RestingCD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

Tcells.matrix <- NA
Tcells.matrix <- cbind("T cells", CSVfile[,c(2,25)])
colnames(Tcells.matrix) <- c("Populations","Genotype","Proportions")

CD4Tcells.matrix <- NA
CD4Tcells.matrix <- cbind("CD4 T cells", CSVfile[,c(2,26)])
colnames(CD4Tcells.matrix) <- c("Populations","Genotype","Proportions")

CD8Tcells.matrix <- NA
CD8Tcells.matrix <- cbind("CD8 T cells", CSVfile[,c(2,27)])
colnames(CD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

Tregs.matrix <- NA
Tregs.matrix <- cbind("Tregs", CSVfile[,c(2,28)])
colnames(Tregs.matrix) <- c("Populations","Genotype","Proportions")

Thelper.matrix <- NA
Thelper.matrix <- cbind("Thelper", CSVfile[,c(2,29)])
colnames(Thelper.matrix) <- c("Populations","Genotype","Proportions")

EffectorTregs.matrix <- NA
EffectorTregs.matrix <- cbind("Effector Tregs", CSVfile[,c(2,30)])
colnames(EffectorTregs.matrix) <- c("Populations","Genotype","Proportions")

RestingTregs.matrix <- NA
RestingTregs.matrix <- cbind("Resting Tregs", CSVfile[,c(2,31)])
colnames(RestingTregs.matrix) <- c("Populations","Genotype","Proportions")

EffectorThelper.matrix <- NA
EffectorThelper.matrix <- cbind("Effector Thelper", CSVfile[,c(2,32)])
colnames(EffectorThelper.matrix) <- c("Populations","Genotype","Proportions")

RestingThelper.matrix <- NA
RestingThelper.matrix <- cbind("Resting Thelper", CSVfile[,c(2,33)])
colnames(RestingThelper.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD8Tcells.matrix <- NA
EffectorCD8Tcells.matrix <- cbind("Effector CD8 Tcells", CSVfile[,c(2,34)])
colnames(EffectorCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

RestingCD8Tcells.matrix <- NA
RestingCD8Tcells.matrix <- cbind("Resting CD8 Tcells", CSVfile[,c(2,35)])
colnames(RestingCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

NaiveCD8Tcells.matrix <- NA
NaiveCD8Tcells.matrix <- cbind("Naive CD8 Tcells", CSVfile[,c(2,36)])
colnames(NaiveCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

main.matrix <- rbind(Live.matrix, FSC.matrix, SSC.matrix, NK.matrix, EffectorNK.matrix, RestingNK.matrix, CD5.matrix,
                     CD161CD8.matrix, CD4posNKT.matrix, CD4negNKT.matrix, EffectorCD4posNKT.matrix, RestingCD4posNKT.matrix,
                     EffectorCD4negNKT.matrix, RestingCD4negNKT.matrix, Tcells.matrix, CD4Tcells.matrix, CD8Tcells.matrix,
                     Tregs.matrix, Thelper.matrix, EffectorThelper.matrix, RestingThelper.matrix, EffectorCD8Tcells.matrix,
                     RestingCD8Tcells.matrix, NaiveCD8Tcells.matrix)


p <- ggplot(data = main.matrix, aes(x=Populations, y=Proportions)) + geom_boxplot(aes(fill=Genotype))
p + facet_wrap(~ Populations, scales = "free")

################################################################
## Plotting bosplots of Proportions calculated based on the Live cell counts

EventsFile <-  read.csv(paste0(results.dir,"/UCD_Panel1_EventCounts2021.csv"))
#EventsFile <- as.matrix(EventsFile)
tempB <- str_replace_all(EventsFile[,c("FCS.files")],fixed("PANEL_A_"), "")
tempC <- str_replace_all(tempB,fixed(".fcs"), "_PANEL_A.fcs")
EventsFile[,c("FCS.files")] <- tempB
tempB <- str_replace_all(EventsFile[1:725,c("FCS.files")],fixed(".fcs"), "_PANEL_A.fcs")
EventsFile[1:725,c("FCS.files")] <- tempB
EventsFile <- EventsFile[,-4:-6]


EventsFile[,c('Assay.Date')] <- EventsFile[,c('Panel.Organ.Folder')]
#EventsFile[,c('Panel.Organ.Folder')] <- EventsFile[,c('Assay.Date')] 
temp <- EventsFile[complete.cases(EventsFile), ]


write.csv(temp, file =  paste0(results.dir, "/Event_Counts_",toupper(centre),"_Panel1.csv"), row.names = FALSE)

EventsFile<-temp

## Creating a dataframe with Proportions based on Live counts

#Proportions.Live.matrix <-  matrix(nrow = nrow(EventsFile), ncol = ncol(EventsFile), data = NA)# matrix for saving the event counts
Proportions.Live.matrix <-  EventsFile


for(i in 1:nrow(EventsFile)){
  Proportions.Live.matrix[i,11] <- (as.numeric(EventsFile[i,11])/as.numeric(EventsFile[i,10]))*100 ## Live
  Proportions.Live.matrix[i,12] <- (as.numeric(EventsFile[i,12])/as.numeric(EventsFile[i,11]))*100 ## FCS.Singlets
  Proportions.Live.matrix[i,13] <- (as.numeric(EventsFile[i,13])/as.numeric(EventsFile[i,11]))*100 ## SCS.Singlets
  Proportions.Live.matrix[i,14] <- (as.numeric(EventsFile[i,14])/as.numeric(EventsFile[i,11]))*100 ## NK cells
  Proportions.Live.matrix[i,15] <- (as.numeric(EventsFile[i,15])/as.numeric(EventsFile[i,11]))*100 ## NOT NK cells
  Proportions.Live.matrix[i,16] <- (as.numeric(EventsFile[i,16])/as.numeric(EventsFile[i,11]))*100 ## Effector NK
  Proportions.Live.matrix[i,17] <- (as.numeric(EventsFile[i,17])/as.numeric(EventsFile[i,11]))*100 ## Resting NK
  Proportions.Live.matrix[i,18] <- (as.numeric(EventsFile[i,18])/as.numeric(EventsFile[i,11]))*100 ## CD5+
  Proportions.Live.matrix[i,19] <- (as.numeric(EventsFile[i,19])/as.numeric(EventsFile[i,11]))*100 ## CD161_CD8-
  Proportions.Live.matrix[i,20] <- (as.numeric(EventsFile[i,20])/as.numeric(EventsFile[i,11]))*100 ## NOT CD161_CD8-
  Proportions.Live.matrix[i,21] <- (as.numeric(EventsFile[i,21])/as.numeric(EventsFile[i,11]))*100 ## CD4+ NKT
  Proportions.Live.matrix[i,22] <- (as.numeric(EventsFile[i,22])/as.numeric(EventsFile[i,11]))*100 ## CD4- NKT
  Proportions.Live.matrix[i,23] <- (as.numeric(EventsFile[i,23])/as.numeric(EventsFile[i,11]))*100 ## Effector CD4+ NKT
  Proportions.Live.matrix[i,24] <- (as.numeric(EventsFile[i,24])/as.numeric(EventsFile[i,11]))*100 ## Resting CD4+ NKT
  Proportions.Live.matrix[i,25] <- (as.numeric(EventsFile[i,25])/as.numeric(EventsFile[i,11]))*100 ## Effector CD4- NKT
  Proportions.Live.matrix[i,26] <- (as.numeric(EventsFile[i,26])/as.numeric(EventsFile[i,11]))*100 ## Resting CD4- NKT
  Proportions.Live.matrix[i,27] <- (as.numeric(EventsFile[i,27])/as.numeric(EventsFile[i,11]))*100 ## T cells
  Proportions.Live.matrix[i,28] <- (as.numeric(EventsFile[i,28])/as.numeric(EventsFile[i,11]))*100 ## CD4 Tcells
  Proportions.Live.matrix[i,29] <- (as.numeric(EventsFile[i,29])/as.numeric(EventsFile[i,11]))*100 ## CD8 Tcells
  Proportions.Live.matrix[i,30] <- (as.numeric(EventsFile[i,30])/as.numeric(EventsFile[i,11]))*100 ## Tregs
  Proportions.Live.matrix[i,31] <- (as.numeric(EventsFile[i,31])/as.numeric(EventsFile[i,11]))*100 ## Thelper
  Proportions.Live.matrix[i,32] <- (as.numeric(EventsFile[i,32])/as.numeric(EventsFile[i,11]))*100 ## Effector Treg
  Proportions.Live.matrix[i,33] <- (as.numeric(EventsFile[i,33])/as.numeric(EventsFile[i,11]))*100 ## Resting Treg
  Proportions.Live.matrix[i,34] <- (as.numeric(EventsFile[i,34])/as.numeric(EventsFile[i,11]))*100 ## Effector Thelper
  Proportions.Live.matrix[i,35] <- (as.numeric(EventsFile[i,35])/as.numeric(EventsFile[i,11]))*100 ## Resting Thelper
  Proportions.Live.matrix[i,36] <- (as.numeric(EventsFile[i,36])/as.numeric(EventsFile[i,11]))*100 ## Effector CD8 Tcells
  Proportions.Live.matrix[i,37] <- (as.numeric(EventsFile[i,37])/as.numeric(EventsFile[i,11]))*100 ## Resting CD8 Tcells
  Proportions.Live.matrix[i,38] <- (as.numeric(EventsFile[i,38])/as.numeric(EventsFile[i,11]))*100 ## Naive CD8 Tcells
  Proportions.Live.matrix[i,39] <- (as.numeric(EventsFile[i,39])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ NK
  Proportions.Live.matrix[i,40] <- (as.numeric(EventsFile[i,40])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ CD4+ NKT
  Proportions.Live.matrix[i,41] <- (as.numeric(EventsFile[i,41])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ CD4- NKT
  Proportions.Live.matrix[i,42] <- (as.numeric(EventsFile[i,42])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ Tregs
  Proportions.Live.matrix[i,43] <- (as.numeric(EventsFile[i,43])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ T helper
  Proportions.Live.matrix[i,44] <- (as.numeric(EventsFile[i,44])/as.numeric(EventsFile[i,11]))*100 ## KLRG1+ CD8 Tcells
  
}

## Boxplot
Live.matrix <- NA
Live.matrix <- cbind("Live", Proportions.Live.matrix[,c(2,11)])
colnames(Live.matrix) <- c("Populations","Genotype","Proportions")

FSC.matrix <- NA
FSC.matrix <- cbind("FSC.Singlets", Proportions.Live.matrix[,c(2,12)])
colnames(FSC.matrix) <- c("Populations","Genotype","Proportions")


SSC.matrix <- NA
SSC.matrix <- cbind("SSC.Singlets", Proportions.Live.matrix[,c(2,13)])
colnames(SSC.matrix) <- c("Populations","Genotype","Proportions")


NK.matrix <- NA
NK.matrix <- cbind("NK", Proportions.Live.matrix[,c(2,14)])
colnames(NK.matrix) <- c("Populations","Genotype","Proportions")

Not.NK.matrix <- NA
Not.NK.matrix <- cbind("Not NK", Proportions.Live.matrix[,c(2,15)])
colnames(Not.NK.matrix) <- c("Populations","Genotype","Proportions")

EffectorNK.matrix <- NA
EffectorNK.matrix <- cbind("Effector NK", Proportions.Live.matrix[,c(2,16)])
colnames(EffectorNK.matrix) <- c("Populations","Genotype","Proportions")

RestingNK.matrix <- NA
RestingNK.matrix <- cbind("Resting NK", Proportions.Live.matrix[,c(2,17)])
colnames(RestingNK.matrix) <- c("Populations","Genotype","Proportions")

CD5.matrix <- NA
CD5.matrix <- cbind("CD5+", Proportions.Live.matrix[,c(2,18)])
colnames(CD5.matrix) <- c("Populations","Genotype","Proportions")

CD161CD8.matrix <- NA
CD161CD8.matrix <- cbind("CD161+CD8-", Proportions.Live.matrix[,c(2,19)])
colnames(CD161CD8.matrix) <- c("Populations","Genotype","Proportions")

CD4posNKT.matrix <- NA
CD4posNKT.matrix <- cbind("CD4+NKT", Proportions.Live.matrix[,c(2,21)])
colnames(CD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

CD4negNKT.matrix <- NA
CD4negNKT.matrix <- cbind("CD4-NKT", Proportions.Live.matrix[,c(2,22)])
colnames(CD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD4posNKT.matrix <- NA
EffectorCD4posNKT.matrix <- cbind("Effector CD4+NKT", Proportions.Live.matrix[,c(2,23)])
colnames(EffectorCD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

RestingCD4posNKT.matrix <- NA
RestingCD4posNKT.matrix <- cbind("Resting CD4+NKT", Proportions.Live.matrix[,c(2,24)])
colnames(RestingCD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD4negNKT.matrix <- NA
EffectorCD4negNKT.matrix <- cbind("Effector CD4-NKT", Proportions.Live.matrix[,c(2,25)])
colnames(EffectorCD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

RestingCD4negNKT.matrix <- NA
RestingCD4negNKT.matrix <- cbind("Resting CD4-NKT", Proportions.Live.matrix[,c(2,26)])
colnames(RestingCD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

Tcells.matrix <- NA
Tcells.matrix <- cbind("T cells", Proportions.Live.matrix[,c(2,27)])
colnames(Tcells.matrix) <- c("Populations","Genotype","Proportions")

CD4Tcells.matrix <- NA
CD4Tcells.matrix <- cbind("CD4 T cells", Proportions.Live.matrix[,c(2,28)])
colnames(CD4Tcells.matrix) <- c("Populations","Genotype","Proportions")

CD8Tcells.matrix <- NA
CD8Tcells.matrix <- cbind("CD8 T cells", Proportions.Live.matrix[,c(2,29)])
colnames(CD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

Tregs.matrix <- NA
Tregs.matrix <- cbind("Tregs", Proportions.Live.matrix[,c(2,30)])
colnames(Tregs.matrix) <- c("Populations","Genotype","Proportions")

Thelper.matrix <- NA
Thelper.matrix <- cbind("Thelper", Proportions.Live.matrix[,c(2,31)])
colnames(Thelper.matrix) <- c("Populations","Genotype","Proportions")

EffectorTregs.matrix <- NA
EffectorTregs.matrix <- cbind("Effector Tregs", Proportions.Live.matrix[,c(2,32)])
colnames(EffectorTregs.matrix) <- c("Populations","Genotype","Proportions")

RestingTregs.matrix <- NA
RestingTregs.matrix <- cbind("Resting Tregs", Proportions.Live.matrix[,c(2,33)])
colnames(RestingTregs.matrix) <- c("Populations","Genotype","Proportions")

EffectorThelper.matrix <- NA
EffectorThelper.matrix <- cbind("Effector Thelper", Proportions.Live.matrix[,c(2,34)])
colnames(EffectorThelper.matrix) <- c("Populations","Genotype","Proportions")

RestingThelper.matrix <- NA
RestingThelper.matrix <- cbind("Resting Thelper", Proportions.Live.matrix[,c(2,35)])
colnames(RestingThelper.matrix) <- c("Populations","Genotype","Proportions")

EffectorCD8Tcells.matrix <- NA
EffectorCD8Tcells.matrix <- cbind("Effector CD8 Tcells", Proportions.Live.matrix[,c(2,36)])
colnames(EffectorCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

RestingCD8Tcells.matrix <- NA
RestingCD8Tcells.matrix <- cbind("Resting CD8 Tcells", Proportions.Live.matrix[,c(2,37)])
colnames(RestingCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

NaiveCD8Tcells.matrix <- NA
NaiveCD8Tcells.matrix <- cbind("Naive CD8 Tcells", Proportions.Live.matrix[,c(2,38)])
colnames(NaiveCD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

klrg1NK.matrix <- NA
klrg1NK.matrix <- cbind("KLRG1+ NK", Proportions.Live.matrix[,c(2,39)])
colnames(klrg1NK.matrix) <- c("Populations","Genotype","Proportions")

klrg1CD4posNKT.matrix <- NA
klrg1CD4posNKT.matrix <- cbind("KLRG1+ CD4+ NKT", Proportions.Live.matrix[,c(2,40)])
colnames(klrg1CD4posNKT.matrix) <- c("Populations","Genotype","Proportions")

klrg1CD4negNKT.matrix <- NA
klrg1CD4negNKT.matrix <- cbind("KLRG1+ CD4- NKT", Proportions.Live.matrix[,c(2,41)])
colnames(klrg1CD4negNKT.matrix) <- c("Populations","Genotype","Proportions")

klrg1Treg.matrix <- NA
klrg1Treg.matrix <- cbind("KLRG1+ Treg", Proportions.Live.matrix[,c(2,42)])
colnames(klrg1Treg.matrix) <- c("Populations","Genotype","Proportions")

klrg1Thelper.matrix <- NA
klrg1Thelper.matrix <- cbind("KLRG1+ T helper", Proportions.Live.matrix[,c(2,43)])
colnames(klrg1Thelper.matrix) <- c("Populations","Genotype","Proportions")

klrg1CD8Tcells.matrix <- NA
klrg1CD8Tcells.matrix <- cbind("KLRG1+ CD8 Tcells", Proportions.Live.matrix[,c(2,44)])
colnames(klrg1CD8Tcells.matrix) <- c("Populations","Genotype","Proportions")

main.matrix <- rbind(Live.matrix, FSC.matrix, SSC.matrix, NK.matrix, Not.NK.matrix, EffectorNK.matrix, RestingNK.matrix, CD5.matrix,
                     CD161CD8.matrix, CD4posNKT.matrix, CD4negNKT.matrix, EffectorCD4posNKT.matrix, RestingCD4posNKT.matrix,
                     EffectorCD4negNKT.matrix, RestingCD4negNKT.matrix, Tcells.matrix, CD4Tcells.matrix, CD8Tcells.matrix,
                     Tregs.matrix, Thelper.matrix, EffectorThelper.matrix, RestingThelper.matrix, EffectorCD8Tcells.matrix,
                     RestingCD8Tcells.matrix, NaiveCD8Tcells.matrix, klrg1NK.matrix, klrg1CD4posNKT.matrix, klrg1CD4negNKT.matrix,
                     klrg1Treg.matrix, klrg1Thelper.matrix, klrg1CD8Tcells.matrix)

main.matrix <- rbind(Live.matrix, FSC.matrix, SSC.matrix, NK.matrix, Not.NK.matrix, EffectorNK.matrix, RestingNK.matrix, CD5.matrix,
                     CD161CD8.matrix, CD4posNKT.matrix, CD4negNKT.matrix, EffectorCD4posNKT.matrix, RestingCD4posNKT.matrix,
                     EffectorCD4negNKT.matrix, RestingCD4negNKT.matrix, Tcells.matrix, CD4Tcells.matrix, CD8Tcells.matrix,
                     Tregs.matrix, Thelper.matrix, EffectorThelper.matrix, RestingThelper.matrix, EffectorCD8Tcells.matrix,
                     RestingCD8Tcells.matrix, NaiveCD8Tcells.matrix)

main.matrix <- Live.matrix


p <- ggplot(data = main.matrix, aes(x=Populations, y=Proportions)) + geom_boxplot(aes(fill=Genotype))
p + facet_wrap(~ Populations, scales = "free")
 
#############################################################################################


## This part is for comparing manual versus automated
results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results-2019")
CSVfile <- read.csv(paste0(results.dir,"/Proportions_UCD_Panel1.csv"))
CSVfile <- as.matrix(CSVfile)

#manualKOfile <- read.csv(paste0(results.dir, "/KOMP Complete Data.csv"))
manualWTfile1 <- read.csv(paste0(results.dir, "/WT.csv"))
manualWTfile1 <- as.matrix(manualWTfile1)

tempA <- str_replace_all(CSVfile[,c("FCS.files")],fixed("_PANEL_A.fcs"), "")

tempM<- toupper(manualWTfile1[,1])
tempM <- str_replace_all(tempM,fixed("_PANEL_A.FCS"), "")
tempM <- str_replace_all(tempM,fixed("_PANEL.FCS"), "")
tempM <- str_replace_all(tempM,fixed(".FCS"), "")
tempM <- str_replace_all(tempM,fixed("F.FCS"), "")
tempM <- str_replace_all(tempM,fixed("M.FCS"), "")
tempM <- str_replace_all(tempM,fixed("WT_"), "")
tempM <- str_replace_all(tempM,fixed("WTS_"), "")
tempM <- str_replace_all(tempM,fixed("WTS "), "")
tempM <- str_replace_all(tempM,fixed("WT`S_"), "")
tempM <- str_replace_all(tempM,fixed("WT TEST_"), "")
tempM <- str_replace_all(tempM,fixed("-"), "_")

matchIndex <- 0
for(i in 1:length(tempA)){
   x<- which(tempA[i] == tempM)
   matchIndex <- c(matchIndex,x)
}
matchIndex <- matchIndex[matchIndex!=0]

# manualCSVfile <- matrix(nrow = 108, ncol = 12, data = NA)# matrix for saving the event counts
# for (i in 1:length(matchIndex)){
#   manualCSVfile[i,] <- manualKOfile[matchIndex[i], c(1,17,19,21,22,24,25,26,28,29,30,39)]
# }


manualWTfile <- matrix(nrow = 1, ncol = 12, data = NA)# matrix for saving the event counts
for (i in 1:length(matchIndex)){
  manualWTfile[i,] <- manualKOfile[matchIndex[i], c(1,17,19,21,22,24,25,26,28,29,30,39)]
}



manualWTfile2 <- read.csv(paste0(results.dir, "/Late-WT.csv"))
manualWTfile2 <- as.matrix(manualWTfile2)

tempM<- toupper(manualWTfile2[,1])
tempM <- str_replace_all(tempM,fixed("_PANEL_A.FCS"), "")
tempM <- str_replace_all(tempM,fixed("_PANEL.FCS"), "")
tempM <- str_replace_all(tempM,fixed(".FCS"), "")
tempM <- str_replace_all(tempM,fixed("F.FCS"), "")
tempM <- str_replace_all(tempM,fixed("M.FCS"), "")
tempM <- str_replace_all(tempM,fixed("WT_"), "")
tempM <- str_replace_all(tempM,fixed("WTS_"), "")
tempM <- str_replace_all(tempM,fixed("WTS "), "")
tempM <- str_replace_all(tempM,fixed("WT`S_"), "")
tempM <- str_replace_all(tempM,fixed("WT TEST_"), "")
tempM <- str_replace_all(tempM,fixed("LWT_"), "")
tempM <- str_replace_all(tempM,fixed("LWTS_"), "")
tempM <- str_replace_all(tempM,fixed("LC"), "C")
tempM <- str_replace_all(tempM,fixed("-"), "_")

matchIndex <- 0
for(i in 1:length(tempA)){
  x<- which(tempA[i] == tempM)
  matchIndex <- c(matchIndex,x)
}
matchIndex <- matchIndex[matchIndex!=0]

# manualCSVfile <- matrix(nrow = 108, ncol = 12, data = NA)# matrix for saving the event counts
# for (i in 1:length(matchIndex)){
#   manualCSVfile[i,] <- manualKOfile[matchIndex[i], c(1,17,19,21,22,24,25,26,28,29,30,39)]
# }


manualWTfile <- matrix(nrow = 1, ncol = 12, data = NA)# matrix for saving the event counts
for (i in 1:length(matchIndex)){
  manualWTfile[i,] <- manualKOfile[matchIndex[i], c(1,17,19,21,22,24,25,26,28,29,30,39)]
}
