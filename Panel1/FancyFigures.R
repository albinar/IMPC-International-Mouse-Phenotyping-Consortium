## Written by Justin Meskas
## Modified by Albina Rahim

setwd("/home/arahim/Desktop/3iTcell/")

library("flowCore")
# library("flowBin")
# library("flowDensity")
# library("flowType")
# library("flowClean")
# library("snowfall")
# library("foreach")

dir.create ( "../3iTcell_Results/Figures/FinalTableOfProportions")

# FinalTableOfProportions <- FinalTable
FinalTableOfProportions <- read.csv("/home/arahim/Desktop/3iTcell_Results/FinalTableOfProportions.csv",header = FALSE)


# PlotNames <- c("UK-Res Prop of Total", "JM-Res Prop of Total", "UK-Res Prop of Parent", "JM-Res Prop of Parent", "UK-Res Cell Count", "JM-Res Cell Count")
# PlotNames <- c("Percentage of Total UK-JM Res", "", "Percentage of Parent UK-JM Res", "", "Cell Count UK-JM Res", "")
PlotNames <- c("Percentage of Events in Each Sub-Population Out of the Total", "", "Percentage of Events in Each Sub-Population Out of the Parent Population", "", "Cell Count of each sub-population", "")

PlotYaxis   <- c("Log10 (Percentage)","","Log10 (Percentage)","","Log10 (Cell Count)")
PlotYaxisNP <- c("Percentage","","Percentage","","Cell Count")


for ( w2 in c(1,3,5) ) {
    dudeTitle <- as.matrix(FinalTableOfProportions[w2,])[2:34]
    dudeTitle <- rbind(paste(dudeTitle, "MA"),paste(dudeTitle,"AA"))
    dudeTitle <- as.vector(dudeTitle)
    dude1 <- as.numeric(as.matrix(FinalTableOfProportions[seq(w2+6,w2+(2658-6), by = 6),2]))
#     dude1[which(is.na(dude1))] <- 0
    dude <- data.frame(log10(dude1))
    NormalPercentagesData <- data.frame(dude1)

        
    dude1b <- as.numeric(as.matrix(FinalTableOfProportions[seq(w2+1+6,w2+1+(2658-6), by = 6),2]))
#     dude1b[which(is.na(dude1b))] <- 0
    dude <- data.frame(dude, log10(dude1b))
    NormalPercentagesData <- data.frame(NormalPercentagesData, dude1b)
    
    for ( w1 in 3:34){

        dude2 <- as.numeric(as.matrix(FinalTableOfProportions[seq(w2+6,w2+(2658-6), by = 6),w1]))
#         dude2[which(is.na(dude2))] <- 0
        dude <- data.frame(dude, log10(dude2))
        NormalPercentagesData <- data.frame(NormalPercentagesData, dude2)
        
        dude3 <- as.numeric(as.matrix(FinalTableOfProportions[seq(w2+1+6,w2+1+(2658-6), by = 6),w1]))
#         dude3[which(is.na(dude3))] <- 0
        dude <- data.frame(dude, log10(dude3))
        NormalPercentagesData <- data.frame(NormalPercentagesData, dude3)
    }
    # Log Percentages
        png ( file = paste("../3iTcell_Results/Figures/FinalTableOfProportions/Box_Plot_", gsub(" ", "_", PlotNames[w2]),".png", sep = "" ), height=1000, width=1600 )
        par(mar=c(35,5,2,2)) # margins
        if (w2 <= 4) {
            boxplot(dude[,-c(1,2)], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab=PlotYaxis[w2], axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
        } else {
            boxplot(dude[,], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab=PlotYaxis[w2], axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
        }
    #     UKoutliers <- sort(c(101, 52, 297, 298, 299, 300, 404, 405, 406, 43, 44, 306, 251, 252, 253, 254, 336))
    #     boxplot(dude[(1:411)[-UKoutliers],], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab="log", axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
    #     axis(1,at=1:33)
        if(w2 <= 4) {
            axis(1,at=1:64,labels=dudeTitle[-c(1,2)],las=2)
        } else {
            axis(1,at=1:66,labels=dudeTitle,las=2)
        }
        axis(2,at=(-17:17)*1)
        dev.off()
    # Regular Percentages
        png ( file = paste("../3iTcell_Results/Figures/FinalTableOfProportions/Box_Plot_", gsub(" ", "_", PlotNames[w2]),"_Normal_Percentage.png", sep = "" ), height=1000, width=1600 )
        par(mar=c(35,5,2,2)) # margins
        if (w2 <= 4) {
            boxplot(NormalPercentagesData[,-c(1,2)], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab=PlotYaxisNP[w2], axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
        } else {
            boxplot(NormalPercentagesData[,], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab=PlotYaxisNP[w2], axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
        }
    #     UKoutliers <- sort(c(101, 52, 297, 298, 299, 300, 404, 405, 406, 43, 44, 306, 251, 252, 253, 254, 336))
    #     boxplot(dude[(1:411)[-UKoutliers],], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab="log", axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
    #     axis(1,at=1:33)
        if(w2 <= 4) {
            axis(1,at=1:64,labels=dudeTitle[-c(1,2)],las=2)
            axis(2,at=seq(0, 100, by= 20))
        } else {
            axis(1,at=1:66,labels=dudeTitle,las=2)
            axis(2,at=seq(0,500000, by=100000))
        }
        dev.off()
}   

    dudeMice <- as.matrix(FinalTableOfProportions[seq(7,1+(2658-6), by = 6),1])

StoreOutliers <- NULL
for ( t1 in seq(1,63, by=2)) {
    walk <- boxplot(dude[,t1], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab="log", axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
    for ( t2 in 1:length(walk$out)){
            StoreOutliers <- c(StoreOutliers, dudeMice[which(walk$out[t2] == dude[,t1])])
    }
}
StoreO <- StoreOutliers
StoreO <- StoreO[duplicated(StoreO)];print(sort(unique(StoreO)))

x <- 1
t.test(dude[,x], dude[,x+1])$p.value; x <- x + 2



