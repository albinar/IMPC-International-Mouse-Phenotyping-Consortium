##################################################################################################
# removes groups of cells that show some statisitcally significant differences from other groups #
##################################################################################################
flowCut <- function(f, segment=500, CleanChan.loc, directory=NULL, FileID=NULL, Plot="None", MaxContin=0.2, MeanOfMeans=0.33, MaxOfMeans=0.5, verbose=F){

    start0 <- Sys.time()

    resTable <- matrix("",15,1)
    rownames(resTable) <- c("Is it monotonically increasing in time",
                         "Largest continuous jump", "Continuous - Pass",
                         "Mean of % of range of means divided by range of data", "Mean of % - Pass",
                         "Max of % of range of means divided by range of data", "Max of % - Pass",
                         "Has a low density section been removed", "% of low density removed",
                         "How many segments have been removed", "% of events removed from segments removed",
                         "Worst Channel", "% of events removed", "FileID", "Has the file passed")

    if (is.null(directory)){ directory <- paste0(getwd(), "/flowCut") }
    if( is.null(FileID) ) { # if no ID is given, generate a random one
        samp <- sample(1:9999, 1)
        FileID <- paste0(rep(0,4-nchar(samp)), samp);
        if (verbose == T){ cat(paste0("The FileID is: ", FileID, "\n"))}
    }
    resTable["FileID",] <- FileID

    if (nrow(f) <= 2*segment){
        message("Either your segment size is too large or your number of cells is too small.")
        to.be.removed <- NULL
        worstChan <- NULL
        return(list(frame=f, ind=to.be.removed, data=resTable, worstChan=worstChan))
    }
    t.name <- f@parameters@data$name
    FSC.loc  <- sort(unique(c(grep("fsc", tolower(t.name)), grep("fs lin", tolower(t.name)))));  names( FSC.loc) <- NULL
    SSC.loc  <- sort(unique(c(grep("ssc", tolower(t.name)), grep("ss lin", tolower(t.name)))));  names( SSC.loc) <- NULL
    # SSC.loc  <- grep("ssc", tolower(f@parameters@data$name)); names( SSC.loc) <- NULL
    Time.loc <- which(tolower(f@parameters@data$name) == "time"); names(Time.loc) <- NULL
    # CleanChan.loc <- (1:ncol(f))[-c(FSC.loc, SSC.loc, Time.loc)]

    if (length(Time.loc) == 0){
        message(paste0("Your data does not have a time channel. flowCut will create one, but now flowCut will not be as fully",
                        " functional as it could be. Consider recording the time for future projects."))
        f@exprs <- cbind(f@exprs, 1:nrow(f))
        colnames(f@exprs)[length(colnames(f))+1] <- "Time"
        # colnames(f)[length(colnames(f))] <- "Time"
        # colnames(f) <- c(colnames(f), "Time")
        f@parameters@data <- rbind(f@parameters@data, c("Time", "Time", 262144, -111, 262143))
        rownames(f@parameters@data)[length(colnames(f))] <- paste0("$P",length(colnames(f)))
        f@description[paste0("P", length(colnames(f))+1, "DISPLAY")] <- "LOG"
        f@description[paste0("P", length(colnames(f))+1, "N")] <- "Time"
        Time.loc <- length(colnames(f))
    }

    ind.removed <- NA
    f.org <- f

    #### Test Monochromatic ###############################################
    if ( all(f@exprs[,Time.loc] == cummax(f@exprs[,Time.loc])) == F ){
        message("The flow frame is not monotonically increasing in time")
        resTable["Is it monotonically increasing in time",] <- "F"
    } else {
        resTable["Is it monotonically increasing in time",] <- "T"
    }

    #### Remove Low Density Sections ######################################
    res.temp <- removeLowDensSections(f)
    f <- res.temp$frame; removeIndLowDens <- res.temp$rem.ind; remove(res.temp)
    ifelse(length(removeIndLowDens) >= 1 , resTable["Has a low density section been removed",] <- "T",
                                           resTable["Has a low density section been removed",] <- "F")
    resTable["% of low density removed",] <- as.character( round(length(removeIndLowDens) / nrow(f.org), digits=4) * 100 )

    #### Calculate means and which bins will be removed ##################
    res.temp <- calcMeansAndBinsRemoved(f, segment, CleanChan.loc, FirstOrSecond="First")
    deletedSegments1 <- res.temp$deletedSegments; quantiles1   <- res.temp$quantiles;   storeMeans1 <- res.temp$storeMeans;
    meanRangePerc1   <- res.temp$meanRangePerc;   timeCentres1 <- res.temp$timeCentres; remove(res.temp)
    #### delete segments that are statistically different ###############
    removed.ind <- NULL
    totalNumSeg <- floor(nrow(f@exprs)/segment)
    if (is.null(deletedSegments1)){ # if nothing is deleted
        if (verbose == T){ cat("None deleted from flowCut segment removal.\n")}
        resTable["How many segments have been removed",] <- as.character(0)
    } else {
        deletedSegments1 <- sort( unique( deletedSegments1 ) , decreasing = T )
        if (verbose == T){cat(paste0("Removing segments ", paste0(sort(deletedSegments1), collapse = ", "), " out of ", totalNumSeg, " sections.\n"))}
        resTable["How many segments have been removed",] <- as.character( length(deletedSegments1) )
        for (n in 1:length(deletedSegments1)){
            if (deletedSegments1[n] == totalNumSeg)
                removed.ind <- c(removed.ind, (segment*(deletedSegments1[n]-1)+1):nrow(f@exprs))
            if (deletedSegments1[n] != totalNumSeg)
                removed.ind <- c(removed.ind, segment*(deletedSegments1[n]-1)+(1:segment))
        }
        f@exprs <- f@exprs[-removed.ind,]
    }
    resTable["% of events removed from segments removed",] <- as.character( round(length(removed.ind) / nrow(f.org), digits=4) * 100 )

    #### 2nd Time: Calculate means and which bins will be removed ########
    res.temp <- calcMeansAndBinsRemoved(f, segment, CleanChan.loc, FirstOrSecond="Second")
    # deletedSegments2 <- res.temp$deletedSegments;
    quantiles2 <- res.temp$quantiles; storeMeans2 <- res.temp$storeMeans;
    meanRangePerc2 <- res.temp$meanRangePerc; timeCentres2 <- res.temp$timeCentres; remove(res.temp)

    #### Continuous Check ################################################
    distJumped <- list()
    maxDistJumped <- rep(0, max(CleanChan.loc))
    for ( j in CleanChan.loc){
        temp.vect <- rep(0,length(storeMeans2[[j]])-1)
        for(l in 1:(length(storeMeans2[[j]])-1) ){
            temp.vect[l] <- abs(storeMeans2[[j]][l]-storeMeans2[[j]][[l+1]]) / (quantiles1[[j]]["98%"]-quantiles1[[j]]["2%"])
        }
        # distJumped[[j]] <- temp.vect
        maxDistJumped[j] <- max(temp.vect)
    }
    resTable["Largest continuous jump",] <- as.character( round(max(maxDistJumped, na.rm = T), digits=3) )
    ifelse ( resTable["Largest continuous jump",] >= as.numeric(MaxContin), resTable["Continuous - Pass",] <- "F", resTable["Continuous - Pass",] <- "T")

    #### Mean of means ###############################################
    resTable["Mean of % of range of means divided by range of data",] <- as.character( round(mean(meanRangePerc2, na.rm = T), digits=3) )
    if (resTable["Mean of % of range of means divided by range of data",] >= as.numeric(MeanOfMeans)) {
        if (verbose == T){ message(paste0("File Failed. The means differ more than ", MeanOfMeans, "% of the range of the 2-98 percentile of the full data."))}
        resTable["Mean of % - Pass",] <- "F"
    } else {
        resTable["Mean of % - Pass",] <- "T"
    }

    #### Max of means ################################################
    resTable["Max of % of range of means divided by range of data",] <- round(max(meanRangePerc2, na.rm = T), digits=3)
    worstChan <- which(meanRangePerc1 == max(meanRangePerc1, na.rm = T)) # use 1 because we want the worst marker before any corrections.
    resTable["Worst Channel",] <- worstChan
    if ( resTable["Max of % of range of means divided by range of data",] >= MaxOfMeans) {
        if (verbose == T){ message("File Failed. The max ranged means differ more than ", MaxOfMeans, "% of the range of the 2-98 percentile of the full data.")}
        resTable["Max of % - Pass",] <- "F"
    } else {
        resTable["Max of % - Pass",] <- "T"
    }

    #### organize the indices that have been removed #################
    if(  is.null(removed.ind ) &&  is.null(removeIndLowDens))
        to.be.removed <- NULL
    if(  is.null(removed.ind ) && !is.null(removeIndLowDens))
        to.be.removed <- removeIndLowDens
    if( !is.null(removed.ind ) &&  is.null(removeIndLowDens))
        to.be.removed <- removed.ind
    if( !is.null(removed.ind ) && !is.null(removeIndLowDens)){
        # to.be.removed <- unique(removed.ind, (1:nrow(f.org))[-removed.ind][removeIndLowDens])
        temp <- setdiff(1:nrow(f.org), removeIndLowDens) # lowDens was removed first
        to.be.kept <- temp[setdiff(1:length(temp), removed.ind)]
        to.be.removed <- setdiff(1:nrow(f.org), to.be.kept)
    }

    resTable["% of events removed",] <- as.character( round(length(to.be.removed) / nrow(f.org), digits=4) * 100 )
    if( "TTTT" == paste0(resTable["Is it monotonically increasing in time",], resTable["Continuous - Pass",],
                         resTable["Mean of % - Pass",],                       resTable["Max of % - Pass",]) ){
        resTable["Has the file passed",] <- "T"
    } else {
        resTable["Has the file passed",] <- "F"
    }

    #################################################################################
    # save the time plots with black points indicating which events will be deleted #
    #################################################################################

    if(resTable["Is it monotonically increasing in time",]   == "T") {PassedMono  <- "T"     } else { PassedMono  <- "F"}
    if(resTable["Continuous - Pass",]   == "T") {PassedCont  <- "T"     } else { PassedCont  <- "F"}
    if(resTable["Mean of % - Pass",]    == "T") {PassedMean  <- "T"     } else { PassedMean  <- "F"}
    if(resTable["Max of % - Pass",]     == "T") {PassedMax   <- "T"     } else { PassedMax   <- "F"}
    if(resTable["Has the file passed",] == "T") {FailedOrNot <- "Passed"} else { FailedOrNot <- "Failed"}

    if (Plot == "All" || (Plot == "Failed Only" && FailedOrNot == "Failed") ){
        z <- ceiling((length(CleanChan.loc)+1)/4)
        suppressWarnings ( dir.create ( paste0(directory), recursive = T) )

        CairoPNG ( file = paste0(directory, "/", FileID, "_", segment, "_", FailedOrNot, "_",
                                PassedMono, PassedCont, PassedMean, PassedMax,".png"), width = 4*400, height = z*400)
        par(mfrow=c(z,4),mar=c(5,5,4,2))
            for (x in CleanChan.loc){
                plotDens(f.org, c(Time.loc,x), cex.main=2, cex.lab=2, cex.axis=2,
                         main=paste0(round(meanRangePerc1[x],digits=3), " / ", round(meanRangePerc2[x],digits=3), " (",
                                     round(maxDistJumped[x], digits=3), ")" ))

                if ( length(to.be.removed) != 0 )
                    points(exprs(f.org)[to.be.removed,c(Time.loc,x)],col=1,pch=".")
                if ( (length(removeIndLowDens) != 0 ) )
                    points(exprs(f.org)[removeIndLowDens,c(Time.loc,x)],pch=".", cex=1, col="grey")
                lines(x=(timeCentres1),y = storeMeans1[[x]], cex.main=2, cex.lab=2, cex.axis=2,  lwd=4, col="deeppink2")
                lines(x=(timeCentres2),y = storeMeans2[[x]], cex.main=2, cex.lab=2, cex.axis=2,  lwd=4, col="brown")
                abline(h=c(quantiles1[[x]]["98%"], quantiles1[[x]]["2%"]), lwd=4, col="chocolate2")
                abline(h=c(quantiles2[[x]]["98%"], quantiles2[[x]]["2%"]), lwd=4, col="chocolate4")
            }
            x <- as.numeric(resTable["Worst Channel",])
            plotDens(f.org, c(Time.loc,x), cex.main=2, cex.lab=2, cex.axis=2, main=paste0("Worst Channel without indices removed"))
        dev.off ( )
    }

if (verbose == T){ if(resTable["Has the file passed",] == "T") {cat("File Passed\n")} else { cat(paste0("File Failed ",
                                PassedMono, PassedCont, PassedMean, PassedMax, "\n"))}}
if (verbose == T){ cat("Cleaning completed in: ",TimeOutput(start0),"\n",sep="")}
return(list(frame=f, ind=to.be.removed, data=resTable, worstChan=worstChan))
}









################################################
# remove sections that have a very low density #
################################################
removeLowDensSections <- function(f){

    Time.loc <- which(tolower(f@parameters@data$name) == "time"); names(Time.loc) <- NULL

    minTime <- min(f@exprs[ ,Time.loc])
    maxTime <- max(f@exprs[ ,Time.loc])
    # Add extra points at the beginning and the end to deal with the fact that density curves always have a zero point at the ends
    # Then when the ranges to be removed contain the left most or right most points they can be removed because they are fake anyways.
    # This method allows the algorithm to remove events that have low density at the endpoints.
    # There are some annoying but required steps here. Mostly because the density functions indices do not match with the flowFrames indices.
    spread <- (max(f@exprs[ ,Time.loc])-min(f@exprs[ ,Time.loc])) / 20
    amount <- nrow(f)*(1000/300000) # most files tested on were 300000 and 1000 was the number I used to use. Useful if very small amount of events
    dens.f <- density(c(rep(minTime-spread, amount), f@exprs[ ,Time.loc], rep(maxTime+spread, amount)), n= nrow(f), adjust=0.1)
    low.dens.removeIndLowDens <- which(dens.f$y <= 0.1*max(dens.f$y))
# plot(dens.f)
    # find all the consecutive sections that have low density
    range.low.dens <- list()
    count <- 1
    low <- 1

    if( length(low.dens.removeIndLowDens) == 1){
        range.low.dens[[1]] <- low.dens.removeIndLowDens
    } else {
        for(b1 in 1:(length(low.dens.removeIndLowDens)-1) ) {
            if( low.dens.removeIndLowDens[b1] != (low.dens.removeIndLowDens[b1+1]-1) ){
                range.low.dens[[count]] <- low.dens.removeIndLowDens[low]:low.dens.removeIndLowDens[b1]
                low <- b1+1
                count <- count + 1
            }
            if( b1 == (length(low.dens.removeIndLowDens)-1) ) {# end point
                if(low.dens.removeIndLowDens[b1] == (low.dens.removeIndLowDens[b1+1]-1)){
                    range.low.dens[[count]] <- low.dens.removeIndLowDens[low] :low.dens.removeIndLowDens[b1+1]
                } else {
                    range.low.dens[[count]] <- low.dens.removeIndLowDens[b1+1]:low.dens.removeIndLowDens[b1+1]
                }
            }
        }
    }

    # does [[1]] or [[N]] contain the first or last point, if so remove them
    if ( length(which(range.low.dens[[length(range.low.dens)]] == nrow(f)) >= 1) )
        range.low.dens[[length(range.low.dens)]] <- NULL
    if ( length(which(range.low.dens[[1]] == 1) >= 1) )
        range.low.dens[[1]] <- NULL

    if (length(range.low.dens) != 0 ){
        # change to range
        range.low.dens <- lapply(1:length(range.low.dens), function(x){ range(range.low.dens[[x]]) })
        # add 2.5% each way to the range.
        if ( length(range.low.dens) >= 3 ) {
            range.low.dens.temp <- range.low.dens
            range.low.dens <- lapply(2:(length(range.low.dens)-1), function(x){print(x)
            # range.low.dens <- lapply(1:(length(range.low.dens)), function(x){print(x)
                range.temp <- range.low.dens[[x]][2]- range.low.dens[[x]][1]
                return( c( max(round(range.low.dens[[x]][1]-0.025*range.temp), 1),
                           min(round(range.low.dens[[x]][2]+0.025*range.temp), length(dens.f$y)) ) ) })
            range.low.dens <- c(list(range.low.dens.temp[[1]]), range.low.dens, list(range.low.dens.temp[[length(range.low.dens.temp)]]))
        }
        # change to time coordinates
        range.low.dens <- lapply(1:length(range.low.dens), function(x){ c(dens.f$x[range.low.dens[[x]][1]], dens.f$x[range.low.dens[[x]][2]]) } )
        # change to indices
        removeIndLowDens <- NULL
        for ( b2 in 1:length(range.low.dens) ){
            removeIndLowDens <- c(removeIndLowDens, intersect(which(f@exprs[,Time.loc] >= range.low.dens[[b2]][1]), which(f@exprs[,Time.loc] <= range.low.dens[[b2]][2])) )
        }

        if (length(removeIndLowDens) == 0 ){
            cat("None deleted from flowCut low dens removal.\n")
        } else {
            temp.text <- range.low.dens
            if( temp.text[[length(temp.text)]][1] > maxTime && temp.text[[length(temp.text)]][2] > maxTime ){ temp.text[[length(temp.text)]] <- NULL}
            if( temp.text[[1]][1] < minTime && temp.text[[1]][2] < minTime ){ temp.text[[1]] <- NULL}
            for( p1 in 1:length(temp.text)){
                temp.text[[p1]] <- paste(round(temp.text[[p1]], digits = 2), collapse = " to ")
            }
            temp.text <- paste(temp.text, collapse = ", ")

            cat(paste0("Removing low density ranges ", temp.text, ".\n"))
            f@exprs <- f@exprs[-removeIndLowDens,]
        }
    } else {
        cat("None deleted from flowCut low dens removal.\n")
         removeIndLowDens <- NULL
    }
    return(list(frame=f, rem.ind=removeIndLowDens))
}

##################################################
# Calculate means and which bins will be removed #
##################################################
calcMeansAndBinsRemoved <- function(f, segment, CleanChan.loc, FirstOrSecond){
    # FirstOrSecond - when we run this the second time we do not need to delete any cells, so we don't calculate that again.
    deletedSegments <- NULL
    meanRangePerc   <- NULL
    storeMeans <- list()
    quantiles <- list()
    totalNumSeg <- floor(nrow(f@exprs)/segment)
    Time.loc <- which(tolower(f@parameters@data$name) == "time"); names(Time.loc) <- NULL
    #For time
    timeCentres <- NULL
    for ( k in 1:(totalNumSeg-1)){
        temp <- f@exprs[segment*(k-1)+(1:segment), Time.loc]
        timeCentres[k] <- mean(temp)
    }
    k <- totalNumSeg
    temp <- f@exprs[(segment*(k-1)+1):nrow(f@exprs), Time.loc]
    timeCentres[k] <- mean(temp)
    for ( j in CleanChan.loc){
    # create matrix from summary for the arrayOutliers function
        segSummary <- matrix(0, totalNumSeg, 8)
        for ( k in 1:(totalNumSeg-1)){
            temp <- f@exprs[segment*(k-1)+(1:segment), c(j)]
            segSummary[k,] <- c( quantile(temp, probs = c(5,20,50,80,95)/100), mean(temp),
                                 sapply(2:3, function(x){ moment(temp, order = x, center = T)}) )
        }
        k <- totalNumSeg
        temp <- f@exprs[(segment*(k-1)+1):nrow(f@exprs), c(j)]
        segSummary[k,] <- c( quantile(temp, probs = c(5,20,50,80,95)/100), mean(temp),
                             sapply(2:3, function(x){ moment(temp, order = x, center = T)}) )

        storeMeans[[j]] <- segSummary[,6] # mean
        quantiles[[j]] <- quantile(f@exprs[,j], probs = c(0,2,98,100)/100)
        meanRangePerc[j] <- (max(storeMeans[[j]]) - min(storeMeans[[j]])) / (quantiles[[j]]["98%"]-quantiles[[j]]["2%"])
        if (FirstOrSecond == "First"){
            cellDelete <- apply(segSummary, 2, remove_outliers, limit = 3)
            cellDelete <- sapply(1:nrow(cellDelete), function(x) { sum(cellDelete[x,]) } )
            cellDelete <- which(cellDelete >= 4)

            # cellDelete <- arrayMvout::ArrayOutliers(data.frame(segSummary), alpha=alpha)

            if ( length(cellDelete)==0 ) {next}
            listDelete <- sort(cellDelete,decreasing=T)
            # listDelete <- sort(as.numeric(rownames(cellDelete$outl)),decreasing=T)
            deletedSegments <- c(deletedSegments, listDelete)
        }

    } # end of for-loop
    if (FirstOrSecond == "Second"){
        return(list(deletedSegments=0, quantiles=quantiles, storeMeans=storeMeans, meanRangePerc=meanRangePerc, timeCentres=timeCentres))
    }
    deletedSegments <- suppressWarnings(sort(unique(deletedSegments)))
    if ( length(deletedSegments) >= 1) {
        # removing single sliver
        deletedSegments2 <- deletedSegments
        for (l in 1:length(deletedSegments) ) {
            if (deletedSegments[l] == 1 || deletedSegments[l] == totalNumSeg ){
            } else {
                if ( length(which((deletedSegments[l] - 1) == deletedSegments)) >= 1 || length(which((deletedSegments[l] + 1) == deletedSegments)) >= 1 ){

                } else {
                    deletedSegments2[l] <- -1234
                }
            }
        }
        if (length(which(deletedSegments2 == -1234)) != 0){
            cat(paste0("Removed ", length(which(deletedSegments2 == -1234)), " slivers.\n"))
            deletedSegments <- deletedSegments2[-which(deletedSegments2 == -1234)]
        }
    }
    if(length(deletedSegments) != 0) {
        # Adding two buffer bins around consectutively removed bins.
        range.seg.rem <- list()
        count <- 1
        low <- 1
        if( length(deletedSegments) == 1){
            range.seg.rem[[1]] <- deletedSegments
        } else {
            for(b1 in 1:(length(deletedSegments)-1) ) {
                if( deletedSegments[b1] != (deletedSegments[b1+1]-1) ){
                    range.seg.rem[[count]] <- deletedSegments[low]:deletedSegments[b1]
                    low <- b1+1
                    count <- count + 1
                }
                if( b1 == (length(deletedSegments)-1) ) {# end point
                    if(deletedSegments[b1] == (deletedSegments[b1+1]-1)){
                        range.seg.rem[[count]] <- deletedSegments[low] :deletedSegments[b1+1]
                    } else {
                        range.seg.rem[[count]] <- deletedSegments[b1+1]:deletedSegments[b1+1]
                    }
                }
            }
        }
        size.in.a.row <- 50
        adding.segments <- NULL
        for ( b2 in 1:length(range.seg.rem) ){
            if(length(range.seg.rem[[b2]]) >= size.in.a.row){
                adding.segments <- c(adding.segments, (min(range.seg.rem[[b2]])-size.in.a.row/2):(min(range.seg.rem[[b2]])-1)
                                                    , (max(range.seg.rem[[b2]])+size.in.a.row/2):(max(range.seg.rem[[b2]])+1))
            }
        }
        if(length(adding.segments) >= 1){
            deletedSegments <- sort(unique(c(deletedSegments, adding.segments[intersect(which(adding.segments >=1), which(adding.segments <= totalNumSeg))])))
        }
    } else {
        deletedSegments <- NULL
    }
    return(list(deletedSegments=deletedSegments, quantiles=quantiles, storeMeans=storeMeans, meanRangePerc=meanRangePerc, timeCentres=timeCentres))
}
##################################################
# Calculate outliers to be removed               #
##################################################
remove_outliers <- function(x, limit = 3) {
    if(all(x == x[1])){
        return(rep(FALSE, length(x)))
    }
    # mn <- mean(x, na.rm = T)
    dens <- density(x)
    maxp <- dens$x[which(dens$y == max(dens$y))][1]
    out <- limit * sd(x, na.rm = T)
    x < (maxp - out) | x > (maxp + out)
}
##################################################
# Time Function                                  #
##################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
