## Albina Rahim
## Last Updated: Feb 22, 2016

## This function swap channels in a flowFrame between 2 or more positions
## This function can be used without knowing which channels needs to be swapped. 
## The function will identify the channels that needs to be swaped and swap it based on the channel positions of the original flowFrame
## @param f contains the flowFrame which requires channels swapping
## @param original.FCS contains the original flowFrame, which can be used as a template for swapping the channels of the flowFrame: f 


swapChannels <- function(f, original.FCS){
	    temp.FCS <- f
	    index.toChangeto <- which(colnames(exprs(original.FCS)) != colnames(exprs(f)), arr.ind = TRUE)
	    name <- colnames(exprs(original.FCS))[index.toChangeto]
	    index.toChangefrom <- sapply(1:length(name), function(i){which(colnames(exprs(f)) == name[i])})
	    storeInfo <- cbind(index.toChangeto, name, index.toChangefrom)
	    
	    # Changing the expression of the temp.FCS flowframe 
	    exprs(temp.FCS)[,as.numeric(storeInfo[,1])] <- exprs(f)[,as.numeric(storeInfo[,3])]
	    colnames(exprs(temp.FCS))[as.numeric(storeInfo[,1])] <- storeInfo[,2]
	    
	    # Changing the spillover matric of the temp.FCS flowframe
	    index.toChangeto.SPILL <- which(colnames(description(original.FCS)$SPILL) != colnames(description(f)$SPILL), arr.ind = TRUE)
	    name.SPILL <- colnames(description(original.FCS)$SPILL)[index.toChangeto.SPILL]
	    index.toChangefrom.SPILL <- sapply(1:length(name), function(i){which(colnames(description(f)$SPILL) == name.SPILL[i])})
	    storeInfo.SPILL <- cbind(index.toChangeto.SPILL, name.SPILL, index.toChangefrom.SPILL)
	    description(temp.FCS)$SPILL[,as.numeric(storeInfo.SPILL[,1])] <- description(f)$SPILL[,as.numeric(storeInfo.SPILL[,3])]  
	    description(temp.FCS)$SPILL[as.numeric(storeInfo.SPILL[,1]),] <- description(temp.FCS)$SPILL[as.numeric(storeInfo.SPILL[,3]),]  
	    colnames(description(temp.FCS)$SPILL)[as.numeric(storeInfo.SPILL[,1])] <- storeInfo.SPILL[,2] 
	   
	    # Changing the parameters of the temp.FCS flowframe
	    parameters(temp.FCS)@data[as.numeric(storeInfo[,1]),] <- parameters(f)@data[as.numeric(storeInfo[,3]),]
	    rownames(parameters(temp.FCS)@data) <- rownames(parameters(f)) 
	    metaData <- varMetadata(parameters(f))
	    
	    f <- new("flowFrame", exprs=exprs(temp.FCS), parameters=new("AnnotatedDataFrame", parameters(temp.FCS)@data, varMetadata=metaData), description=description(temp.FCS))
	
	return(f)
}	    
	
