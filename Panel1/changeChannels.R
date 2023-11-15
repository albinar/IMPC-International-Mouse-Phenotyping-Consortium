## Albina Rahim
## Last Updated: Feb 22, 2016

## This function gets rid of all extra channels in a flowFrame
## The function will identify the extra channels based on the channels of the original flowFrame & delete them
## @param f contains the list of all flowFrames (both the original and the flowFrames that needs to modification)
## @param extra.Markers contains the names of the extra markers that needs to be identified and removed
## @param index.FCS.originalMarkers contains the indices of the original flowFrames, which requires no modification 


library(flowCore)

changeChannels <- function(f, extra.markers, index.FCS.originalMarkers){
  
## Checking the channel positions of the flowFrames and if any swapping is required
## If required then call the swapChannels.R function
  FCS.Markers.position <- sapply(1:length(index.FCS.originalMarkers), function(j){
    length(which(colnames(exprs(f[[index.FCS.originalMarkers[1]]])) != colnames(exprs(f[[index.FCS.originalMarkers[j]]])), arr.ind = TRUE))})
  
  if(any(FCS.Markers.position !=0)){
    FCS.Markers.position <- which(FCS.Markers.position!=0)
    for(k in 1:length(FCS.Markers.position)){
      f[[index.FCS.originalMarkers[FCS.Markers.position[k]]]] <- swapChannels(f=f[[index.FCS.originalMarkers[FCS.Markers.position[k]]]], original.FCS=f[[index.FCS.originalMarkers[1]]])
    }
  }
  

## Getting rid of the extra channels (expres, description, parameter) and re-organizing the channels according to the original ones
f <- sapply(1:length(f), function(i){
    
	if(any(f[[i]]@parameters@data$desc %in% extra.markers == TRUE)){
		index.extraMarkers <- which(f[[i]]@parameters@data$desc %in% extra.markers)
		index.SPILL.extraMarker <- which(colnames(f[[i]]@description$SPILL) %in% extra.markers)
		
		ex <- exprs(f[[i]])[,-index.extraMarkers]
		
		description(f[[i]])$SPILL <- description(f[[i]])$SPILL[-index.SPILL.extraMarker,-index.SPILL.extraMarker]
		for(x in 1:length(index.extraMarkers)) {
	  	f[[i]]@description[grep(paste0("P",index.extraMarkers[x]), names(f[[i]]@description))] <- NULL
    }
		
	
		count <- 0

		for(x in 1:length(index.extraMarkers)){
			if(x==1){
				count <- index.extraMarkers[x]+x
				for(y in count:(index.extraMarkers[x+1]-x)){
					index.Desc <- grep(paste0("P",y), names(f[[i]]@description))
					for(z in 1:length(index.Desc)){
							names(f[[i]]@description)[index.Desc[z]] <- gsub(paste0("P",y), paste0("P",y-x), names(f[[i]]@description)[index.Desc[z]])

						}
				}
			} else {
			  count <- index.extraMarkers[x]+1
				#count <- 15+x
				for(y in count:(index.extraMarkers[x]+x)){
					index.Desc <- grep(paste0("P",y), names(f[[i]]@description))
					for(z in 1:length(index.Desc)){
							names(f[[i]]@description)[index.Desc[z]] <- gsub(paste0("P",y), paste0("P",y-x), names(f[[i]]@description)[index.Desc[z]])

						}
				}
			}

		}

		
		o <- parameters(f[[i]])@data[-index.extraMarkers,]
  	rownames(o) <- rownames(parameters(f[[index.FCS.originalMarkers[1]]])) 
		metaData <- varMetadata(parameters(f[[index.FCS.originalMarkers[1]]]))

		f[[i]] <- new("flowFrame", exprs=ex, parameters=new("AnnotatedDataFrame", o, varMetadata=metaData), description=description(f[[i]])[order(names(f[[i]]@description))])
	
	} else{
	  f[[i]] <- new("flowFrame", exprs=exprs(f[[i]]), parameters=parameters(f[[i]]), description=description(f[[i]])[order(names(f[[i]]@description))])

	}

  })

   	  return(f)
}
