#Removing margin events
#Debris= F default, removes debris events based on minimum value of the channel
#Neg=500 default, removes events based on neg value, when 500 it doesn't remove anything
#library(flowCut)
cleaned.fcs <- function(f, cleaning=T, transformed =T, id, segment=500,failFiles=T,
                        result.folder,comp.matrix=NULL, spill=F,PrintToConsole=F,make.NA=F,return.inds=F,return.all=T)
  #Cleans the FCS files of Anixa data using flowCut algorithm
  #Inputs:
  #f: raw FCS file
  #Cleaning: To perform cleaning or not, default TRUE
  #result.folder: Directory path to save the figures of flowCut.
  #Output
  #f.now: returns the cleaned raw data in flowFrame format
{
  if (cleaning)
  {
    if (spill & class(comp.matrix)=="matrix")
    {
      f1 <- compensate(f, comp.matrix)
    }else if (spill){
     f1<- compensate.flow(f, Manual.gating=NULL,comp=NULL)
        
    }else{
      f1<-f
    }
    if (!transformed){
      f1 <- estimate.logicle(as(f1,"flowSet"),med = F,m = 5,return.raw = F)[[1]]
    }else{
      f1<-f
    }
    res.tvsf <- flowCut(f1, Segment = segment,Directory = result.folder,PrintToConsole=PrintToConsole,
                        FileID = id, Plot ="All",Verbose = T)
  }else{
    res.tvsf <- list(ind=NULL,data= matrix(1,nrow = 16,ncol = 2,byrow = T))
  }
  if (return.all)
     {
       if (failFiles)
       {
         
       if(!any(as.logical( res.tvsf$data[c(3,5,7),1])))
       {
         if(file.exists(paste0(result.folder,"/FlaggedFiles_flowCut.csv")))
         {
           flagged.csv<- as.matrix(read.csv(paste0(result.folder,"/FlaggedFiles_flowCut.csv"),check.names = F))
         }else{
           flagged.csv <-t(as.matrix(c("Name","event","removed events","Flagged")))
         }
           res.tvsf <- c( res.tvsf,"FALSE")
           flagged.csv<- rbind(flagged.csv,cbind(identifier(f),dim(f@exprs)[1],length(res.tvsf$ind),1))
           write.csv(flagged.csv,file=paste0(result.folder,"/FlaggedFiles_flowCut.csv"),row.names=F)
         }else {
         res.tvsf <- c( res.tvsf,"TRUE")
         }
       }
         return(res.tvsf)
       
     }
    if (return.inds)
      return(res.tvsf$ind)
    if(length(res.tvsf$ind)>0){
      if (!make.NA)
      {   f.new <-f[-res.tvsf$ind,]
      }else{
        exprs(f)[res.tvsf$ind,]<-NA
        f.new<-f
      }
    }else{
      f.new <-f
    }
  return(f.new)
}
#################################################################
#sampling from all flowFrames
#################################################################

getGlobalFrame <- function(fs, sample.length=NA, all.cells=F){
  set.seed(123)
  if (is(fs, 'flowFrame')){
    return (frame)
  }
  n <- length(fs)
  sample.n <- ifelse(is.na(sample.length),n,sample.length)
  global.frame <- fsApply(fs[sample(n, sample.n)],
                          function(frame){
                            m <- nrow(frame)
                            sample.size <- ifelse (all.cells,yes = m,no =ceiling(m/sample.n ))
                            
                            exprs(frame )<- exprs(frame)[sample(m, sample.size),]
                            return (frame)
                          })
  global.frame <- as(global.frame, 'flowFrame')
  return (global.frame)
}

#################################################################
#Data rotation
#################################################################
rotate.data <- function(data, chans=NULL, theta=NULL,min.max=F)
{
  if(nrow(data)<3)
  {
    print("Cannot rotate a matrix with less than 3 rows, returning back the original matrix")
    return(list(data=data,theta=theta))
  }else{
    if (class(data)== "flowFrame" & !is.null(chans))
    {
      if(all(is.na(exprs(data)[,1])))
        return("Cannot rotate a flowFrame with all cells NA")
      no.nas <- which(!is.na(exprs(data)[,chans[1]]))
      data.new <- exprs(data)[no.nas,chans]
      if (is.null(theta))
      {
        reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
        slope <-  atan((max(data.new[,2])-min(data.new[,2]))/(max(data.new[,1])-min(data.new[,1])))
        theta <- ifelse(min.max,no = pi/2-reg.slope,yes = slope-pi/2)
      }
      data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
      exprs(data)[no.nas,chans] <- data.new
    }else{
      col.names<- colnames(data)
      data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
      colnames(data)<-col.names
    }
    return(list(data=data,theta=theta))
  }
}

#################################################################
#flowDensity rotation
#################################################################

#Rotating back a flowDensity object
rotate.fd <- function(fd.object, angle)
{
  new.f <-new(Class = "CellPopulation")
  no.na<-which(!is.na(exprs(fd.object@flow.frame)[,1]))
  dat <- fd.object@flow.frame
  temp <-rotate.data(dat[no.na,], fd.object@channels,theta=-angle)$data
  exprs(dat)[no.na,] <-exprs(temp)
  new.f@flow.frame<-dat
  new.f@filter <-rotate.data(fd.object@filter ,theta =-angle)$data
  colnames(new.f@filter)<-colnames(fd.object@filter )
  new.f@channels<-fd.object@channels
  new.f@proportion<-fd.object@proportion
  new.f@cell.count <-fd.object@cell.count
  new.f@index<-fd.object@index
  return(new.f)
}

#################################################################
#Margin removal function
#################################################################
removeMargins<- function(f,chans,sens=1, debris=FALSE,neg=500, verbose = T,return.gate=F)
{
  neg <-rep(neg,length(chans))
  data <- exprs(f)
  margins <- c()
  marg.gates<-c()
  if(!debris)
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.max <-max(data[,chan],na.rm = T)
      margins <- which ( data[, chan] >= stain.max*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
      marg.gates <- append(marg.gates,stain.max*sens-1 )
    }
    if(return.gate)
      return(marg.gates)
  }else
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.min <-min(data[,chan],na.rm=T)
      margins <- which ( data[, chan] <= stain.min*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
  }
  for(i in 1:length(chans))
  {
    if (neg[i]<500)
    {
      ch<-ifelse (is.character(chans[i]),yes = which(colnames(f)==chans[i]),no = chans[i])
      negs <- which ( data[, ch] < neg[i])
      margins <- negs
      if (length(margins)!=0){
        data <- data[ -margins, ]
      }
      if(verbose == T){print(paste(length(margins), "negative events in",colnames(f)[ch], "will be removed.",sep =" "))}
    }
  } 
  exprs(f) <- data
  
  return(f)
}
#################################################################
#Compensation
#################################################################
compensate.flow <- function(data, Manual.gating=NULL,comp=NULL)
  #data: flowSet, flowFrame or a GatingSet
  #Manual.gating: if data is of class GatingSet, you need to provide the manual gatingSet, otherwise the built in Comepnsation matrix will be used.
{
  if (!is.null(comp))
    {
    if (class(data) %in%c("flowFrame" ))
{
      print("Compensating from comp matrix")
      
      return(compensate(data, comp))
    }else if (class(data)=="GatingSet")
    {
      comp <- compensation(comp,compensationId="comp1")
      print("Compensating from comp matrix")
      return(compensate(data, comp))
    }else{
      return(fsApply(data, compensate, comp))
    }
  }else{
  if (class(data)=="GatingSet")
  {
      obj <- getData(data)
      if (!is.null(Manual.gating)){
        comp <- getCompensationMatrices(Manual.gating[[1]])
      }else{
        comp <-lapply(1:length(obj), function(x) {
          sp<-spillover(obj[[x]])
          print("Compensating from FCS file.")
          return(compensation(sp[!sapply(sp,is.null)][[1]]))})
        names(comp)<- sampleNames(data)
      }
      return(compensate(data,comp))
  }else
  {
   
      if(class(data)=="flowFrame")
      {
        sp<-spillover(data)
        comp <- compensation(sp[!sapply(sp,is.null)][[1]])
      }else{
        comp <-lapply(1:length(data), function(x) {
        sp<-spillover(data[[x]])
        return(compensation(sp[!sapply(sp,is.null)][[1]]))})
  }
    if (length(comp)==1)
      return(compensate(data,comp))
    else
    {
      names(comp)<- sampleNames(data)
      return(fsApply(data, function(x) (compensate(x,comp[[identifier(x)]]))))
    }
  }
  }
}


#################################################################
#Creating estimateLogicle transformation filter for flowFrame, flowSet or GatingSet
#################################################################
transform.flow<-function(obj,remove.outliers=T,sd.coeff=4.5,trans.chans=NULL)
{
  #obj: flowFrame, flowSet, or a GatignSet
  #remove.outlier: Default to T, removing far away cells.
  #sd.coeff: Default to 4.5 for removing cells that are below or abour 4.5*sd
  #trans.chans: Channels to be transformed. if it's null, then it tries to find LOG channel in the frame
  if (class(obj)=="flowFrame")
  {
    f <- obj
    f.t <- f
  }else{
    temp <-  obj
    if (class(obj)=="GatingSet")
      temp <- getData(obj, tail(getNodes(obj),1))
    f.t <- temp[[1]]
    f<-getGlobalFrame(temp)
    
  }
  if (is.null(trans.chans))
  {
    log.channels <- paste("P",1:ncol(f.t),"DISPLAY",sep="")
    trans.chans <- which(f.t@description[log.channels]=="LOG")
    if (length(trans.chans)==0)
    {
      warnings("Couldn't find Log channels, all channels except FSC, SSC, and time will be transformed.")
      trans.chans <- 1:ncol(f.t)  
      trans.chans <- trans.chans[-c(grep(colnames(f.t),pattern = "FSC*"),grep(colnames(f.t),pattern = "SSC*"),grep(colnames(f.t),pattern = "Time*"))] 
    }
  }
  if (remove.outliers)
    # Using Justin's code to remove outliers
  {
    low_events <- high_events <- NULL
    time.loc <- which(tolower(f@parameters@data$name) == "time"); names(time.loc) <- NULL
    inds <- c()
    for (p1 in trans.chans){
      
      mean_temp <- mean(f@exprs[,p1])
      sd_temp <- sd(f@exprs[,p1])
      low_events <- c(low_events, which(f@exprs[,p1] <= (mean_temp - sd.coeff*sd_temp)) )
      high_events <- c(high_events, which(f@exprs[,p1] >= (mean_temp + sd.coeff*sd_temp)) )
    }
    ind.marg.neg <- unique(sort(c(low_events, high_events)))
    f.snipped <- f[- ind.marg.neg]
  }else{
    f.snipped <- f
  }
  if (class(obj)=="GatingSet")
  {
    data<- GatingSet(as(f.snipped,Class = "flowSet"))[[1]]
    
  }else
  {
    data <- f.snipped
  }
  #trans <- flowJo_biexp_trans()
  #trans <- transformerList(colnames(f)[trans.chans], trans)
  trans <- tryCatch(estimateLogicle(data, colnames(f)[trans.chans]), error=function(x) {return(1)})
  if (mode(trans)=="numeric")
    trans <- estimateLogicle(data, colnames(f)[trans.chans],type="data")
  return(transform(obj,trans))
}
###########################################################################
# Gating Singlet gate
###########################################################################
singlet.gate <- function(f, channels, angle=-pi/4,alpha=0.05)
{
 rot <- rotate.data(f, channels,theta = angle)$data
  gate.2 <- deGate(rot, channels[1],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=T, alpha=alpha,verbose = F)
  singlets<- rotate.fd(flowDensity(rot,channels,position = c(F,NA),gates=c(gate.2,NA),verbose = F),angle = angle)
  return(singlets)
}


###########################################################################
#
###########################################################################
#' Estimates a common logicle transformation for a flowSet.
#'
#' Of the negative values for each channel specified, the median of the specified
#' quantiles are used.
#'
#' @param flow_set object of class 'flowSet'
#' @param channels character vector of channels to transform
#' @param m TODO -- default value from .lgclTrans
#' @param q quantile
#' @return TODO
estimateMedianLogicle <- function(flow_set, channels, m = 4.5, q = 0.05) {
  if (!is(flow_set, "flowSet")) {
    stop("flow_set has to be an object of class 'flowSet'")
  }
  if (missing(channels)) {
    stop("Please specify the channels to be logicle transformed")
  }
  indx <- channels %in% unname(colnames(exprs(flow_set[[1]])))
  if (!all(indx)) {
    stop(paste("Channels", channels[!indx], "were not found in flow_set "))
  }
  
  neg_marker_quantiles <- fsApply(flow_set, function(sample) {
    apply(exprs(sample), 2, function(markers) {
      quantile(markers[markers < 0], probs = q)
    })
  })
  # Replaces 'r' in flowCore:::.lgclTrans
  neg_marker_quantiles <- apply(neg_marker_quantiles, 2,
                                median, na.rm = TRUE)[channels]
  
  # In the case that no negative markers are present, we set this quantile to the
  # default value of 1/2.
  neg_marker_quantiles <- replace(neg_marker_quantiles,
                                  is.na(neg_marker_quantiles), 0.5)
  
  # Replaces 't' in flowCore:::.lgclTrans
  max_range <- do.call(rbind, lapply(fsApply(flow_set, range), function(x) {
    x[2, channels]
  }))
  max_range <- apply(max_range, 2, max)
  
  # Replaces 'w' in flowCore:::.lgclTrans
  w <- (m - log10(max_range / abs(neg_marker_quantiles))) / 2
  if (any(w<0))
    w[which(w<0)]<-.5
  transformation <- lapply(channels, function(channel) {
    transId <- paste(channel, "medianLogicleTransform", sep = "_")
    
    logicleTransform(transformationId = transId, w = w[channel],
                     t = max_range[channel], m = m, a = 0)
  })
  
  transformList(channels, transformation,
                transformationId = "medianLogicleTransform")
}

###########################################################################
#Transformation for a flowSet
###########################################################################

estimate.logicle<-function(fs.raw,talk= TRUE,return.set=TRUE,med=TRUE, m=NA,trans.chans=NULL,estimate=T,return.raw=F)
  
{
  temp <- fsApply(fs.raw, function(x){
    inds <- which(is.na(exprs(x)[,1]))
    y<-x
    if (length(inds)>0)
      y <- x[-inds,]
    return(list(frame=y, na.inds=inds))
  })
  names(temp) <- sampleNames(fs.raw)
  fs <- as(lapply(temp, function(x)
    return(x$frame)),"flowSet")
  sampleNames(fs ) <- sampleNames(fs.raw)
  f<-fs[[1]]
  if (is.null(trans.chans))
  {
    log.channels <- paste("P",1:ncol(f),"DISPLAY",sep="")
    trans.chans <- which(f@description[log.channels]=="LOG")
    if ( length(trans.chans) ==0){
      trans.chans <- setdiff (1:length(f@parameters@data$desc),
                              c(grep("fsc", tolower(f@parameters@data$desc)),
                                grep("ssc", tolower(f@parameters@data$desc)),
                                grep("time", tolower(f@parameters@data$desc)) ) )
    }
    if(any(is.null(trans.chans)) | any(is.na(trans.chans)) | (length(trans.chans)==0))
    {
      return(cat("You have to find the channels manually, there's no information on FCS file","\n"))
    }
  }
  if( talk)
    cat("Channels to be transformed ", trans.chans, "\n")
  if(estimate)
  {
    if (med==TRUE){
      lgl <-estimateMedianLogicle(flow_set=fs,channels=colnames(f)[trans.chans])
      
    }else if(!is.na(m)){
      
      lgl <- estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans],m=m)
    }else{
      
      lgl<-tryCatch(estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans]), error=function(x) {return(1)})
      if (mode(lgl)=="numeric")
        lgl <- estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans],type="data")
    }
    if(return.set)
    {
      if(!return.raw)
      {
        return(fsApply(fs,function(x) transform(x,lgl)))
      }else{
        fs.temp <-fsApply(fs.raw,function(x) {
          tr <- transform(fs[[identifier(x)]],lgl);
          exprs(x)[-temp[[identifier(x)]]$na.inds,]<-exprs(tr)
          w<-Make.FCS(markers = as.vector(x@parameters@data[,2]),data = exprs(x))
          return(w)
        })
      }
    }
    else
      
      return(lgl)
  }else{
    print("Logicle should of been used but it is commented out")
    library(Logicle)
    trans.fs <-fsApply(fs,function(frame) {
      lg<-Logicle::create(T=262144, W=0.5)
      x<-exprs(frame)[,trans.chans]
      y<-Logicle::scale(lg, x)
      exprs(frame)[,trans.chans]<-y
      # detach(package:Logicle, unload=TRUE)
      return (frame)
    })
  }
}
##########################
##Finds markers in the FCS file
Find.markers <- function(frame,marker.list)
{
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}
########################################################
########################################################
#Outlier detection and replacements after usingdeGate for all FCS files in the flowSet

averageGates <- function(vec, med =T, global.var= NA,sd.coeff = 2, name_of_gate="-"){
  # Args:
  #   vec: a numeric vector of values
  #   med: to use median
  #   global.var: a number used instead of median
  # Value:
  #   vec with outlying values replaced by the median
  # (http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation)
  #s = sd(x)/c4(N) is estimator of SD in vector
  vec[which(is.infinite(vec))]<-NA
  vec <- as.vector(vec,mode = "numeric")
  m <- ifelse(med,median(vec,na.rm=T),global.var)
  if ( any(is.na(vec)) ){
    cat("Changing NAs", which(is.na(vec)), "from", name_of_gate, "\n", sep=" ")
    vec[is.na(vec)] <- m
  }
  
  N <- length(vec)
  c4 <- sqrt(2/(N-1)) * gamma(N/2)/gamma((N-1)/2)
  sdev <- sd(vec,na.rm=T)/c4
  outliers <- which(abs(vec - m)/sdev > sd.coeff)
  if (any(outliers)){
    cat("Changing outliers", outliers, "from", name_of_gate, "\n", sep=" ")
    vec[outliers] <- median(vec[-outliers])
  }
  return (vec)
}
#############################################################
# Using Kde2d density to find shoulders
#############################################################

kde2d.density <- function(fframe, channels, which.dim=2){
  data.new <- na.omit(exprs(fframe)[, channels])
  grid.points<-c(512,50)
  if (which.dim==2)
   grid.points<-rev(grid.points)

  z <- kde2d(data.new[, 1], data.new[, 2], n = grid.points)
  maxDens <- density(fframe@exprs[, channels[which.dim]], na.rm = T)
  maxDens2 <- maxDens
  if(which.dim==1)
  {
  maxDens2$y <- apply(z$z,1,max)
  maxDens2$x <- z$x
  }else{
    maxDens2$x <-z$y
    maxDens2$y <-  apply(z$z,2,max)
  }
  return(maxDens2)
}  

#############################################################
#Returning the coordinates of the outermost contour line
#############################################################
contour.line <- function(frame, channels,which.line=1)
{
  library(MASS)
  new.f <- frame
  dens.2d <- kde2d(exprs(new.f)[,channels[1]],exprs(new.f)[,channels[2]])
  cont <- contourLines(dens.2d)
  cont.coord <- cbind(cont[[which.line]]$x,cont[[which.line]]$y)
  return(cont.coord)
}


#############################################################
#Polygon's differeces
#############################################################
gate.diff <- function(filter1,filter2,subtract=TRUE)
{
  library(rgeos)
  library(sp)
  poly1 <- SpatialPolygons(list(Polygons(list(Polygon(filter1)), ID=c("c"))))                       
  poly2 <- SpatialPolygons(list(Polygons(list(Polygon(filter2)), ID=c("c"))))   
  filt<-filter1
  if(subtract){
    if(rgeos::gIntersects( poly1, poly2)){
      temp <-tryCatch(rgeos::gDifference(poly1,poly2), error=function(ex) return(1))
      if (mode(temp)=="numeric")
      {
        temp <- rgeos::gDifference(poly1,rgeos::gBuffer(poly2))
      }
      filt<-temp@polygons[[1]]@Polygons[[1]]@coords
    }
    return(filt)
  }else{
    
    return(gUnion(poly1,poly2) )}
}

#############################################################
#Showing time processed for analysis
#############################################################
fsApply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- fsApply(X, wrapper, ...)
  close(pb)
  res
}
############################################################
#
############################################################
# res.flowCut <- flowCut(f.trans, Directory = flowCut.fig.dir, LowDensityRemoval=0.05, MaxPercCut = 0.4, 
#                        UseOnlyWorstChannels = T, 
#                        FileID = paste0(basename(dirname(f@description$FILENAME)), '_',basename(f@description$FILENAME)), 
#                        Plot="All", AllowFlaggedRerun = T,
#                        AmountMeanRangeKeep = 0, AmountMeanSDKeep = 4)
# 
# # Save FCS file which contains one extra channel that contains a 0/1 depending on
# # if those cells were removed/kept by flowCut
# ind <- which(fnames ==  f@description$FILENAME)
# f.ind <- getIndices(gs[[ind]], "Margin")
# 
# ind.toremove.flowCut <- setdiff(1:nrow(fs.original[[ind]]), which(f.ind == FALSE))
# ind.toremove.flowCut <- ind.toremove.flowCut[res.flowCut$ind]
# 
# QCvals <- 1:nrow(fs.original[[ind]])
# QCvals <- as.integer(QCvals %in% ind.toremove.flowCut)
# 
# f.QC <- addQC(QCvals, exprs(fs.original[[ind]]), parameters(fs.original[[ind]]), keyword(fs.original[[ind]]))
# f.QC@description$SPILL <- as.matrix(comp) # Put correct compensation matrix into the saved FCS file
# write.FCS(f.QC, paste0(flowCut.dir, basename(dirname(f@description$FILENAME)), '_',basename(f@description$FILENAME)))






## create new flowFrame with the parameter indicating good and bad cells
addQC <- function(QCvals, sub_exprs, params, keyval){
  
  rs <- attr(sub_exprs, "ranges")
  rs <- c(rs, rs[1])
  sub_exprs <- cbind(sub_exprs, QC = QCvals)
  attr(sub_exprs, "ranges") <- rs
  NN <- as.numeric(keyval["$PAR"]) + 1
  names(dimnames(sub_exprs)[[2]]) <- sprintf("$P%sN", 1:NN)
  pnr <- paste0("$P", NN, "R")
  pnb <- paste0("$P", NN, "B")
  pne <- paste0("$P", NN, "E")
  pnn <- paste0("$P", NN, "N")
  pns <- paste0("$P", NN, "S")
  flowCorePnRmax <- paste0("flowCore_$P", NN, "Rmax")
  flowCorePnRmin <- paste0("flowCore_$P", NN, "Rmin")
  o <- params@data
  o[length(o[,1]) + 1,] <- c("QC", "bad = 1", as.numeric(keyval$`$P1R`), 0, 1)
  rownames(o)[length(o[,1])] <- paste("$P", NN, sep = "")
  
  outFCS <- new("flowFrame", exprs=sub_exprs, parameters=new("AnnotatedDataFrame",o), description=keyval)
  description(outFCS)[pnr] <- max(1, description(outFCS)$`$P1R`)
  description(outFCS)[pnb] <- description(outFCS)$`$P1B`
  description(outFCS)[pne] <- "0,0"
  description(outFCS)[pnn] <- "QC"
  description(outFCS)[pns] <- "bad = 1"
  description(outFCS)$`$PAR` <- NN
  description(outFCS)[flowCorePnRmax] <- 1
  description(outFCS)[flowCorePnRmin] <- 0
  
  return(outFCS)
}  


#####################################################################
#####################################################################
#From flowDensity package
.densRange <- function(x, y, gate, pos = FALSE){
  
  ##==================================================================
  ## Plots the output of density-estimate gating method
  ##  Args:
  ##   x: 'x' slot of the density returned by the 'density()' function
  ##   y: 'y' slot of the density returned by the 'density()' function
  ##   gate: a threshold given by the '.densityGating()' function
  ##   pos: refer to the 'pos' in 'deGatePlot()' function
  ## Value:
  ##   (x,y) coordinates
  ##------------------------------------------------------------------
  pts <- list()
  if(is.na(pos))
    return(list(x=c(x,tail(x,1),x[1]), y=c(y,min(y),min(y))))
  if(pos){
    x.pts <- c(x[which(x>=gate)], tail(x[which(x>=gate)],1), x[which(x>=gate)][1])
    y.pts <- c( y[which(x>=gate)], min(y[which(x>=gate)]), min(y[which(x>=gate)]))
  }else{
    x.pts <- c(x[which(x<gate)][1], x[which(x<gate)], gate)
    y.pts <- c(min(y[which(x>=gate)]),y[which(x<gate)], min(y[which(x<gate)]))
  }
  pts$x <- x.pts
  pts$y <- y.pts
  return(pts)
}
##########################################################################3
##########################################################################
Make.FCS<- function(markers, data, f.guid=NULL)
{
  pd <- c()  # 'params' phenoData
  des <- list()  # 'description' list
  
  des[["$DATATYPE"]] <- "F"
  for (c in 1:ncol(data)) {
    c_name <- colnames(data)[c]
    c_marker<-markers[c]
    c_min <- floor(min(data[,c],na.rm = T))
    c_max <- ceiling(max(data[,c],na.rm = T))
    c_rng <- c_max - c_min + 1
    
    pl <- matrix(c(c_name, c_marker, c_rng, c_min, c_max),nrow=1)
    colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
    rownames(pl) <- paste("$P",c,sep="") 
    pd <- rbind(pd, pl)
    
    des[[paste("$P",c,"B",sep="")]] <- "32";      # Number of bits
    des[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
    des[[paste("$P",c,"E",sep="")]] <- "0,0";      # Exponent
    des[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
    des[[paste("$P",c,"S",sep="")]] <- c_marker;	    # Desc	
  }
  frame<-flowFrame(data, as(data.frame(pd), "AnnotatedDataFrame"), description=des)
  if (!is.null(f.guid))
  {
    frame@description$GUID<- f.guid
  }
  return(frame)
}