#####################################################
#Create boxplot: show only boxes using geom_boxplot_interactive (on allprops) (hide outliers with argument 'outlier.shape=NA'),
#and add outliers as separate dataset using geom_point_interactive (plots as red dots with tooltips printing filename)

plot_results.link <- function(alldata,outlierdata,popstoplot,link){
  
    newdata <- alldata[1,]
    for(i in 1:length(popstoplot)){
      newdata <- rbind(newdata,alldata[which(alldata$population==popstoplot[i]),])
    }
    alldata <- unique(newdata)
    lim1 <- floor(min(alldata$proportion,na.rm = T))
    lim2 <- ceiling(max(alldata$proportion,na.rm = T)+1)
    by.factor <- ifelse(lim2>50,yes = 10,no=ifelse(lim2>30,yes = 5,no=ifelse(lim2>10,yes = 4,no = 2)))
  if(length(outlierdata[,1]) > 0){
    outlierdata$file <- unlist(lapply(outlierdata$file, function(file)
      gsub(unlist(strsplit(file , " /"))[1],pattern = ".fcs",replacement = "")))
    #alldata$link<-paste0('window.open("', link , '")')
    ind <- match(outlierdata$file,gsub(link$Title,replacement = "",pattern = ".fcs_Automated.png"))
    outlierdata$link <-  paste0('window.open("',link$URL[ind]
                                    , '")')
    
    r <- ggplot() +
      geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                              "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
      geom_point_interactive(position=position_jitter(width=0.45,height=.45),size=.7,data=outlierdata,
                             aes(x = population, y = proportion,colour=population, onclick=link,
                                 tooltip=file,data_id=file))+
      scale_y_continuous(name="Proportion % (parent population)",limits = c(lim1,lim2),
                         breaks=seq(lim1,lim2,by=by.factor),expand=c(0,0))+
      xlab(label = "Cell Population")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    tooltip_css <- "background-color:transparent;font-style:italic;"
    
    r <- girafe(ggobj=r,options = list(
      opts_tooltip(css = tooltip_css),
      opts_sizing(width = .7) ))
    
    r <- girafe_options(r, opts_toolbar(saveaspng = TRUE),
                        #specifies properties of outlier dots when cursor hovers over them
                        opts_hover(css = "fill:black;stroke:red;stroke-width:5pt;r:3pt"))
  }
  
  #if there are no outlier points... no need to add geom_point_interactive.
  if(length(outlierdata[,1])==0){
    r <- ggplot() +
      geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=factor(population)+1,
                                                tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                              "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
      scale_y_continuous(name="Proportion % (parent population)",limits = c(lim1,lim2),breaks=seq(lim1,lim2,by=by.factor),expand=c(0,0))+
      xlab(label = "Cell Population")+

      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    r <- girafe(ggobj=r)

  }
  return(r)
  
 
}


#####################################################
#Create boxplot: show only boxes using geom_boxplot_interactive (on allprops) (hide outliers with argument 'outlier.shape=NA'),
#and add outliers as separate dataset using geom_point_interactive (plots as red dots with tooltips printing filename)

plot_results<- function(alldata,outlierdata,popstoplot){
  
  newdata <- alldata[1,]
  for(i in 1:length(popstoplot)){
    newdata <- rbind(newdata,alldata[which(alldata$population==popstoplot[i]),])
  }
  alldata <- unique(newdata)
  lim1 <- floor(min(alldata$proportion,na.rm = T))
  lim2 <- ceiling(max(alldata$proportion,na.rm = T)+1)
  by.factor <- ifelse(lim2>50,yes = 10,no=ifelse(lim2>30,yes = 5,no=ifelse(lim2>10,yes = 4,no = 2)))
  if(length(outlierdata[,1]) > 0){
    outlierdata$file <- unlist(lapply(outlierdata$file, function(file) unlist(strsplit(file , " /"))[1]))
    
    r <- ggplot() +
      geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                              "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
      geom_point_interactive(position=position_jitter(width=0.45,height=.45),size=.7,data=outlierdata,
                             aes(x = population, y = proportion,colour=population,
                                 tooltip=paste0(round(proportion,digits=2),"%, ",file),data_id=file))+
      scale_y_continuous(name="Proportion % (parent population)",limits = c(lim1,lim2),
                         breaks=seq(lim1,lim2,by=by.factor),expand=c(0,0))+
      xlab(label = "Cell Population")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    r <- girafe(ggobj=r)
    r <- girafe_options(r, opts_toolbar(saveaspng = TRUE),
                        #specifies properties of outlier dots when cursor hovers over them
                        opts_hover(css = "fill:black;stroke:red;stroke-width:5pt;r:3pt"))
  }
  
  #if there are no outlier points... no need to add geom_point_interactive.
  if(length(outlierdata[,1])==0){
    r <- ggplot() +
      geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=factor(population)+1,
                                                tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                              "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
      scale_y_continuous(name="Proportion % (parent population)",limits = c(lim1,lim2),breaks=seq(lim1,lim2,by=by.factor),expand=c(0,0))+
      xlab(label = "Cell Population")+
      
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
            axis.text.y = element_text(hjust = 1, size=10,color=1),
            legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
    r <- girafe(ggobj=r)
    
  }
  return(r)
  
  
}
#####################################################
#Boxplot Creation

boxplot.svg <- function(gs,sdconst=3,popstoplot,rows.to.avoid=NULL,row.name.no=2,percentage=T,click=F,link=NA)
{
# 1) If gs is a gatingSet, then popstoplot are nodes to be plotted from the gs
# 2) If gs is a matrix of porps or data.frame of props, then columns correspond to samples and rows to populations. 
# And popstoplot corresponds to rownames of gs
# 3) If gs is a vector of path to multiple csv file, then they will be read and merged together, so make sure they have same numer of rows, and are consistent.
# And popstoplot corresponds to rownames of each of these csv files
# And rows.to.avoid, excludes the rows that are not related to samples, and should not be included in the boxplot
# And row.name.no will be passed to read.csv to grab the rownames from the csv file.

if (class(gs)=="GatingSet")
{
   stats <- gs_pop_get_count_fast(gs,format="wide",statistic="freq",path="auto")[popstoplot,]
   tmp <- t(stats)
}else if (class(gs)=="matrix" |class(gs)=="data.frame")
{
  tmp<- gs[,popstoplot]
}else if (class(gs)=="character" & length(grep(gs,pattern = ".csv"))==length(gs))
{
    if (is.null(rows.to.avoid))
        tmp <- do.call(rbind,lapply(gs, function(path) read.csv(path, check.names=F,row.names=row.name.no)))[,popstoplot]
    else
        tmp <- do.call(rbind,lapply(gs, function(path) read.csv(path, check.names=F,row.names=row.name.no)[-rows.to.avoid,popstoplot]))
}else{
  stop("Unsupported input as gs, check the documentation of boxplot.svg")
}

propmatrix <- stack(as.data.frame(tmp))
propmatrix$file <- rownames(tmp)
names(propmatrix) <- c("proportion","population","file")
if (!percentage)
      propmatrix$proportion <- round(propmatrix$proportion*100,digits = 5)
#remove empty rows
propmatrix <-propmatrix[order(propmatrix$population,propmatrix$file,propmatrix$proportion),]
delete <- which(propmatrix$population == "")
if(length(delete) >0){
  propmatrix <- propmatrix[-delete,]
}
populationNumber <- length(unique(propmatrix$population))
propmatrix$median <- ""
propmatrix$mean <- ""

outlierpropmatrix <- propmatrix[1,]
outlierpropmatrix[1,] <- NA
idx <- 1
for(i in unique(propmatrix$population)){
  if(i != ""){
    indices <- which(propmatrix$population==i)
    mean <- mean(x=as.numeric(propmatrix$proportion[indices]))
    propmatrix$mean[indices] <- mean
    median <- median(as.numeric(propmatrix$proportion[indices]))
    propmatrix$median[indices] <- median
    sd <- sd(as.numeric(propmatrix$proportion[indices]))
    upplim <- mean + sd*sdconst
    lowlim <- mean - sd*sdconst
    
    outliers <- which(as.numeric(propmatrix$proportion[indices]) > upplim)
    outliers <- append(outliers,values=which(as.numeric(propmatrix$proportion[indices]) < lowlim))
    if(length(outliers) >0){
      for(k in outliers){
        outlierpropmatrix[idx,] <- propmatrix[indices[k],]
        idx <- idx+1      
      }
    }
  }
}

if(outlierpropmatrix$population[1] == "" || is.na(outlierpropmatrix$population[1])){
  outlierpropmatrix <- outlierpropmatrix[-which(outlierpropmatrix$population==""),]
}

library(ggplot2)
library(ggiraph)
if(click){
  boxplot <- plot_results.link(alldata=propmatrix,outlierdata=outlierpropmatrix,popstoplot=popstoplot,link=link)
  boxplot
}else{
boxplot <- plot_results(alldata=propmatrix,outlierdata=outlierpropmatrix,popstoplot=popstoplot)
boxplot
}
}
