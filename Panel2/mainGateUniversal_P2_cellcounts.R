# Read in all the cell counts and then bind them into one big list.

input.dir <- '/data/projects/Codes/IMPC/Panel2/results'
setwd('/data/projects/Codes/IMPC/Panel2')
fnames <- list.files(input.dir, full.names = T)

load(fnames[1])
cell.counts.table <- cell.counts

for(i in 2:length(fnames)){
  load(fnames[i])
  cell.counts.table <- rbind(cell.counts.table, cell.counts)
}

# This table contains all the cell counts. I want proportions of parent.
population.names <- colnames(cell.counts.table[13:ncol(cell.counts.table)])
parent.pop <- c(1,1,2,3,3,
                5,5,5,8,2,
                2,11,12,2,2,
                15,16,15,18,2,
                2,2,2,2,24)
for(i in 13:ncol(cell.counts.table)){
  cell.counts.table[,i] <- as.integer(as.character(cell.counts.table[,i]))
}
for(i in 13:37){
  new.prop.column <- cell.counts.table[,i]/cell.counts.table[,(12 + parent.pop[i -12])]
  cell.counts.table <- cbind(cell.counts.table, new.prop.column = new.prop.column)
  names(cell.counts.table)[ncol(cell.counts.table)] <- paste0(names(cell.counts.table[i]), " - fraction of ", names(cell.counts.table)[(12 + parent.pop[i -12])])
}

write.table(cell.counts.table, file = paste0('TCP_Panel2_CellCountsProps_', Sys.Date(), '.csv'), quote = F, sep = ',', row.names = F)

library('reshape')
library('ggplot2')
cell.counts.melted <- melt(cell.counts.table, id.vars = c(1:37))
p <- ggplot(data = cell.counts.melted, aes(x = Panel.Organ.Folder, y = value)) + 
  geom_point(alpha = 0.3, color = (as.integer(cell.counts.melted$Genotype == 'WT') + 1)) + theme_bw() + facet_wrap(~variable, scales = "free_y")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(filename = 'TCP_cell_proportions_vs_time.png', p, width = 11*3, height = 7*3)

# Same thing, but labelling only every second time point
my_breaks <- unique(cell.counts.melted$Panel.Organ.Folder)
my_labs <- NULL
my_labs[seq(1,length(unique(cell.counts.melted$Panel.Organ.Folder)), by = 2)] <-  unique(cell.counts.melted$Panel.Organ.Folder)[seq(1,length(unique(cell.counts.melted$Panel.Organ.Folder)), by = 2)]
my_labs[seq(2,length(unique(cell.counts.melted$Panel.Organ.Folder)), by = 2)] <- ''
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(breaks = my_breaks, labels = my_labs)
ggsave(filename = 'TCP_cell_proportions_vs_time.png', p, width = 12*3, height = 7*3)

# print out/get outliers
library('plyr')
outliers <- ldply(colnames(cell.counts.table)[40:ncol(cell.counts.table)], function(x){
  idx <- which(colnames(cell.counts.table) == x)
  median.prop <- median(as.numeric(as.character(cell.counts.table[,idx])))
  sd.prop <- sd(as.numeric(as.character(cell.counts.table[,idx])))
  outlier.FCSfiles <- which(abs(as.numeric(as.numeric(cell.counts.table[, idx])) - median.prop) > 4*sd.prop)
  temp1 <- as.numeric(cell.counts.table[, idx])[outlier.FCSfiles]
  names(temp1) <- outlier.FCSfiles
  
  temp1 <- sort(temp1, decreasing = TRUE)
  outlier.FCSfiles <- as.numeric(names(temp1))
  outlier.FCSfiles <- cell.counts.table[outlier.FCSfiles, c(1,3,5,12)]
  marker <- rep(x, nrow(outlier.FCSfiles))
  outlier.FCSfiles$marker <- marker
  # diff_from_median <- temp1 - median.prop
  # outlier.FCSfiles$diff_from_median <- diff_from_median
  # outlier.FCSfiles$sd.prop <- sd.prop
  return(outlier.FCSfiles)
})

write.table(outliers, file = 'Outlier_cellprops.csv', quote = F, sep = ',', row.names = F)

## copy and paste outlier files scatterplot to a new folder for easy inspections
currentFolder <- "/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel2/Results/Figures/ScatterPlots"
newFolder <- "/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel2/Results/Figures/outlierPlots"

listFiles <- unique(paste0(currentFolder,"/", gsub('.{3}$','png',paste0(outliers$Panel.Organ.Folder,"_",outliers$FCS.files))))
file.copy(listFiles, newFolder)


# Output files which flowCut flagged
flowCut.flagged.idx <- which(cell.counts.table[,"Passed flowCut"] == "F")
flowCut.flagged.files <-cell.counts.table[flowCut.flagged.idx, c(1,3,5,12)]
write.table(flowCut.flagged.files, file = 'flowCut_flagged_files.csv', quote = F, sep = ',', row.names = F)





