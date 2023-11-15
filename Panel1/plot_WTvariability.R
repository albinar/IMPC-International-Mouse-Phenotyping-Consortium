## Developed by Albina Rahim for Panel1 files

rm(list = ls())

library('plyr')
library('ggplot2')


# spreadsheets with event count data
fname.BCM <- '/mnt/f/FCS data/IMPC/IMPC-Results/BCM/Panel1/Results/DCCResults_BCM_Panel1_20170503_1530.csv'
fname.TCP <- '/mnt/f/FCS data/IMPC/IMPC-Results/TCP/Panel1/Results/DCCResults_TCP_Panel1_20170502_1215.csv'
fname.CIPHE <- '/mnt/f/FCS data/IMPC/IMPC-Results/CIPHE/Panel1/Results/DCCResults_WT_CIPHE_Panel1_20170505_2252.csv'

fnames <- c(TCP = fname.TCP, CIPHE = fname.CIPHE, BCM = fname.BCM)

# load IMPReSS ID info 
load('/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/immpressIDconversion.Rdata')
colnames(immpressIDconversion) <- c("IMPReSS.id", "Genotype")

df <- ldply(fnames, function(x){
  df.temp <- read.csv(x, check.names = F, stringsAsFactors = F)
  return(df.temp)
})
WT.idx <- which(df$Genotype %in% c('CIPHE', 'B6', 'CTRL', '+_+', 'Y_+'))
WT.idx <- union(grep('^WT', df$Genotype), WT.idx)
WT.idx <- union(grep('CRL', df$Genotype), WT.idx)
df <- df[WT.idx, ]
df <- df[-which(df$Gender == 'N/A'), ]
df$Gender <- gsub('\\s+','', df$Gender)
df$Gender <- gsub('Female','F', df$Gender)
df$Gender <- gsub('Male','M', df$Gender)
df$Assay.Date <- as.Date(df$`Assay Date`)

# Normalize populations to their parent
dftemp <- df
dftemp[, 'impc_imm_003_001'] <- 100^2*df[, 'impc_imm_003_001']/(df[, 'impc_imm_002_001']*df[, 'impc_imm_026_001'])
dftemp[, 'impc_imm_004_001'] <- 100^2*df[, 'impc_imm_004_001']/(df[, 'impc_imm_002_001']*df[, 'impc_imm_026_001'])
dftemp[, 'impc_imm_005_001'] <- 100^2*df[, 'impc_imm_005_001']/(df[, 'impc_imm_002_001']*df[, 'impc_imm_026_001'])
dftemp[, 'impc_imm_006_001'] <- 100^2*df[, 'impc_imm_006_001']/(df[, 'impc_imm_002_001']*df[, 'impc_imm_026_001'])

dftemp[, 'impc_imm_007_001'] <- 100*df[, 'impc_imm_007_001']/(df[, 'impc_imm_003_001'])
dftemp[, 'impc_imm_008_001'] <- 100*df[, 'impc_imm_008_001']/(df[, 'impc_imm_003_001'])
dftemp[, 'impc_imm_009_001'] <- 100*df[, 'impc_imm_009_001']/(df[, 'impc_imm_003_001'])
dftemp[, 'impc_imm_010_001'] <- 100*df[, 'impc_imm_010_001']/(df[, 'impc_imm_003_001'])

dftemp[, 'impc_imm_011_001'] <- 100*df[, 'impc_imm_011_001']/(df[, 'impc_imm_004_001'])
dftemp[, 'impc_imm_012_001'] <- 100*df[, 'impc_imm_012_001']/(df[, 'impc_imm_004_001'])
dftemp[, 'impc_imm_013_001'] <- 100*df[, 'impc_imm_013_001']/(df[, 'impc_imm_004_001'])

dftemp[, 'impc_imm_014_001'] <- 100*df[, 'impc_imm_014_001']/(df[, 'impc_imm_007_001'])
dftemp[, 'impc_imm_015_001'] <- 100*df[, 'impc_imm_015_001']/(df[, 'impc_imm_007_001'])
dftemp[, 'impc_imm_028_001'] <- 100*df[, 'impc_imm_028_001']/(df[, 'impc_imm_007_001'])
dftemp[, 'impc_imm_029_001'] <- 100*df[, 'impc_imm_029_001']/(df[, 'impc_imm_007_001'])
dftemp[, 'impc_imm_030_001'] <- 100*df[, 'impc_imm_030_001']/(df[, 'impc_imm_007_001'])
dftemp[, 'impc_imm_031_001'] <- 100*df[, 'impc_imm_031_001']/(df[, 'impc_imm_007_001'])

dftemp[, 'impc_imm_016_001'] <- 100*df[, 'impc_imm_016_001']/(df[, 'impc_imm_008_001'])
dftemp[, 'impc_imm_017_001'] <- 100*df[, 'impc_imm_017_001']/(df[, 'impc_imm_008_001'])
dftemp[, 'impc_imm_032_001'] <- 100*df[, 'impc_imm_032_001']/(df[, 'impc_imm_008_001'])
dftemp[, 'impc_imm_033_001'] <- 100*df[, 'impc_imm_033_001']/(df[, 'impc_imm_008_001'])
dftemp[, 'impc_imm_034_001'] <- 100*df[, 'impc_imm_034_001']/(df[, 'impc_imm_008_001'])
dftemp[, 'impc_imm_035_001'] <- 100*df[, 'impc_imm_035_001']/(df[, 'impc_imm_008_001'])

dftemp[, 'impc_imm_018_001'] <- 100*df[, 'impc_imm_018_001']/(df[, 'impc_imm_009_001'])
dftemp[, 'impc_imm_019_001'] <- 100*df[, 'impc_imm_019_001']/(df[, 'impc_imm_009_001'])
dftemp[, 'impc_imm_036_001'] <- 100*df[, 'impc_imm_036_001']/(df[, 'impc_imm_009_001'])
dftemp[, 'impc_imm_037_001'] <- 100*df[, 'impc_imm_037_001']/(df[, 'impc_imm_009_001'])
dftemp[, 'impc_imm_038_001'] <- 100*df[, 'impc_imm_038_001']/(df[, 'impc_imm_009_001'])
dftemp[, 'impc_imm_039_001'] <- 100*df[, 'impc_imm_039_001']/(df[, 'impc_imm_009_001'])

dftemp[, 'impc_imm_020_001'] <- 100*df[, 'impc_imm_020_001']/(df[, 'impc_imm_011_001'])
dftemp[, 'impc_imm_021_001'] <- 100*df[, 'impc_imm_021_001']/(df[, 'impc_imm_011_001'])
dftemp[, 'impc_imm_040_001'] <- 100*df[, 'impc_imm_040_001']/(df[, 'impc_imm_011_001'])
dftemp[, 'impc_imm_041_001'] <- 100*df[, 'impc_imm_041_001']/(df[, 'impc_imm_011_001'])
dftemp[, 'impc_imm_042_001'] <- 100*df[, 'impc_imm_042_001']/(df[, 'impc_imm_011_001'])

dftemp[, 'impc_imm_022_001'] <- 100*df[, 'impc_imm_022_001']/(df[, 'impc_imm_012_001'])
dftemp[, 'impc_imm_023_001'] <- 100*df[, 'impc_imm_023_001']/(df[, 'impc_imm_012_001'])
dftemp[, 'impc_imm_043_001'] <- 100*df[, 'impc_imm_043_001']/(df[, 'impc_imm_012_001'])
dftemp[, 'impc_imm_044_001'] <- 100*df[, 'impc_imm_044_001']/(df[, 'impc_imm_012_001'])
dftemp[, 'impc_imm_045_001'] <- 100*df[, 'impc_imm_045_001']/(df[, 'impc_imm_012_001'])

dftemp[, 'impc_imm_024_001'] <- 100*df[, 'impc_imm_024_001']/(df[, 'impc_imm_013_001'])
dftemp[, 'impc_imm_025_001'] <- 100*df[, 'impc_imm_025_001']/(df[, 'impc_imm_013_001'])
dftemp[, 'impc_imm_046_001'] <- 100*df[, 'impc_imm_046_001']/(df[, 'impc_imm_013_001'])
dftemp[, 'impc_imm_047_001'] <- 100*df[, 'impc_imm_047_001']/(df[, 'impc_imm_013_001'])
dftemp[, 'impc_imm_048_001'] <- 100*df[, 'impc_imm_048_001']/(df[, 'impc_imm_013_001'])

df <- dftemp
rm(dftemp)

df$x <- paste(df$.id, df$Gender, sep = " ")

# Violin plots of  normalized cell counts as a function of (center, sex)

#for(i in c(13:29, 31, 33, 35:36, 38, 39)){
for(i in c(12:55)){
  current.pop <- colnames(df)[i]
  plot.title <- immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == current.pop), 2]
  temp <- df[, c("x", "Assay.Date", "Gender", ".id", "FCS files", colnames(df)[i])]
  colnames(temp)[6] <- 'Value'
  
  p <- ggplot(temp, aes(x = x, y = Value, color = Gender)) + 
    geom_violin(fill = NA) + geom_point(alpha = 0.1) + theme_bw() + scale_color_brewer(palette = 'Dark2') +
    ylab(paste0(current.pop, " (percent parent)")) + ggtitle(plot.title) + labs(color='Sex') + xlab('')
  ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/', current.pop, '.png'), p, height = 5, width = 12)
}


plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_003_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_003_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_003_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_003_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_003_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_003_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_004_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_004_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_004_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_004_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_004_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_004_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_005_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_005_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_005_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_005_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_005_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_005_001_vsTime_Female.png'), p, height = 5, width = 12)


plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_007_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_007_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_007_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_007_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_007_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_007_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_008_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_008_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_008_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_008_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_008_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_008_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_009_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_009_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_009_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_009_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_009_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_009_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_011_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_011_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_011_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_011_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_011_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_011_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_012_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_012_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_012_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_012_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_012_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_012_001_vsTime_Female.png'), p, height = 5, width = 12)



plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_013_001'), 2], " - Male")
p <- ggplot(df[which(df$Gender == 'M'), ], aes(x = Assay.Date, y = impc_imm_013_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_013_001_vsTime_Male.png'), p, height = 5, width = 12)

plot.title <- paste0(immpressIDconversion[which(immpressIDconversion[, 'IMPReSS.id'] == 'impc_imm_013_001'), 2], " - Female")
p <- ggplot(df[which(df$Gender == 'F'), ], aes(x = Assay.Date, y = impc_imm_013_001, color = .id)) + 
  geom_point(alpha = 0.8) + theme_bw() + scale_color_brewer(palette = 'Dark2') + geom_smooth(method = "loess", se = TRUE) + ggtitle(plot.title)
ggsave(paste0('/mnt/f/FCS data/IMPC/IMPC-Results/IMPC_Panel1_WT_comparison/impc_imm_013_001_vsTime_Female.png'), p, height = 5, width = 12)
