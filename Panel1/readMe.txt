Attached are the two scripts: labelFCSMain.R which is the main script that calls the function, labelFCSFunc.R. 
The function, labelFCSFunc.R does the actual work of processing the raw files (.LMD or .FCS) files.

The function, labelFCSFunc.R requires the packages:"flowCore" and "stringR". 
If the packages are not installed the program will automatically install them when running the program.

##############################################################################################

The first step is to save the labelFCSMain.R and labelFCSFunc.R in the same directory.


Next in the main script, labelFCSMain.R, change the following lines:

Line 9:
##This the working directory where the R scripts (both the Main and Function) have been saved.
## Change this to your working directory. 
setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel1")

Line 16:
## Specifying the panel number.
## Change this to panel<-"2" when working on Panel 2
panel <- "1" 

Line 20:
##Change the inputPath to the directory where the raw files have been kept for Panel1
inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel A"
  
Line 23:  
##Change the outputPath to the directory where the processed FCS files will be saved for Panel1
outputPath <- "/home/rstudio/results/IMPC/GMC/Panel1"

Line 33:
##Change the inputPath to the directory where the raw files have been kept for Panel2
inputPath <- "/home/rstudio/data/IMPC/GMC/Original_FCS/Panel B"
  
Line 36:
##Change the outputPath to the directory where the processed FCS files will be saved for Panel2
outputPath <- "/home/rstudio/results/IMPC/GMC/Panel2"

#################################################################################################

Last step: save the main script, labelFCSMain.R and now run the program.
source("labelFCSMain.R")

