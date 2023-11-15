## Date: March 24, 2021
## Created by Albina Rahim
## Script for converting csv spreadsheets to xml format

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
#library("ggplot2")
library("XML")

## This part is for comparing the Proportions distributions between WT versus KO
results.dir <- paste0("/home/rstudio/results/IMPC/", toupper(centre), "/Panel", panel, "/Results")

## Version 1
CSVfile <- read.csv(paste0(results.dir,"/DCCResults_Proportions_UCD_Panel1_20191105_1819.csv"))
csv_df <- as.data.frame(CSVfile)

# create a new xml doc
doc_xml <- newXMLDoc(isHTML = FALSE)

# create a table node
table_node <- newXMLNode("table", doc = doc_xml)

# row data
row_data <- apply(csv_df, 1, function(x) {
  z1 <- newXMLNode('row') # create a new node for each row
  addChildren(z1, lapply(names(x), function(y) newXMLNode(y, x[y])))
})


# add row data to table node
xmlParent(row_data) <- table_node

saveXML(doc_xml, file = paste0(results.dir,"/csv_df.xml"))


#####################################################################

# Version 2
myData <- read.table(text = "Sales       PctSales
Id1     12929.63      0.12278547
Id2     90063.39      0.85528156
Id3      2309.60      0.02193298", header = TRUE, row.names = 1
                     , stringsAsFactors = FALSE)
myData$id <- rownames(myData)
names(myData) <- c("id", "sales", "pctSales")
con <- xmlOutputDOM("salesReport")
for(i in seq(nrow(myData))){
  con$addTag("employee", attrs = myData[i,])
}


con <- xmlOutputDOM("Results")
for(i in seq(nrow(csv_df))){
  con$addTag("FCSfiles", attrs = csv_df[i,])
}
cat(saveXML(con$value()))

saveXML(con$value(), file = paste0(results.dir,"/con_df.xml"))
