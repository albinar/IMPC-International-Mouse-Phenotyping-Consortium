## Developed by Albina Rahim
## Date: December 09, 2016
## This function does the Pre-Processing of the datasets and can be applied on files of all Panels and Centres. 
## As input it requires the paths of the raw FCS files and paths of the metadata spreadsheets and the number of CPUs to use while parallelization
## It will remove FCS files:
## 1. with no Barcodes
## 2. whose Barcodes are not listed in the metadata spreadsheets
## 3. which are Corrupted
## 4. with less than 20,000 cells
## 5. which are Duplicates
## As output it returns a large matrix- store.allFCS, which contains information of the paths, Panel/Organ,
## Genotype of the file,  Names of the FCS files, Label Barcode, Assay Date, Gender, Number of Channels, and Number of Cells.
## It also returns a long list of Genotypes, unique Genotypes, details of FCS files which were sent to us but 
## whose information was missing in the spreadsheet, information of those Barcodes in the spreadsheet for which the corresponding 
## FCS files were not given to us for automated analysis, information on all the Corrupted FCS files so that we can 
## send them to Centre for resending the files through flowRepository, details of files which has less than 20,000 cells,
## and information to determine if all the files have the same number of channels

#############################################################################################################

## preProcessingFunc <- function(inputPath, inputCSV, barcodeExpression){
##            ....
## }
## inputPath is a list which contains paths of all the datasets sent at various times
## inputCSV contains information for all the metadata spreadsheets sent at various times
## barcodeExpression to determine the barcode expression, which are unique for each centre and helps us to compare files against the metadata spreadsheet

##############################################################################################################

preProcessingFunc <- function(inputPath, inputCSV, centre, barcodeExpression){
  library("flowCore")
  library("stringr")
  
  
  ## The following lines are added because the difference spreadsheets of different centres have different column names for Barcodes, Assay Dates, Gender, and Genotype
  barcode.colNames <- c("Label.Barcode", "ID", "Test.Code", "Mouse_Number")
  mouseID.colNames <- c("Mouse_ID", "Experiment_Number", "Animal.ID", "Experiment_ID")
  assayDate.colNames <- c("Assay.Date", "Experiment_Date", "Date.of.Sample.Prep", "Date")
  gender.colNames <- c("Gender", "Sex")
  genotype.colNames <- c("Genotype", "Group", "Gene_ID", "Genotype_Gene")
  colonyID.colNames <- c("Colony_ID", "Line", "Run_Number")
  
  store.allFCS <- NULL
  store.allFMO <- NULL
  NObarcodes.FCS <- 0
  notListed.FCS <- NULL
  corrupted.FCS <- NULL
  lessCells.FCS <- NULL
  duplicate.FCS.temp <- NULL
  Barcodes.NoFCS.temp <- NULL
  Barcodes.NoFCS <- NULL
  Genotype <- c()
  Mouse_Label <- c()
  controlNames <- c("CTRL", "FMO", "UNSTAINEDUNSTAINED", "UNSTAINED", "Beads", "Neg", "DAPI", "SingleStaincells",
                    "_Nikolai", "WT", "WT1", "WT2", "WT3", "WT4", "KO1", "KO2", "KO3", "KO4", "Compensation Controls", "control", "wild type" ) ## Adding the control names since most of centres except for Sanger have control files together with main FCS files
  
  
  
  numPaths <- length(inputPath) # Length of inputPath to determine the number of paths from where the files need to be retrieved
  
  for(i in 1:numPaths){
    # Path to the FCS files
    pathFCS <- unlist(inputPath[i]) 
    # Reads all folders and files in current path folder and makes a list of all of their paths
    allFCS <- dir(pathFCS, full.names=T, recursive=T, pattern = c("*.fcs")) 
    
    
    store.allFCS.temp <- sapply(1:length(allFCS), function(x){pathFCS})
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-1]}))
    if(centre == "bcm"){
      tempCol <- sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-2]})
      store.allFCS.temp[,2] <- paste0(tempCol, "/", store.allFCS.temp[,2])
    }
    store.allFCS.temp <- cbind(store.allFCS.temp, NA) # Column for Genotype
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))]}))
    
    ## Storing information about the FMO controls saving them in a separate matrix
    # if(centre == "gmc"){
    #   allFMO <-  dir(pathFCS, full.names=T, recursive=T, pattern = c("*.FCS")) 
    #   store.allFMO.temp <- sapply(1:length(allFMO), function(x){paste0(pathFCS, "/",unlist(strsplit(allFMO[x], split = "/"))[length(unlist(strsplit(allFMO[x], split = "/")))-1])})
    #   store.allFMO.temp <- cbind(store.allFMO.temp, sapply(1:length(allFMO), function(x){unlist(strsplit(allFMO[x], split = "/"))[length(unlist(strsplit(allFMO[x], split = "/")))]}))
    #   colnames(store.allFMO.temp) <- c("Path", "FMO")
    # }else{
    
    if(centre == "ucd"){
      matches <- grep(paste("FMO",collapse="|"), store.allFCS.temp[,4], value=FALSE)
    }else{
      matches <- grep(paste(controlNames,collapse="|"), store.allFCS.temp[,4], value=FALSE)
    }
    
    if(length(matches) != 0){
      store.allFMO.temp <- store.allFCS.temp[matches,]
      store.allFMO.temp <- store.allFMO.temp[,-3]
      store.allFCS.temp <- store.allFCS.temp[-matches,]
      rownames(store.allFCS.temp) <- 1:nrow(store.allFCS.temp)
    }else{
      store.allFMO.temp <- NA
    }
    #}
    
    
    if(centre == "tcp"| centre == "ucd"){
      FMO.index <- grep("FMO", store.allFMO.temp[,3], value = FALSE)
      store.allFMO.temp <- store.allFMO.temp[FMO.index,]
      if(centre == "tcp"){
        colnames(store.allFMO.temp) <- c("Path", "Assay Date", "FMO")
      }else if(centre == "ucd"){
        colnames(store.allFMO.temp) <- c("Path", "Panel/Organ/Folder", "FMO")
      }
      
    }
    
    
    store.allFCS.temp <- cbind(store.allFCS.temp, str_extract(store.allFCS.temp[,4], barcodeExpression))
    
    if(centre == "ucd" & panel == 1){
      tempA <- str_replace_all(store.allFCS.temp[,5],fixed(c("LWT_"="", "WT_"="", "KO_"="" ,"_CORRECT_PANEL_A"="", "_PANEL_A"="", "_PANEL"="")))
      store.allFCS.temp[,5] <- tempA
      
    }else if(centre == "ucd" & panel == 2){
      tempA <- str_replace_all(store.allFCS.temp[,5],fixed(c("LWT_"="", "WT_"="", "KO_"="" ,"_CORRECT_PANEL_B"="", "_PANEL_B"="", "_PANEL"="")))
      store.allFCS.temp[,5] <- tempA
    }else if(centre == "ccp"){
      tempA <- gsub('_',"", store.allFCS.temp[,5])
      store.allFCS.temp[,5] <- tempA
    }
    store.allFCS.temp <- cbind(store.allFCS.temp, NA) # Column for Mouse ID
    store.allFCS.temp <- cbind(store.allFCS.temp, NA) # Column for Colony ID
    store.allFCS.temp <- cbind(store.allFCS.temp, NA) # Column for Assay date
    store.allFCS.temp <- cbind(store.allFCS.temp, NA) # Column for Gender
    
    ###########################################################################################################
    # 1. Checking for files with NO Barcodes and remove them
    
    if(centre != "ciphe"){ ## Only CIPHE doesn't have Barcode in their file names
      index.Remove <- which(is.na(store.allFCS.temp[,5]))
      
      if(length(index.Remove) != 0){
        NObarcodes.FCS <- NObarcodes.FCS+length(index.Remove)
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
    }
    
    colnames(store.allFCS.temp) <- c("Path", "Panel/Organ/Folder", "Genotype", "FCS files", "Strain Code/Ear Tag/Mouse Number", "Mouse ID/Animal ID", "Colony ID/Line/Run Number", "Assay Date", "Gender")
    
    ################################################################################
    # 2. Files whose Barcodes are not listed in the metadata spreadsheets
    ## Reading the metadata spreadsheet
    CSVfile <- read.csv(unlist(inputCSV[i]))
    if(centre == "tcp" | centre == "ciphe"){
      Label.Barcode <- paste0(CSVfile[,c('Strain_Code')],"_",CSVfile[,c('Ear_Tag')])
      CSVfile <- cbind(Label.Barcode, CSVfile)
    }
    
    
    CSVfile <- as.matrix(CSVfile)
    if(centre == "ucd"){
      Genotype_Gene <- paste0(CSVfile[,c('Genotype')],"_",CSVfile[,c('Gene')])
      CSVfile <- cbind(CSVfile, Genotype_Gene)
      tempB <- str_replace_all(CSVfile[,c("Mouse_Number")],fixed("-"), "_")
      CSVfile[,c("Mouse_Number")] <- tempB
      CSVfile <- CSVfile[,-6:-7]
    }
    
    # write.csv(CSVfile, file = "/home/rstudio/data/IMPC/GMC/GMC_2019-08-28/Panel_A/Metadata_GMC_2013-2018_PanelA.csv", row.names = FALSE)
    
    ## Retrieving Genotype information from the spreadsheet
    Genotype.temp <- CSVfile[,match(genotype.colNames,colnames(CSVfile))]
    Genotype.temp <- Genotype.temp[!is.na(Genotype.temp)]
    #Genotype.temp <- Genotype.temp[,1]
    if(centre == "ccp"){
      Genotype.temp <- sub("control","WT", Genotype.temp)
    }else{
      Genotype.temp <- sub("/","_", Genotype.temp)
    }
    
    
    ## Retrieving Barcode information from the spreadsheet
    Mouse_Label.temp <- CSVfile[,match(barcode.colNames,colnames(CSVfile))]
    #Mouse_Label.temp <- Mouse_Label.temp[,1]
    Mouse_Label.temp <- Mouse_Label.temp[!is.na(Mouse_Label.temp)]
    
    ## Retrieving Mouse ID/Animal ID information from the spreadsheet
    #if(centre != "gmc"){
    #if(centre != "ucd"){
    Mouse_ID.temp <- CSVfile[,match(mouseID.colNames,colnames(CSVfile))]
    Mouse_ID.temp <- Mouse_ID.temp[!is.na(Mouse_ID.temp)]
    #Mouse_ID.temp <- Mouse_ID.temp[,1]
    #}
    
    
    ## Retrieving Colony ID/Line information from the spreadsheet
    Colony_ID.temp <- CSVfile[,match(colonyID.colNames,colnames(CSVfile))]
    if(centre == "tcp" | centre == "jax" | centre == "ccp"){
      #Colony_ID.temp <- Colony_ID.temp[!is.na(Colony_ID.temp)]
      Colony_ID.temp <- Colony_ID.temp[,1]
    }else if(centre == "ucd"){
      Colony_ID.temp <- Colony_ID.temp[,3]
    }
    
    
    ## Retrieving Assay Date information from the spreadsheet
    Assay_Date <- CSVfile[,match(assayDate.colNames,colnames(CSVfile))]
    #Assay_Date <- Assay_Date[,2]
    Assay_Date <- Assay_Date[!is.na(Assay_Date)]
    
    ## Retrieving Gender information from the spreadsheet
    Gender <- CSVfile[,match(gender.colNames,colnames(CSVfile))]
    #Gender <- Gender[,2] 
    Gender <- Gender[!is.na(Gender)]
    
    rownames(store.allFCS.temp) <- 1:nrow(store.allFCS.temp)
    if(centre == "sanger" | centre == "tcp" | centre == "jax" | centre == "gmc"|centre == "ucd" | centre == "ccp"){
      ## Checking for files whose Barcodes are not listed in the spreadsheet
      countX <-0
      index.Remove <- 0
      print("Start checking for files whose Barcodes are NOT listed in the spreadsheet")  
      for(x in 1:nrow(store.allFCS.temp)){
        if(centre == "gmc" & panel == 1){
          temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_A')])
        }else if(centre == "gmc" & panel == 2){
          temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_B')])
        }else{
          temp <- grep(store.allFCS.temp[x,5], Mouse_Label.temp)
        }
        
        if(length(temp) == 1){
          store.allFCS.temp[x,c('Genotype')] <- Genotype.temp[temp]
          store.allFCS.temp[x,c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp]
          store.allFCS.temp[x,c('Colony ID/Line/Run Number')] <- Colony_ID.temp[temp]
          store.allFCS.temp[x,c('Assay Date')] <- Assay_Date[temp]
          store.allFCS.temp[x,c('Gender')] <- Gender[temp]
        }else if(length(temp) > 1){
          for(y in 1:length(temp)){
            store.allFCS.temp[temp[y],c('Genotype')] <- Genotype.temp[temp[y]]
            store.allFCS.temp[temp[y],c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp[y]]
            store.allFCS.temp[temp[y],c('Colony ID/Line/Run Number')] <- Colony_ID.temp[temp[y]]
            store.allFCS.temp[temp[y],c('Assay Date')] <- Assay_Date[temp[y]]
            store.allFCS.temp[temp[y],c('Gender')] <- Gender[temp[y]]
          }
        }else{
          index.Remove <- c(index.Remove,x)
          countX <- countX+1
        }
      } ## end of for loop
      
      index.Remove <-index.Remove[index.Remove !=0]
      
      ## Storing information of all the FCS files which were sent to us but whose information were not listed in the metadata spreadsheet
      # Removing information of the FCS files which were sent to use but whose information were not listed in the metadata spreadsheet from the main storage matrix
      if(length(index.Remove) != 0){
        notListed.FCS <- rbind(notListed.FCS, store.allFCS.temp[index.Remove,])
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
    }else if(centre == "ciphe" & panel == 1){
      countX <-0
      index.Remove <- 0
      print("Start checking for files whose Barcodes are NOT listed in the spreadsheet")  
      for(x in 1:nrow(store.allFCS.temp)){
        temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_A')])
        if(length(temp) == 1){
          store.allFCS.temp[x,c('Strain Code/Ear Tag/Mouse Number')] <- Label.Barcode[temp]
          store.allFCS.temp[x,c('Genotype')] <- Genotype.temp[temp]
          store.allFCS.temp[x,c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp]
          store.allFCS.temp[x,c('Colony ID/Line/Run Number')] <- Colony_ID.temp[temp]
          store.allFCS.temp[x,c('Assay Date')] <- Assay_Date[temp]
          store.allFCS.temp[x,c('Gender')] <- Gender[temp]
        }else{
          index.Remove <- c(index.Remove,x)
          countX <- countX+1
        }
      }
      index.Remove <-index.Remove[index.Remove !=0]
      ## Storing information of all the FCS files which were sent to us but whose information were not listed in the metadata spreadsheet
      # Removing information of the FCS files which were sent to use but whose information were not listed in the metadata spreadsheet from the main storage matrix
      if(length(index.Remove) != 0){
        notListed.FCS <- rbind(notListed.FCS, store.allFCS.temp[index.Remove,])
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
      
      # index <- sapply(1:nrow(store.allFCS.temp), function(j){which(store.allFCS.temp[j,c('FCS files')] == CSVfile[,c('FCS_Files_A')])})
      # store.allFCS.temp[,c('Strain Code/Ear Tag/Mouse Number')] <- Label.Barcode[index]
      # store.allFCS.temp[,c('Mouse ID')] <- Mouse_ID.temp[index]
      # store.allFCS.temp[,c('Colony ID')] <- Colony_ID.temp[index]
      # store.allFCS.temp[,c('Genotype')] <- Genotype.temp[index]
      # store.allFCS.temp[,c('Assay Date')] <- Assay_Date[index]
      # store.allFCS.temp[,c('Gender')] <- Gender[index]
    }else if(centre == "ciphe" & panel == 2){
      countX <-0
      index.Remove <- 0
      print("Start checking for files whose Barcodes are NOT listed in the spreadsheet")  
      for(x in 1:nrow(store.allFCS.temp)){
        temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_B')])
        if(length(temp) == 1){
          store.allFCS.temp[x,c('Strain Code/Ear Tag/Mouse Number')] <- Label.Barcode[temp]
          store.allFCS.temp[x,c('Genotype')] <- Genotype.temp[temp]
          store.allFCS.temp[x,c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp]
          store.allFCS.temp[x,c('Colony ID/Line/Run Number')] <- Colony_ID.temp[temp]
          store.allFCS.temp[x,c('Assay Date')] <- Assay_Date[temp]
          store.allFCS.temp[x,c('Gender')] <- Gender[temp]
        }else{
          index.Remove <- c(index.Remove,x)
          countX <- countX+1
        }
      }
      index.Remove <-index.Remove[index.Remove !=0]
      ## Storing information of all the FCS files which were sent to us but whose information were not listed in the metadata spreadsheet
      # Removing information of the FCS files which were sent to use but whose information were not listed in the metadata spreadsheet from the main storage matrix
      if(length(index.Remove) != 0){
        notListed.FCS <- rbind(notListed.FCS, store.allFCS.temp[index.Remove,])
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
    }else if(centre == "bcm" & panel == 1){
      countX <-0
      index.Remove <- 0
      print("Start checking for files whose Barcodes are NOT listed in the spreadsheet")  
      for(x in 1:nrow(store.allFCS.temp)){
        temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_A')])
        if(length(temp) == 1){
          store.allFCS.temp[x,c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp]
          store.allFCS.temp[x,c('Colony ID/Line')] <- Colony_ID.temp[temp]
          store.allFCS.temp[x,c('Genotype')] <- Genotype.temp[temp]
          store.allFCS.temp[x,c('Assay Date')] <- Assay_Date[temp]
          store.allFCS.temp[x,c('Gender')] <- Gender[temp]
        }else{
          index.Remove <- c(index.Remove,x)
          countX <- countX+1
        }
      }
      index.Remove <-index.Remove[index.Remove !=0]
      
      ## Storing information of all the FCS files which were sent to us but whose information were not listed in the metadata spreadsheet
      # Removing information of the FCS files which were sent to use but whose information were not listed in the metadata spreadsheet from the main storage matrix
      if(length(index.Remove) != 0){
        notListed.allFCS <- store.allFCS.temp[index.Remove,]
        p1.index.Remove <- grep("Panel1", notListed.allFCS[,c('FCS files')])
        notListed.FCS <- rbind(notListed.FCS, notListed.allFCS[p1.index.Remove,])
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
    }else if(centre == "bcm" & panel == 2){
      countX <-0
      index.Remove <- 0
      print("Start checking for files whose Barcodes are NOT listed in the spreadsheet")  
      for(x in 1:nrow(store.allFCS.temp)){
        temp <- grep(store.allFCS.temp[x,c('FCS files')], CSVfile[,c('FCS_Files_B')])
        if(length(temp) == 1){
          store.allFCS.temp[x,c('Mouse ID/Animal ID')] <- Mouse_ID.temp[temp]
          store.allFCS.temp[x,c('Colony ID/Line')] <- Colony_ID.temp[temp]
          store.allFCS.temp[x,c('Genotype')] <- Genotype.temp[temp]
          store.allFCS.temp[x,c('Assay Date')] <- Assay_Date[temp]
          store.allFCS.temp[x,c('Gender')] <- Gender[temp]
        }else{
          index.Remove <- c(index.Remove,x)
          countX <- countX+1
        }
      }
      index.Remove <-index.Remove[index.Remove !=0]
      
      ## Storing information of all the FCS files which were sent to us but whose information were not listed in the metadata spreadsheet
      # Removing information of the FCS files which were sent to use but whose information were not listed in the metadata spreadsheet from the main storage matrix
      if(length(index.Remove) != 0){
        notListed.allFCS <- store.allFCS.temp[index.Remove,]
        p2.index.Remove <- grep("Panel2", notListed.allFCS[,c('FCS files')])
        notListed.FCS <- rbind(notListed.FCS, notListed.allFCS[p2.index.Remove,])
        store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
      }
    }
    
    
    ########################################################################################################
    
    ## 3. Code for removing duplicate FCS files based on their Barcodes 
    if(centre == "sanger"){
      Barcodes <- str_extract(store.allFCS.temp[,c('Strain Code/Ear Tag/Mouse Number')],barcodeExpression)
      duplicate.index <- which(duplicated(Barcodes)==TRUE)
      duplicate.FCS.temp <- rbind(duplicate.FCS.temp, store.allFCS.temp[duplicate.index,])
      store.allFCS.temp <- store.allFCS.temp[!duplicated(store.allFCS.temp[,c('Strain Code/Ear Tag/Mouse Number')]),]
    }else{
      duplicate.index <- which(duplicated(store.allFCS.temp[,c('FCS files')])==TRUE)
      duplicate.FCS.temp <- rbind(duplicate.FCS.temp, store.allFCS.temp[duplicate.index,])
      store.allFCS.temp <- store.allFCS.temp[!duplicated(store.allFCS.temp[,c('FCS files')]),]
    }
    
    
    ###############################################################################################################
    
    ## 4. This part of the script extracts those Barcodes in the spreadsheet for which the corresponding FCS files were not given to us for automated analysis
    if((centre == "ciphe" | centre == "bcm" | centre == "gmc") & panel == 1){
      allBarcodes.metadata <- CSVfile[,c('FCS_Files_A')]  
      Barcodes.NoFCS.temp <- c(Barcodes.NoFCS.temp, setdiff(allBarcodes.metadata, store.allFCS.temp[,c('FCS files')]))
    }else if((centre == "ciphe" | centre == "bcm") & panel == 2){
      allBarcodes.metadata <- CSVfile[,c('FCS_Files_B')]  
      Barcodes.NoFCS.temp <- c(Barcodes.NoFCS.temp, setdiff(allBarcodes.metadata, store.allFCS.temp[,c('FCS files')]))
    }else{
      allBarcodes.metadata <- unlist(strsplit(Mouse_Label.temp, ",", perl = TRUE))
      # To determine if there are any duplicate Barcode entries in the metadata spreadsheet
      duplicate.index <- which(duplicated(allBarcodes.metadata)==TRUE)
      if(length(duplicate.index) != 0){
        allBarcodes.metadata <- allBarcodes.metadata[-duplicate.index]
      }
      
      Barcodes.NoFCS.temp <- c(Barcodes.NoFCS.temp, setdiff(allBarcodes.metadata, store.allFCS.temp[,c('Strain Code/Ear Tag/Mouse Number')]))
    }
    
    ################################################################################################################
    # 5. Checking for files which are Corrupted and removing them. We parallelize this part of the code using ddply()
    # We also record the number of Channels for all the files and the number of Cells.
    
    print("Start finding the Corrupted Files")
    file.names <- data.frame(store.allFCS.temp, stringsAsFactors = F)
    
    corrupted.cols.cells <- ddply(file.names, "FCS.files", function(x){
      set.seed(100)
      index.Corrupted <- matrix(nrow = 1, ncol = 1, data = NA)
      NumberOfCols <- matrix(nrow = 1, ncol = 1, data = NA)
      NumberOfCells <- matrix(nrow = 1, ncol = 1, data = NA)
      #x<- file.names[i,]
      
      if(centre == "sanger" | centre == "ciphe"){
        f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
      }else if(centre == "tcp" | centre == "bcm" | centre == "jax" | centre == "gmc"| centre == "ucd" | centre == "ccp"){
        f <- try(read.FCS(filename = paste0(x$Path, "/", x$Panel.Organ.Folder, "/", x$FCS.files)), silent = TRUE)
      }
      
      if(class(f)=="try-error"){
        index.Corrupted[1] <- "Corrupted"
      }else{
        NumberOfCols[1] <- ncol(f@exprs)
        NumberOfCells[1] <- nrow(f@exprs)
      }
      data.frame(index.Corrupted, NumberOfCols, NumberOfCells)
    }, .parallel = TRUE) # end ddply
    
    ## Arranging the output from ddply() in order of the file.names
    corrupted.cols.cells <- join(file.names, corrupted.cols.cells)
    
    ## Combinging the Number of Channels and Number of Cells with the temporary storage matrix.
    store.allFCS.temp <- cbind(store.allFCS.temp, corrupted.cols.cells[,c('NumberOfCols')], corrupted.cols.cells[,c('NumberOfCells')])
    colnames(store.allFCS.temp) <- c("Path", "Panel/Organ/Folder", "Genotype", "FCS files", "Strain Code/Ear Tag/Mouse Number", "Mouse ID/Animal ID", "Colony ID/Line/Run Number", "Assay Date", "Gender",  "Number of Channels", "Number of Cells")
    
    # Locating the indices of the Corrupted files
    index.Corrupted <- which(!is.na(corrupted.cols.cells[,c('index.Corrupted')]))
    
    
    ## Storing the information for the Corrupted files, so we can send the information to Centre for resending these files through flowRepository
    ## Removing the Corrupted FCS files from store.allFCS.temp
    if(length(index.Corrupted) != 0){
      corrupted.FCS <- rbind(corrupted.FCS, store.allFCS.temp[index.Corrupted,])
      store.allFCS.temp <- store.allFCS.temp[-index.Corrupted,]
    }
    
    print("End of finding the Corrupted Files and removing them from the stored matrix")
    
    
    ##########################################################################################################
    ## 6. Checking for files which has less than 20,000 cells and storing information for such files separately and then removing them from the main storage matrix
    index.lessCells <- 0
    index.lessCells <- which(as.numeric(store.allFCS.temp[,c('Number of Cells')]) < 20000)
    index.lessCells <- index.lessCells[index.lessCells !=0]
    
    if(length(index.lessCells) != 0){
      lessCells.FCS <- rbind(lessCells.FCS, store.allFCS.temp[index.lessCells,])
      store.allFCS.temp <- store.allFCS.temp[-index.lessCells,]
    }
    
    
    ########################################################################################################
    
    
    Genotype <- c(Genotype, Genotype.temp) # Combining all the Genotypes together
    Mouse_Label <- c(Mouse_Label, Mouse_Label.temp) # Combining all the Barcodes together
    
    # if(centre == "bcm"){
    #   index.Remove <- which(is.na(store.allFCS.temp[,c("Assay Date")]))
    #   if(length(index.Remove) != 0){
    #     store.allFCS.temp <- store.allFCS.temp[-index.Remove,]
    #   }
    # }
    ## Combining the Paths, Name of the Files, Genotypes, Barcodes, Assay Dates, Gender, Number of Channels, Number of Cells with each FCS file in a large matrix store.allFCS
    store.allFCS <- rbind(store.allFCS, store.allFCS.temp) 
    
    ## Combining the Paths, Assay Dates, FMO with each FMO file in a large matrix store.allFMO
    if(is.na(store.allFMO.temp)){
      store.allFMO <- store.allFMO
    }else{
      store.allFMO <- rbind(store.allFMO, store.allFMO.temp)
    }
    
    
  } # end of outer for-loop
  
  colnames(store.allFCS) <- c("Path", "Panel/Organ/Folder", "Genotype", "FCS files", "Strain Code/Ear Tag/Mouse Number", "Mouse ID/Animal ID", "Colony ID/Line/Run Number", "Assay Date", "Gender",  "Number of Channels", "Number of Cells")
  
  if(centre == "gmc"| centre == "ucd"){
    colnames(store.allFMO) <- c("Path", "Panel/Organ/Folder","FMO")
  }else if (centre == "ciphe" | centre == "ccp"){
    if(is.null(store.allFMO)){
      store.allFMO <- NA
    }
  }else{
    colnames(store.allFMO) <- c("Path", "Assay Date", "FMO")
  }
  
  
  ########################################################################################################
  
  ## 7. Code for removing duplicate FCS files based on their barcodes, if there are duplicates between the different batches of dataset that were sent to us
  if(centre == "sanger"){
    Barcodes <- str_extract(store.allFCS[,c('Strain Code/Ear Tag/Mouse Number')],"L[0-9]+")
    duplicate.index <- which(duplicated(Barcodes)==TRUE)
    if(length(duplicate.index) != 0){
      duplicate.FCS <- rbind(duplicate.FCS.temp, store.allFCS[duplicate.index,c('Path', 'Panel/Organ/Folder', 'Genotype', 'FCS files', 'Strain Code/Ear Tag/Mouse Number', 'Mouse ID/Animal ID', 'Colony ID/Line', 'Assay Date', 'Gender')])
      store.allFCS <- store.allFCS[!duplicated(store.allFCS[,c('Strain Code/Ear Tag/Mouse Number')]),]
      rownames(store.allFCS) <- 1:length(store.allFCS[,1])
    }else{
      duplicate.FCS <- duplicate.FCS.temp
    }
  }else{
    duplicate.index <- which(duplicated(store.allFCS[,c('FCS files')])==TRUE)
    if(length(duplicate.index) != 0){
      duplicate.FCS <- rbind(duplicate.FCS.temp, store.allFCS[duplicate.index,c('Path', 'Panel/Organ/Folder', 'Genotype', 'FCS files', 'Strain Code/Ear Tag/Mouse Number', 'Mouse ID/Animal ID', 'Colony ID/Line/Run Number', 'Assay Date', 'Gender')])
      store.allFCS <- store.allFCS[!duplicated(store.allFCS[,c('FCS files')]),]
      rownames(store.allFCS) <- 1:length(store.allFCS[,1])
    }else{
      duplicate.FCS <- duplicate.FCS.temp
    }
  }
  
  
  ## Finding duplicates among the Barcodes listed/FCS files(for CIPHE only) for ALL the matadata spreadsheet but whose FCS files were not sent to us
  duplicate.index <- which(duplicated(Barcodes.NoFCS.temp)==TRUE)
  if(length(duplicate.index) != 0){
    Barcodes.NoFCS.temp <- Barcodes.NoFCS.temp[-duplicate.index]
  }
  ## Cross checking between the different batches to determine the Barcodes listed in ALL the metadata spreadsheet but whose FCS files were not sent to us
  if(centre == "ciphe"){
    Barcodes.NoFCS <- setdiff(Barcodes.NoFCS.temp, store.allFCS[,c('FCS files')])
  }else{
    Barcodes.NoFCS <- setdiff(Barcodes.NoFCS.temp, store.allFCS[,c('Strain Code/Ear Tag/Mouse Number')])
  }
  
  
  # Finding the unique Genotypes (KOs + WTs)
  uniqueGT <- unique(Genotype) 
  
  ## Checking if all the FCS files have the same number of Channels
  numChannels <- unique(store.allFCS[,c('Number of Channels')])
  
  return(list(store.allFCS = store.allFCS, store.allFMO = store.allFMO, NObarcodes.FCS = NObarcodes.FCS, Genotype = Genotype, uniqueGT = uniqueGT, notListed.FCS = notListed.FCS, Barcodes.NoFCS = Barcodes.NoFCS, corrupted.FCS = corrupted.FCS, lessCells.FCS = lessCells.FCS, duplicate.FCS = duplicate.FCS, numChannels = numChannels))
  
}

