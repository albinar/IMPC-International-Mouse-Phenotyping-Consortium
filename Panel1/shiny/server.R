

#Change line 21 in server.R, and line 23-24 in shinyHelper if plot and fcs file names matches.
source('shinyHelper.R')

#this section is to get the folder directory, and the reference sample path if provided
server<- function(input, output,session) {
  
  volumes = getVolumes()
  shinyDirChoose(input, "result", roots = volumes, session = session)
  folder.info <- reactiveValues() 
  
  
  values <- reactiveValues(df_data = NULL,df_cols=NULL,df_time=NULL)
  observeEvent(input$result, #priority = 2,
               #input$dir only returns a second variable if a folder is selected
               if(!is.na(input$result[2])){
                 folder.info$path <- shinyFiles::parseDirPath(volumes, input$result) #this location needs to be specified alot of times, so we save it here
                 
                 print(paste0(as.character(folder.info$path),"/",input$panel,"/Results/"))
                 ###Change it to Proportions2021 if you want
                 folder.info$content1 = list.files(paste0(as.character(folder.info$path),"/",input$panel,"/Results/"), "Proportions2021",full.names = T)
                 print(folder.info$content1)
               }
  )
  
 #This setup the values for the main section depending on the QC protocol that is selected by the user
  
  observeEvent(input$protocol, {
    if (is.null(folder.info$content1)) {return(NULL)}
    if (length(folder.info$content1)==0)
    {
      values$df_cols <-NULL
    }else{
       if (input$protocol=="Boxplot_proportion" )
      {
        tmp <- colnames(read.csv(folder.info$content1[1],
                                 check.names = F))
        values$df_cols <- setdiff(grep(tmp,
                                       pattern = "ID*|Strain|Organ|Colony|All|Number|Gender|Assay|Genotype|FCS|Singlets|Live",value = T,invert = T),"")
        
      }else{
        tmp <- colnames(read.csv(folder.info$content1[1],
                                 check.names = F))
        values$df_cols <- grep(tmp,pattern = "G*",value = T)
      }
      
    }
  }) 
 #This shows the population names from the csv reportables in the Shiny app 
  output$checkbox <- renderUI({
    if(is.null(values$df_cols)){return("")}
    tags$div(align = 'left', 
             class = 'multicol', 
    checkboxGroupInput(inputId = "choice",width = '90%' ,
                       label = "Select populations", 
                       choices =values$df_cols,
                       selected =values$df_cols[1:2]))
    
    
  })
 
  
  #This is the main section of the Shiny, where the data is begin generated based on the protocol selected by the user, and generates the boxplot
  observeEvent(input$do, {
    print(input$protocol)
    
  
    if (input$protocol=="Boxplot_proportion")# when time is null this is the problem:| input$time=="Time variable: NA")
    {
      
      
      mat <- c()
     
     
        tmp <- read.table(folder.info$content1[1],check.names = F,header=T,sep=",")[]
        mat <- cbind(paste0(tmp$`FCS files`,":",tmp$Genotype,"_",tmp$Gender),tmp,deparse.level = F)
       colnames(mat)[1] <- "newID"
      }
      dat <<- make.data (mat = mat,cols = input$choice,sdconst = input$slidesd,extra.column=F)
      output$plot <- renderggiraph({
      
        main.name <- ifelse(input$protocol=="Boxplot_QC_Intra",yes ="QC values",
                            no =ifelse((input$protocol=="Boxplot_proportion" |input$protocol=="Boxplot_QC_Inter"),
                                       yes = "Proportion of parent-%",no = "MFIs") )
       
        r <- plot_results(alldata=dat$proportion,outlierdata=dat$outlier,
                          input$choice,name=main.name,sample.ref = NULL)
        box.plot <<-r
        x<-ggiraph(code = print(r),width_svg = 9 )
        x <- girafe_options(x, opts_selection(
          type = "multiple", css = "fill:#FF3333;stroke:black;"),
          opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;"))
        
      })
      
      
      output$down <- downloadHandler(
        filename =  function() {
          paste("boxplot", input$panel, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
          
            pdf(file) # open the pdf device
          # plot(x=x(), y=y(), main = "iris dataset plot", xlab = xl(), ylab = yl()) # draw the plot
          print(box.plot) # for GGPLOT
          dev.off()  # turn the device off
          
        } 
      )
      
      
    
      
      outlier.file <- dat$outlier
      output$download <- downloadHandler(
        filename = function() {
          paste("outliers.csv", sep = "")
        },
        content = function(file) {
          
          if (nrow(outlier.file)==0)
          {
            
            outlier.file <-"no outliers found."
            
          }else{
            outlier.file$intensity <- rep(0,length(outlier.file$proportion))
            outlier.file$intensity <- abs(round(( outlier.file$proportion-as.numeric( outlier.file$mean))/as.numeric( outlier.file$sd)))
            outlier.file$plot <- rep(0,length(outlier.file$proportion))
            outlier.file$plot <- paste0("ScatterPlots/",outlier.file$filesimple,".png")
          }
          write.csv(outlier.file, file)
        }
      )
      #This section is to plot the automated gating strategy if a point is selected from the interactive boxplot
      output$plotdens <- renderImage({
        click.input <- input$plot_selected
        click.pop <-tail(unlist(click.input),1)
        ind.samp <-which(dat$proportion$filesimple==click.pop)[1]
      print(click.input)
         print(click.pop)
         print(dat$proportion$file[1:3])
        samp <- paste0(dat$proportion$filesimple[ind.samp],".png")
       
        outfile <- paste0("~/impc/", input$panel,"/Results/Figures/ScatterPlots/",samp)
        print(outfile)
        if (length(outfile)>0)
        {
          
          list(src = outfile,
               contentType = 'image/png',
               width = 1200,
               height = 1000,
               alt = "-")
        }
      }, deleteFile = FALSE)
   
      
    
    
  })
}


