library(shiny)
library(shinyFiles)
library(flowCore)
library(flowDensity)
library(ncdfFlow)
library(flowWorkspace)
library(gplots)
library(evaluate)
library(e1071)
library(Cairo)
library(CytoML)
library(flowCut)
library(shinythemes)
library(shinyWidgets)
library(ggplot2)
library(ggiraph)
library(shinyjs)
#When running this, please keep in mind to have the floder setting exactly as shown in the ReadMe file,
#as the folder paths point to this setting
source('shinyHelper.R')

ncdf.version <- 3
jsResetCode <- "shinyjs.reset = function() {history.go(0)}" # Define the js method that resets the page
#This is the ui of the Shiny app, with a custom coloring for the drop-down menu for the tube section
ui<-fluidPage( useShinyjs(),theme = shinytheme("lumen"), setBackgroundColor(
  color = c("#FFFFFF", "#3399FF"),
  gradient = "linear",
  direction = "bottom"
), tags$head(
  tags$style(
    HTML('
             .multicol{ 
             -webkit-column-count: 3; /* Chrome, Safari, Opera */ 
             -moz-column-count: 3;    /* Firefox */ 
             column-count: 3; 
             -moz-column-fill: auto;
             -column-fill: auto;
             }
             '))
  ),
title = "IMPC: Automated analysis postprocessing",fluidRow(
  #This is to enter the main directory where the results are
  column(3,offset = 1,
         shinyDirButton("result",label="Choose a parent directory where all saved results are", "Select")
  )
  ,
  column(3,offset = 10,
         img(src = "bccancer_logo_colour.png", height = 70, width = 200)
  )),
fluidRow(
  
  column(3,offset = 10,
         img(src = "impc-logo.png", height = 70, width = 190)))
,fluidRow(
  column(4, 
         #this is the main part that the user selects the QC protocol
         radioButtons(inputId="protocol", 
                      label="Which test:", 
                      choices = c("None selected" = "", 
                                  "Boxplot_proportion"),
                      selected = NULL)),
  column(4, 
         #this is the main part that the user selects the QC protocol
         radioButtons(inputId="panel", 
                      label="Which panel:", 
                      choices = c("Panel1", 
                                  "Panel2"),
                      selected = "Panel1")),
  column(4, sliderInput("slidesd",label = "SD factor",min = 1.5,max = 8,step = .5,value = 3)))

  #This shows the population name extracted from the Shiny app, from the csv reportables
  ,fluidRow(column(8,uiOutput("checkbox"))),
  column(1,offset = 0,actionButton("do", "Perform QC")),
#This is where the boxplot will be shown
fluidRow(column(offset=0, style='padding-top:0px;',width = 12,
                ggiraphOutput("plot"))),

#This will generate the gating plot if a point is selected in the boxplot

#Tese would download the data

fluidRow(column(1,offset=0,downloadButton("down", "Download the boxplot"))),
fluidRow(column(1,offset=0,downloadButton("download", "Download outlier csv"))),


#This will generate the gating plot if a point is selected in the boxplot

div(style = "font-size: 10px; padding: 0px 0px; margin-top:-2em", 
    fluidRow(imageOutput("plotdens")))

)






