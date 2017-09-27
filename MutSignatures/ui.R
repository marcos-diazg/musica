library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(shinythemes)

shinyUI(fluidPage(
   
   #CSS style specification
   theme = shinytheme("united"),
   
   #Activation of Shiny js
   useShinyjs(),
  
   #Title
   titlePanel("MutSignatures"),
  
   sidebarLayout(
      
      sidebarPanel(
         
         #Input format
         radioButtons("datatype", "Input format", c("VCF","MAF","TSV","Excel"),selected = "VCF",inline = TRUE),
         
         #Help menu for format of input file
         actionLink("helpformat","Help with input file format", icon=icon("question-circle-o")),
         bsModal("modal","HELP: Input file format","helpformat", includeHTML("../aux_files/help_with_input.html")),
         
         hr(),
         
         #File uploading
         fileInput("fileinput","Upload your file/s", multiple=TRUE),
#         helpText("Multiple files uploading is allowed"),
         
         hr(),
         
         #Genome selection
         selectInput("genome","Reference Genome",c("UCSC GRCh38/hg38"="hg38","UCSC GRCh37/hg19"="19","1000genomes hs37d5"="37"),selected="hg38"),
         
         #Run button
         actionButton("run","Run"),
         busyIndicator("Running",wait=0)
         
      ),
      


    #      sidebarPanel(
   
     # ),




      #Hidding tabs of mainpanel (results)
      hidden(
    
         mainPanel(id="mainpanel",
                   
            tabsetPanel(type="pills",
                        
               tabPanel("Mutational profile of provided sample/s",
                        br(),
                        plotOutput("prof96")
               ),
               
               tabPanel("Cosmic mutational signatures contributions",
                        br(),
                        downloadButton("download_contr",label="Download table"),
                        plotOutput("heatmap_signatures"),
                        dataTableOutput("contr")
               ),
               
               tabPanel("Comparison with other cancers",
                        br(),
                        downloadButton("download_known",label="Download table"),
                        plotOutput("heatmap_known")
               )
                          
            ) 
         )
      )
   )

#   ,bsTooltip("datatype","Choose the format of the input file")
))
