library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(V8)
library(shinythemes)
library(plotly)

shinyUI(fluidPage(
   
   #CSS style specification
   theme = shinytheme("united"),
   
   #Activation of Shiny js
   useShinyjs(),
   shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
  
   #Title
   titlePanel("MuSiC     -     Mutational Signatures in Cancer"),
   sidebarLayout(
      
      sidebarPanel(
         
         #Input format
         radioButtons("datatype", "Input format", c("VCF","MAF","TSV","Excel"),selected = "VCF",inline = TRUE),
         
         #Help menu for format of input file
         actionLink("helpformat","Help with input file format", icon=icon("question-circle-o")),
         bsModal("modal","HELP  Input file format","helpformat", includeHTML("../aux_files/help_with_input.html")),
         hr(),
         
         #File uploading
         fileInput("fileinput","Upload your file/s", multiple=TRUE),
         hr(),
         
         #Genome selection
         selectInput("genome","Reference Genome",c("UCSC GRCh38/hg38"="hg38","UCSC GRCh37/hg19"="19","1000genomes hs37d5"="37"),selected="hg38"),
         hr(),

         #Run button
         actionButton("run","Run",class = "btn-primary"),
         

         #Busy indicator
         busyIndicator("Running",wait=0),

         
         #Stuff only showed when run button is pressed
         hidden(
            div(id="after_run",
              
              #Sample selection for plots (post push run)         
              uiOutput("selected_samples"),
              
              #Cancer type selection for plots (post push run)
              uiOutput("selected_cancer_types"),
              
              hr(),
              
              #Clear button 
              actionButton("clear","Clear")
            )
         )
      ),
      



      #Hidding tabs of mainpanel (results)
      hidden(
    
         mainPanel(id="mainpanel",
                   
            tabsetPanel(id="tab", type="pills",
               
               #Plot profile 96 changes                  
               tabPanel("Mutational profile of provided sample/s",value="96prof",
                        br(),
                        downloadButton("download_prof96_plot_ID",label="Download plot"),
                        bsModal("modal_prof96","Download plot","download_prof96_plot_ID",
                                radioButtons("type_prof96_plot","Format",c("pdf","png","tiff"),selected="pdf"),
                                downloadButton("download_prof96_plot","OK")),
                        br(),
                        plotOutput("prof96")

               ),
               
               #Contribution of COSMIC mutational signatures (heatmap and table)
               tabPanel("Cosmic mutational signatures contributions",
                        br(),
                        uiOutput("col_dendro_heatmap"),
                        uiOutput("row_dendro_heatmap"),
                         downloadButton("download_signatures_plot_ID",label="Download plot"),
                         bsModal("modal_signatures","Download plot","download_signatures_plot_ID", 
                                 radioButtons("type_signatures_plot","Format",c("pdf","png","tiff")),
                                 downloadButton("download_signatures_plot","OK")),
                        p(),
                        fluidRow(plotlyOutput("heatmap_signatures",width="100%", height="500px")),
                        p(),
                        downloadButton("download_contr",label="Download table"),
                        br(),
                        dataTableOutput("contr")
               ),
               
               #Comparison of COSMIC mutational signatures with cancers
               tabPanel("Comparison with cancers signatures", value="comp_canc_sign",
                        br(),
                        uiOutput("col_dendro_cancers"),
                        uiOutput("row_dendro_cancers"),
                         downloadButton("download_known_plot_ID",label="Download plot"),
                         bsModal("modal_known","Download plot","download_known_plot_ID", 
                                 radioButtons("type_known_plot","Format",c("pdf","png","tiff")),
                                 downloadButton("download_known_plot","OK")),
                        p(),
                        fluidRow(plotlyOutput("heatmap_known",width="100%", height="500px"))
               ),
               
               #Principal Component Analysis (PCA)
               tabPanel("Principal Components Analysis", value="pca",
                        br(),
                        downloadButton("download_pca_ID",label="Download plot"),
                        bsModal("modal_pca","Download plot","download_pca_ID", 
                                radioButtons("type_pca_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_pca_plot","OK")),
                        br(),
                        plotOutput("pca_plot")
               )

            ) 
         )
      )
   )
))
