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
         bsModal("modal","HELP  Input file format","helpformat", includeHTML("../aux_files/help_with_input.html")),
         
         hr(),
         
         #File uploading
         fileInput("fileinput","Upload your file/s", multiple=TRUE),
         helpText("Multiple files uploading is allowed"),
         
         hr(),
         
         #Genome selection
         selectInput("genome","Reference Genome",c("UCSC GRCh38/hg38"="hg38","UCSC GRCh37/hg19"="19","1000genomes hs37d5"="37"),selected="hg38"),

         #Run button
         actionButton("run","Run",class = "btn-primary"),
         
         
         #Stuff only showed when run button is pressed
         hidden(
            div(id="after_run",
              #Adding space between run button and sample selector
              hr(),
              
              #Sample selection for plots (post push run)         
              uiOutput("selected_samples"),
              
              #Clear button and script to hide input file
              actionButton("clear","Clear"),
              
               tags$script('
                 Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {
                 var id = "#" + x + "_progress";
                 var idBar = id + " .bar";
                 $(id).css("visibility", "hidden");
                 $(idBar).css("width", "0%");
                 });
               '), 
              
              #BusyIndicator (post push run)
              busyIndicator("Running",wait=2)
            )
         )
      ),
      



      #Hidding tabs of mainpanel (results)
      hidden(
    
         mainPanel(id="mainpanel",
                   
            tabsetPanel(type="pills",
               
               #Plot profile 96 changes                  
               tabPanel("Mutational profile of provided sample/s",
                        br(),
                        plotOutput("prof96")
               ),
               
               #Contribution of COSMIC mutational signatures (heatmap and table)
               tabPanel("Cosmic mutational signatures contributions",
                        br(),
                        downloadButton("download_signatures_plot_ID",label="Download plot"),
                        bsModal("modal_signatures","Download plot","download_signatures_plot_ID", 
                                radioButtons("type_signatures_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_signatures_plot","OK")),
                        plotOutput("heatmap_signatures"),
                        downloadButton("download_contr",label="Download table"),
                        dataTableOutput("contr")
               ),
               
               #Comparison of COSMIC mutational signatures with cancers
               tabPanel("Comparison with cancers signatures",
                        br(),
                        downloadButton("download_known_plot_ID",label="Download plot"),
                        bsModal("modal_known","Download plot","download_known_plot_ID", 
                                radioButtons("type_known_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_known_plot","OK")),
                        plotOutput("heatmap_known"),
                        downloadButton("download_known",label="Download table"),
                        br(),
                        br(),
                        selectInput("mycancers","Select the cancers to compare", c("All","Adrenocortical.carcinoma","ALL","AML","Bladder","Breast","Cervix","Chondrosarcoma","CLL","Colorectum","Glioblastoma","Glioma.Low.Grade","Head.and.Neck","Kidney.Chromophobe","Kidney.Clear.Cell","Kidney.Papillary","Liver","Lung.Adeno","Lung.Small.Cell","Lung.Squamous","Lymphoma.B.cell","Lymphoma.Hodgkin","Medulloblastoma","Melanoma","Myeloma","Nasopharyngeal.Carcinoma","Neuroblastoma","Oesophagus","Oral.gingivo.buccal.squamous","Osteosarcoma","Ovary","Pancreas","Paraganglioma","Pilocytic.Astrocytoma","Prostate","Stomach","Thyroid","Urothelial.Carcinoma","Uterine.Carcinoma","Uterine.Carcinosarcoma","Uveal.Melanoma"), multiple=TRUE, selectize = FALSE, size=10, selected="All")
               ),
               
               #Principal Component Analysis (PCA)
               tabPanel("Principal Components Analysis",
                        br(),
                        
                        downloadButton("download_pca_ID",label="Download plot"),
                        bsModal("modal_pca","Download plot","download_pca_ID", 
                                radioButtons("type_pca_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_pca_plot","OK")),
                        plotOutput("pca_plot")
               )

            ) 
         )
      )
   )
))
