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

         uiOutput("selected_samples"),

         #Run button
         actionButton("run","Run"),
         busyIndicator("Running",wait=0)
         
      ),
      



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
                        downloadButton("download_signatures_plot_ID",label="Download plot"),
                        bsModal("modal_signatures","Download plot","download_signatures_plot_ID", 
                                radioButtons("type_signatures_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_signatures_plot","OK")),
                  
                        dataTableOutput("contr")
               ),
               
               tabPanel("Comparison with other cancers",
                        br(),
                        downloadButton("download_known",label="Download table"),
                        plotOutput("heatmap_known"),
                        downloadButton("download_known_plot_ID",label="Download plot"),
                        bsModal("modal_known","Download plot","download_known_plot_ID", 
                                radioButtons("type_known_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_known_plot","OK")),
                        selectInput("mycancers","Select the cancers to compare", c("All","Adrenocortical.carcinoma","ALL","AML","Bladder","Breast","Cervix","Chondrosarcoma","CLL","Colorectum","Glioblastoma","Glioma.Low.Grade","Head.and.Neck","Kidney.Chromophobe","Kidney.Clear.Cell","Kidney.Papillary","Liver","Lung.Adeno","Lung.Small.Cell","Lung.Squamous","Lymphoma.B.cell","Lymphoma.Hodgkin","Medulloblastoma","Melanoma","Myeloma","Nasopharyngeal.Carcinoma","Neuroblastoma","Oesophagus","Oral.gingivo.buccal.squamous","Osteosarcoma","Ovary","Pancreas","Paraganglioma","Pilocytic.Astrocytoma","Prostate","Stomach","Thyroid","Urothelial.Carcinoma","Uterine.Carcinoma","Uterine.Carcinosarcoma","Uveal.Melanoma"), multiple=TRUE, selectize = FALSE, size=15, selected="All")
               ),
               tabPanel("Principal Components Analysis",
                        br(),
                        plotOutput("pca_plot"),
                        downloadButton("download_pca_ID",label="Download plot"),
                        bsModal("modal_pca","Download plot","download_pca_ID", 
                                radioButtons("type_pca_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_pca_plot","OK"))
               )

            ) 
         )
      )
   )

#   ,bsTooltip("datatype","Choose the format of the input file")
))
