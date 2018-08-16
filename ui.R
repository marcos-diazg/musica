library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(V8)
library(shinythemes)
library(plotly)

shinyUI(fluidPage(
   
   #Google Analytics tracking
   tags$head(includeHTML(("google-analytics.html"))),   
	
   #CSS style specification
   theme = shinytheme("united"),
   
   tags$style(HTML("
				     .shiny-output-error { visibility: hidden; }
			        .shiny-output-error:before {
			          visibility: visible;
						 color: orangered;
			          content: 'ERROR. An error has occurred. Please contact the app authors at diaz2@clinic.cat or submit your bug-reports at https://github.com/marcos-diazg/musica/issues.'; }
			        }
				")),
   
   #Error management
	uiOutput("custom_error"),

   #Activation of Shiny js
   useShinyjs(),
   shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
  
   #Title
   titlePanel(title=div(img(src="musica.png"),"\t\tMutational Signatures in Cancer"),windowTitle = "MuSiCa - Mutational Signatures in Cancer"),
   sidebarLayout(
      
      sidebarPanel(
         
         #Input format
         radioButtons("datatype", "Input file format", c("VCF","TSV","Excel","MAF"),selected = "VCF",inline = TRUE),
         
         #Help menu for format of input file
         actionLink("helpformat","Help with input file format", icon=icon("question-circle-o")),
         bsModal("modal","HELP  Input file format","helpformat", includeHTML("./aux_files/help_with_input.html")),
         hr(),
         
         #File uploading
         fileInput("fileinput","Upload your file/s", multiple=TRUE),
         hr(),
         
         #Genome selection
         selectInput("genome","Reference genome",c("UCSC GRCh38/hg38"="hg38","UCSC GRCh37/hg19"="19","1000genomes hs37d5"="37"),selected="hg38"),
         hr(),
         
         #Type of Study
         radioButtons("studytype","Type of study", c("Whole Genome Sequencing", "Whole Exome Sequencing", "Targeted Sequencing")),
         

         uiOutput("kb_sequenced"),
         
         hr(),

         uiOutput("run_button"),

         
         #Stuff only showed when run button is pressed
         hidden(
            div(id="after_run",
              
               br(),
               br(),
              #Sample selection for plots (post push run)   
              uiOutput("mean_checkbox"),
              uiOutput("custom_error_style"),
              uiOutput("selected_samples"),
              
              #Cancer type selection for plots (post push run)
              uiOutput("selected_cancer_types"),
              
              hr(),
              
              uiOutput("col_dendro_heatmap"),
              uiOutput("row_dendro_heatmap"),
              uiOutput("col_dendro_cancers"),
              uiOutput("row_dendro_cancers"),
              uiOutput("heatmap_sample_names"),
              
              #uiOutput("error_input"),
              
              #Clear button 
              actionButton("clear","Clear")
            )
         )
      ),
      



      #Hidding tabs of mainpanel (results)
      hidden(
    
         mainPanel(id="mainpanel",
                   
            tabsetPanel(id="tab", type="pills",
                        
               #Somatic Mutation Prevalence
               tabPanel("Somatic mutation prevalence",value="smp",
                        br(),
                        downloadButton("download_smp_plot_ID",label="Download plot"),
                        bsModal("modal_smp","Download plot","download_smp_plot_ID",
                                radioButtons("type_smp_plot","Format",c("pdf","png","tiff"),selected="pdf"),
                                downloadButton("download_smp_plot","OK")),
                        downloadButton("download_smp_table",label="Download table"),
                        p(),
                        plotOutput("smp")
               ),
                        
               #Plot profile 96 changes                  
               tabPanel("Mutational profile",value="96prof",
                        br(),
                        downloadButton("download_prof96_plot_ID",label="Download plot"),
                        bsModal("modal_prof96","Download plot","download_prof96_plot_ID",
                                radioButtons("type_prof96_plot","Format",c("pdf","png","tiff"),selected="pdf"),
                                downloadButton("download_prof96_plot","OK")),
                        downloadButton("download_prof96_table",label="Download table"),
                        p(),
                        plotOutput("prof96")

               ),
               
               #Contribution of COSMIC mutational signatures (heatmap and table)
               tabPanel("COSMIC signatures contributions", value="contrib",
                        br(),
                        downloadButton("download_signatures_plot_ID",label="Download plot"),
                        bsModal("modal_signatures","Download plot","download_signatures_plot_ID", 
                                 radioButtons("type_signatures_plot","Format",c("pdf","png")),
                                 downloadButton("download_signatures_plot","OK")),
                        downloadButton("download_contr",label="Download table"),
                        p(),
                        fluidRow(plotlyOutput("heatmap_signatures",width="100%", height="800px")),
                        p(),
                        dataTableOutput("contr")
               ),
               
               #
               
               
               #Reconstructed mutational profile
               
               tabPanel("Reconstructed mutational profile", value="reconst",
               			tags$style(HTML(".shiny-output-error-reconstructed {visibility: hidden;}
               			.shiny-output-error-reconstructed:before {visibility: visible; color: orangered; content:'Sample selection error, please select just one sample at a time to visualize its reconstructed mutational profile.';}")),
               			
                        br(),
                        downloadButton("download_reconst_plot_ID",label="Download plot"),
                        bsModal("modal_reconst","Download plot","download_reconst_plot_ID",
                                radioButtons("type_reconst_plot","Format",c("pdf","png","tiff"),selected="pdf"),
                                downloadButton("download_reconst_plot","OK")),
                        downloadButton("download_reconst_table",label="Download table"),
                        p(),
                        plotOutput("reconst")
               ),
               
               
               
               #Comparison of COSMIC mutational signatures with cancers
               tabPanel("Comparison with cancers signatures", value="comp_canc_sign",
                        br(),
                         downloadButton("download_known_plot_ID",label="Download plot"),
                         bsModal("modal_known","Download plot","download_known_plot_ID", 
                                 radioButtons("type_known_plot","Format",c("pdf","png")),
                                 downloadButton("download_known_plot","OK")),
                        p(),
                        fluidRow(plotlyOutput("heatmap_known",width="100%", height="800px"),
                        img(src="legend.png",width=400))
               ),
               
               #Principal Component Analysis (PCA)
               tabPanel("Principal components analysis", value="pca",
               			tags$style(HTML(".shiny-output-error-pca {visibility: hidden;}
               			.shiny-output-error-pca:before {visibility: visible; color: orangered; content:'PCA analysis works only with 3 or more samples.';}")),
                        br(),
                        downloadButton("download_pca_ID",label="Download plot"),
                        downloadButton("download_pca_table",label="Download table"),
                        bsModal("modal_pca","Download plot","download_pca_ID", 
                                radioButtons("type_pca_plot","Format",c("pdf","png","tiff")),
                                downloadButton("download_pca_plot","OK")),
                        p(),
                        fluidRow(
                           column(4,tableOutput("pca_plot_table")),
                           column(8,plotOutput("pca_plot",height=650,width=650))
                           
                        )
               )

            ) 
         )
      )
   )
))
