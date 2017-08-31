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
         actionLink("helpformat","Help with input file format"),
         bsModal("modal","HELP: Input file format","helpformat",
                 HTML(
                     "<p>Format of the input file provided to this aplication should be one of the following options (multiple files uploading is allowed in all cases):
                     </p>
                     <p>
                     <u>VCF (Variant Call Format)</u>: Default file format for variant calling and annotation according to the specification provided <a href=https://samtools.github.io/hts-specs/VCFv4.2.pdf>here</a> (optimized for 4.2 version).
                     </p>
                     <p>
                     <u>MAF (Mutation Annotation Format)</u>: File format widely used for variant annotation in TCGA public available data. A complete specification may be consulted <a href=https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification>here</a>.
                     </p>
                     <p>
                     <u>TSV (Tab-Separated Values)</u>: plain text file with four required fields separated by a tab character. A header line is mandatory and the four headers must be provided as listed below:</p>
                     <ol>
                     <li>CHROM - Chromosome, ‘chr’ prefix is optional (e.g. ‘chr10’ and ’10’ are both valid).</li>
                     <li>POS - Chomosomal coordinate of reference allele.</li>
                     <li>REF - Reference allele at position given above.</li>
                     <li>ALT - Observed alternate allele.</li>
                     </ol>
                     </p>
                     <p>
                     <u>Excel</u>: Excel file (xls and xlsx extensions are both valid) with the same structure of four required fields described above (CHROM, POS, REF and ALT).
                     </p>")
         ),
         
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
