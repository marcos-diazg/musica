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
   titlePanel("Mutational Signatures"),
  
   sidebarLayout(
      
      sidebarPanel(

         radioButtons("datatype", "Input format", c("vcf","tsv","Excel"),selected = "vcf",inline = TRUE),
         
         #Help menu for format of input file
         actionLink("helpformat","Help with input file format"),
         bsModal("modal","HELP: Input file format","helpformat",
                 HTML(
                     "<p>Format of the input file provided to this aplication should be one of the following options:
                     </p>
                     <p>
                     <u>vcf</u>: Variant Call Format according to the specification provided <a href=https://samtools.github.io/hts-specs/VCFv4.2.pdf>here</a> (optimized for 4.2 version).
                     </p>
                     <p>
                     <u>tsv</u>: Tab Separated Values plain text file with four required fields. A header line is mandatory and the four headers must be provided as listed below:</p>
                     <ol>
                     <li>CHROM - Chromosome, ‘chr’ prefix is optional (e.g. ‘chr10’ and ’10’ are both valid).</li>
                     <li>POS - Chomosomal coordinate of reference allele.</li>
                     <li>REF - Reference allele at position given above.</li>
                     <li>ALT - Observed alternate allele.</li>
                     </ol>
                     </p>
                     <p>
                     <u>Excel</u>: Excel file (xls and xlsx extensions are both valid) with the same structure of four required fields described above (CHROM, POS, REF and ALT).
                     </p>
                     
                     ")
         ),
         
         hr(),
         
         fileInput("fileinput","Upload your file/s", multiple=TRUE),
         helpText("Multiple files uploading is allowed"),
         
         hr(),
         
         selectInput("genome","Reference Genome",c("NCBI GRCh38"="38","UCSC hg38"="hg38","UCSC hg19"="19","1000genomes hs37d5"="37"),selected="37"),
         
         actionButton("run","Run"),
         busyIndicator("Calculation in progress",wait=0)
         
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
                        dataTableOutput("contr")
               ),
               
               tabPanel("Comparison with other cancers",
                        br(),
                        downloadButton("download_known",label="Download table"),
                        dataTableOutput("known")
               )
                               
            )
         )
      )
   ),

   bsTooltip("datatype","Choose the format of the input file")
))
