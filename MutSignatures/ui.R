library(shiny)

shinyUI(fluidPage(
  
   titlePanel("Mutational Signatures"),
  
   sidebarLayout(
      sidebarPanel(
         helpText("Upload your file:"),
         fileInput("fileinput","", multiple=TRUE),
         radioButtons("datatype", "Formato", c("vcf")),
         selectInput("genome","Reference Genome",c("NCBI GRCh38"="38","UCSC hg19"="19","1000genomes hs37d5"="37"),selected="37")
      ),
    
      mainPanel(
         tabsetPanel(type="tabs",
            tabPanel("Mutational profile of provided sample/s", plotOutput("prof96")),
            tabPanel("Cosmic mutational signatures contributions", downloadButton("download_contr",label="Download table"), dataTableOutput("contr")),
            tabPanel("Comparison with other cancers", downloadButton("download_known",label="Download table"), plotOutput("heatmap_known"), dataTableOutput("known"))
            )
      )
   )
))
