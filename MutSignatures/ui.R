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
       plotOutput("prof96"),
       dataTableOutput("contr")
    )
  )
))
