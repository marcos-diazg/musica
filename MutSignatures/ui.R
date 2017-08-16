library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Mutational Signatures"),
  
  sidebarLayout(
    sidebarPanel(
       fileInput("fileinput",""),
       radioButtons("datatype", "Formato", c("vcf"))
    ),
    
    mainPanel(
       plotOutput("resultplot")
    )
  )
))
