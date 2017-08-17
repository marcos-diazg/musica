library(shiny)

shinyServer(function(input, output) {

   options(shiny.maxRequestSize=300*1024^2)
   
   
   library(MutationalPatterns)
    
   ref_genome<-reactive({
      if (input$genome=="38") return ("BSgenome.Hsapiens.NCBI.GRCh38")
      if (input$genome=="19") return ("BSgenome.Hsapiens.UCSC.hg19")
      if (input$genome=="37") return ("BSgenome.Hsapiens.1000genomes.hs37d5")
   })


   

   
   
   vcfs<-reactive({
      library(ref_genome(), character.only = TRUE)
      inFile<-input$fileinput
      if (input$datatype=="vcf")
         return(read_vcfs_as_granges(inFile$datapath,inFile$dataname,ref_genome))
   })   
   
   
    
   output$result <- renderPrint({
    
      print(summary(vcfs))
    
   })
  
})
