library(shiny)

shinyServer(function(input, output) {

   #Changing maximum file size for uploading
   options(shiny.maxRequestSize=300*1024^2)
   
   
   library(MutationalPatterns)
   library(BSgenome.Hsapiens.NCBI.GRCh38)
   library(BSgenome.Hsapiens.UCSC.hg19)
   library(BSgenome.Hsapiens.1000genomes.hs37d5)
    
   ref_genome<-reactive({
      if (input$genome=="38") return ("BSgenome.Hsapiens.NCBI.GRCh38")
      if (input$genome=="19") return ("BSgenome.Hsapiens.UCSC.hg19")
      if (input$genome=="37") return ("BSgenome.Hsapiens.1000genomes.hs37d5")
   })

 
   vcfs<-reactive({
      inFile<-input$fileinput
      if (input$datatype=="vcf")
         return(read_vcfs_as_granges(inFile$datapath,inFile$name,ref_genome()))
   })
   
   mut_mat <- reactive({ mut_matrix(vcfs(),ref_genome()) })
   
   output$prof96 <- renderPlot({
      plot_96_profile(mut_mat())
   })
   
   sp_url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
   cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
   cancer_signatures <- cancer_signatures[order(cancer_signatures[,1]),]
   cancer_signatures <- as.matrix(cancer_signatures[,4:33])
   
   fit_res <- reactive({ fit_to_signatures(mut_mat(), cancer_signatures) })
   
   proposed_etiology <- c("Age","APOBEC","BRCA1 / BRCA2","Smoking","Unknown (all cancers)","Defective DNA MMR","UV light","Unknown (breast cancer and medulloblastoma)","POLH (CLL, BCL)","POLE","Alkylating agents","Unknown (liver cancer)","APOBEC","Unknown (hypermutation)","Defective DNA MMR","Unknown (liver cancer)","Unknown (different cancers)","Unknown (different cancers)","Unknown (pilocytic astrocytoma)","Defective DNA MMR","Unknown (stomach cancer / MSI tumors)","Aristolochic acid","Unknown (liver cancer)","Aflatoxin","Unknown (Hodgkin lynphoma)","Defective DNA MMR","Unknown (kidney clear cell carcinomas)","Unknown (stomach cancer)","Tobacco chewing","Unknown (breast cancer) / NTHL1")
   

   #divisionRel function creation to print final dataframe
   divisionRel<-function(df){
      sum_df<-sapply(df,sum)
      for (i in 1:ncol(df)){
         df[,i]<-df[,i]/sum_df[i]
      }
      return(df)
   }
   
   
   output$contr <- renderDataTable({
      data.frame(colnames(cancer_signatures), proposed_etiology, divisionRel(as.data.frame(fit_res()$contribution)))
   })


   
})
