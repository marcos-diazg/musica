library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)

shinyServer(function(input, output) {
   
   #Hidding tabs of mainpanel (results)
   observeEvent(input$run,{
      shinyjs::show(id="mainpanel",anim=TRUE,animType="slide")
   })

   #Changing maximum file size for uploading
   options(shiny.maxRequestSize=300*1024^2)
   
   #Library loading (according to genome version)
   library(MutationalPatterns)

   ref_genome<-eventReactive(input$run,{
      if (input$genome=="38"){
         library("BSgenome.Hsapiens.NCBI.GRCh38")
         return ("BSgenome.Hsapiens.NCBI.GRCh38")
      }
      if (input$genome=="19"){
         library("BSgenome.Hsapiens.UCSC.hg19")
         return ("BSgenome.Hsapiens.UCSC.hg19")
      }
      if (input$genome=="37"){
         library ("BSgenome.Hsapiens.1000genomes.hs37d5")
         return ("BSgenome.Hsapiens.1000genomes.hs37d5")
      }
      if (input$genome=="hg38"){
         library("BSgenome.Hsapiens.UCSC.hg38")
         return("BSgenome.Hsapiens.UCSC.hg38")
      }
   })
 
   #Reading vcf files as GRanges objects
   vcfs<-eventReactive(input$run,{
      inFile<-input$fileinput

#      withProgress(message="Running",value=0,{
      
         if (input$datatype=="vcf") return(read_vcfs_as_granges(inFile$datapath,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         
#         incProgress(0.5,detail="GRanges")
         
#      })
      
   })
      
   mut_mat <- reactive({
#      withProgress(message="Running",value=0.5, {
         return(mut_matrix(vcfs(),ref_genome()))
         
#         incProgress(0.5,detail="mut_matrix")
#      })
         
         
   })
      
      
      
   output$prof96 <- renderPlot({
#      withProgress(message="Running",value=0,{
         
#         incProgress(0.5,detail="a")
         
         plot_96_profile(mut_mat())
         
#         incProgress(0.5,detail="b")
 #     })
   })

   
   sp_url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
   cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
   cancer_signatures <- cancer_signatures[order(cancer_signatures[,1]),]
   cancer_signatures <- as.matrix(cancer_signatures[,4:33])
   
   fit_res <- reactive({ fit_to_signatures(mut_mat(), cancer_signatures) })
   
   proposed_etiology <- c("Age","APOBEC","BRCA1 / BRCA2","Smoking","Unknown (all cancers)","Defective DNA MMR","UV light","Unknown (breast cancer and medulloblastoma)","POLH (CLL, BCL)","POLE","Alkylating agents","Unknown (liver cancer)","APOBEC","Unknown (hypermutation)","Defective DNA MMR","Unknown (liver cancer)","Unknown (different cancers)","Unknown (different cancers)","Unknown (pilocytic astrocytoma)","Defective DNA MMR","Unknown (stomach cancer / MSI tumors)","Aristolochic acid","Unknown (liver cancer)","Aflatoxin","Unknown (Hodgkin lynphoma)","Defective DNA MMR","Unknown (kidney clear cell carcinomas)","Unknown (stomach cancer)","Tobacco chewing","Unknown (breast cancer) / NTHL1")
   
   known_cancer_signatures<-read.table("cancermatrix.csv",header=TRUE,sep="\t",row.names=1)
   
   #divisionRel function creation to print final dataframe
   divisionRel<-function(df){
      sum_df<-sapply(df,sum)
      for (i in 1:ncol(df)){
         df[,i]<-round((df[,i]/sum_df[i]),3)
      }
      return(df)
   }
   
   
   output$contr <- renderDataTable({
      data.frame(Signature = 1:30, Proposed_Etiology = proposed_etiology, divisionRel(as.data.frame(fit_res()$contribution)))
   },options = list(lengthChange=FALSE,pageLength=30, paging=FALSE))

   output$download_contr <- downloadHandler( filename="contr.csv", content=function (file){ write.table(x=data.frame(colnames(cancer_signatures), proposed_etiology, fit_res()$contribution), file=file, row.names=FALSE, sep="\t", quote=FALSE) })
   
   output$known<- renderDataTable({
      data.frame(Signature = 1:30, divisionRel(as.data.frame(fit_res()$contribution)), known_cancer_signatures[c(-31),] )
      },options = list(lengthChange=FALSE,pageLength=30, paging=FALSE))

   output$download_known <- downloadHandler( filename="known.csv", content=function (file){ write.table(x=data.frame(colnames(cancer_signatures), fit_res()$contribution, known_cancer_signatures[c(-31),]), file=file, row.names=FALSE, sep="\t", quote=FALSE) })
   
   
})
