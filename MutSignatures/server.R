library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(shinythemes)

shinyServer(function(input, output) {
   
   #Hidding tabs of mainpanel (results)
   observeEvent(input$run,{
      shinyjs::show(id="mainpanel",anim=TRUE,animType="slide")
   })

   #Changing maximum file size for uploading
   options(shiny.maxRequestSize=300*1024^2)
   
   #Library loading (according to genome version)
   library(MutationalPatterns)
   library(reshape2)
   library(ggplot2)
   library(data.table)
   
   ref_genome<-eventReactive(input$run,{
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
 
   #Reading input files as GRanges objects
   vcfs<-eventReactive(input$run,{
      inFile<-input$fileinput

         #VCF
         if (input$datatype=="VCF") return(read_vcfs_as_granges(inFile$datapath,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         
         #MAF
         if (input$datatype=="MAF"){
            aux<-fread(inFile$datapath,header=T,sep="\t",skip="#",data.table=F)
            aux<-aux[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
            colnames(aux)[1:4]<-c("#CHROM","POS","REF","ALT")
            
            #Condition in case "chr" prefix is present at CHROM column in input file
            if (length(grep("chr",aux))>0){
               aux[,c("#CHROM")]<-sapply(strsplit(aux[,c("#CHROM")],"chr"),"[",2)
            }
            
            aux$ID<-"."
            aux$QUAL<-"."
            aux$FILTER<-"PASS"
            aux<-aux[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","Tumor_Sample_Barcode")]
            aux_list<-split(aux,f=aux$Tumor_Sample_Barcode)
            aux_list<-lapply(aux_list,"[",c(1:7))
            ff<-sapply(aux_list,function(i){tempfile(pattern="tp",fileext=".vcf")})
            for (i in 1:length(aux_list)){
               write.table(aux_list[[i]],file=ff[i],row.names=F,quote=F,sep="\t")
            }
            return(read_vcfs_as_granges(ff,names(ff),ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }
      
         #TSV
         if (input$datatype=="TSV"){
            aux<-fread(inFile$datapath,header=T,sep="\t",data.table=F)
            
            #Condition in case "chr" prefix is present at CHROM column in input file
            if (length(grep("chr",aux))>0){
               aux$CHROM<-sapply(strsplit(aux$CHROM,"chr"),"[",2)
            }
            
            colnames(aux)[1]<-"#CHROM"
            aux$ID<-"."
            aux$QUAL<-"."
            aux$FILTER<-"PASS"
            aux<-aux[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")]
            ff<-tempfile("tp",fileext=".vcf")
            write.table(aux,file=ff,row.names=F,quote=F,sep="\t")
            
            return(read_vcfs_as_granges(ff,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }
      
         #Excel
         if (input$datatype=="Excel"){
            library(xlsx)
            aux<-read.xlsx(inFile$datapath,1)
            
            #Condition in case "chr" prefix is present at CHROM column in input file
            if (length(grep("chr",aux))>0){
               aux$CHROM<-sapply(strsplit(aux$CHROM,"chr"),"[",2)
            }
            
            colnames(aux)[1]<-"#CHROM"
            aux$ID<-"."
            aux$QUAL<-"."
            aux$FILTER<-"PASS"
            aux<-aux[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")]
            ff2<-tempfile("tp",fileext=".vcf")
            write.table(aux,file=ff2,row.names=F,quote=F,sep="\t")
            
            return(read_vcfs_as_granges(ff2,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }


   })
      
   mut_mat <- reactive({
         return(mut_matrix(vcfs(),ref_genome()))
   })
      
      
   #Output 96 nucleotide changes profile (samples individually)   
   output$prof96 <- renderPlot({
         plot_96_profile(mut_mat())
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

   output$download_contr <- downloadHandler( filename="contr.csv", content=function (file){ write.table(x=data.frame(colnames(cancer_signatures), proposed_etiology, divisionRel(as.data.frame(fit_res()$contribution))), file=file, row.names=FALSE, sep="\t", quote=FALSE) })
      
   output$heatmap_known <- renderPlot({
      a<-t(data.frame(fit_res()$contribution[30:1,], known_cancer_signatures[30:1,]))
      colnames(a)<-colnames(cancer_signatures)[30:1]
      
      for (i in 1:(nrow(a)-40)) {
         a[i,]<-a[i,]/max(a[i,])
      }
      a.m<-reshape2::melt(as.matrix(a)) 
      a.m$category<-rep(c(rep("Sample",nrow(a)-40),rep("Cancers",40)),30)
      sel<-which(a.m$category=="Cancers")
      a.m[sel,"value"]<-a.m[sel,"value"]+2
      a.m[is.na(a.m)] <- 0
      
      colorends <- c("white","red", "white", "blue")
      
      ggplot(a.m, aes(x=Var1, y=Var2)) + geom_tile(aes(fill = value),
         colour = "white") + theme(axis.text.x=element_text(angle=90)) +
         scale_fill_gradientn(colours = colorends, limits = c(0,3)) + labs(x="",y="")
   })
   
   output$heatmap_signatures <- renderPlot({
      a<-t(divisionRel(data.frame(fit_res()$contribution[30:1,])))
      colnames(a)<-colnames(cancer_signatures)[30:1]
      a.m<-reshape2::melt(as.matrix(a)) 
      colorends <- c("white","red")
      ggplot(a.m, aes(x=Var1, y=Var2)) + geom_tile(aes(fill = value),
       colour = "white") + theme(axis.text.x=element_text(angle=90)) +
       scale_fill_gradientn(colours = colorends, limits = c(0,max(a))) + labs(x="",y="")
   })
   
   
})
