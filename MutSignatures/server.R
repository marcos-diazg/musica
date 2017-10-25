library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(V8)
library(shinythemes)
#webshot::install_phantomjs()  #To plot into pdf Heatmaply

shinyServer(function(input, output,session){
   
   #Setting maximum file size for uploading (500 MB)
   options(shiny.maxRequestSize=500*1024^2)
   
   
   #Resolution of the tiff images
   ppi<-200
   
   
   #Hidding/Showing tabs of mainpanel, run and clear buttons, ...
   observeEvent(input$run,{
      shinyjs::show(id="mainpanel")
      shinyjs::hide(id="run")
      shinyjs::show(id="after_run")
   })
   
   observeEvent(input$clear, {
      shinyjs::js$refresh()
   })
   
   
   
   #Library loading
   library(MutationalPatterns)
   library(reshape2)
   library(ggplot2)
   library(data.table)
   library(VariantAnnotation)
   library(plotly)
   library(heatmaply)
   library(gplots)
   library(xlsx)

   #######################################
   #Reference genome definition and loading [ref_genome]
   #######################################
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
 
   #######################################
   #Reading input files as GRanges objects [vcfs]
   #######################################
   vcfs<-eventReactive(input$run,{
      inFile<-input$fileinput

         #VCF
         if (input$datatype=="VCF"){
            #Warning for file format
            validate(
               need(length(grep(".vcf",inFile$datapath))>0 | length(grep(".txt",inFile$datapath))>0,"File format error, please select the correct input file format before uploading your file/s.")
            )
            
            
            ########################################################################
            #Filtering steps
            #vcfilter<-readVcfAsVRanges(inFile$datapath)
            ########################################################################
            
            
            #Read vcf for MutationalPatterns
            return(read_vcfs_as_granges(inFile$datapath,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }
      
      
         #MAF
         if (input$datatype=="MAF"){
            
            #Warning for file format
            validate(
               need(length(grep(".maf",inFile$datapath))>0 | length(grep(".txt",inFile$datapath))>0,"File format error, please select the correct input file format before uploading your file/s."),
               need(length(inFile$datapath)==1, "Only one multi-sample MAF file is allowed")
            )
            
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
            
            #Warning for file format
            validate(
               need(length(grep(".tsv",inFile$datapath))>0 | length(grep(".txt",inFile$datapath))>0,"File format error, please select the correct input file format before uploading your file/s.")
            )
            
            ff_list<-list()
            for (w in 1:length(inFile$datapath)){
               aux<-fread(inFile$datapath[w],header=T,sep="\t",data.table=F)
               
               #Condition in case "chr" prefix is present at CHROM column in input file
               if (length(grep("chr",aux))>0){
                  aux$CHROM<-sapply(strsplit(aux$CHROM,"chr"),"[",2)
               }
            
               colnames(aux)[1]<-"#CHROM"
               aux$ID<-"."
               aux$QUAL<-"."
               aux$FILTER<-"PASS"
               aux<-aux[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")]
               ff_list[[w]]<-tempfile("tp",fileext=".vcf")
               write.table(aux,file=ff_list[[w]],row.names=F,quote=F,sep="\t")
            }
            
            ff<-do.call("c",ff_list)
            
            return(read_vcfs_as_granges(ff,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }
      
         #Excel
         if (input$datatype=="Excel"){
            
            #Warning for file format
            validate(
               need(length(grep(".xlsx",inFile$datapath))>0 | length(grep(".xls",inFile$datapath))>0,"File format error, please select the correct input file format before uploading your file/s.")
            )
            
            ff_list<-list()
            for (w in 1:length(inFile$datapath)){
               aux<-read.xlsx(inFile$datapath[w],1)
               
               #Condition in case "chr" prefix is present at CHROM column in input file
               if (length(grep("chr",aux))>0){
                  aux$CHROM<-sapply(strsplit(aux$CHROM,"chr"),"[",2)
               }
               
               colnames(aux)[1]<-"#CHROM"
               aux$ID<-"."
               aux$QUAL<-"."
               aux$FILTER<-"PASS"
               aux<-aux[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")]
               ff_list[[w]]<-tempfile("tp",fileext=".vcf")
               write.table(aux,file=ff_list[[w]],row.names=F,quote=F,sep="\t")
            }
            
            ff<-do.call("c",ff_list)
               
            return(read_vcfs_as_granges(ff,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }

   })
   
   
   #######################################
   #Mutation Matrix creation [mut_mat]
   #######################################
   mut_mat <- reactive({
         return(mut_matrix(vcfs(),ref_genome()))
   })
      
   
   #######################################
   #COSMIC Mutational Signatures loading (and adjustment) from COSMIC website [cancer_signatures]
   #######################################
   sp_url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
   cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
   cancer_signatures <- cancer_signatures[order(cancer_signatures[,1]),]
   cancer_signatures <- as.matrix(cancer_signatures[,4:33])
   
   
   #######################################
   #Fitting mutations in samples (mut_mat) to COSMIC signatures [fit_res]
   #######################################
   fit_res <- reactive({ fit_to_signatures(mut_mat(), cancer_signatures) })
         
   
   #Auxiliar files of aetiology and known signatures by cancer type (from COSMIC)
   proposed_etiology <- fread("../aux_files/proposed_etiology_COSMIC_signatures.txt",sep="\t",header=F,data.table=F)[,2]
   known_cancer_signatures<-read.table("../aux_files/cancermatrix.tsv",header=TRUE,sep="\t",row.names=1)
   
   
   #divisionRel function creation to print final dataframe
   divisionRel<-function(df){
      sum_df<-sapply(df,sum)
      for (i in 1:ncol(df)){
         df[,i]<-round((df[,i]/sum_df[i]),3)
      }
      return(df)
   }
   

   #Plot selectize to select samples to plot.
   output$selected_samples<-renderUI({
      
      if (input$tab == "96prof"){
        
         mysamp<-c("All",colnames(as.data.frame(fit_res()$contribution)))
         selectInput("mysamp","Select your samples",mysamp, multiple=TRUE, selectize = FALSE, size=6, selected="All")
         
      } else {
      
         mysamp<-c("All",colnames(as.data.frame(fit_res()$contribution)),"mean")
         selectInput("mysamp","Select your samples",mysamp, multiple=TRUE, selectize = FALSE, size=6, selected="All")
         
      }
      
   })
   
   output$selected_cancer_types<-renderUI({
      
      if (input$tab=="comp_canc_sign"){
         
         selectInput("mycancers","Select the cancers to compare", c("All",colnames(read.table("../aux_files/cancermatrix.tsv",header=TRUE,sep="\t",row.names=1))), multiple=TRUE, selectize=FALSE, size=10, selected="All")
      
      }
      
   })
 
   
   #Select which samples use to plot.
   my_contributions<- reactive({ 
         
      if ("All" %in% input$mysamp) {
         aux<-divisionRel(as.data.frame(fit_res()$contribution))
         con<-data.frame(aux, mean = apply(aux,1,mean))
            
      } else {
            
         if("mean" %in% input$mysamp) {      
            if (length(input$mysamp)>1) {
               aux<-divisionRel(as.data.frame(fit_res()$contribution[,input$mysamp[-length(input$mysamp)]]))
               con<-data.frame(aux, mean = apply(aux,1,mean)) 
               
            } else {
               aux<-divisionRel(as.data.frame(fit_res()$contribution))
               con<-data.frame(mean = apply(aux,1,mean))    
            }
              
         } else {
            con<-data.frame(divisionRel(as.data.frame(fit_res()$contribution[,input$mysamp])))  
         }
      }
         
      
      #Fixing colname of one sample (without mean)
      if (ncol(con)==1 & colnames(con)[1]!="mean"){
         colnames(con)<-setdiff(input$mysamp,c("All","mean"))
      }
         
      #Fixing colname of one sample (with mean)
      if (ncol(con)==2 & colnames(con)[2]=="mean" & input$mysamp!="All"){
         colnames(con)[1]<-setdiff(input$mysamp,c("All","mean"))
      }
         
      #Fixing colname of just mean of samples
      if (ncol(con)==1 & colnames(con)[1]=="mean"){
         colnames(con)<-"mean"
      }
         
      return(con)
      
   })


   #######################################
   #PLOT 96 nucleotide changes profile (samples individually)
   #######################################
   
   #Plot 96 profile
   output$prof96 <- renderPlot({
      aux_96_profile<-as.matrix(mut_mat()[,setdiff(colnames(my_contributions()),c("mean"))])
      colnames(aux_96_profile)<-setdiff(colnames(my_contributions()),c("mean"))
      plot_96_profile(aux_96_profile)
   })
   
   #Download Plot 96 profile 
   output$download_prof96_plot <- downloadHandler (
      filename = function(){
         paste("prof96_plot",input$type_prof96_plot, sep=".")
      },
      content = function(ff) {
         aux_96_profile<-as.matrix(mut_mat()[,setdiff(colnames(my_contributions()),c("mean"))])
         colnames(aux_96_profile)<-setdiff(colnames(my_contributions()),c("mean"))
         plot_96_profile(aux_96_profile)
         ggsave(ff,height=7,width=7,dpi=ppi)
      }
   )
      
      
   #######################################
   ### Plot heatmap with contributions
   #######################################
   
   
   #DataTable
   output$contr <- renderDataTable({
         data.frame(Signature = 1:30, Proposed_Etiology = proposed_etiology, my_contributions())
      },
      options = list(lengthChange=FALSE,pageLength=30, paging=FALSE, searching=FALSE, info=FALSE)
   )
   
   
   #Download Table
   output$download_contr <- downloadHandler( filename="COSMIC_sign_contributions.txt", content=function (file){ write.table(x = data.frame(Signature = 1:30, Proposed_Etiology = proposed_etiology, my_contributions()), file = file, sep = "\t", quote=F, row.names=F) })
   
   
   #check if column or row dendogram is needed
   output$row_dendro_heatmap<-renderUI({
      radioButtons("row_d_heatmap", "Row dendrogram", c("yes","no"),selected = "no",inline = TRUE)
   })
   output$col_dendro_heatmap<-renderUI({
      radioButtons("col_d_heatmap", "Column dendrogram", c("yes","no"),selected = "no",inline = TRUE)
   })
   
   
   #HeatMap
   output$heatmap_signatures <- renderPlotly({
      a<-my_contributions()
      if (ncol(a)==1) colnames(a)<-colnames(my_contributions()) ## fix colnames when there is only one sample
      rownames(a)<-colnames(cancer_signatures)[1:30] 
      colorends <- c("white","red")
      dendro <- "none"
      if (input$row_d_heatmap=="yes") dendro<-"row"
      if (input$col_d_heatmap=="yes") dendro<-"column" 
      if (input$row_d_heatmap=="yes" & input$col_d_heatmap=="yes") dendro<-"both"
   
      heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,1)),
                dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90)
   })
   
   
   #Download HeatMap 
   # output$download_signatures_plot <- downloadHandler (
   #    filename = function(){paste("signatures_plot",input$type_signatures_plot, sep=".")}, 
   #    content = function(ff) {
   #       
   #       a<-my_contributions()
   #       if (ncol(a)==1) colnames(a)<-colnames(my_contributions()) ## fix colnames when there is only one sample
   #       rownames(a)<-colnames(cancer_signatures)[1:30] 
   #       colorends <- c("white","red")
   #       dendro <- "none"
   #       if (input$row_d_heatmap=="yes") dendro<-"row"
   #       if (input$col_d_heatmap=="yes") dendro<-"column" 
   #       if (input$row_d_heatmap=="yes" & input$col_d_heatmap=="yes") dendro<-"both"
   #       
   #       if (input$type_signatures_plot=="pdf") pdf(ff,height=7,width=7)
   #       if (input$type_signatures_plot=="png") png(ff,height=7*ppi,width=7*ppi,res=ppi)
   #       if (input$type_signatures_plot=="tiff") tiff(ff,height=7*ppi,width=7*ppi,res=ppi,compression="lzw")
   # 
   #        heatmap.2(
   #          as.matrix(a),
   #          trace = "none",
   #          Rowv = input$row_d_heatmap=="yes" , 
   #          Colv = input$col_d_heatmap=="yes" ,
   #          col = colorpanel (256, low = "white", high = "red"),
   #          xlab = FALSE
   #       )
   #        
   #       dev.off()
   #       
   #    })
   
   
   
   #######################################
   ### Plot - Comparison with other cancers
   #######################################
   
   #check if column or row dendogram is needed
   output$row_dendro_cancers<-renderUI({
      radioButtons("row_c_heatmap", "Row dendrogram", c("yes","no"),selected = "no",inline = TRUE)
   })
   output$col_dendro_cancers<-renderUI({
      radioButtons("col_c_heatmap", "Column dendrogram", c("yes","no"),selected = "no",inline = TRUE)
   })
   
   
   #HeatMap
   output$heatmap_known <- renderPlotly({
      
      if ("All" %in% input$mycancers) my.sel.cancers<-colnames(known_cancer_signatures)
      else my.sel.cancers<-intersect(input$mycancers,colnames(known_cancer_signatures))
      

      a<-data.frame(my_contributions()[1:30,], known_cancer_signatures[1:30,my.sel.cancers])
      rownames(a)<-colnames(cancer_signatures)[1:30]
      if (ncol(my_contributions())==1) colnames(a)[1]<-colnames(my_contributions()) ## fix colnames when there is only one sample
      if (length(my.sel.cancers)==1) colnames(a)[length(colnames(a))]<-my.sel.cancers ## fix colnames when there is only one cancer type
      
      for (i in 1:(ncol(a)-length(my.sel.cancers))) { 
         #a[,i]<-a[,i]/max(a[,i])  # don't do a rescaling
         a[,i]<-a[,i]/sum(a[,i])   # put the proportions
      }
      for (i in (ncol(a)-length(my.sel.cancers)+1):ncol(a)) { 
         a[,i]<-a[,i]+1.5   # put the proportions   # add 1.5 to cancers
      }

      rownames(a)<-colnames(cancer_signatures)[1:30] 
      colorends <- c("white","red", "white", "blue")
      dendro <- "none"
      if (input$row_c_heatmap=="yes") dendro<-"row"
      if (input$col_c_heatmap=="yes") dendro<-"column" 
      if (input$row_c_heatmap=="yes" & input$col_c_heatmap=="yes") dendro<-"both"

            heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,3)),
                dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90)
      
   })
   

   #  Download HeatMap 
   #  output$download_known_plot <- downloadHandler(filename = function(){paste("comparison_with_other",input$type_known_plot, sep=".")}, content=function (ff) {
   #    
   #    if ("All" %in% input$mycancers) my.sel.cancers<-colnames(known_cancer_signatures)
   #    else my.sel.cancers<-intersect(input$mycancers,colnames(known_cancer_signatures))
   #    
   #    
   #    a<-t(data.frame(my_contributions()[30:1,], known_cancer_signatures[30:1,my.sel.cancers]))
   #    colnames(a)<-colnames(cancer_signatures)[30:1]
   #    if (ncol(my_contributions())==1) rownames(a)[1]<-colnames(my_contributions()) ## fix colnames when there is only one sample
   # 
   #    for (i in 1:(nrow(a)-length(my.sel.cancers))) { 
   #    #   a[i,]<-a[i,]/max(a[i,])   # don't do a rescaling
   #        a[i,]<-a[i,]/sum(a[i,])   # put the proportions
   #    }
   #    a.m<-reshape2::melt(as.matrix(a)) 
   #    a.m$category<-rep(c(rep("Sample",nrow(a)-length(my.sel.cancers)),rep("Cancers",length(my.sel.cancers))),30)
   #    sel<-which(a.m$category=="Cancers")
   #    a.m[sel,"value"]<-a.m[sel,"value"]+1.5
   #    a.m[is.na(a.m)] <- 0
   #    
   #    colorends <- c("white","red", "white", "blue")
   #    
   #    ggplot(a.m, aes(x=Var1, y=Var2)) + geom_tile(aes(fill = value),
   #                                                 colour = "white") + theme(axis.text.x=element_text(angle=90)) +
   #       scale_fill_gradientn(colours = colorends, limits = c(0,3)) + labs(x="",y="")
   #    if (input$type_known_plot=="pdf") ggsave(ff)
   #    if (input$type_known_plot=="png") ggsave(ff)
   #    if (input$type_known_plot=="tiff") ggsave(ff,compression="lzw")
   #    
   # })
   
   
   ###### PCA - Clustering of samples ## only if there are 3 or more samples
   
   output$pca_plot <- renderPlot({
      
      if (ncol(as.data.frame(my_contributions()))>=3) {
         
      a<-t(as.data.frame(my_contributions()[30:1,]))
      for (i in 1:nrow(a)) { 
         a[i,]<-a[i,]/sum(a[i,])   # put the proportions
      }
      a<-a[,which(apply(a,2,sd)>0)] # remove signatures without variation
      pca <- prcomp(a, scale=T)
      plot(pca$x[,1], pca$x[,2],        # x y and z axis
           col="red", pch=19,  
           xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
           ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
           main="PCA")
      text(pca$x[,1], pca$x[,2], rownames(a))
      
      } else {
         par(mar = c(0,0,0,0))
         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
         text(x = 0.5, y = 0.5, paste("PCA analysis works only with >=3 samples"), 
              cex = 1.6, col = "black")
      }

   })
   
   
   output$download_pca_plot <- downloadHandler (
      filename = function(){paste("pca_plot",input$type_pca_plot, sep=".")}, 
      content = function(ff) {
         if (input$type_pca_plot=="pdf") pdf(ff,height=7,width=7)
         if (input$type_pca_plot=="png") png(ff,height=7*ppi,width=7*ppi,res=ppi)
         if (input$type_pca_plot=="tiff") tiff(ff,height=7*ppi,width=7*ppi,res=ppi,compression="lzw")

         if (ncol(as.data.frame(my_contributions()))>=3) {
             a<-t(as.data.frame(my_contributions()[30:1,]))
             for (i in 1:nrow(a)) {
                a[i,]<-a[i,]/sum(a[i,])   # put the proportions
             }
             a<-a[,which(apply(a,2,sd)>0)] # remove signatures without variation
             pca <- prcomp(a, scale=T)
             plot(pca$x[,1], pca$x[,2],        # x y and z axis
                  col="red", pch=19,
                  xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                  ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                  main="PCA")
             text(pca$x[,1], pca$x[,2], rownames(a))
         } else {
             par(mar = c(0,0,0,0))
             plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
             text(x = 0.5, y = 0.5, paste("PCA analysis works only with >=3 samples"),
                  cex = 1.6, col = "black")
         }

         dev.off()
         
      })
   
   
   
})
