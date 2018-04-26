library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(V8)
library(shinythemes)
library(plotly)
#webshot::install_phantomjs()  #To plot into pdf Heatmaply

shinyServer(function(input, output,session){
	

	
   #Setting maximum file size for uploading (1000 MB)
   options(shiny.maxRequestSize=1000*1024^2)
   options(shiny.sanitize.errors = TRUE)
   
   #Resolution of the tiff images
   ppi<-200 
   
   
   #Hidding/Showing tabs of mainpanel, run and clear buttons, ...
   observeEvent(input$run,{
   	shinyjs::hide(id="run")
   	shinyjs::show(id="after_run")
   	
   	
   	if (input$datatype=="TSV" & ((length(grep(".tsv",input[["fileinput"]]$datapath)) + length(grep(".txt",input[["fileinput"]]$datapath)))==length(input[["fileinput"]]$datapath))){
   		
   		mat_list<-list()
   		for (w in 1:length(input[["fileinput"]]$datapath)){
   			mat <- fread(input[["fileinput"]]$datapath[w],header=T,sep="\t",data.table=F)
   			condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)
   			
   			mat_list[[w]] <- condition
   		}
   		mat_vector <- as.vector(do.call(cbind, mat_list))
   	} else {
   		mat_vector<-c("mat_vector")
   	}
   	
   	if (input$datatype=="Excel" & (length(grep(".xls",input[["fileinput"]]$datapath))==length(input[["fileinput"]]$datapath))){
   		mat_list<-list()
   		for (w in 1:length(input[["fileinput"]]$datapath)){
   			
   			if (length(grep(".xlsx",input[["fileinput"]]$datapath[w]))>0){
   				mat<-read.xlsx(input[["fileinput"]]$datapath[w],1)
   			} else {
   				mat<-data.frame(read_xls(input[["fileinput"]]$datapath[w],col_names = TRUE,col_types = "text"))
   			}
   			
   			condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)
   			
   			mat_list[[w]] <- condition
   		}
   		mat_vector_2 <- as.vector(do.call(cbind, mat_list))
   	} else {
   		mat_vector_2<-c("mat_vector_2")
   	}
   		
   	
   	if (length(grep("FALSE", mat_vector))!=0 |
   		 length(grep("FALSE", mat_vector_2))!=0 |
   		(input$datatype=="VCF" & length(grep(".vcf",input[["fileinput"]]$datapath))!=length(input[["fileinput"]]$datapath)) |
		    (input$datatype=="TSV" & ((length(grep(".tsv",input[["fileinput"]]$datapath)) + length(grep(".txt",input[["fileinput"]]$datapath)))!=length(input[["fileinput"]]$datapath))) |
   		 (input$datatype=="Excel" & (length(grep(".xls",input[["fileinput"]]$datapath))!=length(input[["fileinput"]]$datapath))) |
   		 (input$datatype=="MAF" & length(grep(".maf",input[["fileinput"]]$datapath))<=0) |
   		 (input$datatype=="MAF" & length(input[["fileinput"]]$datapath)!=1)
   		 ){}else{
      shinyjs::show(id="mainpanel")
   	}
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
   library(heatmaply)
   library(gplots)
   library(openxlsx)
   library(readxl)

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
 
   
   filedata <- reactive({
      infile <- input$fileinput #file loaded
      if (is.null(infile)) {
         # User has not uploaded a file yet
         return(NULL)
      } else {return(infile)}
   })
   
   output$run_button <- renderUI({
      if (!is.null(filedata())) {
         
         
         div(
            #Run button
            actionButton("run","Run",class = "btn-primary"),
            
            
            #Busy indicator
            busyIndicator("Running",wait=0)
            
         )
         
      } else if (is.null(filedata())) {
         HTML("")
      }
   })

   
   #######WHEN RUN BUTTON##
   observeEvent(input$run,{
     all_ok_for_run <- TRUE
     
     
     if (!is.null(filedata())) {
       output$error_input <- renderUI({
         HTML("All well!")
       })
     } else {
       all_ok_for_run <- FALSE
     }
       



   #######################################
   #Reading input files as GRanges objects [vcfs]
   #######################################
   vcfs<-eventReactive((all_ok_for_run==TRUE),{

      inFile<<-input$fileinput
            
         #VCF
         if (input$datatype=="VCF"){
            
         #Read vcf for MutationalPatterns
            return(read_vcfs_as_granges(inFile$datapath,inFile$name,ref_genome(),group = "auto+sex", check_alleles = TRUE))
         }
      
      
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
            
            ff_list<-list()
            for (w in 1:length(inFile$datapath)){
               
               if (length(grep(".xlsx",inFile$datapath[w]))>0){
                  aux<-read.xlsx(inFile$datapath[w],1)
               } else {
                  aux<-data.frame(read_xls(inFile$datapath[w],col_names = TRUE,col_types = "text"))
               }
               
               
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
   cancer_signatures_aux <- cancer_signatures[order(cancer_signatures[,1]),]
   cancer_signatures <- as.matrix(cancer_signatures_aux[,4:33])
   cancer_signatures_mut_types <- as.matrix(cancer_signatures_aux[,1:3])
   
   
   #######################################
   #Fitting mutations in samples (mut_mat) to COSMIC signatures [fit_res]
   #######################################
   fit_res <- reactive({ fit_to_signatures(mut_mat(), cancer_signatures) })
         
   
   #Auxiliar files of aetiology and known signatures by cancer type (from COSMIC)
   proposed_etiology <- fread("./aux_files/proposed_etiology_COSMIC_signatures.txt",sep="\t",header=F,data.table=F)[,2]
   known_cancer_signatures<-read.table("./aux_files/cancermatrix.tsv",header=TRUE,sep="\t",row.names=1)
   
   
   #divisionRel function creation to print final dataframe
   divisionRel<-function(df){
      sum_df<-sapply(df,sum)
      for (i in 1:ncol(df)){
         df[,i]<-round((df[,i]/sum_df[i]),3)
      }
      return(df)
   }
   

   ################################################################
   output$custom_error<-renderUI({
   	
   	if (input$datatype=="TSV" & ((length(grep(".tsv",input[["fileinput"]]$datapath)) + length(grep(".txt",input[["fileinput"]]$datapath))) == length(input[["fileinput"]]$datapath))){
   		
   		mat_list<-list()
   		for (w in 1:length(input[["fileinput"]]$datapath)){
   			mat <- fread(input[["fileinput"]]$datapath[w],header=T,sep="\t",data.table=F)
   			condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)
   			
   			mat_list[[w]] <- condition
   		}
   		mat_vector <- as.vector(do.call(cbind, mat_list))
   		samples_failed<-grep("FALSE", mat_vector)
   	} else {
   		mat_vector<-c("mat_vector")
   		
   		if (input$datatype=="Excel" & (length(grep(".xls",input[["fileinput"]]$datapath))==length(input[["fileinput"]]$datapath))){
   			mat_list<-list()
   			for (w in 1:length(input[["fileinput"]]$datapath)){
   				
   				if (length(grep(".xlsx",input[["fileinput"]]$datapath[w]))>0){
   					mat<-read.xlsx(input[["fileinput"]]$datapath[w],1)
   				} else {
   					mat<-data.frame(read_xls(input[["fileinput"]]$datapath[w],col_names = TRUE,col_types = "text"))
   				}
   				
   				condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)
   				
   				mat_list[[w]] <- condition
   			}
   			mat_vector_2 <- as.vector(do.call(cbind, mat_list))
   			samples_failed<-grep("FALSE", mat_vector_2)
   		} else {
   			mat_vector_2<-c("mat_vector_2")
   			samples_failed<-NULL
   		}
   		
   	}
   	
   	
   	
   	if (input$datatype=="VCF" & length(grep(".vcf",input[["fileinput"]]$datapath))!=length(input[["fileinput"]]$datapath)){
   		aux<-grep(".vcf",input[["fileinput"]]$datapath,invert=TRUE)
   		samples_bad_format<-paste(input[["fileinput"]]$name[aux], collapse=" | ")
   	} else {
   		
   		if (input$datatype=="TSV" & (length(grep(".tsv",input[["fileinput"]]$datapath)) + length(grep(".txt",input[["fileinput"]]$datapath))) != length(input[["fileinput"]]$datapath)){
   			aux1<-grep(".tsv",input[["fileinput"]]$datapath,invert=TRUE)
   			aux2<-input[["fileinput"]]$name[aux1]
   			aux3<-grep(".txt",aux2,invert=TRUE)
   			samples_bad_format<-paste(aux2[aux3], collapse=" | ")
   		}
   		else{
   			
   			if (input$datatype=="Excel" & length(grep(".xls",input[["fileinput"]]$datapath)) != length(input[["fileinput"]]$datapath)){
   				aux<-grep(".xls",input[["fileinput"]]$datapath,invert=TRUE)
   				samples_bad_format<-paste(input[["fileinput"]]$name[aux], collapse=" | ")
   			} else {
   				if (input$datatype=="MAF" & length(grep(".maf",input[["fileinput"]]$datapath))<=0){
   					aux<-grep(".maf",input[["fileinput"]]$datapath,invert=TRUE)
   					samples_bad_format<-paste(input[["fileinput"]]$name[aux], collapse=" | ")
   				} else {
   					samples_bad_format<-NULL
   				}
   			}
   		}
   	}
   	
   	
   	

   	tags$style(HTML(paste(".shiny-output-error-formats {visibility: hidden;}
   						 .shiny-output-error-formats:before {visibility: visible; color: orangered; content:'Format error in uploaded file/s:  ", samples_bad_format,"  Please select the correct input file format before uploading your file/s.';}
   						 .shiny-output-error-maf {visibility: hidden;}
   						 .shiny-output-error-maf:before {visibility: visible; color: orangered; content:'Only one multi-sample MAF file is allowed.';}
   						 .shiny-output-error-noheader {visibility: hidden;}
   						 .shiny-output-error-noheader:before {visibility: visible; color: orangered; content:'Uploaded file/s:  ", paste(input[["fileinput"]]$name[samples_failed],collapse=" | "),"  do not have the mandatory header with columns CHROM, POS, REF and ALT.';}")
   					)
   	)
   })
   ################################################################

   
   #Plot selectize to select samples to plot.
   output$selected_samples<-renderUI({
      
      #Error managemente for file format
   	error_message<-"File format error, please select the correct input file format before uploading your file/s."
   	
      if (input$datatype=="VCF"){
         validate(
            need(length(grep(".vcf",input[["fileinput"]]$datapath))==length(input[["fileinput"]]$datapath),error_message), errorClass = "formats"
         )
      }
      if (input$datatype=="TSV"){
         validate(
            need((length(grep(".tsv",input[["fileinput"]]$datapath)) + length(grep(".txt",input[["fileinput"]]$datapath))) == length(input[["fileinput"]]$datapath),error_message), errorClass = "formats"
         )
      }
      if (input$datatype=="TSV"){

       mat_list<-list()
       for (w in 1:length(input[["fileinput"]]$datapath)){
         mat <- fread(input[["fileinput"]]$datapath[w],header=T,sep="\t",data.table=F)
         condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)

         mat_list[[w]] <- condition
       }
       mat_vector <- as.vector(do.call(cbind, mat_list))
       validate(
         need(length(grep("FALSE", mat_vector))==0, "Uploaded files header do not have necessary columns CHROM, POS, REF and ALT."), errorClass = "noheader"
       )

      }
      if (input$datatype=="Excel"){
         validate(
            need(length(grep(".xls",input[["fileinput"]]$datapath)) == length(input[["fileinput"]]$datapath),error_message), errorClass = "formats"
         )
      }
   	if (input$datatype=="Excel"){
   		   mat_list<-list()
   		   for (w in 1:length(input[["fileinput"]]$datapath)){

   		      if (length(grep(".xlsx",input[["fileinput"]]$datapath[w]))>0){
   		         mat<-read.xlsx(input[["fileinput"]]$datapath[w],1)
   		      } else {
   		         mat<-data.frame(read_xls(input[["fileinput"]]$datapath[w],col_names = TRUE,col_types = "text"))
   		      }

   		      condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)

   		      mat_list[[w]] <- condition
   		   }
   		   mat_vector <- as.vector(do.call(cbind, mat_list))
   		   validate(
   		      need(length(grep("FALSE", mat_vector))==0, "Uploaded files header do not have necessary columns CHROM, POS, REF and ALT."), errorClass = "noheader"
   		   )

   	}

   	
      if (input$datatype=="MAF"){
         validate(
            need(length(grep(".maf",input[["fileinput"]]$datapath))>0,error_message),errorClass = "formats"
         )
      	validate(
            need(length(input[["fileinput"]]$datapath)==1, "Only one multi-sample MAF file is allowed."),errorClass = "maf"
         )
      }

            
      
      if (length(vcfs())==1){
            mysamp<-colnames(as.data.frame(fit_res()$contribution))
            selectInput("mysamp","Select your samples",mysamp, multiple=TRUE, selectize = FALSE, size=1, selected=colnames(as.data.frame(fit_res()$contribution)))
      } else {
      
         if (input$tab=="reconst"){
            mysamp<-colnames(as.data.frame(fit_res()$contribution))
            selectInput("mysamp","Select your samples",mysamp, multiple=TRUE, selectize = FALSE, size=6, selected = mysamp[1])
         } else {
            mysamp<-c("All samples",colnames(as.data.frame(fit_res()$contribution)))
            selectInput("mysamp","Select your samples",mysamp, multiple=TRUE, selectize = FALSE, size=6, selected="All samples")
         }
      }
   })
   
   output$selected_cancer_types<-renderUI({
      
      if (input$tab=="comp_canc_sign"){
         
         selectInput("mycancers","Select the cancers to compare", c("All cancers",colnames(read.table("./aux_files/cancermatrix.tsv",header=TRUE,sep="\t",row.names=1))), multiple=TRUE, selectize=FALSE, size=10, selected="All cancers")
      
      }
      
   })
 
   #Checkbox to decide if mean is plotted.
   output$mean_checkbox<-renderUI({
      if (input$tab=="smp" | input$tab=="contrib" | input$tab=="comp_canc_sign"){
         checkboxInput("meancheck","Show mean",value=FALSE)
      } else {
         return(invisible(NULL))
      }
   })

   
   
   
   #Select which samples use to plot.
   my_contributions<- reactive({ 
      
      #Error management
      if (length(input$mysamp)==0) return(invisible(NULL))
         

      if ("All samples" %in% input$mysamp){
         if(input$meancheck==TRUE) {
            aux<-divisionRel(as.data.frame(fit_res()$contribution))
            con<-data.frame(aux, mean = apply(aux,1,mean))
            colnames(con)<-c(colnames(aux),"mean")
         } else {
            aux<-divisionRel(as.data.frame(fit_res()$contribution))
            con<-data.frame(aux)
            colnames(con)<-c(colnames(aux))   
         }
         
      } else {
         if(input$meancheck==TRUE) {      
            aux<-divisionRel(as.data.frame(fit_res()$contribution[,input$mysamp]))
            con<-data.frame(aux, mean = apply(aux,1,mean))
            colnames(con)<-c(colnames(aux),"mean")
         } else {
            aux<-divisionRel(as.data.frame(fit_res()$contribution[,input$mysamp]))
            con<-data.frame(aux)
            colnames(con)<-colnames(aux)
         }
      }
      
            
      
      #Fixing colname of one sample (without mean)
      if (ncol(con)==1){
         colnames(con)<-setdiff(input$mysamp,c("All samples"))
      }
         
      #Fixing colname of one sample (with mean)
      if (ncol(con)==1 & input$meancheck==TRUE & input$mysamp!="All samples"){
         colnames(con)[1]<-setdiff(input$mysamp,c("All samples"))
      }
         
         
      return(con)
      
   })


   #######################################
   #PLOT Somatic Mutation Prevalence (number mutations per megabase)
   #######################################
   
   #Selection of type of study and MB affected by it
   output$kb_sequenced<-renderUI({
      
      if (input$studytype == "Targeted Sequencing"){
         numericInput("bases_sequenced","Kilobases sequenced",value = 10,min = 1)
      }
         
   })
   
   #According to Alexandrov et al. 2013
   megabases<-reactive({
      if (input$studytype == "Whole Genome Sequencing"){
         return(2800)
      }
      
      if (input$studytype == "Whole Exome Sequencing"){
         return(30)
      }
      
      if (input$studytype == "Targeted Sequencing"){
         return(input$bases_sequenced/1000)
      }
   })
   
   
   #Selection of samples to plot
   mutation_counts<- reactive({ 
         
      if ("All samples" %in% input$mysamp){
         
         if(input$meancheck==TRUE) {
            mc<-data.frame(samples=c(names(vcfs()),"mean"),smp=(c(sapply(vcfs(),length),mean(sapply(vcfs(),length))))/megabases())
         } else {
            mc<-data.frame(samples=names(vcfs()),smp=(sapply(vcfs(),length))/megabases())
         }
         
      } else {
         
         if(input$meancheck==TRUE) {      
            aux<-input$mysamp
            mc<-data.frame(samples=c(names(vcfs()[aux]),"mean"), smp=(c(sapply(vcfs()[aux],length),mean(sapply(vcfs()[aux],length))))/megabases())
         } else {
            mc<-data.frame(samples=names(vcfs()[input$mysamp]), smp=(sapply(vcfs()[input$mysamp],length))/megabases())
         }
      }
      return(mc)
   })
   
       
#   PLOT somatic mutation prevalence
   output$smp <- renderPlot({

      #Error management
      if (length(input$mysamp)==0) return(invisible(NULL))
      

      mutation_counts_new<-data.frame(samples=mutation_counts()$samples,smp=round(mutation_counts()$smp,1))

      plot_smp<-ggplot(data=mutation_counts_new,aes(x=samples,y=smp))

      plot_smp + geom_bar(stat="identity",fill="orangered2") + theme_minimal() + geom_text(aes(label=smp), size=5, position = position_stack(vjust = 0.5), colour="white") + coord_flip() + labs(x = "", y = "Somatic mutation prevalence (number of mutations per megabase)") + theme(axis.text=element_text(size=12), axis.title = element_text(size = 13, face = "bold"), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

   }, height = function(x=length(colnames(my_contributions())) ){ if (x < 5) {return(250)} 
     if (x >100) {return(10000)}   
     if (x<=100 & x>=5) {return (50*x)}   
   })
   

   #Download Plot somatic mutation prevalence 
   output$download_smp_plot <- downloadHandler (
      filename = function(){
         paste("mut_prevalence_plot",input$type_smp_plot, sep=".")
      },
      content = function(ff) {
         mutation_counts_new<-data.frame(samples=mutation_counts()$samples,smp=round(mutation_counts()$smp,1))
         
         plot_smp<-ggplot(data=mutation_counts_new,aes(x=samples,y=smp))
         
         plot_smp + geom_bar(stat="identity",fill="orangered2") + theme_minimal() + geom_text(aes(label=smp), size=5, position = position_stack(vjust = 0.5), colour="white") + coord_flip() + labs(x = "", y = "Somatic mutation prevalence (number of mutations per megabase)") + theme(axis.text=element_text(size=12), axis.title = element_text(size = 13, face = "bold"), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())       
         
         ggsave(ff,height=min(2*nrow(mutation_counts_new),40),width=25,dpi=ppi,units="cm")

      }
   )
   
   #Download Table somatic mutation prevalence 
   output$download_smp_table <- downloadHandler(
      filename="Somatic_mut_prev.txt",
      content=function (file){
         mutation_counts_new<-data.frame(Sample=mutation_counts()$samples,Somatic_Mutation_Prevalence=round(mutation_counts()$smp,1),Number_of_Samples=length(names(vcfs())))
         
         write.table(x = mutation_counts_new, file = file, sep = "\t", quote=F, row.names=F)
      }
   )
   
   

   #######################################
   #PLOT 96 nucleotide changes profile (samples individually)
   #######################################
   
   #Plot 96 profile
   output$prof96 <- renderPlot({
      aux_96_profile<-as.matrix(mut_mat()[,setdiff(colnames(my_contributions()),c("mean"))])
      colnames(aux_96_profile)<-setdiff(colnames(my_contributions()),c("mean"))
      
      aux_ymax<-as.data.frame(aux_96_profile)
      rownames(aux_ymax)<-1:96
      max_ymax<-max(divisionRel(aux_ymax))
      
      plot_96_profile(aux_96_profile,ymax = max_ymax) + scale_y_continuous(breaks = seq(0, max_ymax, 0.05))
   }, height = function(x=length(colnames(my_contributions())) ){ if (x < 5) {return(400)} 
     if (x >100) {return(10000)}   
     if (x<=100 & x>=5) {return (20+100*x)}   
      }
  )

   
   #Download Plot 96 profile 
   output$download_prof96_plot <- downloadHandler (
      filename = function(){
         paste("prof96_plot",input$type_prof96_plot, sep=".")
      },
      content = function(ff) {
         aux_96_profile<-as.matrix(mut_mat()[,setdiff(colnames(my_contributions()),c("mean"))])
         colnames(aux_96_profile)<-setdiff(colnames(my_contributions()),c("mean"))
         
         aux_ymax<-as.data.frame(aux_96_profile)
         rownames(aux_ymax)<-1:96
         max_ymax<-max(divisionRel(aux_ymax))
         
         plot_96_profile(aux_96_profile,ymax = max_ymax) + scale_y_continuous(breaks = seq(0, max_ymax, 0.05))
         
         ggsave(ff,height=min(4*ncol(aux_96_profile),40),width=25,dpi=ppi,units="cm")
      }
   )
      
   
   #Download Plot 96 TABLE
   
   
   output$download_prof96_table <- downloadHandler(
      filename="Mutational_Profile.txt",
      content=function (file){
         aux_96_profile<-as.matrix(mut_mat()[,setdiff(colnames(my_contributions()),c("mean"))])
         aux_ymax<-divisionRel(as.data.frame(aux_96_profile))
         
         
         write.table(x = data.frame(Substitution_Type = cancer_signatures_mut_types[,1], Trinucleotide = cancer_signatures_mut_types[,2], Somatic_Mutation_Type = cancer_signatures_mut_types[,3], aux_ymax), file = file, sep = "\t", quote=F, row.names=F)
      }
   )
   
   
   
      
   #######################################
   ### Plot heatmap with contributions
   #######################################
   
   
   #DataTable
   output$contr <- renderDataTable({
         data.frame(Signature = 1:30, Proposed_Etiology = proposed_etiology, round(my_contributions(),3))
      },
      options = list(lengthChange=FALSE,pageLength=30, paging=FALSE, searching=FALSE, info=FALSE)
   )
   
   
   #Download Table
   output$download_contr <- downloadHandler( filename="COSMIC_sign_contributions.txt", content=function (file){
      write.table(x = data.frame(Signature = 1:30, Proposed_Etiology = proposed_etiology, round(my_contributions(),3)), file = file, sep = "\t", quote=F, row.names=F)})
   
   
   #check if column or row dendogram is needed
   output$row_dendro_heatmap<-renderUI({
      if (input$tab == "contrib"){
         radioButtons("row_d_heatmap", "Row dendrogram", c("yes","no"),selected = "no",inline = TRUE)
      } else {
         return(invisible(NULL))
      }
   })
   output$col_dendro_heatmap<-renderUI({
      if (input$tab == "contrib"){
         radioButtons("col_d_heatmap", "Column dendrogram", c("yes","no"),selected = "no",inline = TRUE)
      } else {
         return(invisible(NULL))
      }
   })
   
   #check if sample names are needed
   output$heatmap_sample_names<-renderUI({
      if (input$tab == "contrib" | input$tab == "comp_canc_sign"){
         checkboxInput("samplenames","Show sample names",value=TRUE)
      } else {
         return(invisible(NULL))
      }
   })

   #HeatMap
   output$heatmap_signatures <- renderPlotly({
      
      if (length(input$row_d_heatmap)==0) return(invisible(NULL))
      if (length(input$col_d_heatmap)==0) return(invisible(NULL))
       
      a<-my_contributions()
      if (ncol(a)==1) colnames(a)<-colnames(my_contributions()) ## fix colnames when there is only one sample
      
      
      rownames(a)<-colnames(cancer_signatures)[1:30] 
      colorends <- c("white","red")
      dendro <- "none"
      if (input$row_d_heatmap=="yes") dendro<-"row"
      if (input$col_d_heatmap=="yes") dendro<-"column" 
      if (input$row_d_heatmap=="yes" & input$col_d_heatmap=="yes") dendro<-"both"
   
      if (input$samplenames==FALSE){
         colnames(a)<-1:length(colnames(a))
      }
      
      heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,1)),
                dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90, distfun = 'pearson'  )
   })
   
   
   #Download HeatMap 
    output$download_signatures_plot <- downloadHandler (
       filename = function(){paste("signatures_plot",input$type_signatures_plot, sep=".")}, 
       content = function(ff) {
          
          
          a<-my_contributions()
          if (ncol(a)==1) colnames(a)<-colnames(my_contributions()) ## fix colnames when there is only one sample
          rownames(a)<-colnames(cancer_signatures)[1:30] 
          colorends <- c("white","red")
          dendro <- "none"
          if (input$row_d_heatmap=="yes") dendro<-"row"
          if (input$col_d_heatmap=="yes") dendro<-"column" 
          if (input$row_d_heatmap=="yes" & input$col_d_heatmap=="yes") dendro<-"both"
          
          if (input$samplenames==FALSE){
             colnames(a)<-1:length(colnames(a))
          }
          
          heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,1)),
                    dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90, distfun = 'pearson',
                    file = ff)

       })
       
          

    #######################################
    ### PLOT Reconstructed Mutational Profile
    #######################################
    
    #Plot reconstructed profile
    output$reconst <- renderPlot({
      validate(
          need(length(input$mysamp)==1,"Sample selection error, please select just one sample at a time to visualize its reconstructed mutational profile."), errorClass = "reconstructed"
      )
       
      if (input$mysamp=="All samples") return(invisible(NULL))
       
       
      original_prof <- mut_mat()[,input$mysamp]
      reconstructed_prof <- fit_res()$reconstructed[,input$mysamp]
      
      aux_ymax<-data.frame(original_prof,reconstructed_prof)
      max_ymax<-max(divisionRel(aux_ymax))
      
      
      plot_compare_profiles(original_prof,reconstructed_prof,profile_names=c("Original","Reconstructed"),profile_ymax = max_ymax)
      
    })
    
    
    
    #Download Plot reconstructed profile 
    output$download_reconst_plot <- downloadHandler (
       filename = function(){
          paste("reconstructed_plot",input$type_reconst_plot, sep=".")
       },
       content = function(ff) {
          
          
          original_prof <- mut_mat()[,input$mysamp]
          reconstructed_prof <- fit_res()$reconstructed[,input$mysamp]
          
          aux_ymax<-data.frame(original_prof,reconstructed_prof)
          max_ymax<-max(divisionRel(aux_ymax))
          

          
          plot_compare_profiles(original_prof,reconstructed_prof,profile_names=c("Original","Reconstructed"),profile_ymax = max_ymax)
          
          
          
          ggsave(ff,height=6,width=10,dpi=ppi)
       }
    )
    
    #Download reconstructed TABLE
    

    
    output$download_reconst_table <- downloadHandler(
       filename="reconstructed_table.txt",
       content=function (file){
          original_prof <- mut_mat()[,input$mysamp]
          reconstructed_prof <- fit_res()$reconstructed[,input$mysamp]
          
          aux_ymax<-divisionRel(data.frame(Original_Profile = original_prof, Reconstructed_Profile = reconstructed_prof))
          
          
          write.table(x = data.frame(Substitution_Type = cancer_signatures_mut_types[,1], Trinucleotide = cancer_signatures_mut_types[,2], Somatic_Mutation_Type = cancer_signatures_mut_types[,3], aux_ymax), file = file, sep = "\t", quote=F, row.names=F)
       }
    )
   
    

   
   #######################################
   ### Plot - Comparison with other cancers
   #######################################
   
   #check if column or row dendogram is needed
   output$row_dendro_cancers<-renderUI({
      if (input$tab == "comp_canc_sign"){
         radioButtons("row_c_heatmap", "Row dendrogram", c("yes","no"),selected = "no",inline = TRUE)
      } else {
         return(invisible(NULL))
      }
   })
   output$col_dendro_cancers<-renderUI({
      if (input$tab == "comp_canc_sign"){
         radioButtons("col_c_heatmap", "Column dendrogram", c("yes","no"),selected = "no",inline = TRUE)
      } else {
         return(invisible(NULL))
      }
   })
   
   #HeatMap
   output$heatmap_known <- renderPlotly({
      
      if (length(input$mysamp)==0) return(invisible(NULL))
      if (length(input$mycancers)==0) return(invisible(NULL))
      
      if ("All cancers" %in% input$mycancers) my.sel.cancers<-colnames(known_cancer_signatures)
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
         a[,i]<-a[,i]*2.5   # put the proportions   # 1 goes to 2.5 (light blue)

      }

      rownames(a)<-colnames(cancer_signatures)[1:30] 
      colorends <- c("white","red", "white", "blue")
      dendro <- "none"
      if (input$row_c_heatmap=="yes") dendro<-"row"
      if (input$col_c_heatmap=="yes") dendro<-"column" 
      if (input$row_c_heatmap=="yes" & input$col_c_heatmap=="yes") dendro<-"both"

      if (input$samplenames==FALSE){
         colnames(a)<-1:length(colnames(a))
      }
      
      heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,3)),
                dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90, hide_colorbar = TRUE,
                distfun = 'pearson')
      
   })
   

   #  Download HeatMap 
      output$download_known_plot <- downloadHandler(filename = function(){paste("comparison_with_cancers",input$type_known_plot, sep=".")}, content=function (ff) {
        
        
           if ("All cancers" %in% input$mycancers) my.sel.cancers<-colnames(known_cancer_signatures)
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
           
           if (input$samplenames==FALSE){
              colnames(a)<-1:length(colnames(a))
           }
           
           heatmaply(a, scale_fill_gradient_fun = scale_fill_gradientn(colours = colorends, limits = c(0,3)),
                     dendrogram = dendro, k_row = 1, k_col = 1, column_text_angle = 90, hide_colorbar = TRUE,
                     distfun = 'pearson',
                     file = ff)
        
      })
   
      
      
   #######################################
   ###### PCA - Clustering of samples ## only if there are 3 or more samples
   #######################################
   output$pca_plot <- renderPlot({
      
      #Error management
      validate(
         need((length(input$mysamp)>2 | "All samples" %in% input$mysamp ),"PCA analysis works only with 3 or more samples."),errorClass = "pca"
      )
      
      
      #######################################
      my_contributions_mod <- my_contributions()
      #######################################
      
      #Error management
      validate(
         need((length(colnames(my_contributions_mod))>2),"PCA analysis works only with 3 or more samples."),errorClass = "pca"
      )

      if (ncol(as.data.frame(my_contributions_mod))>=3) {
         a<-t(as.data.frame(my_contributions_mod[30:1,]))
         for (i in 1:nrow(a)) { 
           a[i,]<-a[i,]/sum(a[i,])   # put the proportions
         }
         a<-a[,which(apply(a,2,sd)>0)]# remove signatures without variation
         
         samplesnames<-rownames(a)
         rownames(a)<-1:(length(rownames(a)))
         
         pca <- prcomp(a, scale=T)
         plot(pca$x[,1], pca$x[,2],        # x y and z axis
              col="red", pch=19,  
              xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
              ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
              xlim=c(min(pca$x[,1])-0.5*(  max(pca$x[,1])-min(pca$x[,1]) ) ,max(pca$x[,1])+0.5*(  max(pca$x[,1])-min(pca$x[,1]) ) ),
              ylim=c(min(pca$x[,2])-0.5*(  max(pca$x[,2])-min(pca$x[,2]) ) ,max(pca$x[,2])+0.5*(  max(pca$x[,2])-min(pca$x[,2]) ) ),
              main="PCA")
         text(pca$x[,1], pca$x[,2]-0.15, rownames(a))

           
      } else {
         par(mar = c(0,0,0,0))
         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
         text(x = 0.5, y = 0.5, paste("PCA analysis works only with >=3 samples"), 
              cex = 1.6, col = "black")
      }

   })
   
   
   output$pca_plot_table <- renderTable({
      
      my_contributions_mod <- my_contributions()
      
      
      if (ncol(as.data.frame(my_contributions_mod))>=3) {
         a<-t(as.data.frame(my_contributions_mod[30:1,]))
         for (i in 1:nrow(a)) { 
            a[i,]<-a[i,]/sum(a[i,])   # put the proportions
         }
         a<-a[,which(apply(a,2,sd)>0)]# remove signatures without variation
         
         samplesnames<-rownames(a)
         rownames(a)<-1:(length(rownames(a)))
         
         data.frame(ID=rownames(a),Sample=samplesnames)
         
         
      }
      
      
   })
   
   
   
   
   output$download_pca_plot <- downloadHandler (
      filename = function(){paste("pca_plot",input$type_pca_plot, sep=".")}, 
      content = function(ff) {
         if (input$type_pca_plot=="pdf") pdf(ff,height=7,width=7)
         if (input$type_pca_plot=="png") png(ff,height=7*ppi,width=7*ppi,res=ppi)
         if (input$type_pca_plot=="tiff") tiff(ff,height=7*ppi,width=7*ppi,res=ppi,compression="lzw")

         
         #######################################
         my_contributions_mod <- my_contributions()
         #######################################
         
         if (ncol(as.data.frame(my_contributions_mod))>=3) {
            a<-t(as.data.frame(my_contributions_mod[30:1,]))
            for (i in 1:nrow(a)) { 
               a[i,]<-a[i,]/sum(a[i,])   # put the proportions
            }
            a<-a[,which(apply(a,2,sd)>0)]# remove signatures without variation
            
            samplesnames<-rownames(a)
            rownames(a)<-1:(length(rownames(a)))
            
            pca <- prcomp(a, scale=T)
            plot(pca$x[,1], pca$x[,2],        # x y and z axis
                 col="red", pch=19,  
                 xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                 ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                 xlim=c(min(pca$x[,1])-0.5*(  max(pca$x[,1])-min(pca$x[,1]) ) ,max(pca$x[,1])+0.5*(  max(pca$x[,1])-min(pca$x[,1]) ) ),
                 ylim=c(min(pca$x[,2])-0.5*(  max(pca$x[,2])-min(pca$x[,2]) ) ,max(pca$x[,2])+0.5*(  max(pca$x[,2])-min(pca$x[,2]) ) ),
                 main="PCA")
            text(pca$x[,1], pca$x[,2]-0.25, rownames(a))
            
            
         } else {
            par(mar = c(0,0,0,0))
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, paste("PCA analysis works only with >=3 samples"), 
                 cex = 1.6, col = "black")
         }

         dev.off()
         
      })
   
   output$download_pca_table <- downloadHandler(filename="PCA.txt",
                                                content=function (file){
                                                   my_contributions_mod <- my_contributions()
                                                   
                                                   if (ncol(as.data.frame(my_contributions_mod))>=3) {
                                                      a<-t(as.data.frame(my_contributions_mod[30:1,]))
                                                      for (i in 1:nrow(a)) { 
                                                         a[i,]<-a[i,]/sum(a[i,])   # put the proportions
                                                      }
                                                      a<-a[,which(apply(a,2,sd)>0)]# remove signatures without variation
                                                      
                                                      samplesnames<-rownames(a)
                                                      rownames(a)<-1:(length(rownames(a)))
                                                      
                                                      write.table(x = data.frame(ID=rownames(a),Sample=samplesnames), file = file, sep="\t",quote=F,row.names=F)
                                                   }
                                                })
   
   
   })# observeEvent(input$run)

})
