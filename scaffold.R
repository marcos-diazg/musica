

mat <- read.table("aux_files/demo_data/proba.tsv", sep="\t", header=T)

mat <- fread("aux_files/demo_data/proba.tsv",header=T,sep="\t",data.table=F)

         validate(
            need(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4  , error_message)
         )




mat_list<-list()
         for (w in 1:2){
		mat <- fread("aux_files/demo_data/proba.tsv",header=T,sep="\t",data.table=F)
            condition <- length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4
            
            mat_list[[w]] <- condition
         }
         mat_vector <- as.vector(do.call(cbind, mat_list))

length(grep("FALSE", mat_vector))==0


 mat<-read.xlsx("aux_files/demo_data/proba.tsv",1)

mat<-data.frame(read_xls("aux_files/demo_data/PROBA_NO_HEADER.xls", col_names = TRUE,col_types = "text"))

mat_list<-list()
         for (w in 1:2){
            

               mat<-data.frame(read_xls("aux_files/demo_data/PROBA_NO_HEADER.xls",col_names = TRUE,col_types = "text"))

            
            condition <- as.character(length(grep("TRUE", (c("CHROM", "POS", "REF", "ALT") %in% colnames(mat))))==4)
            
            mat_list[[w]] <- condition
         }
         mat_vector <- as.vector(do.call(cbind, mat_list))
