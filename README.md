# MuSiCa - Mutational Signatures in Cancer

<p>
<a href="https://github.com/marcos-diazg/musica/releases" alt="latest release version">
        <img src="https://img.shields.io/github/release/marcos-diazg/musica.svg" /></a>
</p>

MuSiCa (Mutational Signatures in Cancer) is a shiny-based web application aimed to visualize the somatic mutational profile of a series of provided samples (different formats are allowed) and to extract the contribution of the reported mutational signatures ([Alexandrov L.B. et al., Nature (2013)](http://dx.doi.org/10.1038/nature12477), [Catalogue Of Somatic Mutations In Cancer, COSMIC (2017)](http://cancer.sanger.ac.uk/cosmic/signatures)) on their variation profile. It is mainly based on the MutationalPatterns R package ([Blokzijl et al., Genome Medicine (2018)](https://doi.org/10.1186/s13073-018-0539-0)).

Please give credit and cite MuSiCa app when you use it for your genomic analysis ([Díaz-Gay et al., BMC Bioinformatics (2018)](https://doi.org/10.1186/s12859-018-2234-y)).

## Running the app

There are many ways to download and run Mutational Signatures ShinyApp:

First check the dependencies (you only need to do this once):






```R
if(!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if(!require(MutationalPatterns)) {BiocManager::install("MutationalPatterns")}
if(!require(VariantAnnotation)) {BiocManager::install("VariantAnnotation")}
if(!require(BSgenome.Hsapiens.UCSC.hg38)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")}
if(!require(BSgenome.Hsapiens.UCSC.hg19)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)) {BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")}

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(heatmaply)) install.packages("heatmaply")
if(!require(gplots)) install.packages("gplots")
if(!require(reshape2)) install.packages("reshape2")
if(!require(data.table)) install.packages("data.table")
if(!require(readxl)) install.packages("readxl")
if(!require(openxlsx)) install.packages("openxlsx")

if(!require(shiny)) install.packages("shiny")
if(!require(shinyBS)) install.packages("shinyBS")
if(!require(devtools)) install.packages("devtools")
if(!require(shinysky)) {library(devtools); devtools::install_github("AnalytixWare/ShinySky")}
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(V8)) install.packages("V8")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(plotly)) install.packages("plotly")
if(!require(webshot)) install.packages("webshot")
webshot::install_phantomjs()
```

Then load the shiny app:

```R
library(shiny)

# Easiest way is to use runGitHub
runGitHub("musica", "marcos-diazg")
```

Other ways to load the app:

```R
# Run a tar or zip file directly
runUrl("https://github.com/marcos-diazg/musica/archive/master.tar.gz")
runUrl("https://github.com/marcos-diazg/musica/archive/master.zip")

# Using runApp(),  first clone the repository with git. If you have cloned it into
# ~/musica, first go to that directory, then use runApp().
setwd("~/musica")
runApp()
```


## Session information

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.5 LTS

Matrix products: default
BLAS: /opt/R-3.4.2/lib64/R/lib/libRblas.so
LAPACK: /opt/R-3.4.2/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.46.0                  
 [3] rtracklayer_1.38.3                readxl_1.0.0                     
 [5] openxlsx_4.0.17                   gplots_3.0.1                     
 [7] heatmaply_0.14.1                  viridis_0.5.1                    
 [9] viridisLite_0.3.0                 VariantAnnotation_1.24.5         
[11] Rsamtools_1.30.0                  Biostrings_2.46.0                
[13] XVector_0.18.0                    SummarizedExperiment_1.8.1       
[15] DelayedArray_0.4.1                matrixStats_0.53.0               
[17] data.table_1.11.4                 reshape2_1.4.3                   
[19] MutationalPatterns_1.4.2          NMF_0.20.6                       
[21] Biobase_2.38.0                    cluster_2.0.6                    
[23] rngtools_1.2.4                    pkgmaker_0.22                    
[25] registry_0.5                      GenomicRanges_1.30.3             
[27] GenomeInfoDb_1.14.0               IRanges_2.12.0                   
[29] S4Vectors_0.16.0                  BiocGenerics_0.24.0              
[31] plotly_4.7.1                      ggplot2_3.0.0                    
[33] shinythemes_1.1.1                 V8_1.5                           
[35] shinyjs_1.0                       shinysky_0.1.2                   
[37] plyr_1.8.4                        RJSONIO_1.3-0                    
[39] shinyBS_0.61                      shiny_1.0.5                      

loaded via a namespace (and not attached):
 [1] colorspace_1.3-2         class_7.3-14             modeltools_0.2-21       
 [4] mclust_5.4               ggdendro_0.1-20          flexmix_2.3-14          
 [7] bit64_0.9-7              AnnotationDbi_1.40.0     mvtnorm_1.0-7           
[10] codetools_0.2-15         doParallel_1.0.11        robustbase_0.92-8       
[13] jsonlite_1.5             gridBase_0.4-7           kernlab_0.9-25          
[16] compiler_3.4.2           httr_1.3.1               assertthat_0.2.0        
[19] Matrix_1.2-11            lazyeval_0.2.1           htmltools_0.3.6         
[22] prettyunits_1.0.2        tools_3.4.2              bindrcpp_0.2            
[25] gtable_0.2.0             glue_1.2.0               GenomeInfoDbData_1.0.0  
[28] dplyr_0.7.4              Rcpp_0.12.17             cellranger_1.1.0        
[31] trimcluster_0.1-2        gdata_2.18.0             iterators_1.0.9         
[34] fpc_2.1-11               stringr_1.2.0            mime_0.5                
[37] gtools_3.8.1             XML_3.98-1.9             dendextend_1.6.0        
[40] DEoptimR_1.0-8           zlibbioc_1.24.0          MASS_7.3-47             
[43] scales_0.5.0             TSP_1.1-5                RColorBrewer_1.1-2      
[46] yaml_2.1.16              curl_3.1                 memoise_1.1.0           
[49] gridExtra_2.3            biomaRt_2.34.2           stringi_1.1.6           
[52] RSQLite_2.1.1            gclus_1.3.1              foreach_1.4.4           
[55] RMySQL_0.10.15           seriation_1.2-2          caTools_1.17.1          
[58] GenomicFeatures_1.30.3   BiocParallel_1.12.0      rlang_0.2.1             
[61] pkgconfig_2.0.1          prabclus_2.2-6           bitops_1.0-6            
[64] pracma_2.1.4             lattice_0.20-35          purrr_0.2.4             
[67] bindr_0.1                GenomicAlignments_1.14.1 htmlwidgets_1.0         
[70] cowplot_0.9.2            bit_1.1-12               magrittr_1.5            
[73] R6_2.2.2                 DBI_1.0.0                pillar_1.1.0            
[76] whisker_0.3-2            withr_2.1.2              RCurl_1.95-4.11         
[79] nnet_7.3-12              tibble_1.4.2             KernSmooth_2.23-15      
[82] progress_1.1.2           grid_3.4.2               blob_1.1.1              
[85] digest_0.6.15            diptest_0.75-7           webshot_0.5.0           
[88] xtable_1.8-2             tidyr_0.8.0              httpuv_1.3.5            
[91] munsell_0.4.3
```
