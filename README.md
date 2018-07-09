# MuSiCa - Mutational Signatures in Cancer

MuSiCa (Mutational Signatures in Cancer) is a shiny-based web application aimed to visualize the somatic mutational profile of a series of provided samples (different formats are allowed) and to extract the contribution of the reported mutational signatures ([Alexandrov L.B. et al., Nature (2013)](http://dx.doi.org/10.1038/nature12477), [Catalogue Of Somatic Mutations In Cancer, COSMIC (2017)](http://cancer.sanger.ac.uk/cosmic/signatures)) on their variation profile. It is mainly based on the MutationalPatterns R package ([Blokzijl et al., Genome Medicine (2018)](https://doi.org/10.1186/s13073-018-0539-0)).

Please give credit and cite MuSiCa app when you use it for your genomic analysis ([Díaz-Gay et al, BMC Bioinformatics (2018)](https://doi.org/10.1186/s12859-018-2234-y)).

## Running the app

There are many ways to download and run Mutational Signatures ShinyApp:

First check the dependencies (you only need to do this once):




***It is mandatory to use Shiny package version 1.0.5, because the application is not working with the latest version of Shiny (1.1.0). We will work on this soon. Sorry for the inconveniences.***




```R
if(!require(MutationalPatterns)) {source("http://www.bioconductor.org/biocLite.R");biocLite("MutationalPatterns")}
if(!require(VariantAnnotation)) {source("http://www.bioconductor.org/biocLite.R");biocLite("VariantAnnotation")}
if(!require(BSgenome.Hsapiens.UCSC.hg38)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.UCSC.hg38")}
if(!require(BSgenome.Hsapiens.UCSC.hg19)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.UCSC.hg19")}
if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.1000genomes.hs37d5")}

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
# ~/music, first go to that directory, then use runApp().
setwd("~/musica")
runApp()
```
