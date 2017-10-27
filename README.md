# MuSiC - Mutational Signatures in Cancer

MuSiC (Mutational Signatures in Cancer) is a shiny-based web application aimed to visualize the somatic mutational profile of a series of provided samples (different formats are allowed) and to extract the contribution of the reported mutational signatures ([Alexandrov L.B. et al., Nature (2013)](http://dx.doi.org/10.1038/nature12477), [Catalogue Of Somatic Mutations In Cancer, COSMIC (2017)](http://cancer.sanger.ac.uk/cosmic/signatures)) on their variation profile. It is mainly based on the [MutationalPatterns](https://doi.org/10.1101/071761) R package.


## Running the app

There are many ways to download and run Mutational Signatures ShinyApp:

First check the dependencies (you only need to do this once):

```R
if(!require(MutationalPatterns)) {source("http://www.bioconductor.org/biocLite.R");biocLite("MutationalPatterns")}
if(!require(VariantAnnotation)) {source("http://www.bioconductor.org/biocLite.R");biocLite("VariantAnnotation")}
if(!require(BSgenome.Hsapiens.UCSC.hg38)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.UCSC.hg38")}
if(!require(BSgenome.Hsapiens.UCSC.hg19)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.UCSC.hg19")}
if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)) {source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.1000genomes.hs37d5")}

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(plotly)) install.packages("plotly")
if(!require(heatmaply)) install.packages("heatmaply")
if(!require(gplots)) install.packages("gplots")
if(!require(reshape2)) install.packages("reshape2")
if(!require(data.table)) install.packages("data.table")
if(!require(gdata)) install.packages("gdata")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(shiny)) install.packages("shiny")
if(!require(shinyBS)) install.packages("shinyBS")
if(!require(shinysky)) {library(devtools); devtools::install_github("AnalytixWare/ShinySky")}
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(shinySky)) install.packages("shinySky")
if(!require(V8)) install.packages("V8")
```

Then load the shiny app:

```R
library(shiny)

# Easiest way is to use runGitHub
runGitHub("MutSignatures_ShinyApp/MutSignatures", "marcos-diazg")
```

Other ways to load the app:

```R
# Run a tar or zip file directly
runUrl("https://github.com/marcos-diazg/MutSignatures_ShinyApp/MutSignatures/archive/master.tar.gz")
runUrl("https://github.com/marcos-diazg/MutSignatures_ShinyApp/MutSignatures/archive/master.zip")

# Using runApp(),  first clone the repository with git. If you have cloned it into
# ~/MutSignatures_ShinyApp/MutSignatures, first go to that directory, then use runApp().
setwd("~/MutSignatures_ShinyApp/MutSignatures")
runApp()
```

Authors: Marcos Díaz-Gay, Sebastià Franch-Expósito & Maria Vila-Casadesús


