# Mutational Signatures Shiny App

Mutational Signatures ShinyApp is an app aimed to ...


# Running the app

There are many ways to download and run Mutational Signatures ShinyApp:

First check the dependencies (you only need to do this once):

```R
if(!require(MutationalPatterns)) source("http://www.bioconductor.org/biocLite.R");biocLite("MutationalPatterns")
if(!require(BSgenome.Hsapiens.NCBI.GRCh38)) source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.UCSC.hg19")
if(!require(BSgenome.Hsapiens.1000genomes.hs37d)) source("http://www.bioconductor.org/biocLite.R");biocLite("BSgenome.Hsapiens.1000genomes.hs37d")

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(shiny)) install.packages("shiny")
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

Authors: Marcos Díaz-Gay & Maria Vila-Casadesús


