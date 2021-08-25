# Script to install environment for MuSiCa
# Author: Marcos Díaz-Gay
# Date: Aug 22, 2021
# Conda environment: musica

conda env create -f musica_env.yml
conda activate musica

echo 'devtools::install_github("AnalytixWare/ShinySky")' | Rscript -
echo 'webshot::install_phantomjs(force=TRUE)' | Rscript -
echo 'BiocManager::install("MutationalPatterns")' | Rscript -
echo 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")' | Rscript -
echo 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")' | Rscript -
echo 'BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")' | Rscript -
conda install -c conda-forge r-v8 -y
