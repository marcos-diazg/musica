# MuSiCa - Mutational Signatures in Cancer

<p>
<a href="https://github.com/marcos-diazg/musica/releases" alt="latest release version">
  <img src="https://img.shields.io/github/release/marcos-diazg/musica.svg" /></a>
<a href="https://github.com/marcos-diazg/musica/issues">
  <img src="https://img.shields.io/github/issues/marcos-diazg/musica"</a>
<a href="http://bioinfo.ciberehd.org:3838/MuSiCa/">
  <img src="https://img.shields.io/website?down_color=lightgrey&down_message=not%20available%20online%20version&up_color=green&up_message=online%20version&url=http%3A%2F%2Fbioinfo.ciberehd.org%3A3838%2FMuSiCa%2F" /></a>
</p>

MuSiCa (Mutational Signatures in Cancer) is a shiny-based web application aimed to visualize the somatic mutational profile of a series of provided samples (different formats are allowed) and to extract the contribution of the reported mutational signatures ([Alexandrov L.B. et al., Nature (2013)](http://dx.doi.org/10.1038/nature12477), [Catalogue Of Somatic Mutations In Cancer, COSMIC (2020)](http://cancer.sanger.ac.uk/cosmic/signatures)) on their variation profile. It is mainly based on the MutationalPatterns R package ([Blokzijl et al., Genome Medicine (2018)](https://doi.org/10.1186/s13073-018-0539-0)).

Please give credit and cite MuSiCa app when you use it for your genomic analysis ([Díaz-Gay et al., BMC Bioinformatics (2018)](https://doi.org/10.1186/s12859-018-2234-y)).

## Running MuSiCa

There are two ways to use MuSiCa, using the [online web version](http://bioinfo.ciberehd.org:3838/MuSiCa/) or the desktop local version. Please click on the online version link if you prefer the first option or follow the steps below if you would like to install and use the local version.

### Local version installation

First clone the repo on your desired location and access the repo folder.

```shell
cd path/to/desired_location
git clone https://github.com/marcos-diazg/musica.git
cd musica
```

Please execute the following commands to install and activate the required environment to run MuSiCa on your local machine. You will need to have conda installed. A quick way to do this is by installing miniconda3, available for different platforms at https://docs.conda.io/en/latest/miniconda.html.

```shell
source musica_setup.sh
```

### Local version runnning

Once you have your conda environment ready and activated you can now launch the app on your local machine by running the following commands in your terminal.

```shell
conda activate musica
echo 'shiny::runApp()' | Rscript -
```

