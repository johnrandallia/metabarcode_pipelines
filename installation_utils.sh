#!/bin/bash

# Installation of Utilities

conda install -c bioconda fastqc
conda install -c bioconda fastx_toolkit
conda install -c bioconda fastq-pair
conda install -c bioconda trimmomatic
conda install -c "bioconda/label/cf201901" samtools
conda install -c bioconda vsearch
conda install -c conda-forge pigz
conda install biopython

# install cutadapt
sudo apt install cutadapt
pip3 install levenshtein

# install seqkit
conda install bioconda::seqkit

# Install DnoisE
conda install -c adriantich dnoise

# Install R
conda install -c conda-forge r-base=4.0.0

# R packages to execute DADA2
Rscript -e "install.packages(c('dada2','dplyr','stringr','seqinr', 'compiler'), repos = 'https://cloud.r-project.org')"
