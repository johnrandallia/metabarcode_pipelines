# metabarcode_pipelines
This repository contains scripts for multiple (pre)processing approaches of metabarcoding datasets. 
A series of denoising tools are employed, each evaluated across ten stringency levels.

## Installation
To avoid dependency issues and for easier management of the tools we used miniconda with a new environment.
The conda environment was created as follows:

```
> conda create --name metabarcode_data_analysis_env python=3.10
> conda init fish # other shells, such as bash or zsh, can also be used
> conda activate metabarcode_data_analysis_env
```

Now, within the newly created and activated metabarcode_data_analysis_env, the installation of these tools can be done using the **installation_utils.sh** bash script. 
- fastqc 
- fastx_trimmer
- fastq-pair
- trimmomatic
- Samtools pairtosam 
- vsearch


Taxonomic classification (with Kraken2) was performed within a dedicated conda environment, created as follows:
```
> conda create -n kraken2_env -c conda-forge -c bioconda -c defaults --strict-channel-priority kraken2
```
To extract the desired taxa KrakenTools was used, which can be downloaded on github (https://github.com/jenniferlu717/KrakenTools.git)



The binary file for the usearch preprocessing toolkit, as well as the denoising tool UNOISE3 is not included here due to licensing restrictions.
Please download it from the usearch webpage.


LULU can be installed, as described by its author, using the R package devtools (https://github.com/tobiasgf/lulu)
```
> library(devtools)
> install_github("tobiasgf/lulu")  
```


metaMATE was installed in a dedicated conda environment, as recommended in the github repository of metaMATE (https://github.com/tjcreedy/metamate?tab=readme-ov-file#installation)
```
> conda create -n metamate_env metamate -c bioconda -c conda-forge 
```

## Execution of the pipeline
The provided scripts in this github repo are semi-automatic. 
All of them need adjustment of the working directory/ paths.
Once that is done, execute:
```
> ./run_preprocessing_EPA45.sh
> ./run_LULU.sh
> ./run_metamate.sh
```







