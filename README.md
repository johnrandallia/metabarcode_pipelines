# metabarcode_pipelines


The following tools and worksteps describe multiple ways to preprocess metabarcoding datasets.
Here, different (primary) denoising tools are used and each tool is executed with the default denoising stringency, but also 9 different denoising stringencies.

## Installation of working environment and Utilities for quality control and preprocessing of metabarcoding data
To avoid dependency issues and for easier management of the tools we used miniconda with a new environment.

The conda environment was created as follows:

```
conda create --name metabarcode_data_analysis_env python=3.10
conda init fish # other shells, such as bash or zsh, can also be used
conda activate metabarcode_data_analysis_env
```

Now, within the newly created and activated metabarcode_data_analysis_env, the installation of these tools can be done using the **installation_utils.sh** bash script. 
- fastqc 
- fastx_trimmer
- fastq-pair
- trimmomatic
- Samtools pairtosam 
- vsearch





