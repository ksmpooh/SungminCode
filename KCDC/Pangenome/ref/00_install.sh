#!/bin/bash
#wget https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh
#bash Mambaforge-24.3.0-0-Linux-x86_64.sh
#/Users/gyeongheon/Desktop/WORK/kogo_LongRead/00_scripts/00_mamba/mambaforge/bin/mamba init
source ~/.zshrc
#mamba create -n tutorial
#mamba activate tutorial
mamba create -n genomescope2 -c conda-forge -c bioconda genomescope2 kmc
mamba create -n basicGenomics -c conda-forge -c bioconda bioawk samtools seqtk hifiasm busco svim-asm assembly-stats
mamba create -n repeat -c conda-forge -c bioconda repeatmasker
mamba create -n Rplot -c conda-forge -c bioconda -c r r-plyr r-dplyr r-ggplot2