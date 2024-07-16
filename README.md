# scRNAseq-repo

This repository contains the code for the main analyses performed in the ["Comparative Analysis of Single Cell RNA Sequencing Technologies" paper](https://doi.org/10.1101/2024.06.18.599579) (De Simone and Hoover, et. al).

Main editors for this repo are Jonathan Hoover and Daniel Le. For correspondence regarding code please contact Daniel Le (led13@gene.com). For correspondence regarding experimental design please contact Spyros Darmanis (darmanis@gene.com).

## Reproducibility

Analyses were performed using the mamba environment shown in `./mamba_environments/scanpy-default-mamba.yml`. To install this virtual environment using miniforge, make sure miniforge is installed (see https://github.com/conda-forge/miniforge) and then execute the following:

First, clone this directory to your desired local directory (make sure git is installed)

```bash
cd </path/to/clone/location>
git clone https://github.com/danledinh/scRNAseq-repo.git
```

Then, cd into the cloned directory and create the conda repo from the .yml file provided with your desired <envname>

```bash
cd </path/to/scRNAseq-repo>
conda env create -f ./mamba_environments/scanpy-default-mamba.yml -n <envname>
```

## Resources

Pre-print: https://www.biorxiv.org/content/10.1101/2024.06.18.599579v1

SRA BioProject: PRJNA1106903 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1106903)
