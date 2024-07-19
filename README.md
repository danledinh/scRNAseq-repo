# scRNAseq-repo

This repository contains the code for the main analyses performed in the ["Comparative Analysis of Single Cell RNA Sequencing Technologies" paper](https://doi.org/10.1101/2024.06.18.599579) (De Simone and Hoover, et. al).

Main editors for this repo are Jonathan Hoover and Daniel Le. For correspondence regarding code please contact Daniel Le (led13@gene.com). For correspondence regarding experimental design please contact Spyros Darmanis (darmanis@gene.com).

scRNA kits processed in this project include:

- Chromium Single Cell 3' Reagent Kit, v3.1 Chemistry. (10x Genomics): '10X_3'
- Chromium Single Cell 5' Reagent Kit, v2 Chemistry. (10x Genomics): '10X_5'
- Chromium Fixed RNA Profiling Reagent Kit (10x Genomics): '10X_FRP'
- PIPseq T20 3' Single Cell RNA Kit v4.0 (Fluent Biosciences): 'Fluent'
- BD Rhapsody WTA Reagent Kit (Becton Dickinson): 'BD'
- HIVE CLX scRNA Seq Kit v1 (Honeycomb Biotechnologies): 'Honeycomb'
- Evercode WT v2 (Parse Biosciences): 'Parse'
- Single Cell RNA Kit (Scale Biosciences): 'Scale'
- ASTERIA Single-cell RNASeq Kit (Scipio Bioscience): 'Scipio'

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

After setting up the conda environment, some notebooks can be be run directly using data deposited on the CZ Biohub's CellxGene repo (see the **Resources** section below). Some just act as reference for how data was generated. Each notebook will mention whether it can be rerun directly using the resources made available within this repo. 

**For those that are able to be run using the AnnData objects uploaded to CellxGene, please run the `notebooks/CellxGene_Download.ipynb` notebook first. It will generate a `results/anndata_objects/` directory containing all the annotated .h5ad files for each kit along with the harmony integrated .h5ad file.**

## Resources

The following are links to relevant resources related to this project:

- Pre-print: https://www.biorxiv.org/content/10.1101/2024.06.18.599579v1

- SRA BioProject: PRJNA1106903 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1106903)

- CZ Biohub CellxGene Collection: 398e34a9-8736-4b27-a9a7-31a47a67f446 (https://cellxgene.cziscience.com/collections/398e34a9-8736-4b27-a9a7-31a47a67f446)
