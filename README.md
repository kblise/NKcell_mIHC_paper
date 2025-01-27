# Natural killer cells occupy unique spatial neighborhoods in human HER2- and HER2+ breast cancers

This repository contains all code necessary to reproduce the spatial analysis results and figures of the publication "Natural killer cells occupy unique spatial neighborhoods in human HER2- and HER2+ breast cancers," which can be found here: [https://doi.org/10.1186/s13058-025-01964-4](https://doi.org/10.1186/s13058-025-01964-4). All data, including the output of the multiplex immunohistochemistry computational imaging processing workflow for each tissue region and metadata for each region, are available on Zenodo: [DOI: 10.5281/zenodo.10632694](https://doi.org/10.5281/zenodo.10632694).

# Steps to Create Figures

## 1) Clone repository

**a.** Open Terminal and navigate to desired directory: `cd Desired/Directory/Here/`

**b.** Clone repo: `git clone https://github.com/kblise/NKcell_mIHC_paper.git`

## 2) Create new conda environment

**a.** Install Conda if not already installed: [Instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

**b.** Navigate into NKcell_mIHC_paper directory: `cd NKcell_mIHC_paper/`

**b.** Create new conda environment with necessary packages: `conda env create -f nkEnv.yml`

**c.** Activate nkEnv conda environment: `source activate nkEnv` or `conda activate nkEnv` depending on Conda version

## 3) Run bash script to create directories, download data, and run Python script to generate figures

**a.** Run bash script from command line: `bash nkCell_mIHC_paper.sh`

**Note: Csv files created to generate figures will be saved to the 'results/dfCreated' folder and figures will be saved to the 'results/figures' folder.**

This program is intended for Python version 3.
