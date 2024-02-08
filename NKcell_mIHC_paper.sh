#!/usr/bin/env bash

#Author: Katie Blise
#Date: February 2024
#This bash script will: create directories, download data files, and call the nkMakeFigures.py python scprit to generate the results and figures in the manuscript "Natural Killer cells occupy unique spatial neighborhoods in human HER2- and HER2+ breast cancers"
#Call "bash NKcell_mIHC_paper.sh" from the command line to run this script.

#create directory structure
mkdir results
cd results
mkdir {dfCreated,figures} #two subfolders in results
cd dfCreated
mkdir updatedCsvs
cd ../.. #go back to home directory to download data

#download data from Zenodo: DOI: 10.5281/zenodo.10632694
#all data lives in data.zip file, which contains 2 folders: mIHC_files and metadata
wget https://zenodo.org/records/10632694/files/data.zip
#unzip data.zip file to create data folder
unzip data.zip

#run python script to generate results and figures
python nkMakeFigures.py
