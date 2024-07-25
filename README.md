# HELMSLEY_Crohn_scRNA-Seq

This repository contains the code to reproduce the analysis of normal and Crohn's disease scRNA-Seq samples.

[Description](#description)  
[System requirements](#system-requirements)  
[Installation](#installation)  
[Instructions for reproducing the analysis](#instructions-for-reproducing-the-analysis)  
[Contact](#contact) 

## Description

Data, scripts and pipelines are organised in three folders:

* analysis/ is divided into sub-directories, each referring to a specific analysis step. They contain the commands used to run the pipelines and information on the path containing the raw and processed data.
* scripts/ contains the pipeline scripts. They are also grouped by analysis type. 
* data/ contains external data and information required for the analysis.

General-purpose custom code is stored in separate repositories (see below).

## System requirements

To reproduce the analysis, the following software packages are required:

* CellRanger v7.1.0
* scvi vXXX
* R v4.0.3
* python vXXX
* R packages: 
* python packages:

Additionally, the following packages should be downloaded and installed (see [Installation](#installation)):

* https://github.com/fnadalin/r_scripts

The code was executed on a GNU/Linux machine, using 4 cores and 64GB RAM for cell calling and 1 core and ~40GB RAM for downstream analysis.

## Installation

To run the analysis, download and save the present repository and the required packages as follows:

```
$ mkdir <MYDIR>
$ cd <MYDIR>
$ git clone https://github.com/fnadalin/gbc_sum159pt_paper.git
$ git clone https://github.com/fnadalin/r_scripts.git
```

where ```<MYDIR>``` is the directory where analysis and scripts will be stored.

Make sure that both R and Rscript executables are in your PATH and install the required libraries as specified in [System requirements](#system-requirements).

Finally, create the conda environment:

```
$ ...
```

## Instructions for reproducing the analysis

The analysis is organised in XXX steps. There are dependencies between them, so the analysis should be run in the same order as reported below.

### 1. Data download

To run the analysis from the beginning, fastq files should be downloaded and stored at (or linked to) the path as reported in the file "analysis/cell\_calling/sample\_info.tsv".

If starting the analysis from processed data (i.e., count matrices), they should be downloaded and stored as follows:

```
$ wget ...
```

The reference object from the Gut Cell Atlas should be downloaded from the GCA website and processed as follows: 

```
$ cd data/
$ mkdir GCA
$ cd GCA
$ wget https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad
$ cd ..
$ bash EXTRACT_GCA_DATA.sh
```

### 2. Cell calling

### 3. Clustering

### 4. Low-quality filtering and doublet calling

### 5. Automatic annotation with scANVI

### 6. Second round of filtering and annotation curation



