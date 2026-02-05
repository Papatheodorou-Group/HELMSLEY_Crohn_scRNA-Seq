# HELMSLEY Crohn scRNA-Seq

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
* R v4.0.3
* python v3.9
* R packages: ggplot2 v3.5.0.9000, dplyr 1.1.4, cowplot v1.1.3, circlize v0.4.13, SeuratObject v5.0.1, Seurat v4.0.5, Matrix v1.6-5, DoubletFinder v2.0.3, ComplexHeatmap v2.9.4, CellChat v2.1.2.
* python packages: scanpy v1.9.3, scvi-tools v1.0.2, anndata v0.9.2.

The code was executed on a GNU/Linux machine, using 4 cores and 64GB RAM for cell calling and 1 core and ~50GB RAM for downstream analysis.

## Installation

To run the analysis, download and save the present repository and the required packages as follows:

```
$ mkdir <MYDIR>
$ cd <MYDIR>
$ git clone https://github.com/fnadalin/HELMSLEY_Crohn_scRNA-Seq.git
$ git clone https://github.com/fnadalin/r_scripts.git
```

where ```<MYDIR>``` is the directory where analysis and scripts will be stored.

Make sure that both R and Rscript executables are in your PATH and install the required libraries as specified in [System requirements](#system-requirements).

Finally, create the conda environment for the query to reference mapping with scANVI:

```
$ conda create -n scvi-env python=3.9
$ conda activate scvi-env
$ conda install -c conda-forge scvi-tools
$ conda install -c conda-forge scanpy
```

## Instructions for reproducing the analysis

The analysis is organised in 7 steps. There are dependencies between them, so the analysis should be run in the same order as reported below.

### 1. Data download

To run the analysis starting from fastq files, these should be downloaded and stored as follows:

```
$ head -1 analysis/cell_calling/sample_info.tsv
SR1803_20220727_002_N_S_S	PD59664b	fastq/PD59664b/	2000	3151
$ mkdir -p fastq/PD59664b/
$ ls fastq/PD59664b/
PD59664b_S1_L001_R1_001.fastq.gz  PD59664b_S1_L001_R2_001.fastq.gz  PD59664b_S1_L001_singleton.fastq.gz
```

Alternatively, to run the analysis starting from count matrices, these should be downloaded and stored as follows:

```
$ cat HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/matrix_dir_004_N_E_S.tsv 
004_N_E_S	../cell_calling/SR1803_20221111_004_N_E_S/outs/filtered_feature_bc_matrix
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/cell_calling
$ mkdir -p SR1803_20221111_004_N_E_S/outs/filtered_feature_bc_matrix
$ ls SR1803_20221111_004_N_E_S/outs/filtered_feature_bc_matrix
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

The reference object from the Gut Cell Atlas should be downloaded from the GCA website and processed as follows: 

```
$ cd HELMSLEY_Crohn_scRNA-Seq/data/
$ mkdir GCA
$ cd GCA
$ wget https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad
$ cd ..
$ bash EXTRACT_GCA_DATA.sh
```

### 2. Cell calling

Generate gene-cell count matrices with CellRanger:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/cell_calling/
$ bash RUN_cellranger_slurm.sh
```

### 3. Clustering

Process the gene-cell count matrices with Seurat (feature selection, PCA, clustering, differential expression analysis):

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/
$ bash RUN_clustering_ccRegress_slurm.sh
$ bash RUN_DEA_ccRegress_slurm.sh
```

### 4. Low-quality filtering and doublet calling

Remove low-quality (dying) cells, predict doublets and filter them out:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/
$ bash RUN_filtering_slurm.sh
$ bash RUN_filtering_2ndRound.sh
```

### 5. Automatic annotation with scANVI

Generate a latent space with scANVI on the reference Gut Cell Atlas using the cell-type labelling as a covariate; add query cells with scANVI and predict cell-type labels:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/
$ bash RUN_scvi_GCA_slurm.sh
$ bash RUN_scanvi_HELMSLEY_slurm.sh
```

### 6. Second round of filtering and annotation curation

Re-cluster cells after filtering, generate automatic cluster annotation via predicted cell types and run differential expression analysois; use this information to generate a final cell-type annotation and possibly perform furteher quality filtering:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/
$ bash RUN_reclustering_slurm.sh
$ bash RUN_DEA_reclustering_slurm.sh
$ bash RUN_annotation_final_slurm.sh
```

### 7. Cell-cell interaction prediction

Predict cell-cell interactions from each condition separately, using CellChat, and compute the difference between them:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/analysis/gene_expression/
$ bash RUN_cell_cell_interactions_slurm.sh
```

### 7. Compositional modelling

Run scCODA for cell type compositional modelling analysis:

```
$ cd HELMSLEY_Crohn_scRNA-Seq/scripts/pertpy/
$ bash pertpy_slurm.sh HELMSLEY_Crohn_scRNA-Seq/scripts/pertpy/ path/to/files/
```
