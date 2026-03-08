# RNA-Seq Workflow with Snakemake

## Overview

This project implements a reproducible RNA-Seq analysis workflow using Snakemake.  
The workflow processes paired-end RNA-seq data from raw FASTQ files to downstream gene expression analysis.

The pipeline performs quality control, read trimming, alignment, gene quantification, and downstream statistical analysis.

---

## Workflow Steps

### 1. Quality Control
Raw sequencing reads are assessed using **FastQC** to evaluate read quality, adapter contamination, and other sequencing metrics.

### 2. Read Trimming
Low-quality bases and adapter sequences are removed from the raw reads.

### 3. Alignment
Trimmed reads are aligned to the reference genome using **STAR**.

### 4. BAM Processing
Aligned reads are sorted and indexed to produce analysis-ready BAM files.

### 5. Gene Quantification
Reads are assigned to genes using **featureCounts**, producing count tables for each sample.

### 6. Count Aggregation
Individual count files are merged into a combined gene expression count matrix.

### 7. Exploratory Data Analysis
Two analyses are performed:

- **Principal Component Analysis (PCA)** to visualize variation between samples
- **Correlation Heatmap** to evaluate similarity among samples

### 8. DESeq2 Dataset Preparation
A DESeq2 dataset is created and normalized counts are generated for downstream 
TThe workflow successfully processed RNA-seq data from raw reads through alignment, gene quantification, and exploratory analysis. The results demonstrate that the pipeline is functioning correctly and producing the expected outputs.

However, because the dataset used for testing contains highly similar samples, the downstream exploratory analyses (PCA and correlation) show minimal variation between conditions
---