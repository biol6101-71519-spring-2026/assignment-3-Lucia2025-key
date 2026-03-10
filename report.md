# RNA-Seq Workflow with Snakemake

## Overview

This project implements a reproducible RNA-Seq analysis workflow using Snakemake.  
The workflow processes paired-end RNA-seq data from raw FASTQ files to downstream gene expression analysis.

The pipeline performs quality control, read trimming, alignment, gene quantification, and downstream statistical analysis.

---

Quality control of the raw sequencing reads indicated generally good read quality across all samples. FastQC reports showed that most reads had high Phred quality scores across the read length. Adapter contamination and low-quality bases were minimal; trimming had only a modest effect on the total number of reads retained. The trimming step removed a small fraction of low-quality bases and adapter sequences, resulting in clean reads suitable for downstream alignment.

Alignment of the trimmed reads to the reference genome using STAR produced consistent alignment rates across all samples. The majority of reads successfully aligned to the reference genome, indicating good sequencing quality and appropriate reference annotation. Alignment statistics were similar across control and treatment samples, suggesting no major technical issues or sample-specific biases during sequencing or preprocessing.

Principal Component Analysis (PCA) was performed using normalized gene expression counts to explore relationships between samples. The PCA plot showed that samples clustered closely together with minimal separation between control and treatment groups. This pattern suggests limited global differences in gene expression between the two conditions in this dataset. The high similarity among samples is consistent with the high correlation observed across samples and may reflect the use of duplicated or highly similar input reads used for testing the workflow.

However, because the dataset used for testing contains highly similar samples, the downstream exploratory analyses (PCA and correlation) show minimal variation between conditions
---
