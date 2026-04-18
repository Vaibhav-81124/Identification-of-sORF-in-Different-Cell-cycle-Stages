# sORF Discovery Framework

**Author:** Vaibhav D  
**Project Type:** Undergraduate Thesis (B.E. Biotechnology)

---

## What this project does

This project is a lightweight, reproducible framework for discovering and validating novel small open reading frames (sORFs) using RNA-seq and Ribo-seq data.

It standardizes a multi-step workflow—from data preprocessing to translation validation and filtering—into a structured system that can be reused across datasets instead of running one-off analyses.

---

## Why I built this

sORF discovery workflows are often fragmented and difficult to reproduce. Small variations in preprocessing, filtering, or parameter selection can lead to inconsistent results across replicates and datasets.

This framework was built to:
- Standardize the end-to-end workflow  
- Reduce variability across datasets  
- Enable reproducible and consistent outputs  
- Provide a reusable system instead of ad-hoc scripts  

---

## Dataset

- **Source:** PRJNA316618  
- **Data Types:** RNA-seq, Ribo-seq  
- **Conditions:** Multiple cell-cycle stages (Asynchronous, M, S phase) with replicates  
- **Reference Genome:** GRCh38  
- **Annotation:** Ensembl GTF  

---

## Workflow Overview

### 1. sORF Discovery
- Parse transcript sequences and exon structures  
- Map transcript coordinates to genome  
- Identify ORFs:  
  - Start codon: ATG  
  - Length: 10–100 aa (refined ≥15 aa)  
- Remove CDS-overlapping ORFs  

### 2. Filtering & Cleaning
- Remove duplicate peptides and loci  
- Remove nested ORFs (retain longest)  
- Apply minimum length cutoff (≥15 aa)  

### 3. RNA-seq Processing
**Tools:** FastQC, Cutadapt, STAR  
- Quality control and trimming  
- Alignment to genome  
- Generate sorted BAM files  

### 4. Ribo-seq Processing
**Tools:** Cutadapt, SortMeRNA, STAR, Samtools  
- Adapter trimming and rRNA removal  
- Alignment and filtering (RPFs: 26–34 nt)  
- P-site assignment  

### 5. Translation Evidence
- Convert sORFs to BED format  
- Intersect with P-site data (BEDTools)  
- Aggregate ribosome counts  
- Apply threshold (≥10 reads)  

### 6. Periodicity Analysis
- Evaluate 3-nt periodicity  
- Compute frame distribution  
- Apply thresholds:  
  - Minimum reads ≥8  
  - Frame0 fraction ≥0.55  

### 7. Expression Integration
- Quantify RNA expression (featureCounts)  
- Integrate RNA + Ribo signals  

### 8. Annotation & Functional Analysis
- Identify host genes  
- Classify ORFs (intragenic, intergenic, etc.)  
- Perform enrichment analysis (g:Profiler)  
- Analyze protein interactions (STRING)  

---

## Key Features

- End-to-end workflow from discovery to validation  
- Integration of RNA-seq and Ribo-seq data  
- Multi-step filtering for high-confidence ORFs  
- Translation validation using ribosome profiling  
- Designed for reproducibility and reuse  

---

## Outputs

- Novel sORF candidates  
- High-confidence translated ORFs  
- Microprotein candidates (≤100 aa)  
- Host gene annotations  
- Functional enrichment results  

---

## How to use

1. Prepare RNA-seq and Ribo-seq input data  
2. Run scripts in sequence (or via wrapper if available)  
3. Collect outputs from each stage  

> Note: Full wrapper implementation is currently private due to ongoing preprint. This repository provides core scripts and workflow structure.

---

## Summary

This framework enables systematic discovery and validation of novel translated sORFs by combining transcriptomics and ribosome profiling data. It focuses on reproducibility, consistency, and usability across datasets.
