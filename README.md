# ad_oral_microbiome
Pipeline for metagenomic analysis of salivary microbiome alterations in Alzheimer's disease using diversity metrics and machine learning.

# Salivary Microbiome Dysbiosis as a Potential Biomarker for Alzheimer's Disease: A Metagenomics Analysis

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

This repository contains the data, analytical pipeline, and results for the study:

**"Salivary Microbiome Dysbiosis as a Potential Biomarker for Alzheimer's Disease: A Metagenomics Analysis"**  
*Edward Jenner Tettevi, Samuel Armoo, Mike Y. Osei-Atweneboana*

This study investigates the differences in salivary microbiome composition and function between Alzheimer's disease (AD) patients and age-matched healthy controls, using comprehensive metagenomic and machine learning approaches. The findings demonstrate that salivary microbiome dysbiosis is strongly associated with AD and has robust diagnostic potential.

---

## Table of Contents

- [Background](#background)
- [Objectives](#objectives)
- [Methods](#methods)
  - [Data Acquisition](#data-acquisition)
  - [Taxonomic and Functional Analysis](#taxonomic-and-functional-analysis)
  - [Machine Learning](#machine-learning)
- [Key Results](#key-results)
- [Installation & Usage](#installation--usage)
- [Data Availability](#data-availability)
- [Main Findings](#main-findings)
- [License](#license)

---

## Background

Alzheimer’s disease (AD) is the most common cause of dementia, affecting over 55 million people worldwide. Current diagnostic methods are invasive, expensive, and often inaccessible. This study explores the potential of the salivary microbiome as a non-invasive biomarker for early detection and monitoring of AD.

---

## Objectives

1. Characterize the taxonomic composition of the salivary microbiome in AD patients versus controls.
2. Identify differentially abundant microbial taxa associated with AD.
3. Analyze functional differences in microbial communities.
4. Develop machine learning models to predict AD status based on microbial features.

---

## Methods

### Data Acquisition

- **Source:** NCBI Sequence Read Archive (SRA) [BioProject PRJNA770746](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770746)
- **Samples:** Saliva from AD patients and age-matched healthy controls.
- **Preprocessing:** Data processed and normalized via Chan Zuckerberg ID (CZID) and custom Python scripts.

### Taxonomic and Functional Analysis

- **Diversity Metrics:** Shannon and Simpson indices, Bray-Curtis, Canberra, Chebyshev, Hamming, Matching, and Euclidean distances.
- **Statistical Analysis:** Mann-Whitney U test, PERMANOVA, Benjamini-Hochberg correction.
- **Functional Pathways:** Analyzed using Pathway-Tools Python API for pathway reconstruction and compound characterization.
- **Network Analysis:** Microbial co-occurrence and metabolic pathway correlation networks constructed using Spearman rank correlation.

### Machine Learning

- **Algorithms:** Logistic Regression, Support Vector Machine (SVM), Random Forest.
- **Validation:** Stratified 5-fold cross-validation.
- **Metrics:** Accuracy, ROC-AUC, precision, F1-score, Matthews Correlation Coefficient.

---

## Key Results

- **Diversity:** AD patients exhibited significantly higher oral microbial diversity (Shannon p=0.022; Simpson p=0.005).
- **Taxonomic Shifts:** Elevated *Taylorella asinigenitalis*, reduced *Streptococcus suis* in AD.
- **Functional Changes:** Significant downregulation of glycolysis and short-chain fatty acid production pathways in AD.
- **Network Topology:** AD microbiome networks are more fragmented, with reduced integration.
- **Machine Learning:** Classifiers achieved high diagnostic accuracy (AUC = 0.93, accuracy = 86%).

---

## Installation & Usage

1. **Clone this repository:**
   ```bash
   git clone https://github.com/ejtettevi/ad_oral_microbiome.git
   cd ad_oral_microbiome
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Data Preparation:**
   - Download raw data from [NCBI SRA PRJNA770746](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770746).
   - Process data as per the scripts.

4. **Run Analysis:**
   - Use the main pipeline script (o1_automated_analysis_scripts) for taxonomic, diversity, and machine learning analyses.

---

## Data Availability

- **Primary Dataset:** [NCBI SRA PRJNA770746](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770746)
- **Processed Data & Metadata:** See Supplementary files (1A, 1B, 2A, 2B) in main publication.

---

## Main Findings

- Salivary microbiome profiling offers a non-invasive, accessible biomarker for AD.
- Distinct taxonomic (e.g., increased *Taylorella asinigenitalis*) and functional (e.g., reduced glycolysis) signatures differentiate AD from controls.
- Machine learning models based on microbiome features achieve diagnostic accuracy comparable to current clinical standards.
- AD is associated with increased oral microbial diversity, metabolic inflexibility, and fragmented microbial community structure.

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for improvements, bug fixes, or new analyses.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

**Corresponding Author:**  
Edward Jenner Tettevi  
Email: ejtettevi@gmail.com
