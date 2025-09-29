Okay
,


# **GitHub README: Oral Microbiome Dysbiosis in Alzheimer’s Disease: A Comprehensive Analysis**

---

## **1. Introduction**

Oral microbiome dysbiosis is increasingly being linked to the pathogenesis of Alzheimer’s disease (AD). This study investigates the relationship between oral microbiome alterations and AD, using high-resolution 16S rRNA sequencing (HRS) and machine learning to identify dysbiosis-associated metabolic signatures.

The findings reveal that oral microbiome dysbiosis in AD is characterized by:

- **Expansion of Proteobacteria** and **depletion of Actinobacteria** in AD cases.
- **Enrichment of Bacilli** and **Gammaproteobacteria** at the functional level.
- **Coordination of metabolic reprogramming** in AD-associated microbiomes.
- **Dysbiotic signatures** indicating disrupted microbial homeostasis and opportunistic colonization.

These results highlight the potential for microbiome-derived signatures as non-invasive biomarkers for AD diagnosis and monitoring.

---

## **2. Methods**

### **1.1. Study Design and Data Acquisition**

- **Oral Microbiome Composition**: Cases and controls were analyzed using HRS to investigate oral microbiome dysbiosis.
- **Data Sources**: Illumina 16S rRNA sequencing data from Vancouver, Canada (BioProject PRJNA770746).
- **Data Preprocessing**: Raw data and metadata were obtained from the NCBI Sequence Read Archive. The pipeline included: filtering, trimming, and quality control.

### **1.2. Diversity and Statistical Analyses**

- **Alpha Diversity**: Chao1 and Fisher diversity indices were computed.
- **Beta Diversity**: Bray-Curtis and Jaccard distances were used to assess community structure.
- **Functional Profiling**: Compositional differential abundance analysis was performed using ANCOM-BC and ALDEx2.

### **1.3. Machine Learning Classification**

- **Model Selection**: Random Forest, Gradient Boosting, Logistic Regression, and SVM were implemented.
- **Model Evaluation**: Metrics included AUC, balanced accuracy, and ROC-AUC.
- **Visualization**: All plots were generated using ggplot2, seaborn, matplotlib, and pheatmap.

---

## **3. Results**

### **3.1. Microbial Composition and Functional Profiling**

- **Phylum Level**: Cases showed a marked expansion of Proteobacteria and depletion of Actinobacteria.
- **Class Level**: Bacilli were enriched in cases, while Gammaproteobacteria were increased.
- **Order Level**: Lactobacillales were reduced in cases compared to controls.
- **Family Level**: Streptococcus and Haemophilus were reduced in cases compared to controls.
- **Genus Level**: Streptococcus and Veillonella consistently ranked as the most abundant genera in both groups.

### **3.2. Functional Profiling and Metabolic Pathways**

- **Differential Abundance**: Coordinated upregulation of biosynthetic and energy metabolism pathways was observed.
- **Metabolic Pathways**: Pathways involved in phospholipid metabolism and peptidoglycan maturation were significantly upregulated in cases.

### **3.3. Machine Learning and Classification**

- **Model Evaluation**: Random Forest achieved robust case-control discrimination with AUC = 0.85.
- **Confusion Matrix**: Random Forest showed optimal performance with 77.8% sensitivity and 81.8% specificity.

---

## **4. Discussion**

### **4.1. Dysbiotic Signatures and Pathogenesis**

- **Dysbiotic Signature**: Cases exhibited community destabilization marked by opportunistic colonization and pro-inflammatory species.
- **Link to AD**: The observed microbial restructuring aligns with the hypothesis that AD pathogenesis is linked to a loss of oral microbial homeostasis.

---

## **5. Conclusion**

This study integrates high-resolution taxonomic, ecological, and predicted functional analyses to characterize oral microbiome dysbiosis in AD. The results highlight the potential for microbiome-derived signatures as non-invasive biomarkers for AD diagnosis and monitoring.

---

## **6. References**

1. Anderson, M. J. (2017). *Permutational Multivariate Analysis of Variance (PERMANOVA)*. 1–15.
2. Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *The Journal of the Royal Statistical Society, Series B* (Statistical Methodology), 57(1), 289–300.
3. Bradley, A. (1997). *The use of the area under the ROC curve in the evaluation of machine learning algorithms*. *Pattern Recognition*, 30(1145–1159).
4. Callahan, B., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J., & Holmes, S. (2016). *DADA2: High resolution sample inference from Illumina amplicon data*. *Nature Methods*, 13(581–583).
5. Caspi, R., Billington, R., Keseler, I. M., Kothari, A., Krummenacker, M., Midford, P. E., Ong, W. K., Paley, S., Subhraveti, P., & Karp, P. D. (2020). *The MetaCyc database of metabolic pathways and enzymes - a 2019 update*. *Nucleic Acids Research*, 48(D1), D445–D453.
6. Chen, S. (2023). *Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp*. *Imeta*, 2(2), e107.
7. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). *fastp: An ultra-fast all-in-one FASTQ preprocessor*. *Bioinformatics*, 34(17), i884–i890.
8. Dixon, P. (2003). *VEgan: A package of R functions for community ecology*. *14*, 927–930.
9. Long, H., Yan, L., Pu, J., Liu, Y., Zhong, X., & Wang…, H. (2022a). *Multi-omics analysis reveals the effects of microbiota on oral homeostasis*. *Frontiers in …*. https://doi.org/10.3389/fimmu.2022.1005992/pdf
10. Long, H., Yan, L., Pu, J., Liu, Y., Zhong, X., Wang, H., Yang, L., Lou, F., Luo, S., Zhang, Y., Liu, Y., Xie, P., Ji, P., & Jin, X. (2022b). *Multi-omics analysis reveals the effects of microbiota on oral homeostasis*. *Frontiers in Immunology*. https://doi.org/10.3389/fimmu.2022.1005992

---

## **7. Acknowledgments**

- **Authors**: Edward Jenner Tettevi, Samuel Armoo, David Larbi Simpong, Mike Y. Osei - Atweneboana.
- **Data**: Available on GitHub: [https://github.com/ejtettevi/ad_oral_microbiome](https://github.com/ejtettevi/ad_oral_microbiome)
- **Contact**: Correspondence to: Edward Jenner Tettevi, ejtettevi@gmail.com

--- 

This README provides a comprehensive overview of the study, highlighting its methodology, results, and implications for understanding oral microbiome dysbiosis in AD.
