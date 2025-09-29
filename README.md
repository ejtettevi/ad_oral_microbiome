Oral Microbiome Dysbiosis in Alzheimer’s Disease: A Comprehensive Analysis**

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

## **6. Acknowledgments**

- **Authors**: Edward Jenner Tettevi, Samuel Armoo, David Larbi Simpong, Mike Y. Osei - Atweneboana.
- **Contact**: Correspondence to: Edward Jenner Tettevi, ejtettevi@gmail.com

