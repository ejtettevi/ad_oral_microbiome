# Oral Microbiome Dysbiosis Characterizes Alzheimer's Disease: Integrated Taxonomic, Functional, and Machine Learning Insights from a Vancouver Cohort

**Journal:** *The Microbe*  
**Status:** Revised Manuscript Submission  
**BioProject:** PRJNA770746  
**Contact:** Edward Jenner Tettevi (ejtettevi@gmail.com)

---

## Overview

This repository contains the analysis code, data, and manuscript files for a study investigating oral microbiome dysbiosis in Alzheimer's disease (AD) using high-resolution 16S rRNA sequencing and machine learning approaches.

### Key Findings

- **Expansion of Proteobacteria** and **depletion of Actinobacteria** in AD cases
- **Enrichment of Bacilli** and **Gammaproteobacteria** at the functional level
- **Coordination of metabolic reprogramming** in AD-associated microbiomes
- **Dysbiotic signatures** indicating disrupted microbial homeostasis and opportunistic colonization
- Machine learning models (Gradient Boosting AUC = 0.83, Random Forest AUC = 0.56) for AD classification

---

## Repository Structure

```
├── README.md                    # This file
├── LICENSE                      # MIT License
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies
├── sample_metadata.csv          # Sample metadata for all runs
│
├── scripts/                     # Analysis pipeline scripts
│   ├── 1_download_sra_data.sh
│   ├── 2_run_fastp_automatic.sh
│   ├── 3_run_cutadapt_v3v4.sh
│   ├── 4_run_bowtie2_decontam.sh
│   ├── 5_dada2_16S_pipeline_standalone.R
│   ├── 6_abundance_analysis.R
│   ├── 7_alpha_beta_analysis.py
│   ├── 8_microbiome_ml.py
│   └── 9_picrust.py
│
├── dada2/                       # DADA2 pipeline outputs
│   ├── asv_table.tsv
│   ├── asv_taxonomy.tsv
│   ├── phyloseq_singletonfree.rds
│   ├── read_tracking.tsv
│   ├── rep_seqs.fasta
│   ├── seqtab_nochim_singletonfree.rds
│   └── taxa_singletonfree.rds
│
├── manuscript/                  # Manuscript versions
│   ├── Maniuscript_Microbe_v01.docx
│   ├── Maniuscript_Microbe_v1.docx
│   ├── Maniuscript_Microbe_v2.docx
│   ├── Maniuscript_Microbe_v3.docx
│   ├── Maniuscript_Microbe_v3_FINAL.docx
│   ├── Maniuscript_Microbe_v4_FINAL.docx
│   ├── Oral_AD_Manuscript_v3.txt
│   └── Enhanced_Abstract.txt
│
├── figures/                     # All manuscript figures
│   ├── Fig. 2.jpg / Fig. 2.eps          # PCoA + ANCOM-BC2 LFC scores
│   ├── Fig. 3.jpg / Fig. 4.jpg / Fig. 5.jpg  # Alpha/Beta diversity
│   ├── Fig._6A.png / Fig._6B.png        # ML classification results
│   ├── Fig._7.png                        # Feature importance
│   ├── ml_*.png                          # ML model ROC/CM plots
│   ├── oral_pcoa_bray.png               # Beta diversity PCoA
│   ├── rf_feature_importance.png        # RF feature importance
│   └── Supplementary_Fig*.png           # Supplementary figures
│
├── results/                     # Analysis result tables
│   ├── ancombc2_*_results.csv            # ANCOM-BC2 results (5 taxonomic ranks)
│   ├── oral_ancombc2_results.csv         # Combined oral ANCOM-BC2 results
│   ├── combined_beta_pca.csv             # Beta diversity PCA coordinates
│   ├── fecal_alpha_diversity.csv         # Fecal alpha diversity metrics
│   ├── fecal_beta_pca.csv                # Fecal beta diversity PCA
│   ├── fecal_differential_abundance.csv  # Fecal differential abundance
│   └── rf_feature_importance.csv         # RF Gini importance scores
│
├── data/                        # Cohort and sample data
│   ├── bioproject_summary.csv
│   ├── demographics_summary_table.csv
│   ├── complete_runs_inclusion_table.csv
│   └── sample_selection_table.csv
│
├── review/                      # Peer review responses
│   ├── response_to_all_reviewers.docx
│   ├── response_to_all_reviewers.md
│   ├── response_to_reviewers.md
│   ├── reviewer_response_table.docx
│   └── reviewer_response_table.md
│
└── supplementary/               # Supplementary materials
    ├── Supplementary.docx
    └── Supplementary.html
```

---

## Methods Summary

### Study Design
- **Cohort:** 198 sequencing runs from BioProject PRJNA770746 (99 oral, 99 fecal)
- **Oral samples:** 98 passed QC (1 excluded as technical replicate)
- **Fecal samples:** 78 passed QC (21 excluded due to poor quality)
- **Sequencing:** Illumina 16S rRNA (V3-V4 region)

### Analysis Pipeline
1. **Quality Control:** FastP trimming, Cutadapt adapter removal, Bowtie2 decontamination
2. **ASV Processing:** DADA2 with singleton removal
3. **Taxonomic Profiling:** SILVA database (genus-level primary analysis)
4. **Diversity Analysis:** Alpha (Chao1, Fisher, Observed) and Beta (Bray-Curtis, Jaccard)
5. **Differential Abundance:** ANCOM-BC2 with CLR log-fold changes
6. **Functional Prediction:** PICRUSt2 pathway inference
7. **Machine Learning:** Gradient Boosting, Random Forest, Logistic Regression, SVM

### Key Statistical Results
- **Alpha Diversity:** Chao1 genus richness p = 0.003, Fisher p = 0.009
- **Beta Diversity:** PERMANOVA R² = 1.89%, F = 6.24, p = 0.003 (disease effect)
- **ML Best Model:** Gradient Boosting AUC = 0.83, Accuracy = 75.0%
- **Top Biomarkers:** Gemella (Gini = 0.035), Haemophilus (0.030), Neisseria oralis (0.022)

---

## Reproducibility

### Environment Setup
```bash
# Using Conda
conda env create -f environment.yml
conda activate ad_oral_microbiome

# Using pip
pip install -r requirements.txt
```

### Data Access
- Raw sequences: NCBI SRA BioProject PRJNA770746
- Processed ASV table: `dada2/asv_table.tsv`
- Sample metadata: `sample_metadata.csv`

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this code or data, please cite:

```
Tettevi, E.J., Armoo, S., Simpong, D.L., Osei-Atweneboana, M.Y. (2024).
Oral Microbiome Dysbiosis Characterizes Alzheimer's Disease: Integrated Taxonomic,
Functional, and Machine Learning Insights from a Vancouver Cohort.
The Microbe.
```

---

## Contact

Edward Jenner Tettevi  
Email: ejtettevi@gmail.com  
GitHub: [@ejtettevi](https://github.com/ejtettevi)
