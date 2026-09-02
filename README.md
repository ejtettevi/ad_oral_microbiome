<div align="center">

# 🦠 Oral Microbiome Dysbiosis is Associated with Alzheimer's Disease

### Integrated Taxonomic, Functional & Machine-Learning Insights from a Vancouver Cohort

**Edward Jenner Tettevi** · **David Larbi Simpong** · **Robert Peter Biney** · **Mike Y. Osei-Atweneboana**

[![BioProject](https://img.shields.io/badge/BioProject-PRJNA770746-43853D?style=flat-square&logo=dna&logoColor=white)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770746)
[![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=flat-square&logo=python&logoColor=white)]()
[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)]()
[![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)]()
[![16S](https://img.shields.io/badge/Amplicon-16S%20rRNA-blueviolet?style=flat-square)]()
[![ML](https://img.shields.io/badge/AUC-0.85-FF6F00?style=flat-square)]()

_A reproducible 16S rRNA amplicon metagenomics pipeline linking oral microbiome structure,
function, and diagnostic prediction to Alzheimer's disease._

</div>

---

## ✨ Highlights

- 🔬 **Secondary re-analysis** of publicly deposited 16S rRNA sequencing data (BioProject `PRJNA770746`)
- 🧬 **96 oral** (43 AD / 53 control) + **78 fecal** (32 AD / 46 control) = **174 samples analysed**
- 💊 **ANCOM-BC2** compositional differential-abundance testing across 5 taxonomic ranks
- 🧠 **PICRUSt2** predicted functional profiling (biosynthetic pathways enriched in AD)
- 🤖 **Random Forest** classifier reaching **AUC = 0.85** (95% CI 0.72–0.98)
- ⚖️ Fully reproducible end-to-end workflow (Python + R, conda + pip)

---

## 📖 Key Findings

| Finding | Detail | Status |
|---|---|---|
| **Phylum enrichment** | `Bacteroidota` enriched in AD | ✅ significant (LFC +0.75, q = 0.037) |
| **Directional shifts** | `Proteobacteria` ↑, `Actinobacteriota` ↓ | directional, non-significant |
| **Genus markers** | `SUP05`, `Gemella`, `Fusobacterium` ↑; `F0332`, `Selenomonas` ↓ | significant |
| **Alpha diversity** | Increased in AD (community destabilization) | Chao1 p = 0.003 |
| **Functional potential** | Biosynthetic pathways ↑ (peptidoglycan, thiamine, phospholipid) | 58% rise |
| **Diagnostic ML** | Random Forest **AUC 0.85** | sens. 77.8%, spec. 81.8% |

---

## 🩺 Study Design

| Metric | Oral | Fecal | Combined |
|---|---|---|---|
| **Enrolled participants** | 99 | 99 | 99 |
| **Samples analysed** | **96** (43 case / 53 ctrl) | **78** (32 case / 46 ctrl) | **174** |
| **Excluded** | 3 (low depth/QC) | 21 | 24 |
| **Sequencing** | Illumina MiSeq, 16S V3–V4 | Illumina MiSeq, 16S V3–V4 | — |

> ✅ Corrected to the finalized analysis: **96 oral / 78 fecal / 174 combined**.
> Arc-shaped body-site separation in combined PCA driven by sample site (PC1 = 49.6%).

---

## 🔬 Methods Pipeline

```
 raw FASTQ (SRA)
     │  1_download_sra_data.sh
     ▼
 Quality Control ─────────── 2_run_fastp_automatic.sh  (fastp, Q≥20)
     │                       3_run_cutadapt_v3v4.sh     (adapters, primers)
     ▼
 Host Decontamination ────── 4_run_bowtie2_decontam.sh  (hg38; bowtie2, samtools)
     │
     ▼
 ASV Inference ───────────── 5_dada2_16S_pipeline_standalone.R  (DADA2, SILVA v138.1, HOMD)
     │
     ├──► Taxonomy & Abundance ── 6_abundance_analysis.R
     ├──► Alpha/Beta Diversity ── 7_alpha_beta_analysis.py   (vegan, scikit-bio)
     ├──► Function (PICRUSt2) ──── 9_picrust.py
     └──► Machine Learning ─────── 8_microbiome_ml.py        (scikit-learn)
```

---

## 🖼️ Figures

### Combined Oral + Fecal Microbiome (PCA)

<p align="center">
  <img src="figures/Fig_8_combined_pca.png" width="70%">
</p>

> 174 samples, genus-level CLR-transformed abundances. Separation is driven primarily by
> **body site** (PC1 = 49.6%) with a significant independent disease effect (PERMANOVA R² = 1.03%, p = 0.020).

### Machine Learning Classification

<p align="center">
  <img src="figures/ml_multipanel_roc_curve.png" width="55%">
  <img src="figures/ml_RandomForest_confusion_matrix.png" width="38%">
</p>

### Beta Diversity & Feature Importance

<p align="center">
  <img src="figures/oral_pcoa_bray.png" width="45%">
  <img src="figures/rf_feature_importance.png" width="45%">
</p>

---

## 📂 Repository Structure

```
ad_oral_microbiome/
├── README.md                  ← You are here
├── LICENSE                    MIT License
├── environment.yml            Conda environment (Python + R + tools)
├── requirements.txt           Python dependencies
├── sample_metadata.csv        Sample metadata for all runs
│
├── scripts/                   # All 36 analysis scripts used
│   │
│   │  Core pipeline (numbered)
│   ├── 1_download_sra_data.sh
│   ├── 2_run_fastp_automatic.sh
│   ├── 3_run_cutadapt_v3v4.sh
│   ├── 4_run_bowtie2_decontam.sh
│   ├── 5_dada2_16S_pipeline_standalone.R
│   ├── 6_abundance_analysis.R
│   ├── 7_alpha_beta_analysis.py
│   ├── 8_microbiome_ml.py
│   ├── 9_picrust.py
│   ├── 10_pathway-tools.py
│   │
│   │  Fecal processing
│   ├── download_and_process_fecal.sh
│   ├── process_fecal_samples.R
│   ├── process_fecal_batch.R
│   ├── fecal_analysis_status.R
│   ├── fecal_dada2_analysis.R
│   ├── fecal_dada2_simple.R
│   ├── combined_analysis.R
│   │
│   │  (Batch) improved analyses
│   ├── analyze_sample_selection.R
│   ├── identify_excluded_samples.R
│   ├── run_improved_analysis.R
│   ├── improved_analysis_clean.R
│   ├── improved_analysis_final.R
│   ├── improved_analysis_v2.R
│   ├── improved_differential_analysis.R
│   ├── improved_visualization.R
│   ├── full_analysis.R
│   │
│   │  HOMD taxonomy / species-level assignment
│   ├── setup_homd.R
│   ├── homd_fast_taxonomy.R
│   ├── homd_priority_taxonomy.R
│   ├── run_homd_species.R
│   ├── species_level_assignment.R
│   │
│   │  Additional utilities
│   ├── assign_fecal_taxonomy.R
│   ├── generate_manuscript_updates.py
│   ├── run_fecal_and_combined_analysis.R
│   │
│   │  Corrected final-analysis scripts
│   ├── run_ancombc2_phylum.R
│   └── regenerate_combined.R
│
├── dada2/                     # DADA2 pipeline outputs
├── data/                      # Cohort & sample-level data (CSV)
├── results/                   # Analysis result tables (CSV)
├── figures/                   # All figures (PNG / JPG / EPS)
├── manuscript/                # Manuscript versions (.docx)
├── review/                    # Reviewer-response materials
├── supplementary/             # Supplementary materials
└── env/                       # (conda env, created locally)
```

---

## 📊 Data Files (`.csv`)

### Cohort & Sample Metadata — `data/`

| File | Description |
|---|---|
| `bioproject_summary.csv` | BioProject-level summary |
| `sample_selection_table.csv` | **198-run** inclusion / exclusion status |
| `complete_runs_inclusion_table.csv` | Full SRR run table with demographics |
| `demographics_summary_table.csv` | Age / sex / BMI / MoCA by group (96-oral route) |

### Analysis Results — `results/`

| File | Description |
|---|---|
| `ancombc2_phylum_results.csv` | ANCOM-BC2 **phylum-level** results (Bacteroidota q = 0.037) |
| `ancombc2_class_results.csv` | Class-level results |
| `ancombc2_order_results.csv` | Order-level results |
| `ancombc2_family_results.csv` | Family-level results |
| `ancombc2_genus_results.csv` | Genus-level results |
| `oral_ancombc2_results.csv` | Combined oral ANCOM-BC2 output |
| `combined_beta_pca.csv` | **174-sample** combined oral+fecal PCA coordinates |
| `combined_beta_pca_174.csv` | Same as above (explicit label) |
| `fecal_alpha_diversity.csv` | Fecal alpha-diversity metrics |
| `fecal_beta_pca.csv` | Fecal beta-diversity PCA |
| `fecal_differential_abundance.csv` | Fecal differential abundance |
| `rf_feature_importance.csv` | Random Forest Gini importance scores |

---

## 🖌️ Figure Files (`.png` / `.eps` / `.jpg`)

### Main Figures — `figures/`

| File | Figure |
|---|---|
| `Fig. 2.eps` / `Fig. 2.jpg` | PCoA + ANCOM-BC2 log-fold changes (vector EPS) |
| `Fig. 3.jpg` | Chao1 / Fisher / Observed alpha diversity |
| `Fig. 4.jpg` | Beta diversity (NMDS / PCoA) |
| `Fig. 5.jpg` | PICRUSt2 pathway analysis |
| `Fig._6A.png` | ML ROC curves (all models) |
| `Fig._6B.png` | Confusion matrices |
| `Fig._7.png` | Random Forest feature importance |
| `Fig_8_combined_pca.png` | **Combined oral+fecal PCA (174 samples)** |

### Supplementary Figures — `figures/`

| File | Description |
|---|---|
| `Supplementary_Fig_2C.png` | Order-level differential abundance |
| `Supplementary_Fig_2D.png` | Family-level differential abundance |
| `Supplementary_Fig_2E.png` | Genus-level differential abundance |
| `Supplementary_Fig_5C.png` | Pathway volcano plot |

### Model & Utility Plots — `figures/`

| File | Description |
|---|---|
| `ml_multipanel_roc_curve.png` | All-model ROC curves (multipanel) |
| `ml_multipanel_confusion_matrix.png` | All-model confusion matrices |
| `ml_RandomForest_roc_curve.png` / `..._confusion_matrix.png` | RF ROC / confusion matrix |
| `ml_GradientBoosting_roc_curve.png` / `..._confusion_matrix.png` | GBM ROC / confusion matrix |
| `ml_LogisticRegression_roc_curve.png` / `..._confusion_matrix.png` | LogReg ROC / confusion matrix |
| `ml_SVM_roc_curve.png` / `..._confusion_matrix.png` | SVM ROC / confusion matrix |
| `oral_pcoa_bray.png` | Oral beta diversity PCoA |
| `rf_feature_importance.png` | RF feature importance |

---

## 🔁 Reproducibility

### Environment Setup

```bash
# Option A — Conda (recommended)
conda env create -f environment.yml
conda activate oral_microbiome_analysis

# Option B — pip (Python-only)
pip install -r requirements.txt
```

### Data Access

- **Raw reads:** NCBI SRA — BioProject [`PRJNA770746`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770746)
- **ASV table:** `dada2/asv_table.tsv`
- **Taxonomy:** `dada2/asv_taxonomy.tsv`
- **Metadata:** `sample_metadata.csv`

### Run order

```bash
# ── 1. Download ────────────────────────────────────────────
bash scripts/1_download_sra_data.sh
bash scripts/download_and_process_fecal.sh          # fecal FASTQ download

# ── 2. QC & Decontamination ───────────────────────────────
bash scripts/2_run_fastp_automatic.sh
bash scripts/3_run_cutadapt_v3v4.sh
bash scripts/4_run_bowtie2_decontam.sh

# ── 3. DADA2 ASV inference ────────────────────────────────
Rscript scripts/5_dada2_16S_pipeline_standalone.R
Rscript scripts/fecal_dada2_analysis.R              # fecal ASV processing

# ── 4. Taxonomy & Abundance ──────────────────────────────
Rscript scripts/6_abundance_analysis.R
Rscript scripts/setup_homd.R                        # HOMD database setup
Rscript scripts/homd_fast_taxonomy.R                # HOMD-based assignment
Rscript scripts/species_level_assignment.R          # species-level refinement

# ── 5. Sample selection & filtering ───────────────────────
Rscript scripts/analyze_sample_selection.R
Rscript scripts/identify_excluded_samples.R

# ── 6. Diversity & Differential Abundance ─────────────────
python scripts/7_alpha_beta_analysis.py
Rscript scripts/improved_differential_analysis.R

# ── 7. Functional prediction ──────────────────────────────
python scripts/9_picrust.py
python scripts/10_pathway-tools.py

# ── 8. Machine learning ──────────────────────────────────
python scripts/8_microbiome_ml.py

# ── 9. Combined analysis & figures ────────────────────────
Rscript scripts/run_fecal_and_combined_analysis.R
Rscript scripts/regenerate_combined.R               # corrected 174-sample PCA
Rscript scripts/run_ancombc2_phylum.R               # corrected phylum results
```

> **Note:** The 22 additional scripts in `scripts/` include batch-processing helpers,
> HOMD taxonomy pipelines, improved-analysis iterations, and corrected final-analysis
> scripts used to produce the results in this repository.

---

## 📜 Citation

If you use this pipeline, data, or results, please cite:

```bibtex
@article{tettevi2026oral,
  title   = {Oral Microbiome Dysbiosis is Associated with Alzheimer's Disease:
             Integrated Taxonomic, Functional, and Machine-Learning Insights
             from a Vancouver Cohort},
  author  = {Tettevi, Edward Jenner and Simpong, David Larbi and
             Biney, Robert Peter and Osei-Atweneboana, Mike Y.},
  journal = {The Microbe},
  year    = {2026},
  note    = {Revised manuscript submission}
}
```

---

## 📬 Contact

**Edward Jenner Tettevi** — [ejtettevi@gmail.com](mailto:ejtettevi@gmail.com)

<p align="center"><sub>MIT License · Built with reproducibility in mind 🧬</sub></p>
