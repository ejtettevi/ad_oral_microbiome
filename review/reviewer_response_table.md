# Reviewer Response Table — Maniuscript_Microbe_v2.docx

| Comment | Reviewer Request | Manuscript Location | Action Taken |
|---------|-----------------|---------------------|--------------|
| **R1.C1** | Provide metadata table of included SRR runs | p.6, l.23-25 | P49: Table 1 with run inclusions/exclusions and demographics |
| | Perform fecal-only analysis on same cohort | p.17, l.1 | P123–130: Fecal-only microbiome analysis section added |
| | Perform combined oral+fecal analysis | p.18, l.13 | P132–136: Combined multi-niche comparison with PCA |
| **R1.C2** | Provide ASV table from DADA2 pipeline | p.8, l.11-23 | P60: DADA2 methods described; ASV table provided as supplementary material |
| **R1.C3** | Include HOMD 16S rRNA reference for taxonomy assignment | p.8, l.11-23 | P60: Added HOMD v15.2 assignment with priority for oral taxa over SILVA |
| **R1.C4** | Perform species-level taxonomy assignment | p.8, l.11-23 | P60: `addSpecies` function from DADA2 used for species-level calls |
| | | p.21, l.10-18 | P147: Limitations of 16S species-level resolution acknowledged |
| **R1.C5** | Upgrade to ANCOM-BC2; report CLR values and adjusted *q*-values | p.9, l.17-25 | P66: ANCOM-BC2 v2.0.0 specified with CLR-transformed LFC and BH-adjusted *q*-values |
| | | p.11, l.15-21 | P80–82: All results report LFC + *q*-values |
| **R1.C6** | Replace percentage bar charts with CLR/LFC-based reporting | p.11, l.8-13 | P79: Abundance changes reported as LFC + *q*-values, not percentages |
| | | p.12, l.16-18 | P86: Fig. 2 now states ANCOM-BC2 CLR-based analysis |
| **R1.C7** | Correct "species richness" → "genus richness"; state taxonomic level | p.9, l.4-12 | P63: Alpha diversity metrics computed at genus level |
| | | p.12, l.23-29 | P89: "Chao1 estimator showed greater **genus richness**" |
| **R1.C8** | State taxonomic level for beta-diversity analyses | p.9, l.4-12 | P63: Beta diversity at genus level using CLR-transformed abundances |
| | | p.13, l.10-15 | P95: "Beta diversity analyses **at the genus level**" |
| **R1.C9** | Remove "gene expression" language; replace "up/down regulated" with "increased/decreased inferred abundance" | p.9, l.17-25 | P66: "predicted functional potential…not actual gene expression levels" |
| | | p.3, l.1-10 | P19: "coordinated **increase in predicted biosynthetic pathway potential**" |
| | | p.14, l.11-21 | P103: "coordinated **increase in predicted biosynthetic pathway abundance**" |
| | | p.20, l.11-19 | P142: "coordinated **increase in predicted biosynthetic pathway abundance**" |
| | | p.22, l.4-11 | P151: "**increased abundance of predicted** biosynthetic pathways" |
| **R1.C10** | Perform ML on CLR-transformed data; report most influential taxa | p.10, l.4-13 | P69: CLR-transformed ASV abundances used with SelectFromModel |
| | | p.15, l.22-29 | P110: Top features — *Tannerella* (~0.039), *Gemella* (~0.038), *Haemophilus* (~0.029) |
| **R6.C1** | 8 keywords avoiding title-word overlap | p.3, l.20-22 | P22: 8 keywords (16S rRNA sequencing; ANCOM-BC2; biomarker discovery; compositional data analysis; neuroinflammation; PICRUSt2; predicted functional potential; Random Forest) |
| **R6.C2** | Add future research recommendations | p.21, l.20-29 | P148: Longitudinal tracking, multi-omics integration, mechanistic validation |
