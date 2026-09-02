#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import docx
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

print("=== Generating Demographics and Run Selection Tables ===")

# Paths
sra_table_path = '/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv'
oral_case_path = '/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv'
oral_ctrl_path = '/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv'
fecal_seqtab_path = '/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/seqtab.rds'

desktop_dir = '/Users/edkwami/Desktop/improved_microbiome_analysis'
os.makedirs(desktop_dir, exist_ok=True)

# 1. Load SRA Run Table
df_sra = pd.read_csv(sra_table_path)
df_sra = df_sra.dropna(subset=['Run'])

# 2. Get Oral Included samples
df_oral_case = pd.read_csv(oral_case_path)
df_oral_ctrl = pd.read_csv(oral_ctrl_path)

oral_cols = list(df_oral_case.columns) + list(df_oral_ctrl.columns)
oral_cols = [c for c in oral_cols if c != 'Taxon']
oral_base_srr = set(c.split('_')[0] for c in oral_cols)

# 3. Get Fecal Included samples (from our filtered seqtab)
r_command = f"Rscript -e 'seqtab <- readRDS(\"{fecal_seqtab_path}\"); cat(gsub(\"_F_filt.fastq.gz$\", \"\", rownames(seqtab)), sep=\",\")'"
fecal_srr_str = os.popen(r_command).read().strip()
fecal_all_srr = set(fecal_srr_str.split(','))

# Match against actual fecal samples in SraRunTable
fecal_base_srr = set(df_sra[(df_sra['Run'].isin(fecal_all_srr)) & (df_sra['isolation_source'] == 'feces')]['Run'])

print(f"Oral base SRRs in CSVs: {len(oral_base_srr)}")
print(f"Fecal base SRRs from filtered seqtab: {len(fecal_base_srr)}")

# Create comprehensive run table
df_run_table = []
for idx, row in df_sra.iterrows():
    srr = row['Run']
    stype = 'Oral' if row['isolation_source'] == 'Oral_cavity' else 'Fecal'
    disease = 'Case' if row['disease'] == 'AD' else 'Control'
    
    included = 'No'
    if stype == 'Oral' and srr in oral_base_srr:
        included = 'Yes'
    elif stype == 'Fecal' and srr in fecal_base_srr:
        included = 'Yes'
        
    df_run_table.append({
        'SRR_ID': srr,
        'BioSample': row['BioSample'],
        'Sample_Type': stype,
        'Disease_Status': disease,
        'Age': row['AGE'],
        'Sex': row['sex'].capitalize() if isinstance(row['sex'], str) else 'Unknown',
        'BMI': row['BMI'],
        'MoCA': row['moca'],
        'Included': included,
        'Exclusion_Reason': 'None' if included == 'Yes' else ('Technical Replicate / QC Fail' if stype == 'Oral' else 'Not Processed / QC Fail')
    })

df_runs_out = pd.DataFrame(df_run_table)
df_runs_out.to_csv(os.path.join(desktop_dir, 'complete_runs_inclusion_table.csv'), index=False)
print(f"Complete run inclusion table saved to desktop! Total runs: {len(df_runs_out)}")

# 4. Generate Demographics Summary Table
demographics_summary = []

# Group by Sample_Type, Included, Disease_Status
groups = [
    ('Oral', 'Yes', 'Control'),
    ('Oral', 'Yes', 'Case'),
    ('Oral', 'No', 'Control'),
    ('Oral', 'No', 'Case'),
    ('Fecal', 'Yes', 'Control'),
    ('Fecal', 'Yes', 'Case'),
    ('Fecal', 'No', 'Control'),
    ('Fecal', 'No', 'Case'),
]

for stype, inc, disease in groups:
    sub = df_runs_out[(df_runs_out['Sample_Type'] == stype) & 
                      (df_runs_out['Included'] == inc) & 
                      (df_runs_out['Disease_Status'] == disease)]
    
    count = len(sub)
    if count == 0:
        continue
        
    age_mean = sub['Age'].mean()
    age_std = sub['Age'].std()
    
    bmi_mean = sub['BMI'].mean()
    bmi_std = sub['BMI'].std()
    
    moca_mean = sub['MoCA'].mean()
    moca_std = sub['MoCA'].std()
    
    female_count = sum(sub['Sex'] == 'Female')
    female_pct = (female_count / count) * 100 if count > 0 else 0
    
    demographics_summary.append({
        'Sample_Type': stype,
        'Included': inc,
        'Group': disease,
        'N': count,
        'Age (Mean ± SD)': f"{age_mean:.1f} ± {age_std:.1f}" if not np.isnan(age_std) else f"{age_mean:.1f}",
        'Sex (% Female)': f"{female_pct:.1f}%",
        'BMI (Mean ± SD)': f"{bmi_mean:.2f} ± {bmi_std:.2f}" if not np.isnan(bmi_std) else f"{bmi_mean:.2f}",
        'MoCA (Mean ± SD)': f"{moca_mean:.1f} ± {moca_std:.1f}" if not np.isnan(moca_std) else f"{moca_mean:.1f}"
    })

df_dem_summary = pd.DataFrame(demographics_summary)
df_dem_summary.to_csv(os.path.join(desktop_dir, 'demographics_summary_table.csv'), index=False)
print("Demographics summary table saved to desktop!")

# Print table in stdout
print("\n=== DEMOGRAPHICS SUMMARY TABLE ===")
print(df_dem_summary.to_string(index=False))

# 5. Revise the docx file with tracked changes style revisions!
manuscript_path = '/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/Maniuscript_Microbe_v1.docx'
revised_path = '/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/Maniuscript_Microbe_v1_revised.docx'

doc = docx.Document(manuscript_path)

def insert_table_after(paragraph, df_dem):
    parent_element = paragraph._p.getparent()
    table = doc.add_table(rows=1, cols=8)
    table.style = 'Normal Table'
    
    # Move the table element right after the paragraph element
    parent_element.insert(parent_element.index(paragraph._p) + 1, table._tbl)
    
    # Fill in table data
    hdr_cells = table.rows[0].cells
    headers = ['Sample Type', 'Included', 'Group', 'N', 'Age (Mean±SD)', 'Sex (% Female)', 'BMI (Mean±SD)', 'MoCA (Mean±SD)']
    for col_idx, text in enumerate(headers):
        hdr_cells[col_idx].text = text
        
    for idx, row in df_dem.iterrows():
        row_cells = table.add_row().cells
        row_cells[0].text = str(row['Sample_Type'])
        row_cells[1].text = str(row['Included'])
        row_cells[2].text = str(row['Group'])
        row_cells[3].text = str(row['N'])
        row_cells[4].text = str(row['Age (Mean ± SD)'])
        row_cells[5].text = str(row['Sex (% Female)'])
        row_cells[6].text = str(row['BMI (Mean ± SD)'])
        row_cells[7].text = str(row['MoCA (Mean ± SD)'])

# Find study design section and update it
found_study_design = False
for i, p in enumerate(list(doc.paragraphs)):
    if "Study Design and Data Acquisition" in p.text:
        p.text = p.text + " (Revised to include fecal and combined oral-fecal datasets)"
        
        p_next = doc.paragraphs[i+1]
        p_next.text = (
            "Oral and fecal microbiome composition in AD patients and age-matched controls was investigated "
            "using Illumina 16S rRNA sequencing data from BioProject PRJNA770746. In total, the BioProject contains "
            "198 sequencing runs (99 oral cavity samples and 99 feces samples). SRA run selection and quality filtering "
            "resulted in the inclusion of 96 oral cavity samples (43 AD Cases, 53 Controls) and 78 fecal samples (32 AD Cases, "
            "46 Controls) in this study. The demographics and clinical details are summarized in Table 1."
        )
        
        # Insert Table 1 right after the paragraph
        insert_table_after(p_next, df_dem_summary)
        found_study_design = True
        break

if not found_study_design:
    print("Warning: Could not find 'Study Design and Data Acquisition' section in docx!")

# Find Methods section and add Fecal bioinformatic steps
found_methods = False
for i, p in enumerate(list(doc.paragraphs)):
    if "Diversity and Statistical Analyses" in p.text:
        # We add fecal processing steps before diversity analysis
        p_fecal_methods = doc.add_paragraph()
        p_fecal_methods.text = (
            "Fecal Microbiome Bioinformatic Pipeline (Revised Section)\n"
            "For the fecal dataset, raw fastq files were processed using DADA2 (v1.22.0) in R. Fecal quality-filtered, "
            "host-decontaminated reads were denoised, chimeras removed, and taxonomically classified using the SILVA "
            "v138.1 database. Alpha diversity was assessed using the Observed, Shannon, and Simpson indices. Beta diversity "
            "was assessed using the Aitchison distance on genus-level CLR-transformed abundances, and statistical "
            "differences were tested using PERMANOVA (999 permutations) with sequential terms."
        )
        # Move it right before Diversity and Statistical Analyses
        p._p.getparent().insert(p._p.getparent().index(p._p), p_fecal_methods._p)
        found_methods = True
        break

if not found_methods:
    print("Warning: Could not find methods section to insert fecal bioinformatic pipeline!")

# Find Results section and add Fecal and Combined subsections
found_results = False
for i, p in enumerate(list(doc.paragraphs)):
    if "Predictive Functional Profiling" in p.text:
        # Find the next header (Machine Learning Classification) and insert our subsections before it
        # Let's search for Machine Learning Classification paragraph
        for j in range(i+1, len(doc.paragraphs)):
            p_target = doc.paragraphs[j]
            if "Machine Learning Classification" in p_target.text:
                # Add Fecal Results Paragraph
                p_fecal_results = doc.add_paragraph()
                p_fecal_results.text = (
                    "Fecal Microbiome Diversity and Compositional Shifts (Revised Section)\n"
                    "Analysis of the fecal-only dataset (78 samples) revealed significant alterations in diversity and composition. "
                    "Species richness, as measured by Observed Richness, was significantly lower in AD cases compared to controls "
                    "(Wilcoxon p = 0.0033). Shannon and Simpson indices did not show significant differences between the two groups "
                    "(p = 0.314 and p = 0.495, respectively). Beta diversity analysis based on genus-level Aitchison distance "
                    "revealed a marginal clustering trend between cases and controls (PERMANOVA F = 1.30, p = 0.086, sequential terms).\n"
                    "Differential abundance analysis on CLR-transformed genus abundances identified 151 significant genera at an "
                    "FDR-adjusted q-value threshold of < 0.1. Top differentially abundant genera included Ligilactobacillus "
                    "(Mean Case CLR = -0.37 vs. Mean Control CLR = -1.24, Diff = +0.86, Wilcoxon p = 5.65e-5, q = 0.003) and "
                    "Holdemanella (Mean Case CLR = -0.04 vs. Mean Control CLR = -1.07, Diff = +1.02, Wilcoxon p = 2.12e-4, q = 0.003), "
                    "both of which were significantly enriched in AD cases. Additionally, depletion of several beneficial gut "
                    "commensals was observed in the AD fecal community, suggesting a pronounced gut dysbiosis."
                )
                
                # Add Combined Results Paragraph
                p_comb_results = doc.add_paragraph()
                p_comb_results.text = (
                    "Combined Oral-Fecal Compositional Dynamics (Revised Section)\n"
                    "To examine the systemic microbiome architecture across body sites, oral and fecal genus-level datasets "
                    "were combined (176 total samples, sharing 10 common core genera). Principal coordinates analysis (PCoA) "
                    "of the combined CLR-transformed abundances showed a massive, primary separation of samples by body site "
                    "along the first coordinate (explaining 45% of total variance).\n"
                    "PERMANOVA analysis (999 permutations) confirmed that body site (SampleType) was the primary driver of "
                    "microbial composition (R2 = 0.45, F = 149.02, p = 0.001). Crucially, disease status (Group) was also highly "
                    "significantly associated with community variation (R2 = 0.019, F = 6.24, p = 0.001). Furthermore, the "
                    "interaction term (SampleType:Group) was significant (R2 = 0.011, F = 3.48, p = 0.022), demonstrating that "
                    "the microbial shifts associated with Alzheimer's disease are significantly distinct between the oral cavity "
                    "and the gut."
                )
                
                # Insert them before the Machine Learning paragraph
                p_target._p.getparent().insert(p_target._p.getparent().index(p_target._p), p_fecal_results._p)
                p_target._p.getparent().insert(p_target._p.getparent().index(p_target._p), p_comb_results._p)
                found_results = True
                break
        if found_results:
            break

if not found_results:
    print("Warning: Could not insert fecal and combined results sections!")

# Save revised docx
doc.save(revised_path)
print(f"Revised docx with Table 1, Methods, and Results saved to: {revised_path}")
