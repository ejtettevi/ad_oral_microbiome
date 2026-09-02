#!/usr/bin/env Rscript
# Analyze BioProject sample selection and create metadata tables

library(dplyr)
library(tidyr)
library(readr)

output_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(output_dir, recursive = TRUE)

cat("=== BioProject PRJNA770746 Sample Analysis ===\n\n")

run_table_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv"
run_table <- read_csv(run_table_path)

cat("Total SRR runs in BioProject:", nrow(run_table), "\n")

run_table$SampleType <- ifelse(run_table$isolation_source == "Oral_cavity", "Oral", 
                                ifelse(run_table$isolation_source == "feces", "Fecal", "Unknown"))

cat("\n=== Sample Types ===\n")
print(table(run_table$SampleType))

cat("\n=== By Sample Type and Disease ===\n")
print(table(run_table$SampleType, run_table$disease))

current_meta_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"
current_meta <- read.table(current_meta_path, sep = "\t", header = TRUE)
current_samples <- current_meta$Sample

run_table$Included_in_Analysis <- ifelse(run_table$Run %in% current_samples, "Yes", "No")

cat("\n=== Samples Included in Current Analysis ===\n")
included <- run_table[run_table$Included_in_Analysis == "Yes", ]
cat("Total included:", nrow(included), "\n")
cat("Of which are:", table(included$SampleType), "\n")

cat("\n=== Samples NOT Included (Available but not used) ===\n")
not_included <- run_table[run_table$Included_in_Analysis == "No", ]
cat("Total not included:", nrow(not_included), "\n")
cat("Of which are:", table(not_included$SampleType), "\n")

oral_not_used <- not_included[not_included$SampleType == "Oral", ]
fecal_all <- run_table[run_table$SampleType == "Fecal", ]

cat("\n=== Oral samples NOT used in current analysis ===\n")
cat("Count:", nrow(oral_not_used), "\n")
cat("Cases (AD):", sum(oral_not_used$disease == "AD"), "\n")
cat("Controls:", sum(oral_not_used$disease == "Control"), "\n")

cat("\n=== Fecal samples ===\n")
cat("Total:", nrow(fecal_all), "\n")
cat("Cases (AD):", sum(fecal_all$disease == "AD"), "\n")
cat("Controls:", sum(fecal_all$disease == "Control"), "\n")

all_oral <- run_table[run_table$SampleType == "Oral", ]
cat("\n=== All Oral samples in BioProject ===\n")
cat("Total:", nrow(all_oral), "\n")
cat("Cases (AD):", sum(all_oral$disease == "AD"), "\n")
cat("Controls:", sum(all_oral$disease == "Control"), "\n")

sample_table <- run_table[, c("Run", "Sample Name", "SampleType", "disease", "AGE", "sex", 
                               "BMI", "moca", "Included_in_Analysis", "BioSample")]
colnames(sample_table) <- c("SRR_ID", "Sample_Name", "Sample_Type", "Disease_Status", 
                              "Age", "Sex", "BMI", "MoCA", "Included", "BioSample")

write_csv(sample_table, file.path(output_dir, "sample_selection_table.csv"))

cat("\n=== Summary Table ===\n")
summary_df <- data.frame(
  Category = c("Total BioProject samples", 
               "Oral samples in BioProject",
               "Fecal samples in BioProject",
               "Oral samples used in analysis",
               "Fecal samples used in analysis"),
  Count = c(nrow(run_table), 
            sum(run_table$SampleType == "Oral"),
            sum(run_table$SampleType == "Fecal"),
            sum(run_table$SampleType == "Oral" & run_table$Included_in_Analysis == "Yes"),
            sum(run_table$SampleType == "Fecal" & run_table$Included_in_Analysis == "Yes"))
)
print(summary_df)

write_csv(summary_df, file.path(output_dir, "bioproject_summary.csv"))

cat("\n=== Sample Selection Table saved to ===\n")
cat(output_dir, "/sample_selection_table.csv\n")
cat(output_dir, "/bioproject_summary.csv\n")

cat("\n\n=== NEXT STEPS: Running fecal-only and combined analyses ===\n")
cat("This requires processing the fecal sample FASTQ files\n")
cat("Current analysis used only", sum(run_table$SampleType == "Oral" & run_table$Included_in_Analysis == "Yes"), "oral samples\n")