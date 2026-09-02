#!/usr/bin/env Rscript
# Identify excluded samples and create updated metadata

run_table_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv"
current_meta_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"

run_table <- read.csv(run_table_path)
current_meta <- read.table(current_meta_path, sep = "\t", header = TRUE)

all_oral <- run_table[run_table$isolation_source == "Oral_cavity", ]
used_samples <- current_meta$Sample
excluded_oral <- all_oral[!all_oral$Run %in% used_samples, ]

cat("=== Excluded Oral Samples ===\n")
print(excluded_oral[, c("Run", "Sample Name", "disease", "AGE", "sex")])

cat("\nThese 3 oral samples were NOT included in the analysis.\n")
cat("Now creating comprehensive metadata with all samples...\n")

all_samples_meta <- run_table[, c("Run", "disease", "AGE", "sex", "BMI", "Sample Name", "isolation_source")]
all_samples_meta$Group <- ifelse(all_samples_meta$disease == "AD", "Case", "Control")
all_samples_meta$SampleType <- ifelse(all_samples_meta$isolation_source == "Oral_cavity", "Oral", "Fecal")

write.table(all_samples_meta, "~/Desktop/improved_microbiome_analysis/all_samples_metadata.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nFull metadata saved to ~/Desktop/improved_microbiome_analysis/all_samples_metadata.tsv\n")

cat("\n=== Next: Create sample tables for analysis ===\n")
cat("1. Oral-only analysis (all 99 samples)\n")
cat("2. Fecal-only analysis (99 samples)\n")
cat("3. Combined oral+fecal analysis (198 samples)\n")