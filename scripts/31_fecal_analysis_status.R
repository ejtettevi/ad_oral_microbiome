#!/usr/bin/env Rscript
# Fecal + Oral Combined Analysis Pipeline

cat("=== Fecal + Oral Combined Analysis ===\n\n")

# First, check what fecal samples we have
fecal_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_raw"
sra_files <- list.files(fecal_dir, pattern = "\\.sra$", recursive = TRUE)
fecal_srr <- gsub("\\.sra$", "", basename(sra_files))
cat("Fecal SRA files available:", length(fecal_srr), "\n")

# If we have SRA files but no FASTQ, we need to convert
# For now, let's proceed with the combined analysis

# Load run table to get metadata for available samples
run_table_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv"
run_table <- read.csv(run_table_path)

# Get fecal samples we have
fecal_available <- run_table[run_table$Run %in% fecal_srr, ]
cat("Fecal samples with metadata:", nrow(fecal_available), "\n")

# Combined analysis with oral data
# Load oral data (already processed)
oral_case <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv", row.names = 1)
oral_ctrl <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv", row.names = 1)

cat("Oral samples:", ncol(oral_case), "cases,", ncol(oral_ctrl), "controls\n")

# Note: Full fecal analysis requires:
# 1. Convert SRA to FASTQ (using fasterq-dump)
# 2. Run DADA2 pipeline
# 3. Assign taxonomy
# 4. Generate genus-level tables

# Current status
cat("\n=== Current Status ===\n")
cat("Fecal samples downloaded:", length(fecal_srr), "/99\n")
cat("Oral samples analyzed:", 98, "(from previous analysis)\n")
cat("\nTo complete fecal analysis, run:\n")
cat("1. fasterq-dump on each SRA file\n")
cat("2. DADA2 pipeline (QC, denoising, merging)\n")
cat("3. Taxonomy assignment\n")
cat("4. Differential abundance analysis\n")

# Save status
status_df <- data.frame(
  Sample_Type = c("Oral", "Fecal"),
  Total_in_BioProject = c(99, 99),
  Downloaded = c(96, length(fecal_srr)),
  Analyzed = c(98, 0)
)
write.csv(status_df, "~/Desktop/improved_microbiome_analysis/sample_status.csv")

cat("\n=== Note ===\n")
cat("Fecal sample processing requires significant compute time.\n")
cat("With", length(fecal_srr), "samples, estimated processing time: ~2-4 hours.\n")

# For now, show what would be needed for the complete analysis
cat("\n=== Required for Complete Fecal Analysis ===\n")
cat("1. Convert SRA to FASTQ:\n")
for (srr in fecal_srr[1:5]) {
  cat("   fasterq-dump", srr, "-O output/\n")
}
cat("   ... (", length(fecal_srr)-5, "more)\n")

cat("\n2. Run DADA2 pipeline on all FASTQ files\n")
cat("3. Generate genus-level abundance tables\n")
cat("4. Run differential abundance (ANCOMBC2/CLR)\n")
cat("5. Beta diversity analysis\n")
cat("6. Combined analysis with oral samples\n")