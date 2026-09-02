#!/usr/bin/env Rscript
library(dada2)
cat("=== Loading sequence table ===\n")
seqtab <- readRDS("/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/seqtab.rds")
cat("Loaded seqtab with dimensions:", paste(dim(seqtab), collapse="x"), "\n")

cat("=== Assigning taxonomy (this may take a few minutes) ===\n")
taxa <- assignTaxonomy(seqtab, "/Volumes/edKwamiBackUP/analysis_workflow/silva_db/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

cat("=== Saving results ===\n")
saveRDS(taxa, "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/taxa.rds")
write.csv(taxa, "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/taxonomy.csv")
cat("Taxonomy assignment complete!\n")
