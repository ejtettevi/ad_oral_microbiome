#!/usr/bin/env Rscript
library(Biostrings)
library(tidyverse)

cat("=== HOMD Priority Taxonomy Assignment (Fast) ===\n\n")

# Load ASV sequences from fecal data
fecal_seqtab <- readRDS('/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/seqtab.rds')
cat("Fecal ASVs:", ncol(fecal_seqtab), "\n")

asv_seqs <- colnames(fecal_seqtab)
asv_dna <- DNAStringSet(asv_seqs)

# Load HOMD
homd_path <- '/Users/edkwami/Desktop/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta'
homd_seqs <- readDNAStringSet(homd_path)
cat("HOMD reference:", length(homd_seqs), "\n")

# Parse HOMD taxonomy
parse_homd <- function(header) {
  parts <- strsplit(header, ' ')[[1]]
  genus <- parts[3]
  species <- ifelse(length(parts) >= 4 && parts[4] != 'sp.', paste(parts[3], parts[4], sep=' '), NA)
  list(genus=genus, species=species)
}

homd_tax <- data.frame(
  seq_id = names(homd_seqs),
  genus = sapply(names(homd_seqs), parse_homd)$genus,
  species = sapply(names(homd_seqs), parse_homd)$species
)

# Convert to character for faster matching
homd_seq_char <- as.character(homd_seqs)

# Exact matching
cat("\n=== Exact matching ===\n")
exact_matches <- match(asv_seqs, homd_seq_char)
exact_count <- sum(!is.na(exact_matches))
cat("Exact matches:", exact_count, "/", length(asv_seqs), "\n")

# For non-exact, use 100% identity via shorter substring matching (first 100bp)
cat("\n=== Partial matching (first 100bp) ===\n")
asv_prefix <- substr(asv_seqs, 1, 100)
homd_prefix <- substr(homd_seq_char, 1, 100)

prefix_matches <- match(asv_prefix, homd_prefix)
partial_count <- sum(!is.na(prefix_matches) & is.na(exact_matches))
cat("Partial matches (100bp):", partial_count, "\n")

# Use 200bp matching
cat("\n=== Extended matching (200bp) ===\n")
asv_200 <- substr(asv_seqs, 1, 200)
homd_200 <- substr(homd_seq_char, 1, 200)

matches_200 <- match(asv_200, homd_200)
extended_count <- sum(!is.na(matches_200) & is.na(exact_matches))
cat("Extended matches (200bp):", extended_count, "\n")

# Combine results
results <- data.frame(
  ASV = paste0("ASV_", 1:length(asv_seqs)),
  sequence = asv_seqs,
  genus = NA,
  species = NA,
  match_type = NA
)

# Apply exact matches
exact_idx <- which(!is.na(exact_matches))
if(length(exact_idx) > 0) {
  results$genus[exact_idx] <- homd_tax$genus[exact_matches[exact_idx]]
  results$species[exact_idx] <- homd_tax$species[exact_matches[exact_idx]]
  results$match_type[exact_idx] <- "exact"
}

# Apply partial (100bp) matches for remaining
partial_idx <- which(is.na(results$genus) & !is.na(prefix_matches))
if(length(partial_idx) > 0) {
  results$genus[partial_idx] <- homd_tax$genus[prefix_matches[partial_idx]]
  results$species[partial_idx] <- homd_tax$species[prefix_matches[partial_idx]]
  results$match_type[partial_idx] <- "partial_100bp"
}

# Apply extended (200bp) matches for remaining
extended_idx <- which(is.na(results$genus) & !is.na(matches_200))
if(length(extended_idx) > 0) {
  results$genus[extended_idx] <- homd_tax$genus[matches_200[extended_idx]]
  results$species[extended_idx] <- homd_tax$species[matches_200[extended_idx]]
  results$match_type[extended_idx] <- "partial_200bp"
}

cat("\n=== Results Summary ===\n")
cat("Total ASVs:", nrow(results), "\n")
cat("Assigned genus:", sum(!is.na(results$genus)), "\n")
cat("Assigned species:", sum(!is.na(results$species)), "\n")
cat("\nMatch types:\n")
print(table(results$match_type))

# Save
write.csv(results, '/Users/edkwami/Desktop/improved_microbiome_analysis/homd_species_assignment.csv', row.names = FALSE)

cat("\nTop genera:\n")
print(head(sort(table(results$genus[!is.na(results$genus)]), decreasing = TRUE), 15))

cat("\nTop species:\n")
print(head(sort(table(results$species[!is.na(results$species)]), decreasing = TRUE), 15))

cat("\nDone! Results saved.\n")