#!/usr/bin/env Rscript
library(Biostrings)
library(stringr)

cat("=== HOMD Species-Level Taxonomy Assignment ===\n\n")

# Paths
homd_path <- '/Users/edkwami/Desktop/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta'
output_dir <- '/Users/edkwami/Desktop/improved_microbiome_analysis'

# Load ASV sequences from fecal data
fecal_seqtab <- readRDS('/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/seqtab.rds')
cat("Fecal ASVs:", ncol(fecal_seqtab), "\n")

asv_seqs <- colnames(fecal_seqtab)
asv_dna <- DNAStringSet(asv_seqs)
cat("ASV sequences loaded\n")

# Load HOMD reference
homd_seqs <- readDNAStringSet(homd_path)
cat("HOMD reference:", length(homd_seqs), "\n")

# Parse HOMD taxonomy
genus_vec <- sapply(1:length(homd_seqs), function(i) {
  p <- strsplit(names(homd_seqs)[i], ' ')[[1]]
  if(length(p) >= 3) p[3] else NA
})

species_vec <- sapply(1:length(homd_seqs), function(i) {
  p <- strsplit(names(homd_seqs)[i], ' ')[[1]]
  if(length(p) >= 4 && p[4] != 'sp.') paste(p[3], p[4], sep=' ') else NA
})

homd_char <- as.character(homd_seqs)
cat("HOMD parsed: ", length(unique(genus_vec[!is.na(genus_vec)])), " genera, ", sum(!is.na(species_vec)), " species\n", sep="")

# Matching
cat("\n=== Matching ASVs to HOMD ===\n")

# Exact match
exact_match <- match(asv_seqs, homd_char)
exact_count <- sum(!is.na(exact_match))
cat("Exact matches: ", exact_count, "/", length(asv_seqs), "\n", sep="")

# 100bp prefix match for remaining
asv_100 <- substr(asv_seqs, 1, 100)
homd_100 <- substr(homd_char, 1, 100)
prefix_match <- match(asv_100, homd_100)
prefix_count <- sum(!is.na(prefix_match) & is.na(exact_match))
cat("100bp prefix matches: ", prefix_count, "\n", sep="")

# Build results
results <- data.frame(
  ASV = paste0("ASV_", 1:length(asv_seqs)),
  sequence = asv_seqs,
  genus = NA_character_,
  species = NA_character_,
  match_method = NA_character_
)

# Apply exact matches
exact_idx <- which(!is.na(exact_match))
results$genus[exact_idx] <- genus_vec[exact_match[exact_idx]]
results$species[exact_idx] <- species_vec[exact_match[exact_idx]]
results$match_method[exact_idx] <- "exact"

# Apply prefix matches
remaining <- which(is.na(results$genus) & !is.na(prefix_match))
if(length(remaining) > 0) {
  results$genus[remaining] <- genus_vec[prefix_match[remaining]]
  results$species[remaining] <- species_vec[prefix_match[remaining]]
  results$match_method[remaining] <- "100bp_prefix"
}

cat("\n=== Results Summary ===\n")
cat("Total ASVs: ", nrow(results), "\n")
cat("Assigned genus: ", sum(!is.na(results$genus)), " (", round(100*sum(!is.na(results$genus))/nrow(results),1), "%)\n", sep="")
cat("Assigned species: ", sum(!is.na(results$species)), " (", round(100*sum(!is.na(results$species))/nrow(results),1), "%)\n", sep="")

cat("\nMatch method distribution:\n")
print(table(results$match_method, useNA="ifany"))

cat("\nTop genera:\n")
print(head(sort(table(results$genus[!is.na(results$genus)]), decreasing=TRUE), 15))

cat("\nTop species:\n")
print(head(sort(table(results$species[!is.na(results$species)]), decreasing=TRUE), 15))

# Save
write.csv(results, file.path(output_dir, 'homd_species_assignment.csv'), row.names=FALSE)
cat("\nSaved to: homd_species_assignment.csv\n")

# Create combined taxonomy table (for oral data too)
# Need to also add SILVA assignment for non-HOMD matches
cat("\n=== Creating Full Taxonomy Table ===\n")

# For oral ASVs - need to check
if(file.exists('/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv')) {
  cat("Note: Oral ASV sequences not available in this format.\n")
  cat("HOMD species assignment is available for fecal ASVs.\n")
}

cat("\nDone!\n")