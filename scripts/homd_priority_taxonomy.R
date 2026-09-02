#!/usr/bin/env Rscript
library(Biostrings)
library(tidyverse)
library(ShortRead)

cat("=== HOMD Priority Taxonomy Assignment ===\n\n")

# Load ASV sequences
fecal_seqtab <- readRDS('/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed/seqtab.rds')
cat("Fecal ASVs:", ncol(fecal_seqtab), "\n")

# Get unique ASV sequences (column names from seqtab are sequences)
asv_seqs <- colnames(fecal_seqtab)
cat("ASV sequences to classify:", length(asv_seqs), "\n")

# Convert to DNAStringSet
asv_dna <- DNAStringSet(asv_seqs)

# Load HOMD reference sequences
homd_path <- '/Users/edkwami/Desktop/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta'
homd_seqs <- readDNAStringSet(homd_path)
cat("HOMD reference sequences:", length(homd_seqs), "\n")

# Parse HOMD taxonomy from headers
parse_homd <- function(header) {
  parts <- strsplit(header, ' ')[[1]]
  genus <- parts[3]
  species <- ifelse(length(parts) >= 4 && parts[4] != 'sp.', paste(parts[3], parts[4], sep=' '), NA)
  list(genus=genus, species=species)
}

homd_tax <- data.frame(
  seq_id = names(homd_seqs),
  genus = sapply(names(homd_seqs), function(x) parse_homd(x)$genus),
  species = sapply(names(homd_seqs), function(x) parse_homd(x)$species)
)

# For each ASV, find the best HOMD match using exact matching or k-mer similarity
# Using pairwise alignment for better matching
cat("\n=== Performing sequence alignment ===\n")

# Function to find best match
find_best_homd_match <- function(asv_seq, homd_seqs, homd_tax, max_matches = 5) {
  # Try exact match first
  exact_match <- which(as.character(homd_seqs) == asv_seq)
  if(length(exact_match) > 0) {
    return(list(
      genus = homd_tax$genus[exact_match[1]],
      species = homd_tax$species[exact_match[1]],
      identity = 1.0
    ))
  }

  # Use alignment for approximate matching
  aln <- pairwiseAlignment(asv_seq, homd_seqs, type="global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1))
  best_idx <- which.max(score(aln))
  best_score <- score(aln)
  seq_len <- width(asv_seq)

  # Calculate identity
  identity <- best_score / (2 * seq_len)

  if(identity > 0.85) {
    return(list(
      genus = homd_tax$genus[best_idx],
      species = homd_tax$species[best_idx],
      identity = identity
    ))
  }

  return(list(genus = NA, species = NA, identity = identity))
}

# Process in batches
cat("Processing", length(asv_seqs), "ASVs...\n")

results <- list()
for(i in 1:length(asv_seqs)) {
  if(i %% 100 == 0) cat("Progress:", i, "/", length(asv_seqs), "\n")

  asv_seq <- asv_dna[i]
  match <- find_best_homd_match(asv_seq, homd_seqs, homd_tax)

  results[[i]] <- data.frame(
    ASV = paste0("ASV_", i),
    sequence = asv_seqs[i],
    genus = match$genus,
    species = match$species,
    identity = match$identity
  )
}

tax_results <- do.call(rbind, results)
cat("\nAssigned genus:", sum(!is.na(tax_results$genus)), "/", nrow(tax_results), "\n")
cat("Assigned species:", sum(!is.na(tax_results$species)), "/", nrow(tax_results), "\n")

# Save results
write.csv(tax_results, '/Users/edkwami/Desktop/improved_microbiome_analysis/homd_species_assignment.csv', row.names = FALSE)
cat("\nResults saved to homd_species_assignment.csv\n")

# Print summary
cat("\nTop genera found:\n")
print(head(sort(table(tax_results$genus), decreasing = TRUE), 10))

cat("\nTop species found:\n")
print(head(sort(table(tax_results$species[!is.na(tax_results$species)]), decreasing = TRUE), 10))