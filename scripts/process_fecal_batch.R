#!/usr/bin/env Rscript
# Process fecal samples through DADA2 - batch version

suppressPackageStartupMessages({
  library(dada2)
})

cat("=== Starting Fecal Sample Processing ===\n")

run_table_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv"
run_table <- read.csv(run_table_path)

fecal_samples <- run_table[run_table$isolation_source == "feces", ]
fecal_srr <- fecal_samples$Run
cat("Fecal samples in BioProject:", length(fecal_srr), "\n")

raw_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/raw/"
all_files <- list.files(raw_path, pattern = "\\.fastq$")
all_srr <- unique(gsub("_R[12]\\.fastq", "", all_files))

fecal_available <- fecal_srr[fecal_srr %in% all_srr]
cat("Fecal samples available in raw folder:", length(fecal_available), "\n")

if (length(fecal_available) == 0) {
  cat("\nNo fecal samples found in raw folder!\n")
  cat("Fecal SRR IDs in BioProject:\n")
  print(head(fecal_srr))
  cat("\nAvailable SRR IDs in raw:\n")
  print(head(all_srr))
  stop("Cannot proceed without fecal FASTQ files")
}

cat("\nProcessing", length(fecal_available), "fecal samples...\n")

fecal_r1 <- sort(file.path(raw_path, paste0(fecal_available, "_R1.fastq")))
fecal_r2 <- sort(file.path(raw_path, paste0(fecal_available, "_R2.fastq")))

fecal_r1 <- fecal_r1[file.exists(fecal_r1)]
fecal_r2 <- fecal_r2[file.exists(fecal_r2)]

names(fecal_r1) <- gsub("_R1.fastq", "", basename(fecal_r1))
names(fecal_r2) <- gsub("_R2.fastq", "", basename(fecal_r2))

out_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2_fecal"
dir.create(out_dir, recursive = TRUE)

cat("\n1. Quality filtering...\n")
filt_r1 <- file.path(out_dir, paste0(names(fecal_r1), "_R1_filt.fastq.gz"))
filt_r2 <- file.path(out_dir, paste0(names(fecal_r2), "_R2_filt.fastq.gz"))

out <- filterAndTrim(fecal_r1, filt_r1, fecal_r2, filt_r2,
                     truncLen = c(240, 200), minEE = 2, maxN = 0,
                     rm.phix = TRUE, compress = TRUE, multithread = 4)

cat("\nReads remaining after filtering:\n")
print(head(out))
cat("Mean reads per sample:", round(mean(out[,2]), 0), "\n")

cat("\n2. Learning error rates...\n")
err_r1 <- learnErrors(filt_r1, multithread = 4)
err_r2 <- learnErrors(filt_r2, multithread = 4)

cat("\n3. Dereplication...\n")
derep_r1 <- derepFastq(filt_r1, verbose = FALSE)
derep_r2 <- derepFastq(filt_r2, verbose = FALSE)

cat("\n4. Sample inference (DADA2)...\n")
dada_r1 <- dada(derep_r1, err = err_r1, pool = "pseudo", multithread = 4)
dada_r2 <- dada(derep_r2, err = err_r2, pool = "pseudo", multithread = 4)

cat("\n5. Merging paired reads...\n")
mergers <- mergePairs(dada_r1, derep_r1, dada_r2, derep_r2, verbose = FALSE)

cat("\n6. Creating sequence table...\n")
seqtab <- makeSequenceTable(mergers)
cat("Dimensions:", dim(seqtab), "\n")

cat("\n7. Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 4)
cat("After chimera removal:", dim(seqtab_nochim), "\n")

cat("\n8. Assigning taxonomy...\n")
taxa <- assignTaxonomy(seqtab_nochim, 
                       "/Volumes/edKwamiBackUP/analysis_workflow/silva_db/silva_train_set.fa.gz",
                       multithread = 4)

cat("Taxonomy dimensions:", dim(taxa), "\n")

cat("\n9. Saving outputs...\n")
asv_seq <- colnames(seqtab_nochim)
asv_headers <- paste0("ASV", seq_along(asv_seq))
colnames(seqtab_nochim) <- asv_headers

asv_df <- t(seqtab_nochim)
asv_df <- as.data.frame(asv_df)
asv_df <- cbind(SampleID = rownames(asv_df), asv_df)
rownames(asv_df) <- NULL

write.table(asv_df, file.path(out_dir, "asv_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

taxa_df <- as.data.frame(taxa)
taxa_df$ASV <- rownames(taxa_df)
rownames(taxa_df) <- NULL
write.table(taxa_df, file.path(out_dir, "asv_taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_r1, getN), sapply(dada_r2, getN), 
               sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- names(fecal_r1)
write.table(track, file.path(out_dir, "read_tracking.tsv"), sep = "\t", quote = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Total samples processed:", nrow(seqtab_nochim), "\n")
cat("Total ASVs:", ncol(seqtab_nochim), "\n")
cat("Total reads (non-chimeric):", sum(seqtab_nochim), "\n")
cat("\nOutput files saved to:", out_dir, "\n")