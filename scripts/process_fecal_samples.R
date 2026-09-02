#!/usr/bin/env Rscript
# Process fecal samples through DADA2 pipeline

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(ggplot2)
})

cat("=== Fecal Sample Processing ===\n")

run_table_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv"
run_table <- read.csv(run_table_path)

fecal_samples <- run_table[run_table$isolation_source == "feces", ]
fecal_srr <- fecal_samples$Run

cat("Fecal samples to process:", length(fecal_srr), "\n")
cat("AD cases:", sum(fecal_samples$disease == "AD"), "\n")
cat("Controls:", sum(fecal_samples$disease == "Control"), "\n")

raw_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/raw/"
fecal_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_raw"
dir.create(fecal_dir, recursive = TRUE)

existing_files <- list.files(raw_path, pattern = "_R1.fastq$")
existing_srr <- gsub("_R1.fastq", "", existing_files)

fecal_to_download <- fecal_srr[!fecal_srr %in% existing_srr]
cat("\nFecal samples already downloaded:", length(fecal_srr) %in% existing_srr), "\n")
cat("Need to download:", length(fecal_to_download), "\n")

fecal_files_r1 <- list.files(raw_path, pattern = "_R1.fastq$", full.names = TRUE)
fecal_files_r2 <- list.files(raw_path, pattern = "_R2.fastq$", full.names = TRUE)

fecal_srr_in_raw <- unique(gsub("_R1.fastq|_R2.fastq", "", 
                                 c(basename(fecal_files_r1), basename(fecal_files_r2))))

fecal_in_data <- fecal_srr[fecal_srr %in% fecal_srr_in_raw]
cat("\nFecal samples available in raw folder:", length(fecal_in_data), "\n")

fecal_r1 <- file.path(raw_path, paste0(fecal_in_data, "_R1.fastq"))
fecal_r2 <- file.path(raw_path, paste0(fecal_in_data, "_R2.fastq"))

cat("\nProcessing", length(fecal_r1), "fecal samples through DADA2...\n")

fecal_r1 <- fecal_r1[file.exists(fecal_r1)]
fecal_r2 <- fecal_r2[file.exists(fecal_r2)]
names(fecal_r1) <- gsub("_R1.fastq", "", basename(fecal_r1))
names(fecal_r2) <- gsub("_R2.fastq", "", basename(fecal_r2))

cat("Files to process:", length(fecal_r1), "\n")

out_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2_fecal"
dir.create(out_dir, recursive = TRUE)

cat("\n=== Quality filtering ===\n")
filt_r1 <- file.path(out_dir, paste0(names(fecal_r1), "_R1_filt.fastq"))
filt_r2 <- file.path(out_dir, paste0(names(fecal_r2), "_R2_filt.fastq"))

out <- filterAndTrim(fecal_r1, filt_r1, fecal_r2, filt_r2,
                     truncLen = c(240, 200), minEE = 2, maxN = 0,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

cat("Reads after filtering:\n")
print(head(out))

cat("\n=== Learning error rates ===\n")
err_r1 <- learnErrors(filt_r1, multithread = TRUE)
err_r2 <- learnErrors(filt_r2, multithread = TRUE)

cat("\n=== Dereplication ===\n")
derep_r1 <- derepFastq(filt_r1, verbose = TRUE)
derep_r2 <- derepFastq(filt_r2, verbose = TRUE)

cat("\n=== Sample inference ===\n")
dada_r1 <- dada(derep_r1, err = err_r1, pool = "pseudo", multithread = TRUE)
dada_r2 <- dada(derep_r2, err = err_r2, pool = "pseudo", multithread = TRUE)

cat("\n=== Merging pairs ===\n")
mergers <- mergePairs(dada_r1, derep_r1, dada_r2, derep_r2, verbose = TRUE)

cat("\n=== Creating sequence table ===\n")
seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")

cat("\n=== Removing chimeras ===\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
cat("After chimera removal:", dim(seqtab_nochim), "\n")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_r1, getN), sapply(dada_r2, getN), 
               sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- names(fecal_r1)

cat("\n=== Read tracking ===\n")
print(head(track))

write.table(track, file.path(out_dir, "read_tracking.tsv"), sep = "\t", quote = FALSE)

cat("\n=== Assigning taxonomy ===\n")
taxa <- assignTaxonomy(seqtab_nochim, 
                       "/Volumes/edKwamiBackUP/analysis_workflow/silva_db/silva_train_set.fa.gz",
                       multithread = TRUE)

cat("Taxonomy dimensions:", dim(taxa), "\n")

asv_seq <- colnames(seqtab_nochim)
asv_headers <- paste0("ASV", seq_along(asv_seq))
colnames(seqtab_nochim) <- asv_headers

asv_df <- t(seqtab_nochim)
asv_df <- as.data.frame(asv_df)

taxa_df <- as.data.frame(taxa)
taxa_df$ASV <- rownames(taxa_df)

write.table(asv_df, file.path(out_dir, "asv_table.tsv"), sep = "\t", quote = FALSE)
write.table(taxa_df, file.path(out_dir, "asv_taxonomy.tsv"), sep = "\t", quote = FALSE)

cat("\n=== Saving outputs ===\n")
saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_fecal.rds"))
saveRDS(taxa, file.path(out_dir, "taxa_fecal.rds"))

cat("Total ASVs:", ncol(seqtab_nochim), "\n")
cat("Total samples:", nrow(seqtab_nochim), "\n")

cat("\nFecal sample processing complete!\n")
cat("Output directory:", out_dir, "\n")