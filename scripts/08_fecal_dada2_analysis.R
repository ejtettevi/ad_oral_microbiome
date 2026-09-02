#!/usr/bin/env Rscript
library(dada2)
library(phyloseq)
library(tidyverse)

cat("=== Fecal Sample DADA2 Analysis ===\n\n")

fecal_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_raw"
output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed"
dir.create(output_dir, showWarnings = FALSE)

sra_files <- list.files(fecal_dir, pattern = "\\.sra$", recursive = TRUE)
srr_names <- gsub("\\.sra$", "", basename(sra_files))
cat("SRA files found:", length(srr_names), "\n")

cat("\n=== Converting SRA to FASTQ ===\n")
fastq_dir <- file.path(output_dir, "fastq")
dir.create(fastq_dir, showWarnings = FALSE)

for(srr in srr_names[1:min(10, length(srr_names))]) {
  cat("Processing:", srr, "\n")
  sra_path <- file.path(fecal_dir, srr, paste0(srr, ".sra"))
  if(file.exists(sra_path)) {
    out <- tryCatch({
      system2("fasterq-dump", args = c(sra_path, "-O", fastq_dir, "-e", "4"), stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
  }
}

cat("\n=== DADA2 Pipeline ===\n")
fnFs <- sort(list.files(fastq_dir, pattern = "_1.fastq$", full.names = TRUE))
fnRs <- sort(list.files(fastq_dir, pattern = "_2.fastq$", full.names = TRUE))

if(length(fnFs) == 0) {
  cat("No FASTQ files found. Using alternative SRA extraction.\n")
  quit(status = 0)
}

cat("Forward reads:", length(fnFs), "\n")
cat("Reverse reads:", length(fnRs), "\n")

sample_names <- sapply(basename(fnFs), function(x) gsub("_1\\.fastq$", "", x))

plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

filtFs <- file.path(fastq_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(fastq_dir, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(240, 160),
                     maxEE = c(2, 5),
                     trimLeft = 10,
                     compress = TRUE,
                     multithread = TRUE)

cat("\n=== Learning Error Rates ===\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

cat("\n=== Sample Inference ===\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

cat("\n=== Merging Reads ===\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, multithread = TRUE)
cat("After chimera removal:", dim(seqtab.nochim), "\n")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names

cat("\n=== Taxonomy Assignment ===\n")
taxa <- assignTaxonomy(seqtab.nochim, "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/silva_files/silva_train_set_138_1.gz", multithread = TRUE)

cat("\n=== Saving Results ===\n")
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab.rds"))
saveRDS(taxa, file.path(output_dir, "taxa.rds"))
write.csv(track, file.path(output_dir, "tracking.csv"))
write.csv(as.matrix(seqtab.nochim), file.path(output_dir, "asv_table.csv"))
write.csv(taxa, file.path(output_dir, "taxonomy.csv"))

cat("\n=== Complete ===\n")
cat("ASVs:", ncol(seqtab.nochim), "\n")
cat("Samples:", nrow(seqtab.nochim), "\n")