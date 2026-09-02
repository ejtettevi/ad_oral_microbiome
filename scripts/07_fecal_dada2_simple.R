#!/usr/bin/env Rscript
library(dada2)
library(phyloseq)
library(tidyverse)

cat("=== Fecal Sample DADA2 Analysis ===\n\n")

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed"
fastq_dir <- file.path(output_dir, "fastq")

fnFs <- sort(list.files(fastq_dir, pattern = "_1.fastq$", full.names = TRUE))
fnRs <- sort(list.files(fastq_dir, pattern = "_2.fastq$", full.names = TRUE))

cat("Forward reads:", length(fnFs), "\n")
cat("Reverse reads:", length(fnRs), "\n")

sample_names <- sapply(basename(fnFs), function(x) gsub("_1\\.fastq$", "", x))

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
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = FALSE)

seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, multithread = TRUE)
cat("After chimera removal:", dim(seqtab.nochim), "\n")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names

cat("\n=== Saving Results ===\n")
cat("ASVs:", ncol(seqtab.nochim), "\n")
cat("Samples:", nrow(seqtab.nochim), "\n")

saveRDS(seqtab.nochim, file.path(output_dir, "seqtab.rds"))
write.csv(track, file.path(output_dir, "tracking.csv"))
write.csv(as.matrix(seqtab.nochim), file.path(output_dir, "asv_table.csv"))

cat("\nResults saved to:", output_dir, "\n")