# dada2_16S_pipeline_standalone.R

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(tidyverse)
  library(decontam)
  library(Biostrings)
})

# ----------- USER SETTINGS -----------
input_dir   <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/hostfree"
output_dir  <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2"
metadata_file    <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"
silva_train_file <- "/Volumes/edKwamiBackUP/analysis_workflow/silva_db/silva_nr99_v138.1_train_set.fa.gz"
silva_species_file <- "/Volumes/edKwamiBackUP/analysis_workflow/silva_db/silva_species_assignment_v138.1.fa.gz"
# -------------------------------------

# Auto-detect samples
r1_files_all <- list.files(input_dir, pattern="_R1\\.hostfree\\.fastq$", full.names=TRUE)
samples <- sub("_R1\\.hostfree\\.fastq$", "", basename(r1_files_all))

# Check for matching R2 files
r2_files_all <- file.path(input_dir, paste0(samples, "_R2.hostfree.fastq"))
if (!all(file.exists(r2_files_all))) {
  missing <- samples[!file.exists(r2_files_all)]
  stop(paste("Missing R2 files for samples:", paste(missing, collapse=", ")))
}

r1_files <- file.path(input_dir, paste0(samples, "_R1.hostfree.fastq"))
r2_files <- file.path(input_dir, paste0(samples, "_R2.hostfree.fastq"))
names(r1_files) <- samples; names(r2_files) <- samples

cat("Detected samples:\n")
print(samples)

# Filtering and trimming
filt_dir <- file.path(output_dir, "dada2_work")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
filtFs <- file.path(filt_dir, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(samples, "_R_filt.fastq.gz"))

out_f <- filterAndTrim(r1_files, filtFs, r2_files, filtRs,
  truncLen=c(0,0), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  compress=TRUE, multithread=TRUE)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplication, denoising, merging
derepFs <- derepFastq(filtFs); names(derepFs) <- samples
derepRs <- derepFastq(filtRs); names(derepRs) <- samples
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    minFoldParentOverAbundance=1.5,
                                    multithread=TRUE)

# ---------------- REMOVE SINGLETONS ----------------
seqtab.nochim <- seqtab.nochim[, colSums(seqtab.nochim) > 1]
cat("After singleton removal, number of ASVs:", ncol(seqtab.nochim), "\n")
# ---------------------------------------------------

# Taxonomic assignment
taxa <- assignTaxonomy(seqtab.nochim, silva_train_file, multithread=TRUE)
taxa <- addSpecies(taxa, silva_species_file)

# Load metadata
metadata <- read.table(metadata_file, header=TRUE, sep="\t", row.names=1)

# Build phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(taxa))

# ----------- EXPORTS -----------

# 1. ASV Table (samples as rows, ASVs as columns, SampleID as first column, ASV1, ASV2, ...)
asv_seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV", seq_len(length(asv_seqs)))
colnames(seqtab.nochim) <- asv_ids
asv_tab <- as.data.frame(seqtab.nochim)
asv_tab <- tibble::rownames_to_column(asv_tab, var="SampleID")
asv_table_file <- file.path(output_dir, "asv_table.tsv")
tryCatch({
  write.table(asv_tab, asv_table_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote", asv_table_file, "\n")
}, error=function(e) cat("Failed to write asv_table.tsv:", e$message, "\n"))

# 2. ASV Taxonomy Table (ASVs as rows, ASV1, ASV2, ...)
asv_tax_file <- file.path(output_dir, "asv_taxonomy.tsv")
asv_tax <- as.data.frame(taxa)
rownames(asv_tax) <- asv_ids
asv_tax <- tibble::rownames_to_column(asv_tax, var="ASV")
tryCatch({
  write.table(asv_tax, asv_tax_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote", asv_tax_file, "\n")
}, error=function(e) cat("Failed to write asv_taxonomy.tsv:", e$message, "\n"))

# 3. Read Tracking Table
getN <- function(x) sum(getUniques(x))
countFastq <- function(f) {
  con <- if(grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  n <- 0L
  while(length(chunk <- readLines(con, n=4*10000)) > 0) {
    n <- n + length(chunk)/4
  }
  close(con)
  return(as.integer(n))
}
read_tracking <- data.frame(
  Sample = samples,
  input = sapply(r1_files, countFastq),
  filtered = out_f[,1],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)
read_tracking_file <- file.path(output_dir, "read_tracking.tsv")
tryCatch({
  write.table(read_tracking, read_tracking_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote", read_tracking_file, "\n")
}, error=function(e) cat("Failed to write read_tracking.tsv:", e$message, "\n"))

# 4. Representative Sequences FASTA (ASV1, ASV2, ...)
rep_seqs_file <- file.path(output_dir, "rep_seqs.fasta")
names(asv_seqs) <- asv_ids
repseqs <- DNAStringSet(asv_seqs)
names(repseqs) <- asv_ids
tryCatch({
  writeXStringSet(repseqs, rep_seqs_file, format="fasta")
  cat("Wrote", rep_seqs_file, "\n")
}, error=function(e) cat("Failed to write rep_seqs.fasta:", e$message, "\n"))

# Save R objects
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim_singletonfree.rds"))
saveRDS(taxa, file.path(output_dir, "taxa_singletonfree.rds"))
saveRDS(ps, file.path(output_dir, "phyloseq_singletonfree.rds"))

cat("Pipeline complete. Results saved to", output_dir, "\n")
