suppressMessages({
  library(phyloseq)
  library(ANCOMBC)
})

# ---- 1. Load phylum count table (samples x phylum) ----
ptab <- read.delim("/var/folders/45/7jcqmvc92yb4gl2q9jjw6d7h0000gn/T/opencode/phylum_table.tsv",
                   header=TRUE, check.names=FALSE)
rownames(ptab) <- ptab$SampleID
ptab$SampleID <- NULL
counts <- as.matrix(ptab)
counts <- t(counts)  # phyla as rows, samples as cols
counts <- round(counts)
counts[is.na(counts)] <- 0

# ---- 2. Load metadata ----
meta <- read.delim("/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv",
                   header=TRUE, check.names=FALSE)
rownames(meta) <- meta$Sample
meta$Sample <- NULL
meta$Group <- factor(meta$Group, levels=c("Control","Case"))  # control = reference

# align
common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common]
meta <- meta[common, , drop=FALSE]
cat("n samples:", ncol(counts), "\n")
cat("groups:", table(meta$Group), "\n")
cat("n phyla:", nrow(counts), "\n")

# ---- 3. Build phyloseq object ----
OTU <- otu_table(counts, taxa_are_rows=TRUE)
SAM <- sample_data(meta)
ps <- phyloseq(OTU, SAM)

# ---- 4. Run ANCOM-BC2 ----
out <- ancombc2(
  data = ps,
  fix_formula = "Group",
  group = "Group",
  p_adj_method = "BH",
  alpha = 0.05,
  pseudo = 0,
  prv_cut = 0.001,
  lib_cut = 0,
  verbose = TRUE
)

res <- out$res
print(res)

# Save full output
write.csv(res, "/var/folders/45/7jcqmvc92yb4gl2q9jjw6d7h0000gn/T/opencode/phylum_ancombc2_RAW.csv", row.names=FALSE)
cat("\nSaved raw results.\n")
