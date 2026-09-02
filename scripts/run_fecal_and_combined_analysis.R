#!/usr/bin/env Rscript
# Complete fecal and combined oral + fecal analysis script

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(randomForest)
})

cat("=== Starting Fecal-Only and Combined Analysis ===\n\n")

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed"
desktop_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(desktop_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load Fecal Data
# ============================================================================
cat("1. Loading fecal data...\n")
seqtab <- readRDS(file.path(output_dir, "seqtab.rds"))
taxa <- readRDS(file.path(output_dir, "taxa.rds"))

# Clean rownames of seqtab to get SRR IDs
rownames(seqtab) <- gsub("_F_filt\\.fastq\\.gz$", "", rownames(seqtab))

cat("Original Fecal table size: ", nrow(seqtab), " samples, ", ncol(seqtab), " ASVs\n")

# Load SRA Run Table
run_table <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv")

# CRITICAL FIX: Only include runs that are actually fecal samples in the SraRunTable!
fecal_meta <- run_table[run_table$Run %in% rownames(seqtab) & run_table$isolation_source == "feces", ]
rownames(fecal_meta) <- fecal_meta$Run

# Subset seqtab and taxa to keep only the actual fecal samples
seqtab <- seqtab[fecal_meta$Run, ]
cat("Filtered Fecal table size (actual feces only): ", nrow(seqtab), " samples, ", ncol(seqtab), " ASVs\n")

# Ensure metadata Group is factor and disease is aligned
fecal_meta$Group <- ifelse(fecal_meta$disease == "AD", "Case", "Control")
fecal_meta$Group <- factor(fecal_meta$Group, levels = c("Control", "Case"))

# Create fecal phyloseq object
fecal_phy <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(fecal_meta)
)

# ============================================================================
# 2. Fecal-Only Alpha Diversity
# ============================================================================
cat("\n2. Fecal-only alpha diversity...\n")
fecal_alpha <- estimate_richness(fecal_phy, measures = c("Observed", "Shannon", "Simpson"))
fecal_alpha$Sample <- rownames(fecal_alpha)
fecal_alpha$Group <- fecal_meta[fecal_alpha$Sample, "Group"]

# Wilcoxon test on Shannon, Observed, Simpson
shannon_p <- wilcox.test(Shannon ~ Group, data = fecal_alpha)$p.value
observed_p <- wilcox.test(Observed ~ Group, data = fecal_alpha)$p.value
simpson_p <- wilcox.test(Simpson ~ Group, data = fecal_alpha)$p.value

cat("Observed Richness: Control vs Case Wilcoxon p =", observed_p, "\n")
cat("Shannon Index: Control vs Case Wilcoxon p =", shannon_p, "\n")
cat("Simpson Index: Control vs Case Wilcoxon p =", simpson_p, "\n")

write.csv(fecal_alpha, file.path(desktop_dir, "fecal_alpha_diversity.csv"), row.names = FALSE)

# ============================================================================
# 3. Fecal-Only Genus-Level Beta Diversity (Aitchison Distance)
# ============================================================================
cat("\n3. Fecal-only genus-level beta diversity (Aitchison)...\n")
fecal_genus <- tax_glom(fecal_phy, taxrank = "Genus")

# Extract genus abundance matrix
fecal_gen_mat <- as.matrix(otu_table(fecal_genus))
# Map ASV names to Genus names
tax_info <- as.data.frame(tax_table(fecal_genus))
colnames(fecal_gen_mat) <- tax_info$Genus

# Handle any NA genus names (replace with unclassified)
colnames(fecal_gen_mat)[is.na(colnames(fecal_gen_mat))] <- "Unclassified"
# Sum up columns with the same genus name
fecal_gen_mat <- sapply(unique(colnames(fecal_gen_mat)), function(g) {
  rowSums(fecal_gen_mat[, colnames(fecal_gen_mat) == g, drop = FALSE])
})

# CLR transform
clr_transform <- function(x) {
  x[x == 0] <- 0.5
  log(x / exp(mean(log(x))))
}
fecal_clr <- t(apply(fecal_gen_mat, 1, clr_transform))

# PCA
pca <- prcomp(fecal_clr)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# PERMANOVA
dist_aitchison <- dist(fecal_clr, method = "euclidean")
permanova_fecal <- adonis2(dist_aitchison ~ Group, data = fecal_meta[rownames(fecal_clr), ], permutations = 999, by = "terms")
cat("Fecal Beta Diversity PERMANOVA:\n")
print(permanova_fecal)

pca_df <- data.frame(
  Sample = rownames(fecal_clr),
  Group = fecal_meta[rownames(fecal_clr), "Group"],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
write.csv(pca_df, file.path(desktop_dir, "fecal_beta_pca.csv"), row.names = FALSE)

# ============================================================================
# 4. Fecal-Only Differential Abundance
# ============================================================================
cat("\n4. Fecal-only differential abundance (Wilcoxon on CLR)...\n")
fecal_da <- data.frame(Genus = colnames(fecal_clr))
group_labels <- fecal_meta[rownames(fecal_clr), "Group"]

for (i in 1:nrow(fecal_da)) {
  g <- colnames(fecal_clr)[i]
  vals <- fecal_clr[, g]
  
  case_vals <- vals[group_labels == "Case"]
  ctrl_vals <- vals[group_labels == "Control"]
  
  w <- wilcox.test(case_vals, ctrl_vals)
  fecal_da$Mean_Case[i] <- mean(case_vals)
  fecal_da$Mean_Control[i] <- mean(ctrl_vals)
  fecal_da$Diff[i] <- mean(case_vals) - mean(ctrl_vals)
  fecal_da$p_value[i] <- w$p.value
}

fecal_da$q_value <- p.adjust(fecal_da$p_value, method = "fdr")
fecal_da <- fecal_da[order(fecal_da$p_value), ]
sig_fecal_da <- fecal_da[fecal_da$q_value < 0.1, ]
cat("Fecal significant genera (FDR < 0.1):", nrow(sig_fecal_da), "\n")
if (nrow(sig_fecal_da) > 0) {
  print(head(sig_fecal_da))
}

write.csv(fecal_da, file.path(desktop_dir, "fecal_differential_abundance.csv"), row.names = FALSE)

# ============================================================================
# 5. Combined Oral + Fecal Analysis
# ============================================================================
cat("\n5. Loading oral data for combined analysis...\n")
oral_case <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv")
oral_ctrl <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv")

# Set Taxon as rownames
rownames(oral_case) <- oral_case$Taxon
oral_case$Taxon <- NULL
rownames(oral_ctrl) <- oral_ctrl$Taxon
oral_ctrl$Taxon <- NULL

# Merge technical replicates
merge_replicates <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  base_names <- unique(gsub("_1$|_2$", "", colnames(df)))
  result <- matrix(0, nrow = nrow(df), ncol = length(base_names))
  rownames(result) <- rownames(df)
  colnames(result) <- base_names
  for (base in base_names) {
    cols <- grep(paste0("^", base, "_"), colnames(df))
    result[, base] <- rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  return(as.data.frame(result))
}

oral_case_merged <- merge_replicates(oral_case)
oral_ctrl_merged <- merge_replicates(oral_ctrl)

# Find common genera between cases and controls
common_oral_genera <- intersect(rownames(oral_case_merged), rownames(oral_ctrl_merged))
oral_case_aligned <- oral_case_merged[common_oral_genera, ]
oral_ctrl_aligned <- oral_ctrl_merged[common_oral_genera, ]

oral_mat <- cbind(oral_case_aligned, oral_ctrl_aligned)
oral_groups <- c(rep("Case", ncol(oral_case_aligned)), rep("Control", ncol(oral_ctrl_aligned)))

# Find common genera between oral and fecal datasets
common_genera <- intersect(colnames(fecal_gen_mat), rownames(oral_mat))
cat("Common genera between oral and fecal datasets:", length(common_genera), "\n")

# Combine datasets (rows = samples, columns = common genera)
fecal_comb_mat <- fecal_gen_mat[, common_genera]
oral_comb_mat <- t(oral_mat[common_genera, ])

combined_abundance <- rbind(fecal_comb_mat, oral_comb_mat)
combined_meta <- data.frame(
  Sample = rownames(combined_abundance),
  SampleType = c(rep("Fecal", nrow(fecal_comb_mat)), rep("Oral", nrow(oral_comb_mat))),
  Group = c(as.character(group_labels), oral_groups),
  stringsAsFactors = FALSE
)
rownames(combined_meta) <- combined_meta$Sample

cat("Combined Dataset: ", nrow(combined_abundance), " samples\n")
cat("Oral samples: ", sum(combined_meta$SampleType == "Oral"), "\n")
cat("Fecal samples: ", sum(combined_meta$SampleType == "Fecal"), "\n")

# CLR transform combined dataset
combined_clr <- t(apply(combined_abundance, 1, clr_transform))

# PCA on combined dataset
comb_pca <- prcomp(combined_clr)
comb_var_explained <- round(100 * comb_pca$sdev^2 / sum(comb_pca$sdev^2), 1)

# PERMANOVA on combined dataset with sequential term testing!
comb_dist <- dist(combined_clr, method = "euclidean")
permanova_comb <- adonis2(comb_dist ~ SampleType * Group, data = combined_meta, permutations = 999, by = "terms")
cat("\nCombined Dataset PERMANOVA:\n")
print(permanova_comb)

# Save combined results
comb_pca_df <- data.frame(
  Sample = rownames(combined_clr),
  SampleType = combined_meta$SampleType,
  Group = combined_meta$Group,
  PC1 = comb_pca$x[, 1],
  PC2 = comb_pca$x[, 2]
)
write.csv(comb_pca_df, file.path(desktop_dir, "combined_beta_pca.csv"), row.names = FALSE)

# ============================================================================
# 6. Generate Figures
# ============================================================================
cat("\n6. Generating figures...\n")

# Plot Fecal Alpha Diversity
png(file.path(desktop_dir, "fecal_alpha_diversity_plots.png"), width = 10, height = 4, units = "in", res = 700)
par(mfrow = c(1, 3))
boxplot(Observed ~ Group, data = fecal_alpha, main = paste("Observed Richness\nWilcox p =", round(observed_p, 4)), col = c("#4DBBD5", "#E64B35"), ylab = "Observed Genera")
boxplot(Shannon ~ Group, data = fecal_alpha, main = paste("Shannon Index\nWilcox p =", round(shannon_p, 4)), col = c("#4DBBD5", "#E64B35"), ylab = "Shannon Index")
boxplot(Simpson ~ Group, data = fecal_alpha, main = paste("Simpson Index\nWilcox p =", round(simpson_p, 4)), col = c("#4DBBD5", "#E64B35"), ylab = "Simpson Index")
dev.off()

# Plot Fecal Beta Diversity PCA
png(file.path(desktop_dir, "fecal_beta_pca_plot.png"), width = 6, height = 5, units = "in", res = 700)
plot(pca$x[, 1], pca$x[, 2], col = ifelse(group_labels == "Case", "#E64B35", "#4DBBD5"), 
     pch = 19, xlab = paste0("PC1 (", var_explained[1], "%)"), 
     ylab = paste0("PC2 (", var_explained[2], "%)"), 
     main = paste0("Fecal Beta Diversity (Aitchison)\nPERMANOVA p = ", round(permanova_fecal$`Pr(>F)`[1], 4)), cex.lab=1.2, cex.axis=1.1, cex.main=1.3)
legend("topright", legend = c("Control", "Case"), col = c("#4DBBD5", "#E64B35"), pch = 19, cex=1.1)
dev.off()

# Plot Combined PCA (Sample Type and Group)
png(file.path(desktop_dir, "combined_pca_plot.png"), width = 7, height = 6, units = "in", res = 700)
colors <- ifelse(combined_meta$SampleType == "Oral", 
                 ifelse(combined_meta$Group == "Case", "#E64B35", "#4DBBD5"), 
                 ifelse(combined_meta$Group == "Case", "#F39C12", "#2ECC71"))
pch_vals <- ifelse(combined_meta$SampleType == "Oral", 19, 17)

plot(comb_pca$x[, 1], comb_pca$x[, 2], col = colors, pch = pch_vals,
     xlab = paste0("PC1 (", comb_var_explained[1], "%)"), 
     ylab = paste0("PC2 (", comb_var_explained[2], "%)"), 
     main = "Combined PCA (Oral vs Fecal & Case vs Control)", cex.lab=1.2, cex.axis=1.1, cex.main=1.3)
legend("topright", legend = c("Oral Control", "Oral AD Case", "Fecal Control", "Fecal AD Case"), 
       col = c("#4DBBD5", "#E64B35", "#2ECC71", "#F39C12"), pch = c(19, 19, 17, 17), cex=1.1)
dev.off()

# Save a summary file
summary_stats <- paste0("
# Summary of Fecal and Combined Analysis

## 1. Fecal-Only Dataset
- Fecal samples processed: ", nrow(seqtab), "
- Cases (AD): ", sum(group_labels == "Case"), "
- Controls: ", sum(group_labels == "Control"), "
- Significant Genera (FDR < 0.1): ", sum(fecal_da$q_value < 0.1), "

### Alpha Diversity p-values:
- Observed Richness: p = ", round(observed_p, 4), "
- Shannon Index: p = ", round(shannon_p, 4), "
- Simpson Index: p = ", round(simpson_p, 4), "

### Beta Diversity:
- PERMANOVA p-value: p = ", round(permanova_fecal$`Pr(>F)`[1], 4), "

## 2. Combined Dataset (Oral + Fecal)
- Combined samples: ", nrow(combined_abundance), "
- Oral samples: ", sum(combined_meta$SampleType == "Oral"), "
- Fecal samples: ", sum(combined_meta$SampleType == "Fecal"), "
- Common genera analyzed: ", length(common_genera), "

### PERMANOVA (sequential):
- Body Site (SampleType): p = ", round(permanova_comb$`Pr(>F)`[1], 4), "
- Disease (Group): p = ", round(permanova_comb$`Pr(>F)`[2], 4), "
- Interaction (SampleType:Group): p = ", round(permanova_comb$`Pr(>F)`[3], 4), "
")

writeLines(summary_stats, file.path(desktop_dir, "fecal_combined_analysis_summary.md"))

cat("\n=== Analysis Complete! Results saved to ", desktop_dir, " ===\n")
