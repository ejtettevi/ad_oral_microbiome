#!/usr/bin/env Rscript
# Comprehensive Improved Analysis - Addressing All Reviewer Comments
# Version 2: CLR-based analysis, proper reporting, ML improvements

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(pheatmap)
  library(vegan)
  library(randomForest)
  library(lefser)
})

cat("=== Comprehensive Improved Analysis v2 ===\n")
cat("Addressing reviewer comments 3-10\n\n")

output_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(output_dir, recursive = TRUE)

# ============================================================================
# Load Data (using existing oral data)
# ============================================================================
cat("1. Loading data...\n")
case_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv"
ctrl_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv"

case_df <- read.csv(case_path, row.names = 1)
ctrl_df <- read.csv(ctrl_path, row.names = 1)

all_genera <- intersect(rownames(case_df), rownames(ctrl_df))
cat("Genera in dataset:", length(all_genera), "\n")

# ============================================================================
# COMMENT 4: ANCOMBC2 with CLR values instead of percentages
# ============================================================================
cat("\n2. ANCOMBC2 analysis with CLR reporting...\n")

# Prepare data
case_mat <- apply(case_df[all_genera, ], 2, as.numeric)
case_mat[is.na(case_mat)] <- 0
ctrl_mat <- apply(ctrl_df[all_genera, ], 2, as.numeric)
ctrl_mat[is.na(ctrl_mat)] <- 0
combined <- cbind(case_mat, ctrl_mat)

group <- c(rep("Case", ncol(case_mat)), rep("Control", ncol(ctrl_mat)))
metadata <- data.frame(Sample = colnames(combined), Group = group, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Sample

# CLR transformation function
clr_transform <- function(x) {
  x_clr <- x + 1e-6
  gm <- exp(rowMeans(log(x_clr)))
  log(x_clr / gm)
}

# Transform to CLR
clr_data <- clr_transform(combined)
colnames(clr_data) <- colnames(combined)

# Run ANCOMBC2 - make sure Group is proper factor
metadata$Group <- as.factor(metadata$Group)
ps <- phyloseq(otu_table(combined, taxa_are_rows = TRUE),
               sample_data(metadata))

tryCatch({
  ancom_result <- ancombc2(ps, fix_formula = "Group",
                          p_adj_method = "fdr", lib_cut = 0,
                          global = TRUE, n_cl = 1)
  
  res <- ancom_result$res
  res$q_adj <- p.adjust(res$p_val, "fdr")
  res$sig <- res$q_adj < 0.1
  
  # Report CLR values (mean CLR for each group)
  case_clr <- rowMeans(clr_data[, metadata$Group == "Case"])
  ctrl_clr <- rowMeans(clr_data[, metadata$Group == "Control"])
  
  res$CLR_Case <- case_clr
  res$CLR_Control <- ctrl_clr
  res$CLR_diff <- case_clr - ctrl_clr
  
  cat("ANCOMBC2 results:\n")
  sig_res <- res[res$sig, c("taxon", "CLR_Case", "CLR_Control", "CLR_diff", "p_val", "q_adj")]
  print(head(sig_res, 10))
  
  # Save with CLR
  write.csv(res, file.path(output_dir, "ancombc2_clr_results.csv"), row.names = FALSE)
  
}, error = function(e) {
  cat("ANCOMBC2 failed, using Wilcoxon on CLR...\n")
  cat("Error:", e$message, "\n")
  
  results <- data.frame()
  all_genera_clr <- rownames(clr_data)
  
  for (i in seq_along(all_genera_clr)) {
    genus <- all_genera_clr[i]
    case_vals <- clr_data[i, metadata$Group == "Case"]
    ctrl_vals <- clr_data[i, metadata$Group == "Control"]
    
    pval <- tryCatch(wilcox.test(case_vals, ctrl_vals)$p.value, error = function(x) NA)
    if (!is.na(pval)) {
      results <- rbind(results, data.frame(
        Taxon = genus,
        CLR_Case = mean(case_vals),
        CLR_Control = mean(ctrl_vals),
        CLR_diff = mean(case_vals) - mean(ctrl_vals),
        p_value = pval,
        q_value = p.adjust(pval, "fdr")
      ))
    }
  }
  results$sig <- results$q_value < 0.1
  
  write.csv(results, file.path(output_dir, "wilcoxon_clr_results.csv"), row.names = FALSE)
  cat("Saved wilcoxon_clr_results.csv\n")
})

# ============================================================================
# COMMENT 5 & 6: CLR-based visualizations
# ============================================================================
cat("\n3. Creating CLR-based visualizations...\n")

# PCA on CLR data
pca_result <- prcomp(t(clr_data), center = TRUE, scale = FALSE)
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Group = group)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) + stat_ellipse() +
  scale_color_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("PCA on CLR-transformed Genus Abundances") +
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)"))
ggsave(file.path(output_dir, "pca_clr.png"), p1, width = 6, height = 5, dpi = 300)

# ============================================================================
# COMMENT 7: Alpha diversity at multiple levels with proper reporting
# ============================================================================
cat("\n4. Alpha diversity analysis (Observed, Shannon, Simpson)...\n")

# Transform to relative abundance for diversity calc
rel_abund <- sweep(combined, 2, colSums(combined), "/")

# Observed (richness)
observed <- colSums(combined > 0)

# Shannon
shannon <- diversity(t(combined), index = "shannon")

# Simpson
simpson <- diversity(t(combined), index = "simpson")

alpha_df <- data.frame(
  Sample = colnames(combined),
  Group = group,
  Observed = observed,
  Shannon = shannon,
  Simpson = simpson
)

p2a <- ggplot(alpha_df, aes(x = Group, y = Observed, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Observed Genera (Genus Richness)") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(file.path(output_dir, "alpha_observed.png"), p2a, width = 5, height = 4, dpi = 300)

p2b <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Shannon Diversity") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(file.path(output_dir, "alpha_shannon.png"), p2b, width = 5, height = 4, dpi = 300)

p2c <- ggplot(alpha_df, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Simpson Diversity") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(file.path(output_dir, "alpha_simpson.png"), p2c, width = 5, height = 4, dpi = 300)

# ============================================================================
# COMMENT 8: Beta diversity at genus level
# ============================================================================
cat("\n5. Beta diversity at genus level...\n")

# Use Aitchison distance (Euclidean on CLR) instead of Bray-Curtis
# Bray-Curtis doesn't work with CLR (negative values)
dist_aitchison <- dist(t(clr_data), method = "euclidean")

# PCoA
pcoa_res <- cmdscale(dist_aitchison, k = 2)
pcoa_df <- data.frame(PC1 = pcoa_res[,1], PC2 = pcoa_res[,2], Group = group)

p3 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) + stat_ellipse() +
  scale_color_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Beta Diversity (Aitchison Distance, Genus Level, CLR-transformed)") +
  labs(x = "PC1", y = "PC2")
ggsave(file.path(output_dir, "beta_diversity_pcoa.png"), p3, width = 6, height = 5, dpi = 300)

# PERMANOVA
perm_result <- adonis2(dist_aitchison ~ group)
cat("PERMANOVA p-value:", perm_result$Pr[1], "\n")

# ============================================================================
# COMMENT 10: Machine Learning with CLR data
# ============================================================================
cat("\n6. Random Forest classification with CLR data...\n")

# Prepare ML data (use CLR)
ml_data <- t(clr_data)
ml_data[is.na(ml_data)] <- 0

# Make sure group is a factor
group_factor <- as.factor(group)

# Random Forest
set.seed(42)
rf_model <- randomForest(ml_data, group_factor, importance = TRUE, ntree = 500)

cat("OOB error rate:", rf_model$err.rate[nrow(rf_model$err.rate), 1], "\n")

# Feature importance
imp <- importance(rf_model)
imp_df <- data.frame(
  Taxon = rownames(imp),
  Importance = imp[, 1],
  MeanDecreaseAccuracy = imp[, 2],
  MeanDecreaseGini = imp[, 3]
)
imp_df <- imp_df[order(-imp_df$Importance), ]

cat("Top 10 important taxa:\n")
print(head(imp_df, 10))

write.csv(imp_df, file.path(output_dir, "ml_feature_importance.csv"), row.names = FALSE)

# ============================================================================
# COMMENT 5: LDA-like visualization (using effect sizes)
# ============================================================================
cat("\n7. Creating LDA-style effect size plot...\n")

# Get effect sizes from Wilcoxon test
effect_df <- data.frame()
all_genera_clr <- rownames(clr_data)
for (i in seq_along(all_genera_clr)) {
  genus <- all_genera_clr[i]
  case_vals <- clr_data[i, metadata$Group == "Case"]
  ctrl_vals <- clr_data[i, metadata$Group == "Control"]
  
  effect <- mean(case_vals) - mean(ctrl_vals)
  pval <- tryCatch(wilcox.test(case_vals, ctrl_vals)$p.value, error = function(x) NA)
  if (!is.na(pval)) {
    effect_df <- rbind(effect_df, data.frame(
      Taxon = genus,
      CLR_effect = effect,
      p_value = pval
    ))
  }
}
effect_df$q_value <- p.adjust(effect_df$p_value, "fdr")
effect_df$sig <- effect_df$q_value < 0.1
effect_df$CLR_effect <- as.numeric(effect_df$CLR_effect)
effect_df <- effect_df[order(-abs(effect_df$CLR_effect)), ]

# Top 20 by absolute effect
top20 <- head(effect_df, 20)
top20$direction <- ifelse(top20$CLR_effect > 0, "Enriched in AD", "Enriched in Control")

p4 <- ggplot(top20, aes(x = reorder(Taxon, CLR_effect), y = CLR_effect, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Enriched in AD" = "#E64B35", "Enriched in Control" = "#4DBBD5")) +
  theme_bw() + coord_flip() +
  ggtitle("CLR Effect Sizes (Top 20 Genera)") +
  xlab("Genus") + ylab("CLR Difference (Case - Control)")
ggsave(file.path(output_dir, "clr_effect_sizes.png"), p4, width = 8, height = 6, dpi = 300)

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Summary ===\n")
cat("Analysis completed with CLR-based methods\n")
cat("All outputs saved to:", output_dir, "\n")

summary_text <- paste0("
# Analysis Summary - Addressing Reviewer Comments

## Methods Used:
1. **CLR transformation** applied to all abundance data
2. **ANCOMBC2** for differential abundance (or Wilcoxon on CLR)
3. **PCA/PCoA** on CLR-transformed data
4. **Random Forest** on CLR data with importance scores

## Key Changes from Original:
- Report CLR values instead of percentages
- Multiple alpha diversity metrics (Observed, Shannon, Simpson)
- Beta diversity at genus level with CLR transformation
- ML with CLR data, reporting feature importance

## Output Files:
- ancombc2_clr_results.csv or wilcoxon_clr_results.csv
- ml_feature_importance.csv
- pca_clr.png
- alpha_observed.png, alpha_shannon.png, alpha_simpson.png
- beta_diversity_pcoa.png
- clr_effect_sizes.png
")

writeLines(summary_text, file.path(output_dir, "analysis_summary_v2.md"))
cat("\nDone!\n")