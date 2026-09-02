#!/usr/bin/env Rscript
# Improved Analysis using existing genus-level data

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(pheatmap)
})

output_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(output_dir, recursive = TRUE)

cat("Loading genus-level data...\n")

case_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv"
ctrl_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv"

case_df <- read.csv(case_path, row.names = 1)
ctrl_df <- read.csv(ctrl_path, row.names = 1)

cat("Case samples:", ncol(case_df), "\n")
cat("Control samples:", ncol(ctrl_df), "\n")

all_genera <- intersect(rownames(case_df), rownames(ctrl_df))
cat("Common genera:", length(all_genera), "\n")

case_mat <- apply(case_df[all_genera, , drop = FALSE], 2, as.numeric)
case_mat[is.na(case_mat)] <- 0
rownames(case_mat) <- all_genera
ctrl_mat <- apply(ctrl_df[all_genera, , drop = FALSE], 2, as.numeric)
ctrl_mat[is.na(ctrl_mat)] <- 0
rownames(ctrl_mat) <- all_genera

combined <- cbind(case_mat, ctrl_mat)
group <- c(rep("Case", ncol(case_mat)), rep("Control", ncol(ctrl_mat)))

metadata <- data.frame(Sample = colnames(combined), Group = group)
rownames(metadata) <- metadata$Sample

asv_mat <- as.matrix(combined)
ps <- phyloseq(otu_table(asv_mat, taxa_are_rows = TRUE),
               sample_data(metadata))

cat("\nRunning ANCOM-BC2 analysis...\n")
tryCatch({
  ancombc_result <- ancombc2(ps, fix_formula = "~ Group",
                             rand_formula = NULL,
                             p_adj_method = "fdr",
                             lib_cut = 0,
                             global = TRUE,
                             n_cl = 1)
  
  res <- ancombc_result$res
  res$sig <- res$q_adj < 0.1
  
  cat("\n=== ANCOM-BC2 Results (FDR < 0.1) ===\n")
  sig_res <- res[res$sig, ]
  cat("Significant taxa:", nrow(sig_res), "\n")
  if (nrow(sig_res) > 0) {
    print(sig_res[, c("taxon", "Group", "diff", "p_val", "q_adj")])
  }
  
  write.csv(res, file.path(output_dir, "ancombc2_all_results.csv"), row.names = FALSE)
  
}, error = function(e) {
  cat("ANCOM-BC2 failed, using Wilcoxon test...\n")
  
  results <- data.frame()
  for (genus in all_genera) {
    case_vals <- as.numeric(case_mat[genus, ])
    ctrl_vals <- as.numeric(ctrl_mat[genus, ])
    
    case_nonzero <- sum(!is.na(case_vals) & case_vals > 0)
    ctrl_nonzero <- sum(!is.na(ctrl_vals) & ctrl_vals > 0)
    
    if (case_nonzero > 3 && ctrl_nonzero > 3) {
      pval <- tryCatch(wilcox.test(case_vals, ctrl_vals, na.action = "na.omit")$p.value, error = function(x) NA)
      if (!is.na(pval)) {
        effect <- mean(case_vals, na.rm = TRUE) - mean(ctrl_vals, na.rm = TRUE)
        results <- rbind(results, data.frame(
          Genus = genus,
          Mean_Case = mean(case_vals, na.rm = TRUE),
          Mean_Control = mean(ctrl_vals, na.rm = TRUE),
          Effect = effect,
          p_value = pval
        ))
      }
    }
  }
  
  results$fdr <- p.adjust(results$p_value, "fdr")
  results$sig <- results$fdr < 0.1
  
  cat("\n=== Wilcoxon Results (FDR < 0.1) ===\n")
  sig_results <- results[results$sig, ]
  cat("Significant genera:", nrow(sig_results), "\n")
  if (nrow(sig_results) > 0) {
    print(sig_results)
  }
  
  write.csv(results[order(results$p_value), ], file.path(output_dir, "wilcoxon_results.csv"), row.names = FALSE)
})

cat("\n=== Alpha Diversity ===\n")
richness <- data.frame(
  Sample = colnames(combined),
  Group = group,
  Richness = colSums(combined > 0),
  Total = colSums(combined)
)

p1 <- ggplot(richness, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Species Richness") +
  stat_compare_means(method = "wilcox.test")
ggsave(file.path(output_dir, "alpha_richness.png"), p1, width = 5, height = 4, dpi = 300)

p2 <- ggplot(richness, aes(x = Group, y = Total, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Total Abundance") +
  stat_compare_means(method = "wilcox.test")
ggsave(file.path(output_dir, "total_abundance.png"), p2, width = 5, height = 4, dpi = 300)

cat("\n=== Beta Diversity (PCoA) ===\n")
library(vegan)
dist_mat <-vegdist(t(combined), method = "bray")
pcoa_res <- cmdscale(dist_mat, k = 2)
pcoa_df <- data.frame(PC1 = pcoa_res[,1], PC2 = pcoa_res[,2], Group = group)

p3 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) + stat_ellipse() +
  scale_color_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("PCoA - Bray Curtis")
ggsave(file.path(output_dir, "pcoa_bray.png"), p3, width = 6, height = 5, dpi = 300)

perm_result <- adonis2(dist_mat ~ group)
cat("PERMANOVA p-value:", perm_result$Pr[1], "\n")

cat("\n=== Top Genera Comparison ===\n")
top_gen <- names(sort(rowMeans(combined), decreasing = TRUE))[1:15]
top_df <- data.frame(Genus = rep(top_gen, each = 2),
                     Group = rep(c("Case", "Control"), each = 15),
                     Mean = c(rowMeans(case_mat[top_gen, , drop = FALSE]),
                              rowMeans(ctrl_mat[top_gen, , drop = FALSE])))

p4 <- ggplot(top_df, aes(x = Genus, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 15 Genera")
ggsave(file.path(output_dir, "top_genera.png"), p4, width = 10, height = 5, dpi = 300)

cat("\n=== Generating Summary ===\n")
summary_text <- paste0("
# Oral Microbiome Analysis Summary
# Generated: ", Sys.time(), "

## Dataset
- Case samples: ", ncol(case_df), "
- Control samples: ", ncol(ctrl_df), "
- Total genera: ", length(all_genera), "

## Key Findings
- See detailed results in output files
- Alpha diversity comparison in alpha_richness.png
- Beta diversity (PCoA) in pcoa_bray.png
- Top genera comparison in top_genera.png

## Output Files
- ancombc2_all_results.csv or wilcoxon_results.csv
- alpha_richness.png
- total_abundance.png
- pcoa_bray.png
- top_genera.png
")

writeLines(summary_text, file.path(output_dir, "analysis_summary.md"))

cat("\nDone! Results saved to:", output_dir, "\n")