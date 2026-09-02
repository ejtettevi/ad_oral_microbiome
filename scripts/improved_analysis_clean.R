#!/usr/bin/env Rscript
# Comprehensive Improved Analysis - Addressing All Reviewer Comments

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(vegan)
  library(randomForest)
})

cat("=== Comprehensive Improved Analysis v2 ===\n")
cat("Addressing reviewer comments 3-10\n\n")

output_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(output_dir, recursive = TRUE)

# Load Data
cat("1. Loading data...\n")
case_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv"
ctrl_path <- "/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv"

case_df <- read.csv(case_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ctrl_df <- read.csv(ctrl_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Ensure numeric
case_mat <- apply(case_df, 2, function(x) as.numeric(as.character(x)))
case_mat[is.na(case_mat)] <- 0
ctrl_mat <- apply(ctrl_df, 2, function(x) as.numeric(as.character(x)))
ctrl_mat[is.na(ctrl_mat)] <- 0

# Check for valid rownames
if (is.null(rownames(case_mat)) || all(rownames(case_mat) == "")) {
  cat("Row names not set, using first column as taxa names...\n")
  case_mat <- case_df
  case_mat$Taxon <- NULL
  case_mat <- apply(case_mat, 2, as.numeric)
  case_mat[is.na(case_mat)] <- 0
}

if (is.null(rownames(ctrl_mat)) || all(rownames(ctrl_mat) == "")) {
  ctrl_mat <- ctrl_df
  ctrl_mat$Taxon <- NULL
  ctrl_mat <- apply(ctrl_mat, 2, as.numeric)
  ctrl_mat[is.na(ctrl_mat)] <- 0
}

all_genera <- intersect(rownames(case_mat), rownames(ctrl_mat))
cat("Common genera:", length(all_genera), "\n")
combined <- cbind(case_mat[all_genera, ], ctrl_mat[all_genera, ])
group <- c(rep("Case", ncol(case_mat)), rep("Control", ncol(ctrl_mat)))

cat("Samples:", length(group), "(Case:", sum(group=="Case"), ", Control:", sum(group=="Control"), ")\n")
cat("Genera:", length(all_genera), "\n\n")

# CLR transformation
clr_transform <- function(x) {
  x_clr <- x + 1e-6
  gm <- exp(rowMeans(log(x_clr)))
  log(x_clr / gm)
}
clr_data <- clr_transform(combined)
rownames(clr_data) <- all_genera

# Wilcoxon on CLR
cat("2. Differential abundance (CLR-based)...\n")
results <- data.frame()
for (i in 1:nrow(clr_data)) {
  case_vals <- as.numeric(clr_data[i, group == "Case"])
  ctrl_vals <- as.numeric(clr_data[i, group == "Control"])
  
  res <- wilcox.test(case_vals, ctrl_vals)
  results <- rbind(results, data.frame(
    Taxon = all_genera[i],
    CLR_Case = mean(case_vals),
    CLR_Control = mean(ctrl_vals),
    CLR_diff = mean(case_vals) - mean(ctrl_vals),
    p_value = res$p.value
  ))
}
results$q_value <- p.adjust(results$p_value, "fdr")
results$sig <- results$q_value < 0.1
results <- results[order(results$p_value), ]

cat("Significant taxa (FDR<0.1):", sum(results$sig), "\n")
write.csv(results, file.path(output_dir, "differential_abundance_clr.csv"), row.names = FALSE)

# Alpha diversity
cat("\n3. Alpha diversity (Observed, Shannon, Simpson)...\n")
observed <- colSums(combined > 0)
shannon <- diversity(t(combined), index = "shannon")
simpson <- diversity(t(combined), index = "simpson")

alpha_df <- data.frame(Sample = colnames(combined), Group = group,
                      Observed = observed, Shannon = shannon, Simpson = simpson)

for (metric in c("Observed", "Shannon", "Simpson")) {
  p <- ggplot(alpha_df, aes_string(x = "Group", y = metric, fill = "Group")) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
    theme_bw() + ggtitle(paste(metric, "Diversity (Genus Level)")) +
    stat_compare_means(method = "wilcox.test", label = "p.format")
  ggsave(file.path(output_dir, paste0("alpha_", tolower(metric), ".png")), p, width = 5, height = 4, dpi = 300)
}

# Beta diversity - Aitchison (CLR-based Euclidean)
cat("\n4. Beta diversity (Aitchison distance, genus level)...\n")
dist_aitchison <- dist(t(clr_data), method = "euclidean")
pcoa_res <- cmdscale(dist_aitchison, k = 2)
pcoa_df <- data.frame(PC1 = pcoa_res[,1], PC2 = pcoa_res[,2], Group = group)

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) + stat_ellipse() +
  scale_color_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Beta Diversity (Aitchison Distance, Genus Level)") +
  labs(x = "PC1", y = "PC2")
ggsave(file.path(output_dir, "beta_diversity_pcoa.png"), p, width = 6, height = 5, dpi = 300)

perm_result <- adonis2(dist_aitchison ~ group)
cat("PERMANOVA p-value:", perm_result$Pr[1], "\n")

# ML with CLR
cat("\n5. Random Forest classification (CLR-based)...\n")
ml_data <- t(clr_data)
ml_data[is.na(ml_data)] <- 0
group_factor <- as.factor(group)

set.seed(42)
rf_model <- randomForest(ml_data, group_factor, importance = TRUE, ntree = 500)
cat("OOB error rate:", round(rf_model$err.rate[nrow(rf_model$err.rate), 1]*100, 2), "%\n")

imp <- importance(rf_model)
imp_df <- data.frame(Taxon = all_genera, Importance = imp[,1])
imp_df <- imp_df[order(-imp_df$Importance), ]
cat("Top 10 important taxa:\n")
print(head(imp_df, 10))
write.csv(imp_df, file.path(output_dir, "ml_feature_importance.csv"), row.names = FALSE)

# CLR effect size plot
cat("\n6. CLR effect sizes plot...\n")
top20 <- head(results, 20)
top20$direction <- ifelse(top20$CLR_diff > 0, "Enriched in AD", "Enriched in Control")

p <- ggplot(top20, aes(x = reorder(Taxon, CLR_diff), y = CLR_diff, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Enriched in AD" = "#E64B35", "Enriched in Control" = "#4DBBD5")) +
  theme_bw() + coord_flip() +
  ggtitle("CLR Effect Sizes (Top 20 Genera)") +
  xlab("Genus") + ylab("CLR Difference (Case - Control)")
ggsave(file.path(output_dir, "clr_effect_sizes.png"), p, width = 8, height = 6, dpi = 300)

# Summary
cat("\n=== Summary ===\n")
cat("Analysis complete!\n")
cat("Output directory:", output_dir, "\n")

summary_text <- paste0("
# Analysis Summary - Addressing Reviewer Comments

## Methods (Addressing Comments 4-10):
- CLR transformation applied to all abundance data
- Wilcoxon test on CLR-transformed data (replacing ANCOMBC2)
- Report CLR values instead of percentages
- Multiple alpha diversity metrics: Observed, Shannon, Simpson
- Beta diversity at genus level (Aitchison distance on CLR)
- Random Forest on CLR-transformed data with importance scores

## Key Changes from Original:
- Report CLR values instead of percentages
- Specify taxonomy level (genus) throughout
- Use Aitchison distance for beta diversity (compatible with CLR)
- Report feature importance from ML

## Output Files:
- differential_abundance_clr.csv (CLR values + p/q values)
- ml_feature_importance.csv
- alpha_observed.png, alpha_shannon.png, alpha_simpson.png
- beta_diversity_pcoa.png
- clr_effect_sizes.png
")

writeLines(summary_text, file.path(output_dir, "analysis_summary_v2.md"))
cat("\nDone!\n")