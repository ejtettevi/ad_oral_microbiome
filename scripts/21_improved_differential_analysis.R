#!/usr/bin/env Rscript
# Improved Differential Abundance Analysis for Oral Microbiome
# Adds: ANCOM-BC2, effect sizes, better FDR control, comprehensive reporting

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(edgeR)
  library(DESeq2)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(ggpubr)
})

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/improved_analysis"
dir.create(output_dir, recursive = TRUE)

asv_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_table.tsv"
tax_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_taxonomy.tsv"
meta_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"

cat("Loading data...\n")
asv <- read.table(asv_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax <- read.table(tax_path, sep = "\t", header = TRUE, row.names = 1)
meta <- read.table(meta_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meta$Group <- tolower(meta$Group)
rownames(meta) <- meta$Sample

if (grepl("^SRR", rownames(asv)[1])) {
  asv <- t(as.matrix(asv))
}
asv <- as.matrix(asv)
tax <- as.matrix(tax)

ps <- phyloseq(otu_table(asv, taxa_are_rows = TRUE),
               tax_table(tax),
               sample_data(meta))

sample_data(ps)$Group <- factor(sample_data(ps)$Group, levels = c("control", "case"))

cat("Sample sizes:\n")
print(table(sample_data(ps)$Group))

cat("\nRunning ANCOM-BC2 analysis...\n")
ancombc_result <- ancombc2(ps, fix_formula = "Group", 
                           rand_formula = NULL,
                           p_adj_method = "fdr",
                           lib_cut = 0,
                           group = "Group",
                           global = FALSE,
                           pairwise = FALSE,
                           dunnet = FALSE,
                           trend = FALSE,
                           n_cl = 1)

res_ancom <- ancombc_result$res
res_ancom$lfc <- res_ancom$lfc_Groupcase
res_ancom$se <- res_ancom$se_Groupcase
res_ancom$p_val <- res_ancom$p_Groupcase
res_ancom$q_adj <- res_ancom$q_Groupcase
res_ancom$sig <- res_ancom$diff_Groupcase

cat("\n=== ANCOM-BC2 Results ===\n")
sig_taxa <- res_ancom[res_ancom$sig, ]
cat("Significant taxa (q < 0.1):", nrow(sig_taxa), "\n")

if (nrow(sig_taxa) > 0) {
  cat("\nSignificant taxa with direction:\n")
  print(sig_taxa[, c("taxon", "lfc", "p_val", "q_adj")])
}

effect_sizes <- res_ancom %>%
  mutate(
    direction = ifelse(lfc > 0, "Enriched in Case", "Enriched in Control"),
    abs_effect = abs(lfc)
  ) %>%
  arrange(q_adj)

tax_df <- as.data.frame(tax)
tax_df$taxon <- rownames(tax_df)
effect_sizes_annotated <- effect_sizes %>%
  left_join(tax_df, by = "taxon")

write.csv(effect_sizes_annotated, file.path(output_dir, "ancombc2_results.csv"), row.names = FALSE)

cat("\n=== Effect Size Summary ===\n")
cat("Mean absolute effect:", mean(effect_sizes$abs_effect, na.rm = TRUE), "\n")
cat("Max effect:", max(effect_sizes$abs_effect, na.rm = TRUE), "\n")

cat("\n=== Beta-diversity analysis ===\n")
ord <- ordinate(ps, method = "PCoA", distance = "bray")
p_ord <- plot_ordination(ps, ord, color = "Group", shape = "Group") +
  scale_color_manual(values = c("case" = "#E64B35", "control" = "#4DBBD5")) +
  theme_bw() + geom_point(size = 3) + stat_ellipse()
ggsave(file.path(output_dir, "pcoa_bray.png"), p_ord, width = 6, height = 5, dpi = 300)

cat("\nAnalysis complete. Results saved to:", output_dir, "\n")