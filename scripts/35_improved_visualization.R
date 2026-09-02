#!/usr/bin/env Rscript
# Enhanced Visualization Script for Manuscript Figures

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(ggpubr)
  library(pheatmap)
  library(vegan)
})

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/manuscript_figures"
dir.create(output_dir, recursive = TRUE)

asv_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_table.tsv"
tax_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_taxonomy.tsv"
meta_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"

cat("Loading data...\n")
asv <- read.table(asv_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax <- read.table(tax_path, sep = "\t", header = TRUE, row.names = 1)
meta <- read.table(meta_path, sep = "\t", header = TRUE)

if (grepl("^SRR", rownames(asv)[1])) {
  asv <- t(as.matrix(asv))
}
asv <- as.matrix(asv)
tax <- as.matrix(tax)

ps <- phyloseq(otu_table(asv, taxa_are_rows = FALSE),
               tax_table(tax),
               sample_data(meta))

sample_data(ps)$Group <- factor(sample_data(ps)$Group, levels = c("control", "case"))

p1 <- plot_richness(ps, x = "Group", measures = c("Shannon", "Observed")) +
  geom_boxplot(aes(fill = Group), alpha = 0.7) +
  scale_fill_manual(values = c("case" = "#E64B35", "control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Alpha Diversity") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(p1, file.path(output_dir, "alpha_diversity.png"), width = 8, height = 4, dpi = 300)

pcoa <- ordinate(ps, method = "PCoA", distance = "bray")
p2 <- plot_ordination(ps, pcoa, color = "Group", shape = "Group") +
  scale_color_manual(values = c("case" = "#E64B35", "control" = "#4DBBD5")) +
  theme_bw() + ggtitle("Beta Diversity (Bray-Curtis)") +
  geom_point(size = 3) + stat_ellipse(level = 0.95)
ggsave(p2, file.path(output_dir, "beta_diversity_pcoa.png"), width = 6, height = 5, dpi = 300)

ps_genus <- tax_glom(ps, taxrank = "Genus")
top_genera <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:20]
ps_top <- prune_taxa(top_genera, ps_genus)
p3 <- plot_bar(ps_top, x = "Sample", fill = "Genus") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  ggtitle("Top 20 Genera Abundance")
ggsave(p3, file.path(output_dir, "top_genera_barplot.png"), width = 14, height = 6, dpi = 300)

p4 <- plot_heatmap(ps_top, sample.label = "Group", taxa.label = "Genus",
                   method = "bray", distance = "bray",
                   sample.order = "Group")
ggsave(p4, file.path(output_dir, "genera_heatmap.png"), width = 10, height = 8, dpi = 300)

taxa <- tax_table(ps)
genus_names <- taxa[, "Genus"]
abundances <- otu_table(ps)
total_abund <- rowSums(abundances)
genus_abund <- rowsum(abundances, genus_names)
genus_abund <- genus_abund[order(rowSums(genus_abund), decreasing = TRUE), ]
top10_genus <- rownames(genus_abund)[1:10]

case_idx <- which(sample_data(ps)$group == "case")
ctrl_idx <- which(sample_data(ps)$group == "control")

genus_df <- data.frame(
  Genus = rep(top10_genus, each = 2),
  Group = rep(c("Case", "Control"), each = 10),
  Mean = c(
    rowMeans(genus_abund[top10_genus, case_idx]),
    rowMeans(genus_abund[top10_genus, ctrl_idx])
  )
)

p5 <- ggplot(genus_df, aes(x = Genus, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Case" = "#E64B35", "Control" = "#4DBBD5")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 Genera: Case vs Control") + ylab("Mean Abundance")
ggsave(p5, file.path(output_dir, "top10_genera_comparison.png"), width = 10, height = 6, dpi = 300)

cat("Figures saved to:", output_dir, "\n")