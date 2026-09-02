#!/usr/bin/env Rscript
library(phyloseq)
library(tidyverse)
library(vegan)
library(DESeq2)

cat("=== Combined Oral + Fecal Analysis ===\n\n")

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed"
desktop_dir <- "~/Desktop/improved_microbiome_analysis"

fecal_seqtab <- readRDS(file.path(output_dir, "seqtab.rds"))
cat("Fecal samples:", nrow(fecal_seqtab), "ASVs:", ncol(fecal_seqtab), "\n")

run_table <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv")
fecal_meta <- run_table[run_table$Run %in% rownames(fecal_seqtab), ]
cat("Fecal metadata:", nrow(fecal_meta), "samples\n")

oral_asv <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv", row.names = 1)
oral_ctrl <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv", row.names = 1)
oral_otu <- cbind(oral_asv, oral_ctrl)

case_ids <- colnames(oral_asv)
ctrl_ids <- colnames(oral_ctrl)

oral_meta <- data.frame(
  Run = c(case_ids, ctrl_ids),
  Group = c(rep("Case", length(case_ids)), rep("Control", length(ctrl_ids))),
  SampleType = "Oral"
)
fecal_meta$SampleType <- "Fecal"

all_meta <- rbind(oral_meta, fecal_meta[, c("Run", "Group", "SampleType")])

fecal_genus <- tax_glom(fecal_seqtab, taxrank = "Genus")
fecal_mat <- as.matrix(t(fecal_genus))
fecal_mat <- fecal_mat[rownames(fecal_mat) %in% all_meta$Run, ]

oral_mat <- as.matrix(oral_otu)
rownames(oral_mat) <- gsub("_1|_2", "", rownames(oral_mat))

common_genus <- intersect(colnames(fecal_mat), colnames(oral_mat))
cat("Common genera:", length(common_genus), "\n")

combined_mat <- rbind(fecal_mat[, common_genus], oral_mat[, common_genus])
rownames(combined_mat) <- c(rownames(fecal_mat), rownames(oral_mat))

combined_meta <- data.frame(
  Sample = rownames(combined_mat),
  Group = sapply(rownames(combined_mat), function(x) {
    if(x %in% case_ids) return("Case")
    if(x %in% ctrl_ids) return("Control")
    if(x %in% fecal_meta$Run) return(fecal_meta$Group[match(x, fecal_meta$Run)])
    return("Unknown")
  }),
  SampleType = sapply(rownames(combined_mat), function(x) {
    if(x %in% fecal_meta$Run) return("Fecal") else return("Oral")
  })
)

cat("\n=== Combined Dataset ===\n")
cat("Total samples:", nrow(combined_mat), "\n")
cat("Cases:", sum(combined_meta$Group == "Case"), "\n")
cat("Controls:", sum(combined_meta$Group == "Control"), "\n")
cat("Oral:", sum(combined_meta$SampleType == "Oral"), "\n")
cat("Fecal:", sum(combined_meta$SampleType == "Fecal"), "\n")

combined_mat <- combined_mat[combined_meta$Sample, ]

clr <- function(x) {
  x[x == 0] <- 0.5
  log(x / exp(mean(log(x))))
}

clr_mat <- apply(combined_mat + 1, 1, clr)
clr_mat <- t(clr_mat)

cat("\n=== Beta Diversity (Aitchison Distance) ===\n")
pca <- prcomp(clr_mat)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat("PC1:", var_explained[1], "%  PC2:", var_explained[2], "%\n")

dist_mat <- dist(clr_mat)
adonis2(dist_mat ~ Group, data = combined_meta, permutations = 999)

cat("\n=== Differential Abundance (Fecal Only) ===\n")
fecal_clr <- apply(fecal_mat + 1, 1, clr)
fecal_clr <- t(fecal_clr)

fecal_da <- data.frame(Genus = colnames(fecal_clr))
for(g in colnames(fecal_clr)) {
  case_vals <- fecal_clr[fecal_meta$Group[fecal_meta$Run %in% rownames(fecal_clr)] == "Case", g]
  ctrl_vals <- fecal_clr[fecal_meta$Group[fecal_meta$Run %in% rownames(fecal_clr)] == "Control", g]
  w <- wilcox.test(case_vals, ctrl_vals)
  fecal_da$Wilcox_p[fecal_da$Genus == g] <- w$p.value
  fecal_da$Mean_Case[fecal_da$Genus == g] <- mean(case_vals)
  fecal_da$Mean_Control[fecal_da$Genus == g] <- mean(ctrl_vals)
}
fecal_da$Wilcox_padj <- p.adjust(fecal_da$Wilcox_p, "fdr")
fecal_da <- fecal_da[order(fecal_da$Wilcox_padj), ]
cat("Significant genera (FDR < 0.1):", sum(fecal_da$Wilcox_padj < 0.1), "\n")

cat("\n=== Alpha Diversity (Fecal) ===\n")
fecal_otu <- otu_table(fecal_mat, taxa_are_rows = FALSE)
fecal_alpha <- estimate_richness(fecal_otu, measures = c("Observed", "Shannon"))
fecal_alpha$Sample <- rownames(fecal_alpha)
fecal_alpha$Group <- sapply(fecal_alpha$Sample, function(x) {
  fecal_meta$Group[match(x, fecal_meta$Run)]
})

cat("Shannon: Case vs Control p =", wilcox.test(Shannon ~ Group, data = fecal_alpha)$p.value, "\n")

cat("\n=== Saving Results ===\n")
dir.create(desktop_dir, showWarnings = FALSE)
write.csv(fecal_da, file.path(desktop_dir, "fecal_differential_abundance.csv"), row.names = FALSE)
write.csv(pca$x[, 1:3], file.path(desktop_dir, "fecal_pca_scores.csv"))
write.csv(fecal_alpha, file.path(desktop_dir, "fecal_alpha_diversity.csv"), row.names = FALSE)

combined_df <- data.frame(Sample = rownames(combined_mat), combined_meta, clr_mat)
write.csv(combined_df, file.path(desktop_dir, "combined_clr_data.csv"), row.names = FALSE)

cat("Files saved to:", desktop_dir, "\n")
cat("\n=== Complete ===\n")