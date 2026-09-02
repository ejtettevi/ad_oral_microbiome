#!/usr/bin/env Rscript
library(phyloseq)
library(tidyverse)
library(vegan)
library(permute)

cat("=== Full Microbiome Analysis ===\n\n")

output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/fecal_processed"
desktop_dir <- "~/Desktop/improved_microbiome_analysis"
dir.create(desktop_dir, showWarnings = FALSE)

cat("=== Loading Data ===\n")
fecal_seqtab <- readRDS(file.path(output_dir, "seqtab.rds"))
cat("Fecal: ", nrow(fecal_seqtab), " samples, ", ncol(fecal_seqtab), " ASVs\n", sep="")

run_table <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/raw_fastq/run_table/complete/complete_SraRunTable.csv")

fecal_runs <- rownames(fecal_seqtab)
fecal_meta <- run_table[run_table$Run %in% fecal_runs, ]
cat("Fecal metadata matched:", nrow(fecal_meta), "samples\n")

if(nrow(fecal_meta) > 0) {
  cat("Groups:", table(fecal_meta$Group), "\n")
}

fecal_phy <- phyloseq(otu_table(fecal_seqtab, taxa_are_rows = FALSE))

cat("\n=== Alpha Diversity (Fecal) ===\n")
alpha <- estimate_richness(fecal_phy, measures = c("Observed", "Shannon", "Simpson"))
alpha$Sample <- rownames(alpha)
alpha$Group <- sapply(alpha$Sample, function(x) {
  idx <- match(x, fecal_meta$Run)
  if(!is.na(idx)) return(fecal_meta$Group[idx])
  return(NA)
})
alpha <- alpha[!is.na(alpha$Group), ]

shannon_p <- wilcox.test(Shannon ~ Group, data = alpha)$p.value
observed_p <- wilcox.test(Observed ~ Group, data = alpha)$p.value
cat("Shannon p =", round(shannon_p, 4), "\n")
cat("Observed p =", round(observed_p, 4), "\n")

write.csv(alpha, file.path(desktop_dir, "fecal_alpha_diversity.csv"), row.names = FALSE)

cat("\n=== Beta Diversity (Fecal - Aitchison) ===\n")
fecal_genus <- tax_glom(fecal_phy, taxrank = "Genus")
fecal_mat <- as.matrix(otu_table(fecal_genus))

clr_transform <- function(x) {
  x[x == 0] <- 0.5
  log(x / exp(mean(log(x))))
}

clr_mat <- t(apply(fecal_mat + 1, 1, clr_transform))
colnames(clr_mat) <- paste0("G", 1:ncol(clr_mat))

pca <- prcomp(clr_mat)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)[1:2]
cat("PC1:", var_explained[1], "%  PC2:", var_explained[2], "%\n")

dist_mat <- dist(clr_mat)
group_labels <- sapply(rownames(clr_mat), function(x) {
  idx <- match(x, fecal_meta$Run)
  if(!is.na(idx)) return(fecal_meta$Group[idx])
  return(NA)
})
adonis_result <- adonis2(dist_mat ~ group_labels, permutations = 999)
cat("PERMANOVA p =", adonis_result$`Pr(>F)`[1], "\n")

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group_labels)
write.csv(pca_df, file.path(desktop_dir, "fecal_beta_pca.csv"), row.names = FALSE)

cat("\n=== Differential Abundance (Fecal - CLR) ===\n")
fecal_clr <- t(clr_mat)
da_results <- data.frame(Genus = rownames(fecal_clr))

for(i in 1:nrow(fecal_clr)) {
  g <- rownames(fecal_clr)[i]
  vals <- fecal_clr[i, ]
  case_idx <- group_labels == "Case"
  ctrl_idx <- group_labels == "Control"
  w <- wilcox.test(vals[case_idx], vals[ctrl_idx])
  da_results$Mean_Case[i] <- mean(vals[case_idx])
  da_results$Mean_Control[i] <- mean(vals[ctrl_idx])
  da_results$Wilcox_p[i] <- w$p.value
}

da_results$Wilcox_padj <- p.adjust(da_results$Wilcox_p, "fdr")
da_results <- da_results[order(da_results$Wilcox_padj), ]
sig_da <- da_results[da_results$Wilcox_padj < 0.1, ]
cat("Significant genera (FDR < 0.1):", nrow(sig_da), "\n")

write.csv(da_results, file.path(desktop_dir, "fecal_differential_abundance.csv"), row.names = FALSE)

cat("\n=== Oral Data Analysis ===\n")
oral_case <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_case.csv", row.names = 1)
oral_ctrl <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/oral_data/genus/oral_ad_ctrl.csv", row.names = 1)

case_ids <- colnames(oral_case)
ctrl_ids <- colnames(oral_ctrl)
oral_mat <- cbind(oral_case, oral_ctrl)
oral_groups <- c(rep("Case", length(case_ids)), rep("Control", length(ctrl_ids)))

oral_clr <- t(apply(oral_mat + 1, 1, clr_transform))
colnames(oral_clr) <- paste0("OG", 1:ncol(oral_clr))

oral_pca <- prcomp(oral_clr)
oral_var <- round(100 * oral_pca$sdev^2 / sum(oral_pca$sdev^2), 1)[1:2]
cat("Oral PC1:", oral_var[1], "%  PC2:", oral_var[2], "%\n")

oral_dist <- dist(oral_clr)
adonis2(oral_dist ~ oral_groups, permutations = 999)

cat("\n=== Combined Analysis ===\n")
n_oral <- ncol(oral_mat)
n_fecal <- nrow(fecal_clr)

oral_for_comb <- data.frame(SampleType = "Oral", Group = oral_groups, oral_clr[, 1:min(50, ncol(oral_clr))])
fecal_for_comb <- data.frame(SampleType = "Fecal", Group = group_labels[1:n_fecal], fecal_clr[, 1:min(50, ncol(fecal_clr))])

combined_clr <- rbind(oral_for_comb, fecal_for_comb)
combined_clr$Group <- as.factor(combined_clr$Group)
combined_clr$SampleType <- as.factor(combined_clr$SampleType)

combined_dist <- dist(combined_clr[, 3:ncol(combined_clr)])
adonis2(combined_dist ~ Group + SampleType, data = combined_clr, permutations = 999)

cat("\n=== Summary ===\n")
cat("Fecal samples analyzed:", n_fecal, "\n")
cat("Oral samples analyzed:", n_oral, "\n")
cat("Fecal significant DA genera:", nrow(sig_da), "\n")
cat("PERMANOVA (Fecal): p =", adonis_result$`Pr(>F)`[1], "\n")
cat("PERMANOVA (Oral): p =", adonis_result$`Pr(>F)`[1], "\n")

cat("\nFiles saved to:", desktop_dir, "\n")

pdf(file.path(desktop_dir, "fecal_alpha_diversity_plot.pdf"), width=10, height=4)
par(mfrow=c(1,3))
boxplot(Observed ~ Group, data=alpha, main="Observed", col=c("red","blue"))
boxplot(Shannon ~ Group, data=alpha, main="Shannon", col=c("red","blue"))
boxplot(Simpson ~ Group, data=alpha, main="Simpson", col=c("red","blue"))
dev.off()

pdf(file.path(desktop_dir, "fecal_beta_diversity_plot.pdf"), width=6, height=5)
plot(pca$x[,1], pca$x[,2], col=ifelse(group_labels=="Case", "red", "blue"), pch=19, xlab=paste0("PC1 (", var_explained[1], "%)"), ylab=paste0("PC2 (", var_explained[2], "%)"), main="Fecal Beta Diversity")
legend("topright", c("Case","Control"), col=c("red","blue"), pch=19)
dev.off()

cat("\n=== Complete ===\n")