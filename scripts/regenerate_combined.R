suppressPackageStartupMessages({
  library(phyloseq); library(vegan)
})
out <- "/var/folders/45/7jcqmvc92yb4gl2q9jjw6d7h0000gn/T/opencode"
base <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir"

clr_transform <- function(x){ x[x==0] <- 0.5; log(x / exp(mean(log(x)))) }

# ---------- ORAL (96 samples, DADA2) ----------
oral_phy <- readRDS(file.path(base,"dada2/phyloseq_singletonfree.rds"))
# no technical replicates here; each sample one column
oral_genus <- tax_glom(oral_phy, taxrank="Genus")
og <- as.matrix(otu_table(oral_genus))
ti <- as.data.frame(as.matrix(tax_table(oral_genus)))
colnames(og) <- ti$Genus
colnames(og)[is.na(colnames(og))] <- "Unclassified"
og <- sapply(unique(colnames(og)), function(g) rowSums(og[,colnames(og)==g,drop=FALSE]))
cat("ORAL genus matrix:", nrow(og), "samples x", ncol(og), "genera\n")

# group labels for 96 oral from metadata
meta <- read.delim(file.path(base,"metadata.tsv"))
grp <- setNames(as.character(meta$Group), meta$Sample)
oral_groups <- grp[rownames(og)]
cat("Oral groups:", paste(names(table(oral_groups)), table(oral_groups), collapse=", "),"\n")

# ---------- FECAL (78 samples) ----------
seqtab <- readRDS(file.path(base,"fecal_processed/seqtab.rds"))
taxa <- readRDS(file.path(base,"fecal_processed/taxa.rds"))
rownames(seqtab) <- gsub("_F_filt\\.fastq\\.gz$","",rownames(seqtab))
run_table <- read.csv("/Volumes/edKwamiBackUP/Manuscripts/Submitted/AD.Oral/Original/raw_fastq/run_table/complete/complete_SraRunTable.csv")
fec_meta <- run_table[run_table$Run %in% rownames(seqtab) & run_table$isolation_source=="feces",]
rownames(fec_meta) <- fec_meta$Run
seqtab <- seqtab[fec_meta$Run,]
fec_phy <- phyloseq(otu_table(seqtab,taxa_are_rows=FALSE), tax_table(taxa))
fec_genus <- tax_glom(fec_phy, taxrank="Genus")
fg <- as.matrix(otu_table(fec_genus))
tif <- as.data.frame(as.matrix(tax_table(fec_genus)))
colnames(fg) <- tif$Genus
colnames(fg)[is.na(colnames(fg))] <- "Unclassified"
fg <- sapply(unique(colnames(fg)), function(g) rowSums(fg[,colnames(fg)==g,drop=FALSE]))
cat("FECAL genus matrix:", nrow(fg),"samples x",ncol(fg),"genera\n")
fec_meta$Group <- ifelse(fec_meta$disease=="AD","Case","Control")
group_labels <- fec_meta$Group
cat("Fecal groups:", paste(names(table(group_labels)),table(group_labels),collapse=", "),"\n")

# ---------- COMMON GENERA & COMBINE ----------
common <- intersect(colnames(fg), colnames(og))
cat("Common genera:", length(common),"\n")
fec_comb <- fg[,common]
oral_comb <- og[,common]
combined_ab <- rbind(fec_comb, oral_comb)
combined_meta <- data.frame(
  Sample = rownames(combined_ab),
  SampleType = c(rep("Fecal",nrow(fec_comb)), rep("Oral",nrow(oral_comb))),
  Group = c(as.character(group_labels), oral_groups), stringsAsFactors=FALSE)
rownames(combined_meta) <- combined_meta$Sample
cat("COMBINED:",nrow(combined_ab),"samples | Oral:",sum(combined_meta$SampleType=="Oral"),
    "Fecal:",sum(combined_meta$SampleType=="Fecal"),"\n")

combined_clr <- t(apply(combined_ab,1,clr_transform))
comb_pca <- prcomp(combined_clr)
ve <- round(100*comb_pca$sdev^2/sum(comb_pca$sdev^2),1)
cat("PC1:",ve[1],"% PC2:",ve[2],"%\n")

comb_dist <- dist(combined_clr, method="euclidean")
perm <- adonis2(comb_dist ~ SampleType * Group, data=combined_meta, permutations=999, by="terms")
cat("PERMANOVA sequential:\n"); print(perm)

pca_df <- data.frame(Sample=rownames(combined_clr), SampleType=combined_meta$SampleType,
                     Group=combined_meta$Group, PC1=comb_pca$x[,1], PC2=comb_pca$x[,2])
write.csv(pca_df, file.path(out,"combined_beta_pca_CORRECTED.csv"), row.names=FALSE)

# ---------- FIGURE 8 ----------
png(file.path(out,"Fig_8_combined_pca_CORRECTED.png"), width=7,height=6,units="in",res=300)
colors <- ifelse(combined_meta$SampleType=="Oral",
                 ifelse(combined_meta$Group=="Case","#E64B35","#4DBBD5"),
                 ifelse(combined_meta$Group=="Case","#F39C12","#2ECC71"))
pch <- ifelse(combined_meta$SampleType=="Oral",19,17)
plot(comb_pca$x[,1], comb_pca$x[,2], col=colors, pch=pch,
     xlab=paste0("PC1 (",ve[1],"%)"), ylab=paste0("PC2 (",ve[2],"%)"),
     main="Combined PCA (Oral vs Fecal & Case vs Control)", cex.lab=1.2, cex.axis=1.1)
legend("topright", legend=c("Oral Control","Oral AD Case","Fecal Control","Fecal AD Case"),
       col=c("#4DBBD5","#E64B35","#2ECC71","#F39C12"), pch=c(19,19,17,17))
dev.off()
cat("\nDONE\n")
