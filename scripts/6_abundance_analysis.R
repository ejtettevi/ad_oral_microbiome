# comprehensive_oral_microbiome_analysis.R
suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(pheatmap)
  library(pwr)
  library(viridis)
  library(ANCOMBC)   # For compositional DA
  library(ALDEx2)    # For compositional DA
})

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

asv_table_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_table.tsv"
taxonomy_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_taxonomy.tsv"
metadata_path <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"
output_dir <- "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/abundance"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
outpath <- function(filename) file.path(output_dir, filename)

# =============================================================================
# 1. DATA LOADING AND PREPARATION
# =============================================================================

cat("Checking input files...\n")
if (!file.exists(asv_table_path)) stop("ASV table file not found: ", asv_table_path)
if (!file.exists(taxonomy_path)) stop("Taxonomy file not found: ", taxonomy_path)
if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)

cat("Loading ASV table...\n")
asv_table <- read.table(asv_table_path, sep = "\t", header = TRUE, row.names = 1, 
                        check.names = FALSE, stringsAsFactors = FALSE)

# --- Automatically transpose if samples are in rows and ASVs are in columns ---
if (grepl("^SRR", rownames(asv_table)[1])) {
  cat("Transposing ASV table: detected samples in rows and ASVs in columns.\n")
  asv_table <- t(as.matrix(asv_table))
}

cat("Loading taxonomy table...\n")
taxonomy <- read.table(taxonomy_path, sep = "\t", header = TRUE, row.names = 1, 
                       check.names = FALSE, stringsAsFactors = FALSE)

# --- Clean taxonomy table ---
taxonomy[is.na(taxonomy) | taxonomy == "" | taxonomy == " "] <- "Unknown"
col_mapping <- list(
  "kingdom" = "Kingdom", "k__" = "Kingdom",
  "phylum" = "Phylum", "p__" = "Phylum",
  "class" = "Class", "c__" = "Class",
  "order" = "Order", "o__" = "Order",
  "family" = "Family", "f__" = "Family",
  "genus" = "Genus", "g__" = "Genus",
  "species" = "Species", "s__" = "Species"
)
for (old_name in names(col_mapping)) {
  if (old_name %in% colnames(taxonomy)) {
    colnames(taxonomy)[colnames(taxonomy) == old_name] <- col_mapping[[old_name]]
  }
}

cat("Loading metadata...\n")
metadata <- read.table(metadata_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL

# --- Automated sample matching (robust) ---
clean_names <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  return(x)
}

asv_names_raw <- colnames(asv_table)
meta_names_raw <- rownames(metadata)

asv_names_clean <- clean_names(asv_names_raw)
meta_names_clean <- clean_names(meta_names_raw)

common_samples_clean <- intersect(asv_names_clean, meta_names_clean)

if (length(common_samples_clean) == 0) {
  cat("No matching sample names between ASV table and metadata!\n")
  cat("ASV table sample names (cleaned):\n")
  print(asv_names_clean)
  cat("Metadata sample names (cleaned):\n")
  print(meta_names_clean)
  stop("Please ensure sample names match between ASV table and metadata.")
}

asv_match_idx <- which(asv_names_clean %in% common_samples_clean)
meta_match_idx <- which(meta_names_clean %in% common_samples_clean)
asv_match_names <- asv_names_raw[asv_match_idx]
meta_match_names <- meta_names_raw[meta_match_idx]

asv_table <- asv_table[, asv_match_names, drop = FALSE]
metadata <- metadata[meta_match_names, , drop = FALSE]

dropped_asv <- setdiff(asv_names_raw, asv_match_names)
dropped_meta <- setdiff(meta_names_raw, meta_match_names)
if (length(dropped_asv) > 0 || length(dropped_meta) > 0) {
  cat("Warning: The following samples were dropped due to mismatches:\n")
  if (length(dropped_asv) > 0) cat("ASV table:", paste(dropped_asv, collapse = ", "), "\n")
  if (length(dropped_meta) > 0) cat("Metadata:", paste(dropped_meta, collapse = ", "), "\n")
} else {
  cat("All sample names matched between ASV table and metadata.\n")
}

# --- Automated ASV/taxonomy matching ---
common_asvs <- intersect(rownames(asv_table), rownames(taxonomy))
if (length(common_asvs) == 0) {
  cat("No matching ASV names between ASV table and taxonomy!\n")
  cat("ASV table ASV names:\n")
  print(rownames(asv_table))
  cat("Taxonomy ASV names:\n")
  print(rownames(taxonomy))
  stop("Please ensure ASV names match between ASV table and taxonomy.")
}
asv_table <- asv_table[common_asvs, , drop = FALSE]
taxonomy <- taxonomy[common_asvs, , drop = FALSE]

dropped_asvs_table <- setdiff(rownames(asv_table), common_asvs)
dropped_asvs_taxonomy <- setdiff(rownames(taxonomy), common_asvs)
if (length(dropped_asvs_table) > 0 || length(dropped_asvs_taxonomy) > 0) {
  cat("Warning: The following ASVs were dropped due to mismatches:\n")
  if (length(dropped_asvs_table) > 0) cat("ASV table:", paste(dropped_asvs_table, collapse = ", "), "\n")
  if (length(dropped_asvs_taxonomy) > 0) cat("Taxonomy:", paste(dropped_asvs_taxonomy, collapse = ", "), "\n")
} else {
  cat("All ASV names matched between ASV table and taxonomy.\n")
}

# --- Create phyloseq object ---
cat("Creating phyloseq object...\n")
ps <- phyloseq(
  otu_table(as.matrix(asv_table), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)
)

cat("Phyloseq object created with", nsamples(ps), "samples and", ntaxa(ps), "ASVs\n")

# --- Automated group and batch variable detection ---
possible_group_vars <- c("Group", "group", "Treatment", "treatment", "Condition", "condition", 
                        "Type", "type", "Status", "status", "Category", "category")
group_var <- NULL
for (var in possible_group_vars) {
  if (var %in% colnames(sample_data(ps))) {
    group_var <- var
    break
  }
}
if (is.null(group_var)) {
  group_var <- colnames(sample_data(ps))[1]
  cat("No standard group variable found. Using first column:", group_var, "\n")
} else {
  cat("Using group variable:", group_var, "\n")
}

possible_batch_vars <- c("Batch", "batch", "Plate", "plate", "Run", "run", "SeqBatch", "seqbatch")
batch_var <- NULL
for (var in possible_batch_vars) {
  if (var %in% colnames(sample_data(ps))) {
    batch_var <- var
    break
  }
}
if (!is.null(batch_var)) {
  cat("Detected batch variable:", batch_var, "\n")
}

# --- Filter samples and ASVs ---
min_reads <- 1000
ps <- prune_samples(sample_sums(ps) >= min_reads, ps)
min_abundance <- 10
ps <- prune_taxa(taxa_sums(ps) >= min_abundance, ps)
cat("After filtering:", nsamples(ps), "samples and", ntaxa(ps), "ASVs remain\n")

# =============================================================================
# 2. TAXONOMIC COMPOSITION ANALYSIS (Improved for visibility and accessibility)
# =============================================================================

discriminable_palette <- function(n) {
  # Always use viridis for color-blind accessibility
  viridis::viridis(n, option = "D", end = 0.95)
}

create_abundance_plots <- function(ps, group_var = "Group", top_n = 10) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  tax_levels <- intersect(c("Phylum", "Class", "Order", "Family", "Genus"), colnames(tax_table(ps_rel)))
  plot_list <- list()
  for (level in tax_levels) {
    ps_level <- tax_glom(ps_rel, taxrank = level, NArm = TRUE)
    top_taxa <- names(sort(taxa_sums(ps_level), decreasing = TRUE))[1:min(top_n, ntaxa(ps_level))]
    ps_top <- prune_taxa(top_taxa, ps_level)
    plot_data <- psmelt(ps_top)
    plot_data[[level]] <- as.character(plot_data[[level]])
    plot_data[[level]][plot_data[[level]] == "Unknown" | is.na(plot_data[[level]]) | plot_data[[level]] == ""] <- paste("Unknown", level)
    n_colors <- length(unique(plot_data[[level]]))
    p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = .data[[level]])) +
      geom_bar(stat = "identity") +
      facet_wrap(as.formula(paste("~", group_var)), scales = "free_x", nrow = 1) +
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold")
      ) +
      labs(
        title = paste("Relative Abundance -", level, "Level"),
        y = "Relative Abundance",
        x = "Samples"
      ) +
      scale_fill_manual(values = discriminable_palette(n_colors)) +
      guides(fill = guide_legend(ncol = 3, title = level))
    plot_list[[level]] <- p
    avg_data <- plot_data %>%
      group_by(.data[[group_var]], .data[[level]]) %>%
      summarise(Mean_Abundance = mean(Abundance), .groups = 'drop') %>%
      mutate(Group = .data[[group_var]], Taxon = .data[[level]])
    n_colors_avg <- length(unique(avg_data$Taxon))
    p_avg <- ggplot(avg_data, aes(x = Group, y = Mean_Abundance, fill = Taxon)) +
      geom_bar(stat = "identity") +
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      ) +
      labs(
        title = paste("Average Relative Abundance -", level, "Level"),
        y = "Mean Relative Abundance",
        x = "Group"
      ) +
      scale_fill_manual(values = discriminable_palette(n_colors_avg)) +
      guides(fill = guide_legend(ncol = 3, title = level))
    plot_list[[paste0(level, "_avg")]] <- p_avg
  }
  return(plot_list)
}

cat("Generating taxonomic abundance plots...\n")
abundance_plots <- create_abundance_plots(ps, group_var = group_var, top_n = 15)
for (plot_name in names(abundance_plots)) {
  tryCatch({
    ggsave(
      outpath(paste0(plot_name, "_abundance.png")),
      abundance_plots[[plot_name]],
      width = 14, height = 9, dpi = 700
    )
    cat("Saved:", outpath(paste0(plot_name, "_abundance.png")), "\n")
  }, error = function(e) {
    cat("Error creating plot", plot_name, ":", e$message, "\n")
  })
}

# =============================================================================
# 3. COMPOSITIONAL DIFFERENTIAL ABUNDANCE ANALYSIS
# =============================================================================

cat("Running compositional differential abundance analysis (ANCOM-BC)...\n")
ancombc_res <- tryCatch({
  ancombc_out <- ancombc(phyloseq = ps, formula = paste(group_var, if(!is.null(batch_var)) paste("+", batch_var) else ""), p_adj_method = "fdr")
  ancombc_df <- ancombc_out$res$diff_abn
  write.table(ancombc_df, outpath("ancombc_results.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
  cat("ANCOM-BC results saved.\n")
  ancombc_df
}, error = function(e) {
  cat("ANCOM-BC failed:", e$message, "\n")
  NULL
})

cat("Running compositional differential abundance analysis (ALDEx2)...\n")
aldex2_res <- tryCatch({
  # Prepare data for ALDEx2
  asv_counts <- as.data.frame(otu_table(ps))
  group_factor <- as.factor(sample_data(ps)[[group_var]])
  aldex2_out <- aldex(asv_counts, group_factor, test = "t", effect = TRUE, denom = "all", mc.samples = 128)
  write.table(aldex2_out, outpath("aldex2_results.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
  cat("ALDEx2 results saved.\n")
  aldex2_out
}, error = function(e) {
  cat("ALDEx2 failed:", e$message, "\n")
  NULL
})

# =============================================================================
# 4. ADDITIONAL PLOTS AND SESSION INFO (Improved heatmap)
# =============================================================================

cat("Generating heatmap of top ASVs...\n")
top_n_asvs <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:min(30, ntaxa(ps))]
ps_top_asvs <- prune_taxa(top_n_asvs, ps)
asv_mat <- otu_table(ps_top_asvs)
asv_mat <- as.matrix(asv_mat)
rownames(asv_mat) <- paste0("ASV_", rownames(asv_mat))
annotation_col <- data.frame(Group = sample_data(ps_top_asvs)[[group_var]])
rownames(annotation_col) <- sample_names(ps_top_asvs)
pheatmap(
  asv_mat,
  annotation_col = annotation_col,
  color = viridis::viridis(100, option = "D"),
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize_row = 12, fontsize_col = 12,
  main = "Heatmap of Top 30 ASVs",
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  filename = outpath("heatmap_top30_asvs.png"),
  width = 12, height = 10, dpi = 700
)
cat("Saved:", outpath("heatmap_top30_asvs.png"), "\n")

cat("Analysis completed successfully!\n")
cat("All results have been saved to:", output_dir, "\n")
cat("Session info:\n")
writeLines(capture.output(sessionInfo()), outpath("sessionInfo.txt"))
