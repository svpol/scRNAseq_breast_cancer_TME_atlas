library(magick)
library(grid)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# Figure S1a was previously drawn manually.
# Objects merged_all, integrated_all_rpca, integrated_all_harmony, seurat_list_2sd and cell_count_steps are defined in main.R

# Figure S1 a
#---------------------------------------------------------------------------
methods_plot <- image_read("method_steps.png")
methods_plot <- rasterGrob(as.raster(methods_plot), interpolate = TRUE)

# Figure S1 b
#---------------------------------------------------------------------------
# Variance explained - merged
pct_merged <- merged_all@reductions$pca@stdev / sum(merged_all@reductions$pca@stdev) * 100
elbow_df_merged <- data.frame(
  PC = 1:length(pct_merged),
  Variance = pct_merged
)

elbow_merged <- ggplot(elbow_df_merged, aes(x = PC, y = Variance)) +
  geom_point(color = "black") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Merged",
    x = "PC", y = "Variance Explained") + theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# Variance explained - RPCA
pc_matrix_rpca <- Embeddings(integrated_all_rpca, reduction = "rpca_integrated")
pc_stdev_rpca <- apply(pc_matrix_rpca, 2, sd)
pct_var_rpca <- pc_stdev_rpca^2 / sum(pc_stdev_rpca^2) * 100
elbow_df_rpca <- data.frame(
  PC = 1:length(pct_var_rpca),
  Variance = pct_var_rpca
)

elbow_rpca <- ggplot(elbow_df_rpca, aes(x = PC, y = Variance)) +
  geom_point(color = "black") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "RPCA-Integrated",
    x = "PC", y = "Variance Explained") + theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# Variance explained - Harmony
pc_matrix_harmony <- Embeddings(integrated_all_harmony, reduction = "harmony")
pc_stdev_harmony <- apply(pc_matrix_harmony, 2, sd)
pct_var_harmony <- pc_stdev_harmony^2 / sum(pc_stdev_harmony^2) * 100
elbow_df_harmony <- data.frame(
  PC = 1:length(pct_var_harmony),
  Variance = pct_var_harmony
)

elbow_harmony <- ggplot(elbow_df_harmony, aes(x = PC, y = Variance)) +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "black") +
  labs(
    title = "Harmony-Integrated",
    x = "PC", y = "Variance Explained") + theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


# Figure S1 c
#---------------------------------------------------------------------------
dataset_order <- c("GSE161529", "LambrechtPD1", "GSE176078", 
                   "GSE167036", "GSE114727", "GSE148673", "GSE110686")

# Reshape to long format
cell_count_long <- cell_count_steps %>%
  pivot_longer(
    cols = c(initial_cells, cells_min_qc, cells_mito_qc, cells_max_qc, cells_2sd_qc),
    names_to = "qc_step",
    values_to = "cell_count"
  )

# Set QC step order
qc_order <- c("initial_cells", "cells_min_qc", "cells_mito_qc", "cells_max_qc", "cells_2sd_qc")
cell_count_long$qc_step <- factor(cell_count_long$qc_step, levels = qc_order)

# Set dataset order
cell_count_long$dataset_name <- factor(cell_count_long$dataset_name, levels = dataset_order)

# Calculate summary stats
total_before_qc <- sum(cell_count_steps$initial_cells)
total_after_qc <- sum(cell_count_steps$cells_2sd_qc)
avg_after_qc <- mean(cell_count_steps$cells_2sd_qc)

summary_text <- paste0(
  "Total cells before QC: ", format(total_before_qc),
  "\nTotal cells after QC: ", format(total_after_qc),
  "\nAvg cells per dataset after QC: ", round(avg_after_qc, 0)
)

qc_cell_plot <- ggplot(cell_count_long, aes(x = dataset_name, y = cell_count, fill = qc_step)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = cell_count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 4, angle = 45) +
  annotate("text", 
           x = length(dataset_order) + 0.5,
           y = max(cell_count_long$cell_count) * 1.05, 
           label = summary_text, hjust = 1, vjust = 1, size = 4) +
  scale_fill_brewer(palette = "Set2", name = "QC step") +
  labs(x = "Dataset", y = "Number of cells") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


# Figure S1 d
#---------------------------------------------------------------------------
gene_counts_df <- do.call(rbind, lapply(names(seurat_list_2sd), function(nm) {
  data.frame(dataset = nm, n_genes = seurat_list_2sd[[nm]]$nFeature_RNA)
}))

avg_genes_per_dataset <- gene_counts_df %>%
  group_by(dataset) %>%
  summarise(avg_genes = mean(n_genes), .groups = "drop")

# Calculate overall average (mean of dataset means)
overall_avg <- mean(avg_genes_per_dataset$avg_genes)
overall_avg_fmt <- format(round(overall_avg, 0), nsmall = 0)

# Create annotation text
annotation_text <- paste("Average genes expressed by dataset:", overall_avg_fmt)

qc_gene_plot <- ggplot(gene_counts_df, aes(x = dataset, y = n_genes)) +
  geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = dataset), width = 0.15, alpha = 0.02, size = 0.2) +
  annotate("text", 
           x = length(dataset_order) + 0.5, 
           y = max(gene_counts_df$n_genes) * 1.02, 
           label = annotation_text, 
           hjust = 1, vjust = 1, size = 4) +
  labs(x = "Dataset", y = "Number of genes expressed") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )


# Figure S1 grid
#---------------------------------------------------------------------------

elbow_plots <- plot_grid(
  elbow_merged,
  elbow_rpca,
  elbow_harmony,
  ncol = 3,
  labels = c("b", "", "")
)

row1_sup1 <- plot_grid(
  methods_plot,
  elbow_plots,
  ncol = 2,
  labels = c("a", ""),
  label_size = 14,
  rel_widths = c(0.3,0.7)
)

row2_sup1 <- plot_grid(
  qc_cell_plot,
  ncol = 1,
  labels = "c",
  label_size = 14
)

row3_sup1 <- plot_grid(
  qc_gene_plot,
  ncol = 1,
  labels = "d"
)

sup1 <- plot_grid(
  row1_sup1,
  row2_sup1,
  row3_sup1,
  ncol = 1,
  rel_heights = c(0.8, 1, 0.7)  
)


# Save Figure S1
#---------------------------------------------------------------------------
ggsave(filename = "figures/figure_sup1.png", 
       plot = sup1, 
       width = 12, 
       height = 14, 
       dpi=500)
