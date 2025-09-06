library(ggplot2)
library(scales)
library(cowplot)

# Objects merged_counts_dc, rpca_counts_dc and harmony_counts_dc are defined in main.R

# Merged object
#---------------------------------------------------------------------------
merged_prop_heatmap <- ggplot(merged_counts_dc, aes(x = dataset_id, y = seurat_clusters, fill = proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.01)), size = 3.6) +
  scale_fill_gradient(low = "white", high = "darkgreen", labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "Merged",
       x = "Dataset", y = "Cluster", fill = "Percentage")


# RPCA
#---------------------------------------------------------------------------
rpca_prop_heatmap <- ggplot(rpca_counts_dc, aes(x = dataset_id, y = seurat_clusters, fill = proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.01)), size = 3.6) +
  scale_fill_gradient(low = "white", high = "darkgreen", labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "RPCA-Integrated",
       x = "Dataset", y = "Cluster", fill = "Percentage")


# Harmony
#---------------------------------------------------------------------------
harmony_prop_heatmap <- ggplot(harmony_counts_dc, aes(x = dataset_id, y = seurat_clusters, fill = proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.01)), size = 3.6) +
  scale_fill_gradient(low = "white", high = "darkgreen", labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "Harmony-Integrated",
       x = "Dataset", y = "Cluster", fill = "Percentage")


# Figure S2 grid
#---------------------------------------------------------------------------
figure2_sup <- plot_grid(
  merged_prop_heatmap,
  rpca_prop_heatmap,
  harmony_prop_heatmap,
  ncol = 1,
  rel_heights = c(1,1,1)
  # labels = c("a", "b", "c")
)

# Save Figure S2
#---------------------------------------------------------------------------
ggsave(
  filename = "figures/figure2_sup.png",
  plot = figure2_sup,
  width = 9,
  height = 11,
  dpi = 500
)
