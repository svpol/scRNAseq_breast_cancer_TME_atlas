library(Seurat)
library(ggplot2)
library(cowplot)
library(ggalluvial)
library(grid)

# Objects merged_all, integrated_all_rpca, integrated_all_harmony, cluster_dataset_merged_df, 
# cluster_dataset_rpca_df, cluster_dataset_harmony_df, entropy_rpca, entropy_harmony, 
# sankey_all are defined in main.R


# Manually defined colors for better distinction
dataset_colors <- c('blue', 'red', 'yellow', '#FF8DA1', '#4daf4a', '#800080', 'black')
cluster_colors <- c("#3A3136", "#C6121E", "#CE00CA", "#06CF22", "#2263CE", "#CE8F06",
                    "#900058", "yellow", "#708D0C", "#1EB9CF", "blue", "#8A0DCE", '#D8817F')


# Figure 1a
#---------------------------------------------------------------------------
dimplot1 <- DimPlot(merged_all, reduction = "umap", group.by = "dataset_id", 
                    cols = dataset_colors) + theme_void() + NoLegend() + ggtitle('Merged') +
  theme(plot.title = element_text(hjust = 0.1, face = "plain"))
dimplot2 <- DimPlot(integrated_all_rpca, reduction = "umap", group.by = "dataset_id", 
                    cols = dataset_colors) + theme_void() + NoLegend() + ggtitle('RPCA-Integrated') 
dimplot3 <- DimPlot(integrated_all_harmony, reduction = "umap", group.by = "dataset_id", 
                    cols = dataset_colors) + theme_void() + NoLegend() + ggtitle('Harmony-Integrated') 

legend_dimplot <- get_legend(
  DimPlot(merged_all, 
          reduction = "umap", 
          group.by = "dataset_id", 
          cols = dataset_colors) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8.5, margin = margin(l = 1)),
      legend.key.width = unit(0.4, "cm")
    )
)

# Figure 1b
#---------------------------------------------------------------------------
dataset_barplot1 <- ggplot(cluster_dataset_rpca_df, aes(x = Cluster, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity") +
  ylab("Percentage") +
  xlab("Cluster") +
  ggtitle("RPCA") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0)
  ) + NoLegend() +
  scale_fill_manual(values = dataset_colors)

dataset_barplot2 <- ggplot(cluster_dataset_harmony_df, aes(x = Cluster, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity") +
  ylab("Percentage") +
  xlab("Cluster") +
  ggtitle("Harmony") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0)
  ) + NoLegend() +
  scale_fill_manual(values = dataset_colors)

# Figure 1c
#---------------------------------------------------------------------------
entplot_rpca <- ggplot(entropy_rpca, aes(x = as.factor(seurat_clusters), y = entropy)) +
  geom_bar(stat = "identity", fill = "#222222") +
  xlab("Cluster") + ylab("Shannon Entropy") +
  ggtitle("RPCA") + theme_minimal() + NoLegend() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

entplot_harmony <- ggplot(entropy_harmony, aes(x = as.factor(seurat_clusters), y = entropy)) +
  geom_bar(stat = "identity", fill = "#222222") +
  xlab("Cluster") + ylab("Shannon Entropy") +
  ggtitle("Harmony") + theme_minimal() + NoLegend() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

# Figure 1d
#---------------------------------------------------------------------------
sankey_plot <- ggplot(
  sankey_all,
  aes(axis1 = seurat_clusters_rpca, axis2 = seurat_clusters_harmony, y = count)
) +
  geom_alluvium(aes(fill = seurat_clusters_rpca), width = 0.05) + 
  geom_stratum(
    stat = "stratum",
    width = 0.05,                 
    fill = "white",
    color = "black"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3,
    hjust = 0.5
  ) +
  scale_x_discrete(
    limits = c("RPCA", "Harmony"),
    expand = c(.1, .1)
  ) +
  scale_fill_manual(
    values = cluster_colors,
    name = "RPCA Cluster"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = 1, y = max(sankey_all$count) * -0.15,
           label = "RPCA Clusters", size = 5) +
  annotate("text", x = 2, y = max(sankey_all$count) * -0.15,
           label = "Harmony Clusters", size = 5)


# Figure 1 grid
#---------------------------------------------------------------------------
row1_1 <- plot_grid(dimplot1, dimplot2, dimplot3, legend_dimplot, ncol = 4,
                    rel_widths = c(0.3, 0.3, 0.3, 0.1))

row2_1 <- plot_grid(dataset_barplot1, dataset_barplot2, entplot_rpca, entplot_harmony, ncol = 4)

row2_labeled <- ggdraw() +
  draw_plot(row2_1, 0, 0, 1, 1) +
  draw_plot_label(c("b", "c"), 
                  x = c(0.02, 0.52), y = c(1, 1), size = 14)

row3_1 <- plot_grid(
  sankey_plot,
  ncol = 1
)


final_plot_1 <- plot_grid(
  row1_1,
  row2_labeled,
  row3_1,
  ncol = 1,
  rel_heights = c(0.32, 0.32, 0.36),
  labels = c("a", "", "d") 
)

# Save Figure 1
#---------------------------------------------------------------------------
ggsave(
  filename = "figures/figure1.svg",
  plot = final_plot_1,
  width = 12,
  height = 13,
  device = 'svg'
)
