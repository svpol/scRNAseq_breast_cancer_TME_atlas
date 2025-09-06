library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)


# Objects joined_layers_rpca, joined_layers_rpca, gsea_filtered, single_gene_markers1,
# single_gene_markers2, single_gene_markers3 are defined in main.R

# Figure 2a
#---------------------------------------------------------------------------
cluster_levels_rpca <- levels(joined_layers_rpca$seurat_clusters)
cluster_colors <- c("#3A3136", "#C6121E", "#CE00CA", "#06CF22", "#2263CE",
                    "#CE8F06", "#900058", "yellow", "#708D0C", "#1EB9CF",
                    "blue", "#8A0DCE", '#D8817F')

# Make numeric cluster order
orig_clusters <- joined_layers_rpca$seurat_clusters
joined_layers_rpca$seurat_clusters <- as.numeric(orig_clusters)
colors_rpca_clusters <-  setNames(cluster_colors, cluster_levels_rpca)

umap_cluster_plot <- DimPlot(object = joined_layers_rpca, reduction = "umap", 
                             group.by = "seurat_clusters", cols = colors_rpca_clusters) + 
  theme_void() + ggtitle(NULL) + labs(color = "Cluster") +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Make cluster names as before plotting
joined_layers_rpca$seurat_clusters <- orig_clusters


# Figure 2b
#---------------------------------------------------------------------------
umap_cd45_plot <- FeaturePlot(joined_layers_rpca, features = "PTPRC", reduction = "umap", 
                              cols = c('lightgrey', '#36d900')) + theme_void() +
  theme(legend.position = "right") + ggtitle(NULL) + labs(color = "CD45\nNormalized\nExpression") +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Figure 2c
#---------------------------------------------------------------------------
cell_cycle_long_rpca_plotting <- cell_cycle_long_rpca
cell_cycle_long_rpca_plotting$seurat_clusters <- factor(
  cell_cycle_long_rpca_plotting$seurat_clusters,
  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "11", "8", "10", "12")
)

cell_cycle_plot <- ggplot(cell_cycle_long_rpca_plotting, aes(x = seurat_clusters, y = Score, fill = CellCycleScore)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ CellCycleScore, scales = "free_y", ncol = 1) +
  labs(x = "Cluster", y = "Cell Cycle Score") +
  scale_fill_manual(values = c("S Score" = "#6A1B9A", "G2M Score" = "#298c8c")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

# Figure 2d
#---------------------------------------------------------------------------
gsea_filtered$cluster <- as.numeric(gsea_filtered$cluster)

gsea_order <- c(
  "B cells(Bindea)", "Tfh cells(Bindea)", "Cytotoxic cells(Bindea)", "Temra(Schafflick)", 
  "Activated CD8 T cell(Charoentong)", "Tcell G1/S(Chung)", "Tcell G2/M(Chung)",
  "Mast cells(Bindea)", "Neutrophils(Bindea)", "Microglia1", "Microglia2", 
  "Microglia(Patir)", "Macrophage1", "Macrophage2", "Macrophages(Bindea)"   
)

gsea_filtered$Description <- factor(gsea_filtered$Description, levels = gsea_order)

gsea_plot <- ggplot(gsea_filtered, aes(x = factor(cluster), y = Description)) +
  geom_point(aes(size = -log10(p.adjust), color = NES)) +
  scale_color_gradient2(low = "darkgrey", mid = "white", high = "#003a7d", midpoint = 0) +
  labs(
    x = "Cluster",
    y = NULL,
    size = "-log10(adj.p-val)",
    color = "NES"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    plot.title = element_text(face = "bold")
  )


# Figure 2e
#---------------------------------------------------------------------------
single_gene_markers1 <- single_gene_markers1 %>%
  arrange(cell_type, gene) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

markers1_order <- c(
  "CAF", "DC", "DC (more immature)","TAM",  "NK", 
  "Neutrophil", "plasma Blast", "Bcell", "Tcell"
)

# Convert to factor in the specified order
single_gene_markers1 <- single_gene_markers1 %>%
  filter(!is.na(cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels = markers1_order))

markers1_plot <- ggplot(single_gene_markers1, aes(x = factor(cluster), y = gene)) +
  geom_point(aes(size = pct_cells_expressing, color = mean_expr_norm)) +
  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  
  facet_grid(cell_type ~ ., 
             scales = "free_y", 
             space = "free_y", 
             switch = "y") +
  # Horizontal separators between cell types
  geom_hline(
    data = data.frame(
      cell_type = unique(single_gene_markers1$cell_type),
      yintercept = Inf
    ),
    aes(yintercept = yintercept),
    color = "gray60",
    linewidth = 0.5
  ) +
  theme_minimal(base_size = 15) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.text.y = element_text(hjust = 1, colour = 'black'),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cluster") +
  guides(size = guide_legend(title = "% Cells\nExpressing"),
         color = guide_colorbar(title = "Normalized\nExpression"))

# Figure 2f
#---------------------------------------------------------------------------
single_gene_markers2 <- single_gene_markers2 %>%
  arrange(cell_type, gene) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

# Select lymphocyte clusters only
single_gene_markers2 <- single_gene_markers2[single_gene_markers2$cluster %in% c(1, 2, 4, 7, 12),]

markers2_plot <- ggplot(single_gene_markers2, aes(x = factor(cluster), y = gene)) +
  geom_point(aes(size = pct_cells_expressing, color = mean_expr_norm)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  facet_grid(cell_type ~ ., 
             scales = "free_y", 
             space = "free_y", 
             switch = "y") +  
  geom_hline(
    data = data.frame(
      cell_type = unique(single_gene_markers2$cell_type),
      yintercept = Inf
    ),
    aes(yintercept = yintercept),
    color = "gray60",
    linewidth = 0.5
  ) +  
  theme_minimal(base_size = 15) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.text.y = element_text(hjust = 1, colour = 'black'),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cluster") +
  guides(size = guide_legend(title = "% Cells\nExpressing"),
         color = guide_colorbar(title = "Normalized\nExpression"))


# Figure 2g
#---------------------------------------------------------------------------
# Select lymphocyte clusters only
single_gene_markers3 <- single_gene_markers3[single_gene_markers3$cluster %in% c(1, 2, 4, 7, 12),]

single_gene_markers3 <- single_gene_markers3 %>%
  arrange(cell_type, gene) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

markers3_plot <- ggplot(single_gene_markers3, aes(x = factor(cluster), y = gene)) +
  geom_point(aes(size = pct_cells_expressing, color = mean_expr_norm)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
  facet_grid(cell_type ~ ., 
             scales = "free_y", 
             space = "free_y", 
             switch = "y") +
  geom_hline(
    data = data.frame(
      cell_type = unique(single_gene_markers3$cell_type),
      yintercept = Inf
    ),
    aes(yintercept = yintercept),
    color = "gray60",
    linewidth = 0.5
  ) +  
  theme_minimal(base_size = 15) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.text.y = element_text(hjust = 1, colour = 'black'),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cluster") +
  guides(size = guide_legend(title = "% Cells\nExpressing"),
         color = guide_colorbar(title = "Normalized\nExpression"))


# Figure 2 grid
#---------------------------------------------------------------------------
row1_2 <- plot_grid(
  umap_cluster_plot,
  umap_cd45_plot,
  labels = c("a", "b"),
  ncol = 2,
  align = 'h',
  rel_widths = c(1, 1) 
)

row2_2 <- plot_grid(
  cell_cycle_plot,
  gsea_plot,
  labels = c("c", "d"),
  ncol = 2,
  # align = 'h',
  rel_widths = c(1, 1) 
)

left_column3_2 <- plot_grid(
  markers1_plot,
  markers2_plot,
  ncol = 1,
  labels = c("e", "f"),
  label_y = 1.02,
  label_size = 14,
  align = 'v',
  rel_heights = c(0.58, 0.42)
)

row3_2 <- plot_grid(
  left_column3_2,
  markers3_plot,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("", "g"), 
  label_size = 14
)

figure2 <- plot_grid(
  row1_2,
  row2_2,
  row3_2,
  ncol = 1,
  rel_heights = c(0.25, 0.25, 0.5)
)

# Save Figure 2
#---------------------------------------------------------------------------
ggsave(
  filename = "figures/figure2.svg",
  plot = figure2,
  width = 14,
  height = 16.5,
  unit="in",
  device = "svg",
)

