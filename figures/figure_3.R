library(Seurat)
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)


# Objects joined_layers_rpca, leukocyte_summary and leukocyte_data are defined in main.R

# Figure 3a
#---------------------------------------------------------------------------
leuko_cluster_plot <- DimPlot(joined_layers_rpca, group.by = "Leukocyte_cluster", 
                              reduction = "umap") + ggtitle(NULL) +
  theme_void() + theme(legend.position = "right", legend.text = element_text(size = 12)) + 
  scale_color_manual(
    values = c("FALSE" = "lightgrey", "TRUE" = "#3F7CAC"),
    labels = c("FALSE" = "Non-Leukocyte", "TRUE" = "Leukocyte")
  )


# Figure 3b
#---------------------------------------------------------------------------
leuko_pro_dataset <- ggplot(leukocyte_summary, aes(x = dataset_id, y = proportion_leukocytes * 100)) +
  geom_bar(stat = "identity", fill = "#3F7CAC") +
  geom_text(aes(y = proportion_leukocytes * 100, label = total_cells),
            vjust = -0.5, size = 4) +
  geom_text(aes(y = (proportion_leukocytes * 100) / 2,
                label = sprintf("%.2f%%", proportion_leukocytes * 100)),
            color = "white", size = 4) +
  scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45),
    axis.text.y = element_text(size = 12)
  )

# Figure 3c
#---------------------------------------------------------------------------
leuko_by_datset_plot <- ggplot(leukocyte_data, 
                               aes(x = prop, 
                                   y = dataset_id,
                                   fill = cluster_annotation_main)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(
    values = rev(brewer_pal(palette = "Set3")(7)),
    breaks = rev(cell_type_order),  
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))


# Figure 3 grid
#---------------------------------------------------------------------------
row1_3 <- plot_grid(
  leuko_cluster_plot, 
  leuko_pro_dataset,
  labels = c("a", "b"),
  ncol = 2,
  align = 'h',
  rel_widths = c(1, 1) 
)

row2_3 <- plot_grid(
  leuko_by_datset_plot,
  labels = c("c"),
  ncol = 1,
  align = 'h'
)


figure3 <- plot_grid(
  row1_3,
  row2_3,
  ncol = 1,
  rel_heights = c(1,1)
)


# Save Figure 3
ggsave(
  filename = "figures/figure3.svg",
  plot = figure3,
  width = 12,
  height = 10,
  device = 'svg'
)
