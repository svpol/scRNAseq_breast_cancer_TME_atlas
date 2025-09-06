library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(cowplot)
library(rstatix)


# Objects leuko_prop_by_subtype, cell_type_order, leuko_tot, heatmap_data, joined_layers_rpca, 
# patient_fractions and cell_cycle_clinical_patient are defined in main.R


# Figure 4a
#---------------------------------------------------------------------------
leuko_prop_plot <- ggplot(leuko_prop_by_subtype, aes(x = clinical_subtype, y = proportion, fill = leukocyte_status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = scales::percent(proportion, accuracy = 0.01)),  
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 3.5
  ) +
  scale_fill_manual(values = c("#3F7CAC", "#D8BCAB")) +
  labs(
    x = "Clinical Subtype",
    y = NULL,
    fill = NULL
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + theme_minimal() + 
  theme(legend.position = "right",
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 11))


# Figure 4b
#---------------------------------------------------------------------------
original_colors <- brewer_pal(palette = "Set3")(7)
names(original_colors) <- cell_type_order

custom_colors <- c(
  original_colors,
  "All T cells" = "#1b9e77",
  "All B cells" = "#9E1B90"
)

leuko_types_plot <- ggplot(leuko_tot, aes(x = x_pos, y = proportion, fill = cell_group)) +
  
  geom_col(width = 0.35, position = position_stack(reverse = FALSE)) +
  
  geom_text(
    aes(label = ifelse(proportion >= 0.005,
                       scales::percent(proportion, accuracy = 0.01),
                       "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    color = "black",
    size = 3.5
  ) +
  
  scale_x_continuous(
    breaks = 1:length(clinical_subtype_order),
    labels = clinical_subtype_order
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = custom_colors, breaks = stack_order, 
                    drop = FALSE, name = "Cell Type") +
  labs(
    x = "Clinical Subtype",
    y = NULL,
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )


# Figure 4c
#---------------------------------------------------------------------------
clinical_subtype_boxplot <- ggplot(patient_fractions,
                                   aes(x = cluster_annotation_main, y = proportion, fill = clinical_subtype)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(aes(color = clinical_subtype), size = 2, alpha = 0.4,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_fill_manual(values = c("TNBC" = "#E64B35", "HER2BC" = "#4DBBD5", "LBC" = "#00A087")) +
  scale_color_manual(values = c("TNBC" = "#E64B35", "HER2BC" = "#4DBBD5", "LBC" = "#00A087")) +
  labs(x = NULL, y = NULL, fill = "Clinical Subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none") +
  stat_pvalue_manual(stat_tests_detailed, 
                     label = "p.adj.signif", 
                     tip.length = 0.02, 
                     y.position = 1, 
                     hide.ns = TRUE)


# Figure 4d
#---------------------------------------------------------------------------
cell_cycle_clinical_plot <- ggplot(cell_cycle_clinical_patient,
                                   aes(x = Score_Type, y = Score, fill = clinical_subtype)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(color = clinical_subtype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.3) +
  labs(x = NULL, y = "Cell Cycle Score", fill = "Clinical Subtype", color = "Clinical Subtype") +
  theme_minimal(base_size = 14) +
  stat_pvalue_manual(stat_tests_cc, label = "p.adj.signif", tip.length = 0.02, hide.ns = TRUE) +
  scale_fill_manual(values = c("TNBC" = "#E64B35", "HER2BC" = "#4DBBD5", "LBC" = "#00A087")) +
  scale_color_manual(values = c("TNBC" = "#E64B35", "HER2BC" = "#4DBBD5", "LBC" = "#00A087"))


# Figure 4e
#---------------------------------------------------------------------------
heatmap_matrix <- as.matrix(heatmap_data[,-1])
# Prepare clinical subtype annotation
sample_annot <- joined_layers_rpca@meta.data %>%
  select(sample_id, clinical_subtype) %>%
  distinct() %>%
  mutate(clinical_subtype = factor(clinical_subtype, 
                                   levels = c("TNBC", "HER2BC", "LBC")))
subtype_colors <- c("TNBC" = "#E64B35",
                    "HER2BC" = "#4DBBD5",
                    "LBC" = "#00a987")    

# Create heatmap annotation
ha_column <- HeatmapAnnotation(
  Clinical_Subtype = sample_annot$clinical_subtype[
    match(colnames(heatmap_matrix), sample_annot$sample_id)
  ],
  col = list(Clinical_Subtype = subtype_colors),
  annotation_name_side = "left",
  annotation_legend_param = list(
    title = "Clinical Subtype",
    at = names(subtype_colors)
  )
)

ht <- Heatmap(
  heatmap_matrix,
  name = "Proportion",
  col = colorRamp2(c(0, max(heatmap_matrix)), c("white", "steelblue")),
  row_order = rev(cell_type_order),
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  row_names_side = "left",
  show_column_names = F,
  row_names_gp = gpar(fontsize = 12),
  top_annotation = ha_column,
  heatmap_legend_param = list(
    title = "Leukocyte\nProportion",
    title_position = "leftcenter-rot"
  )
)

ht_grob <- grid.grabExpr(draw(ht))

# Figure grid
#---------------------------------------------------------------------------
row1_4 <- plot_grid(
  leuko_prop_plot,
  leuko_types_plot,
  labels = c("a", "b"),
  ncol = 2,
  align = 'h',
  rel_widths = c(0.33, 0.67)
)

row2_4 <- plot_grid(
  clinical_subtype_boxplot,
  cell_cycle_clinical_plot,
  labels = c("c", "d"),
  ncol = 2,
  align = 'h',
  rel_widths = c(0.65, 0.35)
)

row3_4 <- plot_grid(
  ht_grob,
  labels = "e",
  ncol = 1,
  align = 'h',
  rel_widths = 1
)

figure4 <- plot_grid(
  row1_4,
  row2_4,
  row3_4,
  ncol = 1,
  rel_heights = c(0.3, 0.3, 0.4)
)


# Figure save
#---------------------------------------------------------------------------
ggsave(
  filename = "figures/figure4.svg",
  plot = figure4,
  width = 12,height = 14,
  device = 'svg'
)
