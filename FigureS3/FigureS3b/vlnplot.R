setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs")

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/DRG_Neuron_aging_2_24_final.Rdata")

DRG_Neuron_aging$age <- factor(DRG_Neuron_aging$age, levels = c("3 MO","24 MO"))
DefaultAssay(DRG_Neuron_aging) <- "RNA"

age_colors <- c("3 MO" = "#A2BF47", "24 MO" = "#B1672A")

# 选择要展示的5个基因
genes_to_plot <- c("Cdkn1a", "Cdkn2a")

VlnPlot(DRG_Neuron_aging, features = "Cdkn1a",group.by = "Celltype_2_2",split.by = "age",
        cols = age_colors,pt.size = 0,combine = FALSE)

plots <- list()
for (gene in genes_to_plot) {
  p <- VlnPlot(
    DRG_Neuron_aging,
    features = gene,
    group.by = "Celltype_2_2",
    split.by = "age",
    cols = age_colors,
    pt.size = 0,  
    combine = FALSE
  )[[1]] +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = gene)
  
  plots[[gene]] <- p
}

# 使用patchwork包组合所有图形
library(patchwork)
combined_plot <- wrap_plots(plots, ncol = 1) + 
  plot_annotation(title = "Gene Expression Across Cell Types and Ages")

# 显示图形
print(combined_plot)

# 保存图形
ggsave("Cdkn1a_Cdkn2a_vlnplot.pdf", combined_plot, width = 16, height = 10) 


# 选择要展示的5个基因
genes_to_plot <- c("Gabrg3","Trpm8", "Mrgpra3","Mrgprd")

# VlnPlot(DRG_Neuron_aging, features = "Cdkn1a",group.by = "Celltype_2_2",split.by = "age",
#         cols = age_colors,pt.size = 0,combine = FALSE)

plots <- list()
for (gene in genes_to_plot) {
  p <- VlnPlot(
    DRG_Neuron_aging,
    features = gene,
    group.by = "Celltype_2_2",
    split.by = "age",
    cols = age_colors,
    pt.size = 0,  
    combine = FALSE
  )[[1]] +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = gene)
  
  plots[[gene]] <- p
}

# 使用patchwork包组合所有图形
library(patchwork)
combined_plot <- wrap_plots(plots, ncol = 1) + 
  plot_annotation(title = "Gene Expression Across Cell Types and Ages")

# 显示图形
print(combined_plot)

# 保存图形
ggsave("trpm8_mrgpra3_mrgprd_vlnplot.pdf", combined_plot, width = 16, height = 10) 


# 选择要展示的5个基因
genes_to_plot <- c("Ccr2", "Ccr3","Ccr5")
load("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")
DRG_aging_all$age <- factor(DRG_aging_all$age, levels = c("3 MO","24 MO"))
DefaultAssay(DRG_aging_all) <- "RNA"

plots <- list()
for (gene in genes_to_plot) {
  p <- VlnPlot(
    DRG_aging_all,
    features = gene,
    group.by = "Celltype_3",
    split.by = "age",
    cols = age_colors,
    pt.size = 0,  
    combine = FALSE
  )[[1]] +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = gene)
  
  plots[[gene]] <- p
}

# 使用patchwork包组合所有图形
library(patchwork)
combined_plot <- wrap_plots(plots, ncol = 1) + 
  plot_annotation(title = "Gene Expression Across Cell Types and Ages")

# 显示图形
print(combined_plot)

# 保存图形
ggsave("Ccr2_Ccr3_Ccr5_vlnplot.pdf", combined_plot, width = 12, height = 8) 
