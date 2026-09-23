getwd()

cellInfo_2 <- as.data.frame(cellInfo[,"CellType"]) 
rownames(cellInfo_2) <- rownames(cellInfo)
colnames(cellInfo_2) <- "CellType"

library(Seurat)
load("D:/paper_Wang&Chen2025/Rdata/Fig4&5/DRG_Macrophage_aging_2_24.Rdata")
FeaturePlot(DRG_Macrophage_aging,features = c("Nfatc3","Max","Nr1h2","Maf","Zfp384","Ets2","Lgals3","Gngt2","Soat1","Gtf2f1","Dab2","Slc9a9"))


library(ggplot2)
library(cowplot)

library(patchwork)
value_columns <- c("Nfatc3","Max","Nr1h2","Maf","Zfp384","C4b","Ccl8","Cd22")
value_columns <- c("Trem2","C1qa","Cd68","Cyba","Fcer1g","Anxa1","Lrp1","Mif","Pld3","Pparg","Tyrobp","Cd300a")
UMAP <- DRG_Macrophage_aging@reductions$umap
DRG_Macrophage_aging <- ScaleData(DRG_Macrophage_aging,features = value_columns)
log_data <- GetAssayData(DRG_Macrophage_aging,slot = "scale.data")
genes_plot <- as.data.frame(t(log_data[value_columns,]))
UMAP_plot <- data.frame(x = UMAP@cell.embeddings[,1], y = UMAP@cell.embeddings[,2])
index <- match(rownames(genes_plot),rownames(UMAP_plot))
genes_plot <- genes_plot[index,]
head(genes_plot)
UMAP_plot_filter <- UMAP_plot[index,]
UMAP_plot_2 <- cbind(UMAP_plot_filter,genes_plot)
# 使用 lapply 创建图形列表
plot_creator <- function(col_name) {
  data_zero <- UMAP_plot_2[UMAP_plot_2[[col_name]] <= 0, ]
  data_nonzero <- UMAP_plot_2[UMAP_plot_2[[col_name]] > 0, ]
  
  ggplot() +
    geom_point(data = data_zero, 
               aes(x = x, y = y), 
               color = "grey", 
               size = 1.2) +
    geom_point(data = data_nonzero[order(data_nonzero[[col_name]]), ], 
               aes(x = x, y = y, color = .data[[col_name]]), 
               size = 1.2) +
    scale_color_gradient(low = "yellow", high = "red", name = col_name) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    coord_fixed() +
    ggtitle(col_name) +
    labs(x = "UMAP1", y = "UMAP2")
}

# 创建所有图形

plots <- lapply(value_columns, plot_creator)

# 使用 patchwork 组合图形
combined_plot <- wrap_plots(plots, nrow = 4, ncol = 4) + 
  plot_annotation(tag_levels = 'A')

print(combined_plot)
ggsave(filename = "Featureplot_Gpnmbs.pdf",plot = combined_plot)
