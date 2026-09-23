setwd("/paper数据整理/Figure 2/Figure 2b/")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("/paper数据整理/Figure 2/Figure 2b/DRG_Neuron_aging_2_24_final.Rdata")
load("/paper数据整理/Figure 2/Figure 2b/cols_list.Rdata")
pdf(file = "UMAP_DRG_Neuron_aging_subtype_all_split.pdf",width = 17,height = 5)
DimPlot(object = DRG_Neuron_aging, reduction = "umap", label = F, 
        cols = cols_list, split.by = "age", group.by = "Celltype_2_3", repel = TRUE) + coord_fixed()
dev.off()

