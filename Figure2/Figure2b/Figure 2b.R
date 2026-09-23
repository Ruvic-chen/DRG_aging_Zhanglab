setwd("/paper数据整理/Figure 2/Figure 2a/")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("./DRG_aging_all_2_24_final.Rdata")

ident_levels <- c("Neuron","Satellite","Schwann_M","Schwann_N","Fibroblast","VEC","VSMC","Macrophage","Monocyte", "Neutrophil","B Cell", "T Cell")
DRG_aging_all$Celltype_2_2 <- factor(DRG_aging_all$Celltype_2_2,levels = ident_levels)
color_all <- brewer.pal(12,"Paired")
names(color_all) <- c("VEC","Fibroblast","Schwann_M","Schwann_N","Satellite","Neuron","VSMC","Neutrophil","Monocyte","Macrophage","B Cell","T Cell")

pdf(file = "UMAP_DRG_aging_all_celltype_all_split.pdf",width = 17,height = 5)
DimPlot(object = DRG_aging_all, reduction = "umap", label = F, 
        cols = color_all,split.by = "age", 
        group.by = "Celltype_2_2", repel = TRUE) + coord_fixed()
dev.off()