library(Seurat)
library(ggplot2)
library(dplyr)
setwd("E:/web/")

combined_sct <- readRDS(file = "DRG_Neuron_ctrl_SCT_integrated.rds")
metadata <- combined_sct@meta.data[,c("Celltype_integration","celltype_2","Celltype_3","Celltype","ClusterName",
                                      "clusters_2","celltype","Level3","Dataset")]
head(metadata)

DRG_aging <- subset(combined_sct,Dataset %in% "Wang2026,et al.")
P_Wang2026 <- DimPlot(DRG_aging, reduction = "umap", group.by = "Celltype_3",cols = cols_list)+ coord_fixed()
P_Wang2026

DRG_Ginty <- subset(combined_sct,Dataset %in% "Ginty,et al.")
P_Ginty <- DimPlot(DRG_Ginty, reduction = "umap", group.by = "Celltype",cols = cols_list)+ coord_fixed()
P_Ginty

DRG_Woolf <- subset(combined_sct,Dataset %in% "Renthal,et al.")
P_Woolf <- DimPlot(DRG_Woolf, reduction = "umap", group.by = "Celltype",cols = cols_list)+ coord_fixed()
P_Woolf

DRG_mousebrain <- subset(combined_sct,Dataset %in% "Zeisel, et al.")
cols_list <- c(cols_list,"#143268")
P_mousebrain <- DimPlot(DRG_mousebrain, reduction = "umap", group.by = "ClusterName",cols = cols_list)+ coord_fixed()
P_mousebrain

DRG_Li <- subset(combined_sct,Dataset %in% "Li,et al.")
P_Li <- DimPlot(DRG_Li, reduction = "umap", group.by = "clusters_2",cols = cols_list)+ coord_fixed()
P_Li

DRG_Jung <- subset(combined_sct,Dataset %in% "Jung,et al.")
P_Jung <- DimPlot(DRG_Jung, reduction = "umap", group.by = "celltype",cols = cols_list)+ coord_fixed()
P_Jung

DRG_Usoskin <- subset(combined_sct,Dataset %in% "Usoskin,et al.")
P_Usoskin <- DimPlot(DRG_Usoskin, reduction = "umap", group.by = "Level3",cols = cols_list)+ coord_fixed()
P_Usoskin

DRG_Wang <- subset(combined_sct,Dataset %in% "Wang, et al.")
DRG_Wang$celltype_2 <- factor(DRG_Wang$celltype_2,levels = c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
  "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1")) 
P_Wang <- DimPlot(DRG_Wang, reduction = "umap", group.by = "celltype_2",cols = cols_list)+ coord_fixed()
P_Wang

library(patchwork)
combined <- (P_Wang2026 | P_Wang| P_Ginty | P_Woolf) / (P_mousebrain| P_Jung | P_Usoskin | P_Li)
ggsave("plot/UMAP_all_in_one.pdf", combined, width = 25, height = 15)

