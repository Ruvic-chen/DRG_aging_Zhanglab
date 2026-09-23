workingDir <- "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/"
setwd(workingDir)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

load("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/DRG_aging_all_2_24_final.Rdata")

DRG_Macrophage_aging <- subset(DRG_aging_all, Celltype_2_3 %in% c("Macrophage"))
DefaultAssay(object = DRG_Macrophage_aging) <- "RNA"

DRG_Macrophage_aging.list <- SplitObject(object = DRG_Macrophage_aging, split.by = "batch_2")
for (i in 1:length(x = DRG_Macrophage_aging.list)) {
  DRG_Macrophage_aging.list[[i]] <- NormalizeData(object = DRG_Macrophage_aging.list[[i]], scale.factor = 1e6, verbose = FALSE)
  DRG_Macrophage_aging.list[[i]] <- FindVariableFeatures(object = DRG_Macrophage_aging.list[[i]], selection.method = "vst", 
                                                         nfeatures = 2000, verbose = FALSE)
}
reference.list <- DRG_Macrophage_aging.list[c("10X_1_batch3","10X_2_batch1","10X_2_batch2")]
DRG_Macrophage_aging.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
DRG_Macrophage_aging <- IntegrateData(anchorset = DRG_Macrophage_aging.anchors, features.to.integrate = rownames(DRG_Macrophage_aging),k.weight = 80, dims = 1:30)
DefaultAssay(object = DRG_Macrophage_aging) <- "integrated"

# Run the standard workflow for visualization and clustering
DRG_Macrophage_aging <- ScaleData(object = DRG_Macrophage_aging, verbose = FALSE)
DRG_Macrophage_aging <- RunPCA(object = DRG_Macrophage_aging, npcs = 30, verbose = FALSE)
DRG_Macrophage_aging <- RunUMAP(object = DRG_Macrophage_aging, reduction = "pca", dims = 1:30)
DRG_Macrophage_aging <- RunTSNE(object = DRG_Macrophage_aging, reduction = "pca", dims = 1:30)
DRG_Macrophage_aging <- FindNeighbors(object = DRG_Macrophage_aging, dims = 1:30, verbose = FALSE)
DRG_Macrophage_aging <- FindClusters(object = DRG_Macrophage_aging,resolution = 1.2, verbose = FALSE)


pdf(file = "Figure4a.pdf",width = 6,height = 4)
DimPlot(object = DRG_Macrophage_aging, reduction = "umap", label = T, ncol = 3, 
        cols = 
          # c("#FF7F00", "#377EB8", "#4DAF4A", "#984EA3", "#E41A1C"),
          brewer.pal(7,"Accent")[-4],
        group.by = "Celltype_3", split.by = 'age', repel = TRUE) + coord_fixed()
dev.off()

DefaultAssay(DRG_Macrophage_aging) <- "RNA"
features_marker <- c("Mrc1","Isg15","Lyve1","Folr2","Ccr2","Cd163","Ccl8","Gpnmb")
pdf(file = "Figure4b.pdf",width = 4.8,height = 2.5)
DotPlot(DRG_Macrophage_aging, features = features_marker,group.by = "Celltype_3", cols = c("lightgrey", "red"))+theme(axis.text.x = element_text(angle = 90,                                                                                                                                                          vjust = 0.3, hjust=0.3, face = "italic"))
dev.off()





