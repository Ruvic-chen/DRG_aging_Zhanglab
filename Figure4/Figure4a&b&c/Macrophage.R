workingDir <- "E:\\DRG_aging\\paper_Wang&Chen2025\\Figure4\\Fig4a&b"
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

old_idents <- c("Mac_IFN","Mac_TLF","Mac_Ccr2","Mac_Cd163","AAM_Ccl8","AAM_Trem2")
new_idents <- c("Mac_IFN","Mac_TLF","Mac_Ccr2","Mac_Cd163","AAMac_Ccl8","AAMac_Gpnmb")
DRG_Macrophage_aging$Celltype_3 <- plyr::mapvalues(x = DRG_Macrophage_aging$Celltype_3, from = old_idents, to = new_idents)


pdf(file = "UMAP_DRG_Macrophage_aging_celltype_group_split.pdf",width = 6,height = 4)
DimPlot(object = DRG_Macrophage_aging, reduction = "umap", label = T, ncol = 3, 
        cols = 
          # c("#FF7F00", "#377EB8", "#4DAF4A", "#984EA3", "#E41A1C"),
          brewer.pal(7,"Accent")[-4],
        group.by = "Celltype_3", split.by = 'age', repel = TRUE) + coord_fixed()
dev.off()

DefaultAssay(DRG_Macrophage_aging) <- "RNA"


features_marker <- c("Mrc1","Isg15","Lyve1","Folr2","Ccr2","Cd163","Ccl8","Gpnmb")
pdf(file = "Celltype_markers.pdf",width = 4.8,height = 2.5)
DotPlot(DRG_Macrophage_aging, features = features_marker,group.by = "Celltype_3", cols = c("lightgrey", "red"))+theme(axis.text.x = element_text(angle = 90,                                                                                                                                                          vjust = 0.3, hjust=0.3, face = "italic"))
dev.off()

save(DRG_Macrophage_aging,DRG_Macrophage_aging.markers,file = "DRG_Macrophage_aging_2_24.Rdata")


library(sccomp)
library(tidyr)
DefaultAssay(DRG_Macrophage_aging) <- "RNA"
sccomp_result = 
  DRG_Macrophage_aging |>
  sccomp_estimate( 
    formula_composition = ~ age, 
    #formula_variability = ~age,
    sample = "orig.ident", 
    cell_group = "Celltype_3", 
    cores = 1,
    verbose = FALSE
  ) |> 
  sccomp_test()
p <- plot(sccomp_result)
pdf(file = "credible_intervals_1D_nocount.pdf",width = 4,height = 5)
p$credible_intervals_1D
dev.off()
write.csv(sccomp_result,file = "sccomp_result_Mac.csv")


sample_anno_filter <- DRG_Macrophage_aging@meta.data
# sample_anno_filter$Celltype_3 <- factor(sample_anno_filter$Celltype_3,levels = ident_levels)
# table(sample_anno$Celltype_3)

sample_anno_filter$age <- factor(sample_anno_filter$age,levels = c("3 MO","24 MO"))
sample_anno_filter <- sample_anno_filter[,c(c("orig.ident","age","sex","Celltype_3"))]
head(sample_anno_filter)

library(reshape2)
sample_anno_2 <- dcast(sample_anno_filter,sex+age+orig.ident ~ Celltype_3)
head(sample_anno_2)
sample_anno_2 <- dcast(sample_anno_filter,age ~Celltype_3)
head(sample_anno_2)
rownames(sample_anno_2) <- sample_anno_2$age
sample_anno_2 <- sample_anno_2[,-1]
row_num <- rowSums(sample_anno_2)
f1 <- function( x ) x * 10000 / row_num
sample_anno_3 <- as.data.frame(apply(sample_anno_2,2,f1))
sum_num <- colSums(sample_anno_3)
f2 <- function( x ) x / sum_num
sample_anno_percentage <- as.data.frame(apply(sample_anno_3,1,f2))
sample_anno_percentage$celltype <- rownames(sample_anno_percentage)
head(sample_anno_percentage)
sample_anno_percentage_melt <- melt(sample_anno_percentage,id.vars = "celltype")
head(sample_anno_percentage_melt)

p <- ggplot() + geom_bar(data =sample_anno_percentage_melt, 
                         aes(x = factor(celltype,levels = new_idents), y = value, fill = factor(variable)), 
                         stat = "identity", position = "fill")+
  scale_fill_manual(values=c("#a9c844","#b2672a"))+
  geom_hline(yintercept=0.5,linetype='dotted')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
  theme(axis.text = element_text(color="black"),panel.border = element_blank())
p
ggsave(filename = "barplot_subtype_age_ratio_2.pdf",plot = p,width = 5,height = 3)
