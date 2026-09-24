workingDir <- "F:\\paper_data_Wang&Chen_20260919\\Figure5\\Fig5b&c"
setwd(workingDir)

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(Matrix)
rm(list=ls())

##读取单细胞数据
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/DRG_Macrophage_aging_2_24.Rdata")
DRG_Macrophage_aging_ccl8 <- subset(DRG_Macrophage_aging, Celltype_3 %in% c("Mac_Cd163","AAM_Ccl8"))

##创建CDS对象并预处理数据
data <- GetAssayData(DRG_Macrophage_aging_ccl8, assay = 'RNA', slot = 'counts')
cell_metadata <- DRG_Macrophage_aging_ccl8@meta.data[,c(1:5,14,19)]
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$Celltype_3 <- factor(cds@colData$Celltype_3,levels = c("Mac_Cd163","AAM_Ccl8"))

cds@colData$age <- factor(cds@colData$age,levels = c("3 MO","24 MO"))

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds,color_cells_by="Celltype_3")
colnames(colData(cds))
plot_cells(cds, genes=c("Ccl8","Trem2"))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Celltype_3") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(DRG_Macrophage_aging, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Celltype_3") + ggtitle('int.umap')
p2
p = p1|p2
p
ggsave("Reduction_Compare.pdf", plot = p, width = 10, height = 5)

## Monocle3聚类分区
cds <- cluster_cells(cds, reduction_method = "UMAP")
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = p1|p2
p
ggsave("Cell_partition.pdf", plot = p, width = 10, height = 5)

## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "Celltype_3", label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, label_branch_points = FALSE)
p
ggsave("Trajectory_celltype.pdf", plot = p, width = 8, height = 6)

p = plot_cells(cds, color_cells_by = "age", label_cell_groups=FALSE, label_leaves=TRUE,
               label_branch_points=TRUE, graph_label_size=1.5)
p
ggsave("Trajectory_order.pdf", plot = p, width = 8, height = 6)

cds_subset <- choose_cells(cds)
cds_subset<- order_cells(cds_subset)


pdf(file = "Macrophage_pseudotime.pdf",width = 4,height = 3)
plot_cells(cds_subset, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           cell_size = 0.7, label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()

p <- plot_cells(cds_subset,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=2,
                cell_size = 0.7) + coord_fixed(ratio = 1.2)
p
ggsave("Trajectory_Pseudotime_ccl8.pdf", plot = p, width = 4, height = 3)
saveRDS(cds_subset, file = "cds_ccl8.rds")

