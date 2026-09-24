workingDir <- "E:\\DRG_aging\\paper_20251226\\Figure5\\Fig5c&d"
setwd(workingDir)

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(Matrix)
rm(list=ls())

##读取单细胞数据
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/DRG_Macrophage_aging_2_24.Rdata")
DRG_Macrophage_aging_gpnmb <- subset(DRG_Macrophage_aging, Celltype_3 %in% c("Mac_Ccr2","AAM_Trem2"))

##创建CDS对象并预处理数据
data <- GetAssayData(DRG_Macrophage_aging_gpnmb, assay = 'RNA', slot = 'counts')
cell_metadata <- DRG_Macrophage_aging_gpnmb@meta.data[,c(1:5,14,19)]
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$Celltype_3 <- factor(cds@colData$Celltype_3,levels = c("Mac_Ccr2","AAM_Trem2"))

cds@colData$age <- factor(cds@colData$age,levels = c("3 MO","24 MO"))

cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")

## 识别轨迹
cds <- learn_graph(cds)

cds<- order_cells(cds)

pdf(file = "Figure5aleft.pdf",width = 4,height = 3)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           cell_size = 0.7, label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()

saveRDS(cds, file = "cds_gpnmb.rds")
