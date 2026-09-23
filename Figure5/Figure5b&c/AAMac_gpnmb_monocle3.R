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
#cell_metadata <- DRG_Macrophage_aging_2@meta.data
cell_metadata <- DRG_Macrophage_aging_gpnmb@meta.data[,c(1:5,14,19)]
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$Celltype_3 <- factor(cds@colData$Celltype_3,levels = c("Mac_Ccr2","AAM_Trem2"))

cds@colData$age <- factor(cds@colData$age,levels = c("3 MO","24 MO"))

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
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

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes_gpnmb.csv", row.names = F)


pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(DRG_Macrophage_aging_gpnmb@meta.data)]
pseudotime[is.infinite(pseudotime)] <- 20
DRG_Macrophage_aging_gpnmb$pseudotime <- pseudotime
save(DRG_Macrophage_aging_gpnmb,file = "DRG_Macrophage_aging_gpnmb.Rdata")

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
cds <- readRDS("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/cds_gpnmb.rds")
Track_genes <- readRDS("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/all_Track_genes.rds")
load("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/DRG_Macrophage_aging_gpnmb.Rdata")
Track_genes_Gpnmb <- subset(Track_genes,type %in% "AAMac_Gpnmb")
top10 <- Track_genes_Gpnmb %>% top_n(-100, q_value)
genes <- unique(top10$gene_short_name)
pt.matrix <- normalized_counts(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 3 groups
df <- DRG_Macrophage_aging_gpnmb@meta.data[,c(4,19,24)]

df <- df[order(df$pseudotime),]
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

col_pseudotime <- colorRamp2(
  breaks = seq(from = 0, to = 18, by = 18/8), 
  colors = brewer.pal(9,"BuPu")
)
col_con <- brewer.pal(3,"Dark2")[-2]
names(col_con) <- c("3 MO","24 MO")

col_clusters <- brewer.pal(3,"Set1")[-2]
names(col_clusters) <- c("Mac_Ccr2","AAM_Trem2")

ha = HeatmapAnnotation(df = df[,1:3],
                       col = list(age = col_con,
                                  celltype = col_clusters,
                                  pseudotime = col_pseudotime))

htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 3,
  top_annotation = ha,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm
pdf(file = "Figure5a_right.pdf",width = 6,height = 10)
htkm
dev.off()
