workingDir <- "E:\\DRG_aging\\paper_Wang&Chen2025\\Figure5\\Fig5ab"
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
#cell_metadata <- DRG_Macrophage_aging_2@meta.data
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

##细胞按拟时排序
# 解决order_cells(cds)报错"object 'V1' not found"
#rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
#colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds_subset <- choose_cells(cds)
cds_subset<- order_cells(cds_subset)

#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

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

Track_genes <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
#Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes_ccl8.csv", row.names = F)

#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#Track_genes <- ciliated_cds_pr_test_res
#Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
#基因表达趋势图
monocle3::plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by = "age",
                                   min_expr=0.5, ncol = 2)

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(DRG_Macrophage_aging_ccl8@meta.data)]
pseudotime[is.infinite(pseudotime)] <- 20
DRG_Macrophage_aging_ccl8$pseudotime <- pseudotime
save(DRG_Macrophage_aging_ccl8,file = "DRG_Macrophage_aging_ccl8.Rdata")

#plot after Track gene filter
##########################
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
cds <- readRDS("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/cds_ccl8.rds")
Track_genes <- readRDS("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/all_Track_genes.rds")
load("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/DRG_Macrophage_aging_ccl8.Rdata")
Track_genes_ccl8 <- subset(Track_genes,type %in% "AAMac_Ccl8")
top10 <- Track_genes_ccl8 %>% top_n(-100, q_value)
genes <- unique(top10$gene_short_name)
pt.matrix <- normalized_counts(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 3 groups
df <- DRG_Macrophage_aging_ccl8@meta.data[,c(4,19,24)]
df <- df[colnames(cds),]
df <- df[order(df$pseudotime),]
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

col_pseudotime <- colorRamp2(
  breaks = seq(from = 0, to = 18, by = 18/8), 
  colors = brewer.pal(9,"BuPu")
)
col_con <- brewer.pal(3,"Dark2")[-2]
names(col_con) <- c("3 MO","24 MO")

col_clusters <- brewer.pal(3,"Set1")[-2]
names(col_clusters) <- c("Mac_Cd163","AAM_Ccl8")

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
pdf(file = "heatmap_pseudotime_ccl8_top100.pdf",width = 6,height = 10)
htkm
dev.off()
