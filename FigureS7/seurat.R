library(Seurat)
library(dplyr)
setwd(dir = "E:\\DRG_aging\\dataset_Cd163\\GSE195507_RAW")

# 1. 定义文件列表（请按实际文件名调整）
files <- c(
  "GSM5839574_Macrophage_1_Old_filtered_feature_bc_matrix.h5",
  "GSM5839575_Macrophage_1_Young_filtered_feature_bc_matrix.h5",
  "GSM5839576_Macrophage_2_Old_filtered_feature_bc_matrix.h5",
  "GSM5839577_Macrophage_2_Young_filtered_feature_bc_matrix.h5",
  "GSM5839578_Macrophage_3_Old_filtered_feature_bc_matrix.h5",
  "GSM5839579_Macrophage_3_Young_filtered_feature_bc_matrix.h5"
)

# 2. 从文件名解析样本元数据
sample_info <- data.frame(
  file = files,
  sample_id = gsub("_filtered_feature_bc_matrix\\.h5", "", files),  # 如 "GSM5839574_Macrophage_1_Old"
  individual = rep(1:3, each = 2),                                  # 个体编号
  age = rep(c("Old", "Young"), times = 3)                          # 年龄分组
)
# 可自行添加更多信息（如批次等）

# 3. 循环读取每个h5文件，创建Seurat对象，并存入列表
seurat_list <- list()
for (i in 1:nrow(sample_info)) {
  cat("Reading", sample_info$file[i], "...\n")
  
  # 读取稀疏矩阵（基因 × 细胞）
  counts <- Read10X_h5(sample_info$file[i])
  
  # 创建Seurat对象，用样本ID作为项目名
  obj <- CreateSeuratObject(counts = counts, 
                            project = sample_info$sample_id[i],
                            min.cells = 0,      # 不过滤（已有过滤版）
                            min.features = 0)
  
  # 存入样本元数据
  obj$individual <- sample_info$individual[i]
  obj$age <- sample_info$age[i]
  obj$sample <- sample_info$sample_id[i]
  
  seurat_list[[i]] <- obj
}

# 4. 合并所有对象（按细胞合并，基因取交集）
#    add.cell.ids 为每个细胞的barcode加上样本前缀，避免重名
merged_seurat <- merge(x = seurat_list[[1]], 
                       y = seurat_list[-1], 
                       add.cell.ids = sample_info$sample_id,
                       project = "Macrophage_merged")

# 5. 查看合并后的结果
merged_seurat
# 输出示例：An object of class Seurat with N features across M samples within 1 assay.

# 6. （可选）保存合并后的对象，便于后续分析
saveRDS(merged_seurat, file = "Macrophage_merged.rds")

# analysis
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)   # 多图组合
#library(harmony)     # 用于批次校正（可选，强烈推荐）
merged <- readRDS("Macrophage_merged.rds")
head(merged@meta.data)

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

# 可视化QC指标（可选）
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "age", pt.size = 0.1)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# 先ScaleData（Harmony要求）
merged <- ScaleData(merged, vars.to.regress = c("percent.mt"))  # 可选回归线粒体

# 运行PCA
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

# 使用Harmony校正个体（individual）的批次效应
# merged <- RunHarmony(merged, group.by.vars = "individual",
#                      reduction = "pca")

# 此时降维结果存储在 merged[["harmony"]] 中

# 基于Harmony校正后的嵌入进行聚类
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)  # 可调整分辨率
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30, reduction.name = "umap")
# # 运行UMAP（基于harmony）
# merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, 
#                   reduction.name = "umap.harmony")

# 可视化：按cluster和age着色
p1 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Clusters")
p2 <- DimPlot(merged, reduction = "umap", group.by = "age", cols = c("blue", "red")) + 
  ggtitle("Age groups")
p1 + p2
merged <- JoinLayers(merged)
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, test.use = "wilcox")

# 提取每个cluster top5基因
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 查看top5，手动判断
print(top5, n = 50)

top10 <- markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
gene_markers <- c("Mrc1","Lyve1","Folr2","Fcrls","Isg15","Cx3cr1","Bcl2a1b","Cd163","Neat1","Malat1",
                  "Ccr2","Plac8","S100a4","Cd3e","Cd8a","Cd69","Cd79a","Cd79b","S100a9",
                  "S100a8","Naaa","Xcr1","Cd80","Cd86","Il6","Cxcl9","Cxcl10","Ccr7","Cd163","Il10","Ccr2","Ccr5")


gene_heatmap <- c(gene_markers,top10$gene)
#merged <- ScaleData(merged,features = gene_heatmap)
library(ggplot2)
p1 <- DoHeatmap(object = merged, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_merged_muscle.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = merged, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_merged_ctrl_raw.pdf",height = 40, width = 20)
p2
dev.off()

VlnPlot(merged,features = c("nCount_RNA","nFeature_RNA"))


mix <- c(1,4,8,14)
monocytes <- c(3,6)
mac <- c(0,2,5,9,10,11,12)
neutrophil <- c(7,13,16)
B_cell <- c(15)

merged_mac <- subset(merged,seurat_clusters %in% mac)
merged_mac <- NormalizeData(merged_mac, normalization.method = "LogNormalize", scale.factor = 10000)
merged_mac <- FindVariableFeatures(merged_mac, selection.method = "vst", nfeatures = 2000)

# 先ScaleData（Harmony要求）
merged_mac <- ScaleData(merged_mac, vars.to.regress = c("percent.mt"))  # 可选回归线粒体

# 运行PCA
merged_mac <- RunPCA(merged_mac, npcs = 50, verbose = FALSE)

# 使用Harmony校正个体（individual）的批次效应
# merged_mac <- RunHarmony(merged_mac, group.by.vars = "individual",
#                      reduction = "pca")

# 此时降维结果存储在 merged_mac[["harmony"]] 中

# 基于Harmony校正后的嵌入进行聚类
merged_mac <- FindNeighbors(merged_mac, reduction = "pca", dims = 1:30)
merged_mac <- FindClusters(merged_mac, resolution = 0.5)  # 可调整分辨率
merged_mac <- RunUMAP(merged_mac, reduction = "pca", dims = 1:30, reduction.name = "umap")
# # 运行UMAP（基于harmony）
# merged_mac <- RunUMAP(merged_mac, reduction = "harmony", dims = 1:30, 
#                   reduction.name = "umap.harmony")

# 可视化：按cluster和age着色
p1 <- DimPlot(merged_mac, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Clusters")
p2 <- DimPlot(merged_mac, reduction = "umap", group.by = "age", cols = c("blue", "red")) + 
  ggtitle("Age groups")
p1 + p2
#merged_mac <- JoinLayers(merged_mac)
markers <- FindAllMarkers(merged_mac, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, test.use = "wilcox")
saveRDS(markers,file = "mac_markers.rds")
saveRDS(merged_mac,file = "merged_mac.rds")

markers <- readRDS(file = "mac_markers.rds")
merged_mac <- readRDS(file = "merged_mac.rds")

top10 <- markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
gene_markers <- c("Mrc1","Lyve1","Folr2","Fcrls","Isg15","Cx3cr1","Bcl2a1b","Cd163","Neat1","Malat1",
                  "Ccr2","Plac8","S100a4","Cd3e","Cd8a","Cd69","Cd79a","Cd79b","S100a9",
                  "S100a8","Naaa","Xcr1","Cd80","Cd86","Il6","Cxcl9","Cxcl10","Ccr7","Cd163","Il10","Ccr2","Ccr5",
                  "Ccl8","Trem2","Gpnmb","Spp1")


gene_heatmap <- c(gene_markers,top10$gene)
merged_mac <- ScaleData(merged_mac,features = gene_heatmap)
library(ggplot2)
p1 <- DoHeatmap(object = merged_mac, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_merged_mac_muscle.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = merged_mac, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_merged_mac_ctrl_raw.pdf",height = 40, width = 20)
p2
dev.off()

# Endothelial cell 7
# 中性粒细胞 10
merged_mac_2 <- subset(merged_mac,seurat_clusters %in% c(0:6,8,9))
merged_mac_2 <- NormalizeData(merged_mac_2, normalization.method = "LogNormalize", scale.factor = 10000)
merged_mac_2 <- FindVariableFeatures(merged_mac_2, selection.method = "vst", nfeatures = 2000)

# 先ScaleData（Harmony要求）
merged_mac_2 <- ScaleData(merged_mac_2, vars.to.regress = c("percent.mt"))  # 可选回归线粒体

# 运行PCA
merged_mac_2 <- RunPCA(merged_mac_2, npcs = 50, verbose = FALSE)

# 使用Harmony校正个体（individual）的批次效应
# merged_mac_2 <- RunHarmony(merged_mac_2, group.by.vars = "individual",
#                      reduction = "pca")

# 此时降维结果存储在 merged_mac_2[["harmony"]] 中

# 基于Harmony校正后的嵌入进行聚类
merged_mac_2 <- FindNeighbors(merged_mac_2, reduction = "pca", dims = 1:30)
merged_mac_2 <- FindClusters(merged_mac_2, resolution = 0.5)  # 可调整分辨率
merged_mac_2 <- RunUMAP(merged_mac_2, reduction = "pca", dims = 1:30, reduction.name = "umap")
# # 运行UMAP（基于harmony）
# merged_mac_2 <- RunUMAP(merged_mac_2, reduction = "harmony", dims = 1:30, 
#                   reduction.name = "umap.harmony")

# 可视化：按cluster和age着色
p1 <- DimPlot(merged_mac_2, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Clusters")
p2 <- DimPlot(merged_mac_2, reduction = "umap", group.by = "age", cols = c("blue", "red")) + 
  ggtitle("Age groups")
p1 + p2
#merged_mac_2 <- JoinLayers(merged_mac_2)
markers <- FindAllMarkers(merged_mac_2, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, test.use = "wilcox")
saveRDS(markers,file = "mac_markers.rds")
saveRDS(merged_mac_2,file = "merged_mac_2.rds")

markers <- readRDS(file = "mac_markers.rds")
merged_mac_2 <- readRDS(file = "merged_mac_2.rds")

top10 <- markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
gene_markers <- c("Mrc1","Lyve1","Folr2","Fcrls","Isg15","Cx3cr1","Bcl2a1b","Cd163","Neat1","Malat1",
                  "Ccr2","Plac8","S100a4","Cd3e","Cd8a","Cd69","Cd79a","Cd79b","S100a9",
                  "S100a8","Naaa","Xcr1","Cd80","Cd86","Il6","Cxcl9","Cxcl10","Ccr7","Cd163","Il10","Ccr2","Ccr5",
                  "Ccl8","Trem2","Gpnmb","Spp1")


gene_heatmap <- c(gene_markers,top10$gene)
merged_mac_2 <- ScaleData(merged_mac_2,features = gene_heatmap)
library(ggplot2)
p1 <- DoHeatmap(object = merged_mac_2, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_merged_mac_2_muscle.pdf",height = 40, width = 20)
p1
dev.off()
