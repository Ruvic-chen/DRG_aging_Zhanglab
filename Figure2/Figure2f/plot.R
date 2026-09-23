workingDir <- "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/"
setwd(workingDir)
library(ClusterGVis)
library(ggplot2)
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/DRG_Neuron_aging_2_24_final.Rdata")
library(Seurat)
DefaultAssay(DRG_Neuron_aging) <- "RNA"
pdf(file = "Celltype_markers_2.pdf",width = 4,height = 5)
features_marker <- c("Ackr1")
DotPlot(DRG_Neuron_aging, features = features_marker,group.by = "Celltype_3", cols = c("lightgrey", "red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.3, face = "italic"))
dev.off()

Neuron.markers <- read.csv(file = "DRG_neuron_markers.csv")
markGenes <- unique(Neuron.markers$gene_name)
head(Neuron.markers)
#table(DRG_Neuron_aging$Celltype_2_3)
Idents(DRG_Neuron_aging) <- "Celltype_3"
Neuron.markers.all <- Seurat::FindAllMarkers(DRG_Neuron_aging,
                                             only.pos = TRUE,
                                             min.pct = 0.25,
                                             logfc.threshold = 0.25)
Neuron.markers <- subset(Neuron.markers.all, gene %in% Neuron.markers$gene_name) 
Neuron.markers <- Neuron.markers %>% dplyr::group_by(cluster)
Neuron.data <- prepareDataFromscRNA(object = DRG_Neuron_aging,
                                    diffData = Neuron.markers,
                                    group.by = "Celltype_2_3",
                                    showAverage = TRUE)
markGenes <- intersect(markGenes,Neuron.markers$gene)
Neuron.data$wide.res <- Neuron.data$wide.res[
  match(markGenes, Neuron.data$wide.res$gene),   # 按 markGenes 顺序索引
  , drop = FALSE
]


library(ComplexHeatmap)
library(circlize)   # 用于颜色渐变

wide_df <- Neuron.data$wide.res
# 假设第一列名为 "gene"，其余为样本表达值
gene_col <- wide_df[, 1]                 # 基因名
mat_data <- wide_df[, c(-1,-18)]                # 数值矩阵（当前为字符）

# 2. 将整个数值矩阵转为数值型（注意：非数值会变成 NA）
mat <- as.matrix(mat_data)
#mode(mat) <- "numeric"                   # 强制转为数值
# 或使用 apply: mat <- apply(mat_data, 2, as.numeric)

# 3. 恢复行名
rownames(mat) <- gene_col

# 4. 按 markGenes 排序（若需要）
mat <- mat[match(markGenes, rownames(mat)), , drop = FALSE]
mat <- mat[!is.na(rownames(mat)), ]     # 剔除未匹配基因

# 3. 确定列分组（样本属于哪个 cluster）
#    方法：从列名中提取分组信息，或根据您的 cluster.order 映射。
#    示例：假设列名格式为 "clusterX_sampleY"，或您有分组向量。
#    如果列顺序本身已经按 cluster.order 排列，可以直接用 cluster.order 作为分组因子。
#    这里假设您有 14 个 cluster，每个 cluster 有若干样本，总列数 = ncol(mat)
n_samples <- ncol(mat)
# 需根据实际情况构建分组因子，例如：
cluster_labels <- rep(1:16, each = n_samples / 16)  # 仅当每类样本数相等时适用
# 更安全的方式：从原始数据中提取分组信息（如 DRG_Neuron_aging 的 metadata）
# 示例：使用 colnames 中的前缀，请根据您的数据调整。
# 若您有样本分组向量，直接赋值：
# sample_group <- factor(your_group_vector)

# 4. 设置列颜色（使用您的 cols_list，需与分组一一对应）
#    cols_list 应有 14 个颜色，顺序与 cluster.order 一致。
col_anno <- HeatmapAnnotation(
  Cluster = cluster_labels,
  col = list(Cluster = setNames(cols_list, 1:16))   # 将颜色与分组值对应
)

# 5. 设置热图颜色（Z-score 常用蓝-白-红）
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# 6. 绘制热图
pdf('sc1_manual.pdf', height = 10, width = 8)
Heatmap(mat,
        name = "Z-score",
        col = col_fun,
        cluster_rows = FALSE,          # 完全按您的行顺序
        cluster_columns = FALSE,       # 列不聚类（保持原始顺序）
        show_row_names = TRUE,
        show_column_names = TRUE,
        top_annotation = col_anno,     # 添加顶部样本分组注释
        row_title = NULL,
        column_title = NULL,
        heatmap_legend_param = list(title = "Z-score"))
dev.off()

