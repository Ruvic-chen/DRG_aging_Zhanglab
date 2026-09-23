library(Seurat)
library(ggplot2)
library(dplyr)

setwd("E:/DRG_aging/dataset_Cd163/GSE195507_RAW/integrated/")
options(future.globals.maxSize = 18 * 1024^3) 
#krasniewski
Kras_mac <- readRDS(file = "E:/DRG_aging/dataset_Cd163/GSE195507_RAW/merged_mac_2.rds")
head(Kras_mac)
Kras_mac$Dataset <- "krasniewski,et al."

load("E:/DRG_aging/dataset_Cd163/DRG_Macrophage_aging_2_24.Rdata")
Wang_DRG_Mac_count <- DRG_Macrophage_aging@assays$RNA@counts
Wang_mouse_DRG_Mac_meta <- DRG_Macrophage_aging@meta.data[,c(1,4,19)]
colnames(Wang_mouse_DRG_Mac_meta)[3] <- "Cell.types"
Wang_mouse_DRG_Mac <- CreateSeuratObject(Wang_DRG_Mac_count)
Wang_mouse_DRG_Mac@meta.data <- Wang_mouse_DRG_Mac_meta
Wang_mouse_DRG_Mac$Dataset <- "Wang,et al."

head(Wang_mouse_DRG_Mac)

data_all <- merge(x = Kras_mac, y= Wang_mouse_DRG_Mac)
datalist <- SplitObject(data_all, split.by = "Dataset")

# ==================== 步骤1：对每个数据集分别运行SCTransform ====================
# 注意：SCTransform() 内部包含 NormalizeData + FindVariableFeatures + ScaleData
cat("对各数据集独立运行 SCTransform...\n")
for (i in 1:length(datalist)) {
  # 如果需要回归线粒体基因比例等混杂因素，可通过 vars.to.regress 参数添加
  # 例如: datalist[[i]] <- SCTransform(datalist[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
  datalist[[i]] <- SCTransform(datalist[[i]], verbose = FALSE)
  cat("数据集", i, "SCTransform 完成\n")
}

# ==================== 步骤2：选择整合特征 ====================
# SelectIntegrationFeatures 会从每个数据集中挑选前k个高变基因，取交集或排名靠前的基因
cat("\n选择用于整合的特征基因...\n")
integration_features <- SelectIntegrationFeatures(object.list = datalist,
                                                  nfeatures = 3000,
                                                  verbose = FALSE)

#integration_features <- intersect(unique(DRG_Mac.markers$gene),integration_features)
# ==================== 步骤3：PrepSCTIntegration（关键步骤） ====================
# 该函数确保所有 SCT 模型间的必要Pearson残差已被正确计算，为整合做准备
cat("\n准备 SCT 整合...\n")
datalist <- PrepSCTIntegration(object.list = datalist, 
                               anchor.features = integration_features,
                               verbose = FALSE)

# ==================== 步骤4：查找锚点 ====================
# 指定 normalization.method = "SCT" 告知函数使用已生成的SCT残差作为输入[reference:2]
# reduction = "cca" 使用典型相关分析降维（锚点默认方法）
cat("\n查找锚点...\n")
anchors <- FindIntegrationAnchors(object.list = datalist,
                                  anchor.features = integration_features,
                                  normalization.method = "SCT",   # 使用 SCT 残差
                                  reduction = "cca",               # 默认是 cca
                                  dims = 1:30,
                                  verbose = TRUE)

# ==================== 步骤5：整合数据 ====================
# IntegrateData 基于找到的锚点进行整合，同样需要 normalization.method = "SCT"[reference:3]
cat("\n开始整合数据...\n")
combined_sct <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCT",
                              dims = 1:30,
                              verbose = TRUE)

# ==================== 步骤6：整合后处理 ====================
# 默认整合后生成的 assay 为 "integrated"，将其设为默认 assay
DefaultAssay(combined_sct) <- "integrated"

# 运行标准降维与聚类（整合后的数据已经是log归一化形式，无需再次 NormalizeData）
combined_sct <- ScaleData(combined_sct, verbose = FALSE)
combined_sct <- RunPCA(combined_sct, npcs = 50, verbose = FALSE)
combined_sct <- RunUMAP(combined_sct, reduction = "pca", dims = 1:20, verbose = FALSE,seed.use = 50)
combined_sct <- FindNeighbors(combined_sct, reduction = "pca", dims = 1:20)
combined_sct <- FindClusters(combined_sct, resolution = 1)

# ==================== 步骤7：可视化整合效果 ====================
dataset_order <- c("Wang,et al.",  "krasniewski,et al.")
combined_sct$Dataset <- factor(combined_sct$Dataset,levels = rev(dataset_order))
p1 <- DimPlot(combined_sct, reduction = "umap", group.by = "Dataset",split.by = "age") + 
  ggtitle("SCT整合后：各数据集分布") + coord_fixed()
p1
ggsave(plot = p1,filename = "UMPA_integrated.pdf")
p2 <- DimPlot(combined_sct, reduction = "umap", label = TRUE) + 
  ggtitle("SCT整合后：聚类结果") + coord_fixed()
p2
ggsave(plot = p2,filename = "UMPA_cluster.pdf")
p <-p1 + p2
ggsave(filename = "UMAP.pdf",plot = p,width = 10,height = 5)

load("E:/web/cols_list.Rdata")
#DimPlot(combined_sct, reduction = "umap", group.by = "Celltype_3",split.by = "Dataset",cols = cols_list,ncol = 3)+ coord_fixed()
mac_Wang <- subset(combined_sct,Dataset %in% c("Wang,et al."))
DimPlot(mac_Wang, reduction = "umap", group.by = "Cell.types",split.by = "Dataset",cols = cols_list,label = T)+ coord_fixed()

mac_krasniewski <- subset(combined_sct,Dataset %in% c("Wang,et al."))
DimPlot(mac_krasniewski, reduction = "umap", group.by = "Cell.types",split.by = "Dataset",cols = cols_list,label = T)+ coord_fixed()

DimPlot(combined_sct, reduction = "umap", group.by = "Cell.types",split.by = "Dataset",cols = cols_list,label = T)+ coord_fixed()
saveRDS(combined_sct,file = "combined_sct.rds")
combined_sct <- readRDS(file = "combined_sct.rds")

p <- DimPlot(combined_sct, reduction = "umap", split.by = "Dataset",ncol = 2)+ coord_fixed()
p
ggsave(plot = p2,filename = "UMAP_integrated_ctrl.pdf")
DRG_Mac.markers <- FindAllMarkers(object = combined_sct,test.use = "roc",min.pct = 0.2,only.pos = T)
saveRDS(DRG_Mac.markers,file = "DRG_Mac_integrated.markers.rds")
DRG_Mac.markers <- readRDS(file = "DRG_Mac_integrated.markers.rds")
library(dplyr)
top10 <- DRG_Mac.markers %>% group_by(cluster) %>% top_n(20, myAUC)
gene_markers <- c("Mrc1","Lyve1","Folr2","Fcrls","Isg15","Cx3cr1","Bcl2a1b","Cd163","Neat1","Malat1",
                  "Ccr2","Plac8","S100a4","Cd3e","Cd8a","Cd69","Cd79a","Cd79b","S100a9",
                  "S100a8","Naaa","Xcr1","Cd80","Cd86","Il6","Cxcl9","Cxcl10","Ccr7","Cd163","Il10","Ccr2","Ccr5","Ccl8","Gpnmb","Trem2")


gene_heatmap <- c(gene_markers,top10$gene)
combined_sct <- ScaleData(combined_sct,features = gene_heatmap)

p1 <- DoHeatmap(object = combined_sct, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Mac_Ctrl.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = combined_sct, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Mac_ctrl_raw.pdf",height = 40, width = 20)
p2
dev.off()

AAM_Ccl8 <- c(2)
AAM_Gpnmb <- c(3,18)
Mac_Cd163 <- c(0,5,8,9,11,12,15,16,19,20,21,24,25,27)
Mac_Ccr2 <- c(1,4,7,10,14,17,23,26)
Mac_TLF <- c(13)
Mac_IFN <- c(6,22)

combined_sct$group <- paste0(combined_sct$Dataset,"_",combined_sct$age)
p <- FeaturePlot(combined_sct,features = c("Cd163","Ccl8","Trem2","Gpnmb"),split.by = "group")
ggsave(filename = "feature_plot.pdf",plot = p,width = 12,height = 12)
names(AAM_Ccl8) <- rep("AAM_Ccl8",times = length(AAM_Ccl8)) 
names(AAM_Gpnmb) <- rep("AAM_Gpnmb",times = length(AAM_Gpnmb)) 
names(Mac_Cd163) <- rep("Mac_Cd163",times = length(Mac_Cd163)) 
names(Mac_Ccr2) <- rep("Mac_Ccr2",times = length(Mac_Ccr2)) 
names(Mac_TLF) <- rep("Mac_TLF",times = length(Mac_TLF)) 
names(Mac_IFN) <- rep("Mac_IFN",times = length(Mac_IFN)) 

all_idents <- c(Mac_Cd163,Mac_Ccr2,Mac_TLF,Mac_IFN,AAM_Ccl8,AAM_Gpnmb)
all_idents <- sort(all_idents)
old_idents <- c(0:27)
new_idents <- names(all_idents)
combined_sct$Celltype_integration <- combined_sct$seurat_clusters
combined_sct$Celltype_integration <- plyr::mapvalues(x = combined_sct$Celltype_integration, from = old_idents, to = new_idents)

combined_sct$Celltype_integration <- factor(combined_sct$Celltype_integration,
                                            levels = c("Mac_Cd163","Mac_Ccr2","Mac_TLF","Mac_IFN","AAM_Ccl8","AAM_Gpnmb"))

table(combined_sct$Celltype_integration)
DimPlot(combined_sct, reduction = "umap", group.by = "Celltype_integration",split.by = "Dataset",cols = cols_list,label = T)+ coord_fixed()

p <- DimPlot(combined_sct, reduction = "umap", split.by = "Dataset",ncol = 2)+ coord_fixed()
p
ggsave(plot = p2,filename = "UMAP_integrated_ctrl.pdf")
p <- DimPlot(combined_sct, reduction = "umap", split.by = "Dataset",ncol = 2)+ coord_fixed()
p


mac_Wang <- subset(combined_sct,Dataset %in% c("Wang,et al."))
DimPlot(mac_Wang, reduction = "umap", group.by = "Cell.types",split.by = "Dataset",cols = cols_list,label = T)+ coord_fixed()
DimPlot(mac_Wang, reduction = "umap", group.by = "Celltype_integration",split.by = "age",cols = cols_list,label = T)+ coord_fixed()

saveRDS(combined_sct,file = "combined_sct_mac.rds")

mac_Wang <- subset(combined_sct,Dataset %in% c("Wang,et al."))


DefaultAssay(mac_Wang) <- "RNA"
mac_Wang <- NormalizeData(mac_Wang)
mac_Wang <- FindVariableFeatures(mac_Wang, nfeatures = 3000)
mac_Wang <- ScaleData(mac_Wang)
mac_Wang <- RunPCA(mac_Wang, npcs = 50)

mac_krasniewski <- subset(combined_sct,Dataset %in% c("krasniewski,et al."))
DefaultAssay(mac_krasniewski) <- "RNA"
mac_krasniewski <- NormalizeData(mac_krasniewski)
mac_krasniewski <- FindVariableFeatures(mac_krasniewski, nfeatures = 3000)
mac_krasniewski <- ScaleData(mac_krasniewski)

var_genes <- intersect(VariableFeatures(mac_krasniewski), 
                       VariableFeatures(mac_Wang))
length(var_genes)

pancreas.anchors <- FindTransferAnchors(reference = mac_Wang, query = mac_krasniewski, dims = 1:30,
                                        reference.reduction = "pca",features = var_genes)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = mac_Wang$Cell.types, dims = 1:30)
mac_krasniewski <- AddMetaData(mac_krasniewski, metadata = predictions)

mac_Wang <- RunUMAP(mac_Wang, dims = 1:30, reduction = "pca", return.model = TRUE)
mac_krasniewski <- MapQuery(anchorset = pancreas.anchors, reference = mac_Wang, query = mac_krasniewski,
                         refdata = list(celltype = "Cell.types"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(mac_Wang, reduction = "umap", group.by = "Cell.types", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(mac_krasniewski, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p <- p1 + p2

ggsave(filename = "query.pdf",plot = p,width = 10,height = 5)


old_idents <- c("Mac_Cd163","Mac_Ccr2","Mac_TLF","Mac_IFN","AAM_Ccl8","AAM_Gpnmb")
new_idents <- c("Mac_1","Mac_2","Mac_3","Mac_4","AAM_5","AAM_6")
combined_sct$Celltype_integration <- plyr::mapvalues(x = combined_sct$Celltype_integration, from = old_idents, to = new_idents)

#combined_sct$Celltype_integration <- combined_sct$seurat_clusters

metadata <- combined_sct@meta.data[,c("Celltype_integration","Cell.types","Dataset")]
head(metadata)
# 假设您的 Seurat 对象名为 combined
colnames(metadata) <- c("Celltype_integration","original_annotation","Dataset")

# 按数据集、原始注释、整合标签分组计数
plot_sum <- metadata %>%
  group_by(Dataset, original_annotation, Celltype_integration) %>%
  summarise(Freq = n(), .groups = "drop")
#plot_sum_filtered <- plot_sum %>% filter(Freq >= 10)
# 查看统计结果
#head(plot_sum_filtered)

# 1. 数据集的自定义顺序（facet 顺序）
dataset_order <- c("Wang, et al.", "krasniewski,et al.")
plot_sum$Dataset <- factor(plot_sum$Dataset,levels = dataset_order)

integ_levels <- unique(plot_sum$Celltype_integration)

# 如果颜色数量少于整合类别数，使用 colorRampPalette 自动扩展

integ_cols <- c(cols_list,"#143268")
#names(integ_cols) <- integ_levels

library(ggalluvial)
library(ggplot2)
library(dplyr)

# 绘图：fill 映射到 Celltype_integration
p_facet <- ggplot(plot_sum,
                  aes(axis1 = original_annotation, axis2 = Celltype_integration, y = Freq)) +
  # 流线按整合类别着色
  geom_alluvium(aes(fill = Celltype_integration), width = 1/12, alpha = 0.8) +
  # 节点也按整合类别着色（若不想则改为 fill = "grey90"）
  geom_stratum(aes(fill = Celltype_integration), width = 1/12) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  scale_x_discrete(limits = c("Original Annotation", "Integrated Cluster"),
                   expand = c(0.05, 0.05)) +
  facet_wrap(~ Dataset, scales = "free_y", ncol = 2) +
  # 应用自定义颜色
  scale_fill_manual(values = integ_cols) +
  theme_minimal() +
  theme(legend.position = "none",   # 若需图例，改为 "right" 或 "bottom"
        axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Comparison of Original Annotations vs Integrated Clusters",
       y = "Number of cells")

print(p_facet)
ggsave(plot = p_facet,filename = "plot_integrated_clusters.pdf",height = 20,width = 10)

