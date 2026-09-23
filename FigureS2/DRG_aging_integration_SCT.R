library(Seurat)
library(ggplot2)
library(dplyr)
setwd("E:/web/")
options(future.globals.maxSize = 18 * 1024^3) 
datalist <- readRDS(file = "all_DRG_neuron.rds")

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
load("E:/web/DRG_neuron_modify.Rdata")
rm(DRG_Neuron)
DRG_Neuron.markers_ctrl <- subset(DRG_Neuron.markers,cluster %in% c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                                                                    "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1"))
gene_use <- intersect(unique(DRG_Neuron.markers_ctrl$gene),integration_features) 
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
combined_sct <- RunUMAP(combined_sct, reduction = "pca", dims = 1:30, verbose = FALSE,seed.use = 50)
combined_sct <- FindNeighbors(combined_sct, reduction = "pca", dims = 1:30)
combined_sct <- FindClusters(combined_sct, resolution = 0.5)

# ==================== 步骤7：可视化整合效果 ====================
p1 <- DimPlot(combined_sct, reduction = "umap", group.by = "Dataset") + 
  ggtitle("SCT整合后：各数据集分布") + coord_fixed()
p2 <- DimPlot(combined_sct, reduction = "umap", label = TRUE) + 
  ggtitle("SCT整合后：聚类结果") + coord_fixed()
p2
print(p1 + p2)

DimPlot(combined_sct, reduction = "umap", group.by = "celltype_2",split.by = "Dataset",cols = cols_list,ncol = 3)+ coord_fixed()

p <- DimPlot(combined_sct, reduction = "umap", split.by = "Dataset",ncol = 3)+ coord_fixed()
p
ggsave(plot = p2,filename = "UMAP_integrated.pdf")
DRG_Neuron.markers <- FindAllMarkers(object = combined_sct,test.use = "roc",min.pct = 0.2,only.pos = T)
saveRDS(DRG_Neuron.markers,file = "DRG_Neuron_integrated.markers.rds")
DRG_Neuron.markers <- readRDS(file = "DRG_Neuron_integrated.markers.rds")
library(dplyr)
top10 <- DRG_Neuron.markers %>% group_by(cluster) %>% top_n(20, myAUC)
gene_markers <- c("Gfra3","Cldn9","Oprk1","Zcchc12","Sstr2","Dcn","Gabrg3","Trpm8","Rxfp1","Nppb","Th","Tafa4","Mrgpra3","Mrgprb4","Mrgprd","Lpar3","Gm7271","S100b","Pvalb","Wnt7a","Trappc3l","Ikzf1","Calb1","Bmpr1b","Prokr2","Bmpr1b","Smr2","Asic3","Baiap2l1","Colq","Ntrk2")

gene_heatmap <- c(gene_markers,top10$gene)
combined_sct <- ScaleData(combined_sct,features = gene_heatmap)

p1 <- DoHeatmap(object = combined_sct, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = combined_sct, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron_raw.pdf",height = 40, width = 20)
p2
dev.off()


# 保存整合后的对象
saveRDS(combined_sct, file = "DRG_Neuron_SCT_integrated.rds")

DRG_aging <- subset(combined_sct,Dataset %in% "Wang2026, et al")

P_Wang2026 <- DimPlot(DRG_aging, reduction = "umap", group.by = "Celltype_3",cols = cols_list)+ coord_fixed()
P_Wang2026
ggsave(filename = "UMAP_Wang.pdf",plot = P_Wang2026)
DRG_Ginty <- subset(combined_sct,Dataset %in% "Ginty,et al.")
DRG_Ginty$Celltype_Ginty <- rownames(DRG_Ginty@meta.data)
Celltype_Ginty <- gsub("\\d", "", DRG_Ginty$Celltype_Ginty)
Celltype_Ginty <- gsub("\\.$", "", Celltype_Ginty)
DRG_Ginty$Celltype_Ginty <- Celltype_Ginty
P_Ginty <- DimPlot(DRG_Ginty, reduction = "umap", group.by = "Celltype_Ginty",cols = cols_list)+ coord_fixed()
P_Ginty
ggsave(filename = "UMAP_Ginty.pdf",plot = P_Ginty)

DRG_Woolf <- subset(combined_sct,Dataset %in% "Renthal,et al.")
P_Woolf <- DimPlot(DRG_Woolf, reduction = "umap", group.by = "Celltype",cols = cols_list)+ coord_fixed()
P_Woolf

DRG_mousebrain <- subset(combined_sct,Dataset %in% "Zeisel, et al.")
cols_list <- c(cols_list,"#143268")
P_mousebrain <- DimPlot(DRG_mousebrain, reduction = "umap", group.by = "ClusterName",cols = cols_list)+ coord_fixed()
P_mousebrain

DRG_Li <- subset(combined_sct,Dataset %in% "Li,et al.")
DRG_Li$Celltype_Li <- Li_DRG_anno$clusters_2
P_Li <- DimPlot(DRG_Li, reduction = "umap", group.by = "Celltype_Li")+ coord_fixed()
P_Li

library(patchwork)
combined <- (P_Wang2026 | P_Ginty | P_Woolf) / (P_mousebrain | P_Li)
ggsave("UMAP_all_in_one.pdf", combined, width = 14, height = 10)

Inte1 <- c(4)
Inte2 <- c(2,19)
Inte3 <- c(9)
Inte4 <- c(10)
Inte5 <- c(22)
Inte6 <- c(7,25)
Inte7 <- c(1,13)
Inte8 <- c(21)
Inte9 <- c(8)
Inte10 <- c(20)
Inte11 <- c(3,11,17)
Inte12 <- c(0,5)
Inte13 <- c(15)
Inte14 <- c(6,23)
Inte15 <- c(16)
Inte16 <- c(14,24)
Inte17 <- c(12)
Inte18 <- c(18)
names(vGLUT2) <- rep("vGLUT2",times = length(vGLUT2)) 
names(vGAT) <- rep("vGAT",times = length(vGAT))
all_idents <- c(vGLUT2,vGAT)
all_idents <- sort(all_idents)
old_idents <- c(0:22)
new_idents <- names(all_idents)
SC_GFP_neuron_all$Celltype_1 <- SC_GFP_neuron_all$seurat_clusters
SC_GFP_neuron_all$Celltype_1 <- plyr::mapvalues(x = SC_GFP_neuron_all$Celltype_1, from = old_idents, to = new_idents)

DRG_Neuron@active.ident <- DRG_Neuron$Celltype
DRG_Neuron@active.ident <- plyr::mapvalues(x = DRG_Neuron@active.ident, from = ident_levels, to = ident_levels_2)
table(DRG_Neuron@active.ident)
DRG_Neuron@active.ident <- factor(DRG_Neuron@active.ident,levels = (ident))


combined_sct <- readRDS(file = "DRG_Neuron_SCT_integrated.rds")

