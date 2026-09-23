library(Seurat)
library(ggplot2)
library(dplyr)
setwd("E:/web/")
options(future.globals.maxSize = 18 * 1024^3) 
datalist <- readRDS(file = "DRG_neuron_ctrl_all.rds")

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
                               anchor.features = gene_use,
                               verbose = FALSE)

# ==================== 步骤4：查找锚点 ====================
# 指定 normalization.method = "SCT" 告知函数使用已生成的SCT残差作为输入[reference:2]
# reduction = "cca" 使用典型相关分析降维（锚点默认方法）
cat("\n查找锚点...\n")
anchors <- FindIntegrationAnchors(object.list = datalist,
                                  anchor.features = gene_use,
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
combined_sct <- FindClusters(combined_sct, resolution = 0.8)

# ==================== 步骤7：可视化整合效果 ====================
dataset_order <- c("Usoskin,et al.", "Li,et al.", "Zeisel, et al.", "Wang2026,et al.", 
                   "Ginty,et al.", "Renthal,et al.", "Jung,et al.", "Wang, et al.")
combined_sct$Dataset <- factor(combined_sct$Dataset,levels = (dataset_order))
p1 <- DimPlot(combined_sct, reduction = "umap", group.by = "Dataset",order = rev(dataset_order),pt.size = 0.1) + 
  ggtitle("SCT整合后：各数据集分布") + coord_fixed()
p1
ggsave(plot = p1,width = 15,height = 10,filename = "plot/UMPA_integrated_2.pdf")
p2 <- DimPlot(combined_sct, reduction = "umap", label = TRUE) + 
  ggtitle("SCT整合后：聚类结果") + coord_fixed()
p2
print(p1 + p2)

load("E:/web/cols_list.Rdata")
DimPlot(combined_sct, reduction = "umap", group.by = "celltype_2",split.by = "Dataset",cols = cols_list,ncol = 3)+ coord_fixed()

p <- DimPlot(combined_sct, reduction = "umap", split.by = "Dataset",ncol = 3)+ coord_fixed()
p
ggsave(plot = p2,filename = "UMAP_integrated_ctrl.pdf")
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
pdf(file = "heatmap_DRG_Neuron_Ctrl.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = combined_sct, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron_ctrl_raw.pdf",height = 40, width = 20)
p2
dev.off()

Inte1 <- c(3,15)
Inte2 <- c(1,20)
Inte3 <- c(21)
Inte4 <- c(12)
Inte5 <- c(22)
Inte6 <- c(5)
Inte7 <- c(2)
Inte8 <- c(4)
Inte9 <- c(9)
Inte10 <- c(19)
Inte11 <- c(0,7,11,14,23)
Inte12 <- c(8)
Inte13 <- c(16)
Inte14 <- c(6)
Inte15 <- c(17)
Inte16 <- c(13)
Inte17 <- c(10)
Inte18 <- c(18)

names(Inte1) <- rep("Inte1",times = length(Inte1)) 
names(Inte2) <- rep("Inte2",times = length(Inte2)) 
names(Inte3) <- rep("Inte3",times = length(Inte3)) 
names(Inte4) <- rep("Inte4",times = length(Inte4)) 
names(Inte5) <- rep("Inte5",times = length(Inte5)) 
names(Inte6) <- rep("Inte6",times = length(Inte6)) 
names(Inte7) <- rep("Inte7",times = length(Inte7)) 
names(Inte8) <- rep("Inte8",times = length(Inte8)) 
names(Inte9) <- rep("Inte9",times = length(Inte9)) 
names(Inte10) <- rep("Inte10",times = length(Inte10)) 
names(Inte11) <- rep("Inte11",times = length(Inte11)) 
names(Inte12) <- rep("Inte12",times = length(Inte12)) 
names(Inte13) <- rep("Inte13",times = length(Inte13)) 
names(Inte14) <- rep("Inte14",times = length(Inte14)) 
names(Inte15) <- rep("Inte15",times = length(Inte15)) 
names(Inte16) <- rep("Inte16",times = length(Inte16)) 
names(Inte17) <- rep("Inte17",times = length(Inte17)) 
names(Inte18) <- rep("Inte18",times = length(Inte18)) 

all_idents <- c(Inte1,Inte2,Inte3,Inte4,Inte5,Inte6,Inte7,Inte8,Inte9,
                Inte10,Inte11,Inte12,Inte13,Inte14,Inte15,Inte16,Inte17,Inte18)
all_idents <- sort(all_idents)
old_idents <- c(0:23)
new_idents <- names(all_idents)
combined_sct$Celltype_integration <- combined_sct$seurat_clusters
combined_sct$Celltype_integration <- plyr::mapvalues(x = combined_sct$Celltype_integration, from = old_idents, to = new_idents)

combined_sct$Celltype_integration <- factor(combined_sct$Celltype_integration,
                                            levels = c("Inte1","Inte2","Inte3","Inte4"
                                                       ,"Inte5","Inte6","Inte7","Inte8"
                                                       ,"Inte9","Inte10","Inte11","Inte12","Inte13"
                                                       ,"Inte14","Inte15","Inte16","Inte17","Inte18"))

table(combined_sct$Celltype_integration)
# 保存整合后的对象
saveRDS(combined_sct, file = "DRG_Neuron_ctrl_SCT_integrated.rds")
combined_sct <- readRDS(file = "DRG_Neuron_ctrl_SCT_integrated.rds")
metadata <- combined_sct@meta.data[,c("Celltype_integration","celltype_2","Celltype_3","Celltype","ClusterName",
                                      "clusters_2","celltype","Level3","Dataset")]
head(metadata)

DRG_aging <- subset(combined_sct,Dataset %in% "Wang2026,et al.")
P_Wang2026 <- DimPlot(DRG_aging, reduction = "umap", group.by = "Celltype_3",cols = cols_list)+ coord_fixed()
P_Wang2026

DRG_Ginty <- subset(combined_sct,Dataset %in% "Ginty,et al.")
P_Ginty <- DimPlot(DRG_Ginty, reduction = "umap", group.by = "Celltype",cols = cols_list)+ coord_fixed()
P_Ginty

DRG_Woolf <- subset(combined_sct,Dataset %in% "Renthal,et al.")
P_Woolf <- DimPlot(DRG_Woolf, reduction = "umap", group.by = "Celltype",cols = cols_list)+ coord_fixed()
P_Woolf

DRG_mousebrain <- subset(combined_sct,Dataset %in% "Zeisel, et al.")
cols_list <- c(cols_list,"#143268")
P_mousebrain <- DimPlot(DRG_mousebrain, reduction = "umap", group.by = "ClusterName",cols = cols_list)+ coord_fixed()
P_mousebrain

DRG_Li <- subset(combined_sct,Dataset %in% "Li,et al.")
P_Li <- DimPlot(DRG_Li, reduction = "umap", group.by = "clusters_2",cols = cols_list)+ coord_fixed()
P_Li

DRG_Jung <- subset(combined_sct,Dataset %in% "Jung,et al.")
P_Jung <- DimPlot(DRG_Jung, reduction = "umap", group.by = "celltype",cols = cols_list)+ coord_fixed()
P_Jung

DRG_Usoskin <- subset(combined_sct,Dataset %in% "Usoskin,et al.")
P_Usoskin <- DimPlot(DRG_Usoskin, reduction = "umap", group.by = "Level3",cols = cols_list)+ coord_fixed()
P_Usoskin

DRG_Wang <- subset(combined_sct,Dataset %in% "Wang, et al.")
DRG_Wang$celltype_2 <- factor(DRG_Wang$celltype_2,levels = c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
  "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1")) 
P_Wang <- DimPlot(DRG_Wang, reduction = "umap", group.by = "celltype_2",cols = cols_list)+ coord_fixed()
P_Wang

library(patchwork)
combined <- (P_Wang2026 | P_Wang| P_Ginty | P_Woolf) / (P_mousebrain| P_Jung | P_Usoskin | P_Li)
ggsave("plot/UMAP_all_in_one.pdf", combined, width = 25, height = 15)


#combined_sct <- readRDS(file = "DRG_Neuron_SCT_integrated.rds")

