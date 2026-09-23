library(loomR)
library(hdf5r)
library(Seurat)
library(SeuratDisk)
# library(scattermore)
# library(patchwork)
setwd("~/DRG_review/Data_integration_3/")
setwd("E:/web/")
# load("~/DRG_review/Wang/Rdata/DRG_Neuron_all.integrated.Rdata")
# load("~/DRG_review/Wang/Rdata/SNI28d_Neuron_seurat.Rdata")

load("E:/web/DRG_neuron_modify.Rdata")
rm(DRG_Neuron)
DRG_Neuron.markers_ctrl <- subset(DRG_Neuron.markers,cluster %in% c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                                                                    "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1"))
gene_use <- unique(DRG_Neuron.markers_ctrl$gene)

DRG_Neuron_all <- readRDS(file = "Dataset/10x_countdata.rds")
DRG_Neuron_meta <- readRDS(file = "Dataset/10x_meta_data.rds")
#DRG_Neuron_meta <- DRG_Neuron_meta[colnames(DRG_Neuron_all),]       
DRG_Neuron_meta <- subset(DRG_Neuron_meta,celltype_2 %in% c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                                                          "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1"))
DRG_Neuron_meta <- subset(DRG_Neuron_meta,condition %in% "Control")
DRG_Neuron_ctrl <- DRG_Neuron_all[,colnames(DRG_Neuron_all) %in% rownames(DRG_Neuron_meta)]
DRG_Neuron_ctrl_2 <- CreateSeuratObject(DRG_Neuron_ctrl)
DRG_Neuron_ctrl_2@meta.data <- DRG_Neuron_meta
DRG_Neuron_ctrl_2@meta.data$Dataset <- "Wang, et al." 
# DRG_Neuron_ctrl_2 <- NormalizeData(DRG_Neuron_ctrl_2,normalization.method = "LogNormalize",scale.factor = 10000)
# DRG_Neuron_ctrl_2 <- FindVariableFeatures(DRG_Neuron_ctrl_2,verbose = F)
table(DRG_Neuron_ctrl_2$celltype_2)

Jung_DRG_Neuron <- readRDS(file = "Dataset/GSE201654_mousegpcynohumanDRG.rds")
Jung_mouse_DRG_Neuron <- subset(Jung_DRG_Neuron,subset = Species %in% c("Mouse"))
Jung_mouse_DRG_Neuron$Dataset <- "Jung,et al."
Jung_mouse_DRG_Neuron <- UpdateSeuratObject(Jung_mouse_DRG_Neuron)

#Mousebrain_DRG_Neuron 
Mousebrain_DRG.loom <- Connect(filename = "Dataset//l6_r3_peripheral_sensory_neurons.loom", mode = "r")
Mousebrain_DRG_Neuron_2 <- as.Seurat(Mousebrain_DRG.loom)
Mousebrain_DRG_Neuron_2$Dataset <- "Zeisel, et al."
# Mousebrain_DRG_Neuron_2 <- NormalizeData(Mousebrain_DRG_Neuron_2,normalization.method = "LogNormalize",scale.factor = 10000)
# Mousebrain_DRG_Neuron_2 <- FindVariableFeatures(Mousebrain_DRG_Neuron_2,verbose = F)
#woolf

Woolf_DRG_counts <- readRDS("Dataset/GSE154659_C57_Raw_counts_woolf.RDS")
counts_unique <- Woolf_DRG_counts[, !duplicated(colnames(Woolf_DRG_counts))]
Woolf_DRG <- CreateSeuratObject(counts = counts_unique, min.cells = 3, project = "Woolf_DRG")

info <- rownames(Woolf_DRG@meta.data)
library(stringr)
info <- as.data.frame(str_split(info, "_", n=8, simplify = TRUE))
info <- info[,c(1:3,6)]
names(info)<-c("Sex","Strain","Group","Celltype")
Woolf_DRG@meta.data <- cbind(Woolf_DRG@meta.data, info)
Woolf_DRG@meta.data$Celltype[which(Woolf_DRG@meta.data$Celltype =='p')] <- 'p_cLTMR2'
Woolf_DRG_Neuron <- Woolf_DRG[,Woolf_DRG@meta.data[["Group"]] %in% "Naive" 
                              & Woolf_DRG@meta.data[["Celltype"]] %in% c("SST","NP","PEP1","NF2",
                                                                         "NF3","cLTMR1","NF1","p_cLTMR2",
                                                                         "PEP2")]
Woolf_DRG_Neuron@meta.data$Dataset <- "Renthal,et al."

#Usoskin lab
library(readxl)
Usoskin_DRG_counts <- read.table("Dataset/Usoskin_et_al.txt",header = T,row.names = 1,sep = "\t")
Usoskin_info <- as.data.frame(t(Usoskin_DRG_counts[1:5,]))
Usoskin_DRG_counts <- Usoskin_DRG_counts[-c(1:5),]
colnames(Usoskin_info) <- c("Sex","Level1","Level2","Level3","Content")
Usoskin_info <- subset(Usoskin_info,Level3 %in% c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH"))
#rownames(Usoskin_DRG_counts) <- gene
Usoskin_DRG_counts <- Usoskin_DRG_counts[,rownames(Usoskin_info)]
Usoskin_DRG <- CreateSeuratObject(counts = Usoskin_DRG_counts, project = "Usoskin_DRG",min.cells = 3,min.features = 200)
Usoskin_DRG$Dataset <- "Usoskin,et al."
# Usoskin_DRG <- NormalizeData(Usoskin_DRG)
# Usoskin_DRG <- FindVariableFeatures(Usoskin_DRG, selection.method = "vst", nfeatures = 2000)
#Usoskin_DRG <- SCTransform(Usoskin_DRG)

#Ginty lab
Ginty_DRG_counts <- read.csv(file="Dataset/GSM4130750_WT_1.csv", header = T, sep = ",")
cellid <- Ginty_DRG_counts[1,]
cellid <- cellid[,-1]
cellid <- t(cellid)
Ginty_DRG_counts <- Ginty_DRG_counts[-2,]
cellbarcode <- Ginty_DRG_counts[1,]
cellbarcode <- cellbarcode[,-1]
Ginty_DRG_counts <- Ginty_DRG_counts[-1,]
gene <- Ginty_DRG_counts[,1]
Ginty_DRG_counts <- Ginty_DRG_counts[,-1]
gene <- gene[-c(12088,12091)]
Ginty_DRG_counts <- Ginty_DRG_counts[-c(12090,12093),]
rownames(Ginty_DRG_counts) <- gene
colnames(Ginty_DRG_counts) <- colnames(cellbarcode)

Ginty_DRG <- CreateSeuratObject(counts = Ginty_DRG_counts, min.cells = 3, project = "Ginty_DRG")
Ginty_DRG$Dataset <- "Ginty,et al."
Ginty_DRG$Celltype <- cellid
neurons <- unique(Ginty_DRG$Celltype)
neurons <- as.character(neurons)
neurons <- neurons[-10]
Ginty_DRG_Neuron <- Ginty_DRG[,Ginty_DRG@meta.data[["Celltype"]] %in% neurons]

#li_paper
Li_DRG_counts <- read.table(file = "Dataset/FPKM_data_filter.txt", header = T, sep = "\t")
rownames(Li_DRG_counts) <- Li_DRG_counts[,1]
Li_DRG_counts <- Li_DRG_counts[,-1]
Li_DRG_anno <- read.table(file = "Dataset/sample_anno.txt", header = F, sep = "\t")
Li_DRG_anno <- Li_DRG_anno[,-1]
Li_DRG_Neuron <- CreateSeuratObject(counts = Li_DRG_counts, min.cells = 3, project = "Li_DRG")
Li_DRG_Neuron$Dataset <- "Li,et al."

load("E:/web/Dataset/DRG_aging_all_2_24_final.Rdata")
DRG_neuron_wang2026 <- subset(DRG_aging_all,Celltype_2_2 %in% c("Neuron"))
DRG_neuron_wang2026 <- UpdateSeuratObject(DRG_neuron_wang2026)
DRG_neuron_wang2026$Dataset <- "Wang2026, et al"

data_all <- merge(x = DRG_Neuron_ctrl_2, y= c(Jung_mouse_DRG_Neuron, Mousebrain_DRG_Neuron_2, Woolf_DRG_Neuron, Usoskin_DRG, Li_DRG_Neuron, Ginty_DRG_Neuron,DRG_neuron_wang2026))
datalist <- SplitObject(data_all, split.by = "Dataset")
datalist <- datalist[c("Wang, et al.", "Jung,et al.", "Ginty,et al.", "Li,et al.",
                       "Renthal,et al.","Usoskin,et al.","Zeisel, et al.","Wang2026, et al")]
saveRDS(datalist,file = "all_DRG_neuron.rds")
#Data integration use CCA
for (i in 1:length(datalist)) {
  datalist[[i]] <- NormalizeData(datalist[[i]], verbose = FALSE)
  datalist[[i]] <- FindVariableFeatures(datalist[[i]], 
                                        selection.method = "vst", 
                                        nfeatures = 2000,
                                        verbose = FALSE)
  cat("数据集", i, "完成\n")
}
anchors <- FindIntegrationAnchors(object.list = datalist,
                                  anchor.features = gene_use,
                                  normalization.method = "LogNormalize", # 默认可省略
                                  dims = 1:30,
                                  verbose = TRUE)
#anchors <- FindIntegrationAnchors(object.list = datalist, anchor.features = gene_use)

DRG_Neuron.combined_CCA <- IntegrateData(anchorset = anchors,
                                         dims = 1:30,
                                         features.to.integrate = gene_use,
                                         verbose = TRUE)
saveRDS(DRG_Neuron.combined_CCA, file = "DRG_Neuron_integrated.rds")
# DRG_Neuron.combined_CCA <- IntegrateData(anchorset = Anchors)
DefaultAssay(DRG_Neuron.combined_CCA) <- "integrated"

DRG_Neuron.combined_CCA <- ScaleData(DRG_Neuron.combined_CCA, verbose = FALSE)
DRG_Neuron.combined_CCA <- RunPCA(object = DRG_Neuron.combined_CCA, verbose = FALSE)
DRG_Neuron.combined_CCA <- RunUMAP(object = DRG_Neuron.combined_CCA, dims = 1:30, verbose = FALSE)
DRG_Neuron.combined_CCA <- RunTSNE(object = DRG_Neuron.combined_CCA, check_duplicates = FALSE, dims = 1:30, verbose = FALSE)
DRG_Neuron.combined_CCA <- FindNeighbors(object = DRG_Neuron.combined_CCA, dims = 1:30, verbose = FALSE)
DRG_Neuron.combined_CCA <- FindClusters(object = DRG_Neuron.combined_CCA, resolution = 0.6, verbose = FALSE)
library(ggplot2)
pdf(file = "UMAP_DRG_Neuron_combined_CCA.pdf",width = 10,height = 8)
DimPlot(object = DRG_Neuron.combined_CCA, reduction = "umap", label = TRUE,
        repel = TRUE) + coord_fixed()
dev.off()

pdf(file = "UMAP_DRG_Neuron_combined_Dataset_CCA.pdf",width = 10,height = 8)
DimPlot(object = DRG_Neuron.combined_CCA, reduction = "umap", label = TRUE, group.by = "Dataset",
        repel = TRUE) + coord_fixed()
dev.off()

save(DRG_Neuron, Mousebrain_DRG_Neuron, Woolf_DRG_Neuron, Usosokin_DRG_Neuron, Li_DRG_Neuron, Ginty_DRG_Neuron,
     file = "DRG_Neuron_mus_datasets.rds")

#Data integration use SCTransform
datalist <- lapply(X = datalist, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = datalist)
gene_use_2 <- intersect(gene_use,rownames(DRG_Neuron.combined_SCT))
datalist <- PrepSCTIntegration(object.list = datalist, anchor.features = gene_use_2)
Anchors <- FindIntegrationAnchors(object.list = datalist, normalization.method = "SCT",
                                         anchor.features = gene_use_2)
DRG_Neuron.combined_SCT <- IntegrateData(anchorset = Anchors, normalization.method = "SCT")

DRG_Neuron.combined_SCT <- ScaleData(DRG_Neuron.combined_SCT,features = features, verbose = FALSE)
DRG_Neuron.combined_SCT <- RunPCA(DRG_Neuron.combined_SCT, verbose = FALSE)
DRG_Neuron.combined_SCT <- RunUMAP(DRG_Neuron.combined_SCT, reduction = "pca", dims = 1:30)
#DRG_Neuron.combined_SCT <- RunTSNE(object = DRG_Neuron.combined_SCT, check_duplicates = FALSE, dims = 1:30, verbose = FALSE)
DRG_Neuron.combined_SCT <- FindNeighbors(object = DRG_Neuron.combined_SCT, dims = 1:30, verbose = FALSE)
DRG_Neuron.combined_SCT <- FindClusters(object = DRG_Neuron.combined_SCT, resolution = 0.6, verbose = FALSE)
library(ggplot2)
pdf(file = "UMAP_DRG_Neuron_combined_SCT.pdf",width = 10,height = 8)
DimPlot(object = DRG_Neuron.combined_SCT, reduction = "umap", label = TRUE,
        repel = TRUE) + coord_fixed()
dev.off()

order_cell <- c("Usoskin,et al.","Li,et al.","Zeisel, et al.",
                "Ginty,et al.","Renthal,et al.","Wang, et al.","Jung,et al.")
library(RColorBrewer)
pdf(file = "UMAP_DRG_Neuron_combined_Dataset_SCT.pdf",width = 8,height = 6)
DimPlot(object = DRG_Neuron.combined_SCT, reduction = "umap", label = F, group.by = "Dataset",
        repel = TRUE,order = rev(order_cell),cols = rev(brewer.pal(8,"Paired")[c(2,1,3:5,7,8)])) + 
        coord_fixed()
dev.off()
pdf(file = "UMAP_DRG_Neuron_combined_celltype_SCT.pdf",width = 15,height = 6)
DimPlot(object = DRG_Neuron.combined_SCT, reduction = "umap", label = TRUE, group.by = "Celltype",
        repel = TRUE) + coord_fixed()
dev.off()

FeaturePlot(DRG_Neuron.combined_SCT,features = c("Mrgprd","Gm7271","Lpar3"))

DRG_Neuron.markers <- FindAllMarkers(object = DRG_Neuron.combined_SCT,test.use = "roc",min.pct = 0.2,only.pos = T)
library(dplyr)
top10 <- DRG_Neuron.markers %>% group_by(cluster) %>% top_n(20, myAUC)
gene_markers <- c("Gfra3","Cldn9","Oprk1","Zcchc12","Sstr2","Dcn","Gabrg3","Trpm8","Rxfp1","Nppb","Th","Fam19a4","Mrgpra3","Mrgprb4","Mrgprd","S100b","Pvalb","Wnt7a","Trappc3l","Ikzf1","Calb1","Bmpr1b","Prokr2","Bmpr1b","Smr2","Asic3","Baiap2l1","Colq","Ntrk2")

gene_heatmap <- c(gene_markers,top10$gene)
DRG_Neuron.combined_SCT <- ScaleData(DRG_Neuron.combined_SCT,features = gene_heatmap)

p1 <- DoHeatmap(object = DRG_Neuron.combined_SCT, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = DRG_Neuron.combined_SCT, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron_raw.pdf",height = 40, width = 20)
p2
dev.off()


ident_levels <- c("Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Zcchc12/Sstr2", "Zcchc12/Dcn",  "Cldn9", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                  "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Baiap2l1", "S100b/Prokr2", "S100b/Smr2", "S100b/Ntrk3/Gfra1")
ident_levels_2 <- c("Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Zcchc12/Sstr2", "Zcchc12/Dcn",  "Cldn9", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                    "Mrgprd", "Mrgprd", "S100b/Wnt7a", "S100b/Baiap2l1", "S100b/Prokr2", "S100b/Smr2", "S100b/Ntrk3/Gfra1")
ident<- c("Zcchc12/Sstr2", "Zcchc12/Dcn","Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Cldn9", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
          "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Baiap2l1", "S100b/Prokr2", "S100b/Smr2", "S100b/Ntrk3/Gfra1")
DRG_Neuron@active.ident <- DRG_Neuron$Celltype
DRG_Neuron@active.ident <- plyr::mapvalues(x = DRG_Neuron@active.ident, from = ident_levels, to = ident_levels_2)
table(DRG_Neuron@active.ident)
DRG_Neuron@active.ident <- factor(DRG_Neuron@active.ident,levels = (ident))

select_gene <- c("Tac1","Calca","Pvalb","Nefh","S100b","Gfra1","Gfra2","Gfra3","Ntrk1","Ntrk2","Ntrk3","Sstr2","Gabrg3","Trpm8","Rxfp1","Oprk1","Nppb","Fam19a4","Mrgpra3","Mrgprb4","Mrgprd","Lpar3","Gm7271","Wnt7a","Cox6a2","Prokr2","Smr2","Ptgfr")
#select_gene <- unique(subset(DRG_Neuron.markers_ctrl,cluster %in% "Cldn9")$gene)[151:202]
p1 <- DotPlot(object = DRG_Neuron, features = rev(select_gene), cols = c("grey","red"),assay = "RNA")
p1 <- p1 + theme(axis.text.x = element_text(face="italic"),axis.text.y = element_text(face="italic"))
p1 <- p1 + theme(axis.text.x = element_text(angle = 60, hjust = 1))
p1 <- p1 + coord_flip()
pdf(file = "DotPlot_DRG_neuron_ctrl_select_2.pdf", width = 7, height = 7)
p1
dev.off()

library(harmony)
#Data integration use Harmony
DRG_Neuron <- subset(DRG_Neuron_all.integrated,condition %in% "Control")
DefaultAssay(DRG_Neuron) <- "RNA"
DRG_Neuron$Dataset <- "Wang, et al. Cell Res 2021"
DRG_Neuron$Celltype <- DRG_Neuron$celltype
DRG_Neuron_all <- merge(x = DRG_Neuron, y = c(Mousebrain_DRG_Neuron, Woolf_DRG_Neuron, Usosokin_DRG_Neuron, Li_DRG_Neuron, Ginty_DRG_Neuron))
DRG_Neuron_all$Datset2 <- paste0(DRG_Neuron_all$batch,"_",DRG_Neuron_all$Dataset)
DRG_Neuron_all <- NormalizeData(DRG_Neuron_all, normalization.method = "LogNormalize", scale.factor = 10000)
DRG_Neuron_all <- FindVariableFeatures(DRG_Neuron_all, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(DRG_Neuron_all)
DRG_Neuron_all <- ScaleData(DRG_Neuron_all , features = gene_use, vars.to.regress = "nCount_RNA")
#DRG_Neuron_all <- RunPCA(DRG_Neuron_all, features = VariableFeatures(object = DRG_Neuron_all))
DRG_Neuron_all <- RunPCA(DRG_Neuron_all, features = gene_use)
DRG_Neuron.combined_Har <- RunHarmony(DRG_Neuron_all, group.by.vars = "Datset2" , plot_convergence = F,dims.use = 1:20)
#DRG_Neuron.combined_Har = RunTSNE(DRG_Neuron.combined_Har, reduction = "harmony", dims = 1:30)
DRG_Neuron.combined_Har = RunUMAP(DRG_Neuron.combined_Har, reduction = "harmony", dims = 1:30,seed.use = 4)
DRG_Neuron.combined_Har = FindNeighbors(DRG_Neuron.combined_Har, reduction = "harmony",dims = 1:30)
DRG_Neuron.combined_Har = FindClusters(DRG_Neuron.combined_Har, resolution = 0.6)

pdf(file = "UMAP_DRG_Neuron_combined_Har.pdf",width = 10,height = 8)
DimPlot(object = DRG_Neuron.combined_Har, reduction = "umap", label = TRUE,
        repel = TRUE) + coord_fixed()
dev.off()

order_cell <- c("Usoskin, et al. Nat Neurosci 2015","Li, et al. Cell Res. 2016","Zeisel, et al. Cell 2018",
                "Sharma,et al. Nature 2020","Renthal,et al. Neuron 2020","Wang, et al. Cell Res 2021")
pdf(file = "UMAP_DRG_Neuron_combined_Dataset_Har.pdf",width = 10,height = 8)
DimPlot(object = DRG_Neuron.combined_Har, reduction = "umap", label = F, group.by = "Dataset",
        repel = TRUE,order = (order_cell),cols = brewer.pal(7,"Paired")[c(1:5,7)]) + coord_fixed()
dev.off()

pdf(file = "UMAP_DRG_Neuron_combined_Dataset_Har_celltype.pdf",width = 30,height = 8)
DimPlot(object = DRG_Neuron.combined_Har, reduction = "umap", label = TRUE, group.by = "Celltype",
        repel = TRUE) + coord_fixed()
dev.off()

DRG_Neuron.markers <- FindAllMarkers(object = DRG_Neuron.combined_Har,test.use = "roc",min.pct = 0.2,only.pos = T)
top10 <- DRG_Neuron.markers %>% group_by(cluster) %>% top_n(20, myAUC)
gene_markers <- c("Gfra3","Cldn9","Oprk1","Zcchc12","Sstr2","Dcn","Gabrg3","Trpm8","Rxfp1","Nppb","Th","Fam19a4","Mrgpra3","Mrgprb4","Mrgprd","S100b","Pvalb","Wnt7a","Trappc3l","Ikzf1","Calb1","Bmpr1b","Prokr2","Bmpr1b","Smr2","Asic3","Baiap2l1","Colq","Ntrk2")

gene_heatmap <- c(gene_markers,top10$gene)
DRG_Neuron.combined_Har <- ScaleData(DRG_Neuron.combined_Har,features = gene_heatmap)

p1 <- DoHeatmap(object = DRG_Neuron.combined_Har, features = gene_heatmap, slot = "scale.data") 
p1 <- p1 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron.pdf",height = 40, width = 20)
p1
dev.off()

p2 <- DoHeatmap(object = DRG_Neuron.combined_Har, features = gene_heatmap, slot = "data") 
p2 <- p2 + scale_fill_gradient2('legend name', low = 'blue', high = 'red', mid = 'white') 
pdf(file = "heatmap_DRG_Neuron_raw.pdf",height = 40, width = 20)
p2
dev.off()
