library(loomR)
library(hdf5r)
library(Seurat)
library(SeuratDisk)
# library(scattermore)
# library(patchwork)
# setwd("~/DRG_review/Data_integration_3/")
setwd("E:/web/")
# load("~/DRG_review/Wang/Rdata/DRG_Neuron_all.integrated.Rdata")
# load("~/DRG_review/Wang/Rdata/SNI28d_Neuron_seurat.Rdata")

DRG_Neuron_all <- readRDS(file = "Dataset/10x_countdata.rds")
DRG_Neuron_meta <- readRDS(file = "Dataset/10x_meta_data.rds")
DRG_Neuron_meta <- DRG_Neuron_meta[colnames(DRG_Neuron_all),]       
DRG_Neuron_meta <- subset(DRG_Neuron_meta,celltype_2%in% c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                                                          "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1"))
DRG_Neuron_meta <- subset(DRG_Neuron_meta,condition %in% "Control")
DRG_Neuron_ctrl <- DRG_Neuron_all[,colnames(DRG_Neuron_all) %in% rownames(DRG_Neuron_meta)]
DRG_Neuron_ctrl_2 <- CreateSeuratObject(DRG_Neuron_ctrl)
DRG_Neuron_ctrl_2@meta.data <- DRG_Neuron_meta[,c(4,5,7)]
DRG_Neuron_ctrl_2@meta.data$Dataset <- "Wang, et al." 
head(DRG_Neuron_ctrl_2)

table(DRG_Neuron_ctrl_2$celltype_2)
DRG_Neuron_ctrl_2
DRG_Neuron.markers_ctrl <- subset(DRG_Neuron.markers,cluster %in% c("Cldn9", "Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1", "Nppb", "Th/Fam19a4", "Mrgpra3", "Mrgpra3/Mrgprb4",
                                                                    "Mrgprd/Lpar3", "Mrgprd/Gm7271", "S100b/Wnt7a", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Smr2", "S100b/Baiap2l1"))
gene_use <- unique(DRG_Neuron.markers_ctrl$gene)

#Jung
Jung_DRG_Neuron <- readRDS(file = "Dataset/GSE201654_mousegpcynohumanDRG.rds")
Jung_mouse_DRG_Neuron <- subset(Jung_DRG_Neuron,subset = Species %in% c("Mouse"))
Jung_mouse_DRG_Neuron$Dataset <- "Jung,et al."
Jung_DRG_Neuron_count <- Jung_mouse_DRG_Neuron@assays$RNA@counts
Jung_mouse_DRG_Neuron_meta <- Jung_mouse_DRG_Neuron@meta.data[,c(1,7,9)]
Jung_mouse_DRG_Neuron <- CreateSeuratObject(Jung_DRG_Neuron_count)
Jung_mouse_DRG_Neuron@meta.data <- Jung_mouse_DRG_Neuron_meta
head(Jung_mouse_DRG_Neuron)

#Mousebrain_DRG_Neuron 
Mousebrain_DRG.loom <- Connect(filename = "Dataset/l6_r3_peripheral_sensory_neurons.loom", mode = "r")
Mousebrain_DRG_Neuron_2 <- as.Seurat(Mousebrain_DRG.loom)
Mousebrain_DRG_Neuron_2@meta.data <- Mousebrain_DRG_Neuron_2@meta.data[,c("ClusterName","Description")]
Mousebrain_DRG_Neuron_2$Dataset <- "Zeisel, et al."
head(Mousebrain_DRG_Neuron_2)

#woolf
Woolf_DRG_counts <- readRDS("Dataset/GSE154659_C57_Raw_counts_woolf.RDS")
colnames(Woolf_DRG_counts) <- make.unique(colnames(Woolf_DRG_counts))

Woolf_DRG <- CreateSeuratObject(counts = Woolf_DRG_counts)
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
head(Woolf_DRG_Neuron)

#Usoskin lab
library(readxl)
Usoskin_DRG_counts <- read.table("Dataset/Usoskin_et_al.txt",header = T,row.names = 1,sep = "\t")
Usoskin_info <- as.data.frame(t(Usoskin_DRG_counts[1:5,]))
Usoskin_DRG_counts <- Usoskin_DRG_counts[-c(1:5),]
colnames(Usoskin_info) <- c("Sex","Level1","Level2","Level3","Content")
Usoskin_info <- subset(Usoskin_info,Level3 %in% c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH"))
#rownames(Usoskin_DRG_counts) <- gene
Usoskin_DRG_counts <- Usoskin_DRG_counts[,rownames(Usoskin_info)]
Usoskin_DRG <- CreateSeuratObject(counts = Usoskin_DRG_counts)
Usoskin_DRG@meta.data <- Usoskin_info[,c(1,4)]
Usoskin_DRG$Dataset <- "Usoskin,et al."
head(Usoskin_DRG)

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

Ginty_DRG <- CreateSeuratObject(counts = Ginty_DRG_counts)
Ginty_DRG$Dataset <- "Ginty,et al."
Ginty_DRG$Celltype <- rownames(Ginty_DRG@meta.data)
Ginty_DRG$Celltype <- gsub("\\d", "", Ginty_DRG$Celltype)
Ginty_DRG$Celltype <- gsub("\\.$", "", Ginty_DRG$Celltype)
#DRG_Ginty$Celltype_Ginty <- Celltype_Ginty
neurons <- unique(Ginty_DRG$Celltype)
neurons <- as.character(neurons)
neurons <- neurons[-10]
Ginty_DRG_Neuron <- subset(Ginty_DRG, Celltype %in% neurons)

#li_paper
Li_DRG_counts <- read.table(file = "Dataset/FPKM_data_filter.txt", header = T, sep = "\t")
rownames(Li_DRG_counts) <- Li_DRG_counts[,1]
Li_DRG_counts <- Li_DRG_counts[,-1]
Li_DRG_anno <- read.table(file = "Dataset/sample_anno.txt", header = T,row.names = 1, sep = "\t")
Li_DRG_anno <- subset(Li_DRG_anno, Condition %in% "Control")
Li_DRG_ctrl <- Li_DRG_counts[,rownames(Li_DRG_anno)]
Li_DRG_Neuron <- CreateSeuratObject(counts = Li_DRG_ctrl)
Li_DRG_Neuron@meta.data <- Li_DRG_anno
Li_DRG_Neuron@meta.data$Dataset <- "Li,et al."
head(Li_DRG_Neuron)

load("E:/DRG_aging/paper_Wang&Chen2025/Figure2/Figure 2b/DRG_Neuron_aging_2_24_final.Rdata")
DRG_Neuron_aging@meta.data <- DRG_Neuron_aging@meta.data[,c("age","sex","Celltype_3")]
DRG_Neuron_aging@meta.data$Dataset <- "Wang2026,et al."
DRG_Neuron_aging <- UpdateSeuratObject(DRG_Neuron_aging)

data_all <- merge(x = DRG_Neuron_ctrl_2, y= c(DRG_Neuron_aging,Jung_mouse_DRG_Neuron, Mousebrain_DRG_Neuron_2, Woolf_DRG_Neuron, Usoskin_DRG, Li_DRG_Neuron, Ginty_DRG_Neuron))
datalist <- SplitObject(data_all, split.by = "Dataset")
datalist <- datalist[c("Wang, et al.", "Wang2026,et al.", "Jung,et al.", "Ginty,et al.", "Li,et al.",
                       "Renthal,et al.","Usoskin,et al.","Zeisel, et al.")]
saveRDS(datalist,file = "DRG_neuron_ctrl_all.rds")
