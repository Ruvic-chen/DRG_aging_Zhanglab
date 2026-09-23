setwd("E:\\DRG_aging\\paper_Wang&Chen2025\\Figure2\\Figure 2a")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("./DRG_aging_all_2_24_final.Rdata")

rawdata <- DRG_aging_all@assays$RNA@counts
metadata <- DRG_aging_all@meta.data[,c(1:7,18,16,19)]
metadata$barcode <- rownames(metadata)
head(metadata)
metadata$barcode <- substr(metadata$barcode, 1, 18)

table(metadata$Celltype_3)
metadata$Celltype_3 <- as.character(metadata$Celltype_3)
metadata$Celltype_3[metadata$Celltype_3 %in% "AAM_Trem2"] <- "AAM_Gpnmb"
saveRDS(metadata,file = "metadata.rds")
rawdata_neuron <- rawdata[,metadata$Celltype_2_2 %in% "Neuron"]
saveRDS(rawdata_neuron,file = "rawdata_neuron.rds")

saveRDS(rawdata,file = "rawdata.rds")

DRG_aging_all@meta.data <- metadata

ident_levels <- c("Neuron","Satellite","Schwann_M","Schwann_N","Fibroblast","VEC","VSMC","Macrophage","Monocyte", "Neutrophil","B Cell", "T Cell")
DRG_aging_all$Celltype_2_2 <- factor(DRG_aging_all$Celltype_2_2,levels = ident_levels)
color_all <- brewer.pal(12,"Paired")
names(color_all) <- c("VEC","Fibroblast","Schwann_M","Schwann_N","Satellite","Neuron","VSMC","Neutrophil","Monocyte","Macrophage","B Cell","T Cell")

pdf(file = "UMAP_DRG_aging_all_celltype_all_split.pdf",width = 17,height = 5)
DimPlot(object = DRG_aging_all, reduction = "umap", label = F, 
        cols = color_all,split.by = "age", 
        group.by = "Celltype_2_2", repel = TRUE) + coord_fixed()
dev.off()