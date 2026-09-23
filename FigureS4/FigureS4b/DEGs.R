options(stringsAsFactors = F)
setwd(dir = "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/")
library(Seurat)
library(dplyr)
library(stringr)
library(MAST)
library(DESeq2)

load("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")
DefaultAssay(DRG_aging_all) <- "RNA"
table(DRG_aging_all$age)

DRG_aging_all_24MO <- DRG_aging_all
clusters <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                  "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9",
                  "Satellite","Schwann_N","Schwann_M","Fibroblast","VEC","VSMC",
                  "Macrophage", "Monocyte", "Neutrophil","B Cell",  "T Cell")
# clusters <- c("Satellite","Schwann_N","Schwann_M","Fibroblast","VEC","VSMC",
#               "Macrophage", "Monocyte", "Neutrophil","B Cell",  "T Cell")
#MAST
group <- c("3 MO","24 MO")
Idents(DRG_aging_all_24MO) <- "age"
DEGs_all <- as.data.frame(matrix(NA,ncol = 8,nrow = 0))
colnames(DEGs_all) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","group","gene")
for (i in 1:length(clusters)){
  #i=17 Satellite
  subset_seurat <- subset(DRG_aging_all_24MO,subset = Celltype_2_3 %in% clusters[i])
  subset_seurat_1 <- subset(subset_seurat, subset = age %in% c("3 MO","24 MO"))
  print(clusters[i])
  print(group)
  #subset_seurat_1[["RNA"]]@counts<-as.matrix(subset_seurat_1[["RNA"]]@counts)+1
  DEGs_1 <- FindMarkers(subset_seurat_1,ident.1 = group[2],ident.2 = group[1],test.use = "MAST",min.pct = 0,logfc.threshold = 0)
  DEGs_1$clusters <- clusters[i]
  DEGs_1$group <- paste0(group[1],group[2])
  DEGs_1$gene <- rownames(DEGs_1)
  DEGs_all <- rbind(DEGs_all,DEGs_1)
}
saveRDS(DEGs_all,file = "DEGs_24MO_MAST.rds")

#write.csv(DEGs_all,file = "DEGs_all_24mo_MAST.csv")
#DEGs_2_24 <- 
DEGs_2_24 <- DEGs_all
DEGs_2_24 <- subset(DEGs_2_24, pct.1 > 0.2 | pct.2 > 0.2)
DEGs_2_24 <- DEGs_2_24[!grepl("Malat1",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("^mt-",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Rp[^(s)]",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Rp[^(l)]",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Hb[^(P)]",DEGs_2_24$gene),]
#DEGs_all_2_24_MAST <- read.csv(file = "DEGs_all_24mo_MAST.csv")
DEGs_up_2_24_MAST <- DEGs_2_24[which(DEGs_2_24$avg_log2FC >= 0.7 & DEGs_2_24$p_val_adj <= 0.05),]      # 表达量显著上升的基因
DEGs_up_2_24_MAST$regulation <- "up"
DEGs_down_2_24_MAST <- DEGs_2_24[which(DEGs_2_24$avg_log2FC <= -0.7 & DEGs_2_24$p_val_adj <= 0.05),]    # 表达量显著下降的基因
DEGs_down_2_24_MAST$regulation <- "down"
DEGs_filter_2_24_MAST <- rbind(DEGs_up_2_24_MAST, DEGs_down_2_24_MAST)
saveRDS(DEGs_filter_2_24_MAST,file = "DEGs_24MO_filter_2_24_MAST.rds")

