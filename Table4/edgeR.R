library(Libra)
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(Seurat)
setwd("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/edgeR/")

load("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")

DefaultAssay(DRG_aging_all) <- "RNA"

table(DRG_aging_all$Celltype_2_3)

DE <- run_de(DRG_aging_all, replicate_col = "batch_2",cell_type_col = "Celltype_2_3",label_col = "age")

head(DE)

DEGs_2_24 <- DE
DEGs_2_24 <- subset(DEGs_2_24, `3 MO.pct` > 0.2 | `24 MO.pct` > 0.2)
DEGs_2_24 <- DEGs_2_24[!grepl("Malat1",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("^mt-",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Rp[^(s)]",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Rp[^(l)]",DEGs_2_24$gene),]
DEGs_2_24 <- DEGs_2_24[!grepl("Hb[^(P)]",DEGs_2_24$gene),]

DEGs_up_2_24_edgeR <- DEGs_2_24[which(DEGs_2_24$avg_logFC >= 0.2 & DEGs_2_24$p_val_adj <= 0.05),]      # 表达量显著上升的基因
DEGs_up_2_24_edgeR$regulation <- "up"
DEGs_down_2_24_edgeR <- DEGs_2_24[which(DEGs_2_24$avg_logFC <= -0.2 & DEGs_2_24$p_val_adj <= 0.05),]    # 表达量显著下降的基因
DEGs_down_2_24_edgeR$regulation <- "down"
DEGs_filter_2_24_edgeR <- rbind(DEGs_up_2_24_edgeR, DEGs_down_2_24_edgeR)
saveRDS(DEGs_filter_2_24_edgeR,file = "DEGs_24MO_filter_2_24_edgeR.rds")

