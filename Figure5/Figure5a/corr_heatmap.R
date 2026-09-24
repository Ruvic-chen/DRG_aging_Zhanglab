workingDir <- "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/"
setwd(workingDir)
library(Matrix)
library(Seurat)

DRG_Macrophage_aging_anno <- DRG_Macrophage_aging@meta.data
head(DRG_Macrophage_aging_anno)

DRG_Macrophage_aging_data <- DRG_Macrophage_aging@assays$RNA@data

ident_2 <-names(table(DRG_Macrophage_aging_anno$Celltype_3)) 
average_2 <- as.data.frame(matrix(NA,nrow = nrow(DRG_Macrophage_aging_data),ncol = length(ident_2)))
rownames(average_2) <- rownames(DRG_Macrophage_aging_data)
colnames(average_2) <- ident_2
for(i in 1:length(ident_2)){
  ident_data <- DRG_Macrophage_aging_data[,rownames(DRG_Macrophage_aging_anno[DRG_Macrophage_aging_anno$Celltype_3 %in% ident_2[i],])]
  ave_data <- rowMeans(ident_data)
  average_2[,i] <- ave_data
}
cor_matr_2 <- cor(average_2)

write.table(cor_matr_2,file = "cor_matr.txt",sep = "\t",col.names = NA)
library(pheatmap)
library(RColorBrewer)
pdf(file = "pheatmap_corr.pdf", width = 10, height = 10)
pheatmap(cor_matr_2, cellwidth = 20, cellheight = 20, cluster_cols = T, cluster_rows = T,labels_col = F,
         show_colnames = F, display_numbers = T, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
dev.off()
