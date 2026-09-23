setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")

ident_levels <- c("Neuron","Satellite","Schwann_M","Schwann_N","Fibroblast","VEC","VSMC","Macrophage","Monocyte", "Neutrophil","B Cell", "T Cell")
DRG_aging_all$Celltype_2_2 <- factor(DRG_aging_all$Celltype_2_2,levels = ident_levels)
color_all <- brewer.pal(12,"Paired")
names(color_all) <- c("VEC","Fibroblast","Schwann_M","Schwann_N","Satellite","Neuron","VSMC","Neutrophil","Monocyte","Macrophage","B Cell","T Cell")

pdf(file = "Vlnplot_qc_all_1.pdf",width = 10,height = 5)
VlnPlot(object = DRG_aging_all, features = c("nFeature_RNA"),cols = color_all,
        pt.size = 0,group.by = "Celltype_2_2") + 
  geom_boxplot(width=.2,col="black")+ NoLegend()
dev.off()

pdf(file = "Vlnplot_qc_all_2.pdf",width = 10,height = 5)
VlnPlot(object = DRG_aging_all, features = c("nCount_RNA"),cols = color_all,
        pt.size = 0,group.by = "Celltype_2_2") + 
  geom_boxplot(width=.2,col="black")+ NoLegend()
dev.off()


pdf(file = "UMAP_DRG_aging_all_split_age.pdf",width = 7,height = 5)
DimPlot(object = DRG_aging_all, reduction = "umap", label = F,order = c("3 MO","24 MO"), 
        group.by = "age", repel = TRUE) + coord_fixed()
dev.off()

pdf(file = "UMAP_DRG_aging_all_batch.pdf",width = 7,height = 5)
DimPlot(object = DRG_aging_all, reduction = "umap", label = F, order = c("10X_1_batch3","10X_2_batch1","10X_2_batch2"),
        group.by = "batch_2", repel = TRUE) + coord_fixed()
dev.off()


prop.table(table(Idents(DRG_aging_all)))
table(DRG_aging_all$Celltype_2_2, DRG_aging_all$age)#各组不同细胞群细胞数
Cellratio <- prop.table(table(DRG_aging_all$Celltype_2_2, DRG_aging_all$age), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
library(RColorBrewer)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,linewidth = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values=color_all)+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))
ggsave(filename = "cellratio.pdf",width = 5,height = 3)
