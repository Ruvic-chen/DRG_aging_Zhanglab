setwd("/paper????????/Figure 2/Figure 2d/")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("/paper????????/Figure 2/Figure 2a/DRG_aging_all_2_24_final.Rdata")

pdf(file = "Vlnplot_ncount_celltype_2_3.pdf",width = 10,height = 3)
VlnPlot(object = DRG_aging_all, features = c("nFeature_RNA"),
        pt.size = 0,group.by = "Celltype_2_3",cols = rep(x = "white",32)) + 
  geom_boxplot(width=.2,col="black")+ NoLegend()
dev.off()