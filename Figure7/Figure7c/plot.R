workingDir <- "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/"
setwd(workingDir)
library(ClusterGVis)
library(ggplot2)
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/DRG_Neuron_aging_2_24_final.Rdata")
library(Seurat)
DefaultAssay(DRG_Neuron_aging) <- "RNA"
pdf(file = "Celltype_markers_2.pdf",width = 4,height = 5)
features_marker <- c("Ackr1")
DotPlot(DRG_Neuron_aging, features = features_marker,group.by = "Celltype_3", cols = c("lightgrey", "red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.3, face = "italic"))
dev.off()