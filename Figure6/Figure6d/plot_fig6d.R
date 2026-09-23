workingDir <- "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/"
setwd(workingDir)
library(Seurat)
library(ggplot2)

load("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge/DRG_aging_all_2_24_final.Rdata")

p <- FeaturePlot(DRG_Macrophage_aging, features = c("Cd68", "Gpnmb"), blend = TRUE) + coord_fixed()
ggsave(plot= p, filename = "GPNMB_CD68.pdf",width = 12,height = 3)