setwd(dir = "~/pysceinc/analysis/24MO/Macrophage/")
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/DRG_Macrophage_aging_2_24.Rdata")

regulon_target <- read.csv(file = "output/Step2_regulonTargetsInfo.tsv",sep = "\t",header = T)
library(dplyr)
DRG_Macrophage_aging.markers %>%
  group_by(cluster) %>%
  top_n(n = 40, wt = avg_log2FC) -> top10
marker_genes <- subset(top10,cluster %in% c("AAM_Trem2","AAM_Ccl8"))

Regulon <- c("Maf","Zfp384","Max","Nr1h2","Nfatc3")
Genes <- unique(marker_genes$gene)

regulon_AAMac_select <- subset(regulon_target, (TF %in% Regulon) & (gene %in% Genes))

write.csv(regulon_AAMac_select,file = "regulon_AAMac_select.csv",row.names = F)