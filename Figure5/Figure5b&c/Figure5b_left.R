workingDir <- "~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/Monocle3"
setwd(workingDir)

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
cds <- readRDS("cds_gpnmb.rds")
Track_genes <- read.csv(file = "Trajectory_genes_subsampling_freq_gpnmb.csv")
load("E:/DRG_aging/paper_20251226/Rdata/Fig4&5/DRG_Macrophage_aging_gpnmb.Rdata")

top10 <- Track_genes %>% top_n(-100, q_value)
genes <- unique(top10$gene_short_name)
pt.matrix <- normalized_counts(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;

df <- DRG_Macrophage_aging_gpnmb@meta.data[,c(4,19,24)]

df <- df[order(df$pseudotime),]
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

col_pseudotime <- colorRamp2(
  breaks = seq(from = 0, to = 18, by = 18/8), 
  colors = brewer.pal(9,"BuPu")
)
col_con <- brewer.pal(3,"Dark2")[-2]
names(col_con) <- c("3 MO","24 MO")

col_clusters <- brewer.pal(3,"Set1")[-2]
names(col_clusters) <- c("Mac_Ccr2","AAM_Trem2")

ha = HeatmapAnnotation(df = df[,1:3],
                       col = list(age = col_con,
                                  celltype = col_clusters,
                                  pseudotime = col_pseudotime))

htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 3,
  top_annotation = ha,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm
pdf(file = "Figure5a_right.pdf",width = 6,height = 10)
htkm
dev.off()