setwd(dir = "/paper数据整理/Figure 5/Figure 5c/")
library(SeuratObject)

newThresholds <- readRDS(file = "./newThresholds.Rds")
binaryRegulonOrder <- readRDS(file = "./4.4_binaryRegulonOrder.Rds")
binaryRegulonActivity <- readRDS(file = "./4.1_binaryRegulonActivity.Rds")
load("./DRG_Macrophage_aging_2_24.Rdata")

tSNE <- as.data.frame(Embeddings(DRG_Macrophage_aging, reduction = "umap"))
genes_plot <- as.data.frame(t(binaryRegulonActivity[c("Max (89g)","Nfatc3_extended (37g)","Nr1h2_extended (88g)",
                                                      "Maf (12g)","Zfp384 (209g)"),]))
colnames(genes_plot)<-unlist(lapply(colnames(genes_plot), FUN = function(x) {return(strsplit(x, split = " ", fixed=T)[[1]][1])}))

#tsne_plot <- data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2])
tsne_plot <- data.frame(tSNE)
tsne_plot <- tsne_plot[order(rownames(tsne_plot)), ]
genes_plot <- genes_plot[order(rownames(genes_plot)), ]
head(genes_plot)
head(tsne_plot)
#tsne_plot_filter <- tsne_plot[index,]
tsne_plot_2 <- cbind(tsne_plot,genes_plot)

tsne_plot <- tsne_plot_2[rownames(tsne_plot_2) %in% rownames(DRG_Macrophage_aging@meta.data), ]

library(ggplot2)
library(cowplot)

p1 <- ggplot(tsne_plot[order(tsne_plot$Nfatc3_extended),]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Nfatc3_extended)) +
  scale_colour_gradient(low = "lightgrey", high = "#0d4c89") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  coord_fixed()
p1

p2 <- ggplot(tsne_plot[order(tsne_plot$Max),]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Max)) +
  scale_colour_gradient(low = "lightgrey", high = "#0d4c89") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  coord_fixed()
p2

p3 <- ggplot(tsne_plot[order(tsne_plot$Nr1h2_extended),]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Nr1h2_extended)) +
  scale_colour_gradient(low = "lightgrey", high = "#0d4c89") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  coord_fixed()
p3

p4 <- ggplot(tsne_plot[order(tsne_plot$Maf),]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Maf)) +
  scale_colour_gradient(low = "lightgrey", high = "#0d4c89") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  coord_fixed()
p4

p5 <- ggplot(tsne_plot[order(tsne_plot$Zfp384),]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Zfp384)) +
  scale_colour_gradient(low = "lightgrey", high = "#0d4c89") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  coord_fixed()
p5

pdf(file = "DRG_Macrophage_aging_TF_Activity.pdf",width = 16,height = 10)
plot_grid(p1,p2,p3,p4,p5,ncol = 3)
dev.off()




