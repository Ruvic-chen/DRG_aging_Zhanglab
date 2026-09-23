setwd(dir = "C:/Users/Ruvic/Desktop/Figure 6c")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)

load("./Cellchat_DRG_24MO_all.Rdata")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

pathways.show <- c("CCL","IGF","GRN","CD39","PSAP","GAS")
data <- CellChatDB.use$interaction
LRpair <- subset(data, pathway_name %in% pathways.show)
LRpair <- as.data.frame(LRpair[,1])
colnames(LRpair) <- "interaction_name"
LRpair <- LRpair %>% arrange(desc(interaction_name))

pdf("Ccl8-Neuron-bubble.pdf",width = 6,height = 3)
netVisual_bubble(cellchat, sources.use = 27, targets.use = c(1:16),pairLR.use = LRpair,remove.isolate = FALSE)
dev.off()