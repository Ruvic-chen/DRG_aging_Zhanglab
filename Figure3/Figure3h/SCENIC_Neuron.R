setwd(dir = "~/pysceinc/analysis/24MO/Neuron_1210/")
#setwd(dir = "~/SCENIC/analysis/24MO/Neuron/")

library(SCENIC)
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/DRG_Neuron_aging_2_24_final.Rdata")
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron/cols_list.Rdata")
DRG_Neuron_24MO <- subset(DRG_Neuron_aging, age %in% "24 MO")
exprMat <- GetAssayData(DRG_Neuron_24MO, slot = "counts", assay = "RNA")
exprMat <- as.matrix(exprMat)
cellInfo <- DRG_Neuron_24MO@meta.data[,c(1:5,93)]

dim(exprMat)
head(cellInfo)

colnames(cellInfo)[c(3,6)] <- c("nGenes","CellType")
#cbind(table(cellInfo$celltype_2))

dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
idents <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2","C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9")
library(RColorBrewer)

cluster_cols <- cols_list
names(cluster_cols) <- c(idents)
colVars <- list(CellType=setNames(cluster_cols, idents))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

library(SCENIC)
library(RcisTarget)
org <- "mgi" # 
dbDir <- "~/pysceinc/database/mm10/" # RcisTarget databases location
myDatasetTitle <- "SCENIC on Aging DRG Neuron" # choose a name for your analysis
data(defaultDbNames)
dbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather") 
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

#(Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#stop here
#interestingGenes <- union(Trajectory_genes_ccl8$gene_short_name,Trajectory_genes_Gpnmb$gene_short_name)

# any missing?
#interestingGenes[which(!interestingGenes %in% genesKept)]
#genesKept <- union(interestingGenes,genesKept)

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
# rm(exprMat)

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
#runGenie3(exprMat_filtered, scenicOptions)
GRNBoost_linkList <- importArboreto(file.path("~/pysceinc/outputs/adj_24_MO_DRG_aging.tsv"))
head(GRNBoost_linkList)

GRNBoost_linkList <- GRNBoost_linkList[,c(1,2,4)]
colnames(GRNBoost_linkList) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_linkList,file = "int/1.4_GENIE3_linkList.Rds")

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
trace('runSCENIC_2_createRegulons', edit = T, where = asNamespace("SCENIC"))
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
#newThresholds <- savedSelections$thresholds #can't open shiny app
newThresholds <- aucellApp$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

nPcs <- c(10)
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 10
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/DRG_Neuron_24MO_SCENIC.loom"
export2loom(scenicOptions, exprMat)

library(SCopeLoomR)
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

exprMat_log <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
#savedSelections <- shiny::runApp(aucellApp)

print(tsneFileName(scenicOptions))

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
#AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Trnp","Cox5a")],], plots="Expression")

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf(file = "Regulon_activity.pdf",width = 10,height = 50)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pdf(file = "Regulon_binary_activity.pdf",width = 10,height = 40)
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)

scenicOptions <- readRDS(file = "int/scenicOptions.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
cellInfo <- readRDS(file = "int/cellInfo.Rds")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)
p <- plotly::ggplotly(rssPlot$plot)
p
pdf(file = "RSS_plot.pdf",width = 8,height = 20)
plotly::ggplotly(rssPlot$plot)
dev.off()

library(ggplot2)
library(reshape2)

i=16
pdf(file = paste0("./RSS_plot/RSS_plot_",idents[i],".pdf"),width = 8,height = 20)
plotRSS_oneSet(rss, setName = idents[i])
dev.off()

###
regulon_activity_network <- readRDS("int/4.1_binaryRegulonActivity.Rds")
head(regulon_activity_network)[1:5,1:5]

row_names <- rownames(regulon_activity_network)
gene_names<-unlist(lapply(row_names, FUN = function(x) {return(strsplit(x, split = " ", fixed=T)[[1]][1])})) #以.tsv作为分割，截取之前的
gene_names<-unlist(lapply(gene_names, FUN = function(x) {return(strsplit(x, split = "_", fixed=T)[[1]][1])})) #以.tsv作为分割，截取之前的
saveRDS(gene_names,file = "regulons.rds")
#regulon_activity_network <- read.csv(file = "regulon_activity_network.csv",row.names = 1)
cellInfo <- readRDS(file="int/cellInfo.Rds")
#cellInfo$subtype_2 <- factor(cellInfo$subtype_2,levels = ident_levels)
#log_Data <- regulon_activity_network[rownames(regulon_activity_network) %in% regulonSelections$corr,]
log_Data <- regulon_activity_network
#log_Data[is.na(log_Data)] <- 0
head(log_Data)[1:5,1:5]
dim(log_Data)

df <- cellInfo
head(df)
#df <- df[,c(4,2,3)]
#colnames(df) <- c("cluster", "batch", "condition")
groups <- c("3 MO","24 MO")
df_sort <- as.data.frame(matrix(NA,nrow = 0,ncol = ncol(df)))
colnames(df_sort) <- colnames(df)
for(i in 1:length(groups)){
  cells <- df[df$age == groups[i],]
  df_sort <- rbind(df_sort,cells)
}

df_sort_2 <- as.data.frame(matrix(NA,nrow = 0,ncol = ncol(df)))
colnames(df_sort_2) <- colnames(df_sort)
for(i in 1:length(idents)){
  cells <- df_sort[df_sort$CellType == idents[i],]
  df_sort_2 <- rbind(df_sort_2,cells)
}
#df_sort$condition[df_sort$condition %in% "sni0d"] <- "ctrl"
head(df_sort_2)
df_sort_2 <- as.data.frame(df_sort_2[,c(4,6)])
#rownames(df_sort_2) <- rownames(df_sort)
colnames(df_sort_2) <- c("group","CellType")
metainds_cells <- match(rownames(df_sort_2), colnames(log_Data))
missing.meta <- is.na(metainds_cells)
summary(missing.meta)
log_Data <- log_Data[,metainds_cells]
head(log_Data)[1:5,1:5]
df_sort_2$CellType <- factor(df_sort_2$CellType,levels = idents)

#df_sort_3 <- as.data.frame(t(df_sort_2))
cluster_color1 <- cluster_cols
annotation_colors = list(CellType = c(cluster_color1))
library(pheatmap)
library(RColorBrewer)
library(Cairo)
#gene_use <- read.table(file = "genes.txt",header = T,sep = "\t")
#log_Data_2 <- log_Data[rownames(log_Data) %in% gene_use$gene_name,]
pdf(file = "pheatmap_regulon_activity_coor_DRG_24MO_Neuron.pdf", width = 20, height = 80)
pheatmap(log_Data, cluster_cols = F, cluster_rows = T,annotation_col = df_sort_2,  
         labels_col = T, show_colnames = F, color = c("White","Black"))
dev.off()

regulon_target <- read.csv(file = "output/Step2_regulonTargetsInfo.tsv",sep = "\t",header = T)

DEGs_Neuron_celldeath <- read.csv("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge/DEGs/celldeath_analysis/DEGs_celldeath_select.csv")
#Regulon <- unique(TF$gene)

celldeath <- unique(DEGs_Neuron_celldeath$gene)

#Regulon <- celldeath
#Regulon <- c("Mapkapk2","Ern1","Bcl2","Fth1","Bad")
Regulon <- c("Bach1","Ppargc1a","Mxi1","Jun")
#Genes <- unique(c(celldeath,DEGs_Neuron))
Genes <- c("Pdia3","Nfkbia", "Casp1", "Bad","Ube2d1", "Mif", "Gba","Bid","Ppif","Bbc3")

#regulon_Neuron <- subset(regulon_target, (TF %in% Regulon))
#regulon_Neuron <- subset(regulon_target, (gene %in% Genes))
regulon_Neuron_select <- subset(regulon_target, gene %in% Genes)
regulon_Neuron_select <- subset(regulon_target, (TF %in% Regulon) & (gene %in% Genes))

regulon_Neuron_select_up <- subset(regulon_target, TF %in% DEGs_Neuron & gene %in% Regulon )
regulon_Neuron_select_down <- subset(regulon_target, gene %in% Regulon & TF %in% DEGs_Neuron )


#regulon_Neuron_select <- subset(regulon_target, TF %in% Regulon | TF %in% DEGs_Neuron |
#                                      gene %in% Regulon | gene %in% DEGs_Neuron & (highConfAnnot == TRUE))

write.csv(regulon_Neuron_select,file = "regulon_Neuron_select_key_Celldeath_gene.csv",row.names = F)
