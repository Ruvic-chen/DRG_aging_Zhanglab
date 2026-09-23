setwd("/paper????????/Figure 2/Figure 2d/")
library(Seurat)
library(sctransform)
library(ggplot2)
library(RColorBrewer)

load("/paper????????/Figure 2/Figure 2a/DRG_aging_all_2_24_final.Rdata")

sample_anno <- DRG_aging_all@meta.data
head(sample_anno)
sample_anno_filter <- sample_anno
sample_anno_filter$Celltype_2_3 <- factor(sample_anno_filter$Celltype_2_3,levels = ident_levels)
table(sample_anno$Celltype_2_3)

sample_anno_filter$age <- factor(sample_anno_filter$age,levels = c("3 MO","24 MO"))
sample_anno_filter <- sample_anno_filter[,c(c("orig.ident","age","sex","Celltype_2_3"))]
head(sample_anno_filter)

library(reshape2)
sample_anno_2 <- dcast(sample_anno_filter,sex+age+orig.ident ~ Celltype_2_3)
head(sample_anno_2)
sample_anno_2 <- dcast(sample_anno_filter,age ~Celltype_2_3)
head(sample_anno_2)
rownames(sample_anno_2) <- sample_anno_2$age
sample_anno_2 <- sample_anno_2[,-1]
row_num <- rowSums(sample_anno_2)
f1 <- function( x ) x * 10000 / row_num
sample_anno_3 <- as.data.frame(apply(sample_anno_2,2,f1))
sum_num <- colSums(sample_anno_3)
f2 <- function( x ) x / sum_num
sample_anno_percentage <- as.data.frame(apply(sample_anno_3,1,f2))
sample_anno_percentage$celltype <- rownames(sample_anno_percentage)
head(sample_anno_percentage)
sample_anno_percentage_melt <- melt(sample_anno_percentage,id.vars = "celltype")
head(sample_anno_percentage_melt)

p <- ggplot() + geom_bar(data =sample_anno_percentage_melt, 
                         aes(x = factor(celltype,levels = ident_levels), y = value, fill = factor(variable)), 
                         stat = "identity", position = "fill")+
  scale_fill_manual(values=c("#a9c844","#b2672a"))+
  geom_hline(yintercept=0.5,linetype='dotted')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
  theme(axis.text = element_text(color="black"),panel.border = element_blank())
p
ggsave(filename = "barplot_subtype_age_ratio_2.pdf",plot = p,width = 20,height = 3)

