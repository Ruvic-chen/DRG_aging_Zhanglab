workingDir <- "E:\\DRG_aging\\paper_Wang&Chen2025\\Figure4\\Fig4a&b"
setwd(workingDir)
load("DRG_Macrophage_aging_2_24.Rdata")

library(sccomp)
library(tidyr)
DefaultAssay(DRG_Macrophage_aging) <- "RNA"
sccomp_result = 
  DRG_Macrophage_aging |>
  sccomp_estimate( 
    formula_composition = ~ age, 
    #formula_variability = ~age,
    sample = "orig.ident", 
    cell_group = "Celltype_3", 
    cores = 1,
    verbose = FALSE
  ) |> 
  sccomp_test()
p <- plot(sccomp_result)
pdf(file = "credible_intervals_1D_nocount.pdf",width = 4,height = 5)
p$credible_intervals_1D
dev.off()
write.csv(sccomp_result,file = "sccomp_result_Mac.csv")


sample_anno_filter <- DRG_Macrophage_aging@meta.data
# sample_anno_filter$Celltype_3 <- factor(sample_anno_filter$Celltype_3,levels = ident_levels)
# table(sample_anno$Celltype_3)

sample_anno_filter$age <- factor(sample_anno_filter$age,levels = c("3 MO","24 MO"))
sample_anno_filter <- sample_anno_filter[,c(c("orig.ident","age","sex","Celltype_3"))]
head(sample_anno_filter)

library(reshape2)
sample_anno_2 <- dcast(sample_anno_filter,sex+age+orig.ident ~ Celltype_3)
head(sample_anno_2)
sample_anno_2 <- dcast(sample_anno_filter,age ~Celltype_3)
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
                         aes(x = factor(celltype,levels = new_idents), y = value, fill = factor(variable)), 
                         stat = "identity", position = "fill")+
  scale_fill_manual(values=c("#a9c844","#b2672a"))+
  geom_hline(yintercept=0.5,linetype='dotted')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
  theme(axis.text = element_text(color="black"),panel.border = element_blank())
p
ggsave(filename = "barplot_subtype_age_ratio_2.pdf",plot = p,width = 5,height = 3)
