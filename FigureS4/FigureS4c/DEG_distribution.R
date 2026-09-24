setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEG_integrated/")
library(ggplot2)
library(magrittr)
library(dplyr)

df_24MO<- readRDS("robust_ageDEGs.rds")

df <- subset(df_24MO, clusters %in% c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                             "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))
#write.csv(df,file = "DEGs_Neuron_24mo_MAST_0.5.csv")

df<- df[,c("clusters","gene","regulation")]

# 1. 计算每个基因在多少种细胞类群中出现
# 确保每个基因在每个细胞类群只计数一次
unique_pairs <- unique(df[c("gene", "clusters")])

# 计算每个基因的细胞类群计数
gene_counts <- table(unique_pairs$gene)

# 转换为数据框
gene_counts_df <- data.frame(
  gene = names(gene_counts),
  n_celltypes = as.numeric(gene_counts)
)

# 2. 统计每个细胞类群数量对应的基因数量
dist_counts <- as.data.frame(table(gene_counts_df$n_celltypes))
colnames(dist_counts) <- c("n_celltypes", "n_genes")

# 确保因子顺序正确
dist_counts$n_celltypes <- factor(dist_counts$n_celltypes, 
                                  levels = sort(as.numeric(levels(dist_counts$n_celltypes))))

# 3. 打印统计结果
cat("基因分布统计:\n")
print(dist_counts)
dist_counts$n_genes_2 <- log10(dist_counts$n_genes)
# 4. 使用ggplot2绘图
p <- ggplot(dist_counts, aes(x = n_celltypes, y = n_genes)) +
    geom_bar(stat = "identity", fill = "black", width = 0.7) +
    geom_text(aes(label = n_genes), vjust = -0.5, size = 4, color = "black") +
    scale_y_continuous(trans = 'log10')+
    labs(x = "Number of cell types/subclasses", 
         y = "Number of genes",
         title = "Distribution of Age-DE Genes Across Cell Types",
         subtitle = "Histogram of the number of cell types an age-DE gene is significant for") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.6),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits =c(0, 2000) ,expand = c(0,0))+
  coord_flip()
p
ggsave(filename = "DEG_Neuron_distribution_MAST_24MO.pdf",width = 4,height = 8)

