setwd("E:\\DRG_aging\\paper_Wang&Chen2025\\Figure3\\Fig3bc")
library(ggplot2)
library(magrittr)
library(dplyr)

DEGs_filter_2_24_MAST<- readRDS("E:\\DRG_aging\\paper_Wang&Chen2025\\Rdata/Fig3/DEGs_24MO_MAST.rds")
df_24MO <- DEGs_filter_2_24_MAST[which(abs(DEGs_filter_2_24_MAST$avg_log2FC) >= 0.7 & DEGs_filter_2_24_MAST$p_val_adj <= 0.05),] 
gene <- unique(df_24MO$gene)
DEGs_Neuron_filter_2_24 <- subset(df_24MO, clusters %in% c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                             "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))
gene_neuron <- unique(DEGs_Neuron_filter_2_24$gene)
write.csv(DEGs_Neuron_filter_2_24,file = "DEGs_neuron.csv")
df <- DEGs_Neuron_filter_2_24[which(abs(DEGs_Neuron_filter_2_24$avg_log2FC) >= 0.7),]
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
  scale_y_continuous(limits =c(0, 2500) ,expand = c(0,0))+
  coord_flip()
p
ggsave(filename = "DEG_Neuron_distribution_MAST_24MO.pdf",width = 4,height = 8)


multi_cell_genes <- df %>%
  # 确保每个基因-细胞类群组合只计数一次
  distinct(gene, clusters, .keep_all = TRUE) %>%
  
  # 计算每个基因出现的细胞类群数量
  group_by(gene) %>%
  mutate(n_celltypes = n_distinct(clusters)) %>%
  ungroup() %>%
  
  filter(n_celltypes >= 1) %>%
  
  # 按基因分组，汇总细胞类群列表
  group_by(gene, n_celltypes) %>%
  summarise(
    cell_types = paste(sort(unique(clusters)), collapse = ", "),
    .groups = "drop"
  ) %>%
  
  # 按细胞类群数量降序排列
  arrange(desc(n_celltypes))


# 3. 输出详细结果到CSV文件
write.csv(multi_cell_genes, "multi_neuron_genes.csv", row.names = FALSE)
cat("结果已保存到 'multi_cell_genes.csv'\n")
