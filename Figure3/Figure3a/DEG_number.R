setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/")
# 加载必要的包
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
df_24MO <- readRDS(file = "DEGs_24MO_filter_2_24_MAST.rds")
df_24MO_filter <- df_24MO[which(abs(df_24MO$avg_log2FC) >= 0.7 & df_24MO$p_val_adj <= 0.05),]

#df_24MO <- memento_all_celltypes_24MO_results[which(abs(memento_all_celltypes_24MO_results$de_coef) >= 0.5 & memento_all_celltypes_24MO_results$de_pval_adj <= 0.01),] 
#df_32MO <- memento_all_celltypes_32MO_results[which(abs(memento_all_celltypes_32MO_results$de_coef) >= 0.5 & memento_all_celltypes_32MO_results$de_pval_adj <= 0.01),]

df_24MO_filter$regulation <- factor(df_24MO_filter$regulation, levels = c("up","down"))
df_24MO_filter$clusters <- factor(df_24MO_filter$clusters, levels = c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                      "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9",
                                                                      "Satellite","Schwann_N","Schwann_M","Fibroblast","VEC","VSMC",
                                                                      "Macrophage", "Monocyte", "Neutrophil","B Cell",  "T Cell"))

ggplot(df_24MO_filter, aes(x = clusters, fill = regulation)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"))+
  labs(
    title = "Union of DEGs from 24MO by Cluster",
    x = "Clusters",
    y = "Gene Counts",
    fill = "Regulation"
  ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.6),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits =c(0, 1200) ,expand = c(0,0))+
  coord_flip()
ggsave(filename = "DEGs_number_cluster.pdf",width = 5,height = 8)
###

# 获取神经元基因列表
neuron <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
            "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9")
df_Neuron_24MO <- subset(df_24MO_filter, clusters %in% neuron)
up_genes <- unique(df_Neuron_24MO$gene[df_Neuron_24MO$regulation %in% "up"])
down_genes <- unique(df_Neuron_24MO$gene[df_Neuron_24MO$regulation %in% "down"])


# 创建三列数据框（用 NA 填充使长度一致）
max_length <- max(length(up_genes), length(down_genes))

# 填充向量至相同长度
up_genes_filled <- c(up_genes, rep(NA, max_length - length(up_genes)))
down_genes_filled <- c(down_genes, rep(NA, max_length - length(down_genes)))

# 创建数据框
result_df <- data.frame(
  up_genes = up_genes_filled,
  down_genes = down_genes_filled
)

# 查看前几行
head(result_df)

# 保存结果
write.csv(result_df, "DEGs_regulation_results.csv", row.names = FALSE, na = "")

