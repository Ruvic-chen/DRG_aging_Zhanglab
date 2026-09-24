setwd(dir = "~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEG_integrated")

downsample_results <- readRDS(file = "~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/downsample_results.rds")

logfc_threshold <- 0.7
p_val_adj_threshold <- 0.05
n_replicates <- 10 
# 计算每个基因在每个细胞类型中的检出频率
robustness <- downsample_results %>%
  filter(p_val_adj < p_val_adj_threshold & abs(avg_log2FC) > logfc_threshold) %>%
  group_by(clusters, gene) %>%
  summarise(detection_frequency = n() / n_replicates, .groups = "drop")

# 查看最稳健的DEG（频率 >= 0.8）
robust_genes <- robustness %>%
  filter(detection_frequency >= 0.5) %>%
  arrange(clusters, desc(detection_frequency))

# 可选：与原始DEG结果合并，标记稳健性
original_deg <- readRDS("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEGs_24MO_MAST.rds")  # 您之前保存的原始结果

original_deg_with_robustness <- original_deg %>%
  left_join(robustness, by = c("clusters" = "clusters", "gene" = "gene"))

edgeR_res <- readRDS("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/edgeR/DEGs_24MO_filter_2_24_edgeR.rds")
edgeR_sig <- edgeR_res %>%
  select(cell_type, gene, avg_logFC, p_val_adj, regulation) %>%
  rename(clusters = cell_type, logFC_edgeR = avg_logFC, p_adj_edgeR = p_val_adj)

integrated <- original_deg_with_robustness %>%
  full_join(edgeR_sig, by = c("clusters", "gene"))

integrated <- integrated %>%
  mutate(
    detection_frequency = ifelse(is.na(detection_frequency), 0, detection_frequency),
    logFC_edgeR = ifelse(is.na(logFC_edgeR), NA, logFC_edgeR),
    p_adj_edgeR = ifelse(is.na(p_adj_edgeR), NA, p_adj_edgeR),
    # 定义稳健性阈值：检测频率 ≥ 0.7 视为稳健
    robust_downsample = detection_frequency >= 0.5,
    significant_edgeR = !is.na(p_adj_edgeR),
    # 综合分类
    validation_status = case_when(
      significant_edgeR & robust_downsample ~ "High confidence (both methods)",
      significant_edgeR & !robust_downsample ~ "EdgeR only (check robustness)",
      !significant_edgeR & robust_downsample ~ "Downsample only (possible false positive in MAST?)",
      TRUE ~ "Not significant in either"
    )
  )
#write.csv(integrated, file = "integrated_DEGs.csv")
saveRDS(integrated,file = "integrated_DEGs.rds")
integrated <- readRDS(file = "integrated_DEGs.rds")

integrated_filter <- subset(integrated,robust_downsample == TRUE | significant_edgeR == TRUE)
integrated_filter <- integrated_filter %>%
  mutate(
    regulation = case_when(
      avg_log2FC > 0 ~ "up",
      avg_log2FC <0 ~"down"
    )
  )
saveRDS(integrated_filter,file = "robust_ageDEGs.rds")
integrated_filter <- readRDS(file = "robust_ageDEGs.rds")
table(integrated_filter$validation_status)

integrated <- integrated %>%
  arrange(clusters, desc(significant_edgeR), desc(detection_frequency))

plot_data <- integrated %>%
  filter(significant_edgeR) %>%
  mutate(
    detection_frequency = detection_frequency * 100,  # 转为百分比
    cell_type = factor(clusters, levels = sort(unique(clusters)))
  )
plot_data$clusters <- factor(plot_data$clusters, levels = clusters)
plot_data$regulation <- factor(plot_data$regulation, levels = c("up","down"))
if (nrow(plot_data) > 0) {
  p <- ggplot(plot_data, aes(x = detection_frequency, y = logFC_edgeR, color = regulation)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~clusters, scales = "free_y") +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    labs(
      x = "Detection frequency in downsampling MAST (%)",
      y = "log2 Fold Change (edgeR pseudobulk)",
      color = "Regulation (edgeR)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "lightgray"),
      legend.position = "bottom"
    )
  
  ggsave("DEG_consistency_downsample_vs_edgeR.pdf", plot = p, width = 12, height = 8)
  print(p)
} else {
  message("No significant genes in edgeR to plot.")
}


DEGs_robustness_neuron <- subset(integrated_filter, clusters %in% c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                  "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))
DEGs_numbers_neuron <- unique(DEGs_robustness_neuron$gene)

df<- DEGs_robustness_neuron[,c("clusters","gene","regulation")]
#write.csv(df,file = "DEG_for_GO.csv")
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


DEGs_robustness_neuron$regulation <- factor(DEGs_robustness_neuron$regulation, levels = c("up","down"))
DEGs_robustness_neuron$clusters <- factor(DEGs_robustness_neuron$clusters, levels = c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                      "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))

ggplot(DEGs_robustness_neuron, aes(x = clusters, fill = regulation)) +
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
  scale_y_continuous(limits =c(0, 600) ,expand = c(0,0))+
  coord_flip()
ggsave(filename = "DEGs_number_cluster.pdf",width = 5,height = 8)

Up_genes <- unique(subset(df, regulation %in% "up")$gene)
down_genes <- unique(subset(df, regulation %in% "down")$gene)
max_len <- max(length(Up_genes), length(down_genes))

# 构建数据框，较短向量用 NA 补齐
df_GO <- data.frame(
  Up_genes   = c(Up_genes,   rep(NA, max_len - length(Up_genes))),
  down_genes = c(down_genes, rep(NA, max_len - length(down_genes)))
) 
write.csv(df_GO,file = "DEG_for_GO.csv",row.names = F)
