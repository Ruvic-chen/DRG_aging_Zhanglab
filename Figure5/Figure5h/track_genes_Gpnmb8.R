workingDir <- "~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/Monocle3/"
setwd(workingDir)

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(Matrix)
rm(list=ls())

cds_subset <- readRDS(file = "cds_gpnmb.rds")
Track_genes <- read.csv(file = "Trajectory_genes_gpnmb.csv")

# 2. 设置参数
frac <- 0.8                # 每次抽取80%的细胞
n_replicates <- 10         # 重复10次
q_thresh <- 1e-3           # 显著性阈值
total_cells <- ncol(cds_subset)
sample_size <- floor(total_cells * frac)

# 3. 初始化：创建一个空列表，用来存放每次跑出来的显著基因
all_sig_genes <- list()

# 4. 开始循环
for (i in 1:n_replicates) {
  set.seed(123 + i)  # 每次用不同种子，保证结果可重复
  
  # 随机抽取80%的细胞（不放回）
  sampled_cells <- sample(colnames(cds_subset), size = sample_size, replace = FALSE)
  cds_sample <- cds_subset[, sampled_cells]
  
  # 运行拟时序差异分析
  res <- graph_test(cds_sample, neighbor_graph = "principal_graph", cores = 4)
  
  # 筛选显著基因，只保留基因名
  sig_genes <- res %>%
    filter(q_value < q_thresh) %>%
    pull(gene_short_name) %>%
    unique()
  
  # 存到列表里
  all_sig_genes[[i]] <- sig_genes
  
  # 打印进度
  cat("第", i, "次完成，找到显著基因", length(sig_genes), "个\n")
}

# 5. 统计每个基因出现的频率
# 先找出所有被检测到过的基因
all_genes <- unique(unlist(all_sig_genes))

# 计算每个基因在10次里出现了几次
detection_count <- sapply(all_genes, function(g) {
  sum(sapply(all_sig_genes, function(x) g %in% x))
})

# 计算频率（0~1之间）
detection_freq <- detection_count / n_replicates

# 6. 整理成表格，按频率从高到低排序
freq_df <- data.frame(
  gene = all_genes,
  count = detection_count,
  frequency = detection_freq
) %>% arrange(desc(frequency), desc(count))

# 7. 查看结果
print(head(freq_df, 20))  # 打印前20个最稳健的基因

merged <- merge(Track_genes, freq_df, 
                by.x = "gene_short_name", 
                by.y = "gene", 
                all.x = TRUE)

# 把 NA 填充为 0（表示这些基因在10次降采样中从未被检出）
merged$frequency[is.na(merged$frequency)] <- 0
merged$count[is.na(merged$count)] <- 0

load("~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/DRG_Macrophage_aging_2_24.Rdata")
res_gpnmb <- FindMarkers(DRG_Macrophage_aging,features = merged$gene_short_name,ident.1 = "AAM_Trem2",ident.2 = "Mac_Ccr2",test.use = "MAST")
res_gpnmb_1 <- subset(res_gpnmb, pct.1 > 0.2)
res_gpnmb_up <- rownames(res_gpnmb_1[which(res_gpnmb_1$avg_log2FC >= 0),])
merged$regulation <- "NA"
merged$regulation[merged$gene_short_name %in% res_gpnmb_up] <- "up"

res_gpnmb_1 <- subset(res_gpnmb, pct.2 > 0.2)
res_gpnmb_down <- rownames(res_gpnmb_1[which(res_gpnmb_1$avg_log2FC <= 0),])
merged$regulation[merged$gene_short_name %in% res_gpnmb_down] <- "down"


# 8. 保存结果
write.csv(merged, file = "Trajectory_genes_subsampling_freq_gpnmb.csv", row.names = FALSE)

sign_gene <- subset(merged, regulation %in% c("up","down") & frequency >0.4)
table(sign_gene$regulation)
