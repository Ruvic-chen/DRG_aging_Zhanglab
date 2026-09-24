# 加载必要的包
library(Seurat)
library(dplyr)

# 加载您的Seurat对象（已包含Celltype_2_3和age注释）
load("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")
DefaultAssay(DRG_aging_all) <- "RNA"
Idents(DRG_aging_all) <- "age"
# 定义要分析的细胞类型（与原始脚本一致）
clusters <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
              "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9",
              "Satellite","Schwann_N","Schwann_M","Fibroblast","VEC","VSMC",
              "Macrophage", "Monocyte", "Neutrophil", "B Cell", "T Cell")

# 设置降采样参数
downsample_frac <- 0.8      # 每次抽取80%的细胞
n_replicates <- 10          # 重复10次
min_pct_threshold <- 0.2      # 与原始分析一致
logfc_threshold <- 0.7
p_val_adj_threshold <- 0.05

# 存储所有重复的DEG结果
all_downsample_results <- list()
# 可选：设置一个基准种子，确保整体可重复，但每次循环生成不同的子种子
base_seed <- 2026
set.seed(base_seed)

# 生成一组固定的随机种子，供每次循环使用（这样每次循环的结果是可重复的，但循环之间不同）
seed_vec <- sample.int(1e6, size = n_replicates, replace = FALSE)
for (rep in 1:n_replicates) {
  cat("Processing replicate", rep, "of", n_replicates, "\n")
  
  # 对每个年龄组分别降采样
  cells_3mo <- WhichCells(DRG_aging_all, expression = age == "3 MO")
  cells_24mo <- WhichCells(DRG_aging_all, expression = age == "24 MO")
  
  set.seed(seed_vec[rep])   # 每次循环使用不同的种子
  
  sampled_cells_3mo <- sample(cells_3mo, size = floor(length(cells_3mo) * downsample_frac), replace = FALSE)
  sampled_cells_24mo <- sample(cells_24mo, size = floor(length(cells_24mo) * downsample_frac), replace = FALSE)
  
  cells_to_keep <- c(sampled_cells_3mo, sampled_cells_24mo)
  subset_obj <- subset(DRG_aging_all, cells = cells_to_keep)
  
  # 临时存储本次重复的DEG
  deg_list_per_rep <- data.frame()
  
  for (cl in clusters) {
    # 提取该细胞类型
    subset_cl <- subset(subset_obj, subset = Celltype_2_3 == cl)
    if (length(table(subset_cl$age)) < 2) next  # 确保两个年龄组都存在
    
    # 差异分析（使用与原始脚本相同的MAST）
    deg <- FindMarkers(subset_cl, 
                       ident.1 = "24 MO", 
                       ident.2 = "3 MO", 
                       test.use = "MAST",
                       min.pct = min_pct_threshold,
                       logfc.threshold = logfc_threshold)
    
    deg$gene <- rownames(deg)
    deg$clusters <- cl
    deg$replicate <- rep
    deg_list_per_rep <- rbind(deg_list_per_rep, deg)
  }
  
  all_downsample_results[[rep]] <- deg_list_per_rep
}

# 合并所有重复的结果
all_results <- bind_rows(all_downsample_results)
saveRDS(all_results, file = "downsample_results.rds")

