# # Step 1
# devtools::install_github("MangiolaLaboratory/sccomp")
# 
# # Step 2
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
# 
# # Step 3
# cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
# cmdstanr::install_cmdstan()

setwd(dir = "E:/DRG_aging/sccomp")
library(Seurat)
library(dplyr)
library(sccomp)
library(ggplot2)
library(tidyr)
#head(DRG_aging_all[[]])
load("E:/DRG_aging/paper_Wang&Chen2025/Rdata/Fig3/DRG_aging_all_2_24_final.Rdata")

# DRG_aging_all$group <- "Adult"
# DRG_aging_all$group[DRG_aging_all$age %in% "24 MO"] <- "Aged"
# meta <- DRG_aging_all@meta.data[,c("orig.ident","group","Celltype_2_3")]
# head(meta)
# count_matrix <- GetAssayData(DRG_aging_all, assay = "RNA", layer = "counts")

# 1. 按样本 + 细胞类型分组，聚合 counts
agg_data <- AggregateExpression(
  DRG_aging_all,
  group.by = c( "Celltype_2_3","orig.ident"),   # 分组列
  assays = "RNA",
  slot = "counts",                              # 原始计数
  return.seurat = FALSE                         # 返回矩阵列表
)

# agg_data 是一个列表，取 "RNA" 矩阵
count_matrix_agg <- agg_data$RNA   # 行为基因，列为 "orig.ident_Celltype_2_3" 组合

# 2. 计算每个组合的总计数（即所有基因之和）
total_counts <- colSums(count_matrix_agg)

# 3. 转换为数据框，并拆分为样本和细胞类型
counts_obj <- data.frame(
  sample =  gsub("^.*_", "", names(total_counts)),        # 提取 orig.ident
  cell_group = gsub("_.*$", "", names(total_counts)),      # 提取 Celltype_2_3
  count = as.numeric(total_counts),
  stringsAsFactors = FALSE
)

counts_obj <- counts_obj %>%
  mutate(age_group = case_when(
    grepl("2m|3m", sample) ~ "Adult",
    grepl("24m", sample)   ~ "Aged",
    TRUE                   ~ NA_character_
  ))
counts_obj$count <- as.integer(counts_obj$count)

# 在当前工作目录下创建 .sccomp_models 文件夹
cache_dir <- file.path(getwd(), ".sccomp_models")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# 将缓存目录指向这个文件夹（必须在任何 sccomp 调用之前执行）
utils::assignInNamespace("sccomp_stan_models_cache_dir", cache_dir, ns = "sccomp")
get("sccomp_stan_models_cache_dir", envir = asNamespace("sccomp"))

sccomp_result <- counts_obj %>%
  sccomp_estimate(
    formula_composition = ~ age_group,
    #formula_variability = ~age_group,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    mcmc_seed = 12345,
    cores = 4,
    verbose = TRUE
  ) %>%
  sccomp_test()

#sccomp_test(test_composition_above_logit_fold_change = 0.2)

# 箱线图：显示观测比例与后验预测，显著差异标色
sccomp_boxplot(sccomp_result, factor = "age_group")

# 1D 区间图：每个细胞类型的效应大小和置信区间
sccomp_plot_intervals_1D(sccomp_result)

# 综合绘图（包含多个面板）
p <- plot(sccomp_result)
pdf(file = "credible_intervals_1D.pdf",width = 2,height = 5)
p$credible_intervals_1D
dev.off()
p$credible_intervals_2D
write.table(sccomp_result,file = "sccomp_result.csv")
saveRDS(sccomp_result,file = "sccomp_result.rds")

sccomp_result <- readRDS(file = "sccomp_result.rds")


DefaultAssay(DRG_aging_all) <- "RNA"
sccomp_result = 
  DRG_aging_all |>
  sccomp_estimate( 
    formula_composition = ~ age, 
    formula_variability = ~age,
    sample = "orig.ident", 
    cell_group = "Celltype_2_3", 
    cores = 1,
    test_composition_above_logit_fold_change = 0.1,
    verbose = FALSE
  ) |> 
  sccomp_test()
p <- plot(sccomp_result)
pdf(file = "credible_intervals_1D_nocount.pdf",width = 4,height = 5)
p$credible_intervals_1D
dev.off()
p$credible_intervals_2D
write.csv(sccomp_result,file = "sccomp_result_nocount.csv")
sccomp_result %>%
  sccomp_proportional_fold_change(
    formula_composition = ~ age_group,
    from = "Adult",    # 基线组
    to   = "Aged"      # 对比组
  ) %>%
  select(cell_group, statement)
