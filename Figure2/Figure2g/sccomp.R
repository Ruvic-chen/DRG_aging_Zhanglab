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

write.table(sccomp_result,file = "sccomp_result.csv")
saveRDS(sccomp_result,file = "sccomp_result.rds")

sccomp_result <- readRDS(file = "sccomp_result.rds")
