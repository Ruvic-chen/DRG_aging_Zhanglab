# 加载必要的包
setwd(dir = "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/Augur/")
library(tidyverse)
library(ggpmisc) # 用于添加回归方程
library(dplyr)


load("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DRG_aging_all_2_24_final.Rdata")
library(Augur)
library(Seurat)

table(DRG_aging_all$age)
DRG_aging_all$age <- factor(DRG_aging_all$age,levels = c("3 MO","24 MO"))
augur2 <- calculate_auc(DRG_aging_all,cell_type_col = "Celltype_2_3",label_col = "age",augur_mode = "velocity",n_threads = 8)
head(augur2$AUC, 5)

save(augur2,file = "Augur_DRG_aging_all.Rdata")
load("~/DRG_aging_2025/analysis/merged/merged_2_24/merge/Augur/Augur_DRG_aging_all.Rdata")

DEGs_filter_2_24_MAST <- readRDS("~/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEGs_24MO_MAST_filter.rds")

#df_24MO <- DEGs_filter_2_24_MAST[which(abs(DEGs_filter_2_24_MAST$avg_log2FC) >= 0.7),] 

regulation_stats <- df_24MO %>%
  group_by(clusters, regulation) %>%
  summarise(gene_count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = regulation, values_from = gene_count, values_fill = 0)
regulation_stats$sum <- rowSums(regulation_stats[,c(2:3)])

colnames(regulation_stats)[1] <- c("cell_type")
rownames(regulation_stats) <- regulation_stats$cell_type
augur_score <- augur2$AUC
augur_score$cell_type <- as.character(augur_score$cell_type)
rownames(augur_score) <- augur_score$cell_type
augur_score <- augur_score[regulation_stats$cell_type,]

regulation_stats$auc <- augur_score$auc

colnames(regulation_stats) <- c("cell_type","down","up","diff_gene_count","augur_score")
Neuron <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4",
  "C2","C3","C4-1","C4-2","C5-1","C5-2","C7",
  "C8-1","C8-2","C8-3","C9")
regulation_stats <- subset(regulation_stats,cell_type %in% Neuron)

# 执行线性回归分析
model <- lm(diff_gene_count ~ augur_score, data = regulation_stats)
model_summary <- summary(model)

# 提取回归统计信息
r_squared <- round(model_summary$r.squared, 3)
p_value <- ifelse(model_summary$coefficients[2, 4] < 0.001, 
                  "< 0.001", 
                  round(model_summary$coefficients[2, 4], 3))
slope <- round(model_summary$coefficients[2, 1], 3)
intercept <- round(model_summary$coefficients[1, 1], 3)

# 创建回归方程标签
reg_label <- paste0("log10(Gene Count+1) = ", intercept, " + ", slope, " × AUC\n",
                    "R² = ", r_squared, ", p = ", p_value)
saveRDS(reg_label,file = "reg_lable.rds")
# 查看回归结果
print(reg_label)
#"log10(Gene Count+1) = -659.291 + 1212.402 × AUC\nR² = 0.295, p = 0.03"

regulation_stats$cell_type <- factor(regulation_stats$cell_type,levels = Neuron)
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron/cols_list.Rdata")
p <- ggplot(regulation_stats, aes(x = augur_score, y = diff_gene_count)) +
  geom_point(aes(color = cell_type), size = 3) +  # 按细胞类型着色
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey") +  # 添加线性回归线和置信区间
  labs(
    x = "AUC from Augur",
    y = "Age-DE gene count",
    color = "Celltype"
  ) +
  scale_color_manual(values = cols_list) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_blank(),  # 移除背景
    panel.grid = element_blank(),  # 移除网格线
    axis.line = element_line(color = "black", size = 0.5),  # 添加黑色坐标轴
    axis.ticks = element_line(color = "black", size = 0.5),  # 添加刻度线
    legend.position = "right",
    text = element_text(size = 12)
  )
p
ggsave("DEG_Augur_24MO.pdf",p,width = 8,height = 5)
