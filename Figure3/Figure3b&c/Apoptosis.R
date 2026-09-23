#BiocManager::install("escape")
options(stringsAsFactors = F)
setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/apoptosis_analysis/")

# 加载包
library(escape)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(forcats)
apoptotic_genes <- read.csv(file = "apoptotic_genes_20260110.csv")
# 定义凋亡基因集（根据您的研究调整）
pro_apoptotic <- subset(apoptotic_genes,apoptosis_role %in% "pro-apoptotic")$mouse_gene
anti_apoptotic <- subset(apoptotic_genes,apoptosis_role %in% "anti-apoptotic")$mouse_gene

# 检查基因在数据中的存在率
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron_2/DRG_Neuron_aging_2_24_final.Rdata")
DRG_Neuron_aging <- subset(DRG_Neuron_aging, 
                           features = rownames(DRG_Neuron_aging)[Matrix::rowSums(GetAssayData(DRG_Neuron_aging)) > 0])
DRG_Neuron_aging$age <- as.character(DRG_Neuron_aging$age)
DRG_Neuron_aging$age[DRG_Neuron_aging$age %in% "3 MO"] <- "3MO"
DRG_Neuron_aging$age[DRG_Neuron_aging$age %in% "24 MO"] <- "24MO"
table(DRG_Neuron_aging$age)
DRG_Neuron_aging$age <- factor(DRG_Neuron_aging$age,levels = c("3MO","24MO"))
DRG_Neuron_aging$Celltype_2_3 <- factor(DRG_Neuron_aging$Celltype_2_3,levels = c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                                 "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))
DRG_Neuron_aging@meta.data <- DRG_Neuron_aging@meta.data[,c(1:6,87,88,91)]
available_pro <- pro_apoptotic[pro_apoptotic %in% rownames(DRG_Neuron_aging)]
available_anti <- anti_apoptotic[anti_apoptotic %in% rownames(DRG_Neuron_aging)]

cat("促凋亡基因检出率:", length(available_pro)/length(pro_apoptotic)*100, "%\n")
cat("抗凋亡基因检出率:", length(available_anti)/length(anti_apoptotic)*100, "%\n")
apoptosis_genesets <- list(
  `Pro_Apoptosis` = available_pro,
  `Anti_Apoptosis` = available_anti
)
enrichment_scores <- enrichIt(
  obj = DRG_Neuron_aging, 
  gene.sets = apoptosis_genesets, 
  method = "ssGSEA", # 指定方法
  groups = 1000,     # 分组计算以节省内存，可调整
  cores = 4          # 设置并行核数以加速计算
)

DRG_Neuron_aging <- AddMetaData(DRG_Neuron_aging, enrichment_scores)

# 计算平衡指数（Pro - Anti）
DRG_Neuron_aging$Apoptosis_Balance_Index <- 
  DRG_Neuron_aging$Pro_Apoptosis - DRG_Neuron_aging$Anti_Apoptosis
#DRG_Neuron_aging$Apoptosis_Index <- DRG_Neuron_aging$Pro_Apoptosis
DRG_Neuron_aging$barcode <- rownames(DRG_Neuron_aging@meta.data)
plot_data <- DRG_Neuron_aging@meta.data %>%
  select(Cell_Barcode = barcode, 
         Cluster = Celltype_2_3, # 请替换为你的细胞簇列名
         Age = age,           # 请替换为你的年龄列名
         ABI = Apoptosis_Balance_Index,
         pro = Pro_Apoptosis,
         anti = Anti_Apoptosis)
write.csv(plot_data,file = "ABI_data.csv")
# 查看数据结构
head(plot_data)
cat("分析将涉及以下细胞簇：\n", unique(plot_data$Cluster))

# 3. 统计检验：簇水平的年龄组间比较 ----------------------------------------------
# 使用 Wilcoxon 秩和检验（非参数，适用于不服从正态分布的数据）
stat_test_ABI_results <- plot_data %>%
  group_by(Cluster) %>%
  wilcox_test(ABI ~ Age) %>%  # 对每个簇，比较不同年龄组的ABI
  #adjust_pvalue(method = "BH") %>% # 使用Benjamini-Hochberg方法校正p值，控制FDR
  add_significance("p") %>% # 根据校正后p值添加显著性标记（*, **, ***）
  add_xy_position(x = "Cluster", dodge = 0.8) # 为可视化添加位置信息

stat_test_pro_results <- plot_data %>%
  group_by(Cluster) %>%
  wilcox_test(pro ~ Age) %>%  # 对每个簇，比较不同年龄组的ABI
  #adjust_pvalue(method = "BH") %>% # 使用Benjamini-Hochberg方法校正p值，控制FDR
  add_significance("p") %>% # 根据校正后p值添加显著性标记（*, **, ***）
  add_xy_position(x = "Cluster", dodge = 0.8) # 为可视化添加位置信息

stat_test_anti_results <- plot_data %>%
  group_by(Cluster) %>%
  wilcox_test(anti ~ Age) %>%  # 对每个簇，比较不同年龄组的ABI
  #adjust_pvalue(method = "BH") %>% # 使用Benjamini-Hochberg方法校正p值，控制FDR
  add_significance("p") %>% # 根据校正后p值添加显著性标记（*, **, ***）
  add_xy_position(x = "Cluster", dodge = 0.8) # 为可视化添加位置信息


############绘图
library(tidyverse)

# 从Seurat对象元数据中提取关键列，请替换{}内的变量名为你的实际列名
radar_data <- DRG_Neuron_aging@meta.data %>%
  select(
    CellType = Celltype_2_3,  # 例如: "cell_type", "Cluster"
    Age = age,          # 例如: "age", "AgeGroup"
    Bal = Apoptosis_Balance_Index,   # escape包计算的促凋亡分数
    Pro = Pro_Apoptosis, # escape包计算的抗凋亡分数
    Anti = Anti_Apoptosis # escape包计算的抗凋亡分数
  ) %>%
  # 按细胞类型和年龄分组，计算平均富集分数
  group_by(CellType, Age) %>%
  summarise(
    Bal_Mean = mean(Bal, na.rm = TRUE),
    Pro_Mean = mean(Pro, na.rm = TRUE),
    Anti_Mean = mean(Anti, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # 对每个指标，缩放到[0,1]范围以便在同一雷达图上比较
  mutate(
    Bal_Scaled = (Bal_Mean - min(Bal_Mean)) / (max(Bal_Mean) - min(Bal_Mean)),
    Pro_Scaled = (Pro_Mean - min(Pro_Mean)) / (max(Pro_Mean) - min(Pro_Mean)),
    Anti_Scaled = (Anti_Mean - min(Anti_Mean)) / (max(Anti_Mean) - min(Anti_Mean))
  )

# 查看处理后的数据
head(radar_data)

############################
adult_bal_Mean <- subset(radar_data, Age %in% "3MO")
age_bal_Mean <- subset(radar_data, Age %in% "24MO")

OR <- read.csv(file = "/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge/OR.csv")
OR_Neuron <- subset(OR,X %in% c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))

OR_Neuron$Bal_diff <- adult_bal_Mean$Bal_Mean - age_bal_Mean$Bal_Mean
model <- lm(log2OR ~ Bal_diff, data = OR_Neuron)
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
# 查看回归结果
print(reg_label)
#[1] "log10(Gene Count+1) = -724.14 + 1329.708 × AUC\nR² = 0.296, p = 0.029"
print(model_summary)

saveRDS(reg_label,file = "reg_lable.rds")
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron/cols_list.Rdata")
p <- ggplot(OR_Neuron, aes(x = log2OR, y = Bal_diff)) +
  geom_point(aes(color = X), size = 3) +  # 按细胞类型着色
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

#####################
# 1. 将数据转换为长格式
radar_long <- radar_data %>%
  pivot_longer(
    cols = ends_with("_Scaled"),
    names_to = "Signature_Type",
    values_to = "Score"
  ) %>%
  mutate(
    Signature_Type = case_when(
      Signature_Type == "Bal_Scaled" ~ "Bal-apoptotic",
      Signature_Type == "Pro_Scaled" ~ "Pro-apoptotic",
      Signature_Type == "Anti_Scaled" ~ "Anti-apoptotic",
      TRUE ~ Signature_Type
    )
  )
head(radar_long)
# 2. 创建雷达图坐标
# 为每个细胞类型分配一个角度
celltypes <- unique(radar_long$CellType)
n_types <- length(celltypes)
angle_step <- 2 * pi / n_types

coord_df <- radar_long %>%
  mutate(
    angle = (as.numeric(factor(CellType)) - 1) * angle_step, # 计算角度
    x = Score * cos(angle), # 计算x坐标
    y = Score * sin(angle)  # 计算y坐标
  )

# 3. 绘制分面雷达图
library(ggplot2)

p_facet <- ggplot() +
  # 添加网格线
  geom_path(
    data = data.frame(
      x = rep(seq(0, 1, length.out = 5) * cos(seq(0, 2*pi, length.out=100)), each=100),
      y = rep(seq(0, 1, length.out = 5) * sin(seq(0, 2*pi, length.out=100)), each=100),
      group = rep(1:5, each=100)
    ),
    aes(x, y, group = group), color = "gray85", linewidth = 0.3
  ) +
  # 添加轴线
  geom_segment(
    data = data.frame(
      x1 = 0, y1 = 0,
      x2 = cos(seq(0, 2*pi * (n_types-1)/n_types, length.out = n_types)),
      y2 = sin(seq(0, 2*pi * (n_types-1)/n_types, length.out = n_types))
    ),
    aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray70", linewidth = 0.3
  ) +
  # 添加雷达多边形和数据点 (按年龄组)
  geom_polygon(
    data = coord_df,
    aes(x = x, y = y, group = Age, color = Age, fill = Age),
    alpha = 0.15, linewidth = 0.8
  ) +
  geom_point(
    data = coord_df,
    aes(x = x, y = y, color = Age),
    size = 2
  ) +
  # 分面：一行两列，分别显示促凋亡和抗凋亡
  facet_wrap(~ Signature_Type, nrow = 1) +
  # 添加细胞类型标签 (放在最外圈)
  geom_text(
    data = data.frame(
      label = celltypes,
      x = 1.15 * cos(seq(0, 2*pi * (n_types-1)/n_types, length.out = n_types)),
      y = 1.15 * sin(seq(0, 2*pi * (n_types-1)/n_types, length.out = n_types)),
      Signature_Type = "Pro-apoptotic" # 统一标签位置
    ),
    aes(x = x, y = y, label = label), size = 3.5
  ) +
  # 设置颜色和主题
  scale_color_manual(values = c("3MO" = "#2E8B57", "24MO" = "#8B4513")) +
  scale_fill_manual(values = c("3MO" = "#2E8B57", "24MO" = "#8B4513")) +
  coord_equal() +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold", margin = margin(5,0,5,0)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(
    title = "Pro-apoptotic vs Anti-apoptotic Signature Scores by Neuron Type",
    color = "Age Group", fill = "Age Group"
  )

print(p_facet)
# 保存图形
ggsave("Apoptosis_Signatures_Radar_Faceted.pdf", p_facet, width = 12, height = 6, dpi = 300)


heatmap_matrix <- radar_data %>%
  # 创建组合列名，例如 "3MO_Pro", "24MO_Anti"
  pivot_longer(cols = c(Pro_Mean, Anti_Mean), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(Column_Name = paste(Age, gsub("_Mean", "", Metric), sep = "_")) %>%
  select(CellType, Column_Name, Value) %>%
  # 转换为宽格式矩阵
  pivot_wider(names_from = Column_Name, values_from = Value) %>%
  # 将细胞类型设为行名
  column_to_rownames("CellType") %>%
  # 转换为矩阵
  as.matrix()

# 检查矩阵
print("热图数据矩阵 (前几行):")
print(head(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[,c(1,3,2,4)]
# 3. 准备行和列的分组注释信息
# 列注释：区分信号类型 (Pro/Anti) 和年龄组
column_annotation <- data.frame(
  row.names = colnames(heatmap_matrix),
  Age_Group = gsub("_.*", "", colnames(heatmap_matrix)), # 提取年龄部分
  Signature = gsub(".*_", "", colnames(heatmap_matrix))  # 提取信号部分
)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

age_color <- c("3MO" = "#2E8B57", "24MO" = "#8B4513") # 绿色和棕色
signature_color <- c("Pro" = "#d73027", "Anti" = "#4575b4") # 红色和蓝色

annotation_color_list <- list(
  Age_Group = age_color,
  Signature = signature_color
  # 如果有行注释，也在这里添加
  # Brain_Region = region_color
)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
heatmap_plot <- pheatmap(
  mat = heatmap_matrix,
  
  # 颜色与缩放
  color = heatmap_color,
  scale = "row",  # 按行缩放，突出不同细胞类型间的模式差异
  cluster_rows = TRUE,    # 对行（细胞类型）聚类
  cluster_cols = FALSE,   # 不对列聚类，保持Pro/Anti, 3MO/24MO的固定顺序
  
  # 注释信息
  annotation_col = column_annotation,
  # annotation_row = row_annotation, # 如果有行注释则取消注释
  annotation_colors = annotation_color_list,
  
  # 字体与标签
  fontsize_row = 10,
  fontsize_col = 10,
  #angle_col = 45, # 使用数值45而不是字符串"45"
  
  # 图形标题与边距
  main = "Pro- and Anti-apoptotic Signature Scores\nby Neuron Type and Age (Row-scaled)",
  silent = FALSE # 设置为TRUE则不打印图形，直接保存时有用
)
pdf(file = "pheatmap.pdf",width = 5,height = 6)
heatmap_plot
dev.off()
#ggsave(filename = "pheatmap.pdf",plot = heatmap_plot,width = 5,height = 6)
