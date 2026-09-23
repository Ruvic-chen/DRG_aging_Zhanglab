setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEG_integrated/")

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(dplyr)

integrated <- readRDS(file = "integrated_DEGs.rds")
clusters <- c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
             "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9")
DEGs_robustness_neuron <- subset(integrated, clusters %in% clusters)

select_genes <- c("Bad","Cyp1b1","Map3k5","Casp1","Nono", "Pdia3","Nfkbia", "Ube2d1",
                  "Cct2","Cct6a","Cct7","Cct8","Tcp1","Ssb",
                  "Kcnk2","Kcnk4","P2rx4","Asic3","Scn8a","Grin1","Grm3","Penk","Kcnq2")

df_filter <- subset(DEGs_robustness_neuron, abs(avg_log2FC) > 0.3)

df_filter <- df_filter %>%
  mutate(
    is_robust = robust_downsample == TRUE,
    is_selected = gene %in% select_genes & robust_downsample == TRUE,
    # 颜色分组：依据是否稳健
    color_group = ifelse(is_robust, "robust", "nonrobust"),
    # 大小分组：依据是否在 select_genes 中
    size_group = ifelse(is_selected, "selected", "nonselected"),
    color_group = factor(color_group, levels = c("nonrobust", "robust")),
    size_group = factor(size_group, levels = c("nonselected", "selected"))
  )


dfbar <- data.frame(x = clusters, y = rep(6, length(clusters)))
dfbar$x <- factor(dfbar$x, levels = clusters)
dfbar1 <- data.frame(x = clusters, y = rep(-6, length(clusters)))
dfbar1$x <- factor(dfbar1$x, levels = clusters)

# 颜色块数据（用于 cluster 标签）
load("~/DRG_aging_2025/analysis/merged/merged_2_24/Neuron/cols_list.Rdata")

dfcol <- data.frame(x = clusters, 
                    y = 0, 
                    label = clusters,
                    fill_color = cols_list)

# ==================== 5. 绘制基础图层 ====================
p1 <- ggplot() +
  geom_col(data = dfbar, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6)

# ==================== 6. 添加散点图层（颜色：稳健，大小：selected） ====================
p2 <- p1 +
  geom_jitter(
    data = subset(df_filter, color_group == "nonrobust"),
    aes(x = clusters, y = avg_log2FC, 
        color = color_group, 
        size = size_group),
    shape = 16,
    width = 0.4
  ) +
  # 再绘制稳健基因（红色，大小根据是否selected）
  geom_jitter(
    data = subset(df_filter, color_group == "robust"),
    aes(x = clusters, y = avg_log2FC, 
        color = color_group, 
        size = size_group),
    shape = 16,
    width = 0.4
  ) +
  scale_color_manual(
    name = NULL,
    values = c("nonrobust" = "grey60", 
               "robust" = "#E41A1C")
  ) +
  scale_size_manual(
    values = c("nonselected" = 0.5, 
               "selected" = 2),
    guide = "none"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2)   # 图例中点大小统一，避免过小
    )
  )
p2
# ==================== 7. 添加 cluster 标签色块 ====================
p3 <- p2 + 
  geom_tile(data = dfcol,
            aes(x = x, y = y, fill = fill_color),
            height = 0.4,
            color = "black",
            alpha = 0.6,
            show.legend = FALSE) +
  scale_fill_identity()

# ==================== 8. 添加 x 轴标签文本（底部） ====================
p4 <- p3 +
  geom_text(data = dfcol, aes(x = x, y = y, label = label), 
            size = 3, color = "black")

# ==================== 9. 为 selected 基因添加斜体文字标签（置顶，白底） ====================
p5 <- p4 +
  geom_text_repel(
    data = subset(df_filter, is_selected),
    aes(x = clusters, y = avg_log2FC, label = gene),
    color = "black",
    size = 3,
    fontface = "italic",          # 斜体
    segment.color = "black",
    segment.size = 0.5,
    show.legend = FALSE,
    fill = "white",
    label.size = 0.25,
    label.padding = unit(0.2, "lines"),
    box.padding = unit(0.3, "lines"),
    force = 1.5,
    max.overlaps = Inf
  )

# ==================== 10. 设置主题 ====================
p6 <- p5 +
  labs(x = "Cluster", y = "average logFC") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),          # 隐藏默认 x 轴文字，因为已用 geom_text 添加
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 15)
  )

# ==================== 11. 显示并保存 ====================
print(p6)

ggsave("robust_selected_genes_plot.pdf", plot = p6, width = 8, height = 3)
