setwd(dir = "E:\\DRG_aging\\paper_Wang&Chen2026\\Figure3\\Fig3i")
library(ggplot2)
library(dplyr) # 用于数据操作，如排序
library(forcats)
# 2. 读取你的GO分析CSV文件
# 请将 ‘your_go_file.csv‘ 替换为你的实际文件名和路径
go_data <- read.csv("GO_select.csv", stringsAsFactors = FALSE)
go_data <- go_data[,c(4,17,18)]
# 3. 确保group列为因子类型，并可以指定顺序（例如，先up后down）
go_data <- go_data %>%
  mutate(
    # 将 Log.q.value. 转换为正值，通常直接取绝对值
    Abs.Log.q = abs(Log.q.value.)
  ) %>%
  # 按绝对值从大到小排序，这样最显著的项排在前面
  arrange(desc(Abs.Log.q))


go_data$GeneList <- factor(go_data$GeneList, levels = c("Up_genes", "down_genes"))

p1 <- ggplot(go_data, aes(x = fct_reorder(Description, abs(Log.q.value.)), 
                          y = abs(Log.q.value.), 
                          fill = GeneList)) +
  geom_col() +
  scale_fill_manual(values = c("Up_genes" = "#E41A1C", "down_genes" = "#377EB8")) +
  labs(
    title = "GO富集分析 - 上下调基因分开显示",
    x = "GO Term",
    y = "富集显著性 (-log10 Q value)",
    fill = "调控方向"
  ) +
  facet_wrap(~ GeneList, scales = "free_y", ncol = 1) + # 垂直排列两个面板
  coord_flip() + # 翻转坐标轴，便于阅读GO描述
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.6),
        axis.text.x = element_text(angle = 0, hjust = 1))

print(p1)
ggsave(filename = "Go.pdf",plot = p1,width = 8,height = 4)
