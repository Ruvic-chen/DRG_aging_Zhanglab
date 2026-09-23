# 加载包
library(ggplot2)
library(gganimate)
library(gifski)
library(scales)
library(dplyr)
library(ggthemr)
# ---------- 1. 准备数据 ----------
# 假设 daily_df 已存在，包含 date 和 count
setwd("E:\\web")
# 读取文件
lines <- readLines("recode_202212-202602.txt")
file_names <- list.files(path = "record_search/", full.names = FALSE)
times <- c(lines,file_names)
# 提取日期时间
datetime_str <- substr(times, 1, 10)
datetime <- as.POSIXct(datetime_str, format = "%Y-%m-%d", tz = "UTC")
date <- as.Date(datetime)

# 按日统计
daily_count <- table(date)
daily_df <- data.frame(date = as.Date(names(daily_count)), count = as.numeric(daily_count))
daily_df$date <- as.Date(daily_df$date)
daily_df <- daily_df[order(daily_df$date), ]

# 1. 准备数据
daily_df <- daily_df %>%
  arrange(date) %>%
  mutate(cumulative = cumsum(count))

max_cum <- max(daily_df$cumulative)

# 2. 计算每个整万阈值首次达到的日期
thresholds <- seq(10000, max_cum, by = 10000)
milestone <- data.frame(threshold = thresholds, 
                        date = as.Date(NA), 
                        cumulative = thresholds)

for (i in seq_along(thresholds)) {
  idx <- which(daily_df$cumulative >= thresholds[i])[1]
  if (!is.na(idx)) milestone$date[i] <- daily_df$date[idx]
}
milestone <- milestone[!is.na(milestone$date), ]

# 3. 创建动态图（关键：所有 x 保持 Date 类型）
#ggthemr("dust", layout = "clean", spacing = 2)

theme_custom_bg <- function() {
  theme_minimal() %+replace%
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.3),  # 网格线可选白色半透
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      aspect.ratio = 9/16
    )
}

img_path <- "image.png"  # 支持 jpg, png
if (grepl("\\.png$", img_path)) {
  img <- png::readPNG(img_path)
} else {
  img <- jpeg::readJPEG(img_path)
}

# 转换成 rasterGrob 对象
img_grob <- grid::rasterGrob(img, interpolate = F)

p <- ggplot(daily_df, aes(x = date, y = cumulative)) +
  # 添加背景图片（占据整个绘图区域）
  annotation_custom(img_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(color = "forestgreen", size = 4) +
  geom_point(color = "#e74c3c", size = 8, alpha = 0.9) +
  geom_text(aes(x = min(date), y = max(cumulative) * 0.95,
                label = paste0("  累计检索量: ", cumulative)),
            hjust = 0, vjust = 1, size = 16, color = "black", 
            check_overlap = TRUE)+
  # 最终累计标签
  geom_label(data = tail(daily_df, 1),
             aes(label = paste0("总计: ", cumulative)),
             hjust = -0.1, vjust = 0.5,
             #fill = "white", alpha = 0.8, label.size = 0.3, size = 4,
             color = "#2c3e50") +
  scale_x_date(date_breaks = "2 month", date_labels = "%Y-%m") +
  scale_y_continuous(
    breaks = seq(0, max_cum + 10000, by = 10000),
    labels = label_number(scale = 1e-3, suffix = "k"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_custom_bg() +   # 应用透明背景主题
  theme(axis.title = element_text(face = "bold", size = 24),   # 轴标题加粗
        axis.text = element_text(face = "bold", size = 20),    # 刻度标签加粗
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # x轴标签单独指定（覆盖）
        axis.text.y = element_text(face = "bold"),             # y轴标签加粗
        plot.margin = margin(5, 5, 20, 5),
        # 可选：调整其他元素
        plot.title = element_text(face = "bold", size = 32),)+
  
  labs(
    title = "📈 累计 检索量 动态增长",
    x = "",y = "累计 检索量 (千)"
  ) +
  transition_reveal(date) +
  ease_aes('cubic-in-out')
p
# 4. 生成动画（1000+帧，fps 设为 10，播放约 100 秒）
frame_indices <- seq(1, nrow(daily_df), by = 5)
daily_df_sub <- daily_df[frame_indices, ]
animate(p, 
        nframes = nrow(daily_df_sub), 
        fps = 50, 
        width = 1600, height = 900,
        renderer = av_renderer())  # 需要安装 av 包：install.packages("av")
anim_save("cumulative_animation.mp4")
