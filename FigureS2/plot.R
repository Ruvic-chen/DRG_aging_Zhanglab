# 加载必要的包
library(ggplot2)
setwd("E:\\web")
# 读取文件
lines <- readLines("recode_202212-202602.txt")
file_names <- list.files(path = "record_search/", full.names = FALSE)
times <- c(lines,file_names)
# 提取日期时间
datetime_str <- substr(times, 1, 19)
datetime <- as.POSIXct(datetime_str, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
date <- as.Date(datetime)

# 按日统计
daily_count <- table(date)
daily_df <- data.frame(date = as.Date(names(daily_count)), count = as.numeric(daily_count))
daily_df$date <- as.Date(daily_df$date)
daily_df <- daily_df[order(daily_df$date), ]

# 计算累计计数
daily_df$cumulative <- cumsum(daily_df$count)

library(ggplot2)
library(gganimate)
library(gifski)  # 用于生成 GIF，如需其他渲染器也可

p <- ggplot(daily_df, aes(x = date, y = cumulative)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 2) +
  labs(
    title = "累计 counts 数量随时间变化",
    subtitle = "日期: {frame_along}",
    x = "日期", y = "累计 counts"
  ) +
  theme_minimal() +
  transition_reveal(date)  # 沿着日期逐步显示曲线

# 渲染动画（默认生成 GIF）
animate(p, nframes = nrow(daily_df), fps = 5, width = 800, height = 500)

# 绘制日检索图
p_daily <- ggplot(daily_df, aes(x = date, y = count)) +
  geom_line(color = "steelblue") +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(title = "每日检索次数 (2022-12 ~ 2024-12)",
       x = "日期", y = "检索次数") +
  theme_minimal() +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_daily)

# 按月统计
month <- format(date, "%Y-%m")
monthly_count <- table(month)
monthly_df <- data.frame(month = names(monthly_count), count = as.numeric(monthly_count))
monthly_df$month <- factor(monthly_df$month, levels = unique(monthly_df$month))

p_monthly <- ggplot(monthly_df, aes(x = month, y = count, group = 1)) +
  geom_line(color = "darkred") +
  geom_point() +
  labs(title = "月度检索次数", x = "月份", y = "次数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_monthly)


# 提取日期
datetime_str <- substr(lines, 1, 19)
datetime <- as.POSIXct(datetime_str, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
date <- as.Date(datetime)

# 按日统计
daily_count <- table(date)
daily_df <- data.frame(date = as.Date(names(daily_count)), 
                       count = as.numeric(daily_count))

# 排序并计算累计
daily_df <- daily_df[order(daily_df$date), ]
daily_df$cumulative <- cumsum(daily_df$count)

# 绘图
ggplot(daily_df, aes(x = date, y = cumulative)) +
  geom_line(color = "forestgreen", size = 2) +
  labs(title = "累计检索次数 (2022-12 ~ 2024-12)",
       x = "日期", y = "累计检索次数") +
  theme_minimal() +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
