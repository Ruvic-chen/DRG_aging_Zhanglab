setwd(dir = "D:\\paper_Wang&Chen2025\\Figure2")

data_OR <- read.csv(file = "OR_plot.csv",row.names = 1)


library(ggplot2)
library(RColorBrewer)
scale.func <- switch(
  EXPR = 'radius',
  'size' = scale_size,
  'radius' = scale_radius,
  stop("'scale.by' must be either 'size' or 'radius'")
)
data_OR$group <- "OR"
data_OR$cluster <- rownames(data_OR)
data_OR$cluster <- factor(data_OR$cluster, levels = c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                      "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9",
                                                      "Satellite","Schwann_N","Schwann_M","Fibroblast","VEC","VSMC",
                                                      "Macrophage", "Monocyte", "Neutrophil","B Cell",  "T Cell"))
data_OR$log10p <- -log10(data_OR$p_val)
bb <- ggplot(data = data_OR, mapping = aes(x = group, y = cluster)) +
  geom_point(mapping = aes(size = log10p, color = Log2.OR._all)) +
  scale.func(range = c(0, 5), limits = c(NA, NA)) +
  theme_void() + #?հ?????
  #  geom_text(aes(label = paste0("(", avg.exp, ")")), nudge_y = -0.25) +  # 
  coord_fixed(1) + #??????
  #  guides(color = FALSE) + #ȥ??colorͼ??
  theme(axis.title.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "sans",angle = 90, vjust = 0.5,hjust = 1),
        axis.text.y = element_text(family = "sans",color="black",hjust = 1)) + #
  
  scale_color_gradientn(limits=c(-2.5,2.6),breaks = c(-2.5,0,2.5),
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  ) + 
  scale_size_area(max_size = 6, breaks=c(0,1,2),limits = c(0,3)) +
  coord_flip()
bb
ggsave(bb,filename = "OR_dotplot.pdf",width = 5,height = 4)
