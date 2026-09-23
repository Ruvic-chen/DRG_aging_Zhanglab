setwd("/home/rstudio/DRG_aging_2025/analysis/merged/merged_2_24/merge_final/DEGs/DEG_integrated/")
integrated <- readRDS(file = "integrated_DEGs.rds")

DEGs_robustness_neuron <- subset(integrated, clusters %in% c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                                                                    "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9"))

plot_data <- DEGs_robustness_neuron %>%
  # 统一绘图用的 logFC：优先 MAST 的 avg_log2FC，缺失用 edgeR
  mutate(logFC_plot = ifelse(!is.na(avg_log2FC), avg_log2FC, logFC_edgeR))

# 检查行数
nrow(plot_data)

select_genes <- c("Bad","Cyp1b1","Map3k5","Casp1","Nono", "Pdia3","Nfkbia", "Ube2d1",
                  "Cct2","Cct6a","Cct7","Cct8","Tcp1","Ssb",
                  "Kcnk2","Kcnk4","P2rx4","Asic3","Scn8a","Grin1","Grm3","Penk","Kcnq2")

plot_select <- subset(plot_data, gene %in% select_genes)


library(ComplexHeatmap)
library(circlize)

# 筛选用于展示的基因（至少一种方法认为显著）
plot_data <- plot_select %>%
  mutate(
    logFC_plot = ifelse(!is.na(avg_log2FC), avg_log2FC, logFC_edgeR),
    # 用于星号的 p 值：优先 MAST 的 p_val_adj，否则 edgeR 的 p_adj_edgeR
    p_for_star = ifelse(!is.na(p_val_adj), p_val_adj, p_adj_edgeR)
  )

# 显著性星号函数
sig_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# 关键：只有 robust_downsample == TRUE 时才赋予星号，否则为空
plot_data <- plot_data %>%
  mutate(
    star = ifelse(robust_downsample, sig_star(p_for_star), "")
  )

# logFC 矩阵
logfc_mat <- plot_data %>%
  select(clusters, gene, logFC_plot) %>%
  pivot_wider(
    names_from = gene,
    values_from = logFC_plot,
    values_fn = mean
  ) %>%
  column_to_rownames("clusters") %>%
  as.matrix()

logfc_mat <- logfc_mat[rev(c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                         "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9")),select_genes]
# 星号矩阵
star_mat <- plot_data %>%
  select(clusters, gene, star) %>%
  pivot_wider(
    names_from = gene,
    values_from = star,
    values_fn = function(x) paste(x, collapse = "")
  ) %>%
  column_to_rownames("clusters") %>%
  as.matrix()

# 对齐行列
common_clusters <- intersect(rownames(logfc_mat), rownames(star_mat))
common_genes <- intersect(colnames(logfc_mat), colnames(star_mat))
logfc_mat <- logfc_mat[common_clusters, common_genes, drop = FALSE]
star_mat <- star_mat[common_clusters, common_genes, drop = FALSE]

# 颜色映射
col_fun <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick"))

# 绘制热图
ht <- Heatmap(
  logfc_mat,
  name = "log2FC",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,    # 不聚类细胞类型
  cluster_columns = FALSE, # 不聚类基因
  show_row_names = TRUE,   # 显示细胞类型
  show_column_names = TRUE,# 显示基因
  row_names_side = "left",
  column_names_side = "bottom",  # 基因名放底部，避免拥挤
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 6),  # 基因多时字号调小
  width = ncol(logfc_mat) * unit(6, "mm"),   # 列宽（基因数多时自动变宽）
  height = nrow(logfc_mat) * unit(6, "mm"),  # 行高
  cell_fun = function(j, i, x, y, width, height, fill) {
    s <- star_mat[i, j]
    if (!is.na(s) && nzchar(s)) {
      grid.text(s, x = x, y = y, gp = gpar(fontsize = 8, col = "black"))
    }
  }
)

draw(ht, heatmap_legend_side = "right")

pdf("heatmap_select_genes.pdf", width = 12, height = 10)
draw(ht)
dev.off()



##############
Cytokine <- read.table("/home/rstudio/share/gene_anno_cyt_filter.txt", header = T, sep = "\t")
Cytokine_gene <- Cytokine$Symbol

Peptides <-  read.table("/home/rstudio/share/mmus_endogeneous_peptides.txt", header = F, sep = "\t")
Peptides_gene <- Peptides[,1]

GPCR <- read.table("/home/rstudio/share/GPCRs.txt", header = T, sep = "\t")
GPCR_gene <- GPCR$MGI.symbol

Ion_channel <- read.table("/home/rstudio/share/ion_channels.txt", header = T, sep = "\t")
Ion_channel_gene <- Ion_channel$MGI.symbol

robustness_genes <- unique(DEGs_robustness_neuron$gene[DEGs_robustness_neuron$robust_downsample == TRUE])

Cytokine_select <- intersect(Cytokine_gene,robustness_genes)
Peptides_select <- intersect(Peptides_gene,robustness_genes)
GPCR_select <- intersect(GPCR_gene,robustness_genes)
Ion_channel_select <- intersect(Ion_channel_gene,robustness_genes)

plot_data <- DEGs_robustness_neuron %>%
  # 统一绘图用的 logFC：优先 MAST 的 avg_log2FC，缺失用 edgeR
  mutate(logFC_plot = ifelse(!is.na(avg_log2FC), avg_log2FC, logFC_edgeR))

# 检查行数
nrow(plot_data)


plot_select <- subset(plot_data, gene %in% c(Cytokine_select,Peptides_select,GPCR_select,Ion_channel_select))


library(ComplexHeatmap)
library(circlize)

# 筛选用于展示的基因（至少一种方法认为显著）
plot_data <- plot_select %>%
  mutate(
    logFC_plot = ifelse(!is.na(avg_log2FC), avg_log2FC, logFC_edgeR),
    # 用于星号的 p 值：优先 MAST 的 p_val_adj，否则 edgeR 的 p_adj_edgeR
    p_for_star = ifelse(!is.na(p_val_adj), p_val_adj, p_adj_edgeR)
  )

# 显著性星号函数
sig_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# 关键：只有 robust_downsample == TRUE 时才赋予星号，否则为空
plot_data <- plot_data %>%
  mutate(
    star = ifelse(robust_downsample, sig_star(p_for_star), "")
  )

# logFC 矩阵
logfc_mat <- plot_data %>%
  select(clusters, gene, logFC_plot) %>%
  pivot_wider(
    names_from = gene,
    values_from = logFC_plot,
    values_fn = mean
  ) %>%
  column_to_rownames("clusters") %>%
  as.matrix()

logfc_mat <- logfc_mat[rev(c("C1-1","C1-2-1","C1-2-2","C1-2-3","C1-2-4","C2","C3","C4-1","C4-2",
                             "C5-1","C5-2","C7","C8-1","C8-2","C8-3","C9")),c(Cytokine_select,Peptides_select,GPCR_select,Ion_channel_select)]
# 星号矩阵
star_mat <- plot_data %>%
  select(clusters, gene, star) %>%
  pivot_wider(
    names_from = gene,
    values_from = star,
    values_fn = function(x) paste(x, collapse = "")
  ) %>%
  column_to_rownames("clusters") %>%
  as.matrix()

# 对齐行列
common_clusters <- intersect(rownames(logfc_mat), rownames(star_mat))
common_genes <- intersect(colnames(logfc_mat), colnames(star_mat))
logfc_mat <- logfc_mat[common_clusters, common_genes, drop = FALSE]
star_mat <- star_mat[common_clusters, common_genes, drop = FALSE]

# 颜色映射
col_fun <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick"))

# 绘制热图
ht <- Heatmap(
  logfc_mat,
  name = "log2FC",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,    # 不聚类细胞类型
  cluster_columns = FALSE, # 不聚类基因
  show_row_names = TRUE,   # 显示细胞类型
  show_column_names = TRUE,# 显示基因
  row_names_side = "left",
  column_names_side = "bottom",  # 基因名放底部，避免拥挤
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 6),  # 基因多时字号调小
  width = ncol(logfc_mat) * unit(4, "mm"),   # 列宽（基因数多时自动变宽）
  height = nrow(logfc_mat) * unit(4, "mm"),  # 行高
  cell_fun = function(j, i, x, y, width, height, fill) {
    s <- star_mat[i, j]
    if (!is.na(s) && nzchar(s)) {
      grid.text(s, x = x, y = y, gp = gpar(fontsize = 8, col = "black"))
    }
  }
)

draw(ht, heatmap_legend_side = "right")
pdf("heatmap_select_function_genes.pdf", width = 30, height = 10)
draw(ht)
dev.off()
