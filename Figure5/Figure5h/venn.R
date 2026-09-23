workingDir <- "~/DRG_aging_2025/analysis/merged/merged_2_24/Macrophage/Monocle3/"
setwd(workingDir)

merged_ccl8 <- read.csv(file = "Trajectory_genes_subsampling_freq_ccl8.csv")
sign_gene_ccl8 <- subset(merged_ccl8, regulation %in% c("up","down") & frequency >0.4)
sign_gene_ccl8$type <- "AAMac_Ccl8"

merged_Gpnmb <- read.csv(file = "Trajectory_genes_subsampling_freq_gpnmb.csv")
sign_gene_gpnmb <- subset(merged_Gpnmb, regulation %in% c("up","down") & frequency >0.4)
sign_gene_gpnmb$type <- "AAMac_Gpnmb"
all_Track_genes <- rbind(sign_gene_ccl8,sign_gene_gpnmb)
#write.csv(all_Track_genes,file = "all_Track_genes.csv",row.names = F)

sign_gene_ccl8_up <- subset(sign_gene_ccl8, regulation %in% "up")
sign_gene_gpnmb_up <- subset(sign_gene_gpnmb, regulation %in% "up")

library(VennDiagram)
pdf(file = "venn.pdf",width = 5,height = 3)
venn_up <- venn.diagram(
  x = list(Gpnmb = sign_gene_gpnmb_up$gene, Ccl8 = sign_gene_ccl8_up$gene),
  filename = NULL,  # 不直接保存到文件，而是保存为对象
  fill = c("lightblue", "lightpink"),
  alpha = 0.5,
  cat.cex = 1.5,
  cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = c(0.05, 0.05),
  main = "Up-regulated genes",
  main.cex = 1.5
)
grid.newpage()
grid.draw(venn_up)
dev.off()


up_intersect <- intersect(sign_gene_gpnmb_up$gene_short_name, sign_gene_ccl8_up$gene_short_name)

# 获取两个基因列表
gpnmb_genes <- sign_gene_gpnmb_up$gene_short_name
ccl8_genes <- sign_gene_ccl8_up$gene_short_name

# 找出特有和共有的基因
gpnmb_only <- setdiff(gpnmb_genes, ccl8_genes)
ccl8_only <- setdiff(ccl8_genes, gpnmb_genes)
common_genes <- intersect(gpnmb_genes, ccl8_genes)

# 创建三列数据框（用 NA 填充使长度一致）
max_length <- max(length(gpnmb_only), length(ccl8_only), length(common_genes))

# 填充向量至相同长度
gpnmb_only_filled <- c(gpnmb_only, rep(NA, max_length - length(gpnmb_only)))
ccl8_only_filled <- c(ccl8_only, rep(NA, max_length - length(ccl8_only)))
common_filled <- c(common_genes, rep(NA, max_length - length(common_genes)))

# 创建数据框
result_df <- data.frame(
  gpnmb_only = gpnmb_only_filled,
  ccl8_only = ccl8_only_filled,
  common = common_filled
)
head(result_df)

write.csv(result_df, "AAMac_Trackgene_comparison_results_0817.csv", row.names = FALSE, na = "")
