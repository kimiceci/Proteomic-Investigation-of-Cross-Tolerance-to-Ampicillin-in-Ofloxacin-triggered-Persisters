library(DEP2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggvenn)
library(tidyverse)

# 设置工作目录
setwd("D:/Final Report")
filePath <- "proteinGroups_6500_filtered_imputed.csv"
data <- read.csv(filePath, stringsAsFactors = FALSE)

# 假设您的数据框包含蛋白质ID列（通常是"Protein.IDs"或"Majority.protein.IDs" ）
data_unique <- make_unique(
  data,
  names = "name",
  ids = "ID",   # 蛋白质ID列 
  delim = ";"            # 分隔符（用于多ID情况）
)

lfq_col_indices <- grep("^LFQ.intensity.exp", colnames(data_unique))

# OFX: Ofloxacin; AMP: Ampicillin; CTRL: Control
expdesign <- data.frame( 
  label = c("exp7.1", "exp7.2", "exp7.3", "exp1.1", "exp1.2", "exp1.3", "exp6.1", "exp6.2", "exp6.3"),
  condition = c(rep("OFX-AMP",3), rep("OFX-OFX",3),rep("CTRL",3)),
  replicate = c(1, 2, 3,1, 2, 3, 1, 2, 3)
)

data_se <- make_se(data_unique, lfq_col_indices, expdesign)

#使用 limma 进行差异富集测试
data_diff_all_contrasts <- test_diff(data_se, type = "all", fdr.type = "BH")
# Tested contrasts: OFX.OFX_vs_OFX.AMP, OFX.OFX_vs_CTRL, OFX.AMP_vs_CTRL

#标记显著差异表达的蛋白质
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(2))

#得到结果
results <- get_results(dep)

# 手动绘制火山图函数
drawVolcano <- function(df, plotTitle) {
  volcanoPlot <- ggplot(df, aes(x = logFC, y = -log10(p.adj))) +
  geom_point(color = "#BEBEBE", size = 1) +
  geom_point(data = df[df$significant & df$direction == "Upregulated", ], 
             color = "#B2182B", size = 1) +
  geom_point(data = df[df$significant & df$direction == "Downregulated", ], 
             color = "#2166AC", size = 1) +
  scale_x_continuous(
    breaks = c(-10, -5, 0, 5, 10, -log2(2),log2(2)),  # 包含 cutoff 值 
    labels = c(-10, -5, 0, 5, 10, -1, 1)  # 自定义标签 
  ) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6, 8, -log10(0.05)),  # 包含 cutoff 值 
    labels = c(0, 2, 4, 6, 8, expression("-log"[10] ~ "0.05"))  # 自定义标签 
  ) +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05),linetype=2,color = "#BEBEBE")+        # 在图上添加虚线
  geom_vline(xintercept=c(-log2(2),log2(2)),linetype=2,color = "#BEBEBE")+ # 在图上添加虚线
  labs(title = plotTitle, x = expression("log"[2] ~ "Fold Change"), y = expression("-log"[10] ~ "p.adj"))+  # 添加标题 
  theme(
    axis.title  = element_text(size = 14),
    plot.title  = element_text(
      hjust = 0.5,      # 水平居中（0.5表示中心）
      face = "bold",
      size = 16# 字体加粗 
      ),
    panel.border  = element_rect(
      color = "black",    # 边框颜色 
      fill = NA,          # 无填充色
      linewidth = 0.5       # 边框粗细
      )
    ) +
  geom_text_repel(
    data = df[df$significant & df$direction == "Upregulated",],
    aes(label = name),  # 使用基因名作为标签
    box.padding = 0.35,  # 调整标签与点的距离
    point.padding = 0.5,  # 调整标签与点的间距
    size = 3,             # 标签字体大小
    max.overlaps = 15,   # 允许的最大重叠次数，减少未标记的数据点
    segment.size = 0.2,  # 连接线宽度
    segment.color = "#B2182B",  # 连接线颜色
    color = "#B2182B",  # 标签颜色
    nudge_x = 0.1,  # 水平偏移量
    nudge_y = 0.1,  # 垂直偏移量
    direction = "y"  # 标签方向)
  ) +
  geom_text_repel(
    data = df[df$significant & df$direction == "Downregulated",],
    aes(label = name),  # 使用基因名作为标签
    box.padding = 0.35,  # 调整标签与点的距离
    point.padding = 0.5,  # 调整标签与点的间距
    size = 3,             # 标签字体大小
    max.overlaps = 15,   # 允许的最大重叠次数，减少未标记的数据点
    segment.size = 0.2,  # 连接线宽度
    segment.color = "#2166AC",  # 连接线颜色
    color = "#2166AC",  # 标签颜色
    nudge_x = 0.1,  # 水平偏移量
    nudge_y = 0.1,  # 垂直偏移量
    direction = "y"  # 标签方向)
  )
  savePath <- paste0(gsub(": ","_", plotTitle), ".pdf")
  ggsave(savePath, plot = volcanoPlot, device = "pdf", width = 8, height = 6, dpi = 300)
  return(volcanoPlot)
}

# 提取基因显着性列表函数
significantGenes <- function(df, cont) {
  sig_genes <- df %>%
    dplyr::filter(
      !is.na(name)                                     # 确保基因名非空 
    ) %>%
    dplyr::select(
      name, 
      ID, 
      logFC = paste0(cont,"_ratio"),       # 比值列对应log2FC 
      p.val  = paste0(cont,"_p.val"),        # 原始p值 
      p.adj  = paste0(cont, "_p.adj"),         # 校正后p值 
      significant = paste0(cont, "_significant"),  # 自动生成的显著性列
    ) %>%
    dplyr::mutate(
      direction = ifelse(logFC > 0, "Upregulated", ifelse(logFC<0, "Downregulated", "None"))  # 添加上下调信息 
    )
  # 保存为CSV
  savePath <- paste0(gsub("\\.","-",gsub("_"," ",cont)),"_genes.csv")
  write.csv(sig_genes,  savePath,  row.names  = FALSE)
  return(sig_genes)
}

########## OFX.OFX_vs_CTRL ##########
# 提取基因显着性列表（基于你的实际列名）
genes_OFX.OFX_vs_CTRL <- significantGenes(results, "OFX.OFX_vs_CTRL")

# 手动绘制火山图
OFX.OFX_vs_CTRL_volcano <- drawVolcano(genes_OFX.OFX_vs_CTRL, "OFX-OFX vs CTRL: Volcano Plot")

OFX.OFX_vs_CTRL_volcano

########## OFX.AMP_vs_CTRL ##########
# 提取基因显着性列表（基于你的实际列名）
genes_OFX.AMP_vs_CTRL <- significantGenes(results, "OFX.AMP_vs_CTRL")

# 手动绘制火山图
OFX.AMP_vs_CTRL_volcano <- drawVolcano(genes_OFX.AMP_vs_CTRL, "OFX-AMP vs CTRL: Volcano Plot")

OFX.AMP_vs_CTRL_volcano

########## OFX.AMP_vs_OFX.OFX ##########
# 提取基因显着性列表（基于你的实际列名）
genes_OFX.AMP_vs_OFX.OFX <- significantGenes(results, "OFX.AMP_vs_OFX.OFX")

# 手动绘制火山图
OFX.AMP_vs_OFX.OFX_volcano <- drawVolcano(genes_OFX.AMP_vs_OFX.OFX, "OFX-AMP vs OFX-OFX: Volcano Plot")

OFX.AMP_vs_OFX.OFX_volcano

########## Common Genes ##########
OFX.OFX_vs_CTRL_significant <- genes_OFX.OFX_vs_CTRL %>%
  filter(significant == TRUE)
OFX.AMP_vs_CTRL_significant <- genes_OFX.AMP_vs_CTRL %>%
  filter(significant == TRUE)
OFX.OFX_vs_CTRL_up <-  OFX.OFX_vs_CTRL_significant %>%
  filter(direction == "Upregulated")
OFX.OFX_vs_CTRL_down <- OFX.OFX_vs_CTRL_significant %>%
  filter(direction == "Downregulated")
OFX.AMP_vs_CTRL_up <- OFX.AMP_vs_CTRL_significant %>%
  filter(direction == "Upregulated")
OFX.AMP_vs_CTRL_down <- OFX.AMP_vs_CTRL_significant %>%
  filter(direction == "Downregulated")

# 合并两个数据框
common_genes <- inner_join(
  OFX.OFX_vs_CTRL_significant,
  OFX.AMP_vs_CTRL_significant,
  by = c("name", "ID"),
  suffix = c("_OFX.OFX_vs_CTRL", "_OFX.AMP_vs_CTRL")
)

common_up <- inner_join(
  OFX.OFX_vs_CTRL_up,
  OFX.AMP_vs_CTRL_up,
  by = c("name", "ID", "significant", "direction"),
  suffix = c("_OFX.OFX_vs_CTRL", "_OFX.AMP_vs_CTRL")
)
common_down <- inner_join(
  OFX.OFX_vs_CTRL_down,
  OFX.AMP_vs_CTRL_down,
  by = c("name", "ID", "significant", "direction"),
  suffix = c("_OFX.OFX_vs_CTRL", "_OFX.AMP_vs_CTRL")
)

common_cross_1 <- inner_join(
  OFX.OFX_vs_CTRL_up,
  OFX.AMP_vs_CTRL_down,
  by = c("name", "ID", "significant"),
  suffix = c("_OFX.OFX_vs_CTRL", "_OFX.AMP_vs_CTRL")
)

common_cross_2 <- inner_join(
  OFX.OFX_vs_CTRL_down,
  OFX.AMP_vs_CTRL_up,
  by = c("name", "ID", "significant"),
  suffix = c("_OFX.OFX_vs_CTRL", "_OFX.AMP_vs_CTRL")
)

write.csv(common_up, "common_up_genes.csv", row.names = FALSE)
write.csv(common_down, "common_down_genes.csv", row.names = FALSE)
write.csv(common_cross_1, "OFX-OFX vs CTRL_up & OFX-AMP vs CTRL_down_genes.csv", row.names = FALSE)
write.csv(common_cross_2, "OFX-OFX vs CTRL_down & OFX-AMP vs CTRL_up_genes.csv", row.names = FALSE)

# 绘制Venn图
# 方式1：命名列表（最常用）
gene_list <- list(
  `OFX-OFX vs CTRL - Upregulated` = OFX.OFX_vs_CTRL_up$name,
  `OFX-OFX vs CTRL - Downregulated` = OFX.OFX_vs_CTRL_down$name,
  `OFX-AMP vs CTRL - Upregulated` = OFX.AMP_vs_CTRL_up$name,
  `OFX-AMP vs CTRL - Downregulated` = OFX.AMP_vs_CTRL_down$name
)

# 创建基础图形 
venn_up <- ggvenn(
  gene_list, 
  columns = c("OFX-OFX vs CTRL - Upregulated","OFX-AMP vs CTRL - Upregulated"), 
  fill_color = c("#ff5757", "#ffe25f"), 
  stroke_color = "#750000",
  show_percentage = TRUE,
  auto_scale = TRUE,
  set_name_size = 5,
  text_size = 4
)
venn_down <- ggvenn(
  gene_list, 
  columns = c("OFX-OFX vs CTRL - Downregulated","OFX-AMP vs CTRL - Downregulated"), 
  fill_color = c("#7b57ff", "#5fc2ff"), 
  show_percentage = TRUE,
  stroke_color = "#1a0077",
  auto_scale = TRUE,
  set_name_size = 5,
  text_size = 4
)
venn_cross_1 <- ggvenn(
  gene_list, 
  columns = c("OFX-OFX vs CTRL - Upregulated","OFX-AMP vs CTRL - Downregulated"), 
  fill_color = c("#ff5757", "#5fc2ff"), 
  show_percentage = TRUE,
  stroke_color = "#750000",
  auto_scale = TRUE,
  set_name_size = 5,
  text_size = 4
)
venn_cross_2 <- ggvenn(
  gene_list, 
  columns = c("OFX-OFX vs CTRL - Downregulated", "OFX-AMP vs CTRL - Upregulated"), 
  fill_color = c("#7b57ff", "#ffe25f"), 
  show_percentage = TRUE,
  stroke_color = "#1a0077",
  auto_scale = TRUE,
  set_name_size = 5,
  text_size = 4
)
# 保存Venn图
ggsave("venn_up.pdf", plot = venn_up, width = 10, height = 6, device = "pdf", dpi = 300)
ggsave("venn_down.pdf", plot = venn_down, width = 10, height = 6, device = "pdf", dpi = 300)
ggsave("venn_cross_1.pdf", plot = venn_cross_1, width = 10, height = 6, device = "pdf", dpi = 300)
ggsave("venn_cross_2.pdf", plot = venn_cross_2, width = 10, height = 6, device = "pdf", dpi = 300)

# 计算Jaccard系数
jaccard_up <- nrow(common_up) / length(union(OFX.OFX_vs_CTRL_up$name, OFX.AMP_vs_CTRL_up$name))
jaccard_down <- nrow(common_down) / length(union(OFX.OFX_vs_CTRL_down$name, OFX.AMP_vs_CTRL_down$name))
jaccard_cross_1 <- nrow(common_cross_1) / length(union(OFX.OFX_vs_CTRL_up$name, OFX.AMP_vs_CTRL_down$name))
jaccard_cross_2 <- nrow(common_cross_2) / length(union(OFX.OFX_vs_CTRL_down$name, OFX.AMP_vs_CTRL_up$name))

# 将结果组合为数据框
jaccard_results <- data.frame(
  Comparison = c("Upregulated", "Downregulated", "Cross_1(Up+Down)", "Cross_2(Down+Up)"),
  Jaccard_Index = c(jaccard_up, jaccard_down, jaccard_cross_1, jaccard_cross_2),
  Common_Genes = c(nrow(common_up), nrow(common_down), nrow(common_cross_1), nrow(common_cross_2))
)

write.csv(
  jaccard_results,
  file = "jaccard_similarity_results.csv",
  row.names = FALSE,
  quote = FALSE  # 避免字符串带引号
)
