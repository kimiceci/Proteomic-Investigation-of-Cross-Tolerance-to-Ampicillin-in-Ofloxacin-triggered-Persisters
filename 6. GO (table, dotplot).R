# 加载必要的包
library(xml2)
library(ggplot2)
library(ggrepel)  # 智能标签防重叠
library(scales)    # 非线性色阶支持 
library(dplyr)

# 设置工作目录
setwd("D:/Final Report")

contrast <- c("OFX-OFX vs CTRL_up","OFX-OFX vs CTRL_down","OFX-AMP vs CTRL_up","OFX-AMP vs CTRL_down")

annotation <- c("BP","MF","CC")

# 提取xml文件中level0的项并保存为csv
getlv0 <- function(filePath) {
  # 读取 XML 文件
  lines <- readLines(filePath, warn = FALSE)
  xml_data <- read_xml(paste(lines, collapse = "\n"))
  # 提取所有 level=0 的 result 节点
  level0_results <- xml_find_all(xml_data, ".//result[term/level='0']")
  # 提取每个 result 的完整信息
  results_list <- lapply(level0_results, function(node) {
    list(
      term_id = xml_text(xml_find_first(node, ".//term/id")),
      term_label = xml_text(xml_find_first(node, ".//term/label")),
      number_in_reference = as.numeric(xml_text(xml_find_first(node, ".//number_in_reference"))),
      number_in_list = as.numeric(xml_text(xml_find_first(node, ".//input_list/number_in_list"))),
      expected = as.numeric(xml_text(xml_find_first(node, ".//input_list/expected"))),
      pValue = as.numeric(xml_text(xml_find_first(node, ".//input_list/pValue"))),
      fold_enrichment = as.numeric(xml_text(xml_find_first(node, ".//input_list/fold_enrichment"))),
      fdr = as.numeric(xml_text(xml_find_first(node, ".//input_list/fdr")))
    )
  })
  # 转换为数据框
  results_df <- do.call(rbind.data.frame, results_list)
  # 提取 mapped_ids（保持为列表）
  mapped_ids_list <- lapply(level0_results, function(node) {
    xml_text(xml_find_all(node, ".//input_list/mapped_id_list/mapped_id"))
  })
  # 将 mapped_ids 添加到数据框中（可选）
  results_df$mapped_ids <- mapped_ids_list
  # 将列表列转换为逗号分隔的字符串
  results_df_str <- results_df %>%
    mutate(mapped_ids_str = sapply(mapped_ids, paste, collapse = ", "))
  return(results_df_str)
}

drawDotplot <- function(df, plotTitle, plotSubtitle, savePath){
  # 完整绘图代码 
  dot_plot <- ggplot(df,
                    aes(x = fold_enrichment,
                        y = term_label,
                        size = number_in_list,
                        color = -log10(fdr)))  +
    # Core graphical elements 
    geom_point(alpha = 0.85, shape = 19) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey40") +
    annotate("text",
         x = 1.5, y = 0,
         label = "x=1.5", 
         vjust = -1, hjust = -0.2,  # 水平向右偏移
         color = "grey40", size = 3)+
    # Significance markers (non-overlapping)
    ggrepel::geom_text_repel(
      aes(label = ifelse(fdr  < 0.05, 
                        ifelse(fdr  < 0.001, "***", 
                                ifelse(fdr  < 0.01, "**", "*")), "")),
      color = "black",
      fontface = "bold",
      size = 3.5,
      nudge_x = 0.3,
      direction = "y",
      segment.color  = "grey70",
      max.overlaps  = 20,
      point.padding = 0.5,
      box.padding = 0.3
    ) +
    
    # Visual mapping system 
    scale_size_continuous(
      name = "Gene Count",
      range = c(4, 10),
      breaks = seq(5, 30, by = 5)
    ) +
    scale_color_gradientn(
      colors = c("#0f72bd", "#29caff", "#ffd428", "#ff8810", "#d7191c"),
      values = rescale(c(0, 1, 2, 5, 10)),
      name = expression("-log"[10] ~ "p.adj"),
      limits = c(0, max(-log10(df$fdr))) 
    ) +
    
    # Labeling system
    labs(
      x = "Fold Enrichment", 
      y = NULL,
      title = plotTitle,
      subtitle = plotSubtitle
    ) +
    
    # Unified theme settings 
    theme_bw(base_size = 13) +
    theme(
      # Axis customization
      axis.text.y  = element_text(
        face = "italic", 
        size = 10,
        margin = margin(r = 5, unit = "pt")
      ),
      axis.title.x  = element_text(
        face = "bold", 
        size = 12,
        margin = margin(t = 10, unit = "pt")
      ),
      
      # Title customization
      plot.title  = element_text(
        face = "bold", 
        hjust = 0.5,
        size = 16,
        margin = margin(b = 10, unit = "pt")
      ),
      plot.title.position  = "plot",
      plot.subtitle  = element_text(
        hjust = 0.5, 
        color = "grey40",
        size = 11
      ),
      plot.caption  = element_text(
        color = "grey50",
        size = 9 
      ),
      
      # Grid customization 
      panel.grid.minor  = element_blank(),
      panel.grid.major.y  = element_line(
        linetype = "dotted",
        color = "grey90"
      ),
      
      # Legend system (horizontal bottom layout)
      legend.position  = "bottom",
      legend.justification = "right",
      legend.box  = "horizontal",
      legend.direction  = "horizontal",
      legend.title  = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.spacing.x  = unit(0.5, 'cm'),
      
      # Margin system (3cm left margin for pathway names)
      plot.margin  = margin(1, 1, 1, 0.5, "cm")
    ) +
    
    # Legend refinement 
    guides(
      color = guide_colorbar(
        title.position  = "top",
        barwidth = unit(5, "cm"),
        barheight = unit(0.4, "cm"),
        direction = "horizontal",
        ticks.colour  = "black"
      ),
      size = guide_legend(
        title.position  = "top",
        nrow = 1,
        override.aes  = list(alpha = 1)
      )
    )
  save_path <- savePath
  ggsave(savePath, plot = dot_plot, device = "pdf", width = 9, height = 8, dpi = 300)
}

for (cont in contrast){
  geneFileName <- paste0(substr(cont,1,15),"_genes.csv")
  gene_df <- read.csv(geneFileName, header = TRUE)
  reg_direction <- ifelse(substr(cont,17,18)=="up", "Upregulated", ifelse(substr(cont,17,20)=="down", "Downregulated","None"))
  gene_list <- gene_df[gene_df$significant & gene_df$direction==reg_direction, ]$name
  cluster_list <- list()
  for (anno in annotation){
    message(paste0("Processing ", anno, " annotation for ", cont, "..."))
    # 读取GO富集分析结果文件
    go_file_path <- paste0(cont, "_GO_",anno,"_Fisher_FDR.xml")
    go_save_path <- paste0(cont, "_GO_",anno,"_Fisher_FDR_lv0.csv")
    results_df <-getlv0(go_file_path)
    # 保存为 CSV（不保留原始列表列）
    write.csv(select(results_df, -mapped_ids), go_save_path, row.names = FALSE)
    
    results_df <- results_df %>%
      arrange(fdr)  %>%  # 核心修改：按校正p值升序排序 
      head(15)               # 提取Top15最显著通路 
    
    plot_title <- paste0(substr(cont,1,15), " GO ", switch(anno, BP="Biological Process",MF="Molecular Function",CC="Cellular Component"), " Enrichment Analysis")
    plot_subtitle <- paste0(reg_direction," | Significant genes: ", length(gene_list)," | Top ", nrow(results_df), " pathways")
    plot_save_path <- paste0(cont, "_GO_",anno,"_dotplot.pdf")

    drawDotplot(
      df = results_df,
      plotTitle = plot_title,
      plotSubtitle = plot_subtitle,
      savePath = plot_save_path
    )
    message(paste0(cont, " GO ", anno," enrichment analysis completed."))

# 如不需要对Cluster进行GO富集分析，可以注释掉以下代码
    for (i in 1:13){
      message(paste0("Processing ", cont, " - Cluster ", i, " for ", anno, " annotation..."))
      cluster_file_path <- paste0(cont, "_Cluster ", i,"_GO_",anno,"_Fisher_FDR.xml")
      if (!file.exists(cluster_file_path)) {
        message(paste("文件不存在，跳过:", cluster_file_path))
        break  # 跳过余下循环
        }
      cluster_df <-getlv0(cluster_file_path)
      if (nrow(cluster_df) == 0) {
        message(paste("文件为空，跳过:", cluster_file_path))
        next  # 跳过本次循环      
        }
      cluster_list[[paste0(cont, "_Cluster_", i)]][[anno]] <- cluster_df
    }
  }
  clusters_genes <- read.csv(file.path("Hub analysis (Cytoscape)", "MCODE", paste0(cont, "_clusters.csv")), header = TRUE, na.strings = c(""))
  for (i in 1:length(cluster_list)) {
    cluster_length <- length(na.omit(clusters_genes[[paste0("Cluster.",i)]]))
    cluster_combined_df <- bind_rows(cluster_list[[i]], .id = "Annotation")  # 合并所有集群数据
    cluster_save_path <- paste0(cont, "_Cluster ", i,"_GO_Fisher_FDR_lv0.csv")
    # 保存为 CSV（不保留原始列表列）
    write.csv(select(cluster_combined_df, -mapped_ids), cluster_save_path, row.names = FALSE)
    cluster_combined_df <- cluster_combined_df %>%
      mutate(term_label_original = cluster_combined_df$term_label)
    sorted_df <- cluster_combined_df %>%
      mutate(
        Annotation = factor(Annotation, levels = rev(annotation)),  # 再次确保因子顺序
        term_label = paste0(term_label, " (", Annotation, ")")
      ) %>%
      arrange(match(Annotation, rev(annotation)), -fdr) %>%  # 按指定annotation顺序排序
      head(15)
    # 关键步骤：重新生成term_label的因子水平
    sorted_df$term_label <- factor(
      sorted_df$term_label, 
      levels = unique(sorted_df$term_label)  # 保持排序后的顺序
    )
    plot_title <- paste0(substr(cont,1,15), " GO Enrichment Analysis")
    plot_subtitle <- paste0(reg_direction,"-Cluster ", i, " | Significant genes: ", cluster_length," | Top ", nrow(sorted_df), " pathways")
    plot_save_path <- paste0(cont, "_Cluster ", i,"_GO_dotplot.pdf")
    drawDotplot(
        df = sorted_df,
        plotTitle = plot_title,
        plotSubtitle = plot_subtitle,
        savePath = plot_save_path
      )
    message(paste0(cont, " - Cluster ", i, " GO enrichment analysis completed."))
# 如不需要对Cluster进行GO富集分析，可以注释掉以上代码
  }
}
