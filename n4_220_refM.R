# === 自动配置和安装核心依赖包 ===
options(stringsAsFactors = FALSE)

# 2. 确保 BiocManager 可用
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


# 4. 定义你需要的所有包
required_pkgs <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "biomaRt",
  "DOSE",
  "AnnotationDbi",
  "pheatmap",
  "ggplot2",
  "patchwork",
  "ggrepel",
  "tidyverse"
)

# 5. 检查并安装缺失包
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing missing package: ", pkg))
    BiocManager::install(pkg, lib = Sys.getenv("R_LIBS_USER"), ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# 清除可能遗留的锁文件
unlink(paste0(.libPaths()[1], "/00LOCK*"), recursive = TRUE, force = TRUE)


# 设置工作路径
setwd("D:/Ljl_RNA_seq/n4_220_refM")
rm(list = ls())
options(stringsAsFactors = F)

# === 加载包 ===
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(ggrepel)
library(patchwork)

# === 1. 数据读取 ===
count_data <- read.csv("transcript_count_refM.csv", row.names = 1, check.names = FALSE)
group_data <- read.csv("group_refM.csv", row.names = 1)
group_data <- group_data[match(colnames(count_data), rownames(group_data)), , drop = FALSE]
group_data$condition <- relevel(factor(group_data$condition), ref = "Morphine")

# === 2. 构建DESeq2对象并过滤低表达基因 ===
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = group_data,
                              design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# === 3. 导出过滤后的原始计数矩阵 ===
filtered_counts <- counts(dds)
write.csv(filtered_counts, file = "filtered_counts.csv")

# === 4. 差异分析 ===
dds <- DESeq(dds)
print(resultsNames(dds))  # 查看系数名称

# === 5. 设置阈值 ===
fc_cutoff <- log2(1.2)
pval_cutoff <- 0.05

# === 5.1 导出转录本差异分析原始结果（含up/down）===

# Knock vs Morphine
res1 <- lfcShrink(dds, coef = "condition_Knock_vs_Morphine", type = "apeglm")
res1_ordered <- res1[order(res1$padj), ]
res1_df <- as.data.frame(res1_ordered)
res1_df$comparison <- "Knock_vs_Morphine"
res1_df$regulated <- "normal"
res1_df$regulated[res1_df$log2FoldChange > fc_cutoff & res1_df$padj < pval_cutoff] <- "up"
res1_df$regulated[res1_df$log2FoldChange < -fc_cutoff & res1_df$padj < pval_cutoff] <- "down"
res1_df$up_in_group <- "none"
res1_df$up_in_group[res1_df$regulated == "up"] <- "Knock"
res1_df$up_in_group[res1_df$regulated == "down"] <- "Morphine"
write.csv(res1_df, "DE_transcripts_Knock_vs_Morphine.csv", row.names = TRUE)

# Saline vs Morphine
res2 <- lfcShrink(dds, coef = "condition_Saline_vs_Morphine", type = "apeglm")
res2_ordered <- res2[order(res2$padj), ]
res2_df <- as.data.frame(res2_ordered)
res2_df$comparison <- "Saline_vs_Morphine"
res2_df$regulated <- "normal"
res2_df$regulated[res2_df$log2FoldChange > fc_cutoff & res2_df$padj < pval_cutoff] <- "up"
res2_df$regulated[res2_df$log2FoldChange < -fc_cutoff & res2_df$padj < pval_cutoff] <- "down"
res2_df$up_in_group <- "none"
res2_df$up_in_group[res2_df$regulated == "up"] <- "Saline"
res2_df$up_in_group[res2_df$regulated == "down"] <- "Morphine"
write.csv(res2_df, "DE_transcripts_Saline_vs_Morphine.csv", row.names = TRUE)

# === 6. 公用函数定义 ===
annotate_transcripts <- function(df) {
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  id_map <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_transcript_id",
                  values = rownames(df),
                  mart = mart)
  df$transcript_id <- rownames(df)
  df <- merge(df, id_map, by.x = "transcript_id", by.y = "ensembl_transcript_id", all.x = TRUE)
  return(df)
}

run_enrichment_separate <- function(df, prefix) {
  gene_list <- unique(na.omit(df$ensembl_gene_id[df$regulated != "normal"]))
  
  ego_BP <- enrichGO(gene = gene_list, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                     ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  ego_CC <- enrichGO(gene = gene_list, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                     ont = "CC", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  ego_MF <- enrichGO(gene = gene_list, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                     ont = "MF", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  
  write.csv(as.data.frame(ego_BP), paste0(prefix, "_GO_BP.csv"), row.names = FALSE)
  write.csv(as.data.frame(ego_CC), paste0(prefix, "_GO_CC.csv"), row.names = FALSE)
  write.csv(as.data.frame(ego_MF), paste0(prefix, "_GO_MF.csv"), row.names = FALSE)
  
  p_BP <- barplot(ego_BP, showCategory = 10, label_format = 100) + ggtitle("Biological Process")
  p_CC <- barplot(ego_CC, showCategory = 10, label_format = 100) + ggtitle("Cellular Component")
  p_MF <- barplot(ego_MF, showCategory = 10, label_format = 100) + ggtitle("Molecular Function")
  ggsave(paste0(prefix, "_GO_separate.png"), plot = p_BP / p_CC / p_MF, width = 10, height = 16)
  
  entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_list,
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first") %>% na.omit()
  ekegg <- enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)
  write.csv(as.data.frame(ekegg), paste0(prefix, "_KEGG.csv"), row.names = FALSE)
  
  p_kegg1 <- barplot(ekegg, showCategory = 10, label_format = 100) + ggtitle("KEGG Barplot")
  p_kegg2 <- dotplot(ekegg, showCategory = 10, label_format = 100) + ggtitle("KEGG Dotplot")
  ggsave(paste0(prefix, "_KEGG_plot.png"), plot = p_kegg1 / p_kegg2, width = 8, height = 10)
  
  final_plot <- (p_BP / p_CC) | (p_MF / p_kegg1)
  ggsave(paste0(prefix, "_GO_KEGG_combined.png"), plot = final_plot, width = 16, height = 9)
}

plot_volcano <- function(df, prefix, interest_genes = NULL) {
  # 前5个 padj 最小的差异基因（非NA）
  top20 <- df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 2.5) %>%
    arrange(padj) %>%
    slice_head(n = 20) %>%
    pull(transcript_id)
  
  # 感兴趣的基因转录本
  interest_ids <- names(interest_genes)
  interest_ids_in_df <- interest_ids[interest_ids %in% df$ensembl_gene_id]
  interest_transcripts <- df %>%
    filter(ensembl_gene_id %in% interest_ids_in_df) %>%
    pull(transcript_id)
  
  # 合并标注对象
  to_label <- unique(c(top20, interest_transcripts))
  
  # 标注文字：external_gene_name，后续转换为 expression
  df$label_text <- ifelse(df$transcript_id %in% to_label, df$external_gene_name, NA)
  
  # 构造 expression(italic(...)) 格式的字符，再 parse 成表达式向量
  label_exprs <- ifelse(!is.na(df$label_text),
                        paste0("italic('", df$label_text, "')"),
                        NA)
  df$label_expr <- NA
  df$label_expr[!is.na(label_exprs)] <- label_exprs[!is.na(label_exprs)]
  
  # 绘图
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = regulated)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = 2) +
    scale_color_manual(values = c("up" = "#BB0C00", "down" = "#00AFBB", "normal" = "grey")) +
    theme_bw() +
    labs(title = paste(prefix, "Volcano Plot"), x = "log2FoldChange", y = "-log10(padj)") +
    geom_label_repel(
      data = df[!is.na(df$label_expr), ],
      aes(label = label_expr),
      parse = TRUE, size = 5, max.overlaps = 15, show.legend = FALSE
    )
  
  ggsave(paste0(prefix, "_Volcano.png"), plot = p, width = 8, height = 7)
}




plot_heatmap <- function(df, prefix) {
  expr_cpm <- log10(counts(dds, normalized = TRUE) + 1)
  deg_ids <- df$transcript_id[df$regulated != "normal"]
  deg_ids <- deg_ids[deg_ids %in% rownames(expr_cpm)]
  mat <- expr_cpm[deg_ids, ]
  mat <- t(scale(t(mat)))
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  ann_col <- data.frame(group = group_data$condition)
  rownames(ann_col) <- rownames(group_data)
  
  pheatmap(mat,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           annotation_col = ann_col,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           main = paste("Heatmap of DE transcripts:", prefix),
           filename = paste0(prefix, "_Heatmap.png"),
           width = 8, height = 10)
}

# === 7. Knock vs Morphine ===
res1 <- lfcShrink(dds, coef = "condition_Knock_vs_Morphine", type = "apeglm")
res1_df <- as.data.frame(res1[order(res1$padj), ])
res1_df$comparison <- "Knock_vs_Morphine"
res1_df$regulated <- "normal"
res1_df$regulated[res1_df$log2FoldChange > fc_cutoff & res1_df$padj < pval_cutoff] <- "up"
res1_df$regulated[res1_df$log2FoldChange < -fc_cutoff & res1_df$padj < pval_cutoff] <- "down"
res1_df$up_in_group <- ifelse(res1_df$regulated == "up", "Knock",
                              ifelse(res1_df$regulated == "down", "Morphine", "none"))

# 注释
res1_df_annot <- annotate_transcripts(res1_df)

# 🔽 添加：导出转录本与基因的完整对应关系
id_map_res1 <- res1_df_annot[, c("transcript_id", "ensembl_gene_id", "external_gene_name")]
write.csv(id_map_res1, "Transcript_all_to_Gene_Knock_vs_Morphine.csv", row.names = FALSE)

# 去重：每个基因保留padj最小的转录本
#如需改为保留log2FC最大的转录本，只需将slice_min(order_by=padj)改为:slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE)
res1_df_annot_unique <- res1_df_annot %>%
  filter(!is.na(ensembl_gene_id)) %>%
  group_by(ensembl_gene_id) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# 保存去重后结果
write.csv(res1_df_annot_unique, "DE_transcripts_Knock_vs_Morphine_annotated_unique.csv", row.names = FALSE)

# 富集和可视化基于去重后
run_enrichment_separate(res1_df_annot_unique, "Knock_vs_Morphine")
interest_genes <- list(
  "ENSMUSG00000003534" = expression(italic("Ddr1")),
  "ENSMUSG00000020231" = expression(italic("Dip2a")),
  "ENSMUSG00000039191" = expression(italic("Rbpj")),
  "ENSMUSG00000055003" = expression(italic("Lrtm2")),
  "ENSMUSG00000027168" = expression(italic("Pax6")),
  "ENSMUSG00000027457" = expression(italic("Snph")),
  "ENSMUSG00000000223" = expression(italic("Drp2")),
  "ENSMUSG00000023868" = expression(italic("Pde10a")),
  "ENSMUSG00000034135" = expression(italic("Sik3")),
  "ENSMUSG00000028176" = expression(italic("Lrrc7")),
  "ENSMUSG00000053477" = expression(italic("Tcf4")),
  "ENSMUSG00000025485" = expression(italic("Ric8a")),
  "ENSMUSG00000032409" = expression(italic("Atr")),
  "ENSMUSG00000025280" = expression(italic("Polr3a")),
  "ENSMUSG00000037376" = expression(italic("Trmt6")),
  "ENSMUSG00000040055" = expression(italic("Gjb6")),
  "ENSMUSG00000044149" = expression(italic("Nkrf")),
  "ENSMUSG00000062762" = expression(italic("Ei24"))
)

plot_volcano(res1_df_annot_unique, "Knock_vs_Morphine", interest_genes)

plot_heatmap(res1_df_annot_unique, "Knock_vs_Morphine")

# === 8. Saline vs Morphine ===
res2 <- lfcShrink(dds, coef = "condition_Saline_vs_Morphine", type = "apeglm")
res2_df <- as.data.frame(res2[order(res2$padj), ])
res2_df$comparison <- "Saline_vs_Morphine"
res2_df$regulated <- "normal"
res2_df$regulated[res2_df$log2FoldChange > fc_cutoff & res2_df$padj < pval_cutoff] <- "up"
res2_df$regulated[res2_df$log2FoldChange < -fc_cutoff & res2_df$padj < pval_cutoff] <- "down"
res2_df$up_in_group <- ifelse(res2_df$regulated == "up", "Saline",
                              ifelse(res2_df$regulated == "down", "Morphine", "none"))
# 注释
res2_df_annot <- annotate_transcripts(res2_df)

# 🔽 添加：导出转录本与基因的完整对应关系
id_map_res2 <- res2_df_annot[, c("transcript_id", "ensembl_gene_id", "external_gene_name")]
write.csv(id_map_res2, "Transcript_all_to_Gene_Saline_vs_Morphine.csv", row.names = FALSE)

# 🔽 合并两个对应表（保留全部转录本）
#id_map_res1$comparison <- "Knock_vs_Morphine"
#id_map_res2$comparison <- "Saline_vs_Morphine"
#id_map_all <- bind_rows(id_map_res1, id_map_res2)
#id_map_all <- distinct(id_map_all)
#write.csv(id_map_all, "Transcript_to_Gene_All.csv", row.names = FALSE)

# === 去重 ===
res2_df_annot_unique <- res2_df_annot %>%
  filter(!is.na(ensembl_gene_id)) %>%
  group_by(ensembl_gene_id) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(res2_df_annot_unique, "DE_transcripts_Saline_vs_Morphine_annotated_unique.csv", row.names = FALSE)

# 富集和图
run_enrichment_separate(res2_df_annot_unique, "Saline_vs_Morphine")
interest_genes <- list(
  "ENSMUSG00000003534" = expression(italic("Ddr1")),
  "ENSMUSG00000020231" = expression(italic("Dip2a")),
  "ENSMUSG00000039191" = expression(italic("Rbpj")),
  "ENSMUSG00000055003" = expression(italic("Lrtm2")),
  "ENSMUSG00000027168" = expression(italic("Pax6")),
  "ENSMUSG00000027457" = expression(italic("Snph")),
  "ENSMUSG00000000223" = expression(italic("Drp2")),
  "ENSMUSG00000023868" = expression(italic("Pde10a")),
  "ENSMUSG00000034135" = expression(italic("Sik3")),
  "ENSMUSG00000028176" = expression(italic("Lrrc7")),
  "ENSMUSG00000053477" = expression(italic("Tcf4")),
  "ENSMUSG00000025485" = expression(italic("Ric8a")),
  "ENSMUSG00000032409" = expression(italic("Atr")),
  "ENSMUSG00000025280" = expression(italic("Polr3a")),
  "ENSMUSG00000037376" = expression(italic("Trmt6")),
  "ENSMUSG00000040055" = expression(italic("Gjb6")),
  "ENSMUSG00000044149" = expression(italic("Nkrf")),
  "ENSMUSG00000062762" = expression(italic("Ei24"))
)

plot_volcano(res2_df_annot_unique, "Saline_vs_Morphine", interest_genes)

plot_heatmap(res2_df_annot_unique, "Saline_vs_Morphine")



