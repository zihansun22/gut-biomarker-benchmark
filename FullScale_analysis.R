library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 读取数据 ──────────────────────────────────────────────
meta_raw   <- read.csv("valid_matched.csv", sep = ",", stringsAsFactors = FALSE)
counts_df  <- read_tsv("counts_genus.tsv", show_col_types = FALSE)

# ── 整理矩阵：行 = taxa，列 = 样本 ───────────────────────
taxa_names <- counts_df[[1]]
count_mat  <- as.matrix(counts_df[, -1])
rownames(count_mat) <- taxa_names

# 只保留 metadata 中存在的样本
keep_samples <- intersect(colnames(count_mat), meta_raw$run_id)
count_mat    <- count_mat[, keep_samples]

cat("矩阵维度：", nrow(count_mat), "taxa ×", ncol(count_mat), "samples\n")

# ════════════════════════════════════════════════════════════
# 1. 稀疏度汇总
# ════════════════════════════════════════════════════════════

## 1-a  Overall sparsity
overall_sparsity <- mean(count_mat == 0)
cat("\n── 1. 稀疏度汇总 ──────────────────────────────────────\n")
cat(sprintf("Overall sparsity (零值比例): %.4f (%.2f%%)\n",
            overall_sparsity, overall_sparsity * 100))

## 1-b  每个 taxa 的 prevalence（在多少样本中出现）
taxa_prevalence <- rowMeans(count_mat > 0)   # 0~1，出现比例
prev_quantiles  <- quantile(taxa_prevalence, probs = c(0.10, 0.50, 0.90))
cat("\nTaxa prevalence 分位数（出现样本比例）:\n")
cat(sprintf("  10%%: %.4f\n  50%%: %.4f\n  90%%: %.4f\n",
            prev_quantiles[1], prev_quantiles[2], prev_quantiles[3]))

## 1-c  每个 taxa 的 zero_prop（零值比例 = 1 - prevalence）
taxa_zero_prop <- 1 - taxa_prevalence
cat(sprintf("\nTaxa zero_prop 均值  : %.4f\n", mean(taxa_zero_prop)))
cat(sprintf("Taxa zero_prop 中位数: %.4f\n", median(taxa_zero_prop)))

# ════════════════════════════════════════════════════════════
# 2. Library size（测序深度）分布
# ════════════════════════════════════════════════════════════

lib_size <- colSums(count_mat)
cat("\n── 2. Library size 分布 ───────────────────────────────\n")
cat(sprintf("Min    : %s\n", format(min(lib_size),    big.mark = ",")))
cat(sprintf("Q1     : %s\n", format(quantile(lib_size, 0.25), big.mark = ",")))
cat(sprintf("Median : %s\n", format(median(lib_size), big.mark = ",")))
cat(sprintf("Mean   : %s\n", format(round(mean(lib_size)), big.mark = ",")))
cat(sprintf("Q3     : %s\n", format(quantile(lib_size, 0.75), big.mark = ",")))
cat(sprintf("Max    : %s\n", format(max(lib_size),    big.mark = ",")))
cat(sprintf("IQR    : %s\n", format(IQR(lib_size),    big.mark = ",")))

## 直方图
p_lib <- ggplot(data.frame(lib_size = lib_size), aes(x = lib_size)) +
  geom_histogram(bins = 40, fill = "#4E79A7", color = "white", alpha = 0.85) +
  geom_vline(xintercept = median(lib_size), linetype = "dashed",
             color = "#E15759", linewidth = 0.8) +
  annotate("text", x = median(lib_size), y = Inf,
           label = paste0("Median=", format(round(median(lib_size)), big.mark = ",")),
           hjust = -0.1, vjust = 1.5, color = "#E15759", size = 3.5) +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Library Size Distribution",
       x = "Library Size (read counts)", y = "Number of Samples") +
  theme_bw(base_size = 12)

ggsave("real_data_feature/library_size_histogram.pdf", p_lib, width = 7, height = 4)
ggsave("real_data_feature/library_size_histogram.png", p_lib, width = 7, height = 4, dpi = 150)
cat("\n已保存: library_size_histogram.pdf / .png\n")

# ════════════════════════════════════════════════════════════
# 3. 丰度长尾证据
# ════════════════════════════════════════════════════════════

## 相对丰度矩阵（列归一化）
rel_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/")

## 每个 taxa 的平均相对丰度，降序排列
mean_rel_abund <- sort(rowMeans(rel_mat), decreasing = TRUE)

## Top-20 累积占比
top20_cumsum <- cumsum(mean_rel_abund[1:min(20, length(mean_rel_abund))])
cat("\n── 3. 丰度长尾证据 ────────────────────────────────────\n")
cat("Top-N taxa 平均相对丰度累积占比:\n")
for (n in c(1, 5, 10, 20)) {
  if (n <= length(top20_cumsum))
    cat(sprintf("  Top-%2d: %.2f%%\n", n, top20_cumsum[n] * 100))
}

## Rank-abundance 曲线（全部 taxa）
rank_df <- data.frame(
  rank       = seq_along(mean_rel_abund),
  mean_relab = mean_rel_abund
)

p_rank <- ggplot(rank_df, aes(x = rank, y = mean_relab)) +
  geom_line(color = "#4E79A7", linewidth = 0.7) +
  geom_vline(xintercept = 20, linetype = "dashed",
             color = "#E15759", linewidth = 0.7) +
  annotate("text", x = 22, y = max(mean_rel_abund) * 0.9,
           label = sprintf("Top20 = %.1f%%", top20_cumsum[min(20, length(top20_cumsum))] * 100),
           hjust = 0, color = "#E15759", size = 3.5) +
  scale_y_log10(labels = scales::percent_format(accuracy = 0.001)) +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Rank-Abundance Curve (Mean Relative Abundance)",
       x = "Taxa Rank", y = "Mean Relative Abundance (log10)") +
  theme_bw(base_size = 12)

ggsave("real_data_feature/rank_abundance_curve.pdf", p_rank, width = 7, height = 4)
ggsave("real_data_feature/rank_abundance_curve.png", p_rank, width = 7, height = 4, dpi = 150)
cat("已保存: rank_abundance_curve.pdf / .png\n")

## Top-20 条形图（可选，直观展示头部 taxa）
top20_df <- data.frame(
  taxa      = factor(names(mean_rel_abund)[1:min(20, length(mean_rel_abund))],
                     levels = rev(names(mean_rel_abund)[1:min(20, length(mean_rel_abund))])),
  mean_relab = mean_rel_abund[1:min(20, length(mean_rel_abund))]
)

p_top20 <- ggplot(top20_df, aes(x = mean_relab, y = taxa)) +
  geom_col(fill = "#76B7B2", alpha = 0.9) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Top-20 Taxa by Mean Relative Abundance",
       x = "Mean Relative Abundance", y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave("real_data_feature/top20_taxa_barplot.pdf", p_top20, width = 8, height = 6)
ggsave("real_data_feature/top20_taxa_barplot.png", p_top20, width = 8, height = 6, dpi = 150)
cat("已保存: top20_taxa_barplot.pdf / .png\n")



# ══════════════════════════════════════════════════════════════
# 加载包
# ══════════════════════════════════════════════════════════════
library(tidyverse)

QVAL_THR <- 0.05

# ══════════════════════════════════════════════════════════════
# 工具函数：从完整分类路径提取 g__ 属名
# ══════════════════════════════════════════════════════════════
extract_genus_tag <- function(x) {
  # 匹配 g__ 开头的最后一个字段（兼容点、分号、空格分隔）
  # 输入示例:
  #   "d__Bacteria.p__Firmicutes...g__Streptococcus"
  #   "d__Archaea;p__...;g__Methanobrevibacter"
  #   "g__Streptococcus"  （LEfSe已是此格式，直接返回）
  sapply(x, function(s) {
    # 统一分隔符为空格，然后按 g__ 切割
    m <- str_extract(s, "g__[^.;|\\s]+")
    if (is.na(m)) NA_character_ else m
  }, USE.NAMES = FALSE)
}

# ══════════════════════════════════════════════════════════════
# 0. 读入所有原始数据
# ══════════════════════════════════════════════════════════════

# ---------- 元数据 ----------
meta_raw <- read.csv("valid_matched.csv", sep = ",",
                     stringsAsFactors = FALSE)
# 文件是逗号分隔

# ---------- LEfSe ----------
lefse_df <- read_csv(
  "LEfSe/lefse_results_2026_04/lefse_summary_raw.csv",
  show_col_types = FALSE
)

# ---------- MaAsLin2 ----------
maaslin_raw <- read_csv(
  "MaAsLin2/maaslin_summary_raw.csv",
  show_col_types = FALSE
)

# ---------- ANCOM-BC2 ----------
# 文件有列名，直接读取
ancom_raw <- read_csv(
  "ANCOM_BC2/ancombc_summary.csv",
  show_col_types = FALSE
)
# 列: taxon | lfc | se | W | p_val | q_val | diff_abn |
#      comparison | group_a | group_b | ref_group

cat("数据读入完成\n")
cat("LEfSe 行数:", nrow(lefse_df), "\n")
cat("MaAsLin2 行数:", nrow(maaslin_raw), "\n")
cat("ANCOM-BC2 行数:", nrow(ancom_raw), "\n")

# ══════════════════════════════════════════════════════════════
# 1. 分析队列数
# ══════════════════════════════════════════════════════════════

n_cohorts <- meta_raw %>%
  distinct(project_id) %>%
  nrow()

cat("\n① 分析队列数:", n_cohorts, "\n")

# ══════════════════════════════════════════════════════════════
# 2. 覆盖疾病数
# ══════════════════════════════════════════════════════════════

all_phenotypes <- meta_raw %>%
  distinct(phenotype.x) %>%
  pull(phenotype.x)

n_diseases <- length(all_phenotypes)

cat("\n② 覆盖疾病数:", n_diseases, "\n")
cat("   疾病列表:", paste(all_phenotypes, collapse = ", "), "\n")

# ══════════════════════════════════════════════════════════════
# 3. 总样本量
# ══════════════════════════════════════════════════════════════

n_samples <- meta_raw %>%
  distinct(run_id) %>%
  nrow()

cat("\n③ 总样本量:", n_samples, "\n")

# ══════════════════════════════════════════════════════════════
# 4. 提取各方法显著特征（g__ 格式统一）
# ══════════════════════════════════════════════════════════════

# ---------- LEfSe 显著特征 ----------
lefse_sig <- lefse_df %>%
  filter(qval < QVAL_THR) %>%                        # ← padj 改为 qval
  mutate(genus_tag = extract_genus_tag(taxa_full)) %>% # ← 从 taxa_full 提取 g__
  filter(!is.na(genus_tag))

# 每比较对显著特征数
# LEfSe 文件没有 project_id 列，直接用 comparison 分组
lefse_per_comp <- lefse_sig %>%
  group_by(comparison) %>%
  summarise(n_features = n_distinct(genus_tag), .groups = "drop")

lefse_mean  <- mean(lefse_per_comp$n_features)
lefse_total <- n_distinct(lefse_sig$genus_tag)

cat("\n④-LEfSe:\n")
cat("   显著特征总数（去重 genus）:", lefse_total, "\n")
cat("   每比较对平均显著特征数:", round(lefse_mean, 2), "\n")

# ---------- MaAsLin2 显著特征 ----------
# feature 列是完整路径，提取 g__ 部分
maaslin_sig <- maaslin_raw %>%
  filter(qval < QVAL_THR) %>%
  mutate(genus_tag = extract_genus_tag(feature)) %>%
  filter(!is.na(genus_tag))            # 去除无法提取属名的行

maaslin_per_comp <- maaslin_sig %>%
  group_by(comparison) %>%
  summarise(n_features = n_distinct(genus_tag), .groups = "drop")

maaslin_mean  <- if (nrow(maaslin_per_comp) > 0)
  mean(maaslin_per_comp$n_features) else 0
maaslin_total <- n_distinct(maaslin_sig$genus_tag)

cat("\n④-MaAsLin2:\n")
cat("   显著特征总数（去重 genus）:", maaslin_total, "\n")
cat("   每比较对平均显著特征数:", round(maaslin_mean, 2), "\n")

# ---------- ANCOM-BC2 显著特征 ----------
# taxon 列是完整路径，提取 g__ 部分
ancom_sig <- ancom_raw %>%
  filter(q_val < QVAL_THR) %>%
  mutate(genus_tag = extract_genus_tag(taxon)) %>%
  filter(!is.na(genus_tag))

ancom_per_comp <- ancom_sig %>%
  group_by(comparison) %>%
  summarise(n_features = n_distinct(genus_tag), .groups = "drop")

ancom_mean  <- if (nrow(ancom_per_comp) > 0)
  mean(ancom_per_comp$n_features) else 0
ancom_total <- n_distinct(ancom_sig$genus_tag)

cat("\n④-ANCOM-BC2:\n")
cat("   显著特征总数（去重 genus）:", ancom_total, "\n")
cat("   每比较对平均显著特征数:", round(ancom_mean, 2), "\n")

# ══════════════════════════════════════════════════════════════
# 5. 与 LEfSe 的总体重叠率
#    定义：总体重叠 = 各方法与 LEfSe 的 Jaccard 系数
#    Jaccard = |A ∩ B| / |A ∪ B|
# ══════════════════════════════════════════════════════════════

lefse_genera   <- unique(lefse_sig$genus_tag)
maaslin_genera <- unique(maaslin_sig$genus_tag)
ancom_genera   <- unique(ancom_sig$genus_tag)

# ── MaAsLin2 vs LEfSe ──────────────────────────────────────
overlap_ml  <- intersect(maaslin_genera, lefse_genera)
union_ml    <- union(maaslin_genera, lefse_genera)
jaccard_ml  <- length(overlap_ml) / length(union_ml)

# ── ANCOM-BC2 vs LEfSe ─────────────────────────────────────
overlap_al  <- intersect(ancom_genera, lefse_genera)
union_al    <- union(ancom_genera, lefse_genera)
jaccard_al  <- length(overlap_al) / length(union_al)

# ── ANCOM-BC2 vs MaAsLin2 ────────────────────────────────────
overlap_am <- intersect(ancom_genera, maaslin_genera)
union_am   <- union(ancom_genera, maaslin_genera)
jaccard_am <- length(overlap_am) / length(union_am)

cat("   ANCOM-BC2 ∩ MaAsLin2:", length(overlap_am),
    "| Jaccard:", round(jaccard_am, 4),
    paste0("(", round(jaccard_am * 100, 1), "%)"), "\n")
cat("   ANCOM-BC2 ∩ MaAsLin2 重叠属:", paste(overlap_am, collapse = ", "), "\n\n")

# ── 三方法总体重叠（三者交集 / 三者并集）──────────────────
overlap_all <- intersect(intersect(maaslin_genera, ancom_genera),
                         lefse_genera)
union_all   <- union(union(maaslin_genera, ancom_genera),
                     lefse_genera)
jaccard_all <- length(overlap_all) / length(union_all)

cat("\n⑤ 与 LEfSe 的总体重叠率（Jaccard）:\n")
cat("   LEfSe 显著属数:", length(lefse_genera), "\n")
cat("   MaAsLin2 显著属数:", length(maaslin_genera), "\n")
cat("   ANCOM-BC2 显著属数:", length(ancom_genera), "\n\n")

cat("   MaAsLin2 ∩ LEfSe:", length(overlap_ml),
    "| Jaccard:", round(jaccard_ml, 4),
    paste0("(", round(jaccard_ml * 100, 1), "%)"), "\n")
cat("   MaAsLin2 ∩ LEfSe 重叠属:", paste(overlap_ml, collapse = ", "), "\n\n")

cat("   ANCOM-BC2 ∩ LEfSe:", length(overlap_al),
    "| Jaccard:", round(jaccard_al, 4),
    paste0("(", round(jaccard_al * 100, 1), "%)"), "\n")
cat("   ANCOM-BC2 ∩ LEfSe 重叠属:", paste(overlap_al, collapse = ", "), "\n\n")

cat("   三方法总体重叠（交集/并集）:", length(overlap_all), "/", length(union_all),
    "| Jaccard:", round(jaccard_all, 4),
    paste0("(", round(jaccard_all * 100, 1), "%)"), "\n")
cat("   三方法共同重叠属:", paste(overlap_all, collapse = ", "), "\n")

# ══════════════════════════════════════════════════════════════
# 5.5 Top-K（K=50）全局重合率：按属的最小 q 值排序取 Top50
# ══════════════════════════════════════════════════════════════

TOPK <- 50

get_topk_genera <- function(df, genus_col, q_col, k = 50) {
  genus_col <- rlang::ensym(genus_col)
  q_col     <- rlang::ensym(q_col)
  
  df %>%
    filter(!is.na(!!genus_col), !is.na(!!q_col)) %>%
    group_by(!!genus_col) %>%
    summarise(best_q = min(!!q_col), .groups = "drop") %>%
    arrange(best_q) %>%
    slice_head(n = k) %>%
    pull(!!genus_col) %>%
    unique()
}

calc_overlap_metrics <- function(A, B, k = 50) {
  inter <- intersect(A, B)
  uni   <- union(A, B)
  list(
    overlap_n   = length(inter),
    overlap_atk = length(inter) / k,
    jaccard_atk = length(inter) / length(uni)
  )
}

# 各方法 Top50 属（全局）
lefse_top50   <- get_topk_genera(lefse_sig,   genus_tag, qval,  k = TOPK)
maaslin_top50 <- get_topk_genera(maaslin_sig, genus_tag, qval,  k = TOPK)
ancom_top50   <- get_topk_genera(ancom_sig,   genus_tag, q_val, k = TOPK)

# Pairwise overlap@50 和 Jaccard@50
m_ml <- calc_overlap_metrics(maaslin_top50, lefse_top50,   k = TOPK)
m_al <- calc_overlap_metrics(ancom_top50,   lefse_top50,   k = TOPK)
m_am <- calc_overlap_metrics(ancom_top50,   maaslin_top50, k = TOPK)

# 三方法 Top50 交集/并集
top50_overlap_all   <- Reduce(intersect, list(lefse_top50, maaslin_top50, ancom_top50))
top50_union_all     <- Reduce(union,     list(lefse_top50, maaslin_top50, ancom_top50))
top50_jaccard_all   <- length(top50_overlap_all) / length(top50_union_all)

cat("\n5.5 Top50 重合率（全局 genus 最小q 排序）:\n")
cat("   LEfSe Top50 实际属数:", length(lefse_top50), "\n")
cat("   MaAsLin2 Top50 实际属数:", length(maaslin_top50), "\n")
cat("   ANCOM-BC2 Top50 实际属数:", length(ancom_top50), "\n\n")

cat("   MaAsLin2 vs LEfSe:",
    "overlap_n =", m_ml$overlap_n,
    "| overlap@50 =", paste0(round(m_ml$overlap_atk * 100, 2), "%"),
    "| Jaccard@50 =", paste0(round(m_ml$jaccard_atk * 100, 2), "%"), "\n")

cat("   ANCOM-BC2 vs LEfSe:",
    "overlap_n =", m_al$overlap_n,
    "| overlap@50 =", paste0(round(m_al$overlap_atk * 100, 2), "%"),
    "| Jaccard@50 =", paste0(round(m_al$jaccard_atk * 100, 2), "%"), "\n")

cat("   ANCOM-BC2 vs MaAsLin2:",
    "overlap_n =", m_am$overlap_n,
    "| overlap@50 =", paste0(round(m_am$overlap_atk * 100, 2), "%"),
    "| Jaccard@50 =", paste0(round(m_am$jaccard_atk * 100, 2), "%"), "\n")

cat("   三方法 Top50 Jaccard:",
    paste0(round(top50_jaccard_all * 100, 2), "%"),
    "| 三方法共同 Top50 属数:", length(top50_overlap_all), "\n")

# ══════════════════════════════════════════════════════════════
# 6. 汇总报告
# ══════════════════════════════════════════════════════════════

summary_report <- tibble(
  指标 = c(
    "分析队列数",
    "覆盖疾病数",
    "总样本量",
    "LEfSe 显著属总数（去重）",
    "LEfSe 每比较对平均显著属数",
    "MaAsLin2 显著属总数（去重）",
    "MaAsLin2 每比较对平均显著属数",
    "ANCOM-BC2 显著属总数（去重）",
    "ANCOM-BC2 每比较对平均显著属数",
    "MaAsLin2 vs LEfSe Jaccard（%）",
    "ANCOM-BC2 vs LEfSe Jaccard（%）",
    "ANCOM-BC2 vs MaAsLin2 Jaccard（%）",
    "三方法总体 Jaccard（%）",
    "MaAsLin2 vs LEfSe Top50 overlap@50（%）",
    "ANCOM-BC2 vs LEfSe Top50 overlap@50（%）",
    "ANCOM-BC2 vs MaAsLin2 Top50 overlap@50（%）",
    "三方法 Top50 Jaccard（%）"
  ),
  数值 = c(
    n_cohorts,
    n_diseases,
    n_samples,
    lefse_total,
    round(lefse_mean, 2),
    maaslin_total,
    round(maaslin_mean, 2),
    ancom_total,
    round(ancom_mean, 2),
    round(jaccard_ml  * 100, 2),
    round(jaccard_al  * 100, 2),
    round(jaccard_am  * 100, 2),
    round(jaccard_all * 100, 2),
    round(m_ml$overlap_atk * 100, 2),
    round(m_al$overlap_atk * 100, 2),
    round(m_am$overlap_atk * 100, 2),
    round(top50_jaccard_all * 100, 2)
  )
)

print(summary_report)

write.csv(summary_report, "Analysis_summary_report.csv",
          fileEncoding = "GBK",
          row.names = FALSE)

cat("\n汇总报告已保存至: Analysis_summary_report.csv\n")


