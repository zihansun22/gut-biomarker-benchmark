
# 包含有RF的输入函数，统一读入，但RF自己内部做过滤
# 主要输入：计数矩阵 + 元数据
generate_method_inputs <- function(scenario_full, rep_id) {
  counts_df   <- read_tsv(file.path("mock_benchmark_data",
                                    paste0("mock_counts_",   scenario_full, "_rep", rep_id, ".tsv")))
  relab_df    <- read_tsv(file.path("mock_benchmark_data",
                                    paste0("mock_relab_",    scenario_full, "_rep", rep_id, ".tsv")))
  metadata_df <- read_tsv(file.path("mock_benchmark_data",
                                    paste0("mock_metadata_", scenario_full, "_rep", rep_id, ".tsv")))
  
  counts_matrix <- as.matrix(counts_df[, -1]); rownames(counts_matrix) <- counts_df$taxa
  relab_matrix  <- as.matrix(relab_df[, -1]);  rownames(relab_matrix)  <- relab_df$taxa
  
  # ── 统一过滤（用于 LEfSe / ANCOM-BC2 / DESeq2 / MaAsLin2）──
  prevalence <- rowMeans(counts_matrix > 0)
  mean_relab <- rowMeans(relab_matrix)
  keep_taxa  <- prevalence >= 0.1 & mean_relab >= 1e-4
  
  list(
    lefse = list(
      abundance = relab_matrix[keep_taxa, ],
      groups    = metadata_df$phenotype
    ),
    ancombc2 = list(
      counts   = counts_matrix[keep_taxa, ],
      metadata = metadata_df
    ),
    deseq2 = list(
      counts = counts_matrix[keep_taxa, ],
      groups = metadata_df
    ),
    maaslin2 = list(
      abundance = counts_matrix[keep_taxa, ],
      metadata  = metadata_df
    ),
    
    # RF 只做格式化，不过滤、不变换
    # 过滤和 CLR 在训练折内部完成，避免测试集信息泄露
    rf = list(
      counts   = counts_matrix,   # 原始计数，未过滤
      relab    = relab_matrix,    # 原始相对丰度，未过滤
      metadata = metadata_df
    )
  )
}


# 五种方法的结果汇总，macro方式合并
setwd("D:/大四/毕设/data")

library(tidyverse)
library(stringr)

# ----------------------------
# 0) Config
# ----------------------------
RESULT_DIR <- "benchmark/results"
OUT_DIR    <- "benchmark/evaluation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

scenarios <- c(
  "A_no_confound",
  "B_batch_confound",
  "C_sparse_small",
  "D_large_effect"
)

all_comparisons <- tibble::tibble(
  comparison = c("Disease_mild_vs_Disease_severe",
                 "Healthy_vs_Disease_mild",
                 "Healthy_vs_Disease_severe"),
  group_a    = c("Disease_mild", "Healthy", "Healthy"),
  group_b    = c("Disease_severe", "Disease_mild", "Disease_severe")
)

# ----------------------------
# 读入 + 全量修正（method命名 + comparison方向 + 所有关联字段）
# ----------------------------
naive_all <- readr::read_csv(file.path(RESULT_DIR, "naive_summary_all.csv"),
                             show_col_types = FALSE)

# 需要翻转的 comparison（与 all_comparisons 标准方向相反）
FLIP_COMPARISONS <- c(
  "Disease_mild_vs_Healthy",        # ancombc2 / maaslin2 / deseq2
  "Disease_severe_vs_Healthy",      # ancombc2 / maaslin2 / deseq2
  "Disease_severe_vs_Disease_mild"  # deseq2 only
)

# comparison 字符串翻转映射
FLIP_MAP <- c(
  "Disease_mild_vs_Healthy"        = "Healthy_vs_Disease_mild",
  "Disease_severe_vs_Healthy"      = "Healthy_vs_Disease_severe",
  "Disease_severe_vs_Disease_mild" = "Disease_mild_vs_Disease_severe"
)

naive_combined <- naive_all %>%
  mutate(
    # ── 基础类型修正 ──────────────────────────────────────
    significant = as.logical(significant),
    replicate   = as.integer(replicate),
    effect      = as.numeric(effect),
    scenario    = str_trim(as.character(scenario)),
    comparison  = str_trim(as.character(comparison)),
    group_a     = str_trim(as.character(group_a)),
    group_b     = str_trim(as.character(group_b)),
    ref_group   = str_trim(as.character(ref_group)),
    
    # ── taxa_short 清理 ───────────────────────────────────
    taxa_short = str_trim(as.character(taxa_short)),
    taxa_short = na_if(taxa_short, ""),
    taxa_short = na_if(taxa_short, "__"),
    
    # ── ① method + setting → 唯一方法名 ──────────────────
    method = case_when(
      method == "deseq2" & str_trim(setting) == "default"  ~ "deseq2-default",
      method == "deseq2" & str_trim(setting) == "no_cooks" ~ "deseq2-nocook",
      method == "deseq2" & is.na(setting)                  ~ "deseq2-default",
      TRUE ~ as.character(method)
    ),
    method = str_trim(method),
    
    # ── ② 标记需要翻转的行 ────────────────────────────────
    need_flip = comparison %in% FLIP_COMPARISONS,
    
    # ── ③ comparison 字符串翻转 ───────────────────────────
    comparison = ifelse(
      need_flip,
      recode(comparison, !!!FLIP_MAP),
      comparison
    ),
    
    # ── ④ effect（log2FC / coef）取反 ────────────────────
    effect = ifelse(need_flip, -effect, effect),
    
    # ── ⑤ group_a / group_b 互换 ─────────────────────────
    .tmp_a  = group_a,
    group_a = ifelse(need_flip, group_b,  group_a),
    group_b = ifelse(need_flip, .tmp_a,   group_b),
    
    # ── ⑥ ref_group 更新 ──────────────────────────────────
    # 翻转后：ref_group = 新的 group_b（即原来的 group_a）
    ref_group = ifelse(need_flip, group_b, ref_group)
    # 注意：此时 group_b 已经是翻转后的值（即原 group_a）
  ) %>%
  select(-need_flip, -.tmp_a) %>%
  filter(!is.na(taxa_short))

# ============================================================
# 验证
# ============================================================
cat("=== ① Methods（应有5个）===\n")
print(count(naive_combined, method))

cat("\n=== ② Comparisons（应只有3个）===\n")
print(count(naive_combined, comparison))

cat("\n=== ③ group_a / group_b / ref_group 一致性（所有方法应对齐）===\n")
naive_combined %>%
  distinct(method, comparison, group_a, group_b, ref_group) %>%
  arrange(comparison, method) %>%
  print(n = 50)

cat("\n=== ④ 与 all_comparisons 标准对齐检查 ===\n")
align_check <- naive_combined %>%
  distinct(method, comparison, group_a, group_b) %>%
  left_join(all_comparisons, by = "comparison", suffix = c("_data", "_ref")) %>%
  mutate(
    group_a_ok = group_a_data == group_a_ref,
    group_b_ok = group_b_data == group_b_ref,
    all_ok     = group_a_ok & group_b_ok
  )
print(arrange(align_check, comparison, method), n = 50)
cat("全部对齐？", all(align_check$all_ok, na.rm = TRUE), "\n")

cat("\n=== ⑤ effect 符号分布（显著taxa，各方法应方向一致）===\n")
naive_combined %>%
  filter(significant) %>%
  mutate(effect_sign = ifelse(effect > 0, "positive", "negative")) %>%
  count(method, comparison, effect_sign) %>%
  arrange(comparison, method) %>%
  print(n = 60)

# Sanity check
bad_comp <- setdiff(unique(naive_combined$comparison), all_comparisons$comparison)
if (length(bad_comp) > 0) {
  warning("仍有未预期的 comparison: ", paste(bad_comp, collapse = ", "))
} else {
  cat("\n✓ 所有 comparison 均合法\n")
}

# 修改前：naive_all (benchmark/results/naive_summary_all.csv)
# 修改后：naive_combined
# 去除setting（method列已有对deseq-setting的备注）
naive_combined <- naive_combined %>%
  select(-setting)
# 保存在/benchmark/evaluation/
write.csv(naive_combined, 
          file.path("benchmark/evaluation/naive_summary_all_standardized.csv"))


# timing表的合并
library(dplyr)
library(purrr)

# 需要先读入各个方法的timing_raw 和 timing_summary文件。
# 1) 把所有 timing_raw 合并（可选，但建议存档）
timing_raw_all <- bind_rows(
  lefse   %>% mutate(method = as.character(method)),
  maaslin %>% mutate(method = as.character(method)),
  deseq2_default   %>% mutate(method = as.character(method)),
  deseq2_nocooks   %>% mutate(method = as.character(method)),
  ancombc %>% mutate(method = as.character(method)),
  rf   %>% mutate(method = as.character(method))
) %>%
  mutate(
    scenario = as.character(scenario),
    rep = as.integer(rep)
  )

# 2) 把所有 timing_summary 合并（Table 1 的 join 输入）
timing_sum_all <- bind_rows(
  lefse_sum,
  maaslin_sum,
  deseq2_default_sum,
  deseq2_nocooks_sum,
  ancombc_sum,
  rf_sum
) %>%
  mutate(
    method   = as.character(method),
    scenario = as.character(scenario)
  )


# Jaccard表的合并
jaccard_all <- bind_rows(
  ancombc_scn,
  rf_jaccard_summary,
  lefse_jaccard_summary,
  maaslin2_jaccard_summary,
  deseq2_jaccard_summary_default,
  deseq2_jaccard_summary_no_cooks
)

# 有列名不统一的情况
# 直接在jaccard_all上进行修改
rows_to_modify <- 9:nrow(jaccard_all)
# 直接赋值
jaccard_all$n_pairs[rows_to_modify] <- 3
jaccard_all$top_k_total[rows_to_modify] <- 150
jaccard_all$jaccard_mean_across_pairs[rows_to_modify] <- jaccard_all$jaccard_mean_across_reps[rows_to_modify]
jaccard_all$jaccard_median_across_pairs[rows_to_modify] <- jaccard_all$jaccard_median_across_reps[rows_to_modify]
# 重命名：现在的jaccard已经既是across_pairs，又across_reps了
jaccard_all <- jaccard_all %>% rename(
  jaccard_mean = jaccard_mean_across_pairs,
  jaccard_median = jaccard_median_across_pairs
)
# 保存
write.csv(jaccard_all,"benchmark/evaluation/naive_jaccard_all.csv")




## 主表1
# ============================================================
# MAIN TABLE 1 (scenario × method)
# Clean version (RF uses precomputed micro metrics)
# Strict MACRO batch_error
# RF batch_error via Top-K feature importance (from fold files)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(purrr)
})

OUT_DIR <- "benchmark/evaluation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCENARIO_B  <- "B_batch_confound"
RF_RUN_ROOT <- "benchmark/rf/rf_results_naive"
K_RF        <- 50L   # 与 HitRate_at_50 保持一致

# ------------------------------------------------------------
# 0) 统一 method 名称函数
# ------------------------------------------------------------

standardize_method <- function(x){
  x %>%
    str_to_lower() %>%
    str_replace("_naive$", "") %>%
    str_replace("ancom-bc2", "ancombc2") %>%
    str_replace("rf_naive", "rf")
}

# ------------------------------------------------------------
# 1) 读取汇总文件
# ------------------------------------------------------------

diff_all <- read_csv(
  "benchmark/evaluation/naive_summary_all_standardized.csv",
  show_col_types = FALSE
) %>%
  mutate(
    method      = standardize_method(method),
    scenario    = str_trim(as.character(scenario)),
    comparison  = str_trim(as.character(comparison)),
    replicate   = as.integer(replicate),
    taxa_short  = str_trim(as.character(taxa_short)),
    significant = as.logical(significant)
  ) %>%
  filter(!is.na(taxa_short))

rf_metrics <- read.csv(
  "benchmark/results/rf_naive_metrics_summary.csv"
) %>%
  mutate(
    method   = "rf",
    scenario = str_trim(as.character(scenario))
  ) %>%
  select(method, scenario,
         Hit_at_50_mean, Recall_at_50_mean, AUC_mean) %>%
  rename(
    HitRate_at_50 = Hit_at_50_mean,
    Recall_at_50  = Recall_at_50_mean,
    AUC           = AUC_mean
  )

jaccard_all <- read_csv(
  "benchmark/evaluation/naive_jaccard_all.csv",
  show_col_types = FALSE
) %>%
  mutate(
    method   = standardize_method(method),
    scenario = str_trim(as.character(scenario))
  ) %>%
  select(method, scenario, jaccard_mean)

timing_all <- read_csv(
  "benchmark/evaluation/timing_sum_all.csv",
  show_col_types = FALSE
) %>%
  mutate(
    method   = standardize_method(method),
    scenario = str_trim(as.character(scenario))
  ) %>%
  select(method, scenario, wall_median_sec)

# ------------------------------------------------------------
# 2) 构建 truth（仅 disease marker）
# ------------------------------------------------------------

all_comparisons <- tibble(
  comparison = c("Disease_mild_vs_Disease_severe",
                 "Healthy_vs_Disease_mild",
                 "Healthy_vs_Disease_severe"),
  group_a    = c("Disease_mild", "Healthy", "Healthy"),
  group_b    = c("Disease_severe", "Disease_mild", "Disease_severe")
)

truth <- read_tsv(
  "mock_benchmark_data/all_truth_tables.tsv",
  show_col_types = FALSE
)

truth_long <- truth %>%
  filter(marker_type != "batch_effect") %>%
  transmute(
    scenario    = str_trim(as.character(scenario)),
    replicate   = as.integer(replicate),
    taxa,
    groups_list = str_split(present_in_groups, ",") %>%
      map(~ str_trim(.x))
  ) %>%
  crossing(all_comparisons) %>%
  mutate(
    in_a       = pmap_lgl(list(group_a, groups_list), ~ ..1 %in% ..2),
    in_b       = pmap_lgl(list(group_b, groups_list), ~ ..1 %in% ..2),
    is_tp      = xor(in_a, in_b),
    taxa_short = str_extract(taxa, "[^;]+$") %>% str_trim()
  ) %>%
  filter(is_tp) %>%
  distinct(scenario, replicate, comparison, taxa_short)

truth_sets <- truth_long %>%
  group_by(scenario, replicate, comparison) %>%
  summarise(
    truth_set = list(unique(taxa_short)),
    .groups   = "drop"
  )

# ------------------------------------------------------------
# 3) 四种差异方法 MACRO 指标
#    先对每个 (replicate × comparison) 单独计算 precision/recall/F1/FDR，
#    再对 (method × scenario) 取算术均值 → Macro 聚合
# ------------------------------------------------------------

pred_sets <- diff_all %>%
  filter(significant) %>%
  group_by(method, scenario, replicate, comparison) %>%
  summarise(pred_set = list(unique(taxa_short)), .groups = "drop")

# 与 truth 连接：没有显著预测的 (replicate × comparison) 通过 right_join 补全
# 保留所有 truth 存在的组合，pred_set 为 NULL 表示该折无预测
metrics_diff <- truth_sets %>%
  left_join(
    pred_sets,
    by = c("scenario", "replicate", "comparison")
  ) %>%
  # 对于 pred_set 为 NULL 的行（无显著结果），视 pred_set 为空集
  mutate(
    pred_set = map(pred_set, ~ if (is.null(.x)) character(0) else .x),
    TP = map2_int(truth_set, pred_set, ~ length(intersect(.x, .y))),
    FP = map2_int(truth_set, pred_set, ~ length(setdiff(.y, .x))),
    FN = map2_int(truth_set, pred_set, ~ length(setdiff(.x, .y))),
    # 单折指标：分母为 0 时赋 NA（后续 mean 自动 na.rm 处理）
    precision_i = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
    recall_i    = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
    F1_i        = if_else(
      !is.na(precision_i) & !is.na(recall_i) & (precision_i + recall_i) > 0,
      2 * precision_i * recall_i / (precision_i + recall_i),
      NA_real_
    ),
    FDR_i       = if_else(TP + FP > 0, FP / (TP + FP), NA_real_)
  )

# Macro 聚合：对每个 (method × scenario) 取均值
macro_diff <- metrics_diff %>%
  group_by(method, scenario) %>%
  summarise(
    precision = mean(precision_i, na.rm = TRUE),
    recall    = mean(recall_i,    na.rm = TRUE),
    F1        = mean(F1_i,        na.rm = TRUE),
    FDR       = mean(FDR_i,       na.rm = TRUE),
    .groups   = "drop"
  )

# ------------------------------------------------------------
# 4) batch_error（MACRO 聚合）
#    差异方法：每 replicate 算 batch_sig/n_sig，再取均值
#    RF：     每 fold 算 batch_in_TopK/K，再取均值
# ------------------------------------------------------------

# ── batch taxa 列表（场景B）────────────────────────────────
batch_taxa <- truth %>%
  filter(marker_type == "batch_effect",
         scenario    == SCENARIO_B) %>%
  mutate(taxa_short = str_extract(taxa, "[^;]+$")) %>%
  pull(taxa_short) %>%
  unique()

cat("batch_taxa数量:", length(batch_taxa), "\n")

# ── 差异方法：每 (replicate × comparison) 算比率，再 Macro 均值 ──
batch_rep <- diff_all %>%
  filter(significant,
         scenario == SCENARIO_B) %>%
  group_by(method, scenario, replicate, comparison) %>%
  summarise(
    n_sig_fold   = n_distinct(taxa_short),
    n_batch_fold = sum(taxa_short %in% batch_taxa),
    .groups      = "drop"
  ) %>%
  mutate(
    # 单折 batch_error；n_sig_fold == 0 时视为 NA
    batch_error_i = if_else(n_sig_fold > 0,
                            n_batch_fold / n_sig_fold,
                            NA_real_)
  )

batch_rate_diff <- batch_rep %>%
  group_by(method, scenario) %>%
  summarise(
    batch_error = mean(batch_error_i, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  select(method, scenario, batch_error)

# ── RF：读取各fold importance文件，每 fold 算比率，再 Macro 均值 ──

importance_files <- list.files(
  path       = RF_RUN_ROOT,
  pattern    = "_importance\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

cat("找到importance文件数:", length(importance_files), "\n")

if (length(importance_files) == 0) {
  
  warning("未找到RF importance文件，RF batch_error 将设为 NA")
  rf_batch_rate <- tibble(
    method      = "rf",
    scenario    = SCENARIO_B,
    batch_error = NA_real_
  )
  
} else {
  
  rf_importance_all <- importance_files %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE)) %>%
    mutate(
      scenario   = str_trim(as.character(scenario)),
      taxa_short = str_trim(as.character(taxa_short)),
      test_rep   = as.integer(test_rep)
    )
  
  cat("RF importance总行数:", nrow(rf_importance_all), "\n")
  cat("RF importance列名:", paste(colnames(rf_importance_all), collapse = ", "), "\n")
  
  # 每个 (test_rep × comparison) fold 取 Top-K，计算该 fold 的 batch_error_i
  rf_batch_rep <- rf_importance_all %>%
    filter(scenario == SCENARIO_B,
           rank     <= K_RF) %>%
    group_by(scenario, test_rep, comparison) %>%
    summarise(
      k_actual      = n(),
      n_batch_fold  = sum(taxa_short %in% batch_taxa),
      .groups       = "drop"
    ) %>%
    mutate(
      # 单折 batch_error；k_actual == 0 时视为 NA（理论上不会出现）
      batch_error_i = if_else(k_actual > 0,
                              n_batch_fold / k_actual,
                              NA_real_)
    )
  
  cat("\nRF batch_rep 诊断（各fold）:\n")
  print(rf_batch_rep)
  
  # Macro 聚合：对所有 fold 取均值
  rf_batch_rate <- rf_batch_rep %>%
    summarise(
      batch_error = mean(batch_error_i, na.rm = TRUE)
    ) %>%
    mutate(
      method   = "rf",
      scenario = SCENARIO_B
    ) %>%
    select(method, scenario, batch_error)
  
  cat("\nRF batch_error (Macro Top-K):", round(rf_batch_rate$batch_error, 3), "\n")
}

# ── 合并差异方法与RF的batch_error ────────────────────────
batch_rate_all <- bind_rows(batch_rate_diff, rf_batch_rate)

# ------------------------------------------------------------
# 5) 合并所有指标 → Table1
# ------------------------------------------------------------

table1 <- full_join(macro_diff,
                    rf_metrics,
                    by = c("method", "scenario")) %>%
  left_join(jaccard_all,
            by = c("method", "scenario")) %>%
  left_join(timing_all,
            by = c("method", "scenario")) %>%
  left_join(batch_rate_all,
            by = c("method", "scenario")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(scenario)

# ------------------------------------------------------------
# 6) 输出
# ------------------------------------------------------------

write_csv(table1,
          file.path(OUT_DIR, "MAIN_Table1_scenario_method.csv"))

cat("\n================ FINAL TABLE 1 (MACRO) ================\n")
print(table1, n = 100, width = 150)

cat("\n================ batch_error 汇总 (MACRO) ================\n")
print(batch_rate_all)



# table2
# ============================================================
# MAIN TABLE 2 (Scenario B × naive/adjusted 对比)
# batch 混杂场景下的方法鲁棒性
# 四种差异方法 naive/adjusted 的 FDR/F1
# RF naive/adjusted 的 AUC/Hit@K
# 五种方法的 batch误选率 + Jaccard + 时间
# RF batch_error via Top-K importance (from fold files)
# 全部采用 MACRO 聚合
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(purrr)
})

OUT_DIR <- "benchmark/evaluation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCENARIO_B       <- "B_batch_confound"
RF_NAIVE_ROOT    <- "benchmark/rf/rf_results_naive"
RF_ADJUSTED_ROOT <- "benchmark/rf/rf_results_adjusted"
K_RF             <- 50L   # 与 HitRate_at_50 保持一致

# ------------------------------------------------------------
# 0) 统一 method 名称函数
# ------------------------------------------------------------

standardize_method <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace("ancom-bc2", "ancombc2") %>%
    str_replace("ancombc2_naive",    "ancombc2") %>%
    str_replace("ancombc2_adjusted", "ancombc2") %>%
    str_replace("maaslin2_naive",    "maaslin2") %>%
    str_replace("maaslin2_adjusted", "maaslin2") %>%
    str_replace("deseq2-default_naive",    "deseq2-default") %>%
    str_replace("deseq2-default_adjusted", "deseq2-default") %>%
    str_replace("deseq2-nocook_naive",    "deseq2-nocook") %>%
    str_replace("deseq2-nocook_adjusted", "deseq2-nocook") %>%
    str_replace("lefse_naive",    "lefse") %>%
    str_replace("lefse_adjusted", "lefse")
}

# ------------------------------------------------------------
# 1) 读取数据：naive + adjusted 分别读入，打标签后合并
# ------------------------------------------------------------

# ── 差异分析结果 ──────────────────────────────────────────
diff_naive <- read_csv(
  "benchmark/evaluation/naive_summary_all_standardized.csv",
  show_col_types = FALSE
) %>%
  mutate(mode = "naive")

diff_adjusted <- read_csv(
  "benchmark/evaluation/adjusted_summary_all_standardized.csv",
  show_col_types = FALSE
) %>%
  mutate(mode = "adjusted")

diff_all <- bind_rows(diff_naive, diff_adjusted) %>%
  mutate(
    method      = standardize_method(method),
    mode        = str_trim(as.character(mode)),
    scenario    = str_trim(as.character(scenario)),
    comparison  = str_trim(as.character(comparison)),
    replicate   = as.integer(replicate),
    taxa_short  = str_trim(as.character(taxa_short)),
    significant = as.logical(significant)
  ) %>%
  filter(
    !is.na(taxa_short),
    scenario == SCENARIO_B
  )

# ── RF 指标：naive + adjusted ─────────────────────────────
rf_naive_metrics <- read.csv(
  "benchmark/results/rf_naive_metrics_summary.csv"
) %>%
  mutate(mode = "naive", method = "rf",
         scenario = str_trim(as.character(scenario))) %>%
  select(method, mode, scenario,
         Hit_at_50_mean, Recall_at_50_mean, AUC_mean) %>%
  rename(HitRate_at_50 = Hit_at_50_mean,
         Recall_at_50  = Recall_at_50_mean,
         AUC           = AUC_mean)

rf_adjusted_metrics <- read.csv(
  "benchmark/results/rf_adjusted_metrics_summary.csv"
) %>%
  mutate(mode = "adjusted", method = "rf",
         scenario = str_trim(as.character(scenario))) %>%
  select(method, mode, scenario,
         Hit_at_50_mean, Recall_at_50_mean, AUC_mean) %>%
  rename(HitRate_at_50 = Hit_at_50_mean,
         Recall_at_50  = Recall_at_50_mean,
         AUC           = AUC_mean)

rf_metrics <- bind_rows(rf_naive_metrics, rf_adjusted_metrics) %>%
  filter(scenario == SCENARIO_B)

# ── Jaccard：naive + adjusted ─────────────────────────────
jaccard_naive <- read_csv(
  "benchmark/evaluation/naive_jaccard_all.csv",
  show_col_types = FALSE
) %>%
  mutate(mode     = "naive",
         method   = standardize_method(method),
         scenario = str_trim(as.character(scenario)))

jaccard_adjusted <- read_csv(
  "benchmark/evaluation/adjusted_jaccard_all.csv",
  show_col_types = FALSE
) %>%
  mutate(mode     = "adjusted",
         method   = standardize_method(method),
         scenario = str_trim(as.character(scenario)))

jaccard_all <- bind_rows(jaccard_naive, jaccard_adjusted) %>%
  filter(scenario == SCENARIO_B) %>%
  select(method, mode, scenario, jaccard_mean)

# ── Timing：naive + adjusted 在同一文件中 ─────────────────
timing_all <- read_csv(
  "benchmark/evaluation/timing_sum_all.csv",
  show_col_types = FALSE
) %>%
  mutate(
    mode = case_when(
      str_detect(method, "_adjusted") ~ "adjusted",
      str_detect(method, "_naive")    ~ "naive",
      TRUE                            ~ "naive"
    ),
    method   = standardize_method(method),
    scenario = str_trim(as.character(scenario))
  ) %>%
  filter(scenario == SCENARIO_B) %>%
  select(method, mode, scenario, wall_median_sec)

# ------------------------------------------------------------
# 2) 构建 truth（仅 disease marker，场景B）
# ------------------------------------------------------------

all_comparisons <- tibble(
  comparison = c("Disease_mild_vs_Disease_severe",
                 "Healthy_vs_Disease_mild",
                 "Healthy_vs_Disease_severe"),
  group_a    = c("Disease_mild", "Healthy", "Healthy"),
  group_b    = c("Disease_severe", "Disease_mild", "Disease_severe")
)

truth <- read_tsv(
  "mock_benchmark_data/all_truth_tables.tsv",
  show_col_types = FALSE
)

truth_long <- truth %>%
  filter(
    marker_type != "batch_effect",
    scenario    == SCENARIO_B
  ) %>%
  transmute(
    scenario    = str_trim(as.character(scenario)),
    replicate   = as.integer(replicate),
    taxa,
    groups_list = str_split(present_in_groups, ",") %>%
      map(~ str_trim(.x))
  ) %>%
  crossing(all_comparisons) %>%
  mutate(
    in_a       = pmap_lgl(list(group_a, groups_list), ~ ..1 %in% ..2),
    in_b       = pmap_lgl(list(group_b, groups_list), ~ ..1 %in% ..2),
    is_tp      = xor(in_a, in_b),
    taxa_short = str_extract(taxa, "[^;]+$") %>% str_trim()
  ) %>%
  filter(is_tp) %>%
  distinct(scenario, replicate, comparison, taxa_short)

truth_sets <- truth_long %>%
  group_by(scenario, replicate, comparison) %>%
  summarise(
    truth_set = list(unique(taxa_short)),
    .groups   = "drop"
  )

# ------------------------------------------------------------
# 3) 差异方法 MACRO 指标（FDR / F1），按 mode 分开
#    先对每个 (replicate × comparison) 单独计算各指标，
#    再对 (method × mode × scenario) 取算术均值
# ------------------------------------------------------------

pred_sets <- diff_all %>%
  filter(significant) %>%
  group_by(method, mode, scenario, replicate, comparison) %>%
  summarise(pred_set = list(unique(taxa_short)), .groups = "drop")

# truth_sets 为主表，left_join pred_sets
# 确保"无显著预测"的折也参与 Macro 均值计算
metrics_diff <- truth_sets %>%
  # 为每个 (method × mode) 组合都复制一份 truth，再 join pred
  crossing(
    diff_all %>%
      distinct(method, mode)
  ) %>%
  left_join(
    pred_sets,
    by = c("method", "mode", "scenario", "replicate", "comparison")
  ) %>%
  mutate(
    # pred_set 为 NULL 时视为空集（该折无显著结果）
    pred_set    = map(pred_set, ~ if (is.null(.x)) character(0) else .x),
    TP          = map2_int(truth_set, pred_set, ~ length(intersect(.x, .y))),
    FP          = map2_int(truth_set, pred_set, ~ length(setdiff(.y, .x))),
    FN          = map2_int(truth_set, pred_set, ~ length(setdiff(.x, .y))),
    # 单折指标；分母为 0 时赋 NA
    precision_i = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
    recall_i    = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
    F1_i        = if_else(
      !is.na(precision_i) & !is.na(recall_i) & (precision_i + recall_i) > 0,
      2 * precision_i * recall_i / (precision_i + recall_i),
      NA_real_
    ),
    FDR_i       = if_else(TP + FP > 0, FP / (TP + FP), NA_real_)
  )

# Macro 聚合：对每个 (method × mode × scenario) 取均值
macro_diff <- metrics_diff %>%
  group_by(method, mode, scenario) %>%
  summarise(
    precision = mean(precision_i, na.rm = TRUE),
    recall    = mean(recall_i,    na.rm = TRUE),
    F1        = mean(F1_i,        na.rm = TRUE),
    FDR       = mean(FDR_i,       na.rm = TRUE),
    .groups   = "drop"
  )

# ------------------------------------------------------------
# 4) Batch 误选率（MACRO 聚合）
#    差异方法：每 (replicate × comparison) 折算 batch_error_i，再均值
#    RF：     每 fold (test_rep × comparison) 算 batch_error_i，再均值
# ------------------------------------------------------------

# ── batch taxa 列表（场景B）────────────────────────────────
batch_taxa <- truth %>%
  filter(
    marker_type == "batch_effect",
    scenario    == SCENARIO_B
  ) %>%
  mutate(taxa_short = str_extract(taxa, "[^;]+$")) %>%
  pull(taxa_short) %>%
  unique()

cat("batch_taxa 数量:", length(batch_taxa), "\n")

# ── 差异方法：每 (replicate × comparison) 折算比率，再 Macro 均值 ──
batch_rep_diff <- diff_all %>%
  filter(significant) %>%
  group_by(method, mode, scenario, replicate, comparison) %>%
  summarise(
    n_sig_fold   = n_distinct(taxa_short),
    n_batch_fold = sum(taxa_short %in% batch_taxa),
    .groups      = "drop"
  ) %>%
  mutate(
    # 单折 batch_error；n_sig_fold == 0 时赋 NA
    batch_error_i = if_else(n_sig_fold > 0,
                            n_batch_fold / n_sig_fold,
                            NA_real_)
  )

batch_rate_diff <- batch_rep_diff %>%
  group_by(method, mode, scenario) %>%
  summarise(
    batch_error       = mean(batch_error_i, na.rm = TRUE),
    # 绝对数量仍做简单求和，便于参考
    absolute_batch_FP = sum(n_batch_fold),
    .groups           = "drop"
  ) %>%
  select(method, mode, scenario, batch_error, absolute_batch_FP)

# ── RF batch_error：读取 fold importance 文件的通用函数（MACRO）────
compute_rf_batch_error <- function(run_root, mode_label, scenario_name,
                                   batch_taxa_vec, k = K_RF) {
  
  imp_files <- list.files(
    path       = run_root,
    pattern    = "_importance\\.csv$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  cat(sprintf("[RF %s] 找到 importance 文件数: %d\n", mode_label, length(imp_files)))
  
  if (length(imp_files) == 0) {
    warning(sprintf("未找到 RF %s importance 文件，batch_error 设为 NA", mode_label))
    return(
      tibble(
        method            = "rf",
        mode              = mode_label,
        scenario          = scenario_name,
        batch_error       = NA_real_,
        absolute_batch_FP = NA_real_
      )
    )
  }
  
  imp_all <- imp_files %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE)) %>%
    mutate(
      scenario   = str_trim(as.character(scenario)),
      taxa_short = str_trim(as.character(taxa_short)),
      test_rep   = as.integer(test_rep)
    )
  
  cat(sprintf("[RF %s] importance 总行数: %d\n", mode_label, nrow(imp_all)))
  
  imp_scenarioB <- imp_all %>%
    filter(scenario == scenario_name,
           rank     <= k)
  
  if (nrow(imp_scenarioB) == 0) {
    warning(sprintf("RF %s 在场景 %s 无数据，batch_error 设为 NA",
                    mode_label, scenario_name))
    return(
      tibble(
        method            = "rf",
        mode              = mode_label,
        scenario          = scenario_name,
        batch_error       = NA_real_,
        absolute_batch_FP = NA_real_
      )
    )
  }
  
  # 每个 fold（test_rep × comparison）单独计算 batch_error_i
  batch_rep_rf <- imp_scenarioB %>%
    group_by(scenario, test_rep, comparison) %>%
    summarise(
      k_actual      = n(),
      n_batch_fold  = sum(taxa_short %in% batch_taxa_vec),
      .groups       = "drop"
    ) %>%
    mutate(
      # 单折 batch_error；k_actual == 0 时赋 NA（理论上不会出现）
      batch_error_i = if_else(k_actual > 0,
                              n_batch_fold / k_actual,
                              NA_real_)
    )
  
  cat(sprintf("[RF %s] fold 级别诊断:\n", mode_label))
  print(batch_rep_rf)
  
  # Macro 聚合：对所有 fold 取均值
  result <- batch_rep_rf %>%
    summarise(
      batch_error       = mean(batch_error_i, na.rm = TRUE),
      absolute_batch_FP = sum(n_batch_fold)
    ) %>%
    mutate(
      method   = "rf",
      mode     = mode_label,
      scenario = scenario_name
    ) %>%
    select(method, mode, scenario, batch_error, absolute_batch_FP)
  
  cat(sprintf("[RF %s] batch_error Macro (Top-%d): %.3f\n",
              mode_label, k, result$batch_error))
  
  result
}

# ── 分别计算 naive / adjusted 的 RF batch_error ───────────
rf_batch_naive <- compute_rf_batch_error(
  run_root       = RF_NAIVE_ROOT,
  mode_label     = "naive",
  scenario_name  = SCENARIO_B,
  batch_taxa_vec = batch_taxa,
  k              = K_RF
)

rf_batch_adjusted <- compute_rf_batch_error(
  run_root       = RF_ADJUSTED_ROOT,
  mode_label     = "adjusted",
  scenario_name  = SCENARIO_B,
  batch_taxa_vec = batch_taxa,
  k              = K_RF
)

# ── 合并差异方法与 RF 的 batch_error ─────────────────────
batch_rate_all <- bind_rows(
  batch_rate_diff,
  rf_batch_naive,
  rf_batch_adjusted
)

cat("\n====== batch_error 汇总 (MACRO) ======\n")
print(batch_rate_all)

# ------------------------------------------------------------
# 5) 合并所有指标 → Table2 (long)
# ------------------------------------------------------------

# 差异方法部分
diff_part <- macro_diff %>%
  left_join(batch_rate_all,
            by = c("method", "mode", "scenario")) %>%
  left_join(jaccard_all,
            by = c("method", "mode", "scenario")) %>%
  left_join(timing_all,
            by = c("method", "mode", "scenario")) %>%
  mutate(
    AUC           = NA_real_,
    HitRate_at_50 = NA_real_
  ) %>%
  select(method, mode, scenario,
         FDR, F1,
         AUC, HitRate_at_50,
         batch_error, absolute_batch_FP,
         jaccard_mean,
         wall_median_sec)

# RF 部分
rf_part <- rf_metrics %>%
  left_join(batch_rate_all %>% filter(method == "rf"),
            by = c("method", "mode", "scenario")) %>%
  left_join(jaccard_all %>% filter(method == "rf"),
            by = c("method", "mode", "scenario")) %>%
  left_join(timing_all %>% filter(method == "rf"),
            by = c("method", "mode", "scenario")) %>%
  mutate(
    FDR = NA_real_,
    F1  = NA_real_
  ) %>%
  select(method, mode, scenario,
         FDR, F1,
         AUC, HitRate_at_50,
         batch_error, absolute_batch_FP,
         jaccard_mean,
         wall_median_sec)

table2_long <- bind_rows(diff_part, rf_part) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(method, mode)

# ------------------------------------------------------------
# 6) 宽表：naive vs adjusted 并排展示
# ------------------------------------------------------------

table2_wide <- table2_long %>%
  select(-scenario) %>%
  pivot_wider(
    names_from  = mode,
    values_from = c(FDR, F1, AUC, HitRate_at_50,
                    batch_error, absolute_batch_FP,
                    jaccard_mean, wall_median_sec),
    names_glue  = "{.value}_{mode}"
  ) %>%
  select(
    method,
    FDR_naive,               FDR_adjusted,
    F1_naive,                F1_adjusted,
    AUC_naive,               AUC_adjusted,
    HitRate_at_50_naive,     HitRate_at_50_adjusted,
    batch_error_naive,       batch_error_adjusted,
    absolute_batch_FP_naive, absolute_batch_FP_adjusted,
    jaccard_mean_naive,      jaccard_mean_adjusted,
    wall_median_sec_naive,   wall_median_sec_adjusted
  )

# ------------------------------------------------------------
# 7) 输出
# ------------------------------------------------------------

write_csv(table2_long,
          file.path(OUT_DIR, "MAIN_Table2_scenarioB_long.csv"))
write_csv(table2_wide,
          file.path(OUT_DIR, "MAIN_Table2_scenarioB_wide.csv"))

cat("\n================ TABLE 2 (long, MACRO) ================\n")
print(table2_long, n = 50, width = 150)

cat("\n================ TABLE 2 (wide, naive vs adjusted, MACRO) ================\n")
print(table2_wide, n = 20, width = 200)



