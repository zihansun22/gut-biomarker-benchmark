# ══════════════════════════════════════════════════════════════
# Thread lock（建议保留；ranger 也会显式 num.threads=1）
# ══════════════════════════════════════════════════════════════
NTHREAD <- 1L

Sys.setenv(
  OMP_NUM_THREADS        = NTHREAD,
  OPENBLAS_NUM_THREADS   = NTHREAD,
  MKL_NUM_THREADS        = NTHREAD,
  VECLIB_MAXIMUM_THREADS = NTHREAD,
  NUMEXPR_NUM_THREADS    = NTHREAD
)
options(mc.cores = NTHREAD)

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(NTHREAD)
}

cat("Thread lock:", NTHREAD, "\n")
print(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")))

# ══════════════════════════════════════════════════════════════
# Libraries
# ══════════════════════════════════════════════════════════════
library(tidyverse)
library(ranger)
library(pROC)

# ══════════════════════════════════════════════════════════════
# Global configs
# ══════════════════════════════════════════════════════════════
scenarios_to_run <- c(
  "A_no_confound",
  "B_batch_confound",
  "C_sparse_small",
  "D_large_effect"
)
reps <- c(1L, 2L, 3L)

comparisons <- list(
  list(name = "Healthy_vs_Disease_mild",
       groups = c("Healthy", "Disease_mild")),
  list(name = "Healthy_vs_Disease_severe",
       groups = c("Healthy", "Disease_severe")),
  list(name = "Disease_mild_vs_Disease_severe",
       groups = c("Disease_mild", "Disease_severe"))
)

TOP_K     <- 50L
B         <- 30L
SUBSAMPLE <- 0.8

# Output roots (按你指定)
OUT_RUN_ROOT  <- file.path("benchmark", "rf", "rf_results_naive")  # 模型运行结果
OUT_SUMM_ROOT <- file.path("benchmark", "results")                # 汇总指标

dir.create(OUT_RUN_ROOT,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SUMM_ROOT, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════
# Utils
# ══════════════════════════════════════════════════════════════
clr_transform <- function(mat) {
  # mat: samples × taxa
  mat_pseudo <- mat + 1e-6
  t(apply(mat_pseudo, 1, function(x) log(x) - mean(log(x))))
}

jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# ══════════════════════════════════════════════════════════════
# Data IO
#   ✅ RF 只做格式化，不过滤、不变换
#   ✅ 过滤和 CLR 在训练折内部完成，避免测试集信息泄露
# ══════════════════════════════════════════════════════════════
generate_method_inputs <- function(scenario_full, rep_id) {
  
  counts_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_counts_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  relab_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_relab_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  metadata_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_metadata_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  
  counts_matrix <- as.matrix(counts_df[, -1])
  rownames(counts_matrix) <- counts_df$taxa
  
  relab_matrix <- as.matrix(relab_df[, -1])
  rownames(relab_matrix) <- relab_df$taxa
  
  list(
    rf = list(
      counts   = counts_matrix,   # 原始计数，未过滤
      relab    = relab_matrix,    # 原始相对丰度，未过滤
      metadata = metadata_df
    )
  )
}

load_rep_data <- function(scenario_full, rep_id) {
  inputs <- generate_method_inputs(scenario_full, rep_id)
  list(
    relab  = inputs$rf$relab,
    counts = inputs$rf$counts,
    meta   = inputs$rf$metadata
  )
}

# ══════════════════════════════════════════════════════════════
# Truth union (用于 metrics：Hit/Recall)
# ══════════════════════════════════════════════════════════════
truth <- read.table(
  "mock_benchmark_data/all_truth_tables.tsv",
  header = TRUE, sep = "\t"
)

build_truth_union <- function(truth_df, scenario_name) {
  all_comp_defs <- list(
    list(a = "Healthy",      b = "Disease_mild"),
    list(a = "Healthy",      b = "Disease_severe"),
    list(a = "Disease_mild", b = "Disease_severe")
  )
  truth_df %>%
    dplyr::filter(scenario == scenario_name) %>%
    dplyr::mutate(groups_list = stringr::str_split(present_in_groups, ",")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(is_tp_any = any(purrr::map_lgl(all_comp_defs, function(comp) {
      xor(comp$a %in% groups_list, comp$b %in% groups_list)
    }))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_tp_any) %>%
    dplyr::mutate(taxa_short = stringr::str_extract(taxa, "[^;]+$") %>% stringr::str_trim()) %>%
    dplyr::pull(taxa_short) %>%
    unique()
}

# ══════════════════════════════════════════════════════════════
# Fold TRAIN preparation (训练折内部做过滤；无泄露)
# ══════════════════════════════════════════════════════════════
prepare_rf_train_pair <- function(scenario_full, train_reps, grp) {
  
  train_list   <- purrr::map(train_reps, ~load_rep_data(scenario_full, .x))
  train_relab  <- do.call(cbind, purrr::map(train_list, "relab"))
  train_counts <- do.call(cbind, purrr::map(train_list, "counts"))
  train_meta   <- dplyr::bind_rows(purrr::map(train_list, "meta")) %>%
    dplyr::mutate(sample_id = as.character(sample_id),
                  phenotype = as.character(phenotype))
  
  # subset samples for this pair
  tr_idx     <- train_meta$phenotype %in% grp
  train_meta <- train_meta[tr_idx, , drop = FALSE]
  if (nrow(train_meta) == 0) return(NULL)
  
  # align columns to meta sample_id
  common <- intersect(colnames(train_relab), train_meta$sample_id)
  if (length(common) < 4) return(NULL)
  
  train_meta   <- train_meta[match(common, train_meta$sample_id), , drop = FALSE]
  train_relab  <- train_relab[, common, drop = FALSE]
  train_counts <- train_counts[, common, drop = FALSE]
  
  # feasibility
  tab <- table(train_meta$phenotype)
  if (!all(grp %in% names(tab)) || any(tab[grp] < 2)) return(NULL)
  
  # training-only filtering
  prevalence <- rowMeans(train_counts > 0)
  mean_relab <- rowMeans(train_relab)
  keep_taxa  <- prevalence >= 0.1 & mean_relab >= 1e-4
  
  train_relab_f <- train_relab[keep_taxa, , drop = FALSE]
  if (nrow(train_relab_f) == 0) return(NULL)
  
  # build clean_names mapping on training taxa set (this fold)
  taxa_short_names <- rownames(train_relab_f) %>%
    stringr::str_extract("[^;]+$") %>%
    stringr::str_trim()
  clean_names <- make.names(taxa_short_names, unique = TRUE)
  name_map    <- setNames(taxa_short_names, clean_names)  # clean -> taxa_short
  
  list(
    relab_f     = train_relab_f,   # taxa × samples
    meta        = train_meta,      # aligned with columns
    clean_names = clean_names,
    name_map    = name_map
  )
}

# ══════════════════════════════════════════════════════════════
# Train RF on given relab/meta and return top_k taxa_short (Naive)
# ══════════════════════════════════════════════════════════════
train_rf_topk_taxa <- function(relab_mat, meta_df, clean_names, name_map, grp,
                               top_k = TOP_K, seed = 42L) {
  
  common <- intersect(colnames(relab_mat), meta_df$sample_id)
  if (length(common) < 4) return(character(0))
  
  relab_mat <- relab_mat[, common, drop = FALSE]
  meta_df   <- meta_df[match(common, meta_df$sample_id), , drop = FALSE]
  
  tab <- table(meta_df$phenotype)
  if (!all(grp %in% names(tab)) || any(tab[grp] < 2)) return(character(0))
  
  X <- clr_transform(t(relab_mat))    # samples × taxa
  colnames(X) <- clean_names
  
  label_counts  <- table(meta_df$phenotype)
  class_weights <- (1 / label_counts) / sum(1 / label_counts)
  class_weights <- class_weights[grp]
  
  train_df <- as.data.frame(X)
  train_df$label <- factor(meta_df$phenotype, levels = grp)
  
  rf_model <- tryCatch(
    ranger::ranger(
      label ~ .,
      data          = train_df,
      num.trees     = 500,
      importance    = "permutation",
      class.weights = as.numeric(class_weights),
      probability   = TRUE,
      num.threads   = NTHREAD,
      seed          = as.integer(seed)
    ),
    error = function(e) NULL
  )
  if (is.null(rf_model)) return(character(0))
  
  imp <- sort(ranger::importance(rf_model), decreasing = TRUE)
  top_clean <- head(names(imp), top_k)
  unname(name_map[top_clean]) %>% as.character() %>% unique()
}

# ══════════════════════════════════════════════════════════════
# Main LORO fold run (one scenario + one test_rep + one pair)
#   - 写入到：benchmark/rf/rf_results_naive/{scenario}_rep{test_rep}/
#   - 返回结果包含 timing 信息
#   - compute_sec 只包含 ranger() 模型拟合时间
#   - 修改：同时保存比例值和分子分母，支持 Micro/Macro 平均
# ══════════════════════════════════════════════════════════════
run_rf_naive_fold_pair <- function(scenario_full, test_rep, comparison_def, truth_union,
                                   top_k = TOP_K) {
  
  train_reps <- setdiff(reps, test_rep)
  grp        <- comparison_def$groups
  comp_name  <- comparison_def$name
  
  prep <- prepare_rf_train_pair(scenario_full, train_reps, grp)
  if (is.null(prep)) return(NULL)
  
  # test set (only for AUC/metrics), align to training taxa set
  test_data <- load_rep_data(scenario_full, test_rep)
  test_meta <- test_data$meta %>%
    dplyr::mutate(sample_id = as.character(sample_id),
                  phenotype = as.character(phenotype))
  te_idx <- test_meta$phenotype %in% grp
  test_meta <- test_meta[te_idx, , drop = FALSE]
  if (nrow(test_meta) == 0) return(NULL)
  
  # take test relab and align taxa to training relab_f
  common_test <- intersect(colnames(test_data$relab), test_meta$sample_id)
  test_meta <- test_meta[match(common_test, test_meta$sample_id), , drop = FALSE]
  test_relab <- test_data$relab[, common_test, drop = FALSE]
  
  missing_taxa <- setdiff(rownames(prep$relab_f), rownames(test_relab))
  if (length(missing_taxa) > 0) {
    zero_rows <- matrix(0, nrow = length(missing_taxa), ncol = ncol(test_relab),
                        dimnames = list(missing_taxa, colnames(test_relab)))
    test_relab <- rbind(test_relab, zero_rows)
  }
  test_relab_f <- test_relab[rownames(prep$relab_f), , drop = FALSE]
  
  # Prepare data for ranger (这部分不计入 compute_sec)
  X_train <- clr_transform(t(prep$relab_f))
  colnames(X_train) <- prep$clean_names
  X_test  <- clr_transform(t(test_relab_f))
  colnames(X_test) <- prep$clean_names
  
  label_counts  <- table(prep$meta$phenotype)
  class_weights <- (1 / label_counts) / sum(1 / label_counts)
  class_weights <- class_weights[grp]
  
  train_df <- as.data.frame(X_train)
  train_df$label <- factor(prep$meta$phenotype, levels = grp)
  
  # ✅ compute_sec 开始：只计时 ranger() 模型拟合
  t_compute_start <- proc.time()["elapsed"]
  
  rf_model <- tryCatch(
    ranger::ranger(
      label ~ .,
      data          = train_df,
      num.trees     = 500,
      importance    = "permutation",
      class.weights = as.numeric(class_weights),
      probability   = TRUE,
      num.threads   = NTHREAD,
      seed          = 42L
    ),
    error = function(e) NULL
  )
  
  t_compute_end <- proc.time()["elapsed"]
  compute_sec <- as.numeric(t_compute_end - t_compute_start)
  # ✅ compute_sec 结束
  
  if (is.null(rf_model)) return(NULL)
  
  imp <- sort(ranger::importance(rf_model), decreasing = TRUE)
  
  importance_df <- tibble(
    clean_name     = names(imp),
    taxa_short     = prep$name_map[names(imp)],
    importance     = as.numeric(imp),
    rank           = seq_along(imp),
    selected_top_k = rank <= top_k,
    in_truth       = taxa_short %in% truth_union,
    scenario       = scenario_full,
    comparison     = comp_name,
    test_rep       = as.integer(test_rep),
    mode           = "naive"
  )
  
  # AUC (不计入 compute_sec)
  pred_prob <- predict(rf_model, data = as.data.frame(X_test))$predictions
  pos_col   <- which(colnames(pred_prob) == grp[2])
  if (length(pos_col) == 0) pos_col <- 2L
  
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(
      response  = factor(test_meta$phenotype, levels = grp),
      predictor = pred_prob[, pos_col],
      quiet = TRUE
    )))
  }, error = function(e) NA_real_)
  
  # ✅ 修改后的 metrics_df：同时保存比例值和分子分母
  compute_hit_details <- function(k) {
    top_taxa <- importance_df %>% filter(rank <= k) %>% pull(taxa_short)
    if (length(truth_union) == 0) {
      return(list(num = NA_integer_, denom = NA_integer_, ratio = NA_real_))
    }
    num <- length(intersect(top_taxa, truth_union))
    denom <- min(k, length(truth_union))
    ratio <- num / denom
    return(list(num = num, denom = denom, ratio = ratio))
  }
  
  hit_10 <- compute_hit_details(10L)
  hit_20 <- compute_hit_details(20L)
  hit_50 <- compute_hit_details(50L)
  
  recall_50_num <- if(length(truth_union) == 0) NA_integer_ else {
    top_taxa <- importance_df %>% filter(selected_top_k) %>% pull(taxa_short)
    length(intersect(top_taxa, truth_union))
  }
  recall_50_denom <- length(truth_union)
  recall_50_ratio <- if(recall_50_denom > 0) recall_50_num / recall_50_denom else NA_real_
  
  metrics_df <- tibble(
    # 标识信息
    scenario      = scenario_full,
    comparison    = comp_name,
    test_rep      = as.integer(test_rep),
    mode          = "naive",
    
    # AUC
    AUC           = auc_val,
    
    # Hit@10: 比例 + 分子 + 分母
    Hit_at_10     = hit_10$ratio,
    hit_10_num    = hit_10$num,
    hit_10_denom  = hit_10$denom,
    
    # Hit@20
    Hit_at_20     = hit_20$ratio,
    hit_20_num    = hit_20$num,
    hit_20_denom  = hit_20$denom,
    
    # Hit@50
    Hit_at_50     = hit_50$ratio,
    hit_50_num    = hit_50$num,
    hit_50_denom  = hit_50$denom,
    
    # Recall@50
    Recall_at_50  = recall_50_ratio,
    recall_50_num = recall_50_num,
    recall_50_denom = recall_50_denom,
    
    # 其他信息
    n_truth       = length(truth_union),
    n_features    = nrow(prep$relab_f),
    compute_sec   = compute_sec
  )
  
  list(importance_df = importance_df, metrics_df = metrics_df)
}

# ══════════════════════════════════════════════════════════════
# Run main analysis with per-fold saving (按你要求的目录结构)
# 计时到每个场景的每个 rep（不细分 comparison）
# ══════════════════════════════════════════════════════════════
all_importance <- list()
all_metrics    <- list()
timing_records <- list()  # 收集每个 rep 的 timing 信息

for (sc in scenarios_to_run) {
  
  message("\n▶ scenario: ", sc)
  truth_union <- build_truth_union(truth, sc)
  
  for (test_rep in reps) {
    
    # ✅ 每个 rep 的总计时开始
    t_rep_wall_start <- proc.time()["elapsed"]
    t_rep_compute_sum <- 0  # 累加该 rep 下所有 comparison 的 compute_sec
    
    message("  processing test_rep: ", test_rep)
    
    fold_dir <- file.path(OUT_RUN_ROOT, paste0(sc, "_rep", test_rep))
    dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
    
    fold_metrics_list <- list()
    
    for (comp in comparisons) {
      
      comp_name <- comp$name
      message("    comparison: ", comp_name)
      
      res <- tryCatch(
        run_rf_naive_fold_pair(sc, test_rep, comp, truth_union, top_k = TOP_K),
        error = function(e) { message("      !! ERROR: ", e$message); NULL }
      )
      if (is.null(res)) next
      
      # per-fold per-comparison save
      readr::write_csv(
        res$importance_df,
        file.path(fold_dir, paste0(comp_name, "_importance.csv"))
      )
      readr::write_csv(
        res$metrics_df,
        file.path(fold_dir, paste0(comp_name, "_metrics.csv"))
      )
      
      # accumulate
      all_importance <- c(all_importance, list(res$importance_df))
      all_metrics    <- c(all_metrics,    list(res$metrics_df))
      fold_metrics_list <- c(fold_metrics_list, list(res$metrics_df))
      
      # 累加 compute_sec
      t_rep_compute_sum <- t_rep_compute_sum + res$metrics_df$compute_sec
    }
    
    # fold summary inside this folder
    if (length(fold_metrics_list) > 0) {
      fold_metrics_df <- bind_rows(fold_metrics_list)
      readr::write_csv(fold_metrics_df, file.path(fold_dir, "all_folds_summary.csv"))
    } else {
      readr::write_csv(tibble(), file.path(fold_dir, "all_folds_summary.csv"))
    }
    
    # ✅ 每个 rep 的总计时结束
    t_rep_wall_end <- proc.time()["elapsed"]
    rep_wall_sec <- as.numeric(t_rep_wall_end - t_rep_wall_start)
    
    # 记录该 rep 的 timing
    timing_records[[length(timing_records) + 1]] <- tibble(
      method      = "RF_naive",
      scenario    = sc,
      test_rep    = as.integer(test_rep),
      wall_sec    = rep_wall_sec,
      compute_sec = t_rep_compute_sum
    )
  }
}

# ══════════════════════════════════════════════════════════════
# 汇总部分：按照你的要求修改
# - rf_naive_metrics.csv: 粒度到 comparison (Micro Hit/Recall + Macro AUC)
# - rf_naive_metrics_summary.csv: 粒度到 scenario
# ══════════════════════════════════════════════════════════════

# global summaries -> benchmark/results/
rf_importance_all <- bind_rows(all_importance)
rf_metrics_all    <- bind_rows(all_metrics)

# ✅ 1. 粒度到 comparison 的汇总 (rf_naive_metrics.csv)
#    - Hit@K 和 Recall@K: Micro 平均 (总分子/总分母)
#    - AUC: Macro 平均 (因为无法做真正的 Micro)
rf_metrics_comparison <- rf_metrics_all %>%
  group_by(scenario, comparison, mode) %>%
  summarise(
    # Micro 平均 for Hit@K and Recall@K
    Hit_at_10 = sum(hit_10_num, na.rm = TRUE) / sum(hit_10_denom, na.rm = TRUE),
    Hit_at_20 = sum(hit_20_num, na.rm = TRUE) / sum(hit_20_denom, na.rm = TRUE),
    Hit_at_50 = sum(hit_50_num, na.rm = TRUE) / sum(hit_50_denom, na.rm = TRUE),
    Recall_at_50 = sum(recall_50_num, na.rm = TRUE) / sum(recall_50_denom, na.rm = TRUE),
    
    # Macro 平均 for AUC
    AUC = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    
    # 汇总信息
    n_folds = n(),
    total_truth = sum(n_truth, na.rm = TRUE),
    n_features = first(n_features),
    
    .groups = "drop"
  )

# ✅ 2. 粒度到 scenario 的汇总 (rf_naive_metrics_summary.csv)
#    - 对同一 scenario 下的不同 comparison 再平均
rf_metrics_scenario <- rf_metrics_comparison %>%
  group_by(scenario, mode) %>%
  summarise(
    # Hit@K: 对 comparison 的 Micro 平均再做平均
    Hit_at_10_mean = mean(Hit_at_10, na.rm = TRUE),
    Hit_at_10_sd = sd(Hit_at_10, na.rm = TRUE),
    
    Hit_at_20_mean = mean(Hit_at_20, na.rm = TRUE),
    Hit_at_20_sd = sd(Hit_at_20, na.rm = TRUE),
    
    Hit_at_50_mean = mean(Hit_at_50, na.rm = TRUE),
    Hit_at_50_sd = sd(Hit_at_50, na.rm = TRUE),
    
    Recall_at_50_mean = mean(Recall_at_50, na.rm = TRUE),
    Recall_at_50_sd = sd(Recall_at_50, na.rm = TRUE),
    
    # AUC: 对 comparison 的 Macro 平均再做平均
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    
    # 汇总信息
    n_comparisons = n(),
    total_folds = sum(n_folds, na.rm = TRUE),
    n_features = first(n_features),
    
    .groups = "drop"
  )

# 保存文件
readr::write_csv(rf_metrics_all,
                 file.path(OUT_SUMM_ROOT, "rf_naive_metrics_all.csv"))
readr::write_csv(rf_metrics_comparison,
                 file.path(OUT_SUMM_ROOT, "rf_naive_metrics.csv"))
readr::write_csv(rf_metrics_scenario,
                 file.path(OUT_SUMM_ROOT, "rf_naive_metrics_summary.csv"))

# 新增：timing 汇总（按 scenario + rep 粒度）
timing_df <- bind_rows(timing_records)

# 按 scenario 汇总
timing_summary <- timing_df %>%
  group_by(method, scenario) %>%
  summarise(
    wall_mean_sec    = mean(wall_sec,      na.rm = TRUE),
    wall_median_sec  = median(wall_sec,    na.rm = TRUE),
    wall_min_sec     = min(wall_sec,       na.rm = TRUE),
    wall_max_sec     = max(wall_sec,       na.rm = TRUE),
    compute_mean_sec = mean(compute_sec,   na.rm = TRUE),
    compute_median_sec = median(compute_sec, na.rm = TRUE),
    compute_min_sec  = min(compute_sec,    na.rm = TRUE),
    compute_max_sec  = max(compute_sec,    na.rm = TRUE),
    n_reps           = n(),
    .groups = "drop"
  )

cat("\n── RF Naive 运行时间汇总（按场景）──\n")
print(timing_summary)

cat("\n── RF Naive 运行时间明细（每个 rep）──\n")
print(timing_df)

# 保存 timing 结果
readr::write_csv(timing_df,      file.path(OUT_SUMM_ROOT, "rf_naive_timing_raw.csv"))
readr::write_csv(timing_summary, file.path(OUT_SUMM_ROOT, "rf_naive_timing_summary.csv"))

cat("\n", strrep("=", 62), "\n")
cat("RF Naive 主分析完成！\n")
cat("Timing 汇总已保存至:", OUT_SUMM_ROOT, "\n")
cat("  - rf_naive_timing_raw.csv (每个 rep 的详细时间)\n")
cat("  - rf_naive_timing_summary.csv (按场景汇总)\n")
cat("\n新增的指标汇总文件：\n")
cat("  - rf_naive_metrics_all.csv (原始数据，含分子分母)\n")
cat("  - rf_naive_metrics.csv (按 comparison 汇总: Micro Hit/Recall + Macro AUC)\n")
cat("  - rf_naive_metrics_summary.csv (按 scenario 汇总)\n")
cat(strrep("=", 62), "\n")




# ══════════════════════════════════════════════════════════════
# RF-adjusted (Scenario B only)
# - LORO CV (rep as fold): train on 2 reps, test on held-out rep
# - In-fold filtering + CLR
# - Batch correction on feature matrix BEFORE RF (no leakage)
# - Jaccard stability aligned to fold; subsampling WITHOUT replacement
#   stratified by phenotype × batch
# Outputs:
#   benchmark/rf/rf_results_adjusted/{scenario}_rep{f}/...
#   benchmark/results/rf_adjusted_*.csv
# ══════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════
# Thread lock
# ══════════════════════════════════════════════════════════════
NTHREAD <- 1L

Sys.setenv(
  OMP_NUM_THREADS        = NTHREAD,
  OPENBLAS_NUM_THREADS   = NTHREAD,
  MKL_NUM_THREADS        = NTHREAD,
  VECLIB_MAXIMUM_THREADS = NTHREAD,
  NUMEXPR_NUM_THREADS    = NTHREAD
)
options(mc.cores = NTHREAD)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(NTHREAD)

cat("Thread lock:", NTHREAD, "\n")

# ══════════════════════════════════════════════════════════════
# Libraries
# ══════════════════════════════════════════════════════════════
library(tidyverse)
library(ranger)
library(pROC)

# ══════════════════════════════════════════════════════════════
# Configs
# ══════════════════════════════════════════════════════════════
scenario_B <- "B_batch_confound"
reps <- c(1L, 2L, 3L)

comparisons <- list(
  list(name = "Healthy_vs_Disease_mild",
       groups = c("Healthy", "Disease_mild")),
  list(name = "Healthy_vs_Disease_severe",
       groups = c("Healthy", "Disease_severe")),
  list(name = "Disease_mild_vs_Disease_severe",
       groups = c("Disease_mild", "Disease_severe"))
)

TOP_K     <- 50L
B         <- 30L
SUBSAMPLE <- 0.8

# Output roots
OUT_RUN_ROOT  <- file.path("benchmark", "rf", "rf_results_adjusted")
OUT_SUMM_ROOT <- file.path("benchmark", "results")
dir.create(OUT_RUN_ROOT,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SUMM_ROOT, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════
# Utils
# ══════════════════════════════════════════════════════════════
clr_transform <- function(mat) {
  # mat: samples × taxa
  mat_pseudo <- mat + 1e-6
  t(apply(mat_pseudo, 1, function(x) log(x) - mean(log(x))))
}

jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# Find batch column in metadata (robust)
get_batch_col <- function(meta_df) {
  candidates <- c("batch", "Batch", "batch_id", "BatchID", "run", "plate", "sequencing_batch")
  hit <- candidates[candidates %in% colnames(meta_df)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

# ══════════════════════════════════════════════════════════════
# Batch adjustment (location-only centering in CLR space)
# Fit on training only; apply to train/test
# X: samples × features
# ══════════════════════════════════════════════════════════════
fit_batch_centering <- function(X_train, batch_train) {
  X_train <- as.matrix(X_train)
  batch_train <- factor(batch_train)
  lvls <- levels(batch_train)
  
  overall_mean <- colMeans(X_train)
  
  # levels × features
  batch_means <- vapply(lvls, function(bl) {
    idx <- which(batch_train == bl)
    if (length(idx) == 0) rep(NA_real_, ncol(X_train)) else colMeans(X_train[idx, , drop = FALSE])
  }, FUN.VALUE = numeric(ncol(X_train))) |>
    t()
  
  rownames(batch_means) <- lvls
  colnames(batch_means) <- colnames(X_train)
  
  # shifts: levels × features
  shifts <- sweep(batch_means, 2, overall_mean, FUN = "-")
  colnames(shifts) <- colnames(X_train)
  
  list(levels = lvls, overall_mean = overall_mean, shifts = shifts)
}

apply_batch_centering <- function(X, batch, fit) {
  X <- as.matrix(X)
  batch_chr <- as.character(batch)
  idx <- match(batch_chr, fit$levels)
  
  # unseen batch -> shift=0
  shift_mat <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  colnames(shift_mat) <- colnames(X)
  
  ok <- !is.na(idx)
  if (any(ok)) {
    shift_mat[ok, ] <- fit$shifts[idx[ok], , drop = FALSE]
  }
  
  X - shift_mat
}

# ══════════════════════════════════════════════════════════════
# Data IO
# ✅ RF 输入只格式化，不过滤、不变换
# ══════════════════════════════════════════════════════════════
generate_method_inputs <- function(scenario_full, rep_id) {
  
  counts_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_counts_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  relab_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_relab_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  metadata_df <- readr::read_tsv(
    file.path("mock_benchmark_data",
              paste0("mock_metadata_", scenario_full, "_rep", rep_id, ".tsv")),
    show_col_types = FALSE
  )
  
  counts_matrix <- as.matrix(counts_df[, -1])
  rownames(counts_matrix) <- counts_df$taxa
  
  relab_matrix <- as.matrix(relab_df[, -1])
  rownames(relab_matrix) <- relab_df$taxa
  
  list(
    rf = list(
      counts   = counts_matrix,   # 原始计数，未过滤
      relab    = relab_matrix,    # 原始相对丰度，未过滤
      metadata = metadata_df
    )
  )
}

load_rep_data <- function(scenario_full, rep_id) {
  inputs <- generate_method_inputs(scenario_full, rep_id)
  list(
    relab  = inputs$rf$relab,
    counts = inputs$rf$counts,
    meta   = inputs$rf$metadata
  )
}

# ══════════════════════════════════════════════════════════════
# Truth union (用于 Hit/Recall)
# ══════════════════════════════════════════════════════════════
truth <- read.table(
  "mock_benchmark_data/all_truth_tables.tsv",
  header = TRUE, sep = "\t"
)

build_truth_union <- function(truth_df, scenario_name) {
  all_comp_defs <- list(
    list(a = "Healthy",      b = "Disease_mild"),
    list(a = "Healthy",      b = "Disease_severe"),
    list(a = "Disease_mild", b = "Disease_severe")
  )
  truth_df %>%
    dplyr::filter(scenario == scenario_name) %>%
    dplyr::mutate(groups_list = stringr::str_split(present_in_groups, ",")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(is_tp_any = any(purrr::map_lgl(all_comp_defs, function(comp) {
      xor(comp$a %in% groups_list, comp$b %in% groups_list)
    }))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_tp_any) %>%
    dplyr::mutate(taxa_short = stringr::str_extract(taxa, "[^;]+$") %>% stringr::str_trim()) %>%
    dplyr::pull(taxa_short) %>%
    unique()
}

# ══════════════════════════════════════════════════════════════
# Prepare training data for a fold+pair
# - training-only filtering (counts prevalence + relab mean)
# - keep metadata aligned
# - extract batch column
# ══════════════════════════════════════════════════════════════
prepare_train_pair_adjusted <- function(scenario_full, train_reps, grp) {
  
  train_list   <- purrr::map(train_reps, ~load_rep_data(scenario_full, .x))
  train_relab  <- do.call(cbind, purrr::map(train_list, "relab"))
  train_counts <- do.call(cbind, purrr::map(train_list, "counts"))
  train_meta   <- dplyr::bind_rows(purrr::map(train_list, "meta")) %>%
    dplyr::mutate(sample_id = as.character(sample_id),
                  phenotype = as.character(phenotype))
  
  batch_col <- get_batch_col(train_meta)
  if (is.null(batch_col)) {
    stop("No batch column found in metadata. Please add a 'batch' column (or common alias).")
  }
  
  # subset samples for this pair
  tr_idx     <- train_meta$phenotype %in% grp
  train_meta <- train_meta[tr_idx, , drop = FALSE]
  if (nrow(train_meta) == 0) return(NULL)
  
  # align columns to meta sample_id
  common <- intersect(colnames(train_relab), train_meta$sample_id)
  if (length(common) < 4) return(NULL)
  
  train_meta   <- train_meta[match(common, train_meta$sample_id), , drop = FALSE]
  train_relab  <- train_relab[, common, drop = FALSE]
  train_counts <- train_counts[, common, drop = FALSE]
  
  # feasibility (each class enough)
  tab <- table(train_meta$phenotype)
  if (!all(grp %in% names(tab)) || any(tab[grp] < 2)) return(NULL)
  
  # training-only filtering
  prevalence <- rowMeans(train_counts > 0)
  mean_relab <- rowMeans(train_relab)
  keep_taxa  <- prevalence >= 0.1 & mean_relab >= 1e-4
  
  train_relab_f <- train_relab[keep_taxa, , drop = FALSE]
  if (nrow(train_relab_f) == 0) return(NULL)
  
  taxa_short_names <- rownames(train_relab_f) %>%
    stringr::str_extract("[^;]+$") %>%
    stringr::str_trim()
  clean_names <- make.names(taxa_short_names, unique = TRUE)
  name_map    <- setNames(taxa_short_names, clean_names)  # clean -> taxa_short
  
  list(
    relab_f     = train_relab_f,            # taxa × samples
    meta        = train_meta,               # aligned
    batch_col   = batch_col,
    clean_names = clean_names,
    name_map    = name_map
  )
}

# ══════════════════════════════════════════════════════════════
# Train adjusted RF and return importance + model + batch-fit
# - Steps inside training fold:
#   CLR -> batch centering fit -> apply -> RF
# - Returns compute_sec (ranger fitting time only)
# ══════════════════════════════════════════════════════════════
fit_rf_adjusted <- function(relab_train_taxa_by_samp, meta_train,
                            clean_names, grp, batch_col, seed = 42L) {
  
  # Build CLR
  X_train_raw <- clr_transform(t(relab_train_taxa_by_samp))  # samples × taxa
  colnames(X_train_raw) <- clean_names
  
  batch_train <- meta_train[[batch_col]] %>% as.character()
  fit <- fit_batch_centering(X_train_raw, batch_train)
  X_train_adj <- apply_batch_centering(X_train_raw, batch_train, fit)
  
  label_counts  <- table(meta_train$phenotype)
  class_weights <- (1 / label_counts) / sum(1 / label_counts)
  class_weights <- class_weights[grp]
  
  train_df <- as.data.frame(X_train_adj)
  train_df$label <- factor(meta_train$phenotype, levels = grp)
  
  # ✅ compute_sec 开始：只计时 ranger() 模型拟合
  t_compute_start <- proc.time()["elapsed"]
  
  rf_model <- tryCatch(
    ranger::ranger(
      label ~ .,
      data          = train_df,
      num.trees     = 500,
      importance    = "permutation",
      class.weights = as.numeric(class_weights),
      probability   = TRUE,
      num.threads   = NTHREAD,
      seed          = as.integer(seed)
    ),
    error = function(e) NULL
  )
  
  t_compute_end <- proc.time()["elapsed"]
  compute_sec <- as.numeric(t_compute_end - t_compute_start)
  # ✅ compute_sec 结束
  
  list(model = rf_model, batch_fit = fit, compute_sec = compute_sec)
}

# ══════════════════════════════════════════════════════════════
# One fold + one pair: train adjusted RF, test, output dfs
# ══════════════════════════════════════════════════════════════
run_rf_adjusted_fold_pair <- function(scenario_full, test_rep, comparison_def, truth_union,
                                      top_k = TOP_K) {
  
  train_reps <- setdiff(reps, test_rep)
  grp        <- comparison_def$groups
  comp_name  <- comparison_def$name
  
  prep <- prepare_train_pair_adjusted(scenario_full, train_reps, grp)
  if (is.null(prep)) return(NULL)
  
  # Test set
  test_data <- load_rep_data(scenario_full, test_rep)
  test_meta <- test_data$meta %>%
    mutate(sample_id = as.character(sample_id),
           phenotype = as.character(phenotype))
  
  batch_col_test <- get_batch_col(test_meta)
  if (is.null(batch_col_test)) stop("No batch column found in test metadata.")
  # 强制同名列（用于 apply）
  if (batch_col_test != prep$batch_col) {
    # 若名字不一样但都存在，你可以改成重命名；这里简单处理为重命名 test 到 train 的 batch_col
    test_meta[[prep$batch_col]] <- test_meta[[batch_col_test]]
  }
  
  te_idx <- test_meta$phenotype %in% grp
  test_meta <- test_meta[te_idx, , drop = FALSE]
  if (nrow(test_meta) == 0) return(NULL)
  
  common_test <- intersect(colnames(test_data$relab), test_meta$sample_id)
  test_meta <- test_meta[match(common_test, test_meta$sample_id), , drop = FALSE]
  test_relab <- test_data$relab[, common_test, drop = FALSE]
  
  # Align taxa to training taxa set (missing -> 0)
  missing_taxa <- setdiff(rownames(prep$relab_f), rownames(test_relab))
  if (length(missing_taxa) > 0) {
    zero_rows <- matrix(0, nrow = length(missing_taxa), ncol = ncol(test_relab),
                        dimnames = list(missing_taxa, colnames(test_relab)))
    test_relab <- rbind(test_relab, zero_rows)
  }
  test_relab_f <- test_relab[rownames(prep$relab_f), , drop = FALSE]
  
  # Fit adjusted RF on training (includes compute_sec)
  fit <- fit_rf_adjusted(prep$relab_f, prep$meta, prep$clean_names, grp, prep$batch_col, seed = 42L)
  if (is.null(fit$model)) return(NULL)
  
  compute_sec <- fit$compute_sec
  
  # Importance
  imp <- sort(ranger::importance(fit$model), decreasing = TRUE)
  importance_df <- tibble(
    clean_name     = names(imp),
    taxa_short     = prep$name_map[names(imp)],
    importance     = as.numeric(imp),
    rank           = seq_along(imp),
    selected_top_k = rank <= top_k,
    in_truth       = taxa_short %in% truth_union,
    scenario       = scenario_full,
    comparison     = comp_name,
    test_rep       = as.integer(test_rep),
    mode           = "adjusted_batch_centering"
  )
  
  # Build adjusted X_test using training-fitted batch correction
  X_test_raw <- clr_transform(t(test_relab_f))  # samples × taxa
  colnames(X_test_raw) <- prep$clean_names
  batch_test <- test_meta[[prep$batch_col]] %>% as.character()
  X_test_adj <- apply_batch_centering(X_test_raw, batch_test, fit$batch_fit)
  
  pred_prob <- predict(fit$model, data = as.data.frame(X_test_adj))$predictions
  pos_col   <- which(colnames(pred_prob) == grp[2])
  if (length(pos_col) == 0) pos_col <- 2L
  
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(
      response  = factor(test_meta$phenotype, levels = grp),
      predictor = pred_prob[, pos_col],
      quiet = TRUE
    )))
  }, error = function(e) NA_real_)
  
  # ✅ 修改后的 metrics_df：同时保存比例值和分子分母
  compute_hit_details <- function(k) {
    top_taxa <- importance_df %>% filter(rank <= k) %>% pull(taxa_short)
    if (length(truth_union) == 0) {
      return(list(num = NA_integer_, denom = NA_integer_, ratio = NA_real_))
    }
    num <- length(intersect(top_taxa, truth_union))
    denom <- min(k, length(truth_union))
    ratio <- num / denom
    return(list(num = num, denom = denom, ratio = ratio))
  }
  
  hit_10 <- compute_hit_details(10L)
  hit_20 <- compute_hit_details(20L)
  hit_50 <- compute_hit_details(50L)
  
  recall_50_num <- if(length(truth_union) == 0) NA_integer_ else {
    top_taxa <- importance_df %>% filter(selected_top_k) %>% pull(taxa_short)
    length(intersect(top_taxa, truth_union))
  }
  recall_50_denom <- length(truth_union)
  recall_50_ratio <- if(recall_50_denom > 0) recall_50_num / recall_50_denom else NA_real_
  
  metrics_df <- tibble(
    scenario      = scenario_full,
    comparison    = comp_name,
    test_rep      = as.integer(test_rep),
    mode          = "adjusted_batch_centering",
    
    # AUC
    AUC           = auc_val,
    
    # Hit@10: 比例 + 分子 + 分母
    Hit_at_10     = hit_10$ratio,
    hit_10_num    = hit_10$num,
    hit_10_denom  = hit_10$denom,
    
    # Hit@20
    Hit_at_20     = hit_20$ratio,
    hit_20_num    = hit_20$num,
    hit_20_denom  = hit_20$denom,
    
    # Hit@50
    Hit_at_50     = hit_50$ratio,
    hit_50_num    = hit_50$num,
    hit_50_denom  = hit_50$denom,
    
    # Recall@50
    Recall_at_50  = recall_50_ratio,
    recall_50_num = recall_50_num,
    recall_50_denom = recall_50_denom,
    
    # 其他信息
    n_truth       = length(truth_union),
    n_features    = nrow(prep$relab_f),
    compute_sec   = compute_sec
  )
  
  list(importance_df = importance_df, metrics_df = metrics_df)
}

# ══════════════════════════════════════════════════════════════
# Run main LORO analysis (Scenario B only) with per-fold saving
# ✅ 计时模块：每个 rep 记录 wall_sec 和 compute_sec
# ══════════════════════════════════════════════════════════════
truth_union_B <- build_truth_union(truth, scenario_B)

all_importance <- list()
all_metrics    <- list()
timing_records <- list()  # 收集每个 rep 的 timing 信息

for (test_rep in reps) {
  
  # ✅ 每个 rep 的总计时开始
  t_rep_wall_start <- proc.time()["elapsed"]
  t_rep_compute_sum <- 0  # 累加该 rep 下所有 comparison 的 compute_sec
  
  fold_dir <- file.path(OUT_RUN_ROOT, paste0(scenario_B, "_rep", test_rep))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  
  fold_metrics_list <- list()
  message("\n▶ scenario: ", scenario_B, " | fold(test_rep): ", test_rep)
  
  for (comp in comparisons) {
    
    comp_name <- comp$name
    message("  pair: ", comp_name)
    
    res <- tryCatch(
      run_rf_adjusted_fold_pair(scenario_B, test_rep, comp, truth_union_B, top_k = TOP_K),
      error = function(e) { message("    !! ERROR: ", e$message); NULL }
    )
    if (is.null(res)) next
    
    readr::write_csv(res$importance_df, file.path(fold_dir, paste0(comp_name, "_importance.csv")))
    readr::write_csv(res$metrics_df,    file.path(fold_dir, paste0(comp_name, "_metrics.csv")))
    
    all_importance <- c(all_importance, list(res$importance_df))
    all_metrics    <- c(all_metrics,    list(res$metrics_df))
    fold_metrics_list <- c(fold_metrics_list, list(res$metrics_df))
    
    # 累加 compute_sec
    t_rep_compute_sum <- t_rep_compute_sum + res$metrics_df$compute_sec
  }
  
  # same filename convention as your naive script
  if (length(fold_metrics_list) > 0) {
    readr::write_csv(bind_rows(fold_metrics_list),
                     file.path(fold_dir, "all_folds_summary.csv"))
  } else {
    readr::write_csv(tibble(), file.path(fold_dir, "all_folds_summary.csv"))
  }
  
  # ✅ 每个 rep 的总计时结束
  t_rep_wall_end <- proc.time()["elapsed"]
  rep_wall_sec <- as.numeric(t_rep_wall_end - t_rep_wall_start)
  
  # 记录该 rep 的 timing
  timing_records[[length(timing_records) + 1]] <- tibble(
    method      = "RF_adjusted",
    scenario    = scenario_B,
    test_rep    = as.integer(test_rep),
    wall_sec    = rep_wall_sec,
    compute_sec = t_rep_compute_sum
  )
}

rf_importance_all <- bind_rows(all_importance)
rf_metrics_all    <- bind_rows(all_metrics)

# ✅ 1. 粒度到 comparison 的汇总 (rf_adjusted_metrics.csv)
#    - Hit@K 和 Recall@K: Micro 平均 (总分子/总分母)
#    - AUC: Macro 平均
rf_metrics_comparison <- rf_metrics_all %>%
  group_by(scenario, comparison, mode) %>%
  summarise(
    # Micro 平均 for Hit@K and Recall@K
    Hit_at_10 = sum(hit_10_num, na.rm = TRUE) / sum(hit_10_denom, na.rm = TRUE),
    Hit_at_20 = sum(hit_20_num, na.rm = TRUE) / sum(hit_20_denom, na.rm = TRUE),
    Hit_at_50 = sum(hit_50_num, na.rm = TRUE) / sum(hit_50_denom, na.rm = TRUE),
    Recall_at_50 = sum(recall_50_num, na.rm = TRUE) / sum(recall_50_denom, na.rm = TRUE),
    
    # Macro 平均 for AUC
    AUC = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    
    # 汇总信息
    n_folds = n(),
    total_truth = sum(n_truth, na.rm = TRUE),
    n_features = first(n_features),
    
    .groups = "drop"
  )

# ✅ 2. 粒度到 scenario 的汇总 (rf_adjusted_metrics_summary.csv)
#    - 对同一 scenario 下的不同 comparison 再平均
rf_metrics_scenario <- rf_metrics_comparison %>%
  group_by(scenario, mode) %>%
  summarise(
    # Hit@K: 对 comparison 的 Micro 平均再做平均
    Hit_at_10_mean = mean(Hit_at_10, na.rm = TRUE),
    Hit_at_10_sd = sd(Hit_at_10, na.rm = TRUE),
    
    Hit_at_20_mean = mean(Hit_at_20, na.rm = TRUE),
    Hit_at_20_sd = sd(Hit_at_20, na.rm = TRUE),
    
    Hit_at_50_mean = mean(Hit_at_50, na.rm = TRUE),
    Hit_at_50_sd = sd(Hit_at_50, na.rm = TRUE),
    
    Recall_at_50_mean = mean(Recall_at_50, na.rm = TRUE),
    Recall_at_50_sd = sd(Recall_at_50, na.rm = TRUE),
    
    # AUC: 对 comparison 的 Macro 平均再做平均
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    
    # 汇总信息
    n_comparisons = n(),
    total_folds = sum(n_folds, na.rm = TRUE),
    n_features = first(n_features),
    
    .groups = "drop"
  )

# 保存文件
readr::write_csv(rf_metrics_all,
                 file.path(OUT_SUMM_ROOT, "rf_adjusted_metrics_all.csv"))
readr::write_csv(rf_metrics_comparison,
                 file.path(OUT_SUMM_ROOT, "rf_adjusted_metrics.csv"))
readr::write_csv(rf_metrics_scenario,
                 file.path(OUT_SUMM_ROOT, "rf_adjusted_metrics_summary.csv"))

# ✅ 新增：timing 汇总
timing_df <- bind_rows(timing_records)

# 按 scenario 汇总
timing_summary <- timing_df %>%
  group_by(method, scenario) %>%
  summarise(
    wall_mean_sec    = mean(wall_sec,      na.rm = TRUE),
    wall_median_sec  = median(wall_sec,    na.rm = TRUE),
    wall_min_sec     = min(wall_sec,       na.rm = TRUE),
    wall_max_sec     = max(wall_sec,       na.rm = TRUE),
    compute_mean_sec = mean(compute_sec,   na.rm = TRUE),
    compute_median_sec = median(compute_sec, na.rm = TRUE),
    compute_min_sec  = min(compute_sec,    na.rm = TRUE),
    compute_max_sec  = max(compute_sec,    na.rm = TRUE),
    n_reps           = n(),
    .groups = "drop"
  )

cat("\n── RF Adjusted 运行时间汇总（按场景）──\n")
print(timing_summary)

cat("\n── RF Adjusted 运行时间明细（每个 rep）──\n")
print(timing_df)

# 保存 timing 结果
readr::write_csv(timing_df,      file.path(OUT_SUMM_ROOT, "rf_adjusted_timing_raw.csv"))
readr::write_csv(timing_summary, file.path(OUT_SUMM_ROOT, "rf_adjusted_timing_summary.csv"))

cat("\n", strrep("=", 62), "\n")
cat("RF Adjusted 主分析完成！\n")
cat("Timing 汇总已保存至:", OUT_SUMM_ROOT, "\n")
cat("  - rf_adjusted_timing_raw.csv (每个 rep 的详细时间)\n")
cat("  - rf_adjusted_timing_summary.csv (按场景汇总)\n")
cat("\n指标汇总文件：\n")
cat("  - rf_adjusted_metrics_all.csv (原始数据，含分子分母)\n")
cat("  - rf_adjusted_metrics.csv (按 comparison 汇总: Micro Hit/Recall + Macro AUC)\n")
cat("  - rf_adjusted_metrics_summary.csv (按 scenario 汇总)\n")
cat(strrep("=", 62), "\n")