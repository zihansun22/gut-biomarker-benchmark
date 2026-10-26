
# ══════════════════════════════════════════════════════════════
# Benchmark thread lock（必须在在可能触发并行计算之前设置）
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

# ══════════════════════════════════════════════════════════════
# 加载包
# ══════════════════════════════════════════════════════════════
library(Maaslin2)
library(tidyverse)

# ══════════════════════════════════════════════════════════════
# Thread Lock（library() 之后补充，控制已加载的 BLAS）
# ══════════════════════════════════════════════════════════════
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(NTHREAD)
  RhpcBLASctl::omp_set_num_threads(NTHREAD)
}

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(NTHREAD)
}

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

cat("Thread lock:", NTHREAD, "\n")
print(Sys.getenv(c("OMP_NUM_THREADS",
                   "OPENBLAS_NUM_THREADS",
                   "MKL_NUM_THREADS")))

# ══════════════════════════════════════════════════════════════
# 全局参数（只定义一次）
# ══════════════════════════════════════════════════════════════
scenarios_to_run <- c(
  "A_no_confound",
  "B_batch_confound",
  "C_sparse_small",
  "D_large_effect"
)
reps <- as.character(1:3)

QVAL_THR  <- 0.05
LFC_THR   <- 1        # MaAsLin2 不强制 LFC 过滤，保留供扩展
TOP_K     <- 50
B         <- 30
SUBSAMPLE <- 0.8

# ══════════════════════════════════════════════════════════════
# 工具函数（只定义一次）
# ══════════════════════════════════════════════════════════════
jaccard <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

ensure_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

empty_jaccard_tibble <- function() {
  tibble(
    method            = character(0),
    scenario          = character(0),
    rep               = integer(0),
    comparison        = character(0),
    comparison_set    = character(0),
    top_k_per_pair    = integer(0),
    n_pairs           = integer(0),
    top_k_total       = integer(0),
    jaccard_mean      = numeric(0),
    jaccard_median    = numeric(0),
    jaccard_sd        = numeric(0),
    jaccard_min       = numeric(0),
    jaccard_max       = numeric(0),
    n_boots           = integer(0),
    n_boots_valid     = integer(0),
    n_features_full   = integer(0)
  )
}

# ══════════════════════════════════════════════════════════════
# 核心函数：从 feature_table + meta 运行 MaAsLin2，返回显著特征名
# feature_table : 行 = 样本，列 = taxa（与主函数一致）
# meta_df       : 含 rownames = sample_id，含 phenotype 列
# debug         : TRUE 时打印维度信息
# ══════════════════════════════════════════════════════════════
run_maaslin2_core <- function(feature_table,
                              meta_df,
                              qval_cutoff = 0.05,
                              top_k       = NULL,
                              debug       = FALSE) {
  
  if (debug) {
    cat("  [DEBUG] feature_table dim:",
        paste(dim(feature_table), collapse = "x"), "\n")
  }
  
  # ── 先对齐样本，再过滤分组 ────────────────────────────────
  common_samples <- intersect(rownames(feature_table), rownames(meta_df))
  if (length(common_samples) < 4) return(character(0))
  
  feature_table <- feature_table[common_samples, , drop = FALSE]
  meta_df       <- meta_df[common_samples,       , drop = FALSE]
  
  # ── 过滤有效分组 ───────────────────────────────────────────
  group_counts <- table(meta_df$phenotype)
  valid_groups <- names(group_counts[group_counts >= 2])
  if (length(valid_groups) < 2) return(character(0))
  
  meta_df       <- meta_df[meta_df$phenotype %in% valid_groups, , drop = FALSE]
  feature_table <- feature_table[rownames(meta_df), , drop = FALSE]
  
  disease_pairs <- combn(valid_groups, 2, simplify = FALSE)
  all_features  <- character(0)
  
  for (pair in disease_pairs) {
    
    pair_samples <- rownames(meta_df)[meta_df$phenotype %in% pair]
    if (length(pair_samples) < 4) next
    
    pair_meta    <- meta_df[pair_samples,    , drop = FALSE]
    pair_feature <- feature_table[pair_samples, , drop = FALSE]
    
    # 去全 0 列
    nonzero_taxa <- colSums(pair_feature) > 0
    pair_feature <- pair_feature[, nonzero_taxa, drop = FALSE]
    if (ncol(pair_feature) == 0) next
    
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    
    tmp_dir <- tempfile(pattern = "maaslin2_core_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    
    fit <- tryCatch({
      capture.output(
        fit_obj <- Maaslin2(
          input_data       = pair_feature,
          input_metadata   = pair_meta,
          output           = tmp_dir,
          fixed_effects    = "phenotype",  # Naive: 只考虑phenotype
          normalization    = "TSS",
          transform        = "LOG",
          analysis_method  = "LM",
          correction       = "BH",
          max_significance = 0.25,
          min_abundance    = 0.0,
          min_prevalence   = 0.0,
          min_variance     = 0.0,
          standardize      = FALSE,
          cores            = NTHREAD, 
          plot_heatmap     = FALSE,
          plot_scatter     = FALSE
        ),
        type = "output"
      )
      fit_obj
    }, error = function(e) {
      if (debug) cat("  [core] MaAsLin2 失败:", e$message, "\n")
      NULL
    })
    
    if (is.null(fit)) next
    
    res_file <- file.path(tmp_dir, "all_results.tsv")
    if (!file.exists(res_file)) next
    
    res_df <- tryCatch(
      read_tsv(res_file, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(res_df) || nrow(res_df) == 0) next
    
    feats <- res_df %>%
      filter(!is.na(qval), qval < qval_cutoff) %>%
      arrange(desc(abs(coef)))
    
    if (!is.null(top_k)) feats <- slice_head(feats, n = top_k)
    
    all_features <- c(all_features, feats$feature)
  }
  
  unique(all_features)
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 MaAsLin2，含完整保存逻辑
# wall_sec   ：总耗时（含 I/O）
# compute_sec：纯模型拟合耗时（每个 pair 的 Maaslin2() 调用时间之和）
# ══════════════════════════════════════════════════════════════
run_mock_maaslin2_genus <- function(scenario_full, rep_id) {
  
  set.seed(42)
  
  cat("开始处理:", scenario_full, "rep", rep_id, "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  # 1. 读取数据 ────────────────────────────────────────────
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$maaslin2$abundance
  metadata_df   <- inputs$maaslin2$metadata
  
  cat("特征矩阵维度:", dim(counts_matrix), "\n")
  cat("样本数:", nrow(metadata_df), "\n")
  flush.console()
  
  # 2. 转置为 MaAsLin2 要求的格式：行 = 样本，列 = taxa ───
  counts_matrix <- ensure_matrix(counts_matrix)
  feature_table <- as.data.frame(t(counts_matrix))
  
  # 3. 整理 metadata ───────────────────────────────────────
  # Naive: 只保留 sample_id 和 phenotype，丢弃 batch
  meta <- metadata_df %>%
    dplyr::select(sample_id, phenotype) %>%
    tibble::column_to_rownames("sample_id")
  
  # 先对齐样本
  common_samples <- intersect(rownames(meta), rownames(feature_table))
  meta           <- meta[common_samples,          , drop = FALSE]
  feature_table  <- feature_table[common_samples, , drop = FALSE]
  
  cat("对齐后样本数:", length(common_samples), "\n")
  cat("分组分布:\n")
  print(table(meta$phenotype))
  flush.console()
  
  # 4. 枚举所有两两比较对 ──────────────────────────────────
  all_groups    <- sort(unique(as.character(meta$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n\n")
  flush.console()
  
  # ✅ 修改：输出目录改为 maaslin2_results_genus_naive
  base_result_dir <- file.path(
    "benchmark", "maaslin2", "maaslin2_results_genus_naive",
    paste0(scenario_full, "_rep", rep_id)
  )
  dir.create(base_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  compute_sec  <- 0  # 累加每个 pair 的 Maaslin2() 耗时
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]
    g2        <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    
    cat("比较", i, "/", length(disease_pairs), ":", comp_name, "\n")
    flush.console()
    
    pair_samples <- rownames(meta)[meta$phenotype %in% pair]
    if (length(pair_samples) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_meta    <- meta[pair_samples,          , drop = FALSE]
    pair_feature <- feature_table[pair_samples, , drop = FALSE]
    
    # 去全 0 列
    nonzero_taxa <- colSums(pair_feature) > 0
    pair_feature <- pair_feature[, nonzero_taxa, drop = FALSE]
    
    cat("  样本数:", nrow(pair_feature),
        "| 特征数:", ncol(pair_feature), "\n")
    cat("  分组分布:\n")
    print(table(pair_meta$phenotype))
    flush.console()
    
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    
    cat("  参考组:", ref_group, "\n")
    flush.console()
    
    comp_output_dir <- file.path(base_result_dir, comp_name)
    dir.create(comp_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ── 纯模型拟合计时（每个 pair 单独计，累加）────────────
    t_compute_start <- proc.time()["elapsed"]
    
    maaslin_res <- tryCatch({
      capture.output(
        fit <- Maaslin2(
          input_data       = pair_feature,
          input_metadata   = pair_meta,
          output           = comp_output_dir,
          fixed_effects    = "phenotype",  # Naive: 只考虑phenotype
          normalization    = "TSS",
          transform        = "LOG",
          analysis_method  = "LM",
          correction       = "BH",
          max_significance = 0.25,
          min_abundance    = 0.0,
          min_prevalence   = 0.0,
          min_variance     = 0.0,
          standardize      = FALSE,
          plot_heatmap     = FALSE,
          plot_scatter     = FALSE
        ),
        type = "output"
      )
      fit
    }, error = function(e) {
      cat("  MaAsLin2 失败:", e$message, "\n"); NULL
    })
    
    t_compute_end <- proc.time()["elapsed"]
    compute_sec   <- compute_sec +
      as.numeric(t_compute_end - t_compute_start)
    
    if (is.null(maaslin_res)) next
    
    res_file <- file.path(comp_output_dir, "all_results.tsv")
    if (!file.exists(res_file)) {
      cat("  结果文件不存在，跳过\n"); next
    }
    
    res_df <- read_tsv(res_file, show_col_types = FALSE) %>%
      mutate(
        comparison = comp_name,
        group_a    = g1,
        group_b    = g2,
        ref_group  = ref_group,
        scenario   = scenario_full,
        replicate  = rep_id
      )
    
    cat("  发现", nrow(res_df), "个特征结果，",
        sum(res_df$qval < 0.05, na.rm = TRUE), "个 FDR<0.05\n")
    flush.console()
    
    # I/O（不计入 compute_sec）
    write_csv(res_df,
              file.path(base_result_dir,
                        paste0(comp_name, "_results.csv")))
    results_list[[comp_name]] <- res_df
  }
  
  # 合并所有比较对结果（I/O，不计入 compute_sec）
  if (length(results_list) > 0) {
    all_results_df <- bind_rows(results_list)
    write_csv(all_results_df,
              file.path(base_result_dir, "all_comparisons_results.csv"))
    cat("结果已保存至:", base_result_dir, "\n")
  } else {
    cat("无有效比较结果\n")
  }
  
  t_wall_end <- proc.time()["elapsed"]
  wall_sec   <- as.numeric(t_wall_end - t_wall_start)
  
  cat("\n函数执行完成\n")
  cat("  总耗时   (wall_sec)   :", round(wall_sec,    2), "秒\n")
  cat("  纯模型   (compute_sec):", round(compute_sec, 2), "秒\n")
  cat("  共完成", length(results_list), "个比较\n")
  flush.console()
  
  return(list(
    results     = results_list,
    wall_sec    = wall_sec,
    compute_sec = compute_sec
  ))
}

# ══════════════════════════════════════════════════════════════
# 批量运行所有场景 + 收集 timing
# ══════════════════════════════════════════════════════════════
run_all <- function() {
  
  cat("\n=============================\n")
  cat("开始运行 MaAsLin2 Naive\n")
  cat("=============================\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n--- 场景:", scenario, "Rep:", rep, "---\n")
      flush.console()
      
      key <- paste0(scenario, "_rep", rep)
      
      out <- tryCatch(
        run_mock_maaslin2_genus(scenario, rep),
        error = function(e) { cat("失败:", e$message, "\n"); NULL }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out$results
        timing_records[[length(timing_records) + 1]] <- tibble(
          method      = "MaAsLin2",
          scenario    = scenario,
          rep         = as.integer(rep),
          wall_sec    = out$wall_sec,
          compute_sec = out$compute_sec
        )
      } else {
        results_all[[key]] <- NULL
      }
    }
  }
  
  cat("\nMaAsLin2 完成，成功:",
      sum(!sapply(results_all, is.null)), "/", length(results_all), "\n")
  
  timing_df <- bind_rows(timing_records)
  
  timing_summary <- timing_df %>%
    group_by(method, scenario) %>%
    summarise(
      wall_mean_sec    = mean(wall_sec,      na.rm = TRUE),
      wall_median_sec  = median(wall_sec,    na.rm = TRUE),
      wall_min_sec     = min(wall_sec,       na.rm = TRUE),
      wall_max_sec     = max(wall_sec,       na.rm = TRUE),
      compute_mean_sec = mean(compute_sec,   na.rm = TRUE),
      compute_med_sec  = median(compute_sec, na.rm = TRUE),
      compute_min_sec  = min(compute_sec,    na.rm = TRUE),
      compute_max_sec  = max(compute_sec,    na.rm = TRUE),
      n_runs           = n(),
      .groups = "drop"
    )
  
  cat("\n── MaAsLin2 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  # ✅ 修改：保存为 maaslin2_naive_timing_*
  saveRDS(timing_df,      "benchmark/results/maaslin2_naive_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/maaslin2_naive_timing_summary.rds")
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

out_maaslin2 <- run_all()

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV
# ══════════════════════════════════════════════════════════════
collect_maaslin2_results <- function(
    base_dir  = "benchmark/maaslin2/maaslin2_results_genus_naive",
    scenarios, reps) {
  
  all_rows <- list()
  
  for (scenario in scenarios) {
    for (rep in reps) {
      result_dir <- file.path(base_dir,
                              paste0(scenario, "_rep", rep))
      if (!dir.exists(result_dir)) next
      
      csv_files <- list.files(result_dir,
                              pattern    = "_results\\.csv$",
                              full.names = TRUE)
      csv_files <- csv_files[
        !grepl("all_comparisons_results\\.csv$", csv_files)
      ]
      
      for (f in csv_files) {
        df <- tryCatch(
          read_csv(f, show_col_types = FALSE),
          error = function(e) NULL
        )
        if (is.null(df) || nrow(df) == 0) next
        all_rows[[length(all_rows) + 1]] <- df
      }
    }
  }
  
  if (length(all_rows) == 0) { cat("未找到任何结果文件\n"); return(NULL) }
  bind_rows(all_rows)
}

maaslin2_summary_raw <- collect_maaslin2_results(
  scenarios = scenarios_to_run,
  reps      = reps
)

if (!is.null(maaslin2_summary_raw)) {
  cat("汇总行数:", nrow(maaslin2_summary_raw), "\n")
  print(head(maaslin2_summary_raw))
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  # ✅ 修改：保存为 maaslin2_naive_summary_*
  saveRDS(maaslin2_summary_raw,  "benchmark/results/maaslin2_naive_summary_raw.rds")
  write_csv(maaslin2_summary_raw, "benchmark/results/maaslin2_naive_summary_raw.csv")
  cat("\n汇总结果已保存\n")
}

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估
# ══════════════════════════════════════════════════════════════
run_maaslin2_jaccard <- function() {
  
  BASE_SEED <- 42L
  jaccard_records <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("MaAsLin2 Naive Jaccard 稳定性评估（pairwise + 场景内pair均值汇总）\n")
  cat("══════════════════════════════════════════════════════════════\n")
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n[Jaccard] 场景:", scenario, "Rep:", rep, "\n")
      flush.console()
      
      inputs <- tryCatch(
        generate_method_inputs(scenario, rep),
        error = function(e) {
          cat("数据读取失败:", e$message, "\n"); NULL
        }
      )
      if (is.null(inputs)) next
      
      counts_matrix <- inputs$maaslin2$abundance
      metadata_df   <- inputs$maaslin2$metadata
      
      # ── 硬防护 + 基本清洗 ─────────────────────────────
      counts_matrix <- ensure_matrix(counts_matrix)
      counts_matrix <- round(pmax(counts_matrix, 0))
      
      # 去全 0 样本
      col_keep      <- colSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # 去全 0 特征
      row_keep      <- rowSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # 转置为 samples × taxa
      feature_table <- as.data.frame(t(counts_matrix))
      
      # Naive: 只保留 phenotype
      meta <- metadata_df %>%
        dplyr::select(sample_id, phenotype) %>%
        tibble::column_to_rownames("sample_id")
      
      # ── 对齐样本 ─────────────────────────────────────
      common_samples <- intersect(rownames(meta), rownames(feature_table))
      meta          <- meta[common_samples,          , drop = FALSE]
      feature_table <- feature_table[common_samples, , drop = FALSE]
      
      if (nrow(feature_table) == 0 || ncol(feature_table) == 0) {
        cat("  过滤后矩阵为空，跳过\n"); next
      }
      
      cat("  [预处理后] feature_table dim:",
          paste(dim(feature_table), collapse = "x"), "\n")
      
      # ── 过滤有效分组 ──────────────────────────────────
      group_counts <- table(meta$phenotype)
      valid_groups <- names(group_counts[group_counts >= 2])
      if (length(valid_groups) < 2) next
      
      meta          <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
      feature_table <- feature_table[rownames(meta), , drop = FALSE]
      
      # ── 枚举所有比较对（3个pair）───────────────────────
      all_groups    <- sort(unique(as.character(meta$phenotype)))
      disease_pairs <- combn(all_groups, 2, simplify = FALSE)
      n_pairs_total <- length(disease_pairs)
      
      scenario_idx <- match(scenario, scenarios_to_run)
      rep_idx      <- as.integer(rep)
      
      # 用于在raw里也保存“本rep下的3个pair标签集合”
      comparison_set_label <- paste(
        vapply(disease_pairs, \(p) paste0(p[1], "_vs_", p[2]), ""),
        collapse = ";"
      )
      
      for (i in seq_along(disease_pairs)) {
        
        pair      <- disease_pairs[[i]]
        g1        <- pair[1]
        g2        <- pair[2]
        comp_name <- paste0(g1, "_vs_", g2)
        
        pair_seed <- BASE_SEED +
          1000000L * scenario_idx +
          10000L   * rep_idx +
          100L     * i
        
        cat("  比较:", comp_name, "\n"); flush.console()
        
        samples_all <- rownames(meta)[meta$phenotype %in% pair]
        if (length(samples_all) < 4) {
          cat("  样本数不足，跳过\n"); next
        }
        
        feat_all <- feature_table[samples_all, , drop = FALSE]
        feat_all <- feat_all[, colSums(feat_all) > 0, drop = FALSE]
        if (ncol(feat_all) == 0) next
        
        meta_all <- meta[samples_all, , drop = FALSE]
        meta_all$phenotype <- as.character(meta_all$phenotype)
        
        # ── Step 1：全量运行（该pair）────────────────────
        S_all <- tryCatch(
          run_maaslin2_core(
            feature_table = feat_all,
            meta_df       = meta_all,
            qval_cutoff   = QVAL_THR,
            top_k         = TOP_K,     # 每pair top-k
            debug         = FALSE
          ),
          error = function(e) {
            cat("  全量MaAsLin2失败:", e$message, "\n")
            character(0)
          }
        )
        cat("  全量特征数:", length(S_all), "\n"); flush.console()
        
        if (length(S_all) == 0) {
          jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
            method            = "MaAsLin2_naive",
            scenario          = scenario,
            rep               = rep_idx,
            comparison        = comp_name,
            comparison_set    = comparison_set_label,
            top_k_per_pair    = as.integer(TOP_K),
            n_pairs           = as.integer(n_pairs_total),
            top_k_total       = as.integer(TOP_K * n_pairs_total),  # 通常150
            jaccard_mean      = NA_real_,
            jaccard_median    = NA_real_,
            jaccard_sd        = NA_real_,
            jaccard_min       = NA_real_,
            jaccard_max       = NA_real_,
            n_boots           = as.integer(B),
            n_boots_valid     = 0L,
            n_features_full   = 0L
          )
          next
        }
        
        # ── Step 2：B 次 subsampling（pair内按phenotype分层）─
        jac_vals <- rep(NA_real_, B)
        
        for (b in seq_len(B)) {
          
          set.seed(pair_seed + b)
          
          sub_ids <- meta_all %>%
            tibble::rownames_to_column("sample_id") %>%
            dplyr::group_by(phenotype) %>%
            dplyr::slice_sample(prop = SUBSAMPLE) %>%
            dplyr::ungroup() %>%
            dplyr::pull(sample_id) %>%
            as.character()
          
          if (length(sub_ids) < 4) next
          
          feat_sub <- feat_all[sub_ids, , drop = FALSE]
          feat_sub <- feat_sub[, colSums(feat_sub) > 0, drop = FALSE]
          if (ncol(feat_sub) == 0) next
          
          meta_sub <- meta_all[sub_ids, , drop = FALSE]
          
          S_b <- tryCatch(
            run_maaslin2_core(
              feature_table = feat_sub,
              meta_df       = meta_sub,
              qval_cutoff   = QVAL_THR,
              top_k         = TOP_K,
              debug         = FALSE
            ),
            error = function(e) character(0)
          )
          
          jac_vals[b] <- jaccard(S_b, S_all)
        }
        
        jac_clean     <- jac_vals[!is.na(jac_vals)]
        n_boots_valid <- length(jac_clean)
        
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method            = "MaAsLin2_naive",
          scenario          = scenario,
          rep               = rep_idx,
          comparison        = comp_name,
          comparison_set    = comparison_set_label,
          top_k_per_pair    = as.integer(TOP_K),
          n_pairs           = as.integer(n_pairs_total),
          top_k_total       = as.integer(TOP_K * n_pairs_total),
          jaccard_mean      = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
          jaccard_median    = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
          jaccard_sd        = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
          jaccard_min       = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
          jaccard_max       = if (n_boots_valid == 0) NA_real_ else max(jac_clean),
          n_boots           = as.integer(B),
          n_boots_valid     = as.integer(n_boots_valid),
          n_features_full   = as.integer(length(S_all))
        )
      }
    }
  }
  
  # ── 汇总并保存 ────────────────────────────────────────────
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    empty <- empty_jaccard_tibble()
    saveRDS(empty, "benchmark/results/maaslin2_naive_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/maaslin2_naive_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/maaslin2_naive_jaccard_summary.rds")
    return(list(
      jaccard_raw = empty,
      jaccard_pair_summary = empty,
      jaccard_scenario_summary = empty
    ))
  }
  
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
  # ── 1) Pair 级汇总（保留 comparison）──────────────────────────
  jaccard_pair_summary <- jaccard_raw %>%
    dplyr::group_by(method, scenario, comparison, top_k_per_pair) %>%
    dplyr::summarise(
      jaccard_mean_across_reps   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_across_reps     = mean(jaccard_sd,     na.rm = TRUE),
      n_boots                    = max(n_boots,          na.rm = TRUE),
      n_boots_valid_mean         = mean(n_boots_valid,   na.rm = TRUE),
      n_pairs_mean               = mean(n_pairs,         na.rm = TRUE),
      top_k_total_mean           = mean(top_k_total,     na.rm = TRUE),
      .groups = "drop"
    )
  
  # ── 2) 场景级汇总：每个 rep 内先对 pair 取均值，再 across reps ──
  jaccard_rep_pairavg <- jaccard_raw %>%
    dplyr::group_by(method, scenario, rep, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set         = paste(sort(unique(comparison)), collapse = ";"),
      n_pairs                = dplyr::n_distinct(comparison),
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      top_k_total = top_k_per_pair * n_pairs  # ✅ 移至 mutate，避免 summarise 内长度歧义
    )
  
  jaccard_scenario_summary <- jaccard_rep_pairavg %>%
    dplyr::group_by(method, scenario, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set             = paste(sort(unique(comparison_set)), collapse = " | "),
      n_pairs_mean               = mean(n_pairs,     na.rm = TRUE),
      top_k_total_mean           = mean(top_k_total, na.rm = TRUE),
      jaccard_mean_across_reps   = mean(jaccard_mean_pairavg,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median_pairavg, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n── MaAsLin2 Naive Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  cat("\n── MaAsLin2 Naive Jaccard 稳定性汇总（场景级：3个pair取平均） ──\n")
  print(jaccard_scenario_summary)
  
  # ── 保存：raw + pair summary + scenario summary ───────────────
  saveRDS(jaccard_raw,              "benchmark/results/maaslin2_naive_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/maaslin2_naive_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/maaslin2_naive_jaccard_summary.rds")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# 执行
jaccard_maaslin2 <- run_maaslin2_jaccard()





# ══════════════════════════════════════════════════════════════
# Benchmark thread lock（必须在在可能触发并行计算之前设置）
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

# ══════════════════════════════════════════════════════════════
# 加载包
# ══════════════════════════════════════════════════════════════
library(Maaslin2)
library(tidyverse)

# ══════════════════════════════════════════════════════════════
# Thread Lock（library() 之后补充，控制已加载的 BLAS）
# ══════════════════════════════════════════════════════════════
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(NTHREAD)
  RhpcBLASctl::omp_set_num_threads(NTHREAD)
}

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(NTHREAD)
}

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

cat("Thread lock:", NTHREAD, "\n")
print(Sys.getenv(c("OMP_NUM_THREADS",
                   "OPENBLAS_NUM_THREADS",
                   "MKL_NUM_THREADS")))

# ══════════════════════════════════════════════════════════════
# 全局参数（只定义一次）
# ══════════════════════════════════════════════════════════════
scenarios_to_run <- c("B_batch_confound")  # ✅ 只跑场景B
reps <- as.character(1:3)

QVAL_THR  <- 0.05
LFC_THR   <- 1
TOP_K     <- 50
B         <- 30
SUBSAMPLE <- 0.8

# 确保NTHREAD存在
NTHREAD <- getOption("mc.cores", 1)

# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════
ensure_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

jaccard <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

empty_jaccard_tibble <- function() {
  tibble(
    method            = character(0),
    scenario          = character(0),
    rep               = integer(0),
    comparison        = character(0),
    comparison_set    = character(0),
    top_k_per_pair    = integer(0),
    n_pairs           = integer(0),
    top_k_total       = integer(0),
    jaccard_mean      = numeric(0),
    jaccard_median    = numeric(0),
    jaccard_sd        = numeric(0),
    jaccard_min       = numeric(0),
    jaccard_max       = numeric(0),
    n_boots           = integer(0),
    n_boots_valid     = integer(0),
    n_features_full   = integer(0)
  )
}

# ══════════════════════════════════════════════════════════════
# core：给定 feature_table + meta_df 跑 MaAsLin2（两两比较）
# ✅ 修改：fixed_effects 包含 batch，只返回 phenotype 效应的显著特征
# ══════════════════════════════════════════════════════════════
run_maaslin2_core <- function(feature_table,
                              meta_df,
                              qval_cutoff = 0.05,
                              top_k       = NULL,
                              debug       = FALSE) {
  
  if (debug) {
    cat("  [DEBUG] feature_table dim:",
        paste(dim(feature_table), collapse = "x"), "\n")
  }
  
  # ── 对齐样本 ────────────────────────────────────────────────
  common_samples <- intersect(rownames(feature_table), rownames(meta_df))
  if (length(common_samples) < 4) return(character(0))
  
  feature_table <- feature_table[common_samples, , drop = FALSE]
  meta_df       <- meta_df[common_samples,       , drop = FALSE]
  
  # ── 过滤有效分组 ───────────────────────────────────────────
  group_counts <- table(meta_df$phenotype)
  valid_groups <- names(group_counts[group_counts >= 2])
  if (length(valid_groups) < 2) return(character(0))
  
  meta_df       <- meta_df[meta_df$phenotype %in% valid_groups, , drop = FALSE]
  feature_table <- feature_table[rownames(meta_df), , drop = FALSE]
  
  # ── batch 列检查（必须有）──────────────────────────────────
  if (!("batch" %in% colnames(meta_df))) {
    stop("run_maaslin2_core: meta_df 缺少 batch 列，无法进行 batch-adjusted MaAsLin2")
  }
  
  disease_pairs <- combn(valid_groups, 2, simplify = FALSE)
  all_features  <- character(0)
  
  for (pair in disease_pairs) {
    
    pair_samples <- rownames(meta_df)[meta_df$phenotype %in% pair]
    if (length(pair_samples) < 4) next
    
    pair_meta    <- meta_df[pair_samples,    , drop = FALSE]
    pair_feature <- feature_table[pair_samples, , drop = FALSE]
    
    nonzero_taxa <- colSums(pair_feature) > 0
    pair_feature <- pair_feature[, nonzero_taxa, drop = FALSE]
    if (ncol(pair_feature) == 0) next
    
    # phenotype 参考组
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    
    # batch 作为固定效应：建议显式设为 factor
    pair_meta$batch <- factor(pair_meta$batch)
    
    tmp_dir <- tempfile(pattern = "maaslin2_core_")
    dir.create(tmp_dir, recursive = TRUE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    
    fit <- tryCatch({
      capture.output(
        fit_obj <- Maaslin2(
          input_data       = pair_feature,
          input_metadata   = pair_meta,
          output           = tmp_dir,
          fixed_effects    = c("phenotype", "batch"),  # ✅ batch-adjusted
          normalization    = "TSS",
          transform        = "LOG",
          analysis_method  = "LM",
          correction       = "BH",
          max_significance = 0.25,
          min_abundance    = 0.0,
          min_prevalence   = 0.0,
          min_variance     = 0.0,
          standardize      = FALSE,
          cores            = NTHREAD,
          plot_heatmap     = FALSE,
          plot_scatter     = FALSE
        ),
        type = "output"
      )
      fit_obj
    }, error = function(e) {
      if (debug) cat("  [core] MaAsLin2 失败:", e$message, "\n")
      NULL
    })
    
    if (is.null(fit)) next
    
    res_file <- file.path(tmp_dir, "all_results.tsv")
    if (!file.exists(res_file)) next
    
    res_df <- tryCatch(
      readr::read_tsv(res_file, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(res_df) || nrow(res_df) == 0) next
    
    # ✅ 关键修改：只取 metadata == "phenotype" 的特征（疾病效应）
    # 保留 batch 效应的结果在完整输出中，但Jaccard只评估phenotype特征
    feats <- res_df %>%
      dplyr::filter(metadata == "phenotype") %>%  # ✅ 只保留phenotype效应
      dplyr::filter(!is.na(qval), qval < qval_cutoff) %>%
      dplyr::arrange(desc(abs(coef)))
    
    if (!is.null(top_k)) feats <- dplyr::slice_head(feats, n = top_k)
    
    all_features <- c(all_features, feats$feature)
  }
  
  unique(all_features)
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 MaAsLin2（加入 batch 固定效应）
# 输出目录/文件名：全部改为 maaslin2_results_genus_adjusted / maaslin2_adjusted_*
# ══════════════════════════════════════════════════════════════
run_mock_maaslin2_genus <- function(scenario_full, rep_id) {
  
  # ✅ 避免rep间随机序列复制（使用42+rep_id，与lefse-stratified一致）
  set.seed(42 + as.integer(rep_id))
  
  cat("开始处理:", scenario_full, "rep", rep_id, "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  # 1. 读取数据 ────────────────────────────────────────────
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$maaslin2$abundance
  metadata_df   <- inputs$maaslin2$metadata
  
  cat("特征矩阵维度:", dim(counts_matrix), "\n")
  cat("样本数:", nrow(metadata_df), "\n")
  flush.console()
  
  # 2. 转置：行=样本，列=taxa
  counts_matrix <- ensure_matrix(counts_matrix)
  feature_table <- as.data.frame(t(counts_matrix))
  
  # 3. 整理 metadata：✅保留 batch
  if (!("batch" %in% colnames(metadata_df))) {
    stop("metadata_df 缺少 batch 列：无法进行 batch-adjusted 分析")
  }
  
  meta <- metadata_df %>%
    dplyr::select(sample_id, phenotype, batch) %>%
    tibble::column_to_rownames("sample_id")
  
  # 对齐样本
  common_samples <- intersect(rownames(meta), rownames(feature_table))
  meta          <- meta[common_samples, , drop = FALSE]
  feature_table <- feature_table[common_samples, , drop = FALSE]
  
  cat("对齐后样本数:", length(common_samples), "\n")
  cat("phenotype分布:\n"); print(table(meta$phenotype))
  cat("batch分布:\n");     print(table(meta$batch))
  cat("phenotype x batch:\n"); print(table(meta$phenotype, meta$batch))
  flush.console()
  
  # 4. 枚举所有两两比较对
  all_groups    <- sort(unique(as.character(meta$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n\n")
  flush.console()
  
  # ✅ 输出目录：adjusted
  base_result_dir <- file.path(
    "benchmark", "maaslin2", "maaslin2_results_genus_adjusted",
    paste0(scenario_full, "_rep", rep_id)
  )
  dir.create(base_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  compute_sec  <- 0
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]
    g2        <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    
    cat("比较", i, "/", length(disease_pairs), ":", comp_name, "\n")
    flush.console()
    
    pair_samples <- rownames(meta)[meta$phenotype %in% pair]
    if (length(pair_samples) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_meta    <- meta[pair_samples, , drop = FALSE]
    pair_feature <- feature_table[pair_samples, , drop = FALSE]
    
    # 去全0列
    nonzero_taxa <- colSums(pair_feature) > 0
    pair_feature <- pair_feature[, nonzero_taxa, drop = FALSE]
    
    cat("  样本数:", nrow(pair_feature),
        "| 特征数:", ncol(pair_feature), "\n")
    cat("  phenotype分布:\n"); print(table(pair_meta$phenotype))
    cat("  batch分布:\n");     print(table(pair_meta$batch))
    flush.console()
    
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    pair_meta$batch <- factor(pair_meta$batch)
    
    cat("  参考组:", ref_group, "\n")
    flush.console()
    
    comp_output_dir <- file.path(base_result_dir, comp_name)
    dir.create(comp_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ── 纯模型拟合计时
    t_compute_start <- proc.time()["elapsed"]
    
    maaslin_res <- tryCatch({
      capture.output(
        fit <- Maaslin2(
          input_data       = pair_feature,
          input_metadata   = pair_meta,
          output           = comp_output_dir,
          fixed_effects    = c("phenotype", "batch"),  # ✅ batch-adjusted
          normalization    = "TSS",
          transform        = "LOG",
          analysis_method  = "LM",
          correction       = "BH",
          max_significance = 0.25,
          min_abundance    = 0.0,
          min_prevalence   = 0.0,
          min_variance     = 0.0,
          standardize      = FALSE,
          plot_heatmap     = FALSE,
          plot_scatter     = FALSE
        ),
        type = "output"
      )
      fit
    }, error = function(e) {
      cat("  MaAsLin2 失败:", e$message, "\n"); NULL
    })
    
    t_compute_end <- proc.time()["elapsed"]
    compute_sec   <- compute_sec + as.numeric(t_compute_end - t_compute_start)
    
    if (is.null(maaslin_res)) next
    
    res_file <- file.path(comp_output_dir, "all_results.tsv")
    if (!file.exists(res_file)) {
      cat("  结果文件不存在，跳过\n"); next
    }
    
    res_df <- readr::read_tsv(res_file, show_col_types = FALSE) %>%
      dplyr::mutate(
        comparison = comp_name,
        group_a    = g1,
        group_b    = g2,
        ref_group  = ref_group,
        scenario   = scenario_full,
        replicate  = rep_id
      )
    
    # 统计：phenotype效应和batch效应分别有多少显著
    n_sig_pheno <- res_df %>%
      dplyr::filter(metadata == "phenotype", !is.na(qval), qval < QVAL_THR) %>%
      nrow()
    n_sig_batch <- res_df %>%
      dplyr::filter(metadata == "batch", !is.na(qval), qval < QVAL_THR) %>%
      nrow()
    
    cat("  结果行数:", nrow(res_df),
        "| phenotype显著:", n_sig_pheno,
        "| batch显著:", n_sig_batch, "\n")
    flush.console()
    
    readr::write_csv(
      res_df,
      file.path(base_result_dir, paste0(comp_name, "_results.csv"))
    )
    results_list[[comp_name]] <- res_df
  }
  
  if (length(results_list) > 0) {
    all_results_df <- dplyr::bind_rows(results_list)
    readr::write_csv(all_results_df,
                     file.path(base_result_dir, "all_comparisons_results.csv"))
    cat("结果已保存至:", base_result_dir, "\n")
  } else {
    cat("无有效比较结果\n")
  }
  
  t_wall_end <- proc.time()["elapsed"]
  wall_sec   <- as.numeric(t_wall_end - t_wall_start)
  
  cat("\n函数执行完成\n")
  cat("  总耗时   (wall_sec)   :", round(wall_sec,    2), "秒\n")
  cat("  纯模型   (compute_sec):", round(compute_sec, 2), "秒\n")
  cat("  共完成", length(results_list), "个比较\n")
  flush.console()
  
  list(
    results     = results_list,
    wall_sec    = wall_sec,
    compute_sec = compute_sec
  )
}

# ══════════════════════════════════════════════════════════════
# 只跑场景B + 3个rep，并保存 timing（全部改为 adjusted 后缀）
# ══════════════════════════════════════════════════════════════
run_all <- function() {
  
  cat("\n=============================\n")
  cat("开始运行 MaAsLin2 Adjusted (phenotype + batch)\n")
  cat("=============================\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n--- 场景:", scenario, "Rep:", rep, "---\n")
      flush.console()
      
      key <- paste0(scenario, "_rep", rep)
      
      out <- tryCatch(
        run_mock_maaslin2_genus(scenario, rep),
        error = function(e) { cat("失败:", e$message, "\n"); NULL }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out$results
        timing_records[[length(timing_records) + 1]] <- tibble::tibble(
          method      = "MaAsLin2_adjusted",
          scenario    = scenario,
          rep         = as.integer(rep),
          wall_sec    = out$wall_sec,
          compute_sec = out$compute_sec
        )
      } else {
        results_all[[key]] <- NULL
      }
    }
  }
  
  cat("\nMaAsLin2 Adjusted 完成，成功:",
      sum(!sapply(results_all, is.null)), "/", length(results_all), "\n")
  
  timing_df <- dplyr::bind_rows(timing_records)
  
  timing_summary <- timing_df %>%
    dplyr::group_by(method, scenario) %>%
    dplyr::summarise(
      wall_mean_sec    = mean(wall_sec,      na.rm = TRUE),
      wall_median_sec  = median(wall_sec,    na.rm = TRUE),
      wall_min_sec     = min(wall_sec,       na.rm = TRUE),
      wall_max_sec     = max(wall_sec,       na.rm = TRUE),
      compute_mean_sec = mean(compute_sec,   na.rm = TRUE),
      compute_med_sec  = median(compute_sec, na.rm = TRUE),
      compute_min_sec  = min(compute_sec,    na.rm = TRUE),
      compute_max_sec  = max(compute_sec,    na.rm = TRUE),
      n_runs           = dplyr::n(),
      .groups = "drop"
    )
  
  cat("\n── MaAsLin2 Adjusted 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(timing_df,      "benchmark/results/maaslin2_adjusted_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/maaslin2_adjusted_timing_summary.rds")
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV（MaAsLin2 Adjusted 版本）
# ══════════════════════════════════════════════════════════════
collect_maaslin2_adjusted_results <- function(
    base_dir  = "benchmark/maaslin2/maaslin2_results_genus_adjusted",
    scenarios, reps) {
  
  all_rows <- list()
  
  for (scenario in scenarios) {
    for (rep in reps) {
      result_dir <- file.path(base_dir,
                              paste0(scenario, "_rep", rep))
      if (!dir.exists(result_dir)) next
      
      csv_files <- list.files(result_dir,
                              pattern    = "_results\\.csv$",
                              full.names = TRUE)
      csv_files <- csv_files[
        !grepl("all_comparisons_results\\.csv$", csv_files)
      ]
      
      for (f in csv_files) {
        df <- tryCatch(
          read_csv(f, show_col_types = FALSE),
          error = function(e) NULL
        )
        if (!is.null(df) && nrow(df) > 0)
          all_rows[[length(all_rows) + 1]] <- df
      }
    }
  }
  
  if (length(all_rows) == 0) { cat("未找到任何结果文件\n"); return(NULL) }
  bind_rows(all_rows)
}

# 注意：scenarios_to_run 应该只包含 "B_batch_confound"
maaslin2_adjusted_summary_raw <- collect_maaslin2_adjusted_results(
  scenarios = scenarios_to_run,  # c("B_batch_confound")
  reps      = reps
)

if (!is.null(maaslin2_adjusted_summary_raw)) {
  cat("MaAsLin2 Adjusted 汇总行数:", nrow(maaslin2_adjusted_summary_raw), "\n")
  print(head(maaslin2_adjusted_summary_raw))
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(maaslin2_adjusted_summary_raw, 
          "benchmark/results/maaslin2_adjusted_summary_raw.rds")
  write_csv(maaslin2_adjusted_summary_raw,
            "benchmark/results/maaslin2_adjusted_summary_raw.csv")
  cat("\nMaAsLin2 Adjusted 汇总结果已保存\n")
}

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估（Adjusted版）
# 参考 lefse_stratified_jaccard 的设计
# ══════════════════════════════════════════════════════════════
run_maaslin2_adjusted_jaccard <- function() {
  
  BASE_SEED <- 42L
  SCENARIO  <- "B_batch_confound"
  
  jaccard_records <- list()
  
  cat("\n", strrep("=", 62), "\n")
  cat("MaAsLin2 Adjusted Jaccard 稳定性评估 — 场景 B（pairwise + pair均值汇总）\n")
  cat(strrep("=", 62), "\n")
  
  # ── 内部小工具：对“一个pair的数据”跑一次 batch-adjusted MaAsLin2 并取特征 ──
  run_one_pair_adjusted <- function(pair_feature, pair_meta, pair_levels, top_k) {
    
    # 去全0列
    pair_feature <- pair_feature[, colSums(pair_feature) > 0, drop = FALSE]
    if (ncol(pair_feature) == 0) return(character(0))
    
    # 设置参考组（两组比较只会有一个 phenotype 系数）
    ref_group <- if ("Healthy" %in% pair_levels) "Healthy" else sort(pair_levels)[1]
    
    pair_meta <- pair_meta[, c("phenotype", "batch"), drop = FALSE]
    pair_meta$phenotype <- stats::relevel(factor(pair_meta$phenotype, levels = pair_levels), ref = ref_group)
    pair_meta$batch     <- factor(pair_meta$batch)
    
    tmp_dir <- tempfile(pattern = "maaslin2_jaccard_adj_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    
    fit <- tryCatch({
      capture.output(
        fit_obj <- Maaslin2::Maaslin2(
          input_data       = pair_feature,
          input_metadata   = pair_meta,
          output           = tmp_dir,
          fixed_effects    = c("phenotype", "batch"),
          normalization    = "TSS",
          transform        = "LOG",
          analysis_method  = "LM",
          correction       = "BH",
          max_significance = 0.25,
          min_abundance    = 0.0,
          min_prevalence   = 0.0,
          min_variance     = 0.0,
          standardize      = FALSE,
          cores            = NTHREAD,
          plot_heatmap     = FALSE,
          plot_scatter     = FALSE
        ),
        type = "output"
      )
      fit_obj
    }, error = function(e) NULL)
    
    if (is.null(fit)) return(character(0))
    
    res_file <- file.path(tmp_dir, "all_results.tsv")
    if (!file.exists(res_file)) return(character(0))
    
    res_df <- tryCatch(
      readr::read_tsv(res_file, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
    
    feats <- res_df %>%
      dplyr::filter(metadata == "phenotype", !is.na(qval), qval < QVAL_THR) %>%
      dplyr::arrange(dplyr::desc(abs(coef)))
    
    if (!is.null(top_k) && top_k > 0) {
      feats <- dplyr::slice_head(feats, n = top_k)
    }
    
    unique(feats$feature)
  }
  
  for (rep in reps) {
    
    cat("\n[Jaccard-Adjusted] Rep:", rep, "\n"); flush.console()
    
    inputs <- tryCatch(
      generate_method_inputs(SCENARIO, rep),
      error = function(e) {
        cat("数据读取失败:", e$message, "\n"); NULL
      }
    )
    if (is.null(inputs)) next
    
    counts_matrix <- inputs$maaslin2$abundance
    metadata_df   <- inputs$maaslin2$metadata
    
    # ── 硬防护 + 基本清洗 ─────────────────────────────────────
    counts_matrix <- ensure_matrix(counts_matrix)
    counts_matrix <- round(pmax(counts_matrix, 0))
    
    col_keep      <- colSums(counts_matrix) > 0
    counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
    
    row_keep      <- rowSums(counts_matrix) > 0
    counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
    
    feature_table <- as.data.frame(t(counts_matrix))  # samples × taxa
    
    if (!("batch" %in% colnames(metadata_df))) {
      cat("  metadata 缺少 batch 列，跳过\n"); next
    }
    
    meta <- metadata_df %>%
      dplyr::select(sample_id, phenotype, batch) %>%
      tibble::column_to_rownames("sample_id")
    
    # ── 对齐样本 ─────────────────────────────────────────────
    common_samples <- intersect(rownames(meta), rownames(feature_table))
    meta          <- meta[common_samples,          , drop = FALSE]
    feature_table <- feature_table[common_samples, , drop = FALSE]
    
    if (nrow(feature_table) == 0 || ncol(feature_table) == 0) {
      cat("  过滤后矩阵为空，跳过\n"); next
    }
    
    cat("  [预处理后] feature_table dim:",
        paste(dim(feature_table), collapse = "x"), "\n")
    
    # ── 过滤有效分组 ─────────────────────────────────────────
    group_counts <- table(meta$phenotype)
    valid_groups <- names(group_counts[group_counts >= 2])
    if (length(valid_groups) < 2) next
    
    meta          <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
    feature_table <- feature_table[rownames(meta), , drop = FALSE]
    
    cat("  batch分布:\n")
    print(table(meta$phenotype, meta$batch, dnn = c("phenotype", "batch")))
    flush.console()
    
    # ── 枚举 pairs（3个两两比较）─────────────────────────────
    all_groups    <- sort(unique(as.character(meta$phenotype)))
    disease_pairs <- combn(all_groups, 2, simplify = FALSE)
    n_pairs_total <- length(disease_pairs)
    
    comparison_set_label <- paste(
      vapply(disease_pairs, \(p) paste0(p[1], "_vs_", p[2]), ""),
      collapse = ";"
    )
    
    rep_idx <- as.integer(rep)
    
    for (i in seq_along(disease_pairs)) {
      
      pair      <- disease_pairs[[i]]
      g1        <- pair[1]
      g2        <- pair[2]
      comp_name <- paste0(g1, "_vs_", g2)
      
      # seed：rep + pair + boot
      scenario_idx <- 1  # 场景 B 的索引
      pair_seed <- BASE_SEED + 1000000L*scenario_idx + 10000L * rep_idx + 100L * i
      
      cat("  比较:", comp_name, "\n"); flush.console()
      
      samples_all <- rownames(meta)[meta$phenotype %in% pair]
      if (length(samples_all) < 4) {
        cat("    样本数不足，跳过\n"); next
      }
      
      feat_all <- feature_table[samples_all, , drop = FALSE]
      meta_all <- meta[samples_all, , drop = FALSE]
      
      # ── ① 全量（该pair）→ S_all ───────────────────────────
      cat("    计算该pair全量参考集 S_all...\n"); flush.console()
      S_all <- tryCatch(
        run_one_pair_adjusted(
          pair_feature = feat_all,
          pair_meta    = meta_all,
          pair_levels  = pair,
          top_k        = TOP_K
        ),
        error = function(e) character(0)
      )
      cat("    S_all 特征数:", length(S_all), "\n"); flush.console()
      
      if (length(S_all) == 0) {
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method            = "MaAsLin2_adjusted",
          scenario          = SCENARIO,
          rep               = rep_idx,
          comparison        = comp_name,
          comparison_set    = comparison_set_label,
          top_k_per_pair    = as.integer(TOP_K),
          n_pairs           = as.integer(n_pairs_total),
          top_k_total       = as.integer(TOP_K * n_pairs_total),  # 通常150
          jaccard_mean      = NA_real_,
          jaccard_median    = NA_real_,
          jaccard_sd        = NA_real_,
          jaccard_min       = NA_real_,
          jaccard_max       = NA_real_,
          n_boots           = as.integer(B),
          n_boots_valid     = 0L,
          n_features_full   = 0L
        )
        next
      }
      
      # ── ② Bootstrap：pair内 batch × phenotype 分层抽样 ─────
      jac_vals <- rep(NA_real_, B)
      
      # 为分层抽样准备 rownames
      meta_all2 <- meta_all
      rownames(meta_all2) <- samples_all
      
      for (b in seq_len(B)) {
        
        set.seed(pair_seed + b)
        
        sub_ids <- meta_all2 %>%
          tibble::rownames_to_column("sample_id") %>%
          dplyr::group_by(batch, phenotype) %>%
          dplyr::slice_sample(prop = SUBSAMPLE) %>%
          dplyr::ungroup() %>%
          dplyr::pull(sample_id) %>%
          as.character()
        
        if (length(sub_ids) < 4) next
        
        feat_sub <- feat_all[sub_ids, , drop = FALSE]
        feat_sub <- feat_sub[, colSums(feat_sub) > 0, drop = FALSE]
        if (ncol(feat_sub) == 0) next
        
        meta_sub <- meta_all2[sub_ids, , drop = FALSE]
        
        # 基本可估性检查：两组都在 + 至少两个batch（更稳一些）
        if (dplyr::n_distinct(meta_sub$phenotype) < 2) next
        if (dplyr::n_distinct(meta_sub$batch) < 2) next
        
        S_b <- tryCatch(
          run_one_pair_adjusted(
            pair_feature = feat_sub,
            pair_meta    = meta_sub,
            pair_levels  = pair,
            top_k        = TOP_K
          ),
          error = function(e) character(0)
        )
        
        jac_vals[b] <- jaccard(S_b, S_all)
        
        if (b %% 10 == 0) {
          cat(sprintf("    boot %d/%d | S_b=%d | jac=%.3f\n",
                      b, B, length(S_b),
                      ifelse(is.na(jac_vals[b]), NA, jac_vals[b])))
          flush.console()
        }
      }
      
      jac_clean     <- jac_vals[!is.na(jac_vals)]
      n_boots_valid <- length(jac_clean)
      
      jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
        method            = "MaAsLin2_adjusted",
        scenario          = SCENARIO,
        rep               = rep_idx,
        comparison        = comp_name,
        comparison_set    = comparison_set_label,
        top_k_per_pair    = as.integer(TOP_K),
        n_pairs           = as.integer(n_pairs_total),
        top_k_total       = as.integer(TOP_K * n_pairs_total),
        jaccard_mean      = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
        jaccard_median    = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
        jaccard_sd        = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
        jaccard_min       = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
        jaccard_max       = if (n_boots_valid == 0) NA_real_ else max(jac_clean),
        n_boots           = as.integer(B),
        n_boots_valid     = as.integer(n_boots_valid),
        n_features_full   = as.integer(length(S_all))
      )
      
      cat(sprintf(
        "    [%s | %s | rep%s] Jaccard均值=%.3f (SD=%.3f, 有效boots=%d/%d, S_all=%d)\n",
        SCENARIO, comp_name, rep_idx,
        if (n_boots_valid == 0) NA else mean(jac_clean),
        if (n_boots_valid == 0) NA else sd(jac_clean),
        n_boots_valid, B, length(S_all)
      ))
      flush.console()
    }
  }
  
  # ── 汇总 & 保存 ────────────────────────────────────────────
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    cat("\n警告：所有 rep 均无有效记录\n")
    empty <- empty_jaccard_tibble()
    saveRDS(empty, "benchmark/results/maaslin2_adjusted_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/maaslin2_adjusted_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/maaslin2_adjusted_jaccard_summary.rds")
    return(list(
      jaccard_raw = empty,
      jaccard_pair_summary = empty,
      jaccard_scenario_summary = empty
    ))
  }
  
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
  # 1) Pair级汇总（保留 comparison）
  jaccard_pair_summary <- jaccard_raw %>%
    dplyr::group_by(method, scenario, comparison, top_k_per_pair) %>%
    dplyr::summarise(
      jaccard_mean_across_reps   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_across_reps     = mean(jaccard_sd,     na.rm = TRUE),
      n_boots                    = max(n_boots, na.rm = TRUE),
      n_boots_valid_mean         = mean(n_boots_valid,  na.rm = TRUE),
      n_pairs_mean               = mean(n_pairs,        na.rm = TRUE),
      top_k_total_mean           = mean(top_k_total,    na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2) 场景级汇总：每个rep内先对3个pair取均值，再 across reps
  jaccard_rep_pairavg <- jaccard_raw %>%
    dplyr::group_by(method, scenario, rep, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set = paste(sort(unique(comparison)), collapse = ";"),
      n_pairs        = dplyr::n_distinct(comparison),
      top_k_total    = top_k_per_pair * n_pairs,
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      .groups = "drop"
    )%>%
    dplyr::mutate(
      top_k_total = top_k_per_pair * n_pairs  # ✅ 移至 mutate，避免 summarise 内长度歧义
    )
  
  jaccard_scenario_summary <- jaccard_rep_pairavg %>%
    dplyr::group_by(method, scenario, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set      = paste(sort(unique(comparison_set)), collapse = " | "),
      n_pairs_mean        = mean(n_pairs, na.rm = TRUE),
      top_k_total_mean    = mean(top_k_total, na.rm = TRUE),
      jaccard_mean_across_reps   = mean(jaccard_mean_pairavg,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median_pairavg, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n── MaAsLin2 Adjusted Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  
  cat("\n── MaAsLin2 Adjusted Jaccard 稳定性汇总（场景级：3个pair取平均） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,              "benchmark/results/maaslin2_adjusted_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/maaslin2_adjusted_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/maaslin2_adjusted_jaccard_summary.rds")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# 执行入口
out_maaslin2_adjusted <- run_all()
jaccard_maaslin2_adjusted <- run_maaslin2_adjusted_jaccard()