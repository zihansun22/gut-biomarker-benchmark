
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

library(DESeq2)
library(tidyverse)

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
LFC_THR   <- 1
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
# 核心函数：针对单个 pair 运行 DESeq2，返回显著特征名
# ══════════════════════════════════════════════════════════════
run_deseq2_core <- function(counts_matrix,
                            metadata_df,
                            test_group,
                            ref_group,
                            use_cooks        = TRUE,
                            qval_cutoff      = 0.05,
                            lfc_cutoff       = 1,
                            top_k            = NULL,
                            debug            = FALSE) {
  
  # ── 硬防护 ────────────────────────────────────────────────
  counts_matrix <- ensure_matrix(counts_matrix)
  
  if (debug) {
    cat("  [DEBUG] dim:", paste(dim(counts_matrix), collapse = "x"), "\n")
  }
  
  # ── 标量保护 ──────────────────────────────────────────────
  test_group <- as.character(test_group)[1]
  ref_group  <- as.character(ref_group)[1]
  
  # ── 确保含 sample_id 列 ───────────────────────────────────
  if (!"sample_id" %in% colnames(metadata_df)) {
    metadata_df <- tibble::rownames_to_column(
      as.data.frame(metadata_df), "sample_id"
    )
  }
  rownames(metadata_df) <- NULL
  
  # ── 只保留该 pair 的样本 ───────────────────────────────────
  pair        <- c(test_group, ref_group)
  metadata_df <- metadata_df[metadata_df$phenotype %in% pair, , drop = FALSE]
  rownames(metadata_df) <- NULL
  
  # ── 先对齐样本，再做所有过滤 ──────────────────────────────
  common_samples <- intersect(colnames(counts_matrix),
                              metadata_df$sample_id)
  if (length(common_samples) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  metadata_df <- metadata_df[
    match(common_samples, metadata_df$sample_id), , drop = FALSE
  ]
  rownames(metadata_df) <- NULL
  
  # ── 基本过滤（对齐后） ────────────────────────────────────
  counts_matrix <- round(pmax(counts_matrix, 0))
  
  # 去全 0 行
  row_keep <- rowSums(counts_matrix) > 0
  counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  # 与主分析一致的 keep 过滤（出现率 ≥ 20% 且均值 > 1）
  keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
    rowMeans(counts_matrix) > 1
  counts_matrix <- counts_matrix[keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  # ── 验证两组各有足够样本 ───────────────────────────────────
  group_counts <- table(metadata_df$phenotype)
  if (!all(pair %in% names(group_counts))) return(character(0))
  if (any(group_counts[pair] < 2))         return(character(0))
  
  # ── 构建 colData（Naive: 只保留 phenotype）────────────────
  col_data           <- as.data.frame(metadata_df)
  rownames(col_data) <- col_data$sample_id
  col_data$sample_id <- NULL
  col_data$phenotype <- factor(col_data$phenotype,
                               levels = c(ref_group, test_group))
  
  # ── design（Naive: 只用 phenotype）────────────────────────
  design_formula <- ~ phenotype
  
  # ── 构建并运行 DESeq2 ─────────────────────────────────────
  dds <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = counts_matrix,
      colData   = col_data,
      design    = design_formula
    )
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- estimateDispersions(dds, fitType = "local")
    dds <- nbinomWaldTest(dds)
    dds
  }, error = function(e) NULL)
  
  if (is.null(dds)) return(character(0))
  
  # ── 提取结果 ──────────────────────────────────────────────
  res <- tryCatch(
    results(
      dds,
      contrast      = c("phenotype", test_group, ref_group),
      alpha         = qval_cutoff,
      pAdjustMethod = "BH",
      cooksCutoff   = use_cooks
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(character(0))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("taxon") %>%
    filter(
      !is.na(padj),
      padj < qval_cutoff,
      abs(log2FoldChange) >= lfc_cutoff
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  if (nrow(res_df) == 0) return(character(0))
  if (!is.null(top_k)) res_df <- slice_head(res_df, n = top_k)
  
  res_df$taxon
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 DESeq2（Naive版本）
# ══════════════════════════════════════════════════════════════
run_mock_deseq2_genus <- function(scenario_full, rep_id, use_cooks = TRUE) {
  
  # ✅ 主函数 seed 固定为 42（与 lefse/maaslin naive 一致）
  set.seed(42)
  
  setting_tag <- if (use_cooks) "default" else "no_cooks"
  
  cat("开始处理:", scenario_full, "rep", rep_id, "[", setting_tag, "]\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  # 1. 读取数据 ────────────────────────────────────────────
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$deseq2$counts
  metadata_df   <- inputs$deseq2$groups
  
  cat("特征矩阵维度:", dim(counts_matrix), "\n")
  cat("样本数:", nrow(metadata_df), "\n")
  flush.console()
  
  # 2. 数据清洗 ────────────────────────────────────────────
  counts_matrix <- ensure_matrix(counts_matrix)
  counts_matrix <- round(pmax(counts_matrix, 0))
  
  # 去全 0 样本
  sample_sums <- colSums(counts_matrix)
  if (any(sample_sums == 0)) {
    cat("警告: 发现", sum(sample_sums == 0), "个全 0 样本，已删除\n")
    counts_matrix <- counts_matrix[, sample_sums > 0, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
  }
  
  # 去全 0 特征
  feature_sums <- rowSums(counts_matrix)
  if (any(feature_sums == 0)) {
    cat("警告: 发现", sum(feature_sums == 0), "个全 0 特征，已删除\n")
    counts_matrix <- counts_matrix[feature_sums > 0, , drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
  }
  
  # ✅ Naive: 只保留 sample_id 和 phenotype，丢弃 batch
  meta <- metadata_df %>%
    dplyr::select(sample_id, phenotype) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      phenotype = as.character(phenotype)
    )
  
  # 先对齐样本
  common_samples <- intersect(colnames(counts_matrix), meta$sample_id)
  counts_matrix  <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix  <- ensure_matrix(counts_matrix)
  meta           <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
  
  keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
    rowMeans(counts_matrix) > 1
  cat("低丰度过滤: 保留", sum(keep), "/", nrow(counts_matrix), "个特征\n")
  counts_matrix <- counts_matrix[keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  cat("对齐后样本数:", ncol(counts_matrix), "\n")
  cat("对齐后特征数:", nrow(counts_matrix), "\n")
  cat("分组分布:\n")
  print(table(meta$phenotype))
  flush.console()
  
  # 3. 枚举所有两两比较对 ──────────────────────────────────
  all_groups    <- sort(unique(as.character(meta$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n")
  cat("分组顺序:", paste(all_groups, collapse = ", "), "\n\n")
  flush.console()
  
  # ✅ 输出目录：deseq2_results_genus_naive_{default/no_cooks}
  base_result_dir <- file.path(
    "benchmark", "deseq2",
    paste0("deseq2_results_genus_naive_", setting_tag),
    paste0(scenario_full, "_rep", rep_id)
  )
  dir.create(base_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  compute_sec  <- 0
  
  for (i in seq_along(disease_pairs)) {
    
    pair       <- disease_pairs[[i]]
    ref_group  <- if ("Healthy" %in% pair) "Healthy" else pair[1]
    ref_group  <- as.character(ref_group)[1]
    test_group <- as.character(setdiff(pair, ref_group))[1]
    comp_name  <- paste0(test_group, "_vs_", ref_group)
    pair_id    <- paste(sort(pair), collapse = "_")
    
    cat("比较", i, "/", length(disease_pairs), ":", comp_name, "\n")
    cat("  配对:", pair_id, "| 参考组:", ref_group,
        "| 测试组:", test_group, "\n")
    flush.console()
    
    # 子集数据
    pair_idx <- meta$phenotype %in% pair
    if (sum(pair_idx) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_counts <- counts_matrix[, pair_idx, drop = FALSE]
    pair_counts <- ensure_matrix(pair_counts)
    pair_meta   <- meta[pair_idx, , drop = FALSE]
    rownames(pair_meta) <- NULL
    
    # 去全 0 行
    nonzero_taxa <- rowSums(pair_counts) > 0
    pair_counts  <- pair_counts[nonzero_taxa, , drop = FALSE]
    pair_counts  <- ensure_matrix(pair_counts)
    
    cat("  样本数:", ncol(pair_counts),
        "| 特征数:", nrow(pair_counts), "\n")
    cat("  分组分布:\n")
    print(table(pair_meta$phenotype))
    flush.console()
    
    # ── 纯模型拟合计时 ────────────────────────────────────
    t_compute_start <- proc.time()["elapsed"]
    
    # 构建 DESeq2 对象
    col_data           <- as.data.frame(pair_meta)
    rownames(col_data) <- col_data$sample_id
    col_data$sample_id <- NULL
    col_data$phenotype <- factor(col_data$phenotype,
                                 levels = c(ref_group, test_group))
    
    dds <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = pair_counts,
        colData   = col_data,
        design    = ~ phenotype
      )
      dds <- estimateSizeFactors(dds, type = "poscounts")
      dds <- estimateDispersions(dds, fitType = "local")
      dds <- nbinomWaldTest(dds)
      dds
    }, error = function(e) {
      cat("  DESeq2 运行失败:", e$message, "\n"); NULL
    })
    
    t_compute_end <- proc.time()["elapsed"]
    compute_sec   <- compute_sec + as.numeric(t_compute_end - t_compute_start)
    
    if (is.null(dds)) next
    
    res <- tryCatch({
      results(
        dds,
        contrast      = c("phenotype", test_group, ref_group),
        alpha         = 0.05,
        pAdjustMethod = "BH",
        cooksCutoff   = use_cooks
      )
    }, error = function(e) {
      cat("  提取结果失败:", e$message, "\n"); NULL
    })
    
    if (is.null(res)) next
    
    result_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("taxon") %>%
      dplyr::rename(
        base_mean = baseMean,
        log2fc    = log2FoldChange,
        lfc_se    = lfcSE,
        stat      = stat,
        p_val     = pvalue,
        q_val     = padj
      ) %>%
      dplyr::mutate(
        sig_fdr      = !is.na(q_val) & q_val < 0.05,
        sig_fdr_eff  = sig_fdr & abs(log2fc) >= 1,
        comparison   = comp_name,
        pair_id      = pair_id,
        test_group   = test_group,
        ref_group    = ref_group,
        scenario     = scenario_full,
        replicate    = rep_id,
        setting      = setting_tag
      ) %>%
      dplyr::select(
        taxon, base_mean, log2fc, lfc_se, stat, p_val, q_val,
        sig_fdr, sig_fdr_eff,
        comparison, pair_id, test_group, ref_group,
        scenario, replicate, setting
      )
    
    n_p_sig <- sum(result_df$p_val < 0.05, na.rm = TRUE)
    n_q_sig <- sum(result_df$sig_fdr,      na.rm = TRUE)
    n_q_eff <- sum(result_df$sig_fdr_eff,  na.rm = TRUE)
    
    cat("  发现", nrow(result_df), "个特征结果\n")
    cat("  p < 0.05:", n_p_sig, "个\n")
    cat("  q < 0.05 (FDR):", n_q_sig, "个\n")
    cat("  q < 0.05 + |log2FC| >= 1:", n_q_eff, "个\n")
    flush.console()
    
    # I/O（不计入 compute_sec）
    comp_output_dir <- file.path(base_result_dir, pair_id)
    dir.create(comp_output_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv(result_df, file.path(comp_output_dir, "results.csv"))
    
    results_list[[comp_name]] <- result_df
  }
  
  # 合并所有比较对结果
  if (length(results_list) > 0) {
    all_results_df <- bind_rows(results_list)
    write_csv(all_results_df,
              file.path(base_result_dir, "all_comparisons_results.csv"))
    cat("结果已保存至:", base_result_dir, "\n")
    cat("\n=== 汇总统计 ===\n")
    cat("sig_fdr (q < 0.05):",
        sum(all_results_df$sig_fdr, na.rm = TRUE), "个\n")
    cat("sig_fdr_eff (q < 0.05 + |log2FC| >= 1):",
        sum(all_results_df$sig_fdr_eff, na.rm = TRUE), "个\n")
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
# 批量运行 + timing（Naive版本）
# ══════════════════════════════════════════════════════════════
run_all_naive <- function(use_cooks) {
  
  tag <- if (use_cooks) "default" else "no_cooks"
  cat("\n=============================\n")
  cat("开始运行 DESeq2 Naive [", tag, "]\n")
  cat("=============================\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n--- 场景:", scenario, "Rep:", rep, "---\n")
      flush.console()
      
      key <- paste0(scenario, "_rep", rep)
      
      out <- tryCatch(
        run_mock_deseq2_genus(scenario, rep, use_cooks = use_cooks),
        error = function(e) { cat("失败:", e$message, "\n"); NULL }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out$results
        timing_records[[length(timing_records) + 1]] <- tibble(
          method      = paste0("DESeq2_naive_", tag),
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
  
  cat("\nDESeq2 Naive [", tag, "] 完成，成功:",
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
  
  cat("\n── DESeq2 Naive [", tag, "] 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(timing_df,
          paste0("benchmark/results/deseq2_naive_timing_raw_",     tag, ".rds"))
  saveRDS(timing_summary,
          paste0("benchmark/results/deseq2_naive_timing_summary_", tag, ".rds"))
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

out_naive_default  <- run_all_naive(use_cooks = TRUE)
out_naive_no_cooks <- run_all_naive(use_cooks = FALSE)

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV（Naive版本）
# ══════════════════════════════════════════════════════════════
collect_deseq2_results_naive <- function(setting_tag = "default",
                                         scenarios, reps) {
  base_dir <- file.path(
    "benchmark", "deseq2",
    paste0("deseq2_results_genus_naive_", setting_tag)
  )
  
  all_rows <- list()
  for (scenario in scenarios) {
    for (rep in reps) {
      result_dir <- file.path(base_dir,
                              paste0(scenario, "_rep", rep))
      if (!dir.exists(result_dir)) next
      
      for (comp_dir in list.dirs(result_dir, recursive = FALSE)) {
        res_file <- file.path(comp_dir, "results.csv")
        if (!file.exists(res_file)) next
        df <- tryCatch(
          read_csv(res_file, show_col_types = FALSE),
          error = function(e) NULL
        )
        if (!is.null(df) && nrow(df) > 0)
          all_rows[[length(all_rows) + 1]] <- df
      }
    }
  }
  
  if (length(all_rows) == 0) { cat("未找到结果\n"); return(NULL) }
  bind_rows(all_rows)
}

summary_naive_default  <- collect_deseq2_results_naive("default",  scenarios_to_run, reps)
summary_naive_no_cooks <- collect_deseq2_results_naive("no_cooks", scenarios_to_run, reps)

dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(summary_naive_default,  "benchmark/results/deseq2_naive_summary_default.rds")
saveRDS(summary_naive_no_cooks, "benchmark/results/deseq2_naive_summary_no_cooks.rds")
write_csv(summary_naive_default,  "benchmark/results/deseq2_naive_summary_default.csv")
write_csv(summary_naive_no_cooks, "benchmark/results/deseq2_naive_summary_no_cooks.csv")

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估（Naive版本）
# ══════════════════════════════════════════════════════════════
run_deseq2_naive_jaccard <- function(use_cooks = TRUE) {
  
  setting_tag     <- if (use_cooks) "default" else "no_cooks"
  jaccard_records <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("DESeq2 Naive Jaccard 稳定性评估 [", setting_tag, "]\n")
  cat("══════════════════════════════════════════════════════════════\n")
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n[Jaccard] 场景:", scenario, "Rep:", rep,
          "[", setting_tag, "]\n")
      flush.console()
      
      inputs <- tryCatch(
        generate_method_inputs(scenario, rep),
        error = function(e) {
          cat("数据读取失败:", e$message, "\n"); NULL
        }
      )
      if (is.null(inputs)) next
      
      counts_matrix <- inputs$deseq2$counts
      metadata_df   <- inputs$deseq2$groups
      
      # ── 硬防护 + 基本清洗 ─────────────────────────────
      counts_matrix <- ensure_matrix(counts_matrix)
      counts_matrix <- round(pmax(counts_matrix, 0))
      
      # 去全 0 列
      col_keep      <- colSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # 去全 0 行
      row_keep      <- rowSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # ✅ Naive: metadata 只保留 phenotype
      meta <- metadata_df %>%
        dplyr::select(sample_id, phenotype) %>%
        dplyr::mutate(
          sample_id = as.character(sample_id),
          phenotype = as.character(phenotype)
        )
      
      # 对齐样本
      common <- intersect(colnames(counts_matrix), meta$sample_id)
      counts_matrix <- counts_matrix[, common, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      meta <- meta[match(common, meta$sample_id), , drop = FALSE]
      
      # keep 过滤（与主分析保持一致）
      keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
        rowMeans(counts_matrix) > 1
      counts_matrix <- counts_matrix[keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
        cat("  过滤后矩阵为空，跳过\n"); next
      }
      
      cat("  [预处理后] dim:",
          paste(dim(counts_matrix), collapse = "x"), "\n")
      
      # ── 过滤有效分组 ──────────────────────────────────
      group_counts <- table(meta$phenotype)
      valid_groups <- names(group_counts[group_counts >= 2])
      if (length(valid_groups) < 2) next
      
      meta <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
      counts_matrix <- counts_matrix[, meta$sample_id, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # ── 枚举所有比较对 ────────────────────────────────
      all_groups    <- sort(unique(as.character(meta$phenotype)))
      disease_pairs <- combn(all_groups, 2, simplify = FALSE)
      n_pairs_total <- length(disease_pairs)
      
      scenario_idx <- match(scenario, scenarios_to_run)
      rep_idx      <- as.integer(rep)
      
      # 用于在raw里保存"本rep下的3个pair标签集合"
      comparison_set_label <- paste(
        vapply(disease_pairs, \(p) paste0(p[1], "_vs_", p[2]), ""),
        collapse = ";"
      )
      
      for (i in seq_along(disease_pairs)) {
        
        pair       <- disease_pairs[[i]]
        ref_group  <- if ("Healthy" %in% pair) "Healthy" else pair[1]
        ref_group  <- as.character(ref_group)[1]
        test_group <- as.character(setdiff(pair, ref_group))[1]
        comp_name  <- paste0(test_group, "_vs_", ref_group)
        
        # ✅ Jaccard seed 公式（与 lefse/maaslin naive 一致）
        pair_seed <- 42L +
          1000000L * scenario_idx +
          10000L   * rep_idx +
          100L     * i
        
        cat("  比较:", comp_name, "\n"); flush.console()
        
        samples_all <- meta$sample_id[meta$phenotype %in% pair]
        if (length(samples_all) < 4) {
          cat("  样本数不足，跳过\n"); next
        }
        
        counts_all <- counts_matrix[, samples_all, drop = FALSE]
        counts_all <- ensure_matrix(counts_all)
        counts_all <- counts_all[
          rowSums(counts_all) > 0, , drop = FALSE
        ]
        counts_all <- ensure_matrix(counts_all)
        if (nrow(counts_all) == 0) next
        
        meta_all <- meta[meta$sample_id %in% samples_all, , drop = FALSE]
        
        # ── Step 1：全量运行 ──────────────────────────
        S_all <- tryCatch(
          run_deseq2_core(
            counts_matrix    = counts_all,
            metadata_df      = meta_all,
            test_group       = test_group,
            ref_group        = ref_group,
            use_cooks        = use_cooks,
            qval_cutoff      = QVAL_THR,
            lfc_cutoff       = LFC_THR,
            top_k            = TOP_K,
            debug            = FALSE
          ),
          error = function(e) {
            cat("  全量DESeq2失败:", e$message, "\n")
            character(0)
          }
        )
        cat("  全量特征数:", length(S_all), "\n"); flush.console()
        
        if (length(S_all) == 0) {
          jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
            method            = paste0("DESeq2_naive_", setting_tag),
            scenario          = scenario,
            rep               = rep_idx,
            comparison        = comp_name,
            comparison_set    = comparison_set_label,
            top_k_per_pair    = as.integer(TOP_K),
            n_pairs           = as.integer(n_pairs_total),
            top_k_total       = as.integer(TOP_K * n_pairs_total),
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
        
        # ── Step 2：B 次 subsampling（按 phenotype 分层）──
        jac_vals <- rep(NA_real_, B)
        
        for (b in seq_len(B)) {
          
          set.seed(pair_seed + b)
          
          # ✅ Naive: 只按 phenotype 分层抽样
          sub_ids <- meta_all %>%
            dplyr::group_by(phenotype) %>%
            dplyr::slice_sample(prop = SUBSAMPLE) %>%
            dplyr::ungroup() %>%
            dplyr::pull(sample_id) %>%
            as.character()
          
          if (length(sub_ids) < 4) next
          
          counts_sub <- counts_all[, sub_ids, drop = FALSE]
          counts_sub <- ensure_matrix(counts_sub)
          counts_sub <- counts_sub[
            rowSums(counts_sub) > 0, , drop = FALSE
          ]
          counts_sub <- ensure_matrix(counts_sub)
          if (nrow(counts_sub) == 0 || ncol(counts_sub) == 0) next
          
          meta_sub <- meta_all[
            meta_all$sample_id %in% sub_ids, , drop = FALSE
          ]
          rownames(meta_sub) <- NULL
          
          S_b <- tryCatch(
            run_deseq2_core(
              counts_matrix    = counts_sub,
              metadata_df      = meta_sub,
              test_group       = test_group,
              ref_group        = ref_group,
              use_cooks        = use_cooks,
              qval_cutoff      = QVAL_THR,
              lfc_cutoff       = LFC_THR,
              top_k            = TOP_K,
              debug            = FALSE
            ),
            error = function(e) character(0)
          )
          
          jac_vals[b] <- jaccard(S_b, S_all)
        }
        
        jac_clean     <- jac_vals[!is.na(jac_vals)]
        n_boots_valid <- length(jac_clean)
        
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method            = paste0("DESeq2_naive_", setting_tag),
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
        
        cat(sprintf(
          "  [%s | %s | rep%s | %s] Jaccard均值=%.3f (SD=%.3f, 有效boots=%d/%d)\n",
          scenario, comp_name, rep, setting_tag,
          if (n_boots_valid == 0) NA else mean(jac_clean),
          if (n_boots_valid == 0) NA else sd(jac_clean),
          n_boots_valid, B
        ))
        flush.console()
      }
    }
  }
  
  # ── 汇总并保存 ────────────────────────────────────────────
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    empty <- empty_jaccard_tibble()
    saveRDS(empty, paste0("benchmark/results/deseq2_naive_jaccard_raw_", setting_tag, ".rds"))
    saveRDS(empty, paste0("benchmark/results/deseq2_naive_jaccard_pair_summary_", setting_tag, ".rds"))
    saveRDS(empty, paste0("benchmark/results/deseq2_naive_jaccard_summary_", setting_tag, ".rds"))
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
  
  # 2) 场景级汇总：每rep先对3个pair取均值，再 across reps
  jaccard_rep_pairavg <- jaccard_raw %>%
    dplyr::group_by(method, scenario, rep, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set = paste(sort(unique(comparison)), collapse = ";"),
      n_pairs_distinct = dplyr::n_distinct(comparison),
      
      # 防止同组内 top_k_per_pair 不一致/非标量
      top_k_scalar = dplyr::first(top_k_per_pair),
      
      top_k_total = as.integer(top_k_scalar * n_pairs_distinct),
      
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    dplyr::select(-top_k_scalar) %>%              # 可选：不想保留就删掉
    dplyr::rename(n_pairs = n_pairs_distinct)
  
  
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
  
  cat("\n── DESeq2 Naive [", setting_tag, "] Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  
  cat("\n── DESeq2 Naive [", setting_tag, "] Jaccard 稳定性汇总（场景级：3个pair取平均） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,              paste0("benchmark/results/deseq2_naive_jaccard_raw_", setting_tag, ".rds"))
  saveRDS(jaccard_pair_summary,     paste0("benchmark/results/deseq2_naive_jaccard_pair_summary_", setting_tag, ".rds"))
  saveRDS(jaccard_scenario_summary, paste0("benchmark/results/deseq2_naive_jaccard_summary_", setting_tag, ".rds"))
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# ══════════════════════════════════════════════════════════════
# 执行
# ══════════════════════════════════════════════════════════════
jaccard_deseq2_naive_default  <- run_deseq2_naive_jaccard(use_cooks = TRUE)
jaccard_deseq2_naive_no_cooks <- run_deseq2_naive_jaccard(use_cooks = FALSE)



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

library(DESeq2)
library(tidyverse)

# ══════════════════════════════════════════════════════════════
# 全局参数（只定义一次）
# ══════════════════════════════════════════════════════════════
scenarios_to_run <- "B_batch_confound"   # 只做场景B

reps <- as.character(1:3)                # 只做3个rep

QVAL_THR  <- 0.05
LFC_THR   <- 1      # 保留参数但不用于significant判定
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
# 核心函数：针对单个 pair 运行 DESeq2（含batch协变量），返回显著特征名
# 修改要点：
#   ① 加入 batch 参数，设计公式改为 ~ batch + phenotype
#   ② significant 仅由 q值 < qval_cutoff 决定，不加LFC限制
# ══════════════════════════════════════════════════════════════
run_deseq2_core <- function(counts_matrix,
                            metadata_df,
                            test_group,
                            ref_group,
                            use_batch        = TRUE,   # [修改] 是否加入batch协变量
                            use_cooks        = TRUE,
                            qval_cutoff      = 0.05,
                            lfc_cutoff       = 1,      # 保留但不用于显著性判定
                            top_k            = NULL,
                            debug            = FALSE) {
  
  counts_matrix <- ensure_matrix(counts_matrix)
  
  if (debug) {
    cat("  [DEBUG] dim:", paste(dim(counts_matrix), collapse = "x"), "\n")
  }
  
  test_group <- as.character(test_group)[1]
  ref_group  <- as.character(ref_group)[1]
  
  if (!"sample_id" %in% colnames(metadata_df)) {
    metadata_df <- tibble::rownames_to_column(
      as.data.frame(metadata_df), "sample_id"
    )
  }
  rownames(metadata_df) <- NULL
  
  pair        <- c(test_group, ref_group)
  metadata_df <- metadata_df[metadata_df$phenotype %in% pair, , drop = FALSE]
  rownames(metadata_df) <- NULL
  
  common_samples <- intersect(colnames(counts_matrix),
                              metadata_df$sample_id)
  if (length(common_samples) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  metadata_df <- metadata_df[
    match(common_samples, metadata_df$sample_id), , drop = FALSE
  ]
  rownames(metadata_df) <- NULL
  
  counts_matrix <- round(pmax(counts_matrix, 0))
  
  row_keep <- rowSums(counts_matrix) > 0
  counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
    rowMeans(counts_matrix) > 1
  counts_matrix <- counts_matrix[keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  group_counts <- table(metadata_df$phenotype)
  if (!all(pair %in% names(group_counts))) return(character(0))
  if (any(group_counts[pair] < 2))         return(character(0))
  
  col_data           <- as.data.frame(metadata_df)
  rownames(col_data) <- col_data$sample_id
  col_data$sample_id <- NULL
  col_data$phenotype <- factor(col_data$phenotype,
                               levels = c(ref_group, test_group))
  
  # [修改] 根据 use_batch 决定设计公式
  # 同时检查：若使用batch，需要batch列存在且有多于1个水平
  use_batch_actual <- use_batch &&
    "batch" %in% colnames(col_data) &&
    length(unique(col_data$batch)) > 1
  
  # [修改] 额外检查：batch与phenotype不能完全共线
  # 若每个batch内只有一种phenotype，则batch与phenotype完全混淆，
  # 无法同时估计，退回 ~ phenotype
  if (use_batch_actual) {
    cross_tab <- table(col_data$batch, col_data$phenotype)
    # 每行（batch）至少有两个非零cell才算可以分离
    rows_ok <- apply(cross_tab, 1, function(r) sum(r > 0) >= 2)
    if (!any(rows_ok)) {
      use_batch_actual <- FALSE
      if (debug) cat("  [DEBUG] batch与phenotype完全共线，退回~phenotype\n")
    }
  }
  
  design_formula <- if (use_batch_actual) ~ batch + phenotype else ~ phenotype
  
  if (debug) cat("  [DEBUG] design:", deparse(design_formula), "\n")
  
  dds <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = counts_matrix,
      colData   = col_data,
      design    = design_formula
    )
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- estimateDispersions(dds, fitType = "local")
    dds <- nbinomWaldTest(dds)
    dds
  }, error = function(e) {
    if (debug) cat("  [DEBUG] DESeq2失败:", e$message, "\n")
    NULL
  })
  
  if (is.null(dds)) return(character(0))
  
  res <- tryCatch(
    results(
      dds,
      contrast      = c("phenotype", test_group, ref_group),
      alpha         = qval_cutoff,
      pAdjustMethod = "BH",
      cooksCutoff   = use_cooks
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(character(0))
  
  # [确认②] significant 仅由 q值 < qval_cutoff 决定，不加LFC限制
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("taxon") %>%
    filter(
      !is.na(padj),
      padj < qval_cutoff      # 仅q值阈值，无LFC过滤
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  if (nrow(res_df) == 0) return(character(0))
  if (!is.null(top_k)) res_df <- slice_head(res_df, n = top_k)
  
  res_df$taxon
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 DESeq2（Adjusted版本）
# ══════════════════════════════════════════════════════════════
run_mock_deseq2_genus_adjusted <- function(scenario_full, rep_id, use_cooks = TRUE) {
  
  set.seed(42)
  
  setting_tag <- if (use_cooks) "default" else "no_cooks"
  
  cat("开始处理:", scenario_full, "rep", rep_id, "[", setting_tag, "]\n")
  flush.console()
  
  t_wall_start <- as.numeric(proc.time()["elapsed"])
  
  # 1. 读取数据
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$deseq2$counts
  metadata_df   <- inputs$deseq2$groups
  
  cat("特征矩阵维度:", dim(counts_matrix), "\n")
  cat("样本数:", nrow(metadata_df), "\n")
  flush.console()
  
  # 2. 数据清洗
  counts_matrix <- ensure_matrix(counts_matrix)
  counts_matrix <- round(pmax(counts_matrix, 0))
  
  sample_sums <- colSums(counts_matrix)
  if (any(sample_sums == 0)) {
    cat("警告: 发现", sum(sample_sums == 0), "个全 0 样本，已删除\n")
    counts_matrix <- counts_matrix[, sample_sums > 0, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
  }
  
  feature_sums <- rowSums(counts_matrix)
  if (any(feature_sums == 0)) {
    cat("警告: 发现", sum(feature_sums == 0), "个全 0 特征，已删除\n")
    counts_matrix <- counts_matrix[feature_sums > 0, , drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
  }
  
  # Adjusted: 保留 sample_id, phenotype, batch
  meta <- metadata_df %>%
    dplyr::select(sample_id, phenotype, batch) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      phenotype = as.character(phenotype),
      batch     = as.factor(batch)
    )
  
  common_samples <- intersect(colnames(counts_matrix), meta$sample_id)
  counts_matrix  <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix  <- ensure_matrix(counts_matrix)
  meta           <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
  
  keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
    rowMeans(counts_matrix) > 1
  cat("低丰度过滤: 保留", sum(keep), "/", nrow(counts_matrix), "个特征\n")
  counts_matrix <- counts_matrix[keep, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  cat("对齐后样本数:", ncol(counts_matrix), "\n")
  cat("对齐后特征数:", nrow(counts_matrix), "\n")
  cat("分组分布:\n")
  print(table(meta$phenotype))
  cat("批次分布:\n")
  print(table(meta$batch))
  flush.console()
  
  # 3. 枚举所有两两比较对
  all_groups    <- sort(unique(as.character(meta$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n")
  cat("分组顺序:", paste(all_groups, collapse = ", "), "\n\n")
  flush.console()
  
  base_result_dir <- file.path(
    "benchmark", "deseq2",
    paste0("deseq2_results_genus_adjusted_", setting_tag),
    paste0(scenario_full, "_rep", rep_id)
  )
  dir.create(base_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  compute_sec  <- 0
  
  for (i in seq_along(disease_pairs)) {
    
    pair       <- disease_pairs[[i]]
    ref_group  <- if ("Healthy" %in% pair) "Healthy" else pair[1]
    ref_group  <- as.character(ref_group)[1]
    test_group <- as.character(setdiff(pair, ref_group))[1]
    comp_name  <- paste0(test_group, "_vs_", ref_group)
    pair_id    <- paste(sort(pair), collapse = "_")
    
    cat("比较", i, "/", length(disease_pairs), ":", comp_name, "\n")
    cat("  配对:", pair_id, "| 参考组:", ref_group,
        "| 测试组:", test_group, "\n")
    flush.console()
    
    pair_idx <- meta$phenotype %in% pair
    if (sum(pair_idx) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_counts <- counts_matrix[, pair_idx, drop = FALSE]
    pair_counts <- ensure_matrix(pair_counts)
    pair_meta   <- meta[pair_idx, , drop = FALSE]
    rownames(pair_meta) <- NULL
    
    nonzero_taxa <- rowSums(pair_counts) > 0
    pair_counts  <- pair_counts[nonzero_taxa, , drop = FALSE]
    pair_counts  <- ensure_matrix(pair_counts)
    
    cat("  样本数:", ncol(pair_counts),
        "| 特征数:", nrow(pair_counts), "\n")
    cat("  分组分布:\n")
    print(table(pair_meta$phenotype))
    cat("  批次分布:\n")
    print(table(pair_meta$batch))
    flush.console()
    
    t_compute_start <- as.numeric(proc.time()["elapsed"])
    
    col_data           <- as.data.frame(pair_meta)
    rownames(col_data) <- col_data$sample_id
    col_data$sample_id <- NULL
    col_data$phenotype <- factor(col_data$phenotype,
                                 levels = c(ref_group, test_group))
    # batch已经是factor，保持不变
    
    # [确认④] Adjusted: 使用 ~ batch + phenotype 设计
    # 检查batch与phenotype是否完全共线
    use_batch_actual <- length(unique(col_data$batch)) > 1
    if (use_batch_actual) {
      cross_tab <- table(col_data$batch, col_data$phenotype)
      rows_ok <- apply(cross_tab, 1, function(r) sum(r > 0) >= 2)
      if (!any(rows_ok)) {
        use_batch_actual <- FALSE
        cat("  警告: batch与phenotype完全共线，退回~phenotype\n")
      }
    }
    design_formula <- if (use_batch_actual) ~ batch + phenotype else ~ phenotype
    cat("  设计公式:", deparse(design_formula), "\n")
    
    dds <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = pair_counts,
        colData   = col_data,
        design    = design_formula
      )
      dds <- estimateSizeFactors(dds, type = "poscounts")
      dds <- estimateDispersions(dds, fitType = "local")
      dds <- nbinomWaldTest(dds)
      dds
    }, error = function(e) {
      cat("  DESeq2 运行失败:", e$message, "\n"); NULL
    })
    
    t_compute_end <- as.numeric(proc.time()["elapsed"])
    compute_sec   <- compute_sec + (t_compute_end - t_compute_start)
    
    if (is.null(dds)) next
    
    res <- tryCatch({
      results(
        dds,
        contrast      = c("phenotype", test_group, ref_group),
        alpha         = 0.05,
        pAdjustMethod = "BH",
        cooksCutoff   = use_cooks
      )
    }, error = function(e) {
      cat("  提取结果失败:", e$message, "\n"); NULL
    })
    
    if (is.null(res)) next
    
    result_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("taxon") %>%
      dplyr::rename(
        base_mean = baseMean,
        log2fc    = log2FoldChange,
        lfc_se    = lfcSE,
        stat      = stat,
        p_val     = pvalue,
        q_val     = padj
      ) %>%
      dplyr::mutate(
        # [确认②] sig_fdr: 仅 q < 0.05，无LFC限制
        sig_fdr      = !is.na(q_val) & q_val < 0.05,
        # sig_fdr_eff 与 sig_fdr 相同（均不加LFC限制）
        sig_fdr_eff  = !is.na(q_val) & q_val < 0.05,
        comparison   = comp_name,
        pair_id      = pair_id,
        test_group   = test_group,
        ref_group    = ref_group,
        scenario     = scenario_full,
        replicate    = rep_id,
        setting      = setting_tag
      ) %>%
      dplyr::select(
        taxon, base_mean, log2fc, lfc_se, stat, p_val, q_val,
        sig_fdr, sig_fdr_eff,
        comparison, pair_id, test_group, ref_group,
        scenario, replicate, setting
      )
    
    n_p_sig <- sum(result_df$p_val < 0.05, na.rm = TRUE)
    n_q_sig <- sum(result_df$sig_fdr,      na.rm = TRUE)
    
    cat("  发现", nrow(result_df), "个特征结果\n")
    cat("  p < 0.05:", n_p_sig, "个\n")
    cat("  q < 0.05 (FDR，无LFC限制):", n_q_sig, "个\n")
    flush.console()
    
    comp_output_dir <- file.path(base_result_dir, pair_id)
    dir.create(comp_output_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv(result_df, file.path(comp_output_dir, "results.csv"))
    
    results_list[[comp_name]] <- result_df
  }
  
  if (length(results_list) > 0) {
    all_results_df <- bind_rows(results_list)
    write_csv(all_results_df,
              file.path(base_result_dir, "all_comparisons_results.csv"))
    cat("结果已保存至:", base_result_dir, "\n")
    cat("\n=== 汇总统计 ===\n")
    cat("sig_fdr (q < 0.05，无LFC限制):",
        sum(all_results_df$sig_fdr, na.rm = TRUE), "个\n")
  } else {
    cat("无有效比较结果\n")
  }
  
  t_wall_end <- as.numeric(proc.time()["elapsed"])
  wall_sec   <- t_wall_end - t_wall_start
  
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
# 批量运行 + timing（仅场景B，3个rep）
# [确认③] scenarios_to_run = "B_batch_confound", reps = 1:3
# ══════════════════════════════════════════════════════════════
run_all_adjusted <- function(use_cooks) {
  
  tag <- if (use_cooks) "default" else "no_cooks"
  cat("\n=============================\n")
  cat("开始运行 DESeq2 Adjusted [", tag, "]\n")
  cat("=============================\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (scenario in scenarios_to_run) {      # 只有 "B_batch_confound"
    for (rep in reps) {                   # 只有 "1","2","3"
      
      cat("\n--- 场景:", scenario, "Rep:", rep, "---\n")
      flush.console()
      
      key <- paste0(scenario, "_rep", rep)
      
      out <- tryCatch(
        run_mock_deseq2_genus_adjusted(scenario, rep, use_cooks = use_cooks),
        error = function(e) { cat("失败:", e$message, "\n"); NULL }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out$results
        timing_records[[length(timing_records) + 1]] <- tibble(
          method      = paste0("DESeq2_adjusted_", tag),
          scenario    = scenario,
          rep         = as.integer(rep),
          wall_sec    = as.numeric(out$wall_sec),
          compute_sec = as.numeric(out$compute_sec)
        )
      } else {
        results_all[[key]] <- NULL
      }
    }
  }
  
  cat("\nDESeq2 Adjusted [", tag, "] 完成，成功:",
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
  
  cat("\n── DESeq2 Adjusted [", tag, "] 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(timing_df,
          paste0("benchmark/results/deseq2_adjusted_timing_raw_",     tag, ".rds"))
  saveRDS(timing_summary,
          paste0("benchmark/results/deseq2_adjusted_timing_summary_", tag, ".rds"))
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

out_adjusted_default  <- run_all_adjusted(use_cooks = TRUE)
out_adjusted_no_cooks <- run_all_adjusted(use_cooks = FALSE)

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV（Adjusted版本）
# ══════════════════════════════════════════════════════════════
collect_deseq2_results_adjusted <- function(setting_tag = "default",
                                            scenarios, reps) {
  base_dir <- file.path(
    "benchmark", "deseq2",
    paste0("deseq2_results_genus_adjusted_", setting_tag)
  )
  
  all_rows <- list()
  for (scenario in scenarios) {
    for (rep in reps) {
      result_dir <- file.path(base_dir,
                              paste0(scenario, "_rep", rep))
      if (!dir.exists(result_dir)) next
      
      for (comp_dir in list.dirs(result_dir, recursive = FALSE)) {
        res_file <- file.path(comp_dir, "results.csv")
        if (!file.exists(res_file)) next
        df <- tryCatch(
          read_csv(res_file, show_col_types = FALSE),
          error = function(e) NULL
        )
        if (!is.null(df) && nrow(df) > 0)
          all_rows[[length(all_rows) + 1]] <- df
      }
    }
  }
  
  if (length(all_rows) == 0) { cat("未找到结果\n"); return(NULL) }
  bind_rows(all_rows)
}

summary_adjusted_default  <- collect_deseq2_results_adjusted(
  "default",  scenarios_to_run, reps)
summary_adjusted_no_cooks <- collect_deseq2_results_adjusted(
  "no_cooks", scenarios_to_run, reps)

dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(summary_adjusted_default,  "benchmark/results/deseq2_adjusted_summary_default.rds")
saveRDS(summary_adjusted_no_cooks, "benchmark/results/deseq2_adjusted_summary_no_cooks.rds")
write_csv(summary_adjusted_default,  "benchmark/results/deseq2_adjusted_summary_default.csv")
write_csv(summary_adjusted_no_cooks, "benchmark/results/deseq2_adjusted_summary_no_cooks.csv")

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估（Adjusted版本）
# 修改要点：
#   ④ run_deseq2_core 加入 batch 协变量（use_batch=TRUE）
#   ④ subsampling 改为按 phenotype × batch 联合分层
#   ⑤ 每个pair单独记录，综合时用3对pair的Jaccard均值
# ══════════════════════════════════════════════════════════════
run_deseq2_adjusted_jaccard <- function(use_cooks = TRUE) {
  
  setting_tag     <- if (use_cooks) "default" else "no_cooks"
  jaccard_records <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("DESeq2 Adjusted Jaccard 稳定性评估 [", setting_tag, "]\n")
  cat("══════════════════════════════════════════════════════════════\n")
  
  for (scenario in scenarios_to_run) {          # 只有 "B_batch_confound"
    for (rep in reps) {                       # 只有 "1","2","3"
      
      cat("\n[Jaccard] 场景:", scenario, "Rep:", rep,
          "[", setting_tag, "]\n")
      flush.console()
      
      inputs <- tryCatch(
        generate_method_inputs(scenario, rep),
        error = function(e) {
          cat("数据读取失败:", e$message, "\n"); NULL
        }
      )
      if (is.null(inputs)) next
      
      counts_matrix <- inputs$deseq2$counts
      metadata_df   <- inputs$deseq2$groups
      
      counts_matrix <- ensure_matrix(counts_matrix)
      counts_matrix <- round(pmax(counts_matrix, 0))
      
      col_keep      <- colSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      row_keep      <- rowSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # [确认④] Adjusted: metadata 保留 phenotype 和 batch
      meta <- metadata_df %>%
        dplyr::select(sample_id, phenotype, batch) %>%
        dplyr::mutate(
          sample_id = as.character(sample_id),
          phenotype = as.character(phenotype),
          batch     = as.factor(batch)
        )
      
      common <- intersect(colnames(counts_matrix), meta$sample_id)
      counts_matrix <- counts_matrix[, common, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      meta <- meta[match(common, meta$sample_id), , drop = FALSE]
      
      keep <- rowSums(counts_matrix > 0) >= (ncol(counts_matrix) * 0.2) &
        rowMeans(counts_matrix) > 1
      counts_matrix <- counts_matrix[keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
        cat("  过滤后矩阵为空，跳过\n"); next
      }
      
      cat("  [预处理后] dim:",
          paste(dim(counts_matrix), collapse = "x"), "\n")
      cat("  批次分布:\n")
      print(table(meta$batch, meta$phenotype))
      
      group_counts <- table(meta$phenotype)
      valid_groups <- names(group_counts[group_counts >= 2])
      if (length(valid_groups) < 2) next
      
      meta <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
      counts_matrix <- counts_matrix[, meta$sample_id, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      all_groups    <- sort(unique(as.character(meta$phenotype)))
      disease_pairs <- combn(all_groups, 2, simplify = FALSE)
      n_pairs_total <- length(disease_pairs)
      
      scenario_idx <- match(scenario, scenarios_to_run)
      rep_idx      <- as.integer(rep)
      
      # [确认⑤] comparison_set：列出所有pair的名称
      comparison_set_label <- paste(
        vapply(disease_pairs,
               function(p) {
                 rg <- if ("Healthy" %in% p) "Healthy" else p[1]
                 tg <- setdiff(p, rg)[1]
                 paste0(tg, "_vs_", rg)
               }, ""),
        collapse = ";"
      )
      
      # [确认⑤] 每个pair单独循环记录
      for (i in seq_along(disease_pairs)) {
        
        pair       <- disease_pairs[[i]]
        ref_group  <- if ("Healthy" %in% pair) "Healthy" else pair[1]
        ref_group  <- as.character(ref_group)[1]
        test_group <- as.character(setdiff(pair, ref_group))[1]
        comp_name  <- paste0(test_group, "_vs_", ref_group)
        
        pair_seed <- 42L +
          1000000L * scenario_idx +
          10000L   * rep_idx +
          100L     * i
        
        cat("  比较:", comp_name, "\n"); flush.console()
        
        samples_all <- meta$sample_id[meta$phenotype %in% pair]
        if (length(samples_all) < 4) {
          cat("  样本数不足，跳过\n"); next
        }
        
        counts_all <- counts_matrix[, samples_all, drop = FALSE]
        counts_all <- ensure_matrix(counts_all)
        counts_all <- counts_all[
          rowSums(counts_all) > 0, , drop = FALSE
        ]
        counts_all <- ensure_matrix(counts_all)
        if (nrow(counts_all) == 0) next
        
        meta_all <- meta[meta$sample_id %in% samples_all, , drop = FALSE]
        rownames(meta_all) <- NULL
        
        # Step 1：全量运行（含batch，仅qval阈值）
        # [确认④] 全量运行时使用 use_batch=TRUE
        S_all <- tryCatch(
          run_deseq2_core(
            counts_matrix = counts_all,
            metadata_df   = meta_all,
            test_group    = test_group,
            ref_group     = ref_group,
            use_batch     = TRUE,      # [修改] 加入batch
            use_cooks     = use_cooks,
            qval_cutoff   = QVAL_THR,
            lfc_cutoff    = LFC_THR,   # 不参与过滤
            top_k         = TOP_K,
            debug         = FALSE
          ),
          error = function(e) {
            cat("  全量DESeq2失败:", e$message, "\n")
            character(0)
          }
        )
        cat("  全量特征数:", length(S_all), "\n"); flush.console()
        
        if (length(S_all) == 0) {
          jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
            method            = paste0("DESeq2_adjusted_", setting_tag),
            scenario          = scenario,
            rep               = rep_idx,
            comparison        = comp_name,
            comparison_set    = comparison_set_label,
            top_k_per_pair    = as.integer(TOP_K),
            n_pairs           = as.integer(n_pairs_total),
            top_k_total       = as.integer(TOP_K * n_pairs_total),
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
        
        # Step 2：B 次 subsampling
        # [修改④] 按 phenotype × batch 联合分层抽样
        jac_vals <- rep(NA_real_, B)
        
        for (b in seq_len(B)) {
          
          set.seed(pair_seed + b)
          
          # [修改] 按 phenotype × batch 联合分层
          # 每个 phenotype×batch 组合内按 SUBSAMPLE 比例抽取
          # 若某组合只有1个样本，仍保留（至少取1个）
          sub_ids <- meta_all %>%
            dplyr::group_by(phenotype, batch) %>%
            dplyr::slice_sample(
              prop = SUBSAMPLE,
              replace = FALSE
            ) %>%
            dplyr::ungroup() %>%
            dplyr::pull(sample_id) %>%
            as.character()
          
          # 确保每个phenotype至少2个样本
          pheno_check <- meta_all %>%
            dplyr::filter(sample_id %in% sub_ids) %>%
            dplyr::count(phenotype)
          if (any(pheno_check$n < 2)) next
          if (length(sub_ids) < 4) next
          
          counts_sub <- counts_all[, sub_ids, drop = FALSE]
          counts_sub <- ensure_matrix(counts_sub)
          counts_sub <- counts_sub[
            rowSums(counts_sub) > 0, , drop = FALSE
          ]
          counts_sub <- ensure_matrix(counts_sub)
          if (nrow(counts_sub) == 0 || ncol(counts_sub) == 0) next
          
          meta_sub <- meta_all[
            meta_all$sample_id %in% sub_ids, , drop = FALSE
          ]
          rownames(meta_sub) <- NULL
          
          # [修改④] subsampling 运行时也使用 use_batch=TRUE
          S_b <- tryCatch(
            run_deseq2_core(
              counts_matrix = counts_sub,
              metadata_df   = meta_sub,
              test_group    = test_group,
              ref_group     = ref_group,
              use_batch     = TRUE,      # [修改] 加入batch
              use_cooks     = use_cooks,
              qval_cutoff   = QVAL_THR,
              lfc_cutoff    = LFC_THR,   # 不参与过滤
              top_k         = TOP_K,
              debug         = FALSE
            ),
            error = function(e) character(0)
          )
          
          jac_vals[b] <- jaccard(S_b, S_all)
        }
        
        jac_clean     <- jac_vals[!is.na(jac_vals)]
        n_boots_valid <- length(jac_clean)
        
        # [确认⑤] 每个pair单独存一行记录
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method            = paste0("DESeq2_adjusted_", setting_tag),
          scenario          = scenario,
          rep               = rep_idx,
          comparison        = comp_name,           # 当前pair
          comparison_set    = comparison_set_label, # 全部pair集合
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
          "  [%s | %s | rep%s | %s] Jaccard均值=%.3f (SD=%.3f, 有效boots=%d/%d)\n",
          scenario, comp_name, rep, setting_tag,
          if (n_boots_valid == 0) NA else mean(jac_clean),
          if (n_boots_valid == 0) NA else sd(jac_clean),
          n_boots_valid, B
        ))
        flush.console()
      }  # end pair loop
    }  # end rep loop
  }  # end scenario loop
  
  # ══════════════════════════════════════════════════════════════
  # 汇总并保存
  # ══════════════════════════════════════════════════════════════
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    empty <- empty_jaccard_tibble()
    saveRDS(empty, paste0("benchmark/results/deseq2_adjusted_jaccard_raw_",
                          setting_tag, ".rds"))
    saveRDS(empty, paste0("benchmark/results/deseq2_adjusted_jaccard_pair_summary_",
                          setting_tag, ".rds"))
    saveRDS(empty, paste0("benchmark/results/deseq2_adjusted_jaccard_summary_",
                          setting_tag, ".rds"))
    return(list(
      jaccard_raw              = empty,
      jaccard_pair_summary     = empty,
      jaccard_scenario_summary = empty
    ))
  }
  
  # raw：每行 = 一个 scenario × rep × pair 的Jaccard统计
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
  # [确认⑤] Pair级汇总：每个 scenario × comparison 跨3个rep取均值
  jaccard_pair_summary <- jaccard_raw %>%
    dplyr::group_by(method, scenario, comparison, top_k_per_pair) %>%
    dplyr::summarise(
      jaccard_mean_across_reps   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_across_reps     = mean(jaccard_sd,     na.rm = TRUE),
      n_boots                    = max(n_boots,         na.rm = TRUE),
      n_boots_valid_mean         = mean(n_boots_valid,  na.rm = TRUE),
      n_features_full_mean       = mean(n_features_full,na.rm = TRUE),
      n_reps                     = dplyr::n(),
      .groups = "drop"
    )
  
  # [确认⑤] 场景级汇总：每个场景一行
  # 先在每个 scenario × rep 内对所有pair取均值，再跨rep取均值
  jaccard_rep_pairavg <- jaccard_raw %>%
    dplyr::group_by(method, scenario, rep, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set         = dplyr::first(comparison_set),
      n_pairs                = dplyr::n_distinct(comparison),  # 实际pair数
      top_k_total            = as.integer(top_k_per_pair[1] * dplyr::n_distinct(comparison)),
      # [⑤] 用该场景×rep下全部pair的Jaccard均值作为该rep的汇总值
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      n_boots_valid_sum      = sum(n_boots_valid,   na.rm = TRUE),
      n_features_full_mean   = mean(n_features_full,na.rm = TRUE),
      .groups = "drop"
    )
  
  # 场景级：每个场景一行（跨3个rep再平均）
  jaccard_scenario_summary <- jaccard_rep_pairavg %>%
    dplyr::group_by(method, scenario, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set               = dplyr::first(comparison_set),
      n_pairs_mean                 = mean(n_pairs,    na.rm = TRUE),
      top_k_total_mean             = mean(top_k_total,na.rm = TRUE),
      # 每个场景一个数：3rep × 3pair 的Jaccard总平均
      jaccard_mean_across_reps     = mean(jaccard_mean_pairavg,   na.rm = TRUE),
      jaccard_median_across_reps   = mean(jaccard_median_pairavg, na.rm = TRUE),
      jaccard_sd_across_reps       = sd(jaccard_mean_pairavg,     na.rm = TRUE),
      n_reps                       = dplyr::n(),
      n_boots_valid_mean           = mean(n_boots_valid_sum,      na.rm = TRUE),
      n_features_full_mean         = mean(n_features_full_mean,   na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n── DESeq2 Adjusted [", setting_tag,
      "] Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  
  cat("\n── DESeq2 Adjusted [", setting_tag,
      "] Jaccard 稳定性汇总（场景级，每场景一行） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,
          paste0("benchmark/results/deseq2_adjusted_jaccard_raw_",
                 setting_tag, ".rds"))
  saveRDS(jaccard_pair_summary,
          paste0("benchmark/results/deseq2_adjusted_jaccard_pair_summary_",
                 setting_tag, ".rds"))
  saveRDS(jaccard_scenario_summary,
          paste0("benchmark/results/deseq2_adjusted_jaccard_summary_",
                 setting_tag, ".rds"))
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# ══════════════════════════════════════════════════════════════
# 执行
# ══════════════════════════════════════════════════════════════
jaccard_deseq2_adjusted_default  <- run_deseq2_adjusted_jaccard(use_cooks = TRUE)
jaccard_deseq2_adjusted_no_cooks <- run_deseq2_adjusted_jaccard(use_cooks = FALSE)
