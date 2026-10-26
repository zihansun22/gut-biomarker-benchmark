
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

library(ANCOMBC)
library(phyloseq)
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
LFC_THR   <- 1        # ANCOMBC不强制LFC过滤，保留参数供扩展
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
# 核心函数：从 counts 矩阵 + meta 运行 ANCOM-BC2，返回显著特征名
# debug = TRUE 时打印维度信息
# ══════════════════════════════════════════════════════════════
run_ancombc2_core <- function(counts_matrix,
                              metadata_df,
                              pair_levels,
                              qval_cutoff = 0.05,
                              top_k = NULL,
                              debug = FALSE) {
  
  # counts_matrix: taxa × samples（numeric matrix）
  # metadata_df  : data.frame，至少含 sample_id, phenotype（不含batch）
  # pair_levels  : c("G1","G2")
  
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # metadata 规范化（确保 sample_id 列存在）
  if (!"sample_id" %in% colnames(metadata_df)) {
    metadata_df <- tibble::rownames_to_column(as.data.frame(metadata_df), "sample_id")
  }
  
  metadata_df <- metadata_df %>%
    dplyr::select(sample_id, phenotype) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      phenotype = as.character(phenotype)
    )
  
  # 对齐样本
  common <- intersect(colnames(counts_matrix), metadata_df$sample_id)
  if (length(common) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, common, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  metadata_df   <- metadata_df[match(common, metadata_df$sample_id), , drop = FALSE]
  
  # 只保留 pair
  metadata_df <- metadata_df[metadata_df$phenotype %in% pair_levels, , drop = FALSE]
  if (nrow(metadata_df) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, metadata_df$sample_id, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # 至少两组
  if (dplyr::n_distinct(metadata_df$phenotype) < 2) return(character(0))
  
  # 去全0特征
  counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  # 设置参考组
  ref_group <- if ("Healthy" %in% pair_levels) "Healthy" else sort(pair_levels)[1]
  metadata_df$phenotype <- stats::relevel(factor(metadata_df$phenotype, levels = pair_levels),
                                          ref = ref_group)
  
  # 构建 phyloseq
  sample_data_df <- tibble::column_to_rownames(as.data.frame(metadata_df), "sample_id")
  physeq <- tryCatch(
    phyloseq::phyloseq(
      phyloseq::otu_table(counts_matrix, taxa_are_rows = TRUE),
      phyloseq::sample_data(sample_data_df)
    ),
    error = function(e) NULL
  )
  if (is.null(physeq)) return(character(0))
  
  # ANCOM-BC2（显式单线程）
  fit <- tryCatch(
    ANCOMBC::ancombc2(
      data         = physeq,
      fix_formula  = "phenotype",   # ✅ naive：不含batch
      p_adj_method = "BH",
      prv_cut      = 0,
      lib_cut      = 0,
      s0_perc      = 0.05,
      struc_zero   = FALSE,
      neg_lb       = FALSE,
      global       = FALSE,
      pairwise     = FALSE,
      dunnet       = FALSE,
      trend        = FALSE,
      alpha        = 0.05,
      n_cl         = 1,             # ✅ 单线程
      verbose      = FALSE
    ),
    error = function(e) {
      if (debug) cat("  [core] ANCOM-BC2失败:", e$message, "\n")
      NULL
    }
  )
  if (is.null(fit)) return(character(0))
  
  res_df <- fit$res
  if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
  
  # phenotype 的非参考水平（两组比较只有1列）
  nonref <- setdiff(pair_levels, ref_group)
  if (length(nonref) != 1) return(character(0))
  
  pheno_col <- paste0("phenotype", nonref)
  q_col     <- paste0("q_",   pheno_col)
  lfc_col   <- paste0("lfc_", pheno_col)
  
  if (!all(c("taxon", q_col, lfc_col) %in% colnames(res_df))) return(character(0))
  
  feats <- res_df %>%
    dplyr::filter(!is.na(.data[[q_col]]), .data[[q_col]] < qval_cutoff) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[lfc_col]])))
  
  if (!is.null(top_k) && top_k > 0) {
    feats <- dplyr::slice_head(feats, n = top_k)
  }
  
  unique(as.character(feats$taxon))
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 ANCOM-BC2，含完整保存逻辑
# wall_sec   ：总耗时（含 I/O）
# compute_sec：纯模型拟合耗时（仅 ancombc2() 调用时间之和）
# ══════════════════════════════════════════════════════════════
run_mock_ancombc2_genus <- function(scenario_full, rep_id) {
  
  set.seed(42)
  
  cat("开始处理:", scenario_full, "rep", rep_id, "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  # 1. 读取数据 ────────────────────────────────────────────
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$ancombc2$counts
  metadata_df   <- inputs$ancombc2$metadata
  
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
  
  # 先对齐样本，再做后续过滤
  common_samples <- intersect(colnames(counts_matrix),
                              metadata_df$sample_id)
  counts_matrix  <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix  <- ensure_matrix(counts_matrix)
  metadata_df    <- metadata_df[
    match(common_samples, metadata_df$sample_id), , drop = FALSE
  ]
  rownames(metadata_df) <- NULL
  
  cat("对齐后样本数:", ncol(counts_matrix), "\n")
  cat("分组分布:\n")
  print(table(metadata_df$phenotype))
  flush.console()
  
  # 3. 枚举所有两两比较对 ──────────────────────────────────
  all_groups    <- sort(unique(as.character(metadata_df$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n\n")
  flush.console()
  
  # ✅ 修改：输出目录改为 ancombc2_results_genus_naive
  base_result_dir <- file.path(
    "benchmark", "ancombc2", "ancombc2_results_genus_naive",
    paste0(scenario_full, "_rep", rep_id)
  )
  dir.create(base_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  compute_sec  <- 0  # 累加每次 ancombc2() 的耗时
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]
    g2        <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    
    cat("比较", i, "/", length(disease_pairs), ":", comp_name, "\n")
    flush.console()
    
    pair_idx <- metadata_df$phenotype %in% pair
    if (sum(pair_idx) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_counts <- counts_matrix[, pair_idx, drop = FALSE]
    pair_counts <- ensure_matrix(pair_counts)
    pair_meta   <- metadata_df[pair_idx, , drop = FALSE]
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
    
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    
    cat("  参考组:", ref_group, "\n")
    flush.console()
    
    sample_data_df  <- tibble::column_to_rownames(
      as.data.frame(pair_meta), "sample_id"
    )
    physeq <- tryCatch(
      phyloseq(
        otu_table(pair_counts, taxa_are_rows = TRUE),
        sample_data(sample_data_df)
      ),
      error = function(e) {
        cat("  phyloseq 构建失败:", e$message, "\n"); NULL
      }
    )
    if (is.null(physeq)) next
    
    cat("  phyloseq 对象构建完成\n")
    flush.console()
    
    # ── 纯模型拟合计时（每个 pair 单独计，累加）────────────
    t_compute_start <- proc.time()["elapsed"]
    
    ancombc_res <- tryCatch(
      ancombc2(
        data         = physeq,
        fix_formula  = "phenotype",
        p_adj_method = "BH",
        prv_cut      = 0,
        lib_cut      = 0,
        s0_perc      = 0.05,
        struc_zero   = FALSE,
        neg_lb       = FALSE,
        global       = FALSE,
        pairwise     = FALSE,
        dunnet       = FALSE,
        trend        = FALSE,
        alpha        = 0.05,
        n_cl         = 1,
        verbose      = FALSE
      ),
      error = function(e) {
        cat("  ANCOM-BC2 失败:", e$message, "\n"); NULL
      }
    )
    
    t_compute_end <- proc.time()["elapsed"]
    compute_sec   <- compute_sec + as.numeric(t_compute_end - t_compute_start)
    
    if (is.null(ancombc_res)) next
    
    res_df <- ancombc_res$res
    if (is.null(res_df) || nrow(res_df) == 0) {
      cat("  无有效结果，跳过\n"); next
    }
    
    pheno_col <- paste0("phenotype", setdiff(pair, ref_group))
    
    result_df <- data.frame(
      taxon      = res_df$taxon,
      lfc        = res_df[[paste0("lfc_",  pheno_col)]],
      se         = res_df[[paste0("se_",   pheno_col)]],
      W          = res_df[[paste0("W_",    pheno_col)]],
      p_val      = res_df[[paste0("p_",    pheno_col)]],
      q_val      = res_df[[paste0("q_",    pheno_col)]],
      diff_abn   = res_df[[paste0("diff_", pheno_col)]],
      comparison = comp_name,
      group_a    = g1,
      group_b    = g2,
      ref_group  = ref_group,
      scenario   = scenario_full,
      replicate  = rep_id,
      stringsAsFactors = FALSE
    )
    
    cat("  发现", nrow(result_df), "个特征结果，",
        sum(result_df$q_val < 0.05, na.rm = TRUE), "个 FDR<0.05\n")
    flush.console()
    
    # I/O（不计入 compute_sec）
    comp_output_dir <- file.path(base_result_dir, comp_name)
    dir.create(comp_output_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv(result_df, file.path(comp_output_dir, "results.csv"))
    
    results_list[[comp_name]] <- result_df
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
  cat("开始运行 ANCOM-BC2 Naive\n")
  cat("=============================\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n--- 场景:", scenario, "Rep:", rep, "---\n")
      flush.console()
      
      key <- paste0(scenario, "_rep", rep)
      
      out <- tryCatch(
        run_mock_ancombc2_genus(scenario, rep),
        error = function(e) { cat("失败:", e$message, "\n"); NULL }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out$results
        timing_records[[length(timing_records) + 1]] <- tibble(
          method      = "ANCOM-BC2_naive",
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
  
  cat("\nANCOM-BC2 Naive 完成，成功:",
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
  
  cat("\n── ANCOM-BC2 Naive 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  # ✅ 修改：保存为 ancombc2_naive_timing_*
  saveRDS(timing_df,      "benchmark/results/ancombc2_naive_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/ancombc2_naive_timing_summary.rds")
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

out_ancombc2 <- run_all()

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV
# ══════════════════════════════════════════════════════════════
collect_ancombc2_results <- function(
    base_dir  = "benchmark/ancombc2/ancombc2_results_genus_naive",
    scenarios, reps) {
  
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
  
  if (length(all_rows) == 0) { cat("未找到任何结果文件\n"); return(NULL) }
  bind_rows(all_rows)
}

ancombc2_summary_raw <- collect_ancombc2_results(
  scenarios = scenarios_to_run,
  reps      = reps
)

cat("汇总行数:", nrow(ancombc2_summary_raw), "\n")
print(head(ancombc2_summary_raw))

dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
# ✅ 修改：保存为 ancombc2_naive_summary_*
saveRDS(ancombc2_summary_raw, "benchmark/results/ancombc2_naive_summary_raw.rds")
write_csv(ancombc2_summary_raw, "benchmark/results/ancombc2_naive_summary_raw.csv")
cat("\n汇总结果已保存\n")

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估
# ══════════════════════════════════════════════════════════════
run_ancombc2_naive_jaccard <- function() {
  
  BASE_SEED <- 42L
  jaccard_records <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("ANCOM-BC2 Naive Jaccard 稳定性评估（pairwise + 场景内pair均值汇总）\n")
  cat("══════════════════════════════════════════════════════════════\n")
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n[Jaccard] 场景:", scenario, "Rep:", rep, "\n")
      flush.console()
      
      inputs <- tryCatch(
        generate_method_inputs(scenario, rep),
        error = function(e) { cat("数据读取失败:", e$message, "\n"); NULL }
      )
      if (is.null(inputs)) next
      
      counts_matrix <- inputs$ancombc2$counts
      metadata_df   <- inputs$ancombc2$metadata
      
      # ── 基本清洗 ─────────────────────────────────────────
      counts_matrix <- ensure_matrix(counts_matrix)
      counts_matrix <- round(pmax(counts_matrix, 0))
      
      # 去全0样本
      col_keep      <- colSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # 去全0特征
      row_keep      <- rowSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # ✅ naive：metadata 只保留 phenotype
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
      
      if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
        cat("  过滤后矩阵为空，跳过\n"); next
      }
      
      cat("  [预处理后] dim:", paste(dim(counts_matrix), collapse = "x"), "\n")
      
      # 过滤有效分组
      group_counts <- table(meta$phenotype)
      valid_groups <- names(group_counts[group_counts >= 2])
      if (length(valid_groups) < 2) next
      
      meta <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
      counts_matrix <- counts_matrix[, meta$sample_id, drop = FALSE]
      counts_matrix <- ensure_matrix(counts_matrix)
      
      # 枚举 pair
      all_groups    <- sort(unique(as.character(meta$phenotype)))
      disease_pairs <- combn(all_groups, 2, simplify = FALSE)
      n_pairs_total <- length(disease_pairs)
      
      scenario_idx <- match(scenario, scenarios_to_run)
      rep_idx      <- as.integer(rep)
      
      comparison_set_label <- paste(
        vapply(disease_pairs, \(p) paste0(p[1], "_vs_", p[2]), ""),
        collapse = ";"
      )
      
      for (i in seq_along(disease_pairs)) {
        
        pair      <- disease_pairs[[i]]
        g1        <- pair[1]
        g2        <- pair[2]
        comp_name <- paste0(g1, "_vs_", g2)
        
        cat("  比较:", comp_name, "\n"); flush.console()
        
        # pair 子集
        samples_all <- meta$sample_id[meta$phenotype %in% pair]
        if (length(samples_all) < 4) {
          cat("    样本数不足，跳过\n"); next
        }
        
        counts_all <- counts_matrix[, samples_all, drop = FALSE]
        counts_all <- ensure_matrix(counts_all)
        counts_all <- counts_all[rowSums(counts_all) > 0, , drop = FALSE]
        counts_all <- ensure_matrix(counts_all)
        if (nrow(counts_all) == 0) next
        
        meta_all <- meta[meta$sample_id %in% samples_all, , drop = FALSE]
        
        # ── Step 1：全量（该pair）──────────────────────────
        S_all <- tryCatch(
          run_ancombc2_core(
            counts_matrix = counts_all,
            metadata_df   = meta_all,
            pair_levels   = pair,
            qval_cutoff   = QVAL_THR,
            top_k         = TOP_K,
            debug         = FALSE
          ),
          error = function(e) {
            cat("    全量ANCOM-BC2失败:", e$message, "\n")
            character(0)
          }
        )
        cat("    全量特征数:", length(S_all), "\n"); flush.console()
        
        if (length(S_all) == 0) {
          jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
            method            = "ANCOM-BC2_naive",
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
        
        # ── Step 2：Bootstrap（pair内按 phenotype 分层抽样）──
        jac_vals <- rep(NA_real_, B)
        
        # ✅ seed 公式（与 maaslin2 naive 一致）
        base_seed <- BASE_SEED +
          1000000L * as.integer(scenario_idx) +
          10000L   * as.integer(rep_idx) +
          100L     * as.integer(i)
        
        for (b in seq_len(B)) {
          
          set.seed(base_seed + b)
          
          sub_ids <- meta_all %>%
            dplyr::group_by(phenotype) %>%
            dplyr::slice_sample(prop = SUBSAMPLE) %>%
            dplyr::ungroup() %>%
            dplyr::pull(sample_id) %>%
            as.character()
          
          if (length(sub_ids) < 4) next
          
          counts_sub <- counts_all[, sub_ids, drop = FALSE]
          counts_sub <- ensure_matrix(counts_sub)
          counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
          counts_sub <- ensure_matrix(counts_sub)
          if (nrow(counts_sub) == 0 || ncol(counts_sub) == 0) next
          
          meta_sub <- meta_all[meta_all$sample_id %in% sub_ids, , drop = FALSE]
          if (dplyr::n_distinct(meta_sub$phenotype) < 2) next
          
          # ✅ 修复：使用正确的函数名 run_ancombc2_core（已传递 pair_levels）
          S_b <- tryCatch(
            run_ancombc2_core(
              counts_matrix = counts_sub,
              metadata_df   = meta_sub,
              pair_levels   = pair,
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
          method            = "ANCOM-BC2_naive",
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
          "    [%s | %s | rep%s] Jaccard均值=%.3f (SD=%.3f, 有效boots=%d/%d)\n",
          scenario, comp_name, rep,
          if (n_boots_valid == 0) NA else mean(jac_clean),
          if (n_boots_valid == 0) NA else sd(jac_clean),
          n_boots_valid, B
        ))
        flush.console()
      }
    }
  }
  
  # ── 汇总 & 保存 ────────────────────────────────────────────
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    empty <- empty_jaccard_tibble()
    saveRDS(empty, "benchmark/results/ancombc2_naive_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/ancombc2_naive_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/ancombc2_naive_jaccard_summary.rds")
    return(list(
      jaccard_raw = empty,
      jaccard_pair_summary = empty,
      jaccard_scenario_summary = empty
    ))
  }
  
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
  # 1) Pair级汇总
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
      n_pairs        = dplyr::n_distinct(comparison),
      top_k_total    = top_k_per_pair * n_pairs,
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      .groups = "drop"
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
  
  cat("\n── ANCOM-BC2 Naive Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  
  cat("\n── ANCOM-BC2 Naive Jaccard 稳定性汇总（场景级：3个pair取平均） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,              "benchmark/results/ancombc2_naive_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/ancombc2_naive_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/ancombc2_naive_jaccard_summary.rds")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# 执行
jaccard_ancombc2 <- run_ancombc2_naive_jaccard()




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

library(ANCOMBC)
library(phyloseq)
library(tidyverse)

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

# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
# core：给定 counts_matrix + meta_df 跑 ANCOM-BC2（单个 pair）
# ✅ 返回显著特征名（只取 phenotype 效应）
# ══════════════════════════════════════════════════════════════
run_ancombc2_core_adjusted <- function(counts_matrix,
                                       metadata_df,
                                       pair_levels,
                                       qval_cutoff = 0.05,
                                       top_k = NULL,
                                       debug = FALSE) {
  
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # 统一处理 metadata
  metadata_df <- as.data.frame(metadata_df)
  rownames(metadata_df) <- NULL
  
  if (!"sample_id" %in% colnames(metadata_df)) {
    metadata_df <- tibble::rownames_to_column(metadata_df, "sample_id")
  }
  
  # adjusted：必须有 batch
  if (!"batch" %in% colnames(metadata_df)) {
    if (debug) cat("  [core-adjusted] metadata 缺少 batch，返回空\n")
    return(character(0))
  }
  
  metadata_df <- metadata_df %>%
    dplyr::select(sample_id, phenotype, batch) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      phenotype = as.character(phenotype),
      batch     = as.character(batch)
    )
  
  # 对齐样本
  common <- intersect(colnames(counts_matrix), metadata_df$sample_id)
  if (length(common) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, common, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  metadata_df   <- metadata_df[match(common, metadata_df$sample_id), , drop = FALSE]
  
  # 只保留 pair
  metadata_df <- metadata_df[metadata_df$phenotype %in% pair_levels, , drop = FALSE]
  if (nrow(metadata_df) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, metadata_df$sample_id, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # 至少两组
  if (dplyr::n_distinct(metadata_df$phenotype) < 2) return(character(0))
  
  # 去全 0 特征
  counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  # 设置参考组
  ref_group <- if ("Healthy" %in% pair_levels) "Healthy" else sort(pair_levels)[1]
  metadata_df$phenotype <- stats::relevel(
    factor(metadata_df$phenotype, levels = pair_levels),
    ref = ref_group
  )
  
  metadata_df$batch <- factor(metadata_df$batch)
  
  # 构建 phyloseq
  sample_data_df <- as.data.frame(metadata_df)
  rownames(sample_data_df) <- NULL
  sample_data_df <- tibble::column_to_rownames(sample_data_df, "sample_id")
  
  physeq <- tryCatch(
    phyloseq::phyloseq(
      phyloseq::otu_table(counts_matrix, taxa_are_rows = TRUE),
      phyloseq::sample_data(sample_data_df)
    ),
    error = function(e) NULL
  )
  if (is.null(physeq)) return(character(0))
  
  # fix_formula：phenotype + batch
  fix_formula_use <- if (dplyr::n_distinct(metadata_df$batch) >= 2) {
    "phenotype + batch"
  } else {
    if (debug) cat("  [core-adjusted] batch 仅1个水平，回退 fix_formula = phenotype\n")
    "phenotype"
  }
  
  fit <- tryCatch(
    ANCOMBC::ancombc2(
      data         = physeq,
      fix_formula  = fix_formula_use,
      p_adj_method = "BH",
      prv_cut      = 0,
      lib_cut      = 0,
      s0_perc      = 0.05,
      struc_zero   = FALSE,
      neg_lb       = FALSE,
      global       = FALSE,
      pairwise     = FALSE,
      dunnet       = FALSE,
      trend        = FALSE,
      alpha        = 0.05,
      n_cl         = 1,
      verbose      = FALSE
    ),
    error = function(e) {
      if (debug) cat("  [core-adjusted] ANCOM-BC2失败:", e$message, "\n")
      NULL
    }
  )
  if (is.null(fit)) return(character(0))
  
  res_df <- fit$res
  if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
  
  # phenotype 的非参考水平（两组比较只有1列）
  nonref <- setdiff(pair_levels, ref_group)
  if (length(nonref) != 1) return(character(0))
  
  pheno_col <- paste0("phenotype", nonref)
  q_col     <- paste0("q_",   pheno_col)
  lfc_col   <- paste0("lfc_", pheno_col)
  
  if (!all(c("taxon", q_col, lfc_col) %in% colnames(res_df))) return(character(0))
  
  feats <- res_df %>%
    dplyr::filter(!is.na(.data[[q_col]]), .data[[q_col]] < qval_cutoff) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[lfc_col]])))
  
  if (!is.null(top_k) && top_k > 0) {
    feats <- dplyr::slice_head(feats, n = top_k)
  }
  
  unique(as.character(feats$taxon))
}

# ══════════════════════════════════════════════════════════════
# 主函数：运行单个 scenario+rep 的 ANCOM-BC2 adjusted
# 输出目录：ancombc2_results_genus_adjusted
# ══════════════════════════════════════════════════════════════
run_ancombc2_adjusted_B <- function(scenario_full, rep_id) {
  
  set.seed(42 + as.integer(rep_id))
  
  cat("\n", strrep("=", 60), "\n")
  cat("场景:", scenario_full, " | Rep:", rep_id, "\n")
  cat(strrep("=", 60), "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  # 1. 读取数据
  inputs        <- generate_method_inputs(scenario_full, rep_id)
  counts_matrix <- inputs$ancombc2$counts
  metadata_df   <- inputs$ancombc2$metadata
  
  cat("特征矩阵维度:", dim(counts_matrix), "\n")
  cat("样本数:", nrow(metadata_df), "\n")
  flush.console()
  
  # 2. 数据清洗
  counts_matrix <- ensure_matrix(counts_matrix)
  counts_matrix <- round(pmax(counts_matrix, 0))
  
  # 去全 0 样本
  sample_sums <- colSums(counts_matrix)
  if (any(sample_sums == 0)) {
    cat("警告: 发现", sum(sample_sums == 0), "个全 0 样本，已删除\n")
    counts_matrix <- counts_matrix[, sample_sums > 0, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
  }
  
  # 整理 metadata：保留 batch
  if (!("batch" %in% colnames(metadata_df))) {
    stop("metadata_df 缺少 batch 列：无法进行 batch-adjusted 分析")
  }
  
  meta <- metadata_df %>%
    dplyr::select(sample_id, phenotype, batch) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      phenotype = as.character(phenotype),
      batch     = as.character(batch)
    )
  
  # 对齐样本
  common_samples <- intersect(colnames(counts_matrix), meta$sample_id)
  counts_matrix  <- counts_matrix[, common_samples, drop = FALSE]
  counts_matrix  <- ensure_matrix(counts_matrix)
  meta           <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
  
  cat("对齐后样本数:", ncol(counts_matrix), "\n")
  cat("分组分布:\n")
  print(table(meta$phenotype))
  cat("batch分布:\n")
  print(table(meta$batch))
  cat("phenotype x batch:\n")
  print(table(meta$phenotype, meta$batch))
  flush.console()
  
  # 3. 枚举所有两两比较对
  all_groups    <- sort(unique(as.character(meta$phenotype)))
  disease_pairs <- combn(all_groups, 2, simplify = FALSE)
  
  cat("将进行", length(disease_pairs), "个两两比较\n\n")
  flush.console()
  
  # 输出目录
  base_result_dir <- file.path(
    "benchmark", "ancombc2", "ancombc2_results_genus_adjusted",
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
    
    pair_idx <- meta$phenotype %in% pair
    if (sum(pair_idx) < 4) {
      cat("  样本数不足，跳过\n"); next
    }
    
    pair_counts <- counts_matrix[, pair_idx, drop = FALSE]
    pair_counts <- ensure_matrix(pair_counts)
    pair_meta   <- meta[pair_idx, , drop = FALSE]
    
    # 去全 0 行
    nonzero_taxa <- rowSums(pair_counts) > 0
    pair_counts  <- pair_counts[nonzero_taxa, , drop = FALSE]
    pair_counts  <- ensure_matrix(pair_counts)
    
    cat("  样本数:", ncol(pair_counts),
        "| 特征数:", nrow(pair_counts), "\n")
    cat("  分组分布:\n")
    print(table(pair_meta$phenotype))
    cat("  batch分布:\n")
    print(table(pair_meta$batch))
    flush.console()
    
    ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
    pair_meta$phenotype <- relevel(
      factor(pair_meta$phenotype, levels = pair), ref = ref_group
    )
    pair_meta$batch <- factor(pair_meta$batch)
    
    cat("  参考组:", ref_group, "\n")
    flush.console()
    
    # 构建 physeq
    sample_data_df <- tibble::column_to_rownames(
      as.data.frame(pair_meta), "sample_id"
    )
    physeq <- tryCatch(
      phyloseq::phyloseq(
        phyloseq::otu_table(pair_counts, taxa_are_rows = TRUE),
        phyloseq::sample_data(sample_data_df)
      ),
      error = function(e) {
        cat("  phyloseq 构建失败:", e$message, "\n"); NULL
      }
    )
    if (is.null(physeq)) next
    
    cat("  phyloseq 对象构建完成\n")
    flush.console()
    
    # 纯模型拟合计时
    t_compute_start <- proc.time()["elapsed"]
    
    ancombc_res <- tryCatch(
      ANCOMBC::ancombc2(
        data         = physeq,
        fix_formula  = "phenotype + batch",
        p_adj_method = "BH",
        prv_cut      = 0,
        lib_cut      = 0,
        s0_perc      = 0.05,
        struc_zero   = FALSE,
        neg_lb       = FALSE,
        global       = FALSE,
        pairwise     = FALSE,
        dunnet       = FALSE,
        trend        = FALSE,
        alpha        = 0.05,
        n_cl         = 1,
        verbose      = FALSE
      ),
      error = function(e) {
        cat("  ANCOM-BC2 失败:", e$message, "\n"); NULL
      }
    )
    
    t_compute_end <- proc.time()["elapsed"]
    compute_sec   <- compute_sec + as.numeric(t_compute_end - t_compute_start)
    
    if (is.null(ancombc_res)) next
    
    res_df <- ancombc_res$res
    if (is.null(res_df) || nrow(res_df) == 0) {
      cat("  无有效结果，跳过\n"); next
    }
    
    pheno_col <- paste0("phenotype", setdiff(pair, ref_group))
    
    result_df <- data.frame(
      taxon      = res_df$taxon,
      lfc        = res_df[[paste0("lfc_",  pheno_col)]],
      se         = res_df[[paste0("se_",   pheno_col)]],
      W          = res_df[[paste0("W_",    pheno_col)]],
      p_val      = res_df[[paste0("p_",    pheno_col)]],
      q_val      = res_df[[paste0("q_",    pheno_col)]],
      diff_abn   = res_df[[paste0("diff_", pheno_col)]],
      comparison = comp_name,
      group_a    = g1,
      group_b    = g2,
      ref_group  = ref_group,
      scenario   = scenario_full,
      replicate  = rep_id,
      stringsAsFactors = FALSE
    )
    
    cat("  发现", nrow(result_df), "个特征结果，",
        sum(result_df$q_val < 0.05, na.rm = TRUE), "个 FDR<0.05\n")
    flush.console()
    
    # I/O
    comp_output_dir <- file.path(base_result_dir, comp_name)
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
# 批量运行（场景B，rep1~3）
# ══════════════════════════════════════════════════════════════
run_all_adjusted_B <- function() {
  
  cat("\n", strrep("=", 60), "\n")
  cat("ANCOM-BC2 Adjusted - 场景 B (batch confound)\n")
  cat(strrep("=", 60), "\n\n")
  
  results_all    <- list()
  timing_records <- list()
  
  for (rep in reps) {
    
    key <- paste0(SCENARIO_B, "_rep", rep)
    
    out <- tryCatch(
      run_ancombc2_adjusted_B(SCENARIO_B, rep),
      error = function(e) {
        cat("Rep", rep, "失败:", e$message, "\n"); NULL
      }
    )
    
    if (!is.null(out)) {
      results_all[[key]] <- out
      timing_records[[length(timing_records) + 1]] <- tibble(
        method      = "ANCOM-BC2_adjusted",
        scenario    = SCENARIO_B,
        rep         = as.integer(rep),
        wall_sec    = out$wall_sec,
        compute_sec = out$compute_sec
      )
    }
  }
  
  cat("\n", strrep("=", 60), "\n")
  cat("全部 Rep 完成，成功:",
      sum(!sapply(results_all, is.null)), "/", length(reps), "\n")
  
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
  
  cat("\n── 运行时间汇总 ──\n")
  print(timing_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(timing_df,      "benchmark/results/ancombc2_adjusted_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/ancombc2_adjusted_timing_summary.rds")
  
  list(
    results        = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary
  )
}

# ══════════════════════════════════════════════════════════════
# 汇总结果 CSV（ANCOM-BC2 Adjusted 版本）
# ══════════════════════════════════════════════════════════════
collect_ancombc2_adjusted_results <- function(
    base_dir  = "benchmark/ancombc2/ancombc2_results_genus_adjusted",
    scenarios, reps) {
  
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
  
  if (length(all_rows) == 0) { cat("未找到任何结果文件\n"); return(NULL) }
  bind_rows(all_rows)
}

# 注意：scenarios_to_run 应该只包含 "B_batch_confound"
ancombc2_adjusted_summary_raw <- collect_ancombc2_adjusted_results(
  scenarios = scenarios_to_run,  # 此时 scenarios_to_run = c("B_batch_confound")
  reps      = reps
)

cat("汇总行数:", nrow(ancombc2_adjusted_summary_raw), "\n")
print(head(ancombc2_adjusted_summary_raw))

dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(ancombc2_adjusted_summary_raw, "benchmark/results/ancombc2_adjusted_summary_raw.rds")
write_csv(ancombc2_adjusted_summary_raw, "benchmark/results/ancombc2_adjusted_summary_raw.csv")
cat("\n汇总结果已保存\n")

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估（Adjusted版）
# 参考代码④的架构：run_one_pair_adjusted 内部函数 + 分层抽样
# ══════════════════════════════════════════════════════════════
run_ancombc2_adjusted_jaccard <- function() {
  
  BASE_SEED <- 42L
  SCENARIO  <- "B_batch_confound"
  
  jaccard_records <- list()
  
  cat("\n", strrep("=", 62), "\n")
  cat("ANCOM-BC2 Adjusted Jaccard 稳定性评估 — 场景 B（pairwise + pair均值汇总）\n")
  cat(strrep("=", 62), "\n")
  
  # ── 内部小工具：对“一个pair的数据”跑一次 batch-adjusted ANCOM-BC2 并取特征 ──
  run_one_pair_adjusted <- function(pair_counts, pair_meta, pair_levels, top_k) {
    
    # 去全0行
    pair_counts <- pair_counts[rowSums(pair_counts) > 0, , drop = FALSE]
    if (nrow(pair_counts) == 0) return(character(0))
    
    # 设置参考组
    ref_group <- if ("Healthy" %in% pair_levels) "Healthy" else sort(pair_levels)[1]
    
    # ✅ 保留 sample_id，先转 data.frame
    pair_meta <- as.data.frame(pair_meta)
    pair_meta <- pair_meta[, c("sample_id", "phenotype", "batch"), drop = FALSE]
    
    # ✅ 用 sample_id 设置 rownames，再删除该列
    rownames(pair_meta) <- pair_meta$sample_id
    pair_meta$sample_id <- NULL
    
    # 设置因子
    pair_meta$phenotype <- stats::relevel(
      factor(pair_meta$phenotype, levels = pair_levels),
      ref = ref_group
    )
    pair_meta$batch <- factor(pair_meta$batch)
    
    # ✅ 对齐 counts 列与 meta 行
    common_ids <- intersect(colnames(pair_counts), rownames(pair_meta))
    if (length(common_ids) == 0) return(character(0))
    
    pair_counts <- pair_counts[, common_ids, drop = FALSE]
    pair_meta   <- pair_meta[common_ids, , drop = FALSE]
    
    physeq <- tryCatch(
      phyloseq::phyloseq(
        phyloseq::otu_table(pair_counts, taxa_are_rows = TRUE),
        phyloseq::sample_data(pair_meta)
      ),
      error = function(e) NULL
    )
    if (is.null(physeq)) return(character(0))
    
    # fix_formula
    fix_formula_use <- if (dplyr::n_distinct(pair_meta$batch) >= 2) {
      "phenotype + batch"
    } else {
      "phenotype"
    }
    
    fit <- tryCatch(
      ANCOMBC::ancombc2(
        data         = physeq,
        fix_formula  = fix_formula_use,
        p_adj_method = "BH",
        prv_cut      = 0,
        lib_cut      = 0,
        s0_perc      = 0.05,
        struc_zero   = FALSE,
        neg_lb       = FALSE,
        global       = FALSE,
        pairwise     = FALSE,
        dunnet       = FALSE,
        trend        = FALSE,
        alpha        = 0.05,
        n_cl         = 1,
        verbose      = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(character(0))
    
    res_df <- fit$res
    if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
    
    nonref <- setdiff(pair_levels, ref_group)
    if (length(nonref) != 1) return(character(0))
    
    pheno_col <- paste0("phenotype", nonref)
    q_col     <- paste0("q_", pheno_col)
    lfc_col   <- paste0("lfc_", pheno_col)
    
    if (!all(c("taxon", q_col, lfc_col) %in% colnames(res_df))) return(character(0))
    
    feats <- res_df %>%
      dplyr::filter(!is.na(.data[[q_col]]), .data[[q_col]] < QVAL_THR) %>%
      dplyr::arrange(dplyr::desc(abs(.data[[lfc_col]])))
    
    if (!is.null(top_k) && top_k > 0) {
      feats <- dplyr::slice_head(feats, n = top_k)
    }
    
    unique(as.character(feats$taxon))
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
    
    counts_matrix <- inputs$ancombc2$counts
    metadata_df   <- inputs$ancombc2$metadata
    
    # ── 硬防护 + 基本清洗 ─────────────────────────────────────
    counts_matrix <- ensure_matrix(counts_matrix)
    counts_matrix <- round(pmax(counts_matrix, 0))
    
    col_keep      <- colSums(counts_matrix) > 0
    counts_matrix <- counts_matrix[, col_keep, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
    
    row_keep      <- rowSums(counts_matrix) > 0
    counts_matrix <- counts_matrix[row_keep, , drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
    
    if (!("batch" %in% colnames(metadata_df))) {
      cat("  metadata 缺少 batch 列，跳过\n"); next
    }
    
    meta <- metadata_df %>%
      dplyr::select(sample_id, phenotype, batch) %>%
      dplyr::mutate(
        sample_id = as.character(sample_id),
        phenotype = as.character(phenotype),
        batch     = as.character(batch)
      )
    
    # 对齐样本
    common_samples <- intersect(colnames(counts_matrix), meta$sample_id)
    counts_matrix  <- counts_matrix[, common_samples, drop = FALSE]
    counts_matrix  <- ensure_matrix(counts_matrix)
    meta           <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
    
    if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
      cat("  过滤后矩阵为空，跳过\n"); next
    }
    
    cat("  [预处理后] counts_matrix dim:",
        paste(dim(counts_matrix), collapse = "x"), "\n")
    
    # ── 过滤有效分组 ─────────────────────────────────────────
    group_counts <- table(meta$phenotype)
    valid_groups <- names(group_counts[group_counts >= 2])
    if (length(valid_groups) < 2) next
    
    meta          <- meta[meta$phenotype %in% valid_groups, , drop = FALSE]
    counts_matrix <- counts_matrix[, meta$sample_id, drop = FALSE]
    counts_matrix <- ensure_matrix(counts_matrix)
    
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
    scenario_idx <- 1L  # 场景 B 的索引
    
    for (i in seq_along(disease_pairs)) {
      
      pair      <- disease_pairs[[i]]
      g1        <- pair[1]
      g2        <- pair[2]
      comp_name <- paste0(g1, "_vs_", g2)
      
      # seed：rep + pair + boot
      pair_seed <- BASE_SEED + 1000000L * scenario_idx + 10000L * rep_idx + 100L * i
      
      cat("  比较:", comp_name, "\n"); flush.console()
      
      samples_all <- meta$sample_id[meta$phenotype %in% pair]
      if (length(samples_all) < 4) {
        cat("    样本数不足，跳过\n"); next
      }
      
      counts_all <- counts_matrix[, samples_all, drop = FALSE]
      counts_all <- ensure_matrix(counts_all)
      meta_all   <- meta[meta$sample_id %in% samples_all, , drop = FALSE]
      
      # ── ① 全量（该pair）→ S_all ───────────────────────────
      cat("    计算该pair全量参考集 S_all...\n"); flush.console()
      S_all <- tryCatch(
        run_one_pair_adjusted(
          pair_counts = counts_all,
          pair_meta   = meta_all,
          pair_levels = pair,
          top_k       = TOP_K
        ),
        error = function(e) character(0)
      )
      cat("    S_all 特征数:", length(S_all), "\n"); flush.console()
      
      if (length(S_all) == 0) {
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method            = "ANCOM-BC2_adjusted",
          scenario          = SCENARIO,
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
      
      # ── ② Bootstrap：pair内 batch × phenotype 分层抽样 ─────
      jac_vals <- rep(NA_real_, B)
      
      # 准备分层抽样用的数据框
      # ✅ 直接用 meta_all，sample_id 已是普通列，不需要 rownames 操作
      meta_all_df <- as.data.frame(meta_all)
      
      for (b in seq_len(B)) {
        
        set.seed(pair_seed + b)
        
        sub_ids <- meta_all_df %>%
          dplyr::group_by(batch, phenotype) %>%
          dplyr::slice_sample(prop = SUBSAMPLE) %>%
          dplyr::ungroup() %>%
          dplyr::pull(sample_id) %>%
          as.character()
        
        if (length(sub_ids) < 4) next
        
        counts_sub <- counts_all[, sub_ids, drop = FALSE]
        counts_sub <- ensure_matrix(counts_sub)
        counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
        counts_sub <- ensure_matrix(counts_sub)
        if (nrow(counts_sub) == 0 || ncol(counts_sub) == 0) next
        
        # ✅ filter 而非 rownames 索引
        meta_sub <- meta_all_df %>% 
          dplyr::filter(sample_id %in% sub_ids)
        
        # 检查可估性：两组都在 + 至少两个batch
        if (dplyr::n_distinct(meta_sub$phenotype) < 2) next
        if (dplyr::n_distinct(meta_sub$batch) < 2) next
        
        S_b <- tryCatch(
          run_one_pair_adjusted(
            pair_counts = counts_sub,
            pair_meta   = meta_sub,
            pair_levels = pair,
            top_k       = TOP_K
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
        method            = "ANCOM-BC2_adjusted",
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
    saveRDS(empty, "benchmark/results/ancombc2_adjusted_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/ancombc2_adjusted_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/ancombc2_adjusted_jaccard_summary.rds")
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
    ) %>%
    dplyr::mutate(
      top_k_total = top_k_per_pair * n_pairs
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
  
  cat("\n── ANCOM-BC2 Adjusted Jaccard 稳定性汇总（Pair级） ──\n")
  print(jaccard_pair_summary)
  
  cat("\n── ANCOM-BC2 Adjusted Jaccard 稳定性汇总（场景级：3个pair取平均） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,              "benchmark/results/ancombc2_adjusted_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/ancombc2_adjusted_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/ancombc2_adjusted_jaccard_summary.rds")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# ══════════════════════════════════════════════════════════════
# 执行入口
# ══════════════════════════════════════════════════════════════
out_ancombc2_adjusted <- run_all_adjusted_B()
jaccard_ancombc2_adjusted <- run_ancombc2_adjusted_jaccard()