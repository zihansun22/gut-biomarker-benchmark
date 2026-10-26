
# ══════════════════════════════════════════════════════════════
# LEfSe (Naive) - 修改版
# 修改①: BH 校正范围 = 该次运行被检验的全部 taxa（m 正确）
# 修改②: 断点续跑（主分析 + Jaccard 均支持）
# ══════════════════════════════════════════════════════════════

# ── Thread lock ───────────────────────────────────────────────
NTHREAD <- 1L
options(mc.cores = NTHREAD)
Sys.setenv(
  OMP_NUM_THREADS        = NTHREAD,
  OPENBLAS_NUM_THREADS   = NTHREAD,
  MKL_NUM_THREADS        = NTHREAD,
  VECLIB_MAXIMUM_THREADS = NTHREAD,
  NUMEXPR_NUM_THREADS    = NTHREAD
)
if (requireNamespace("data.table",   quietly = TRUE)) data.table::setDTthreads(NTHREAD)
if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::register(BiocParallel::SerialParam())

# ── Libraries ─────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(microbiomeMarker)
})

# ── Global parameters ─────────────────────────────────────────
scenarios_to_run <- c(
  "A_no_confound",
  "B_batch_confound",
  "C_sparse_small",
  "D_large_effect"
)
reps      <- as.character(1:3)
QVAL_THR  <- 0.05
LDA_THR   <- 2
TOP_K     <- 50
B         <- 30
SUBSAMPLE <- 0.8
BASE_SEED <- 42L

TARGET_PAIRS <- list(
  c("Disease_mild",  "Disease_severe"),
  c("Healthy",       "Disease_mild"),
  c("Healthy",       "Disease_severe")
)

# ── 断点续跑：checkpoint 文件路径 ────────────────────────────
# 主分析：每个 (scenario, rep) 完成后写一条记录
CKPT_MAIN   <- "benchmark/results/lefse_naive_checkpoint_main.rds"
# Jaccard：每个 (scenario, rep, comparison, boot) 完成后写
CKPT_JAC    <- "benchmark/results/lefse_naive_checkpoint_jaccard.rds"

# ── Checkpoint 读写工具 ───────────────────────────────────────

# 读取 checkpoint（不存在则返回空 tibble）
load_checkpoint <- function(path, empty_fn) {
  if (file.exists(path)) {
    ckpt <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(ckpt)) {
      cat("[Checkpoint] 读取:", path, "->", nrow(ckpt), "条记录\n")
      return(ckpt)
    }
  }
  empty_fn()
}

# 向 checkpoint 追加一行并保存
append_checkpoint <- function(path, new_row, empty_fn) {
  ckpt <- load_checkpoint(path, empty_fn)
  ckpt <- dplyr::bind_rows(ckpt, new_row)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ckpt, path)
  invisible(ckpt)
}

# 主分析 checkpoint 的空结构
empty_ckpt_main <- function() {
  tibble(
    scenario  = character(0),
    rep       = integer(0),
    n_sig     = integer(0),
    done      = logical(0)
  )
}

# Jaccard checkpoint 的空结构（每个 boot 一行）
empty_ckpt_jac <- function() {
  tibble(
    scenario    = character(0),
    rep         = integer(0),
    comparison  = character(0),
    boot        = integer(0),
    jac_val     = numeric(0),
    n_full      = integer(0),
    m_tested    = integer(0)
  )
}

# ── Utils ─────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

as_int_rep <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  d <- stringr::str_extract(as.character(x), "\\d+")
  if (is.na(d)) stop("rep_id 无法解析: ", x)
  as.integer(d)
}

align_groups <- function(groups, sample_names) {
  if (!is.null(names(groups))) {
    return(groups[sample_names])
  } else {
    if (length(groups) != length(sample_names))
      stop("align_groups: 长度不一致")
    names(groups) <- sample_names
    groups
  }
}

jaccard <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
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
    n_features_full   = integer(0),
    m_tested          = integer(0)
  )
}

# ── Taxonomy parsing ──────────────────────────────────────────
split_tax_string <- function(x) {
  if (stringr::str_detect(x, ";"))   return(stringr::str_split(x, ";",   simplify = TRUE)[1, ])
  if (stringr::str_detect(x, "\\|")) return(stringr::str_split(x, "\\|", simplify = TRUE)[1, ])
  if (stringr::str_detect(x, "\\.")) return(stringr::str_split(x, "\\.", simplify = TRUE)[1, ])
  return(x)
}

parse_taxonomy2 <- function(taxa_strings) {
  ranks  <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_df <- matrix("", nrow = length(taxa_strings), ncol = length(ranks),
                   dimnames = list(taxa_strings, ranks))
  for (i in seq_along(taxa_strings)) {
    parts <- as.character(split_tax_string(taxa_strings[i]))
    parts <- stringr::str_trim(parts[!is.na(parts) & parts != ""])
    clean <- purrr::map_chr(seq_along(parts), function(j) {
      p <- stringr::str_replace(parts[j], "^[a-z]__", "")
      if (is.na(p) || p == "" || p == "__")
        p <- paste0("Unknown_", ranks[pmin(j, 7)], "_", i)
      p
    })
    for (j in seq_len(min(length(clean), length(ranks)))) tax_df[i, j] <- clean[j]
    if (length(clean) < length(ranks))
      for (j in (length(clean) + 1):length(ranks))
        tax_df[i, j] <- paste0("Unknown_", ranks[j], "_", i)
  }
  as.data.frame(tax_df, stringsAsFactors = FALSE)
}

# ── Taxa 映射 ─────────────────────────────────────────────────
build_taxa_lookup_from_ps <- function(ps) {
  tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  original_names <- rownames(tt)
  
  rank_prefixes <- c(
    Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__",
    Family = "f__", Genus = "g__", Species = "s__"
  )
  avail_ranks <- intersect(names(rank_prefixes), colnames(tt))
  lookup <- list()
  
  for (i in seq_along(original_names)) {
    orig <- original_names[i]
    row  <- as.character(tt[i, avail_ranks, drop = TRUE])
    for (depth in seq_along(avail_ranks)) {
      parts <- paste0(rank_prefixes[avail_ranks[seq_len(depth)]],
                      row[seq_len(depth)])
      key           <- paste(parts, collapse = "|")
      lookup[[key]] <- orig
    }
  }
  
  keys   <- names(lookup)
  values <- unlist(lookup, use.names = FALSE)
  setNames(values, keys)
}

restore_taxa_name <- function(lefse_feature_names, taxa_lookup) {
  result  <- taxa_lookup[lefse_feature_names]
  is_miss <- is.na(result)
  if (any(is_miss)) {
    warning(sum(is_miss), " 个 taxa 映射失败，保留原始格式", call. = FALSE)
    result[is_miss] <- lefse_feature_names[is_miss]
  }
  unname(result)
}

extract_taxa_short_g <- function(taxa_full) {
  g <- stringr::str_trim(stringr::str_extract(taxa_full, "g__[^;|\\. ]+"))
  ifelse(is.na(g) | g == "g__", NA_character_, g)
}

# ══════════════════════════════════════════════════════════════
# ① 全 taxa BH 校正核心辅助函数
#
# 输入：
#   marker_df   : 含 taxa_full + pval 列的 data.frame（LEfSe marker_table 输出）
#   tested_taxa : 本次运行被检验的全部 taxa（去全0后的 rownames）
#   q_thr       : 显著性阈值
#
# 逻辑：
#   - 未出现在 marker_df 的 taxa → 补 p=1（不显著占位）
#   - p.adjust(p_all, "BH") 的长度 = length(tested_taxa) = 真实 m
#   - 取 marker_df 中每个 taxa 对应的 q 值
# ══════════════════════════════════════════════════════════════
bh_fullm <- function(marker_taxa, marker_pval, tested_taxa, q_thr = QVAL_THR) {
  # marker_taxa : character vector，marker_df 中的 taxa_full
  # marker_pval : numeric vector，与 marker_taxa 等长
  # tested_taxa : 本次运行全部 tested taxa（定义 m）
  # 返回：与 marker_taxa 等长的 q 值 vector
  
  tested_taxa <- unique(as.character(tested_taxa))
  if (length(tested_taxa) == 0) return(rep(NA_real_, length(marker_taxa)))
  
  # 建 p 全集（默认 1.0）
  p_full <- rep(1.0, length(tested_taxa))
  names(p_full) <- tested_taxa
  
  # 覆盖 marker 中有 p 值的 taxa（重复取 min p）
  valid <- !is.na(marker_pval) & marker_taxa %in% tested_taxa
  if (any(valid)) {
    p_min <- tapply(
      marker_pval[valid],
      marker_taxa[valid],
      function(x) min(suppressWarnings(as.numeric(x)), na.rm = TRUE)
    )
    hit <- intersect(names(p_min), names(p_full))
    p_full[hit] <- p_min[hit]
  }
  
  # BH 校正，m = length(tested_taxa)
  q_full <- p.adjust(p_full, method = "BH")
  
  # 返回 marker_taxa 对应的 q（不在 tested_taxa 内的返回 NA）
  q_out <- q_full[as.character(marker_taxa)]
  unname(q_out)
}

# ══════════════════════════════════════════════════════════════
# Jaccard 核心：单次运行 LEfSe，返回显著 taxa 名
# ① BH 使用全 taxa 作为 m
# ══════════════════════════════════════════════════════════════
run_lefse_core_jaccard <- function(abundance_matrix,
                                   metadata_df,
                                   taxonomy_df,
                                   top_k  = NULL,
                                   debug  = FALSE) {
  
  abundance_matrix <- ensure_matrix(abundance_matrix)
  
  common_samples <- intersect(colnames(abundance_matrix), rownames(metadata_df))
  if (length(common_samples) < 4) return(list(taxa = character(0), m_tested = 0L))
  
  abundance_matrix <- abundance_matrix[, common_samples, drop = FALSE]
  metadata_df      <- metadata_df[common_samples, , drop = FALSE]
  abundance_matrix <- ensure_matrix(abundance_matrix)
  
  group_counts <- table(metadata_df$phenotype)
  valid_groups <- names(group_counts[group_counts >= 2])
  if (length(valid_groups) < 2) return(list(taxa = character(0), m_tested = 0L))
  
  # 去全0行 → 定义 tested_taxa（BH 的 m）
  abundance_matrix <- abundance_matrix[rowSums(abundance_matrix) > 0, , drop = FALSE]
  abundance_matrix <- ensure_matrix(abundance_matrix)
  if (nrow(abundance_matrix) == 0) return(list(taxa = character(0), m_tested = 0L))
  
  tested_taxa <- rownames(abundance_matrix)          # ← m 的来源
  m_tested    <- length(tested_taxa)
  
  tax_sub <- taxonomy_df[tested_taxa, , drop = FALSE]
  
  ps <- tryCatch(
    phyloseq(
      otu_table(abundance_matrix, taxa_are_rows = TRUE),
      sample_data(metadata_df),
      tax_table(as.matrix(tax_sub))
    ),
    error = function(e) NULL
  )
  if (is.null(ps)) return(list(taxa = character(0), m_tested = m_tested))
  
  taxa_lookup <- build_taxa_lookup_from_ps(ps)
  
  lefse_res <- tryCatch(
    run_lefse(
      ps,
      wilcoxon_cutoff = 0.05,
      kw_cutoff       = 0.05,
      lda_cutoff      = 0,           # 全量输出，不预过滤
      group           = "phenotype",
      multigrp_strat  = FALSE,
      taxa_rank       = "all"
    ),
    error = function(e) {
      if (debug) cat("  [core] LEfSe 失败:", e$message, "\n")
      NULL
    }
  )
  if (is.null(lefse_res)) return(list(taxa = character(0), m_tested = m_tested))
  
  markers <- tryCatch({
    mt <- lefse_res@marker_table
    df <- data.frame(mt@.Data, stringsAsFactors = FALSE)
    names(df) <- mt@names
    df
  }, error = function(e) NULL)
  
  if (is.null(markers) || nrow(markers) == 0)
    return(list(taxa = character(0), m_tested = m_tested))
  
  get_col <- function(df, cands) {
    nm <- intersect(cands, names(df))
    if (!length(nm)) NULL else nm[1]
  }
  col_feature <- get_col(markers, c("feature","Feature","taxon","Taxon"))
  col_lda     <- get_col(markers, c("ef_lda","lda","lda_score","LDA","ef_lda_mean"))
  col_p       <- get_col(markers, c("pvalue","p_val","pvalue_kw","pvalue_wilcoxon",
                                    "p.value","p_value"))
  if (is.null(col_feature)) return(list(taxa = character(0), m_tested = m_tested))
  
  markers$feature_raw <- as.character(markers[[col_feature]])
  markers$p_raw       <- if (!is.null(col_p))
    suppressWarnings(as.numeric(markers[[col_p]])) else NA_real_
  
  # ① 全 taxa BH：m = length(tested_taxa)，未出现的 taxa 补 p=1
  # ✅ 修复（先 restore，再传给 bh_fullm）
  markers$feature_restored <- restore_taxa_name(markers$feature_raw, taxa_lookup)
  
  markers$q_val <- bh_fullm(
    marker_taxa  = markers$feature_restored,  # ← 已还原为原始分号格式
    marker_pval  = markers$p_raw,
    tested_taxa  = tested_taxa,               # ← 同为分号格式，可匹配
    q_thr        = QVAL_THR
  )
  # 后续也用 restored 名字筛选
  feats <- markers %>%
    dplyr::filter(!is.na(feature_restored), !is.na(q_val), q_val <= QVAL_THR)
  
  if (!is.null(col_lda)) {
    feats <- feats %>%
      dplyr::mutate(lda_raw = suppressWarnings(as.numeric(.data[[col_lda]]))) %>%
      dplyr::arrange(q_val, dplyr::desc(lda_raw))
  } else {
    feats <- feats %>% dplyr::arrange(q_val)
  }
  
  if (nrow(feats) == 0) return(list(taxa = character(0), m_tested = m_tested))
  if (!is.null(top_k)) feats <- dplyr::slice_head(feats, n = top_k)
  
  list(taxa = unique(feats$feature_restored), m_tested = m_tested)
}

# ══════════════════════════════════════════════════════════════
# 主分析：单个 pair
# ① BH 使用全 taxa 作为 m
# ══════════════════════════════════════════════════════════════
run_lefse_naive_one_pair <- function(abundance_matrix,
                                     metadata_df,
                                     scenario,
                                     rep_id,
                                     pair,
                                     lda_threshold = LDA_THR,
                                     debug = FALSE) {
  g1 <- pair[1]; g2 <- pair[2]
  
  keep_samples <- rownames(metadata_df)[metadata_df$phenotype %in% c(g1, g2)]
  if (length(keep_samples) < 4) return(NULL)
  
  meta    <- metadata_df[keep_samples, , drop = FALSE]
  grp_cnt <- table(meta$phenotype)
  if (sum(grp_cnt >= 2) < 2) return(NULL)
  
  ab <- abundance_matrix[, keep_samples, drop = FALSE]
  ab <- ensure_matrix(ab)
  
  # 去全0行 → tested_taxa（BH 的 m）
  ab          <- ab[rowSums(ab) > 0, , drop = FALSE]
  ab          <- ensure_matrix(ab)
  if (nrow(ab) == 0) return(NULL)
  
  tested_taxa         <- rownames(ab)      # ← BH m 的来源
  original_taxa_names <- tested_taxa
  taxonomy_df         <- parse_taxonomy2(original_taxa_names)
  meta$phenotype      <- factor(meta$phenotype, levels = c(g1, g2))
  
  ps <- tryCatch(
    phyloseq(
      otu_table(ab, taxa_are_rows = TRUE),
      sample_data(meta),
      tax_table(as.matrix(taxonomy_df))
    ),
    error = function(e) {
      if (debug) cat("  phyloseq 失败:", e$message, "\n"); NULL
    }
  )
  if (is.null(ps)) return(NULL)
  
  taxa_lookup <- build_taxa_lookup_from_ps(ps)
  
  lefse_res <- tryCatch(
    run_lefse(
      ps,
      wilcoxon_cutoff = 0.05,
      kw_cutoff       = 0.05,
      lda_cutoff      = 0,
      group           = "phenotype",
      multigrp_strat  = FALSE,
      taxa_rank       = "all"
    ),
    error = function(e) {
      if (debug) cat("  LEfSe 失败:", e$message, "\n"); NULL
    }
  )
  if (is.null(lefse_res)) return(NULL)
  
  markers <- tryCatch({
    mt <- lefse_res@marker_table
    df <- data.frame(mt@.Data, stringsAsFactors = FALSE)
    names(df) <- mt@names
    df
  }, error = function(e) NULL)
  
  if (is.null(markers) || nrow(markers) == 0) return(NULL)
  
  get_col <- function(df, cands) {
    nm <- intersect(cands, names(df)); if (!length(nm)) NULL else nm[1]
  }
  col_feature <- get_col(markers, c("feature","Feature","taxon","Taxon"))
  col_lda     <- get_col(markers, c("ef_lda","lda","lda_score","LDA","ef_lda_mean"))
  col_enrich  <- get_col(markers, c("enrich_group","enriched_group","group","class"))
  col_p       <- get_col(markers, c("pvalue","p_val","pvalue_kw","pvalue_wilcoxon",
                                    "p.value","p_value"))
  
  comp_name <- paste0(g1, "_vs_", g2)
  rep_int   <- as_int_rep(rep_id)
  
  out <- markers %>%
    dplyr::mutate(
      feature_raw = if (!is.null(col_feature)) as.character(.data[[col_feature]]) else NA_character_,
      ef_lda_raw  = if (!is.null(col_lda))     suppressWarnings(as.numeric(.data[[col_lda]])) else NA_real_,
      enrich_raw  = if (!is.null(col_enrich))  as.character(.data[[col_enrich]]) else NA_character_,
      p_raw       = if (!is.null(col_p))       suppressWarnings(as.numeric(.data[[col_p]])) else NA_real_
    ) %>%
    dplyr::filter(!is.na(feature_raw)) %>%
    dplyr::mutate(
      taxa_full_restored = restore_taxa_name(feature_raw, taxa_lookup)
    ) %>%
    dplyr::transmute(
      method       = "LEfSe_naive",
      mode         = "naive",
      scenario     = as.character(scenario),
      replicate    = rep_int,
      comparison   = comp_name,
      group_a      = g1,
      group_b      = g2,
      ref_group    = g1,
      batch        = "all",
      taxa_full    = taxa_full_restored,
      taxa_short   = extract_taxa_short_g(taxa_full_restored),
      effect       = ef_lda_raw,
      pval         = p_raw,
      qval         = NA_real_,    # 下面补充
      significant  = FALSE,       # 下面补充
      enrich_group = enrich_raw
    )
  
  # ① 全 taxa BH：m = length(tested_taxa)
  out$qval <- bh_fullm(
    marker_taxa  = out$taxa_full,
    marker_pval  = out$pval,
    tested_taxa  = tested_taxa,
    q_thr        = QVAL_THR
  )
  out$significant <- !is.na(out$qval) & out$qval <= QVAL_THR
  
  attr(out, "m_tested") <- length(tested_taxa)
  out
}

# ══════════════════════════════════════════════════════════════
# Main: one scenario + one rep
# ② 断点续跑：跳过已完成的 (scenario, rep)
# ══════════════════════════════════════════════════════════════
run_lefse_naive_one_rep <- function(scenario_full, rep_id, debug = TRUE) {
  
  rep_int <- as_int_rep(rep_id)
  set.seed(42)
  
  cat("\n", strrep("=", 60), "\n")
  cat("LEfSe NAIVE | scenario:", scenario_full,
      "| Rep:", rep_id, "(int:", rep_int, ")\n")
  cat(strrep("=", 60), "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  inputs           <- generate_method_inputs(scenario_full, rep_id)
  abundance_matrix <- ensure_matrix(inputs$lefse$abundance)
  groups           <- inputs$lefse$groups
  
  phenotype_vec <- align_groups(groups, colnames(abundance_matrix))
  
  metadata_df <- data.frame(
    sample_id = colnames(abundance_matrix),
    phenotype = as.character(phenotype_vec),
    row.names = colnames(abundance_matrix),
    stringsAsFactors = FALSE
  )
  
  common_samples   <- intersect(colnames(abundance_matrix), rownames(metadata_df))
  abundance_matrix <- abundance_matrix[, common_samples, drop = FALSE]
  metadata_df      <- metadata_df[common_samples, , drop = FALSE]
  
  cat("\n表型分布:\n"); print(table(metadata_df$phenotype)); flush.console()
  
  result_dir <- file.path("benchmark","lefse","lefse_results_genus_naive",
                          paste0(scenario_full, "_rep", rep_int))
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_groups <- sort(unique(metadata_df$phenotype))
  pairs      <- purrr::keep(TARGET_PAIRS, ~all(.x %in% all_groups))
  
  compute_sec <- 0
  results_all <- list()
  
  for (pair in pairs) {
    g1        <- pair[1]; g2 <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    
    # ② 断点续跑：检查单个 pair 的 csv 是否已存在
    pair_csv <- file.path(result_dir, paste0(comp_name, "_naive_markers.csv"))
    if (file.exists(pair_csv)) {
      cat("  [跳过] 已完成 (从磁盘读取):", comp_name, "\n")
      res_cached <- tryCatch(readr::read_csv(pair_csv, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(res_cached) && nrow(res_cached) > 0) {
        results_all[[comp_name]] <- res_cached
      }
      next
    }
    
    cat("\n", strrep("-", 50), "\n", "Comparison:", comp_name, "\n")
    flush.console()
    
    t0  <- proc.time()["elapsed"]
    res <- run_lefse_naive_one_pair(
      abundance_matrix = abundance_matrix,
      metadata_df      = metadata_df,
      scenario         = scenario_full,
      rep_id           = rep_int,
      pair             = pair,
      lda_threshold    = LDA_THR,
      debug            = debug
    )
    compute_sec <- compute_sec + as.numeric(proc.time()["elapsed"] - t0)
    
    if (!is.null(res) && nrow(res) > 0) {
      results_all[[comp_name]] <- res
      readr::write_csv(res, pair_csv)
      
      cat("  m_tested:", attr(res, "m_tested") %||% NA, "\n")
      cat("  n_total :", nrow(res), "\n")
      cat("  n_sig   :", sum(res$significant, na.rm = TRUE), "\n")
    } else {
      cat("  no markers\n")
    }
  }
  
  final_df <- dplyr::bind_rows(results_all)
  
  all_csv <- file.path(result_dir, "all_comparisons_naive_markers.csv")
  if (nrow(final_df) > 0) readr::write_csv(final_df, all_csv)
  
  wall_sec <- as.numeric(proc.time()["elapsed"] - t_wall_start)
  
  # ② 写主分析 checkpoint
  append_checkpoint(
    CKPT_MAIN,
    tibble(scenario = scenario_full, rep = rep_int,
           n_sig = sum(final_df$significant, na.rm = TRUE),
           done  = TRUE),
    empty_ckpt_main
  )
  
  cat("\n── Timing ──\nwall_sec:", round(wall_sec, 2),
      "| compute_sec:", round(compute_sec, 2), "\n")
  flush.console()
  
  list(
    results_df  = final_df,
    wall_sec    = wall_sec,
    compute_sec = compute_sec,
    n_pairs     = length(pairs)
  )
}

# ══════════════════════════════════════════════════════════════
# Run all scenarios × all reps
# ② 断点续跑：自动跳过已完成的 (scenario, rep)
# ══════════════════════════════════════════════════════════════
run_all_naive <- function() {
  
  cat("\n", strrep("=", 60), "\n")
  cat("LEfSe NAIVE - all scenarios, all reps\n")
  cat(strrep("=", 60), "\n\n")
  
  # ② 读取已完成记录
  ckpt_main <- load_checkpoint(CKPT_MAIN, empty_ckpt_main)
  done_keys <- if (nrow(ckpt_main) > 0 && "done" %in% names(ckpt_main)) {
    ckpt_main %>%
      dplyr::filter(done) %>%
      dplyr::mutate(key = paste0(scenario, "_rep", rep)) %>%
      dplyr::pull(key)
  } else character(0)
  
  cat("[断点续跑] 已完成:", length(done_keys), "个 (scenario×rep)\n")
  if (length(done_keys) > 0) cat("  ", paste(done_keys, collapse = ", "), "\n")
  
  results_all    <- list()
  timing_records <- list()
  all_results    <- list()
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      key <- paste0(scenario, "_rep", rep)
      
      # ② 跳过已完成
      if (key %in% done_keys) {
        cat("\n[跳过]", key, "\n")
        # 尝试从磁盘读取结果以合入汇总
        result_dir <- file.path("benchmark","lefse","lefse_results_genus_naive",
                                paste0(scenario,"_rep",as_int_rep(rep)))
        cached_csv <- file.path(result_dir, "all_comparisons_naive_markers.csv")
        if (file.exists(cached_csv)) {
          cached <- tryCatch(readr::read_csv(cached_csv, show_col_types = FALSE), error = function(e) NULL)
          if (!is.null(cached)) all_results[[key]] <- cached
        }
        next
      }
      
      cat("\n--- Scenario:", scenario, "| Rep:", rep, "---\n"); flush.console()
      
      out <- tryCatch(
        run_lefse_naive_one_rep(scenario_full = scenario, rep_id = rep, debug = TRUE),
        error = function(e) {
          cat("Rep", rep, "失败:", conditionMessage(e), "\n"); NULL
        }
      )
      
      if (!is.null(out)) {
        results_all[[key]] <- out
        
        timing_records[[length(timing_records) + 1]] <- tibble(
          method      = "LEfSe_naive",
          mode        = "naive",
          scenario    = scenario,
          replicate   = as_int_rep(rep),
          wall_sec    = out$wall_sec,
          compute_sec = out$compute_sec,
          n_pairs     = out$n_pairs
        )
        
        if (!is.null(out$results_df) && nrow(out$results_df) > 0)
          all_results[[key]] <- out$results_df
      }
    }
  }
  
  timing_df <- if (length(timing_records) > 0)
    dplyr::bind_rows(timing_records) else
      tibble(method=character(), mode=character(), scenario=character(),
             replicate=integer(), wall_sec=numeric(),
             compute_sec=numeric(), n_pairs=integer())
  
  timing_summary <- timing_df %>%
    dplyr::group_by(method, mode, scenario) %>%
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
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(timing_df,      "benchmark/results/lefse_naive_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/lefse_naive_timing_summary.rds")
  
  final_all <- dplyr::bind_rows(all_results)
  
  if (nrow(final_all) > 0) {
    readr::write_csv(final_all, "benchmark/results/lefse_naive_summary_raw.csv")
    saveRDS(final_all,          "benchmark/results/lefse_naive_summary_raw.rds")
  } else {
    warning("所有 scenario/rep 均无结果")
  }
  
  cat("\n── LEfSe naive timing summary ──\n")
  print(timing_summary)
  
  list(
    results_all    = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary,
    final_all      = final_all
  )
}

# ══════════════════════════════════════════════════════════════
# Jaccard 稳定性评估
# ① BH 使用全 taxa 作为 m
# ② 断点续跑：以 (scenario, rep, comparison, boot) 为粒度
# ══════════════════════════════════════════════════════════════
run_lefse_naive_jaccard <- function() {
  
  jaccard_records <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("LEfSe Naive Jaccard 稳定性评估\n")
  cat("① BH m = 全部 tested taxa  ② 断点续跑\n")
  cat("══════════════════════════════════════════════════════════════\n")
  
  # ② 读取已完成的 boot 记录
  ckpt_jac <- load_checkpoint(CKPT_JAC, empty_ckpt_jac)
  done_boots <- if (nrow(ckpt_jac) > 0) {
    ckpt_jac %>%
      dplyr::mutate(
        boot_key = paste0(scenario,"_rep",rep,"_",comparison,"_boot",boot)
      ) %>%
      dplyr::pull(boot_key)
  } else character(0)
  
  cat("[断点续跑] 已完成 boot 数:", length(done_boots), "\n\n")
  
  for (scenario in scenarios_to_run) {
    for (rep in reps) {
      
      cat("\n[Jaccard] 场景:", scenario, "Rep:", rep, "\n"); flush.console()
      
      inputs <- tryCatch(
        generate_method_inputs(scenario, rep),
        error = function(e) { cat("数据读取失败:", e$message, "\n"); NULL }
      )
      if (is.null(inputs)) next
      
      abundance_matrix <- ensure_matrix(inputs$lefse$abundance)
      groups           <- inputs$lefse$groups
      
      col_keep <- colSums(abundance_matrix) > 0
      abundance_matrix <- abundance_matrix[, col_keep, drop = FALSE]
      row_keep <- rowSums(abundance_matrix) > 0
      abundance_matrix <- abundance_matrix[row_keep, , drop = FALSE]
      abundance_matrix <- ensure_matrix(abundance_matrix)
      
      phenotype_vec <- align_groups(groups, colnames(abundance_matrix))
      
      metadata_df <- data.frame(
        sample_id = colnames(abundance_matrix),
        phenotype = as.character(phenotype_vec),
        stringsAsFactors = FALSE
      )
      rownames(metadata_df) <- metadata_df$sample_id
      
      common_samples   <- intersect(colnames(abundance_matrix), rownames(metadata_df))
      abundance_matrix <- abundance_matrix[, common_samples, drop = FALSE]
      abundance_matrix <- ensure_matrix(abundance_matrix)
      metadata_df      <- metadata_df[common_samples, , drop = FALSE]
      
      if (nrow(abundance_matrix) == 0 || ncol(abundance_matrix) == 0) next
      
      group_counts <- table(metadata_df$phenotype)
      valid_groups <- names(group_counts[group_counts >= 2])
      if (length(valid_groups) < 2) next
      
      metadata_df      <- metadata_df[metadata_df$phenotype %in% valid_groups, , drop = FALSE]
      abundance_matrix <- abundance_matrix[, rownames(metadata_df), drop = FALSE]
      abundance_matrix <- ensure_matrix(abundance_matrix)
      
      taxonomy_df_full <- parse_taxonomy2(rownames(abundance_matrix))
      all_groups       <- sort(unique(as.character(metadata_df$phenotype)))
      disease_pairs    <- combn(all_groups, 2, simplify = FALSE)
      n_pairs_total    <- length(disease_pairs)
      
      scenario_idx   <- match(scenario, scenarios_to_run)
      rep_idx        <- as.integer(rep)
      comp_set_label <- paste(
        vapply(disease_pairs, \(p) paste0(p[1],"_vs_",p[2]), ""),
        collapse = ";"
      )
      
      for (i in seq_along(disease_pairs)) {
        
        pair      <- disease_pairs[[i]]
        g1        <- pair[1]; g2 <- pair[2]
        comp_name <- paste0(g1, "_vs_", g2)
        pair_seed <- BASE_SEED + 1000000L*scenario_idx + 10000L*rep_idx + 100L*i
        
        cat("  比较:", comp_name, "\n"); flush.console()
        
        # ② 检查该 pair 的全量结果是否已有（从 checkpoint 读取）
        full_run_key <- paste0(scenario,"_rep",rep_idx,"_",comp_name,"_boot0")
        cached_full <- ckpt_jac %>%
          dplyr::filter(scenario == !!scenario, rep == rep_idx,
                        comparison == comp_name, boot == 0L)
        
        samples_all <- rownames(metadata_df)[metadata_df$phenotype %in% pair]
        if (length(samples_all) < 4) next
        
        ab_all  <- abundance_matrix[, samples_all, drop = FALSE]
        ab_all  <- ensure_matrix(ab_all)
        ab_all  <- ab_all[rowSums(ab_all) > 0, , drop = FALSE]
        ab_all  <- ensure_matrix(ab_all)
        if (nrow(ab_all) == 0 || ncol(ab_all) == 0) next
        
        meta_all           <- metadata_df[samples_all, , drop = FALSE]
        meta_all$phenotype <- factor(meta_all$phenotype)
        tax_all            <- taxonomy_df_full[rownames(ab_all), , drop = FALSE]
        
        # Step 1：全量运行（只要未缓存就跑）
        if (nrow(cached_full) > 0) {
          # 从 checkpoint 恢复全量结果
          S_all    <- character(0)   # 全量 taxa 不存在 checkpoint，需重跑
          # 注意：全量结果（S_all）不保存到 checkpoint，因为字符串集合不便
          # 因此全量每次都重跑（快，通常几秒），只有 boot 才做粒度续跑
          cat("    [全量] 重新计算 S_all（全量运行快）\n")
        }
        
        full_res <- tryCatch(
          run_lefse_core_jaccard(
            abundance_matrix = ab_all,
            metadata_df      = meta_all,
            taxonomy_df      = tax_all,
            top_k            = TOP_K,
            debug            = FALSE
          ),
          error = function(e) { cat("  全量失败:", e$message, "\n"); list(taxa=character(0), m_tested=0L) }
        )
        
        S_all    <- full_res$taxa
        m_tested <- full_res$m_tested
        
        cat("  全量特征数:", length(S_all), " | m_tested:", m_tested, "\n")
        
        if (length(S_all) == 0) {
          jaccard_records[[length(jaccard_records) + 1]] <- tibble(
            method=        "LEfSe_naive", scenario=scenario, rep=rep_idx,
            comparison=    comp_name, comparison_set=comp_set_label,
            top_k_per_pair=as.integer(TOP_K), n_pairs=as.integer(n_pairs_total),
            top_k_total=   as.integer(TOP_K*n_pairs_total),
            jaccard_mean=  NA_real_, jaccard_median=NA_real_, jaccard_sd=NA_real_,
            jaccard_min=   NA_real_, jaccard_max=NA_real_,
            n_boots=       as.integer(B), n_boots_valid=0L,
            n_features_full=0L, m_tested=as.integer(m_tested)
          )
          next
        }
        
        # Step 2：B 次 subsampling
        jac_vals <- rep(NA_real_, B)
        
        for (b in seq_len(B)) {
          
          # ② 断点续跑：检查该 boot 是否已完成
          bkey <- paste0(scenario,"_rep",rep_idx,"_",comp_name,"_boot",b)
          if (bkey %in% done_boots) {
            cached_b <- ckpt_jac %>%
              dplyr::filter(scenario==!!scenario, rep==rep_idx,
                            comparison==comp_name, boot==b)
            if (nrow(cached_b) > 0) {
              jac_vals[b] <- cached_b$jac_val[1]
              next
            }
          }
          
          set.seed(pair_seed + b)
          
          sub_ids <- meta_all %>%
            as.data.frame() %>%
            tibble::rownames_to_column("tmp_rowname") %>%
            dplyr::group_by(phenotype) %>%
            dplyr::slice_sample(prop = SUBSAMPLE) %>%
            dplyr::ungroup() %>%
            dplyr::pull(tmp_rowname) %>%
            as.character()
          
          if (length(sub_ids) < 4) next
          
          ab_sub  <- ab_all[, sub_ids, drop = FALSE]
          ab_sub  <- ensure_matrix(ab_sub)
          ab_sub  <- ab_sub[rowSums(ab_sub) > 0, , drop = FALSE]
          ab_sub  <- ensure_matrix(ab_sub)
          if (nrow(ab_sub) == 0 || ncol(ab_sub) == 0) next
          
          meta_sub           <- meta_all[sub_ids, , drop = FALSE]
          meta_sub$phenotype <- factor(meta_sub$phenotype)
          tax_sub            <- taxonomy_df_full[rownames(ab_sub), , drop = FALSE]
          
          sub_res <- tryCatch(
            run_lefse_core_jaccard(
              abundance_matrix = ab_sub,
              metadata_df      = meta_sub,
              taxonomy_df      = tax_sub,
              top_k            = TOP_K,
              debug            = FALSE
            ),
            error = function(e) list(taxa=character(0), m_tested=0L)
          )
          
          S_b         <- sub_res$taxa
          jac_b       <- jaccard(S_b, S_all)
          jac_vals[b] <- jac_b
          
          # ② 写 boot checkpoint
          append_checkpoint(
            CKPT_JAC,
            tibble(scenario=scenario, rep=as.integer(rep_idx),
                   comparison=comp_name, boot=as.integer(b),
                   jac_val=ifelse(is.na(jac_b), NA_real_, jac_b),
                   n_full=as.integer(length(S_all)),
                   m_tested=as.integer(sub_res$m_tested)),
            empty_ckpt_jac
          )
          # 更新 done_boots（当前会话内也跳过）
          done_boots <- c(done_boots, bkey)
        }
        
        jac_clean     <- jac_vals[!is.na(jac_vals)]
        n_boots_valid <- length(jac_clean)
        
        jaccard_records[[length(jaccard_records) + 1]] <- tibble(
          method=         "LEfSe_naive", scenario=scenario, rep=rep_idx,
          comparison=     comp_name, comparison_set=comp_set_label,
          top_k_per_pair= as.integer(TOP_K), n_pairs=as.integer(n_pairs_total),
          top_k_total=    as.integer(TOP_K*n_pairs_total),
          jaccard_mean=   if (n_boots_valid==0) NA_real_ else mean(jac_clean),
          jaccard_median= if (n_boots_valid==0) NA_real_ else median(jac_clean),
          jaccard_sd=     if (n_boots_valid==0) NA_real_ else sd(jac_clean),
          jaccard_min=    if (n_boots_valid==0) NA_real_ else min(jac_clean),
          jaccard_max=    if (n_boots_valid==0) NA_real_ else max(jac_clean),
          n_boots=        as.integer(B), n_boots_valid=as.integer(n_boots_valid),
          n_features_full=as.integer(length(S_all)),
          m_tested=       as.integer(m_tested)
        )
        
        cat(sprintf(
          "  Jaccard=%.3f (SD=%.3f, 有效boots=%d/%d)\n",
          if (n_boots_valid==0) NA else mean(jac_clean),
          if (n_boots_valid==0) NA else sd(jac_clean),
          n_boots_valid, B
        ))
        flush.console()
      }
    }
  }
  
  # ── 汇总 ────────────────────────────────────────────────────
  if (length(jaccard_records) == 0) {
    cat("\n警告：所有比较均无显著特征，Jaccard 结果为空\n")
    empty <- empty_jaccard_tibble()
    saveRDS(empty, "benchmark/results/lefse_naive_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/lefse_naive_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/lefse_naive_jaccard_summary.rds")
    return(list(jaccard_raw=empty, jaccard_pair_summary=empty, jaccard_scenario_summary=empty))
  }
  
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"), ~ifelse(is.nan(.), NA_real_, .)))
  
  jaccard_pair_summary <- jaccard_raw %>%
    dplyr::group_by(method, scenario, comparison, top_k_per_pair) %>%
    dplyr::summarise(
      jaccard_mean_across_reps   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_across_reps     = mean(jaccard_sd,     na.rm = TRUE),
      jaccard_min_across_reps    = min(jaccard_min,     na.rm = TRUE),
      jaccard_max_across_reps    = max(jaccard_max,     na.rm = TRUE),
      n_boots                    = max(n_boots,         na.rm = TRUE),
      n_boots_valid_mean         = mean(n_boots_valid,  na.rm = TRUE),
      m_tested_mean              = mean(m_tested,       na.rm = TRUE),
      .groups = "drop"
    )
  
  jaccard_rep_pairavg <- jaccard_raw %>%
    dplyr::group_by(method, scenario, rep, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set         = paste(sort(unique(comparison)), collapse=";"),
      n_pairs                = dplyr::n_distinct(comparison),
      jaccard_mean_pairavg   = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median_pairavg = mean(jaccard_median, na.rm = TRUE),
      m_tested_mean          = mean(m_tested,       na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(top_k_total = as.integer(top_k_per_pair * n_pairs))
  
  jaccard_scenario_summary <- jaccard_rep_pairavg %>%
    dplyr::group_by(method, scenario, top_k_per_pair) %>%
    dplyr::summarise(
      comparison_set             = paste(sort(unique(comparison_set)), collapse=" | "),
      n_pairs_mean               = mean(n_pairs,    na.rm = TRUE),
      top_k_total_mean           = mean(top_k_total, na.rm = TRUE),
      jaccard_mean_across_reps   = mean(jaccard_mean_pairavg,   na.rm = TRUE),
      jaccard_median_across_reps = mean(jaccard_median_pairavg, na.rm = TRUE),
      m_tested_mean              = mean(m_tested_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"), ~ifelse(is.nan(.), NA_real_, .)))
  
  cat("\n── Pair 级汇总 ──\n"); print(jaccard_pair_summary)
  cat("\n── 场景级汇总 ──\n"); print(jaccard_scenario_summary)
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(jaccard_raw,              "benchmark/results/lefse_naive_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/lefse_naive_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/lefse_naive_jaccard_summary.rds")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# ══════════════════════════════════════════════════════════════
# Execute
# ══════════════════════════════════════════════════════════════
out_lefse_naive    <- run_all_naive()
jaccard_lefse_naive <- run_lefse_naive_jaccard()



# ══════════════════════════════════════════════════════════════
# LEfSe (Adjusted / Stratified) - Benchmark-aligned version
# FIX: BH correction is applied to ALL tested taxa (within this analysis),
#      not only LDA-filtered taxa.
# Includes:
#   - Stratified per-batch LEfSe for each TARGET_PAIRS comparison
#   - Meta-analysis across batches (Fisher) + BH on full tested set
#   - Timing summary
#   - Jaccard bootstrap stability evaluation (pairwise, stratified resampling)
# Outputs:
#   benchmark/results/lefse_stratified_summary_raw.{csv,rds}
#   benchmark/results/lefse_stratified_batch_raw.{csv,rds}
#   benchmark/results/lefse_stratified_timing_raw.rds
#   benchmark/results/lefse_stratified_timing_summary.rds
#   benchmark/results/lefse_stratified_jaccard_*.rds
# ══════════════════════════════════════════════════════════════

# ── Thread lock (before any parallel work) ────────────────────
NTHREAD <- 1L
options(mc.cores = NTHREAD)
Sys.setenv(
  OMP_NUM_THREADS        = NTHREAD,
  OPENBLAS_NUM_THREADS   = NTHREAD,
  MKL_NUM_THREADS        = NTHREAD,
  VECLIB_MAXIMUM_THREADS = NTHREAD,
  NUMEXPR_NUM_THREADS    = NTHREAD
)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(NTHREAD)
if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::register(BiocParallel::SerialParam())

cat("Thread lock:", NTHREAD, "\n")
print(Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")))

# ── Libraries ────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(microbiomeMarker)
})

# ── Global parameters ────────────────────────────────────────
SCENARIO_B <- "B_batch_confound"
reps       <- as.character(1:3)

QVAL_THR      <- 0.05
LDA_THR       <- 2
DIR_CONSENSUS <- 1.0

# Jaccard params
TOP_K      <- 50          # kept for compatibility; stratified meta returns all significant features
B          <- 30
SUBSAMPLE  <- 0.8

# For some microbiomeMarker versions, lda_cutoff=0 may error; set to 1e-6 if needed
LEFSE_LDA_CUTOFF_FOR_OUTPUT <- 0

# Fixed benchmark comparison set (direction must match truth)
TARGET_PAIRS <- list(
  c("Disease_mild",   "Disease_severe"),
  c("Healthy",        "Disease_mild"),
  c("Healthy",        "Disease_severe")
)

# ── Small utils ───────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

as_int_rep <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  s <- as.character(x)
  d <- stringr::str_extract(s, "\\d+")
  if (is.na(d)) stop("rep_id 无法解析为整数: ", s)
  as.integer(d)
}

align_groups <- function(groups, sample_names) {
  if (!is.null(names(groups))) {
    missing_samples <- setdiff(sample_names, names(groups))
    if (length(missing_samples) > 0) {
      stop("align_groups: 以下样本在 groups 中找不到对应名称: ",
           paste(missing_samples, collapse = ", "))
    }
    return(groups[sample_names])
  } else {
    if (length(groups) != length(sample_names)) {
      stop("align_groups: groups 长度 (", length(groups),
           ") 与样本数 (", length(sample_names), ") 不一致，且 groups 无名称，无法安全对齐。")
    }
    warning("align_groups: groups 无名称，按位置对齐，请确认顺序与 colnames(abundance_matrix) 完全一致。")
    out <- groups
    names(out) <- sample_names
    out
  }
}

# ── Taxonomy parsing: supports | ; . ──────────────────────────
split_tax_string <- function(x) {
  if (stringr::str_detect(x, "\\|")) return(stringr::str_split(x, "\\|", simplify = TRUE)[1, ])
  if (stringr::str_detect(x, ";"))  return(stringr::str_split(x, ";",  simplify = TRUE)[1, ])
  if (stringr::str_detect(x, "\\."))return(stringr::str_split(x, "\\.", simplify = TRUE)[1, ])
  return(x)
}

parse_taxonomy2 <- function(taxa_strings) {
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_df <- matrix("", nrow = length(taxa_strings), ncol = length(ranks),
                   dimnames = list(taxa_strings, ranks))
  
  for (i in seq_along(taxa_strings)) {
    parts <- split_tax_string(taxa_strings[i])
    parts <- as.character(parts)
    parts <- parts[!is.na(parts) & parts != ""]
    parts <- stringr::str_trim(parts)
    
    clean <- purrr::map_chr(seq_along(parts), function(j) {
      part <- parts[j]
      part <- stringr::str_replace(part, "^[a-z]__", "")
      if (is.na(part) || part == "" || part == "__") part <- paste0("Unknown_", ranks[pmin(j,7)], "_", i)
      part
    })
    
    for (j in seq_len(min(length(clean), length(ranks)))) tax_df[i, j] <- clean[j]
    if (length(clean) < length(ranks)) {
      for (j in (length(clean)+1):length(ranks)) tax_df[i, j] <- paste0("Unknown_", ranks[j], "_", i)
    }
  }
  as.data.frame(tax_df, stringsAsFactors = FALSE)
}

extract_taxa_short_g <- function(taxa_full) {
  g <- stringr::str_extract(taxa_full, "g__[^|;\\.]+")
  g <- stringr::str_trim(g)
  ifelse(is.na(g) | g == "g__", NA_character_, g)
}

# ── Combine p-values (Fisher) ─────────────────────────────────
fisher_combine_p <- function(p_vec) {
  p_vec <- p_vec[!is.na(p_vec)]
  p_vec <- pmax(p_vec, 1e-300)
  if (length(p_vec) == 0) return(NA_real_)
  chi2 <- -2 * sum(log(p_vec))
  df   <- 2 * length(p_vec)
  stats::pchisq(chi2, df = df, lower.tail = FALSE)
}

# ── Standardize marker_table across microbiomeMarker versions ─
# NOTE: keep ALL taxa (no LDA filtering); add lda_pass flag.
standardize_single_batch_markers <- function(markers,
                                             scenario, rep_id,
                                             comparison, group_a, group_b,
                                             batch_label,
                                             lda_thr = LDA_THR) {
  rep_int <- as_int_rep(rep_id)
  
  get_col <- function(df, candidates) {
    nm <- intersect(candidates, names(df))
    if (length(nm) == 0) return(NULL)
    nm[1]
  }
  
  col_feature <- get_col(markers, c("feature","Feature","taxon","Taxon"))
  col_lda     <- get_col(markers, c("ef_lda","lda","lda_score","LDA","ef_lda_mean"))
  col_enrich  <- get_col(markers, c("enrich_group","enriched_group","group","class"))
  col_p       <- get_col(markers, c("pvalue","p_val","pvalue_kw","pvalue_wilcoxon","p.value","p_value"))
  
  out <- markers %>%
    dplyr::mutate(
      feature_raw = if (!is.null(col_feature)) as.character(.data[[col_feature]]) else NA_character_,
      ef_lda_raw   = if (!is.null(col_lda))    suppressWarnings(as.numeric(.data[[col_lda]])) else NA_real_,
      enrich_raw   = if (!is.null(col_enrich)) as.character(.data[[col_enrich]]) else NA_character_,
      p_raw        = if (!is.null(col_p))      suppressWarnings(as.numeric(.data[[col_p]])) else NA_real_
    ) %>%
    dplyr::filter(!is.na(ef_lda_raw)) %>%  # keep rows with LDA
    dplyr::transmute(
      method     = "LEfSe_stratified",
      mode       = "adjusted",
      scenario   = as.character(scenario),
      replicate  = rep_int,
      comparison = as.character(comparison),
      group_a    = as.character(group_a),
      group_b    = as.character(group_b),
      ref_group  = as.character(group_a),
      
      batch      = as.character(batch_label),
      
      taxa_full  = feature_raw,
      taxa_short = extract_taxa_short_g(feature_raw),
      
      effect     = ef_lda_raw,
      pval       = p_raw,
      qval       = NA_real_,  # batch-level BH not used; meta-level BH is the key
      
      lda_pass   = (!is.na(ef_lda_raw) & ef_lda_raw >= lda_thr),
      
      # informational
      significant = (!is.na(p_raw) & p_raw <= 0.05) & lda_pass,
      
      enrich_group = enrich_raw
    )
  
  out
}

# ── Run LEfSe for ONE batch, ONE comparison ───────────────────
run_lefse_one_batch_one_pair <- function(abundance_matrix,
                                         metadata_df,
                                         scenario,
                                         rep_id,
                                         batch_label,
                                         pair,
                                         lda_threshold = LDA_THR,
                                         debug = FALSE) {
  
  abundance_matrix <- ensure_matrix(abundance_matrix)
  
  g1 <- pair[1]; g2 <- pair[2]
  keep_samples <- rownames(metadata_df)[metadata_df$batch == batch_label &
                                          metadata_df$phenotype %in% c(g1, g2)]
  if (length(keep_samples) < 4) {
    if (debug) cat("  [batch", batch_label, "|", g1, "vs", g2, "] 样本不足4，跳过\n")
    return(NULL)
  }
  
  meta <- metadata_df[keep_samples, , drop = FALSE]
  group_counts <- table(meta$phenotype)
  valid_groups <- names(group_counts[group_counts >= 2])
  if (length(valid_groups) < 2) {
    if (debug) cat("  [batch", batch_label, "|", g1, "vs", g2, "] 每组样本不足2，跳过\n")
    return(NULL)
  }
  
  ab <- abundance_matrix[, keep_samples, drop = FALSE]
  ab <- ab[rowSums(ab) > 0, , drop = FALSE]
  if (nrow(ab) == 0) return(NULL)
  
  taxonomy_df <- parse_taxonomy2(rownames(ab))
  meta$phenotype <- factor(meta$phenotype, levels = c(g1, g2))
  
  ps <- tryCatch(
    phyloseq(
      otu_table(ab, taxa_are_rows = TRUE),
      sample_data(meta),
      tax_table(as.matrix(taxonomy_df))
    ),
    error = function(e) {
      if (debug) cat("  [batch", batch_label, "|", g1, "vs", g2, "] phyloseq失败:", e$message, "\n")
      NULL
    }
  )
  if (is.null(ps)) return(NULL)
  
  lefse_res <- tryCatch(
    run_lefse(
      ps,
      wilcoxon_cutoff = 0.05,
      kw_cutoff       = 0.05,
      # key: do not LDA-filter here; we need full tested set for BH downstream
      lda_cutoff      = LEFSE_LDA_CUTOFF_FOR_OUTPUT,
      group           = "phenotype",
      multigrp_strat  = FALSE,
      taxa_rank       = "all"
    ),
    error = function(e) {
      if (debug) cat("  [batch", batch_label, "|", g1, "vs", g2, "] LEfSe失败:", e$message, "\n")
      NULL
    }
  )
  if (is.null(lefse_res)) return(NULL)
  
  markers <- tryCatch({
    mt <- lefse_res@marker_table
    df <- data.frame(mt@.Data, stringsAsFactors = FALSE)
    names(df) <- mt@names
    df
  }, error = function(e) NULL)
  
  if (is.null(markers) || nrow(markers) == 0) return(NULL)
  
  comp_name <- paste0(g1, "_vs_", g2)
  
  markers_std <- standardize_single_batch_markers(
    markers     = markers,
    scenario    = scenario,
    rep_id      = rep_id,
    comparison  = comp_name,
    group_a     = g1,
    group_b     = g2,
    batch_label = batch_label,
    lda_thr     = lda_threshold
  )
  
  attr(markers_std, "n_samples") <- ncol(ab)
  markers_std
}

# ── Meta-analysis across batches for ONE comparison ───────────
# BH is applied to ALL taxa in this comparison-level meta table.
meta_analysis_lefse_pair <- function(batch_markers_list,
                                     qval_thr = QVAL_THR,
                                     lda_thr = LDA_THR,
                                     dir_consensus = DIR_CONSENSUS,
                                     min_batch_hit = NULL) {
  
  batch_markers_list <- Filter(Negate(is.null), batch_markers_list)
  n_batch <- length(batch_markers_list)
  if (n_batch == 0) return(NULL)
  
  if (is.null(min_batch_hit)) min_batch_hit <- ceiling(n_batch / 2)
  
  n_samples_vec <- sapply(batch_markers_list, function(x) attr(x, "n_samples") %||% NA_integer_)
  
  combined <- dplyr::bind_rows(batch_markers_list)
  if (nrow(combined) == 0) return(NULL)
  
  meta_df <- combined %>%
    dplyr::group_by(taxa_full) %>%
    dplyr::summarise(
      n_batch_hit = dplyr::n_distinct(batch),
      n_batch_lda = sum(lda_pass, na.rm = TRUE),
      
      direction_ratio = {
        eg <- enrich_group[lda_pass]
        tbl <- table(eg)
        if (sum(tbl) == 0) NA_real_ else max(tbl) / sum(tbl)
      },
      enrich_group_consensus = {
        eg <- enrich_group[lda_pass]
        tbl <- table(eg)
        if (sum(tbl) == 0) NA_character_ else names(which.max(tbl))
      },
      
      fisher_p = fisher_combine_p(pval),
      
      lda_weighted_mean = {
        batch_labels <- batch
        batch_ids <- match(batch_labels, names(batch_markers_list))
        w_raw <- n_samples_vec[batch_ids]
        w <- sqrt(replace(w_raw, is.na(w_raw), 1))
        stats::weighted.mean(effect, w = w, na.rm = TRUE)
      },
      
      lda_mean = mean(effect, na.rm = TRUE),
      lda_min  = min(effect, na.rm = TRUE),
      lda_max  = max(effect, na.rm = TRUE),
      batches_present = paste(sort(unique(batch)), collapse = ";"),
      
      method     = dplyr::first(method),
      mode       = dplyr::first(mode),
      scenario   = dplyr::first(scenario),
      replicate  = dplyr::first(replicate),
      comparison = dplyr::first(comparison),
      group_a    = dplyr::first(group_a),
      group_b    = dplyr::first(group_b),
      ref_group  = dplyr::first(ref_group),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # key: BH on full tested taxa set
      q_meta = p.adjust(fisher_p, method = "BH"),
      taxa_short = extract_taxa_short_g(taxa_full),
      
      lda_pass_meta = (!is.na(lda_weighted_mean) & lda_weighted_mean >= lda_thr),
      
      # keep your "meta 三重标准"（这里是：n_batch_hit + direction + q_meta）
      # 并额外保留 LDA 阈值（如果你不想要，把 lda_pass_meta 这一行从 significant 中去掉）
      significant =
        (n_batch_hit >= min_batch_hit) &
        lda_pass_meta &
        (!is.na(direction_ratio) & direction_ratio >= dir_consensus) &
        (!is.na(q_meta) & q_meta <= qval_thr)
    ) %>%
    dplyr::arrange(q_meta, desc(lda_weighted_mean))
  
  meta_df
}

# ── Helper: run stratified meta features for Jaccard ──────────
# Returns the set of meta-level significant taxa_full for one pair dataset.
run_stratified_meta_features_one_pair <- function(abundance_matrix,
                                                  metadata_df,
                                                  pair,
                                                  lda_thr = LDA_THR,
                                                  qval_thr = QVAL_THR,
                                                  dir_consensus = DIR_CONSENSUS,
                                                  debug = FALSE) {
  
  abundance_matrix <- ensure_matrix(abundance_matrix)
  
  g1 <- pair[1]; g2 <- pair[2]
  
  keep_samples <- rownames(metadata_df)[metadata_df$phenotype %in% c(g1, g2)]
  if (length(keep_samples) < 4) return(character(0))
  
  ab_all <- abundance_matrix[, keep_samples, drop = FALSE]
  ab_all <- ab_all[rowSums(ab_all) > 0, , drop = FALSE]
  if (nrow(ab_all) == 0) return(character(0))
  
  meta_all <- metadata_df[keep_samples, , drop = FALSE]
  meta_all$phenotype <- factor(meta_all$phenotype, levels = c(g1, g2))
  
  batches <- sort(unique(as.character(meta_all$batch)))
  batch_markers_list <- list()
  
  for (bat in batches) {
    res <- run_lefse_one_batch_one_pair(
      abundance_matrix = ab_all,
      metadata_df      = meta_all,
      scenario         = SCENARIO_B,
      rep_id           = 1,               # placeholder; not used for feature set
      batch_label      = bat,
      pair             = pair,
      lda_threshold    = lda_thr,
      debug            = debug
    )
    if (!is.null(res) && nrow(res) > 0) {
      batch_markers_list[[as.character(bat)]] <- res
    }
  }
  
  meta_res <- meta_analysis_lefse_pair(
    batch_markers_list = batch_markers_list,
    qval_thr      = qval_thr,
    lda_thr       = lda_thr,
    dir_consensus = dir_consensus,
    min_batch_hit = NULL
  )
  
  if (is.null(meta_res) || nrow(meta_res) == 0) return(character(0))
  
  meta_sig <- meta_res %>% dplyr::filter(significant)
  if (nrow(meta_sig) == 0) return(character(0))
  as.character(meta_sig$taxa_full)
}

# ── Main: run stratified LEfSe for scenario B and one rep ─────
run_lefse_stratified_B <- function(rep_id, debug = TRUE) {
  
  scenario_full <- SCENARIO_B
  set.seed(42)
  
  rep_int <- as_int_rep(rep_id)
  
  cat("\n", strrep("=", 60), "\n")
  cat("LEfSe stratified (adjusted) | scenario:", scenario_full, "| Rep:", rep_id, "(int:", rep_int, ")\n")
  cat(strrep("=", 60), "\n")
  flush.console()
  
  t_wall_start <- proc.time()["elapsed"]
  
  inputs <- generate_method_inputs(scenario_full, rep_id)
  abundance_matrix <- ensure_matrix(inputs$lefse$abundance)
  groups <- inputs$lefse$groups
  metadata_source <- inputs$maaslin2$metadata
  
  phenotype_vec <- align_groups(groups, colnames(abundance_matrix))
  
  metadata_df <- tibble(
    sample_id = colnames(abundance_matrix),
    phenotype = as.character(phenotype_vec)
  )
  
  if ("batch" %in% colnames(metadata_source)) {
    batch_map <- metadata_source %>% dplyr::select(sample_id, batch) %>% dplyr::distinct()
    metadata_df <- metadata_df %>% dplyr::left_join(batch_map, by = "sample_id")
  } else {
    warning("metadata 中无 batch 列，所有样本视为 batch1")
    metadata_df$batch <- "batch1"
  }
  
  metadata_df <- as.data.frame(metadata_df)
  rownames(metadata_df) <- metadata_df$sample_id
  
  common_samples <- intersect(colnames(abundance_matrix), rownames(metadata_df))
  abundance_matrix <- abundance_matrix[, common_samples, drop = FALSE]
  metadata_df <- metadata_df[common_samples, , drop = FALSE]
  
  cat("\n总体分布(phenotype x batch):\n")
  print(table(metadata_df$phenotype, metadata_df$batch, dnn = c("phenotype","batch")))
  flush.console()
  
  result_dir <- file.path("benchmark","lefse","lefse_results_genus_stratified",
                          paste0(scenario_full, "_rep", rep_int))
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  batches <- sort(unique(as.character(metadata_df$batch)))
  all_groups <- sort(unique(as.character(metadata_df$phenotype)))
  pairs <- purrr::keep(TARGET_PAIRS, ~all(.x %in% all_groups))
  
  compute_sec <- 0
  batch_raw_all <- list()
  meta_all <- list()
  
  for (pair in pairs) {
    g1 <- pair[1]; g2 <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    
    cat("\n", strrep("-", 50), "\n")
    cat("Comparison:", comp_name, "\n")
    flush.console()
    
    batch_markers_list <- list()
    
    for (bat in batches) {
      cat("  Batch:", bat, "\n"); flush.console()
      
      t0 <- proc.time()["elapsed"]
      res <- run_lefse_one_batch_one_pair(
        abundance_matrix = abundance_matrix,
        metadata_df      = metadata_df,
        scenario         = scenario_full,
        rep_id           = rep_int,
        batch_label      = bat,
        pair             = pair,
        lda_threshold    = LDA_THR,
        debug            = debug
      )
      compute_sec <- compute_sec + as.numeric(proc.time()["elapsed"] - t0)
      
      if (!is.null(res) && nrow(res) > 0) {
        batch_markers_list[[as.character(bat)]] <- res
        batch_raw_all[[length(batch_raw_all) + 1]] <- res
        
        readr::write_csv(res, file.path(result_dir, paste0("batch_", bat, "_", comp_name, "_markers_std.csv")))
        cat("    -> n_taxa_tested:", nrow(res),
            "| n_lda_pass:", sum(res$lda_pass, na.rm = TRUE), "\n")
      } else {
        cat("    -> no taxa\n")
      }
    }
    
    meta_res <- meta_analysis_lefse_pair(
      batch_markers_list = batch_markers_list,
      qval_thr = QVAL_THR,
      lda_thr = LDA_THR,
      dir_consensus = DIR_CONSENSUS,
      min_batch_hit = NULL
    )
    
    if (!is.null(meta_res) && nrow(meta_res) > 0) {
      meta_all[[comp_name]] <- meta_res
      readr::write_csv(meta_res, file.path(result_dir, paste0(comp_name, "_meta_summary_std.csv")))
      cat("  Meta significant:", sum(meta_res$significant, na.rm = TRUE),
          " / total tested:", nrow(meta_res), "\n")
    } else {
      cat("  Meta: no features\n")
    }
  }
  
  batch_raw_df <- dplyr::bind_rows(batch_raw_all)
  meta_df <- dplyr::bind_rows(meta_all)
  
  if (nrow(batch_raw_df) > 0) {
    readr::write_csv(batch_raw_df, file.path(result_dir, "all_batches_all_comparisons_markers_std.csv"))
  }
  if (nrow(meta_df) > 0) {
    readr::write_csv(meta_df, file.path(result_dir, "meta_all_comparisons_summary_std.csv"))
  }
  
  t_wall_end <- proc.time()["elapsed"]
  wall_sec <- as.numeric(t_wall_end - t_wall_start)
  
  cat("\n── Timing ──\n")
  cat("wall_sec   :", round(wall_sec, 2), "\n")
  cat("compute_sec:", round(compute_sec, 2), "\n")
  flush.console()
  
  list(
    batch_raw   = batch_raw_df,
    meta_df     = meta_df,
    wall_sec    = wall_sec,
    compute_sec = compute_sec,
    n_batches   = length(batches),
    n_pairs     = length(pairs)
  )
}

# ── Run all reps (scenario B) ─────────────────────────────────
run_all_stratified_B <- function() {
  
  cat("\n", strrep("=", 60), "\n")
  cat("LEfSe stratified (adjusted) - Scenario B\n")
  cat(strrep("=", 60), "\n\n")
  
  results_all <- list()
  timing_records <- list()
  meta_all_reps <- list()
  batch_all_reps <- list()
  
  timing_df_template <- tibble(
    method      = character(),
    mode        = character(),
    scenario    = character(),
    replicate   = integer(),
    wall_sec    = numeric(),
    compute_sec = numeric(),
    n_batches   = integer(),
    n_pairs     = integer()
  )
  
  for (rep in reps) {
    key <- paste0(SCENARIO_B, "_rep", rep)
    
    out <- tryCatch(
      run_lefse_stratified_B(rep_id = rep, debug = TRUE),
      error = function(e) { cat("Rep", rep, "失败:", e$message, "\n"); NULL }
    )
    
    if (!is.null(out)) {
      results_all[[key]] <- out
      
      timing_records[[length(timing_records) + 1]] <- tibble(
        method      = "LEfSe_stratified",
        mode        = "adjusted",
        scenario    = SCENARIO_B,
        replicate   = as_int_rep(rep),
        wall_sec    = out$wall_sec,
        compute_sec = out$compute_sec,
        n_batches   = out$n_batches,
        n_pairs     = out$n_pairs
      )
      
      if (!is.null(out$meta_df) && nrow(out$meta_df) > 0)  meta_all_reps[[key]]  <- out$meta_df
      if (!is.null(out$batch_raw) && nrow(out$batch_raw) > 0) batch_all_reps[[key]] <- out$batch_raw
    }
  }
  
  timing_df <- if (length(timing_records) > 0) dplyr::bind_rows(timing_records) else timing_df_template
  
  timing_summary <- timing_df %>%
    dplyr::group_by(method, mode, scenario) %>%
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
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(timing_df,      "benchmark/results/lefse_stratified_timing_raw.rds")
  saveRDS(timing_summary, "benchmark/results/lefse_stratified_timing_summary.rds")
  
  meta_all <- dplyr::bind_rows(meta_all_reps)
  batch_all <- dplyr::bind_rows(batch_all_reps)
  
  if (nrow(meta_all) > 0) {
    readr::write_csv(meta_all, "benchmark/results/lefse_stratified_summary_raw.csv")
    saveRDS(meta_all,          "benchmark/results/lefse_stratified_summary_raw.rds")
  } else {
    warning("meta_all 为空：未生成 lefse_stratified_summary_raw")
  }
  
  if (nrow(batch_all) > 0) {
    readr::write_csv(batch_all, "benchmark/results/lefse_stratified_batch_raw.csv")
    saveRDS(batch_all,          "benchmark/results/lefse_stratified_batch_raw.rds")
  }
  
  cat("\n── LEfSe stratified timing summary ──\n")
  print(timing_summary)
  
  list(
    results_all    = results_all,
    timing_raw     = timing_df,
    timing_summary = timing_summary,
    meta_all       = meta_all,
    batch_all      = batch_all
  )
}

# ══════════════════════════════════════════════════════════════
# Jaccard utilities + bootstrap (stratified)
# BH range is automatically correct because it reuses the same meta_analysis logic.
# ══════════════════════════════════════════════════════════════
jaccard <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
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

run_lefse_stratified_jaccard <- function() {
  
  BASE_SEED <- 42
  jaccard_records <- list()
  
  cat("\n", strrep("=", 62), "\n")
  cat("LEfSe stratified Jaccard stability — Scenario B (pairwise)\n")
  cat(strrep("=", 62), "\n")
  
  for (rep in reps) {
    
    cat("\n[Jaccard-Stratified] Rep:", rep, "\n"); flush.console()
    
    inputs <- tryCatch(
      generate_method_inputs(SCENARIO_B, rep),
      error = function(e) {
        cat("数据读取失败:", e$message, "\n"); NULL
      }
    )
    if (is.null(inputs)) next
    
    abundance_matrix <- ensure_matrix(inputs$lefse$abundance)
    groups           <- inputs$lefse$groups
    metadata_source  <- inputs$maaslin2$metadata
    
    phenotype_vec <- align_groups(groups, colnames(abundance_matrix))
    
    metadata_df <- data.frame(
      sample_id        = colnames(abundance_matrix),
      phenotype        = as.character(phenotype_vec),
      stringsAsFactors = FALSE
    )
    rownames(metadata_df) <- metadata_df$sample_id
    
    if ("batch" %in% colnames(metadata_source)) {
      batch_map <- metadata_source %>%
        dplyr::select(sample_id, batch) %>%
        dplyr::distinct()
      metadata_df <- metadata_df %>%
        dplyr::left_join(batch_map, by = "sample_id")
    } else {
      warning("无 batch 列，所有样本视为 batch1")
      metadata_df$batch <- "batch1"
    }
    rownames(metadata_df) <- metadata_df$sample_id
    
    # filter all-0
    col_keep         <- colSums(abundance_matrix) > 0
    abundance_matrix <- abundance_matrix[, col_keep, drop = FALSE]
    row_keep         <- rowSums(abundance_matrix) > 0
    abundance_matrix <- abundance_matrix[row_keep, , drop = FALSE]
    abundance_matrix <- ensure_matrix(abundance_matrix)
    
    # align samples
    common           <- intersect(colnames(abundance_matrix), metadata_df$sample_id)
    abundance_matrix <- abundance_matrix[, common, drop = FALSE]
    metadata_df      <- metadata_df[match(common, metadata_df$sample_id), , drop = FALSE]
    rownames(metadata_df) <- metadata_df$sample_id
    
    if (ncol(abundance_matrix) < 4) {
      cat("  样本过少，跳过\n"); next
    }
    
    cat("  [预处理后] dim:", paste(dim(abundance_matrix), collapse = "x"), "\n")
    cat("  phenotype x batch 分布:\n")
    print(table(metadata_df$phenotype, metadata_df$batch, dnn = c("phenotype", "batch")))
    flush.console()
    
    # pairs: only those existing in data
    all_groups <- sort(unique(as.character(metadata_df$phenotype)))
    disease_pairs <- purrr::keep(TARGET_PAIRS, ~all(.x %in% all_groups))
    n_pairs_total <- length(disease_pairs)
    
    if (n_pairs_total == 0) next
    
    comparison_set_label <- paste(
      vapply(disease_pairs, \(p) paste0(p[1], "_vs_", p[2]), ""),
      collapse = ";"
    )
    
    rep_idx <- as_int_rep(rep)
    
    for (i in seq_along(disease_pairs)) {
      
      pair      <- disease_pairs[[i]]
      g1        <- pair[1]; g2 <- pair[2]
      comp_name <- paste0(g1, "_vs_", g2)
      
      scenario_idx <- 1L
      pair_seed <- BASE_SEED +
        1000000L * scenario_idx +
        10000L   * rep_idx +
        100L     * i
      
      cat("  比较:", comp_name, "\n"); flush.console()
      
      # subset samples for this pair
      samples_all <- rownames(metadata_df)[metadata_df$phenotype %in% pair]
      if (length(samples_all) < 4) {
        cat("    样本数不足，跳过\n"); next
      }
      
      ab_all <- abundance_matrix[, samples_all, drop = FALSE]
      ab_all <- ensure_matrix(ab_all)
      ab_all <- ab_all[rowSums(ab_all) > 0, , drop = FALSE]
      if (nrow(ab_all) == 0 || ncol(ab_all) == 0) next
      
      meta_all <- metadata_df[samples_all, , drop = FALSE]
      meta_all$phenotype <- factor(meta_all$phenotype, levels = c(g1, g2))
      meta_all$batch     <- factor(meta_all$batch)
      
      # S_all (full reference set)
      cat("    计算该pair全量参考集 S_all...\n"); flush.console()
      S_all <- tryCatch(
        run_stratified_meta_features_one_pair(
          abundance_matrix = ab_all,
          metadata_df      = meta_all,
          pair             = pair,
          lda_thr          = LDA_THR,
          qval_thr         = QVAL_THR,
          dir_consensus    = DIR_CONSENSUS,
          debug            = FALSE
        ),
        error = function(e) {
          cat("    S_all 计算失败:", e$message, "\n")
          character(0)
        }
      )
      cat("    S_all 特征数:", length(S_all), "\n"); flush.console()
      
      if (length(S_all) == 0) {
        cat("    全量无显著特征，记录 NA\n")
        jaccard_records[[length(jaccard_records) + 1]] <- tibble::tibble(
          method          = "LEfSe_stratified",
          scenario        = SCENARIO_B,
          rep             = rep_idx,
          comparison      = comp_name,
          comparison_set  = comparison_set_label,
          top_k_per_pair  = as.integer(TOP_K),
          n_pairs         = as.integer(n_pairs_total),
          top_k_total     = as.integer(TOP_K * n_pairs_total),
          jaccard_mean    = NA_real_,
          jaccard_median  = NA_real_,
          jaccard_sd      = NA_real_,
          jaccard_min     = NA_real_,
          jaccard_max     = NA_real_,
          n_boots         = as.integer(B),
          n_boots_valid   = 0L,
          n_features_full = 0L
        )
        next
      }
      
      # bootstrap
      jac_vals <- rep(NA_real_, B)
      
      for (b in seq_len(B)) {
        
        set.seed(pair_seed + b)
        
        sub_ids <- meta_all %>%
          tibble::rownames_to_column("tmp_row") %>%
          dplyr::group_by(batch, phenotype) %>%
          dplyr::slice_sample(prop = SUBSAMPLE) %>%
          dplyr::ungroup() %>%
          dplyr::pull(tmp_row) %>%
          as.character()
        
        if (length(sub_ids) < 4) next
        
        ab_sub <- ab_all[, sub_ids, drop = FALSE]
        ab_sub <- ensure_matrix(ab_sub)
        ab_sub <- ab_sub[rowSums(ab_sub) > 0, , drop = FALSE]
        if (ncol(ab_sub) == 0 || nrow(ab_sub) == 0) next
        
        meta_sub <- meta_all[sub_ids, , drop = FALSE]
        
        batch_ok <- meta_sub %>%
          dplyr::group_by(batch) %>%
          dplyr::summarise(n_grp = dplyr::n_distinct(phenotype), .groups = "drop") %>%
          dplyr::filter(n_grp >= 2) %>%
          nrow()
        if (batch_ok == 0) next
        
        S_b <- tryCatch(
          run_stratified_meta_features_one_pair(
            abundance_matrix = ab_sub,
            metadata_df      = meta_sub,
            pair             = pair,
            lda_thr          = LDA_THR,
            qval_thr         = QVAL_THR,
            dir_consensus    = DIR_CONSENSUS,
            debug            = FALSE
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
        method          = "LEfSe_stratified",
        scenario        = SCENARIO_B,
        rep             = rep_idx,
        comparison      = comp_name,
        comparison_set  = comparison_set_label,
        top_k_per_pair  = as.integer(TOP_K),
        n_pairs         = as.integer(n_pairs_total),
        top_k_total     = as.integer(TOP_K * n_pairs_total),
        jaccard_mean    = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
        jaccard_median  = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
        jaccard_sd      = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
        jaccard_min     = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
        jaccard_max     = if (n_boots_valid == 0) NA_real_ else max(jac_clean),
        n_boots         = as.integer(B),
        n_boots_valid   = as.integer(n_boots_valid),
        n_features_full = as.integer(length(S_all))
      )
      
      cat(sprintf(
        "    [%s | %s | rep%s] Jaccard均值=%.3f (SD=%.3f, 有效boots=%d/%d, S_all=%d)\n",
        SCENARIO_B, comp_name, rep,
        if (n_boots_valid == 0) NA else mean(jac_clean),
        if (n_boots_valid == 0) NA else sd(jac_clean),
        n_boots_valid, B, length(S_all)
      ))
      flush.console()
    }
  }
  
  dir.create("benchmark/results", recursive = TRUE, showWarnings = FALSE)
  
  if (length(jaccard_records) == 0) {
    cat("\n警告：所有 rep 均无有效记录\n")
    empty <- empty_jaccard_tibble()
    saveRDS(empty, "benchmark/results/lefse_stratified_jaccard_raw.rds")
    saveRDS(empty, "benchmark/results/lefse_stratified_jaccard_pair_summary.rds")
    saveRDS(empty, "benchmark/results/lefse_stratified_jaccard_summary.rds")
    return(list(
      jaccard_raw              = empty,
      jaccard_pair_summary     = empty,
      jaccard_scenario_summary = empty
    ))
  }
  
  jaccard_raw <- dplyr::bind_rows(jaccard_records) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
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
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
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
      top_k_total = top_k_per_pair * n_pairs
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
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("jaccard_"),
                                ~ ifelse(is.nan(.), NA_real_, .)))
  
  cat("\n── LEfSe stratified Jaccard summary（Pair级） ──\n")
  print(jaccard_pair_summary)
  cat("\n── LEfSe stratified Jaccard summary（场景级：pair均值） ──\n")
  print(jaccard_scenario_summary)
  
  saveRDS(jaccard_raw,              "benchmark/results/lefse_stratified_jaccard_raw.rds")
  saveRDS(jaccard_pair_summary,     "benchmark/results/lefse_stratified_jaccard_pair_summary.rds")
  saveRDS(jaccard_scenario_summary, "benchmark/results/lefse_stratified_jaccard_summary.rds")
  cat("Jaccard 文件已保存至 benchmark/results/\n")
  
  list(
    jaccard_raw              = jaccard_raw,
    jaccard_pair_summary     = jaccard_pair_summary,
    jaccard_scenario_summary = jaccard_scenario_summary
  )
}

# ── Execute ───────────────────────────────────────────────────
# out_lefse_stratified     <- run_all_stratified_B()
jaccard_lefse_stratified <- run_lefse_stratified_jaccard()



# ══════════════════════════════════════════════════════════════
# LEfSe Adjusted 结果去重代码
# 策略：
# 1. 真重复：完全相同的行 (n_distinct_full == 1) -> 直接去重
# 2. 截断差异：不同分类深度 (n_distinct_full > 1) -> 保留 taxa_full 最长的行
# ══════════════════════════════════════════════════════════════

# ── 1. 加载必要包 ─────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ── 2. 读取数据 ───────────────────────────────────────────────
if (!exists("out_lefse_adjusted")) {
  cat("从文件读取 lefse_adjusted_summary_raw...\n")
  lefse_raw <- readr::read_csv("benchmark/results/lefse_stratified_summary_raw.csv")
} else {
  lefse_raw <- out_lefse_adjusted$final_all
}

cat("原始数据行数:", nrow(lefse_raw), "\n")
cat("数据列名:", paste(names(lefse_raw), collapse = ", "), "\n")

# ── 3. 检查重复情况 ──────────────────────────────────────────
dup_check <- lefse_raw %>%
  dplyr::group_by(replicate, comparison, taxa_short) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup()

cat("\n重复组数（按 replicate + comparison + taxa_short）:",
    length(unique(paste(dup_check$replicate, dup_check$comparison, dup_check$taxa_short))), "\n")
cat("涉及重复的行数:", nrow(dup_check), "\n")

if (nrow(dup_check) > 0) {
  # 分析重复类型
  dup_summary <- dup_check %>%
    dplyr::group_by(replicate, comparison, taxa_short) %>%
    dplyr::summarise(
      n_rows          = dplyr::n(),
      n_distinct_full = dplyr::n_distinct(taxa_full),
      taxa_full_vals  = paste(unique(taxa_full), collapse = " ||| "),
      max_len         = max(nchar(taxa_full)),
      min_len         = min(nchar(taxa_full)),
      .groups = "drop"
    )
  
  cat("\n=== 重复类型分析 ===\n")
  cat("真重复 (taxa_full 完全相同):",
      sum(dup_summary$n_distinct_full == 1), "组\n")
  cat("截断差异 (taxa_full 不同):",
      sum(dup_summary$n_distinct_full > 1), "组\n")
  
  # 展示截断差异的示例
  if (sum(dup_summary$n_distinct_full > 1) > 0) {
    cat("\n截断差异示例（前5组）:\n")
    dup_summary %>%
      dplyr::filter(n_distinct_full > 1) %>%
      dplyr::select(replicate, comparison, taxa_short,
                    n_rows, max_len, min_len, taxa_full_vals) %>%
      head(5) %>%
      print()
  }
  
  # ── 4. 去重策略 ──────────────────────────────────────────
  lefse_dedup <- lefse_raw %>%
    dplyr::mutate(
      taxa_len   = nchar(taxa_full),
      taxa_depth = stringr::str_count(taxa_full, "\\|") + 1
    ) %>%
    dplyr::group_by(replicate, comparison, taxa_short) %>%
    dplyr::arrange(desc(taxa_len), desc(taxa_depth), dplyr::desc(effect)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_len, -taxa_depth)
  
  cat("\n=== 去重结果 ===\n")
  cat("去重前总行数:", nrow(lefse_raw), "\n")
  cat("去重后总行数:", nrow(lefse_dedup), "\n")
  cat("删除重复行数:", nrow(lefse_raw) - nrow(lefse_dedup), "\n")
  
  # 验证去重效果
  dup_after <- lefse_dedup %>%
    dplyr::group_by(replicate, comparison, taxa_short) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  if (nrow(dup_after) == 0) {
    cat("✓ 去重成功：每组 (replicate, comparison, taxa_short) 现在都是唯一的\n")
  } else {
    cat("⚠ 警告：仍有", nrow(dup_after), "行重复，可能需要更细致的处理\n")
  }
  
} else {
  cat("\n✓ 无重复行，无需去重\n")
  lefse_dedup <- lefse_raw
}

# ── 5. 保存去重结果 ──────────────────────────────────────────
output_dir <- "benchmark/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(lefse_dedup,
                 file.path(output_dir, "lefse_stratified_summary_dedup.csv"))
