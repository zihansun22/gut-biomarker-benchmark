
# ============================================================
# 代表性验证 —— 抽样数据准备
# ============================================================

# ── 0. 依赖包 ────────────────────────────────────────────────
library(tidyverse)
library(jsonlite)
library(glue)

# ============================================================
# 1. 用户配置区（按实际情况修改）
# ============================================================
COUNTS_FILE     <- "counts_genus.tsv"
META_FILE       <- "valid_matched.csv"
OUTPUT_DIR      <- "representative_validation"

PROJECT_COL     <- "project_id"   # ← metadata 中 project 列名
SAMPLE_COL      <- "run_id"    # ← metadata 中样本 ID 列名
GROUP_COL       <- "phenotype.x"        # ← metadata 中分组列名（用于两两比较）

SEED            <- 42
SUBSAMPLE_RATIO <- 0.70           # 每个 project 内抽 70%

# 分层抽样数量
N_SMALL         <- 2              # 小档抽 2 个 project
N_MEDIUM        <- 2              # 中档抽 2 个 project
N_LARGE         <- 1              # 大档抽 1 个 project
N_TOP           <- 5              # Top N 最大样本量 project

# ============================================================
# 工具函数
# ============================================================
SEP_LINE <- paste(rep("=", 60), collapse = "")

log_step <- function(step, msg) {
  cat("\n", SEP_LINE, "\n", sep = "")
  cat("STEP ", step, ": ", msg, "\n", sep = "")
  cat(SEP_LINE, "\n", sep = "")
}

log_info <- function(msg) {
  cat("  ", msg, "\n", sep = "")
}

# ============================================================
# STEP 1. 读取原始数据
# ============================================================
log_step(1, "读取原始数据")

counts_raw <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

meta_raw <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)

log_info("counts_genus 维度: {nrow(counts_raw)} features × {ncol(counts_raw)} samples")
log_info("metadata     维度: {nrow(meta_raw)} 行 × {ncol(meta_raw)} 列")
log_info("metadata 列名: {paste(colnames(meta_raw), collapse = ', ')}")

# 列名校验
stopifnot(
  "PROJECT_COL 在 metadata 中不存在，请检查配置" = PROJECT_COL %in% colnames(meta_raw),
  "SAMPLE_COL  在 metadata 中不存在，请检查配置" = SAMPLE_COL  %in% colnames(meta_raw),
  "GROUP_COL   在 metadata 中不存在，请检查配置" = GROUP_COL   %in% colnames(meta_raw)
)

# ============================================================
# STEP 2. 对齐 counts 与 metadata（取交集）
# ============================================================
log_step(2, "对齐 counts 与 metadata")

common_samples <- intersect(colnames(counts_raw), meta_raw[[SAMPLE_COL]])
log_info("counts  样本数: {ncol(counts_raw)}")
log_info("metadata样本数: {nrow(meta_raw)}")
log_info("交集样本数:     {length(common_samples)}")

if (length(common_samples) == 0) {
  stop("交集样本数为 0，请检查 counts 列名与 metadata sample 列是否一致")
}

counts <- counts_raw[, common_samples, drop = FALSE]
meta   <- meta_raw[meta_raw[[SAMPLE_COL]] %in% common_samples, ]

# ============================================================
# STEP 3. 按 project 统计样本量 → 三档划分（tertiles）
# ============================================================
log_step(3, "统计样本量 → 三档划分（tertiles）")

proj_size <- meta |>
  group_by(across(all_of(PROJECT_COL))) |>
  summarise(n_samples = n(), .groups = "drop") |>
  arrange(n_samples)

# 计算 tertile 边界
t33 <- quantile(proj_size$n_samples, 1/3)
t67 <- quantile(proj_size$n_samples, 2/3)
log_info("Tertile 边界 ── 33%: {t33},  67%: {t67}")

proj_size <- proj_size |>
  mutate(
    tier = case_when(
      n_samples <= t33 ~ "small",
      n_samples <= t67 ~ "medium",
      TRUE             ~ "large"
    ),
    tier = factor(tier, levels = c("small", "medium", "large"))
  )

cat("\n  各档 project 数量:\n")
print(table(proj_size$tier))

cat("\n  project 样本量分布:\n")
print(as.data.frame(proj_size), row.names = FALSE)

# ============================================================
# STEP 4. 分层随机抽样 + Top5 合并
# ============================================================
log_step(4, "分层抽样 + Top5 合并")

set.seed(SEED)

safe_sample <- function(pool, n) {
  # 防止池子不足时报错
  n_actual <- min(n, length(pool))
  if (n_actual < n) {
    warning(glue("该档只有 {length(pool)} 个 project，少于目标 {n}，全部选取"))
  }
  sample(pool, size = n_actual, replace = FALSE)
}

proj_small  <- proj_size |> filter(tier == "small")  |> pull(!!sym(PROJECT_COL))
proj_medium <- proj_size |> filter(tier == "medium") |> pull(!!sym(PROJECT_COL))
proj_large  <- proj_size |> filter(tier == "large")  |> pull(!!sym(PROJECT_COL))

drawn_small  <- safe_sample(proj_small,  N_SMALL)
drawn_medium <- safe_sample(proj_medium, N_MEDIUM)
drawn_large  <- safe_sample(proj_large,  N_LARGE)

log_info("[small ] 可选 {length(proj_small )} 个，抽取 {length(drawn_small )}: {paste(drawn_small,  collapse=', ')}")
log_info("[medium] 可选 {length(proj_medium)} 个，抽取 {length(drawn_medium)}: {paste(drawn_medium, collapse=', ')}")
log_info("[large ] 可选 {length(proj_large )} 个，抽取 {length(drawn_large )}: {paste(drawn_large,  collapse=', ')}")

# Top N 最大样本量 project
top_projects <- proj_size |>
  slice_max(order_by = n_samples, n = N_TOP, with_ties = FALSE) |>
  pull(!!sym(PROJECT_COL))
log_info("[Top{N_TOP}] 最大样本 project: {paste(top_projects, collapse=', ')}")

# 合并去重（保留顺序）
all_projects <- unique(c(drawn_small, drawn_medium, drawn_large, top_projects))
log_info("\n✅ 最终选定 {length(all_projects)} 个 project:")

proj_info <- proj_size |>
  filter(!!sym(PROJECT_COL) %in% all_projects) |>
  mutate(in_top5 = !!sym(PROJECT_COL) %in% top_projects)

print(as.data.frame(proj_info), row.names = FALSE)

# ============================================================
# STEP 5. 对每个 project 内随机抽 70% 样本（固定 seed）
# ============================================================
log_step(5, "对每个 project 内随机抽 70% 样本")

sample_draws <- map(seq_along(all_projects), function(i) {
  
  proj <- all_projects[i]
  
  proj_samples <- meta |>
    filter(!!sym(PROJECT_COL) == proj) |>
    pull(!!sym(SAMPLE_COL))
  
  n_total <- length(proj_samples)
  n_pick  <- max(1L, floor(n_total * SUBSAMPLE_RATIO))
  
  # ✅ 方案B：seed = 全局seed + 索引*1000，简单可靠
  proj_seed <- SEED + i * 1000L
  
  set.seed(proj_seed)
  picked <- sample(proj_samples, size = n_pick, replace = FALSE)
  
  log_info("{proj}[{i}]: 总样本 {n_total} → 抽取 {n_pick} 个 (seed={proj_seed})")
  
  list(project = proj, samples = picked, seed = proj_seed,
       n_total = n_total, n_picked = n_pick)
})

names(sample_draws) <- all_projects

# ============================================================
# STEP 6. 构建最终数据集
# ============================================================
log_step(6, "构建最终数据集")

# 汇总所有选中样本
all_selected_samples <- map(sample_draws, "samples") |> unlist() |> unname()

final_meta <- meta |>
  filter(!!sym(SAMPLE_COL) %in% all_selected_samples) |>
  left_join(
    tibble(!!SAMPLE_COL := all_selected_samples,
           source_project = rep(all_projects,
                                times = map_int(sample_draws, ~ length(.x$samples)))),
    by = SAMPLE_COL
  )

final_counts <- counts[, all_selected_samples, drop = FALSE]

log_info("最终 metadata 样本总数: {nrow(final_meta)}")
log_info("最终 counts 维度: {nrow(final_counts)} features × {ncol(final_counts)} samples")

# ============================================================
# STEP 7. 生成 project 摘要与两两比较组合
# ============================================================
log_step(7, "生成 project 摘要 & 两两比较组合")

project_summary <- map_dfr(all_projects, function(proj) {
  
  proj_meta <- final_meta |> filter(!!sym(PROJECT_COL) == proj)
  groups    <- unique(proj_meta[[GROUP_COL]])
  
  # 该 project 内所有两两组合
  if (length(groups) >= 2) {
    combos <- combn(sort(groups), 2, simplify = FALSE) |>
      map_chr(~ paste(.x, collapse = " vs "))
    combo_str <- paste(combos, collapse = " | ")
    n_combos  <- length(combos)
  } else {
    combo_str <- "单组，无法比较"
    n_combos  <- 0
  }
  
  tibble(
    !!PROJECT_COL   := proj,
    tier             = proj_size |> filter(!!sym(PROJECT_COL) == proj) |> pull(tier) |> as.character(),
    n_total          = sample_draws[[proj]]$n_total,
    n_selected       = nrow(proj_meta),
    subsample_pct    = glue("{round(nrow(proj_meta)/sample_draws[[proj]]$n_total*100,1)}%"),
    in_top5          = proj %in% top_projects,
    groups           = paste(sort(groups), collapse = ", "),
    n_pairwise       = n_combos,
    pairwise_combos  = combo_str
  )
})

cat("\n  Project 摘要:\n")
print(as.data.frame(project_summary), row.names = FALSE)

# ============================================================
# STEP 8. 保存所有输出文件
# ============================================================
log_step(8, "保存文件")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 8.1 全局 counts 子集
counts_out <- file.path(OUTPUT_DIR, "counts_genus_repval.tsv")
write.table(final_counts, counts_out, sep = "\t", quote = FALSE, col.names = NA)
log_info("✅ counts 已保存:         {counts_out}")

# 8.2 全局 metadata 子集
meta_out <- file.path(OUTPUT_DIR, "metadata_repval.csv")
write.csv(final_meta, meta_out, row.names = FALSE, quote = TRUE)
log_info("✅ metadata 已保存:       {meta_out}")

# 8.3 project 摘要表
summary_out <- file.path(OUTPUT_DIR, "project_summary.csv")
write.csv(project_summary, summary_out, row.names = FALSE, quote = TRUE)
log_info("✅ project摘要 已保存:    {summary_out}")

# 8.4 每个 project 独立子文件夹
per_proj_dir <- file.path(OUTPUT_DIR, "per_project")
dir.create(per_proj_dir, recursive = TRUE, showWarnings = FALSE)

walk(all_projects, function(proj) {
  
  proj_dir <- file.path(per_proj_dir, proj)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  
  proj_samples <- sample_draws[[proj]]$samples
  proj_meta    <- final_meta |> filter(!!sym(PROJECT_COL) == proj)
  proj_counts  <- final_counts[, proj_samples, drop = FALSE]
  
  write.table(proj_counts,
              file.path(proj_dir, "counts.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  
  write.csv(proj_meta,
            file.path(proj_dir, "metadata.csv"),
            row.names = FALSE, quote = TRUE)
  
  log_info("{proj}: 已保存到 per_project/{proj}/")
})

# 8.5 抽样配置 JSON（复现用）
config <- list(
  seed               = SEED,
  subsample_ratio    = SUBSAMPLE_RATIO,
  tertile_lower      = unname(t33),
  tertile_upper      = unname(t67),
  n_small_drawn      = N_SMALL,
  n_medium_drawn     = N_MEDIUM,
  n_large_drawn      = N_LARGE,
  n_top_added        = N_TOP,
  selected_projects  = all_projects,
  top5_projects      = top_projects,
  total_samples      = nrow(final_meta),
  per_project_seeds  = map(sample_draws, "seed")
)

config_out <- file.path(OUTPUT_DIR, "sampling_config.json")
write_json(config, config_out, pretty = TRUE, auto_unbox = TRUE)
log_info("✅ 抽样配置 已保存:       {config_out}")

# ============================================================
# STEP 9. 最终汇总
# ============================================================

cat("完成！输出文件清单：\n")

all_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in all_files) {
  size_kb <- round(file.size(f) / 1024, 1)
  cat(sprintf("  %-45s %6.1f KB\n", basename(f), size_kb))
}



# ============================================================
# 代表性验证 —— 生成相对丰度文件
# 输入：representative_validation/ 下的 counts 文件
# 输出：与 counts 同目录，文件名 counts → relab
# ============================================================

library(readr)
library(purrr)

OUTPUT_DIR <- "representative_validation"

# ============================================================
# 统一的相对丰度生成函数
# ============================================================
generate_relab_from_counts <- function(counts_path, relab_path) {
  
  # 读取 counts（第一列为 feature/taxa 名）
  counts_df <- read_tsv(counts_path, show_col_types = FALSE)
  
  # 第一列作为行名（feature 列）
  feature_col  <- colnames(counts_df)[1]
  feature_names <- counts_df[[feature_col]]
  counts_matrix <- as.matrix(counts_df[, -1])
  rownames(counts_matrix) <- feature_names
  
  # 列和校验
  col_sums <- colSums(counts_matrix)
  zero_cols <- names(col_sums[col_sums == 0])
  if (length(zero_cols) > 0) {
    warning(sprintf(
      "⚠️  %s 中以下样本列和为0，相对丰度将为NaN: %s",
      basename(counts_path),
      paste(zero_cols, collapse = ", ")
    ))
  }
  
  # 计算相对丰度（列方向归一化）
  relab_matrix <- sweep(counts_matrix, 2, col_sums, "/")
  
  # 转回数据框
  relab_df <- data.frame(
    feature_names,
    relab_matrix,
    check.names = FALSE
  )
  colnames(relab_df)[1] <- feature_col  # 保持原列名
  
  # 保存
  write_tsv(relab_df, relab_path)
  
  cat(sprintf(
    "  ✅ %-55s → %s\n",
    basename(counts_path),
    basename(relab_path)
  ))
}

# ============================================================
# PART A：全局 counts（根目录下）
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("PART A: 全局 counts → relab\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

global_counts_path <- file.path(OUTPUT_DIR, "counts_genus_repval.tsv")
global_relab_path  <- file.path(OUTPUT_DIR, "relab_genus_repval.tsv")

if (file.exists(global_counts_path)) {
  generate_relab_from_counts(global_counts_path, global_relab_path)
} else {
  warning("全局 counts 文件不存在: ", global_counts_path)
}

# ============================================================
# PART B：per_project 子文件夹（遍历每个 project）
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("PART B: per_project/ 下各 project counts → relab\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

per_proj_dir <- file.path(OUTPUT_DIR, "per_project")

# 找到所有 project 子目录下的 counts.tsv
proj_counts_files <- list.files(
  per_proj_dir,
  pattern    = "^counts\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(proj_counts_files) == 0) {
  warning("per_project/ 下未找到任何 counts.tsv 文件，请检查路径")
} else {
  cat(sprintf("  共找到 %d 个 per-project counts 文件\n\n", length(proj_counts_files)))
  
  walk(proj_counts_files, function(counts_path) {
    proj_dir  <- dirname(counts_path)
    relab_path <- file.path(proj_dir, "relab.tsv")   # 同目录下保存 relab.tsv
    generate_relab_from_counts(counts_path, relab_path)
  })
}

# ============================================================
# 最终汇总
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("完成！输出文件清单：\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")

relab_files <- list.files(OUTPUT_DIR, pattern = "relab", recursive = TRUE, full.names = TRUE)

for (f in relab_files) {
  size_kb <- round(file.size(f) / 1024, 1)
  cat(sprintf("  %-55s %7.1f KB\n",
              sub(paste0(OUTPUT_DIR, "/"), "", f),   # 显示相对路径，更清晰
              size_kb))
}

cat(sprintf("\n   共生成 %d 个相对丰度文件\n", length(relab_files)))
cat(sprintf("     输出根目录: %s\n\n", normalizePath(OUTPUT_DIR)))



# ══════════════════════════════════════════════════════════════
# 代表性验证 —— MaAsLin2 Jaccard 稳定性分析
# 修复：max_significance 1.0 确保所有结果写入，再手动筛 q<0.05
# ══════════════════════════════════════════════════════════════

library(Maaslin2)
library(tidyverse)

INPUT_DIR  <- "representative_validation/per_project"
OUTPUT_DIR <- "representative_validation/rep_val_maaslin"

SAMPLE_COL  <- "run_id"
GROUP_COL   <- "phenotype.x"

BASE_SEED    <- 42L
B            <- 30L   # 修改：10 → 30
SUBSAMPLE    <- 0.80
QVAL_THR     <- 0.05
TOP_K        <- 50L

SEP <- paste(rep("=", 60), collapse = "")
log_step <- function(step, msg) {
  cat("\n", SEP, "\nSTEP ", step, ": ", msg, "\n", SEP, "\n", sep = "")
}
log_info <- function(msg) cat("  ", msg, "\n", sep = "")

jaccard <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

ensure_matrix <- function(x) {
  rn <- rownames(x); cn <- colnames(x)
  x  <- as.matrix(x)
  storage.mode(x) <- "numeric"
  rownames(x) <- rn; colnames(x) <- cn
  x
}

empty_jaccard_row <- function(proj, comp_name, n_pairs_total) {
  tibble(
    project         = proj,
    comparison      = comp_name,
    n_pairs         = as.integer(n_pairs_total),
    top_k           = as.integer(TOP_K),
    n_features_full = 0L,
    n_boots         = as.integer(B),
    n_boots_valid   = 0L,
    jaccard_mean    = NA_real_,
    jaccard_median  = NA_real_,
    jaccard_sd      = NA_real_,
    jaccard_min     = NA_real_,
    jaccard_max     = NA_real_
  )
}

# ══════════════════════════════════════════════════════════════
# 核心函数：单次 MaAsLin2
# 修复要点：
#   1. max_significance = 1.0，让所有特征结果都写入 all_results.tsv
#   2. 读取结果后手动用 qval < QVAL_THR 筛选
#   3. 对 qval = NA 的行（数值问题导致）做显式过滤
# ══════════════════════════════════════════════════════════════
run_maaslin2_core <- function(feature_table,
                              meta_df,
                              pair,
                              output_dir,
                              qval_cutoff = QVAL_THR,
                              top_k       = TOP_K) {
  
  common <- intersect(rownames(feature_table), rownames(meta_df))
  if (length(common) < 4) return(character(0))
  
  feature_table <- feature_table[common, , drop = FALSE]
  meta_df       <- meta_df[common,       , drop = FALSE]
  
  pair_samples <- rownames(meta_df)[meta_df$phenotype %in% pair]
  if (length(pair_samples) < 4) return(character(0))
  
  pair_feat <- feature_table[pair_samples, , drop = FALSE]
  pair_meta <- meta_df[pair_samples,       , drop = FALSE]
  
  # 去全 0 列
  pair_feat <- pair_feat[, colSums(pair_feat) > 0, drop = FALSE]
  if (ncol(pair_feat) == 0) return(character(0))
  
  grp_counts <- table(pair_meta$phenotype)
  if (any(grp_counts < 2)) return(character(0))
  ref_group <- if ("Healthy" %in% pair) "Healthy" else sort(pair)[1]
  pair_meta$phenotype <- relevel(
    factor(pair_meta$phenotype, levels = pair), ref = ref_group
  )
  
  use_tmp <- is.null(output_dir)
  if (use_tmp) {
    output_dir <- tempfile(pattern = "maaslin2_tmp_")
    dir.create(output_dir, recursive = TRUE)
  }
  on.exit({
    if (use_tmp) unlink(output_dir, recursive = TRUE)
  }, add = TRUE)
  
  fit <- tryCatch({
    capture.output(
      fit_obj <- Maaslin2(
        input_data       = pair_feat,
        input_metadata   = pair_meta,
        output           = output_dir,
        fixed_effects    = "phenotype",
        normalization    = "TSS",
        transform        = "LOG",
        analysis_method  = "LM",
        correction       = "BH",
        # ★ 修复：改为 1.0，确保所有特征结果都写入 all_results.tsv
        # 原来 0.25 会导致 q>0.25 的特征根本不出现在结果文件里
        max_significance = 1.0,
        min_abundance    = 0.0,
        min_prevalence   = 0.0,
        min_variance     = 0.0,
        standardize      = FALSE,
        plot_heatmap     = FALSE,
        plot_scatter     = FALSE
      ),
      type = "output"
    )
    fit_obj
  }, error = function(e) {
    cat("    [MaAsLin2 失败]", e$message, "\n")
    NULL
  })
  
  if (is.null(fit)) return(character(0))
  
  res_file <- file.path(output_dir, "all_results.tsv")
  if (!file.exists(res_file)) return(character(0))
  
  res_df <- tryCatch(
    read_tsv(res_file, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
  
  # ★ 修复：手动筛 qval < QVAL_THR（原来依赖 max_significance 隐式过滤）
  # 同时过滤 qval = NA（TSS+LOG 归一化后数值问题导致的行）
  feats <- res_df %>%
    filter(!is.na(qval), !is.nan(qval), qval < qval_cutoff) %>%
    arrange(desc(abs(coef))) %>%
    slice_head(n = top_k)
  
  unique(feats$feature)
}

# ══════════════════════════════════════════════════════════════
# 单个 project 的完整 MaAsLin2 + Jaccard 流程
# ══════════════════════════════════════════════════════════════
run_project_jaccard <- function(proj, proj_idx) {
  
  cat("\n", SEP, "\n[Project ", proj_idx, "] ", proj, "\n", SEP, "\n", sep = "")
  
  proj_dir  <- file.path(INPUT_DIR, proj)
  out_proj  <- file.path(OUTPUT_DIR, proj)
  dir.create(out_proj, recursive = TRUE, showWarnings = FALSE)
  
  counts_path <- file.path(proj_dir, "counts.tsv")
  meta_path   <- file.path(proj_dir, "metadata.csv")
  
  if (!file.exists(counts_path) || !file.exists(meta_path)) {
    cat("  ⚠️  缺少 counts.tsv 或 metadata.csv，跳过\n")
    return(NULL)
  }
  
  counts_raw <- read.delim(counts_path, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE)
  meta_raw   <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  log_info(paste0("counts 维度: ", nrow(counts_raw), " features × ",
                  ncol(counts_raw), " samples"))
  log_info(paste0("metadata 行数: ", nrow(meta_raw)))
  
  counts_mat <- ensure_matrix(counts_raw)
  counts_mat <- round(pmax(counts_mat, 0))
  counts_mat <- counts_mat[, colSums(counts_mat) > 0, drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  counts_mat <- counts_mat[rowSums(counts_mat) > 0, , drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  
  # 行=样本，列=taxa
  feature_table <- as.data.frame(t(counts_mat))
  
  meta <- meta_raw %>%
    select(all_of(c(SAMPLE_COL, GROUP_COL))) %>%
    rename(phenotype = all_of(GROUP_COL)) %>%
    mutate(phenotype = as.character(phenotype)) %>%
    filter(!is.na(phenotype), phenotype != "") %>%
    column_to_rownames(SAMPLE_COL)
  
  common <- intersect(rownames(meta), rownames(feature_table))
  cat("  对齐后样本数:", length(common), "\n")
  
  if (length(common) < 4) {
    cat("  ⚠️  对齐后样本数不足 4，跳过\n"); return(NULL)
  }
  meta          <- meta[common, , drop = FALSE]
  feature_table <- feature_table[common, , drop = FALSE]
  
  grp_counts <- table(meta$phenotype)
  valid_grps <- names(grp_counts[grp_counts >= 2])
  if (length(valid_grps) < 2) {
    cat("  ⚠️  有效分组数不足 2，跳过\n"); return(NULL)
  }
  
  meta          <- meta[meta$phenotype %in% valid_grps, , drop = FALSE]
  feature_table <- feature_table[rownames(meta), , drop = FALSE]
  
  cat("  分组分布:\n"); print(table(meta$phenotype))
  
  all_grps      <- sort(unique(meta$phenotype))
  disease_pairs <- combn(all_grps, 2, simplify = FALSE)
  n_pairs_total <- length(disease_pairs)
  cat("  两两比较对数:", n_pairs_total, "\n")
  
  jaccard_records <- list()
  maaslin_results <- list()
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]; g2 <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    cat("\n  [Pair", i, "/", n_pairs_total, "]", comp_name, "\n")
    
    pair_seed <- BASE_SEED + proj_idx * 1000000L + i * 100L
    
    full_out_dir <- file.path(out_proj, "fulldata", comp_name)
    dir.create(full_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    S_all <- tryCatch(
      run_maaslin2_core(
        feature_table = feature_table,
        meta_df       = meta,
        pair          = pair,
        output_dir    = full_out_dir
      ),
      error = function(e) {
        cat("    全量运行失败:", e$message, "\n"); character(0)
      }
    )
    cat("    全量显著特征数:", length(S_all), "\n")
    
    res_file <- file.path(full_out_dir, "all_results.tsv")
    if (file.exists(res_file)) {
      res_df <- read_tsv(res_file, show_col_types = FALSE) %>%
        mutate(project = proj, comparison = comp_name,
               group_a = g1, group_b = g2)
      maaslin_results[[comp_name]] <- res_df
      write_csv(res_df, file.path(out_proj, paste0(comp_name, "_fulldata_results.csv")))
    }
    
    if (length(S_all) == 0) {
      jaccard_records[[length(jaccard_records) + 1]] <-
        empty_jaccard_row(proj, comp_name, n_pairs_total)
      next
    }
    
    jac_vals <- rep(NA_real_, B)
    
    for (b in seq_len(B)) {
      
      set.seed(pair_seed + b)
      
      pair_samples_all <- rownames(meta)[meta$phenotype %in% pair]
      meta_pair        <- meta[pair_samples_all, , drop = FALSE]
      
      sub_ids <- meta_pair %>%
        rownames_to_column("sid") %>%
        group_by(phenotype) %>%
        slice_sample(prop = SUBSAMPLE) %>%
        ungroup() %>%
        pull(sid) %>%
        as.character()
      
      if (length(sub_ids) < 4) next
      
      feat_sub <- feature_table[sub_ids, , drop = FALSE]
      feat_sub <- feat_sub[, colSums(feat_sub) > 0, drop = FALSE]
      if (ncol(feat_sub) == 0) next
      
      meta_sub <- meta[sub_ids, , drop = FALSE]
      if (n_distinct(meta_sub$phenotype) < 2) next
      
      S_b <- tryCatch(
        run_maaslin2_core(
          feature_table = feat_sub,
          meta_df       = meta_sub,
          pair          = pair,
          output_dir    = NULL
        ),
        error = function(e) character(0)
      )
      
      jac_vals[b] <- jaccard(S_b, S_all)
    }
    
    jac_clean     <- jac_vals[!is.na(jac_vals)]
    n_boots_valid <- length(jac_clean)
    cat("    有效 bootstrap 次数:", n_boots_valid, "/", B, "\n")
    if (n_boots_valid > 0)
      cat(sprintf("    Jaccard 均值: %.3f  中位数: %.3f\n",
                  mean(jac_clean), median(jac_clean)))
    
    jaccard_records[[length(jaccard_records) + 1]] <- tibble(
      project         = proj,
      comparison      = comp_name,
      n_pairs         = as.integer(n_pairs_total),
      top_k           = as.integer(TOP_K),
      n_features_full = as.integer(length(S_all)),
      n_boots         = as.integer(B),
      n_boots_valid   = as.integer(n_boots_valid),
      jaccard_mean    = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
      jaccard_median  = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
      jaccard_sd      = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
      jaccard_min     = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
      jaccard_max     = if (n_boots_valid == 0) NA_real_ else max(jac_clean)
    )
  }
  
  proj_jaccard <- bind_rows(jaccard_records)
  write_csv(proj_jaccard, file.path(out_proj, "jaccard_results.csv"))
  
  if (length(maaslin_results) > 0) {
    write_csv(bind_rows(maaslin_results),
              file.path(out_proj, "all_comparisons_results.csv"))
  }
  
  log_info(paste0("✅ ", proj, " 完成，结果保存至 ", out_proj))
  proj_jaccard
}

# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════
log_step("MAIN", "代表性验证 MaAsLin2 Jaccard 分析")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

all_projects <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
all_projects <- all_projects[all_projects != ""]
cat("检测到", length(all_projects), "个 project:",
    paste(all_projects, collapse = ", "), "\n")

all_jaccard  <- list()
failed_projs <- character(0)

for (idx in seq_along(all_projects)) {
  proj <- all_projects[idx]
  result <- tryCatch(
    run_project_jaccard(proj, idx),
    error = function(e) {
      cat("  ❌ Project", proj, "失败:", e$message, "\n")
      failed_projs <<- c(failed_projs, proj)
      NULL
    }
  )
  if (!is.null(result)) all_jaccard[[proj]] <- result
}

# ══════════════════════════════════════════════════════════════
# 全局汇总
# ══════════════════════════════════════════════════════════════
log_step("SUMMARY", "全局 Jaccard 汇总")

if (length(all_jaccard) > 0) {
  
  jaccard_raw <- bind_rows(all_jaccard) %>%
    mutate(across(starts_with("jaccard_"), ~ ifelse(is.nan(.), NA_real_, .)))
  
  write_csv(jaccard_raw, file.path(OUTPUT_DIR, "jaccard_raw.csv"))
  saveRDS(jaccard_raw,   file.path(OUTPUT_DIR, "jaccard_raw.rds"))
  
  cat("\n── Raw Jaccard 明细（前10行）──\n")
  print(head(jaccard_raw, 10))
  
  jaccard_proj_summary <- jaccard_raw %>%
    group_by(project) %>%
    summarise(
      n_pairs_run          = sum(!is.na(jaccard_mean)),
      n_pairs_total        = max(n_pairs),
      n_features_full_mean = mean(n_features_full, na.rm = TRUE),
      n_boots_valid_mean   = mean(n_boots_valid,   na.rm = TRUE),
      jaccard_mean         = mean(jaccard_mean,    na.rm = TRUE),
      jaccard_median       = mean(jaccard_median,  na.rm = TRUE),
      jaccard_sd_mean      = mean(jaccard_sd,      na.rm = TRUE),
      jaccard_min          = min(jaccard_min,       na.rm = TRUE),
      jaccard_max          = max(jaccard_max,       na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(jaccard_proj_summary, file.path(OUTPUT_DIR, "jaccard_project_summary.csv"))
  saveRDS(jaccard_proj_summary,   file.path(OUTPUT_DIR, "jaccard_project_summary.rds"))
  
  cat("\n── Project 级 Jaccard 汇总 ──\n")
  print(jaccard_proj_summary)
  
  jaccard_comp_summary <- jaccard_raw %>%
    group_by(comparison) %>%
    summarise(
      n_projects         = n(),
      jaccard_mean       = mean(jaccard_mean,  na.rm = TRUE),
      jaccard_median     = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_between = sd(jaccard_mean,    na.rm = TRUE),
      jaccard_min        = min(jaccard_min,     na.rm = TRUE),
      jaccard_max        = max(jaccard_max,     na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(jaccard_mean))
  
  write_csv(jaccard_comp_summary, file.path(OUTPUT_DIR, "jaccard_comparison_summary.csv"))
  saveRDS(jaccard_comp_summary,   file.path(OUTPUT_DIR, "jaccard_comparison_summary.rds"))
  
  cat("\n── Comparison 级 Jaccard 汇总 ──\n")
  print(jaccard_comp_summary)
  
} else {
  cat("⚠️  所有 project 均未产生有效 Jaccard 结果\n")
}

cat("\n", SEP, "\n完成！输出文件清单：\n", SEP, "\n", sep = "")
all_out_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in all_out_files) {
  size_kb  <- round(file.size(f) / 1024, 1)
  rel_path <- sub(paste0(normalizePath(OUTPUT_DIR, winslash = "/"), "/"), "",
                  normalizePath(f, winslash = "/"))
  cat(sprintf("  %-60s %7.1f KB\n", rel_path, size_kb))
}
cat(sprintf("\n  输出根目录:  %s\n",  normalizePath(OUTPUT_DIR)))
cat(sprintf("  成功 project: %d / %d\n", length(all_jaccard), length(all_projects)))
if (length(failed_projs) > 0)
  cat(sprintf("  失败 project: %s\n", paste(failed_projs, collapse = ", ")))



# ══════════════════════════════════════════════════════════════
# 代表性验证 —— ANCOM-BC2 Jaccard 稳定性分析
# 输入：representative_validation/per_project/ 下各 project
# 输出：representative_validation/rep_val_ancombc2/
# ══════════════════════════════════════════════════════════════

setwd("D:/大四/毕设/data")

library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# ══════════════════════════════════════════════════════════════
# 全局参数
# ══════════════════════════════════════════════════════════════
INPUT_DIR  <- "representative_validation/per_project"
OUTPUT_DIR <- "representative_validation/rep_val_ancombc2"

SAMPLE_COL <- "run_id"
GROUP_COL  <- "phenotype.x"

BASE_SEED  <- 42L
B          <- 10L
SUBSAMPLE  <- 0.80
QVAL_THR   <- 0.05
TOP_K      <- 50L

# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════
SEP <- paste(rep("=", 60), collapse = "")

log_step <- function(step, msg) {
  cat("\n", SEP, "\nSTEP ", step, ": ", msg, "\n", SEP, "\n", sep = "")
}
log_info <- function(msg) cat("  ", msg, "\n", sep = "")

jaccard <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# 修复：保留行列名
ensure_matrix <- function(x) {
  rn <- rownames(x)
  cn <- colnames(x)
  x  <- as.matrix(x)
  storage.mode(x) <- "numeric"
  rownames(x) <- rn
  colnames(x) <- cn
  x
}

empty_jaccard_row <- function(proj, comp_name, n_pairs_total) {
  tibble(
    project         = proj,
    comparison      = comp_name,
    n_pairs         = as.integer(n_pairs_total),
    top_k           = as.integer(TOP_K),
    n_features_full = 0L,
    n_boots         = as.integer(B),
    n_boots_valid   = 0L,
    jaccard_mean    = NA_real_,
    jaccard_median  = NA_real_,
    jaccard_sd      = NA_real_,
    jaccard_min     = NA_real_,
    jaccard_max     = NA_real_
  )
}

# ══════════════════════════════════════════════════════════════
# 核心函数：单次 ANCOM-BC2，返回显著特征名向量
# counts_matrix : taxa × samples（numeric，整数计数）
# metadata_df   : data.frame，含 sample_id 列和 phenotype 列
# pair_levels   : c("G1", "G2")
# ══════════════════════════════════════════════════════════════
run_ancombc2_core <- function(counts_matrix,
                              metadata_df,
                              pair_levels,
                              output_dir  = NULL,
                              qval_cutoff = QVAL_THR,
                              top_k       = TOP_K) {
  
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # 规范化 metadata：确保有 sample_id 列
  if (!"sample_id" %in% colnames(metadata_df)) {
    metadata_df <- tibble::rownames_to_column(as.data.frame(metadata_df), "sample_id")
  }
  metadata_df <- metadata_df %>%
    select(sample_id, phenotype) %>%
    mutate(sample_id = as.character(sample_id),
           phenotype  = as.character(phenotype))
  
  # 对齐样本
  common <- intersect(colnames(counts_matrix), metadata_df$sample_id)
  if (length(common) < 4) return(character(0))
  
  counts_matrix <- counts_matrix[, common, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  metadata_df   <- metadata_df[match(common, metadata_df$sample_id), , drop = FALSE]
  
  # 过滤 pair 内样本
  metadata_df <- metadata_df[metadata_df$phenotype %in% pair_levels, , drop = FALSE]
  if (nrow(metadata_df) < 4) return(character(0))
  if (n_distinct(metadata_df$phenotype) < 2) return(character(0))
  
  counts_matrix <- counts_matrix[, metadata_df$sample_id, drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  
  # 去全 0 特征
  counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, , drop = FALSE]
  counts_matrix <- ensure_matrix(counts_matrix)
  if (nrow(counts_matrix) == 0) return(character(0))
  
  # 各组样本数校验
  grp_counts <- table(metadata_df$phenotype)
  if (any(grp_counts < 2)) return(character(0))
  
  # 参考组
  ref_group <- if ("Healthy" %in% pair_levels) "Healthy" else sort(pair_levels)[1]
  metadata_df$phenotype <- relevel(
    factor(metadata_df$phenotype, levels = pair_levels), ref = ref_group
  )
  
  # 构建 phyloseq
  # 修复：重置行索引，避免残留数字 rownames 导致 column_to_rownames 报错
  sample_data_df <- metadata_df %>%
    as.data.frame() %>%
    `rownames<-`(NULL) %>%
    column_to_rownames("sample_id")
  
  physeq <- tryCatch(
    phyloseq(
      otu_table(counts_matrix, taxa_are_rows = TRUE),
      sample_data(sample_data_df)
    ),
    error = function(e) {
      cat("    [phyloseq 构建失败]", e$message, "\n"); NULL
    }
  )
  if (is.null(physeq)) return(character(0))
  
  # 运行 ANCOM-BC2
  fit <- tryCatch(
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
      cat("    [ANCOM-BC2 失败]", e$message, "\n"); NULL
    }
  )
  if (is.null(fit)) return(character(0))
  
  res_df <- fit$res
  if (is.null(res_df) || nrow(res_df) == 0) return(character(0))
  
  # 定位结果列
  nonref    <- setdiff(pair_levels, ref_group)
  if (length(nonref) != 1) return(character(0))
  pheno_col <- paste0("phenotype", nonref)
  q_col     <- paste0("q_",   pheno_col)
  lfc_col   <- paste0("lfc_", pheno_col)
  
  if (!all(c("taxon", q_col, lfc_col) %in% colnames(res_df))) return(character(0))
  
  feats <- res_df %>%
    filter(!is.na(.data[[q_col]]), .data[[q_col]] < qval_cutoff) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    slice_head(n = top_k)
  
  # 落盘
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    result_full <- res_df %>%
      mutate(ref_group = ref_group,
             pair      = paste(pair_levels, collapse = "_vs_"))
    write_csv(result_full, file.path(output_dir, "results.csv"))
  }
  
  unique(as.character(feats$taxon))
}

# ══════════════════════════════════════════════════════════════
# 单个 project 的完整 ANCOM-BC2 + Jaccard 流程
# ══════════════════════════════════════════════════════════════
run_project_jaccard <- function(proj, proj_idx) {
  
  cat("\n", SEP, "\n[Project ", proj_idx, "] ", proj, "\n", SEP, "\n", sep = "")
  
  proj_dir <- file.path(INPUT_DIR,  proj)
  out_proj <- file.path(OUTPUT_DIR, proj)
  dir.create(out_proj, recursive = TRUE, showWarnings = FALSE)
  
  # ── 读取数据 ───────────────────────────────────────────────
  counts_path <- file.path(proj_dir, "counts.tsv")
  meta_path   <- file.path(proj_dir, "metadata.csv")
  
  if (!file.exists(counts_path) || !file.exists(meta_path)) {
    cat("  ⚠️  缺少 counts.tsv 或 metadata.csv，跳过\n")
    return(NULL)
  }
  
  counts_raw <- read.delim(counts_path, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE)
  meta_raw   <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  log_info(paste0("counts 维度: ", nrow(counts_raw),
                  " features × ", ncol(counts_raw), " samples"))
  log_info(paste0("metadata 行数: ", nrow(meta_raw)))
  
  # ── 诊断：确认列名与 sample_id 能对齐 ─────────────────────
  cat("  counts 列名示例:      ", paste(head(colnames(counts_raw), 3), collapse = ", "), "\n")
  cat("  metadata sample_id 示例:", paste(head(meta_raw[[SAMPLE_COL]], 3), collapse = ", "), "\n")
  
  # ── 清洗 counts（taxa × samples，保留行列名）──────────────
  counts_mat <- ensure_matrix(counts_raw)
  counts_mat <- round(pmax(counts_mat, 0))
  # 去全 0 列，显式保留列名
  keep_cols  <- colSums(counts_mat) > 0
  counts_mat <- counts_mat[, keep_cols, drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  # 去全 0 行，显式保留行名
  keep_rows  <- rowSums(counts_mat) > 0
  counts_mat <- counts_mat[keep_rows, , drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  
  # ── 整理 metadata ──────────────────────────────────────────
  meta <- meta_raw %>%
    select(all_of(c(SAMPLE_COL, GROUP_COL))) %>%
    rename(sample_id = all_of(SAMPLE_COL),
           phenotype  = all_of(GROUP_COL)) %>%
    mutate(sample_id = as.character(sample_id),
           phenotype  = as.character(phenotype)) %>%
    filter(!is.na(phenotype), phenotype != "")
  
  # ── 对齐样本 ───────────────────────────────────────────────
  common <- intersect(colnames(counts_mat), meta$sample_id)
  cat("  对齐后样本数:", length(common), "\n")
  
  if (length(common) < 4) {
    cat("  ⚠️  对齐后样本数不足 4，跳过\n")
    return(NULL)
  }
  
  counts_mat <- counts_mat[, common, drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  meta       <- meta[match(common, meta$sample_id), , drop = FALSE]
  
  # ── 过滤有效分组（每组 ≥ 2）──────────────────────────────
  grp_counts <- table(meta$phenotype)
  valid_grps <- names(grp_counts[grp_counts >= 2])
  if (length(valid_grps) < 2) {
    cat("  ⚠️  有效分组数不足 2，跳过\n")
    return(NULL)
  }
  
  meta       <- meta[meta$phenotype %in% valid_grps, , drop = FALSE]
  counts_mat <- counts_mat[, meta$sample_id, drop = FALSE]
  counts_mat <- ensure_matrix(counts_mat)
  
  cat("  分组分布:\n"); print(table(meta$phenotype))
  
  # ── 枚举两两比较对 ─────────────────────────────────────────
  all_grps      <- sort(unique(meta$phenotype))
  disease_pairs <- combn(all_grps, 2, simplify = FALSE)
  n_pairs_total <- length(disease_pairs)
  cat("  两两比较对数:", n_pairs_total, "\n")
  
  jaccard_records <- list()
  maaslin_results <- list()
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]; g2 <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    cat("\n  [Pair", i, "/", n_pairs_total, "]", comp_name, "\n")
    
    pair_seed <- BASE_SEED + proj_idx * 1000000L + i * 100L
    
    # pair 子集
    pair_samples <- meta$sample_id[meta$phenotype %in% pair]
    if (length(pair_samples) < 4) {
      cat("    样本数不足，跳过\n")
      jaccard_records[[length(jaccard_records) + 1]] <-
        empty_jaccard_row(proj, comp_name, n_pairs_total)
      next
    }
    
    counts_pair <- counts_mat[, pair_samples, drop = FALSE]
    counts_pair <- ensure_matrix(counts_pair)
    counts_pair <- counts_pair[rowSums(counts_pair) > 0, , drop = FALSE]
    counts_pair <- ensure_matrix(counts_pair)
    meta_pair   <- meta[meta$sample_id %in% pair_samples, , drop = FALSE]
    
    if (nrow(counts_pair) == 0) {
      jaccard_records[[length(jaccard_records) + 1]] <-
        empty_jaccard_row(proj, comp_name, n_pairs_total)
      next
    }
    
    # ── Step 1：全量运行（落盘供复查）─────────────────────
    full_out_dir <- file.path(out_proj, "fulldata", comp_name)
    
    S_all <- tryCatch(
      run_ancombc2_core(
        counts_matrix = counts_pair,
        metadata_df   = meta_pair,
        pair_levels   = pair,
        output_dir    = full_out_dir
      ),
      error = function(e) {
        cat("    全量运行失败:", e$message, "\n"); character(0)
      }
    )
    cat("    全量显著特征数:", length(S_all), "\n")
    
    res_file <- file.path(full_out_dir, "results.csv")
    if (file.exists(res_file)) {
      res_df <- read_csv(res_file, show_col_types = FALSE) %>%
        mutate(project = proj, comparison = comp_name,
               group_a = g1,  group_b    = g2)
      maaslin_results[[comp_name]] <- res_df
    }
    
    if (length(S_all) == 0) {
      jaccard_records[[length(jaccard_records) + 1]] <-
        empty_jaccard_row(proj, comp_name, n_pairs_total)
      next
    }
    
    # ── Step 2：B 次 bootstrap Jaccard ──────────────────
    jac_vals <- rep(NA_real_, B)
    
    for (b in seq_len(B)) {
      
      set.seed(pair_seed + b)
      
      sub_ids <- meta_pair %>%
        group_by(phenotype) %>%
        slice_sample(prop = SUBSAMPLE) %>%
        ungroup() %>%
        pull(sample_id) %>%
        as.character()
      
      if (length(sub_ids) < 4) next
      
      counts_sub <- counts_pair[, sub_ids, drop = FALSE]
      counts_sub <- ensure_matrix(counts_sub)
      counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
      counts_sub <- ensure_matrix(counts_sub)
      meta_sub   <- meta_pair[meta_pair$sample_id %in% sub_ids, , drop = FALSE]
      
      if (nrow(counts_sub) == 0 || ncol(counts_sub) == 0) next
      if (n_distinct(meta_sub$phenotype) < 2) next
      
      S_b <- tryCatch(
        run_ancombc2_core(
          counts_matrix = counts_sub,
          metadata_df   = meta_sub,
          pair_levels   = pair,
          output_dir    = NULL
        ),
        error = function(e) character(0)
      )
      
      jac_vals[b] <- jaccard(S_b, S_all)
    }
    
    jac_clean     <- jac_vals[!is.na(jac_vals)]
    n_boots_valid <- length(jac_clean)
    cat("    有效 bootstrap 次数:", n_boots_valid, "/", B, "\n")
    if (n_boots_valid > 0)
      cat(sprintf("    Jaccard 均值: %.3f  中位数: %.3f\n",
                  mean(jac_clean), median(jac_clean)))
    
    jaccard_records[[length(jaccard_records) + 1]] <- tibble(
      project         = proj,
      comparison      = comp_name,
      n_pairs         = as.integer(n_pairs_total),
      top_k           = as.integer(TOP_K),
      n_features_full = as.integer(length(S_all)),
      n_boots         = as.integer(B),
      n_boots_valid   = as.integer(n_boots_valid),
      jaccard_mean    = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
      jaccard_median  = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
      jaccard_sd      = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
      jaccard_min     = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
      jaccard_max     = if (n_boots_valid == 0) NA_real_ else max(jac_clean)
    )
  }
  
  # ── 保存当前 project 结果 ─────────────────────────────────
  proj_jaccard <- bind_rows(jaccard_records)
  write_csv(proj_jaccard, file.path(out_proj, "jaccard_results.csv"))
  
  if (length(maaslin_results) > 0) {
    write_csv(bind_rows(maaslin_results),
              file.path(out_proj, "all_comparisons_results.csv"))
  }
  
  log_info(paste0("✅ ", proj, " 完成，结果保存至 ", out_proj))
  proj_jaccard
}

# ══════════════════════════════════════════════════════════════
# 主流程：遍历所有 project
# ══════════════════════════════════════════════════════════════
log_step("MAIN", "代表性验证 ANCOM-BC2 Jaccard 分析")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

all_projects <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
all_projects <- all_projects[all_projects != ""]
cat("检测到", length(all_projects), "个 project:",
    paste(all_projects, collapse = ", "), "\n")

all_jaccard  <- list()
failed_projs <- character(0)

for (idx in seq_along(all_projects)) {
  proj <- all_projects[idx]
  result <- tryCatch(
    run_project_jaccard(proj, idx),
    error = function(e) {
      cat("  ❌ Project", proj, "失败:", e$message, "\n")
      failed_projs <<- c(failed_projs, proj)
      NULL
    }
  )
  if (!is.null(result)) all_jaccard[[proj]] <- result
}

# ══════════════════════════════════════════════════════════════
# 全局汇总
# ══════════════════════════════════════════════════════════════
log_step("SUMMARY", "全局 Jaccard 汇总")

if (length(all_jaccard) > 0) {
  
  jaccard_raw <- bind_rows(all_jaccard) %>%
    mutate(across(starts_with("jaccard_"),
                  ~ ifelse(is.nan(.), NA_real_, .)))
  
  write_csv(jaccard_raw, file.path(OUTPUT_DIR, "jaccard_raw.csv"))
  saveRDS(jaccard_raw,   file.path(OUTPUT_DIR, "jaccard_raw.rds"))
  
  cat("\n── Raw Jaccard 明细（前10行）──\n")
  print(head(jaccard_raw, 10))
  
  jaccard_proj_summary <- jaccard_raw %>%
    group_by(project) %>%
    summarise(
      n_pairs_run          = sum(!is.na(jaccard_mean)),
      n_pairs_total        = max(n_pairs),
      n_features_full_mean = mean(n_features_full, na.rm = TRUE),
      n_boots_valid_mean   = mean(n_boots_valid,   na.rm = TRUE),
      jaccard_mean         = mean(jaccard_mean,    na.rm = TRUE),
      jaccard_median       = mean(jaccard_median,  na.rm = TRUE),
      jaccard_sd_mean      = mean(jaccard_sd,      na.rm = TRUE),
      jaccard_min          = min(jaccard_min,       na.rm = TRUE),
      jaccard_max          = max(jaccard_max,       na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(jaccard_proj_summary, file.path(OUTPUT_DIR, "jaccard_project_summary.csv"))
  saveRDS(jaccard_proj_summary,   file.path(OUTPUT_DIR, "jaccard_project_summary.rds"))
  
  cat("\n── Project 级 Jaccard 汇总 ──\n")
  print(jaccard_proj_summary)
  
  jaccard_comp_summary <- jaccard_raw %>%
    group_by(comparison) %>%
    summarise(
      n_projects         = n(),
      jaccard_mean       = mean(jaccard_mean,  na.rm = TRUE),
      jaccard_median     = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_between = sd(jaccard_mean,    na.rm = TRUE),
      jaccard_min        = min(jaccard_min,     na.rm = TRUE),
      jaccard_max        = max(jaccard_max,     na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(jaccard_mean))
  
  write_csv(jaccard_comp_summary, file.path(OUTPUT_DIR, "jaccard_comparison_summary.csv"))
  saveRDS(jaccard_comp_summary,   file.path(OUTPUT_DIR, "jaccard_comparison_summary.rds"))
  
  cat("\n── Comparison 级 Jaccard 汇总 ──\n")
  print(jaccard_comp_summary)
  
} else {
  cat("⚠️  所有 project 均未产生有效 Jaccard 结果\n")
}

# ══════════════════════════════════════════════════════════════
# 最终文件清单
# ══════════════════════════════════════════════════════════════
cat("\n", SEP, "\n完成！输出文件清单：\n", SEP, "\n", sep = "")

all_out_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in all_out_files) {
  size_kb  <- round(file.size(f) / 1024, 1)
  rel_path <- sub(paste0(normalizePath(OUTPUT_DIR), .Platform$file.sep), "", normalizePath(f))
  cat(sprintf("  %-60s %7.1f KB\n", rel_path, size_kb))
}

cat(sprintf("\n  输出根目录:   %s\n",  normalizePath(OUTPUT_DIR)))
cat(sprintf("  成功 project: %d / %d\n", length(all_jaccard), length(all_projects)))
if (length(failed_projs) > 0)
  cat(sprintf("  失败 project: %s\n", paste(failed_projs, collapse = ", ")))




# ══════════════════════════════════════════════════════════════
# 代表性验证 —— LEfSe Jaccard 稳定性分析
# 输入：representative_validation/per_project/ 下各 project
# 输出：representative_validation/rep_val_lefse/
# 数据读取逻辑与 MaAsLin2 / ANCOMBC 版完全一致
# ══════════════════════════════════════════════════════════════

setwd("D:/大四/毕设/data")

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(microbiomeMarker)
})

# ══════════════════════════════════════════════════════════════
# 全局参数
# ══════════════════════════════════════════════════════════════
INPUT_DIR  <- "representative_validation/per_project"
OUTPUT_DIR <- "representative_validation/rep_val_lefse"

SAMPLE_COL <- "run_id"
GROUP_COL  <- "phenotype.x"

BASE_SEED  <- 42L
B          <- 10L
SUBSAMPLE  <- 0.80
QVAL_THR   <- 0.05
TOP_K      <- 50L

# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════
SEP <- paste(rep("=", 60), collapse = "")

log_step <- function(step, msg) {
  cat("\n", SEP, "\nSTEP ", step, ": ", msg, "\n", SEP, "\n", sep = "")
}
log_info <- function(msg) cat("  ", msg, "\n", sep = "")

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

empty_jaccard_row <- function(proj, comp_name, n_pairs_total) {
  tibble(
    project         = proj,
    comparison      = comp_name,
    n_pairs         = as.integer(n_pairs_total),
    top_k           = as.integer(TOP_K),
    n_features_full = 0L,
    m_tested        = 0L,
    n_boots         = as.integer(B),
    n_boots_valid   = 0L,
    jaccard_mean    = NA_real_,
    jaccard_median  = NA_real_,
    jaccard_sd      = NA_real_,
    jaccard_min     = NA_real_,
    jaccard_max     = NA_real_
  )
}

# ══════════════════════════════════════════════════════════════
# Taxonomy 解析（与原 LEfSe 代码完全一致）
# ══════════════════════════════════════════════════════════════
split_tax_string <- function(x) {
  if (stringr::str_detect(x, ";"))   return(stringr::str_split(x, ";",   simplify = TRUE)[1, ])
  if (stringr::str_detect(x, "\\|")) return(stringr::str_split(x, "\\|", simplify = TRUE)[1, ])
  if (stringr::str_detect(x, "\\.")) return(stringr::str_split(x, "\\.", simplify = TRUE)[1, ])
  return(x)
}

parse_taxonomy2 <- function(taxa_strings) {
  ranks  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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

# ══════════════════════════════════════════════════════════════
# Taxa 映射（与原 LEfSe 代码完全一致）
# ══════════════════════════════════════════════════════════════
build_taxa_lookup_from_ps <- function(ps) {
  tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  original_names <- rownames(tt)
  
  rank_prefixes <- c(
    Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__",
    Family  = "f__", Genus  = "g__", Species = "s__"
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
    result[is_miss] <- lefse_feature_names[is_miss]
  }
  unname(result)
}

# ══════════════════════════════════════════════════════════════
# 全 taxa BH 校正（与原 LEfSe 代码完全一致）
# ══════════════════════════════════════════════════════════════
bh_fullm <- function(marker_taxa, marker_pval, tested_taxa, q_thr = QVAL_THR) {
  tested_taxa <- unique(as.character(tested_taxa))
  if (length(tested_taxa) == 0) return(rep(NA_real_, length(marker_taxa)))
  
  p_full <- rep(1.0, length(tested_taxa))
  names(p_full) <- tested_taxa
  
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
  
  q_full  <- p.adjust(p_full, method = "BH")
  q_out   <- q_full[as.character(marker_taxa)]
  unname(q_out)
}

# ══════════════════════════════════════════════════════════════
# 核心函数：单次 LEfSe，返回显著 taxa 名向量 + m_tested
# 入参格式与 MaAsLin2/ANCOMBC 版一致：
#   feature_table : 行=样本，列=taxa（data.frame）
#   meta_df       : rownames=sample_id，含 phenotype 列
#   pair          : 长度2字符向量
# ══════════════════════════════════════════════════════════════
run_lefse_core <- function(feature_table,
                           meta_df,
                           pair,
                           output_dir = NULL,
                           top_k      = TOP_K) {
  
  # ── 对齐（与 MaAsLin2/ANCOMBC 完全相同）─────────────────
  common <- intersect(rownames(feature_table), rownames(meta_df))
  if (length(common) < 4) return(list(taxa = character(0), m_tested = 0L))
  
  feature_table <- feature_table[common, , drop = FALSE]
  meta_df       <- meta_df[common,       , drop = FALSE]
  
  # 只保留 pair 内样本
  pair_samples <- rownames(meta_df)[meta_df$phenotype %in% pair]
  if (length(pair_samples) < 4) return(list(taxa = character(0), m_tested = 0L))
  
  pair_feat <- feature_table[pair_samples, , drop = FALSE]
  pair_meta <- meta_df[pair_samples,       , drop = FALSE]
  
  # 各组样本数校验
  grp_counts <- table(pair_meta$phenotype)
  if (any(grp_counts < 2)) return(list(taxa = character(0), m_tested = 0L))
  
  # ── 转置为 taxa×samples，去全 0 行 ───────────────────────
  ab_mat      <- t(ensure_matrix(as.matrix(pair_feat)))
  ab_mat      <- ab_mat[rowSums(ab_mat) > 0, , drop = FALSE]
  ab_mat      <- ensure_matrix(ab_mat)
  if (nrow(ab_mat) == 0) return(list(taxa = character(0), m_tested = 0L))
  
  tested_taxa <- rownames(ab_mat)   # BH 的 m
  m_tested    <- length(tested_taxa)
  
  # 整理 metadata（phyloseq 要求 rownames = sample_id）
  pair_meta$phenotype <- factor(pair_meta$phenotype, levels = pair)
  
  # 解析 taxonomy
  taxonomy_df <- parse_taxonomy2(tested_taxa)
  
  # 构建 phyloseq
  ps <- tryCatch(
    phyloseq(
      otu_table(ab_mat, taxa_are_rows = TRUE),
      sample_data(pair_meta),
      tax_table(as.matrix(taxonomy_df))
    ),
    error = function(e) { cat("    [phyloseq 构建失败]", e$message, "\n"); NULL }
  )
  if (is.null(ps)) return(list(taxa = character(0), m_tested = m_tested))
  
  taxa_lookup <- build_taxa_lookup_from_ps(ps)
  
  # 运行 LEfSe
  # wilcoxon_cutoff / kw_cutoff 均设为 1.0，让所有特征进入 marker_table，
  # 事后由 bh_fullm() 对全部 m_tested 个特征做 BH 校正，再筛 q < QVAL_THR。
  # 小样本时 KW p 值普遍 > 0.05，预过滤会导致 marker_table 为空，
  # BH 校正也无从施展，因此必须先放宽阈值至 1.0。
  lefse_res <- tryCatch(
    suppressMessages(
      run_lefse(
        ps,
        wilcoxon_cutoff = 1.0,         # ← 改动：原 0.05，放宽为 1.0（预过滤关闭）
        kw_cutoff       = 1.0,         # ← 改动：原 0.05，放宽为 1.0（预过滤关闭）
        lda_cutoff      = 0,
        group           = "phenotype",
        multigrp_strat  = FALSE,
        taxa_rank       = "all"
      )
    ),
    error = function(e) { cat("    [LEfSe 失败]", e$message, "\n"); NULL }
  )
  if (is.null(lefse_res)) return(list(taxa = character(0), m_tested = m_tested))
  
  # 提取 marker_table
  markers <- tryCatch({
    mt <- lefse_res@marker_table
    df <- data.frame(mt@.Data, stringsAsFactors = FALSE)
    names(df) <- mt@names
    df
  }, error = function(e) NULL)
  
  if (is.null(markers) || nrow(markers) == 0)
    return(list(taxa = character(0), m_tested = m_tested))
  
  # 定位列名（与原 LEfSe 代码完全一致）
  get_col <- function(df, cands) {
    nm <- intersect(cands, names(df))
    if (!length(nm)) NULL else nm[1]
  }
  col_feature <- get_col(markers, c("feature", "Feature", "taxon", "Taxon"))
  col_lda     <- get_col(markers, c("ef_lda", "lda", "lda_score", "LDA", "ef_lda_mean"))
  col_p       <- get_col(markers, c("pvalue", "p_val", "pvalue_kw", "pvalue_wilcoxon",
                                    "p.value", "p_value"))
  
  if (is.null(col_feature)) return(list(taxa = character(0), m_tested = m_tested))
  
  markers$feature_raw <- as.character(markers[[col_feature]])
  markers$p_raw       <- if (!is.null(col_p))
    suppressWarnings(as.numeric(markers[[col_p]])) else NA_real_
  
  # 还原 taxa 名 + BH 全 m 校正
  markers$feature_restored <- restore_taxa_name(markers$feature_raw, taxa_lookup)
  markers$q_val <- bh_fullm(
    marker_taxa  = markers$feature_restored,
    marker_pval  = markers$p_raw,
    tested_taxa  = tested_taxa,
    q_thr        = QVAL_THR
  )
  
  feats <- markers %>%
    filter(!is.na(feature_restored), !is.na(q_val), q_val <= QVAL_THR)
  
  if (!is.null(col_lda)) {
    feats <- feats %>%
      mutate(lda_raw = suppressWarnings(as.numeric(.data[[col_lda]]))) %>%
      arrange(q_val, desc(lda_raw))
  } else {
    feats <- feats %>% arrange(q_val)
  }
  
  if (nrow(feats) == 0) return(list(taxa = character(0), m_tested = m_tested))
  feats <- slice_head(feats, n = top_k)
  
  # 落盘（仅全量运行时传入 output_dir）
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ref_group <- sort(pair)[1]
    nonref    <- setdiff(pair, ref_group)[1]
    result_full <- markers %>%
      mutate(ref_group = ref_group,
             pair      = paste(pair, collapse = "_vs_"))
    write_csv(result_full, file.path(output_dir, "results.csv"))
  }
  
  list(taxa = unique(as.character(feats$feature_restored)), m_tested = m_tested)
}

# ══════════════════════════════════════════════════════════════
# 单个 project 的完整 LEfSe + Jaccard 流程
# 数据读取 & 对齐逻辑与 MaAsLin2/ANCOMBC 版逐行对应
# ══════════════════════════════════════════════════════════════
run_project_jaccard <- function(proj, proj_idx) {
  
  cat("\n", SEP, "\n[Project ", proj_idx, "] ", proj, "\n", SEP, "\n", sep = "")
  
  proj_dir <- file.path(INPUT_DIR,  proj)
  out_proj <- file.path(OUTPUT_DIR, proj)
  dir.create(out_proj, recursive = TRUE, showWarnings = FALSE)
  
  counts_path <- file.path(proj_dir, "counts.tsv")
  meta_path   <- file.path(proj_dir, "metadata.csv")
  
  if (!file.exists(counts_path) || !file.exists(meta_path)) {
    cat("  ⚠️  缺少 counts.tsv 或 metadata.csv，跳过\n")
    return(NULL)
  }
  
  # ── 读取（与 MaAsLin2/ANCOMBC 完全相同）─────────────────
  counts_raw <- read.delim(counts_path, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE)
  meta_raw   <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  log_info(paste0("counts 维度: ", nrow(counts_raw), " features × ",
                  ncol(counts_raw), " samples"))
  log_info(paste0("metadata 行数: ", nrow(meta_raw)))
  
  # ── 清洗 counts（与 MaAsLin2/ANCOMBC 完全相同）──────────
  counts_mat <- ensure_matrix(counts_raw)
  counts_mat <- round(pmax(counts_mat, 0))
  counts_mat <- counts_mat[, colSums(counts_mat) > 0, drop = FALSE]
  counts_mat <- counts_mat[rowSums(counts_mat) > 0,   , drop = FALSE]
  
  # ── ★ 转置为行=样本，与 MaAsLin2/ANCOMBC 一致 ───────────
  feature_table <- as.data.frame(t(counts_mat))
  # feature_table: 行=samples，列=taxa
  
  # ── 整理 metadata（与 MaAsLin2/ANCOMBC 完全相同）────────
  meta <- meta_raw %>%
    select(all_of(c(SAMPLE_COL, GROUP_COL))) %>%
    rename(phenotype = all_of(GROUP_COL)) %>%
    mutate(phenotype = as.character(phenotype)) %>%
    filter(!is.na(phenotype), phenotype != "") %>%
    column_to_rownames(SAMPLE_COL)    # rownames = sample_id
  
  # ── 对齐（与 MaAsLin2/ANCOMBC 完全相同）─────────────────
  common <- intersect(rownames(meta), rownames(feature_table))
  if (length(common) < 4) {
    cat("  ⚠️  对齐后样本数不足 4，跳过\n"); return(NULL)
  }
  meta          <- meta[common, , drop = FALSE]
  feature_table <- feature_table[common, , drop = FALSE]
  
  # ── 过滤有效分组（与 MaAsLin2/ANCOMBC 完全相同）─────────
  grp_counts <- table(meta$phenotype)
  valid_grps <- names(grp_counts[grp_counts >= 2])
  if (length(valid_grps) < 2) {
    cat("  ⚠️  有效分组数不足 2，跳过\n"); return(NULL)
  }
  
  meta          <- meta[meta$phenotype %in% valid_grps, , drop = FALSE]
  feature_table <- feature_table[rownames(meta), , drop = FALSE]
  
  cat("  分组分布:\n"); print(table(meta$phenotype))
  
  # ── 枚举两两比较对 ────────────────────────────────────────
  all_grps      <- sort(unique(meta$phenotype))
  disease_pairs <- combn(all_grps, 2, simplify = FALSE)
  n_pairs_total <- length(disease_pairs)
  cat("  两两比较对数:", n_pairs_total, "\n")
  
  jaccard_records <- list()
  lefse_results   <- list()
  
  for (i in seq_along(disease_pairs)) {
    
    pair      <- disease_pairs[[i]]
    g1        <- pair[1]; g2 <- pair[2]
    comp_name <- paste0(g1, "_vs_", g2)
    cat("\n  [Pair", i, "/", n_pairs_total, "]", comp_name, "\n")
    
    pair_seed <- BASE_SEED + proj_idx * 1000000L + i * 100L
    
    # ── Step 1：全量运行（落盘供复查）────────────────────
    full_out_dir <- file.path(out_proj, "fulldata", comp_name)
    
    full_res <- tryCatch(
      run_lefse_core(
        feature_table = feature_table,
        meta_df       = meta,
        pair          = pair,
        output_dir    = full_out_dir    # 落盘
      ),
      error = function(e) {
        cat("    全量运行失败:", e$message, "\n")
        list(taxa = character(0), m_tested = 0L)
      }
    )
    
    S_all    <- full_res$taxa
    m_tested <- full_res$m_tested
    cat("    全量显著特征数:", length(S_all),
        " | m_tested:", m_tested, "\n")
    
    # 读取落盘结果汇总
    res_file <- file.path(full_out_dir, "results.csv")
    if (file.exists(res_file)) {
      res_df <- read_csv(res_file, show_col_types = FALSE) %>%
        mutate(project = proj, comparison = comp_name,
               group_a = g1,  group_b    = g2)
      lefse_results[[comp_name]] <- res_df
    }
    
    if (length(S_all) == 0) {
      jaccard_records[[length(jaccard_records) + 1]] <-
        empty_jaccard_row(proj, comp_name, n_pairs_total)
      next
    }
    
    # ── Step 2：B 次 bootstrap Jaccard ─────────────────
    jac_vals <- rep(NA_real_, B)
    
    for (b in seq_len(B)) {
      
      set.seed(pair_seed + b)
      
      # 分层抽样（与 MaAsLin2/ANCOMBC 完全相同）
      pair_samples_all <- rownames(meta)[meta$phenotype %in% pair]
      meta_pair        <- meta[pair_samples_all, , drop = FALSE]
      
      sub_ids <- meta_pair %>%
        rownames_to_column("sid") %>%
        group_by(phenotype) %>%
        slice_sample(prop = SUBSAMPLE) %>%
        ungroup() %>%
        pull(sid) %>%
        as.character()
      
      if (length(sub_ids) < 4) next
      
      feat_sub <- feature_table[sub_ids, , drop = FALSE]
      feat_sub <- feat_sub[, colSums(feat_sub) > 0, drop = FALSE]
      meta_sub <- meta[sub_ids, , drop = FALSE]
      
      if (ncol(feat_sub) == 0) next
      if (n_distinct(meta_sub$phenotype) < 2) next
      
      sub_res <- tryCatch(
        run_lefse_core(
          feature_table = feat_sub,
          meta_df       = meta_sub,
          pair          = pair,
          output_dir    = NULL          # 临时，不落盘
        ),
        error = function(e) list(taxa = character(0), m_tested = 0L)
      )
      
      jac_vals[b] <- jaccard(sub_res$taxa, S_all)
    }
    
    jac_clean     <- jac_vals[!is.na(jac_vals)]
    n_boots_valid <- length(jac_clean)
    cat("    有效 bootstrap 次数:", n_boots_valid, "/", B, "\n")
    if (n_boots_valid > 0)
      cat(sprintf("    Jaccard 均值: %.3f  中位数: %.3f\n",
                  mean(jac_clean), median(jac_clean)))
    
    jaccard_records[[length(jaccard_records) + 1]] <- tibble(
      project         = proj,
      comparison      = comp_name,
      n_pairs         = as.integer(n_pairs_total),
      top_k           = as.integer(TOP_K),
      n_features_full = as.integer(length(S_all)),
      m_tested        = as.integer(m_tested),
      n_boots         = as.integer(B),
      n_boots_valid   = as.integer(n_boots_valid),
      jaccard_mean    = if (n_boots_valid == 0) NA_real_ else mean(jac_clean),
      jaccard_median  = if (n_boots_valid == 0) NA_real_ else median(jac_clean),
      jaccard_sd      = if (n_boots_valid == 0) NA_real_ else sd(jac_clean),
      jaccard_min     = if (n_boots_valid == 0) NA_real_ else min(jac_clean),
      jaccard_max     = if (n_boots_valid == 0) NA_real_ else max(jac_clean)
    )
  }
  
  # ── 保存 project 结果 ─────────────────────────────────────
  proj_jaccard <- bind_rows(jaccard_records)
  write_csv(proj_jaccard, file.path(out_proj, "jaccard_results.csv"))
  
  if (length(lefse_results) > 0) {
    write_csv(bind_rows(lefse_results),
              file.path(out_proj, "all_comparisons_results.csv"))
  }
  
  log_info(paste0("✅ ", proj, " 完成，结果保存至 ", out_proj))
  proj_jaccard
}

# ══════════════════════════════════════════════════════════════
# 主流程：遍历所有 project
# ══════════════════════════════════════════════════════════════
log_step("MAIN", "代表性验证 LEfSe Jaccard 分析")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

all_projects <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
all_projects <- all_projects[all_projects != ""]
cat("检测到", length(all_projects), "个 project:",
    paste(all_projects, collapse = ", "), "\n")

all_jaccard  <- list()
failed_projs <- character(0)

for (idx in seq_along(all_projects)) {
  proj <- all_projects[idx]
  result <- tryCatch(
    run_project_jaccard(proj, idx),
    error = function(e) {
      cat("  ❌ Project", proj, "失败:", e$message, "\n")
      failed_projs <<- c(failed_projs, proj)
      NULL
    }
  )
  if (!is.null(result)) all_jaccard[[proj]] <- result
}

# ══════════════════════════════════════════════════════════════
# 全局汇总（结构与 MaAsLin2/ANCOMBC 版完全一致）
# ══════════════════════════════════════════════════════════════
log_step("SUMMARY", "全局 Jaccard 汇总")

if (length(all_jaccard) > 0) {
  
  # ── Raw 明细 ─────────────────────────────────────────────
  jaccard_raw <- bind_rows(all_jaccard) %>%
    mutate(across(starts_with("jaccard_"),
                  ~ ifelse(is.nan(.), NA_real_, .)))
  
  write_csv(jaccard_raw, file.path(OUTPUT_DIR, "jaccard_raw.csv"))
  saveRDS(jaccard_raw,   file.path(OUTPUT_DIR, "jaccard_raw.rds"))
  
  cat("\n── Raw Jaccard 明细（前10行）──\n")
  print(head(jaccard_raw, 10))
  
  # ── Project 级汇总 ───────────────────────────────────────
  jaccard_proj_summary <- jaccard_raw %>%
    group_by(project) %>%
    summarise(
      n_pairs_run          = sum(!is.na(jaccard_mean)),
      n_pairs_total        = max(n_pairs),
      n_features_full_mean = mean(n_features_full, na.rm = TRUE),
      m_tested_mean        = mean(m_tested,        na.rm = TRUE),
      n_boots_valid_mean   = mean(n_boots_valid,   na.rm = TRUE),
      jaccard_mean         = mean(jaccard_mean,    na.rm = TRUE),
      jaccard_median       = mean(jaccard_median,  na.rm = TRUE),
      jaccard_sd_mean      = mean(jaccard_sd,      na.rm = TRUE),
      jaccard_min          = min(jaccard_min,       na.rm = TRUE),
      jaccard_max          = max(jaccard_max,       na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(jaccard_proj_summary, file.path(OUTPUT_DIR, "jaccard_project_summary.csv"))
  saveRDS(jaccard_proj_summary,   file.path(OUTPUT_DIR, "jaccard_project_summary.rds"))
  
  cat("\n── Project 级 Jaccard 汇总 ──\n")
  print(jaccard_proj_summary)
  
  # ── Comparison 级汇总（跨 project）─────────────────────
  jaccard_comp_summary <- jaccard_raw %>%
    group_by(comparison) %>%
    summarise(
      n_projects         = n(),
      jaccard_mean       = mean(jaccard_mean,   na.rm = TRUE),
      jaccard_median     = mean(jaccard_median, na.rm = TRUE),
      jaccard_sd_between = sd(jaccard_mean,     na.rm = TRUE),
      jaccard_min        = min(jaccard_min,      na.rm = TRUE),
      jaccard_max        = max(jaccard_max,      na.rm = TRUE),
      m_tested_mean      = mean(m_tested,        na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(jaccard_mean))
  
  write_csv(jaccard_comp_summary, file.path(OUTPUT_DIR, "jaccard_comparison_summary.csv"))
  saveRDS(jaccard_comp_summary,   file.path(OUTPUT_DIR, "jaccard_comparison_summary.rds"))
  
  cat("\n── Comparison 级 Jaccard 汇总 ──\n")
  print(jaccard_comp_summary)
  
} else {
  cat("⚠️  所有 project 均未产生有效 Jaccard 结果\n")
}

# ══════════════════════════════════════════════════════════════
# 最终文件清单
# ══════════════════════════════════════════════════════════════
cat("\n", SEP, "\n完成！输出文件清单：\n", SEP, "\n", sep = "")

all_out_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in all_out_files) {
  size_kb  <- round(file.size(f) / 1024, 1)
  rel_path <- sub(paste0(normalizePath(OUTPUT_DIR, winslash = "/"), "/"), "",
                  normalizePath(f, winslash = "/"))
  cat(sprintf("  %-60s %7.1f KB\n", rel_path, size_kb))
}

cat(sprintf("\n  输出根目录:   %s\n",  normalizePath(OUTPUT_DIR)))
cat(sprintf("  成功 project: %d / %d\n", length(all_jaccard), length(all_projects)))
