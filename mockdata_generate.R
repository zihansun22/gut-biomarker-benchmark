
# 微生物组Mock数据生成脚本
# 专注于：Step 2 特征提取 + Step 3-4 差异标志物植入

library(dplyr)
library(readr)
library(MCMCpack)  # for rdirichlet

# =============================================================================
# 配置参数
# =============================================================================
CONFIG <- list(
  # 输入文件
  metadata_file = "metadata_genus.tsv",
  counts_file = "counts_genus.tsv",
  
  # Mock参数
  n_mock_samples = 200,  # 每组样本数
  n_replicates = 3,      # 每个场景重复次数
  
  # 差异设置
  n_diff_taxa = c(150, 300, 500),      # 差异taxa数量
  
  # 输出
  output_dir = "mock_benchmark_data"
)

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir,
                                               recursive = TRUE)



# =============================================================================
# Step 2: 提取真实数据特征
# =============================================================================

extract_real_features <- function() {
  cat("Step 2: 提取真实数据特征...\n")
  
  # 读取数据
  metadata <- read_tsv(CONFIG$metadata_file)
  counts_raw <- read_tsv(CONFIG$counts_file)
  
  # 转换为样本×taxa矩阵
  count_matrix <- as.matrix(counts_raw[, -1])
  rownames(count_matrix) <- counts_raw[[1]]
  count_matrix <- t(count_matrix)  # 样本为行
  
  # 1. 测序深度分布
  library_sizes <- rowSums(count_matrix)
  
  # 2. 计算相对丰度
  rel_abundance <- count_matrix / rowSums(count_matrix)     # rel = relative
  
  # 3. Taxa级特征
  taxa_features <- data.frame(
    taxa = colnames(count_matrix),
    prevalence = colMeans(count_matrix > 0),           # p_i
    mean_rel_abund = colMeans(rel_abundance),          # mu_i (相对丰度)
    zero_prop = colMeans(count_matrix == 0),
    stringsAsFactors = FALSE
  )
  
  # 4. 估计Dirichlet-Multinomial参数
  # 使用矩方法估计总alpha参数
  overall_mean_rel <- colMeans(rel_abundance)       # 每个物种的平均相对丰度
  overall_var_rel <- apply(rel_abundance, 2, var)   # 每个物种的方差
  
  # 计算离散参数（计算均值-方差比率）
  mean_variance_ratio <- mean(overall_var_rel / (overall_mean_rel * (1 - overall_mean_rel)), na.rm = TRUE)
  alpha_sum <- max(1, 1 / (mean_variance_ratio - 1/nrow(count_matrix)))  # 防止过小
  
  # 个体alpha参数
  alpha_vec <- overall_mean_rel * alpha_sum
  alpha_vec[alpha_vec < 0.001] <- 0.001  # 防止过小
  
  features <- list(
    count_matrix = count_matrix,
    metadata = metadata,
    library_sizes = library_sizes,
    taxa_features = taxa_features,
    alpha_vec = alpha_vec,
    alpha_sum = alpha_sum
  )
  
  cat("  样本数:", nrow(count_matrix), "| Taxa数:", ncol(count_matrix), "\n")
  cat("  平均测序深度:", round(mean(library_sizes)), "\n")
  cat("  平均零比例:", round(mean(taxa_features$zero_prop), 3), "\n")
  cat("  Alpha总和:", round(alpha_sum, 2), "\n\n")
  
  return(features)
}


# =============================================================================
# Step 3: 定义Benchmark场景（修订版）
# =============================================================================

define_scenarios <- function(features) {
  cat("Step 3: 定义benchmark场景...\n")
  
  n_taxa <- length(features$taxa_features$taxa)
  
  scenarios <- list()
  
  # 场景A: Realistic baseline - 无混杂、中等效应
  scenarios$A_no_confound <- list(
    name = "A_no_confound",
    description = "Realistic baseline: 97-98% 稀疏，中等效应",
    n_diff_taxa = round(n_taxa * 0.09),      # ~9% taxa 真差异（~300个）
    log2_fold_change = 1.2,                  # log2FC = 1.2（对应 2.3倍）
    confound_batch = FALSE,
    target_sparsity = 0.975,                 # 目标97.5%零值率
    target_median_depth = 12000,
    high_var_taxa_fraction = 0.27            # 27%高变异taxa
  )
  
  # 场景B: Batch confound - 中等效应 + batch偏移
  scenarios$B_batch_confound <- list(
    name = "B_batch_confound",
    description = "Batch confound: 97-98% 稀疏 + batch偏移",
    n_diff_taxa = round(n_taxa * 0.09),      # ~9% taxa 真差异
    log2_fold_change = 1.2,
    confound_batch = TRUE,
    target_sparsity = 0.975,
    target_median_depth = 12000,
    high_var_taxa_fraction = 0.27,
    batch_effect = list(
      n_batch_taxa = round(n_taxa * 0.05),   # 5% taxa受batch影响
      batch_depth_multiplier = c(1.3, 0.7)   # batch1: ×1.3, batch2: ×0.7
    )
  )
  
  # 场景C: 极端稀疏 - 小效应
  scenarios$C_sparse_small <- list(
    name = "C_sparse_small",
    description = "极端稀疏: 99-99.5% 稀疏，小效应",
    n_diff_taxa = round(n_taxa * 0.045),     # ~4.5% taxa 真差异（~150个）
    log2_fold_change = 0.6,                  # log2FC = 0.6（对应 1.5倍）
    confound_batch = FALSE,
    target_sparsity = 0.992,                 # 目标99.2%零值率
    target_median_depth = 12000,
    high_var_taxa_fraction = 0.15,           # 15%高变异taxa
    extreme_prevalence_cutoff = 0.05         # 多数taxa <5%样本出现
  )
  
  # 场景D: 大效应 - 保持稀疏但效应明显
  scenarios$D_large_effect <- list(
    name = "D_large_effect",
    description = "大效应: 97-98% 稀疏，logFC ≥ 2",
    n_diff_taxa = round(n_taxa * 0.12),      # ~12% taxa 真差异（~390个）
    log2_fold_change = 2.0,                  # log2FC = 2.0（对应 4倍）
    confound_batch = FALSE,
    target_sparsity = 0.975,
    target_median_depth = 12000,
    high_var_taxa_fraction = 0.35            # 35%高变异taxa
  )
  
  cat("  定义了", length(scenarios), "个benchmark场景\n\n")
  return(scenarios)
}


# =============================================================================
# Step 4: 生成Mock数据核心函数（修订版）
# =============================================================================

generate_mock_counts <- function(features, scenario, replicate_id = 1) {
  
  n_samples_per_group <- CONFIG$n_mock_samples
  n_taxa <- length(features$taxa_features$taxa)
  group_names <- c("Healthy", "Disease_mild", "Disease_severe")
  total_samples <- n_samples_per_group * 3
  
  # ========== 1. 选择差异taxa ==========
  taxa_by_abund <- order(features$taxa_features$mean_rel_abund, decreasing = TRUE)
  # 从中等丰度范围选择（避免极稀有和极丰富的）
  candidate_taxa <- taxa_by_abund[round(n_taxa*0.05):round(n_taxa*0.4)]
  
  n_total <- scenario$n_diff_taxa
  n_shared <- round(n_total * 0.4)
  n_prog <- round(n_total * 0.3)
  n_stage <- n_total - n_shared - n_prog
  
  sampled_all <- sample(candidate_taxa, n_shared + n_prog + n_stage)
  shared_idx <- sampled_all[1:n_shared]
  prog_idx <- sampled_all[(n_shared+1):(n_shared+n_prog)]
  stage_idx <- sampled_all[(n_shared+n_prog+1):length(sampled_all)]
  
  shared_names <- features$taxa_features$taxa[shared_idx]
  prog_names <- features$taxa_features$taxa[prog_idx]
  stage_names <- features$taxa_features$taxa[stage_idx]
  
  # ========== 2. 效应方向 ==========
  make_direction <- function(n) {
    n_up <- ceiling(n / 2)
    sample(c(rep("up", n_up), rep("down", n - n_up)))
  }
  shared_dir <- make_direction(n_shared)
  prog_dir <- make_direction(n_prog)
  stage_dir <- make_direction(n_stage)
  
  # ========== 3. 转换log2FC为fold change ==========
  fold_change <- 2^scenario$log2_fold_change
  
  # ========== 4. 基线alpha参数 ==========
  alpha_healthy <- pmax(features$alpha_vec + 0.05, 0.01)
  
  # ========== 5. 应用效应 ==========
  apply_effect <- function(alpha_base, taxa_indices, directions, fc) {
    alpha_new <- alpha_base
    for (i in seq_along(taxa_indices)) {
      idx <- taxa_indices[i]
      if (directions[i] == "up") {
        alpha_new[idx] <- alpha_new[idx] * fc
      } else {
        alpha_new[idx] <- alpha_new[idx] / fc
      }
    }
    alpha_new / sum(alpha_new) * sum(alpha_base)
  }
  
  alpha_mild <- apply_effect(
    alpha_healthy,
    c(shared_idx, stage_idx),
    c(shared_dir, stage_dir),
    fold_change
  )
  
  alpha_severe <- apply_effect(
    alpha_healthy,
    c(shared_idx, prog_idx),
    c(shared_dir, prog_dir),
    fold_change
  )
  
  alpha_list <- list(
    Healthy = alpha_healthy,
    Disease_mild = alpha_mild,
    Disease_severe = alpha_severe
  )
  
  # ========== 6. 调整库大小以匹配目标深度 ==========
  target_depth <- scenario$target_median_depth
  lib_sizes_scaled <- features$library_sizes * (target_depth / median(features$library_sizes))
  good_lib_sizes <- lib_sizes_scaled[lib_sizes_scaled > quantile(lib_sizes_scaled, 0.3)]
  
  # ========== 7. 生成计数矩阵 ==========
  mock_counts <- matrix(0, nrow = total_samples, ncol = n_taxa)
  colnames(mock_counts) <- features$taxa_features$taxa
  
  sample_groups <- c()
  sample_names <- c()
  sample_idx <- 1
  
  for (group_name in names(alpha_list)) {
    alpha_params <- alpha_list[[group_name]]
    
    for (i in 1:n_samples_per_group) {
      lib_size <- sample(good_lib_sizes, 1)
      # 先从 Dirichlet 采样得到这个样本的组成比例
      p <- rdirichlet(1, alpha_params)
      # 再用这个 p 生成计数
      mock_counts[sample_idx, ] <- rmultinom(1, lib_size, p)
      
      sample_groups <- c(sample_groups, group_name)
      sample_names <- c(sample_names, 
                        paste0("Mock_", scenario$name, "_R", replicate_id, "_", group_name, "_S", i))
      sample_idx <- sample_idx + 1
    }
  }
  
  # ========== 8. 应用零膨胀以达到目标稀疏性 ==========
  all_diff_idx <- c(shared_idx, prog_idx, stage_idx)
  mock_counts <- apply_zero_inflation_v2(
    mock_counts, 
    features$taxa_features, 
    scenario, 
    all_diff_idx, 
    sample_groups
  )
  
  # ========== 9. 生成metadata ==========
  mock_metadata <- data.frame(
    sample_id = sample_names,
    phenotype = sample_groups,
    scenario = scenario$name,
    replicate = replicate_id,
    library_size = rowSums(mock_counts),
    stringsAsFactors = FALSE
  )
  
  # ========== 10. 添加batch混杂 ==========
  if (scenario$confound_batch) {
    mock_metadata$batch <- "batch1"
    
    healthy_idx <- which(mock_metadata$phenotype == "Healthy")
    healthy_batch2_idx <- sample(healthy_idx, round(length(healthy_idx) * 0.3))
    mock_metadata$batch[healthy_batch2_idx] <- "batch2"
    
    disease_idx <- which(mock_metadata$phenotype != "Healthy")
    disease_batch2_idx <- sample(disease_idx, round(length(disease_idx) * 0.7))
    mock_metadata$batch[disease_batch2_idx] <- "batch2"
    
    # 对batch2样本应用batch effect
    batch2_idx <- which(mock_metadata$batch == "batch2")
    batch_taxa_idx <- sample(setdiff(1:n_taxa, all_diff_idx), 
                             scenario$batch_effect$n_batch_taxa)
    
    for (idx in batch2_idx) {
      group_name <- mock_metadata$phenotype[idx]
      alpha_params <- alpha_list[[group_name]]
      
      # 应用batch effect：改变特定taxa的丰度
      alpha_batch <- alpha_params
      for (b_idx in batch_taxa_idx) {
        alpha_batch[b_idx] <- alpha_batch[b_idx] * scenario$batch_effect$batch_depth_multiplier[1]
      }
      alpha_batch <- alpha_batch / sum(alpha_batch) * sum(alpha_params)
      
      # 使用 Dirichlet-Multinomial 生成 batch effect 样本
      p_batch <- rdirichlet(1, alpha_batch)
      lib_size <- mock_metadata$library_size[idx]
      mock_counts[idx, ] <- rmultinom(1, lib_size, p_batch)
    }
    
    mock_metadata$library_size <- rowSums(mock_counts)
    
    # 对batch2样本也应用零膨胀
    mock_counts <- apply_zero_inflation_v2(
      mock_counts, 
      features$taxa_features, 
      scenario, 
      all_diff_idx, 
      sample_groups
    )
  }
  
  # ========== 11. 生成truth table ==========
  truth_table <- bind_rows(
    data.frame(
      taxa = shared_names,
      direction = shared_dir,
      log2_fold_change = scenario$log2_fold_change,
      fold_change = ifelse(shared_dir == "up", fold_change, 1/fold_change),
      marker_type = "shared",
      present_in_groups = "Disease_mild,Disease_severe",
      scenario = scenario$name,
      replicate = replicate_id,
      stringsAsFactors = FALSE
    ),
    data.frame(
      taxa = prog_names,
      direction = prog_dir,
      log2_fold_change = scenario$log2_fold_change,
      fold_change = ifelse(prog_dir == "up", fold_change, 1/fold_change),
      marker_type = "progression",
      present_in_groups = "Disease_severe",
      scenario = scenario$name,
      replicate = replicate_id,
      stringsAsFactors = FALSE
    ),
    data.frame(
      taxa = stage_names,
      direction = stage_dir,
      log2_fold_change = scenario$log2_fold_change,
      fold_change = ifelse(stage_dir == "up", fold_change, 1/fold_change),
      marker_type = "stage_specific",
      present_in_groups = "Disease_mild",
      scenario = scenario$name,
      replicate = replicate_id,
      stringsAsFactors = FALSE
    )
  )
  
  if (scenario$confound_batch) {
    batch_truth <- data.frame(
      taxa = features$taxa_features$taxa[batch_taxa_idx],
      direction = rep("batch_effect", length(batch_taxa_idx)),
      log2_fold_change = NA,
      fold_change = rep(scenario$batch_effect$batch_depth_multiplier[1], length(batch_taxa_idx)),
      marker_type = "batch_effect",
      present_in_groups = "batch2",
      scenario = scenario$name,
      replicate = replicate_id,
      stringsAsFactors = FALSE
    )
    truth_table <- bind_rows(truth_table, batch_truth)
  }
  
  return(list(counts = mock_counts, metadata = mock_metadata, truth = truth_table))
}


# =============================================================================
# 辅助函数：应用零膨胀 v2（达到目标稀疏性）
# =============================================================================

apply_zero_inflation_v2 <- function(counts, taxa_features, scenario, diff_taxa_idx, sample_groups) {
  
  n_samples <- nrow(counts)
  n_taxa <- ncol(counts)
  
  # 计算当前稀疏性
  current_sparsity <- sum(counts == 0) / (n_samples * n_taxa)
  target_sparsity <- scenario$target_sparsity
  
  # 如果已经达到目标，直接返回
  if (current_sparsity >= target_sparsity) {
    return(counts)
  }
  
  # 需要增加的零值数量
  zeros_needed <- round((target_sparsity - current_sparsity) * n_samples * n_taxa)
  
  # 优先在非差异taxa中增加零值
  non_diff_taxa <- setdiff(1:n_taxa, diff_taxa_idx)
  
  # 按prevalence排序，优先处理低prevalence的taxa
  prevalence <- colMeans(counts > 0)
  non_diff_prevalence <- prevalence[non_diff_taxa]
  sorted_non_diff <- non_diff_taxa[order(non_diff_prevalence)]
  
  zeros_added <- 0
  
  for (taxa_idx in sorted_non_diff) {
    if (zeros_added >= zeros_needed) break
    
    nonzero_idx <- which(counts[, taxa_idx] > 0)
    if (length(nonzero_idx) == 0) next
    
    # 计算该taxa需要增加的零值
    current_prev <- length(nonzero_idx) / n_samples
    target_prev <- taxa_features$prevalence[taxa_idx] * 0.5  # 降低到原prevalence的50%
    
    n_zeros_for_taxa <- round(n_samples * (current_prev - target_prev))
    n_zeros_for_taxa <- min(n_zeros_for_taxa, length(nonzero_idx) - 1)  # 至少保留1个非零
    
    if (n_zeros_for_taxa > 0) {
      # 优先将小的count变为0
      small_count_idx <- nonzero_idx[order(counts[nonzero_idx, taxa_idx])][1:n_zeros_for_taxa]
      counts[small_count_idx, taxa_idx] <- 0
      zeros_added <- zeros_added + n_zeros_for_taxa
    }
  }
  
  # 场景C：极端稀疏处理
  if (!is.null(scenario$extreme_prevalence_cutoff)) {
    cutoff <- scenario$extreme_prevalence_cutoff
    for (taxa_idx in non_diff_taxa) {
      current_prev <- mean(counts[, taxa_idx] > 0)
      if (current_prev > cutoff) {
        # 强制降低到cutoff以下
        nonzero_idx <- which(counts[, taxa_idx] > 0)
        target_n_nonzero <- round(n_samples * cutoff * 0.8)
        if (length(nonzero_idx) > target_n_nonzero) {
          zeros_idx <- sample(nonzero_idx, length(nonzero_idx) - target_n_nonzero)
          counts[zeros_idx, taxa_idx] <- 0
        }
      }
    }
  }
  
  return(counts)
}


# =============================================================================
# 主执行函数
# =============================================================================
generate_all_mock_data <- function(seed = 123) {
  # 在函数开头设置种子
  set.seed(seed)
  cat("=== 微生物组Mock数据生成流水线 ===\n\n")
  
  # Step 2: 提取特征
  features <- extract_real_features()
  
  # Step 3: 定义场景
  scenarios <- define_scenarios(features)
  
  # Step 4: 生成Mock数据
  cat("Step 4: 生成Mock数据...\n")
  
  all_results <- list()
  
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    cat("  生成场景:", scenario$name, "\n")
    
    scenario_results <- list()
    
    for (rep in 1:CONFIG$n_replicates) {
      cat("    重复", rep, "/", CONFIG$n_replicates, "\n")
      
      mock_data <- generate_mock_counts(features, scenario, rep)
      
      rep_name <- paste0(scenario$name, "_rep", rep)
      scenario_results[[rep_name]] <- mock_data
      
      # 保存文件
      counts_file <- file.path(CONFIG$output_dir, paste0("mock_counts_", rep_name, ".tsv"))
      metadata_file <- file.path(CONFIG$output_dir, paste0("mock_metadata_", rep_name, ".tsv"))
      truth_file <- file.path(CONFIG$output_dir, paste0("truth_table_", rep_name, ".tsv"))
      
      counts_df <- data.frame(taxa = features$taxa_features$taxa, 
                              t(mock_data$counts), check.names = FALSE)
      colnames(counts_df)[-1] <- mock_data$metadata$sample_id
      
      write_tsv(counts_df, counts_file)
      write_tsv(mock_data$metadata, metadata_file)  
      write_tsv(mock_data$truth, truth_file)
    }
    
    all_results[[scenario_name]] <- scenario_results
  }
  
  # 合并所有truth tables
  all_truth <- do.call(rbind, lapply(all_results, function(scenario) {
    do.call(rbind, lapply(scenario, function(rep) rep$truth))
  }))
  
  write_tsv(all_truth, file.path(CONFIG$output_dir, "all_truth_tables.tsv"))
  
  # 保存总结信息
  summary_info <- list(
    config = CONFIG,
    scenarios = scenarios,
    group_design = "Healthy / Disease_mild / Disease_severe",  # 新增说明
    marker_types = c("shared", "progression", "stage_specific"),
    random_seed = seed,  # 记录使用的种子
    features_summary = list(
      n_samples = nrow(features$count_matrix),
      n_taxa = ncol(features$count_matrix),
      mean_lib_size = mean(features$library_sizes),
      alpha_sum = features$alpha_sum
    ),
    generation_time = Sys.time()
  )
  
  saveRDS(summary_info, file.path(CONFIG$output_dir, "mock_generation_summary.rds"))
  
  cat("\n=== Mock数据生成完成! ===\n")
  cat("输出目录:", CONFIG$output_dir, "\n")
  cat("生成场景数:", length(scenarios), "x", CONFIG$n_replicates, "重复\n")
  cat("总差异taxa记录数:", nrow(all_truth), "\n")
  cat("组别设计: Healthy → Disease_mild → Disease_severe\n")
  cat("Marker类型: shared / progression / stage_specific\n")
  
  cat("主要输出文件:\n")
  cat("- mock_counts_[scenario]_rep[N].tsv: Mock计数矩阵\n")  
  cat("- mock_metadata_[scenario]_rep[N].tsv: 样本元数据（phenotype列含Healthy和疾病名）\n")
  cat("- truth_table_[scenario]_rep[N].tsv: 真实差异taxa\n")
  cat("- all_truth_tables.tsv: 所有truth table合集\n")
  cat("- mock_generation_summary.rds: 生成参数和总结\n\n")
  
  return(all_results)
}

if (!exists("skip_run")) {
  results <- generate_all_mock_data()
}



# 加载包
library(readr)

# 统一的相对丰度生成函数
generate_relative_abundance <- function(output_dir = "mock_benchmark_data") {
  
  # 获取所有counts文件
  counts_files <- list.files(output_dir, pattern = "mock_counts_.*\\.tsv", full.names = TRUE)
  
  for (counts_file in counts_files) {
    # 读取counts
    counts_df <- read_tsv(counts_file, show_col_types = FALSE)
    counts_matrix <- as.matrix(counts_df[, -1])
    rownames(counts_matrix) <- counts_df$taxa
    
    # 计算相对丰度
    relab_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
    
    # 转换回数据框格式
    relab_df <- data.frame(
      taxa = counts_df$taxa,
      relab_matrix,
      check.names = FALSE
    )
    
    # 保存相对丰度文件
    relab_file <- gsub("mock_counts_", "mock_relab_", basename(counts_file))
    write_tsv(relab_df, file.path(output_dir, relab_file))
  }
  
  cat("生成了", length(counts_files), "个相对丰度文件\n")
}

# 运行生成相对丰度
generate_relative_abundance()

# 生成了 12 个相对丰度文件 


# 真值标准 (Ground Truth)
truth_files <- list(
  individual = "truth_table_[scenario]_rep[N].tsv",  # 单个实验真值
  combined = "all_truth_tables.tsv"                  # 所有实验真值合并
)

# 评估
# 修改：新增 comparison_pair 参数，用于按比较对筛选真值
evaluate_method_performance <- function(predictions, scenario, rep_id,
                                        comparison_pair,   # 新增，如 c("Healthy","Disease_mild")
                                        all_taxa) {
  
  truth_file <- paste0("truth_table_", scenario, "_rep", rep_id, ".tsv")
  truth_df   <- read_tsv(file.path("mock_benchmark_data", truth_file), show_col_types = FALSE)
  
  # 构建真值向量
  # 修改：按比较对筛选真值
  # 原来：true_positive_taxa <- truth_df$taxa（所有差异taxa都算真值）
  # 现在：只有在这个比较对中"应该有差异"的taxa才算真值
  group_a <- comparison_pair[1]
  group_b <- comparison_pair[2]
  
  # 根据 present_in_groups 判断该taxa在这个比较对中是否应该有差异
  true_positive_taxa <- truth_df %>%
    filter(
      # 两组中至少一组在 present_in_groups 里，且另一组不在（或是Healthy）
      sapply(present_in_groups, function(pg) {
        groups_with_effect <- unlist(strsplit(pg, ","))
        # 比较对中有一组有效应、另一组没有 → 这个taxa在此比较对中是真阳性
        xor(group_a %in% groups_with_effect, group_b %in% groups_with_effect) |
          # 两组都有效应但方向可能不同（shared taxa在两个疾病组比较时不是marker）
          # Healthy vs 任何疾病组：只要疾病组有效应就算
          (group_a == "Healthy" & group_b %in% groups_with_effect) |
          (group_b == "Healthy" & group_a %in% groups_with_effect)
      })
    ) %>%
    pull(taxa)
  
  true_labels <- ifelse(all_taxa %in% true_positive_taxa, 1, 0)
  
  # 计算性能指标
  performance <- calculate_metrics(
    true_labels = true_labels,
    predictions = predictions,  # 方法输出的p值或重要性分数
    method_name = "LefSe"      # 或其他方法名
  )
  
  # 新增：记录比较对信息
  performance$comparison <- paste(comparison_pair, collapse = "_vs_")
  
  return(performance)
}

calculate_metrics <- function(true_labels, predictions, method_name) {
  # ROC/AUPRC
  roc_obj <- roc(true_labels, -log10(predictions$p_value))  # p值越小排名越高
  auc_roc <- auc(roc_obj)
  auc_prc <- pr.curve(predictions$p_value[true_labels == 1], 
                      predictions$p_value[true_labels == 0])$auc.integral
  
  # 固定阈值性能
  sig_taxa <- predictions$p_value < 0.05
  tp <- sum(true_labels == 1 & sig_taxa)
  fp <- sum(true_labels == 0 & sig_taxa)
  fn <- sum(true_labels == 1 & !sig_taxa)
  tn <- sum(true_labels == 0 & !sig_taxa)
  
  metrics <- list(
    method = method_name,
    auroc = auc_roc,
    auprc = auc_prc,
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    precision = tp / (tp + fp),
    f1 = 2 * tp / (2 * tp + fp + fn),
    fdr = fp / (tp + fp)
  )
  
  return(metrics)
}



