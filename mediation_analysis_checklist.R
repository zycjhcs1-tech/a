# NHANES中介效应分析问题诊断检查清单
# 快速诊断工具

# =============================================================================
# 快速诊断函数
# =============================================================================

diagnose_mediation_analysis <- function(data, exposure_var, mediator_var, outcome_var, 
                                       covariates = NULL, weights = NULL) {
  
  cat("=== NHANES中介效应分析诊断报告 ===\n")
  cat("分析时间:", Sys.time(), "\n\n")
  
  # 1. 数据基本信息
  cat("1. 数据基本信息:\n")
  cat("   总样本量:", nrow(data), "\n")
  
  # 检查缺失值
  missing_exposure <- sum(is.na(data[[exposure_var]]))
  missing_mediator <- sum(is.na(data[[mediator_var]]))
  missing_outcome <- sum(is.na(data[[outcome_var]]))
  
  cat("   暴露变量缺失:", missing_exposure, "(", round(missing_exposure/nrow(data)*100, 1), "%)\n")
  cat("   中介变量缺失:", missing_mediator, "(", round(missing_mediator/nrow(data)*100, 1), "%)\n")
  cat("   结果变量缺失:", missing_outcome, "(", round(missing_outcome/nrow(data)*100, 1), "%)\n")
  
  # 完整案例数
  complete_cases <- sum(complete.cases(data[c(exposure_var, mediator_var, outcome_var)]))
  cat("   完整案例数:", complete_cases, "(", round(complete_cases/nrow(data)*100, 1), "%)\n")
  
  # 2. 变量分布检查
  cat("\n2. 变量分布检查:\n")
  
  # 暴露变量
  exp_summary <- summary(data[[exposure_var]])
  cat("   暴露变量 (", exposure_var, "):\n")
  print(exp_summary)
  
  # 检查偏态
  if(require(moments, quietly = TRUE)) {
    exp_skew <- skewness(data[[exposure_var]], na.rm = TRUE)
    cat("   偏度:", round(exp_skew, 3))
    if(abs(exp_skew) > 2) cat(" [严重偏态，建议转换]")
    cat("\n")
  }
  
  # 中介变量
  med_summary <- summary(data[[mediator_var]])
  cat("   中介变量 (", mediator_var, "):\n")
  print(med_summary)
  
  if(require(moments, quietly = TRUE)) {
    med_skew <- skewness(data[[mediator_var]], na.rm = TRUE)
    cat("   偏度:", round(med_skew, 3))
    if(abs(med_skew) > 2) cat(" [严重偏态，建议转换]")
    cat("\n")
  }
  
  # 结果变量
  out_summary <- summary(data[[outcome_var]])
  cat("   结果变量 (", outcome_var, "):\n")
  print(out_summary)
  
  # 3. 相关性检查
  cat("\n3. 关键关联检查:\n")
  
  # 暴露-结果关联
  cor_exp_out <- cor(data[[exposure_var]], data[[outcome_var]], use = "complete.obs")
  cat("   暴露-结果相关系数:", round(cor_exp_out, 4))
  if(abs(cor_exp_out) < 0.05) cat(" [关联很弱，可能无总效应]")
  cat("\n")
  
  # 暴露-中介关联
  cor_exp_med <- cor(data[[exposure_var]], data[[mediator_var]], use = "complete.obs")
  cat("   暴露-中介相关系数:", round(cor_exp_med, 4))
  if(abs(cor_exp_med) < 0.05) cat(" [关联很弱，a路径可能不显著]")
  cat("\n")
  
  # 中介-结果关联
  cor_med_out <- cor(data[[mediator_var]], data[[outcome_var]], use = "complete.obs")
  cat("   中介-结果相关系数:", round(cor_med_out, 4))
  if(abs(cor_med_out) < 0.05) cat(" [关联很弱，b路径可能不显著]")
  cat("\n")
  
  # 4. 统计效能评估
  cat("\n4. 统计效能评估:\n")
  
  # 基于相关系数的效能计算
  if(require(pwr, quietly = TRUE)) {
    power_exp_out <- pwr.r.test(n = complete_cases, r = abs(cor_exp_out), sig.level = 0.05)$power
    power_exp_med <- pwr.r.test(n = complete_cases, r = abs(cor_exp_med), sig.level = 0.05)$power
    power_med_out <- pwr.r.test(n = complete_cases, r = abs(cor_med_out), sig.level = 0.05)$power
    
    cat("   检测暴露-结果关联的效能:", round(power_exp_out, 3))
    if(power_exp_out < 0.8) cat(" [效能不足]")
    cat("\n")
    
    cat("   检测暴露-中介关联的效能:", round(power_exp_med, 3))
    if(power_exp_med < 0.8) cat(" [效能不足]")
    cat("\n")
    
    cat("   检测中介-结果关联的效能:", round(power_med_out, 3))
    if(power_med_out < 0.8) cat(" [效能不足]")
    cat("\n")
  }
  
  # 5. 快速中介效应检验
  cat("\n5. 快速中介效应检验:\n")
  
  # 准备数据
  analysis_data <- data[complete.cases(data[c(exposure_var, mediator_var, outcome_var)]), ]
  
  if(!is.null(covariates)) {
    covariate_complete <- complete.cases(analysis_data[covariates])
    analysis_data <- analysis_data[covariate_complete, ]
    cat("   调整协变量后样本量:", nrow(analysis_data), "\n")
  }
  
  # Baron & Kenny步骤
  if(nrow(analysis_data) > 50) {
    
    # 构建公式
    if(is.null(covariates)) {
      formula_total <- as.formula(paste(outcome_var, "~", exposure_var))
      formula_a <- as.formula(paste(mediator_var, "~", exposure_var))
      formula_direct <- as.formula(paste(outcome_var, "~", exposure_var, "+", mediator_var))
    } else {
      cov_string <- paste(covariates, collapse = " + ")
      formula_total <- as.formula(paste(outcome_var, "~", exposure_var, "+", cov_string))
      formula_a <- as.formula(paste(mediator_var, "~", exposure_var, "+", cov_string))
      formula_direct <- as.formula(paste(outcome_var, "~", exposure_var, "+", mediator_var, "+", cov_string))
    }
    
    # 拟合模型
    if(is.null(weights)) {
      model_total <- lm(formula_total, data = analysis_data)
      model_a <- lm(formula_a, data = analysis_data)
      model_direct <- lm(formula_direct, data = analysis_data)
    } else {
      model_total <- lm(formula_total, data = analysis_data, weights = analysis_data[[weights]])
      model_a <- lm(formula_a, data = analysis_data, weights = analysis_data[[weights]])
      model_direct <- lm(formula_direct, data = analysis_data, weights = analysis_data[[weights]])
    }
    
    # 提取系数和p值
    total_coef <- summary(model_total)$coefficients[exposure_var, "Estimate"]
    total_p <- summary(model_total)$coefficients[exposure_var, "Pr(>|t|)"]
    
    a_coef <- summary(model_a)$coefficients[exposure_var, "Estimate"]
    a_p <- summary(model_a)$coefficients[exposure_var, "Pr(>|t|)"]
    
    direct_coef <- summary(model_direct)$coefficients[exposure_var, "Estimate"]
    direct_p <- summary(model_direct)$coefficients[exposure_var, "Pr(>|t|)"]
    
    b_coef <- summary(model_direct)$coefficients[mediator_var, "Estimate"]
    b_p <- summary(model_direct)$coefficients[mediator_var, "Pr(>|t|)"]
    
    # 间接效应
    indirect_effect <- a_coef * b_coef
    
    cat("   总效应 (c路径):", round(total_coef, 4), ", p =", round(total_p, 4))
    if(total_p > 0.05) cat(" [不显著]")
    cat("\n")
    
    cat("   a路径 (暴露→中介):", round(a_coef, 4), ", p =", round(a_p, 4))
    if(a_p > 0.05) cat(" [不显著]")
    cat("\n")
    
    cat("   b路径 (中介→结果):", round(b_coef, 4), ", p =", round(b_p, 4))
    if(b_p > 0.05) cat(" [不显著]")
    cat("\n")
    
    cat("   直接效应 (c'路径):", round(direct_coef, 4), ", p =", round(direct_p, 4))
    if(direct_p > 0.05) cat(" [不显著]")
    cat("\n")
    
    cat("   间接效应 (a×b):", round(indirect_effect, 4), "\n")
    
    # 中介效应判断
    cat("\n   中介效应判断:\n")
    if(total_p > 0.05) {
      cat("   - 总效应不显著，可能无中介效应或存在抑制效应\n")
    }
    if(a_p > 0.05) {
      cat("   - a路径不显著，中介效应不成立\n")
    }
    if(b_p > 0.05) {
      cat("   - b路径不显著，中介效应不成立\n")
    }
    if(a_p <= 0.05 && b_p <= 0.05) {
      if(direct_p > 0.05) {
        cat("   - 可能存在完全中介效应\n")
      } else {
        cat("   - 可能存在部分中介效应\n")
      }
    }
  }
  
  # 6. 问题诊断和建议
  cat("\n6. 问题诊断和改进建议:\n")
  
  issues <- c()
  suggestions <- c()
  
  # 样本量问题
  if(complete_cases < 500) {
    issues <- c(issues, "样本量偏小")
    suggestions <- c(suggestions, "考虑合并多个NHANES周期数据")
  }
  
  # 缺失值问题
  if(missing_exposure/nrow(data) > 0.2 || missing_mediator/nrow(data) > 0.2 || missing_outcome/nrow(data) > 0.2) {
    issues <- c(issues, "缺失值比例过高")
    suggestions <- c(suggestions, "使用多重插补处理缺失值")
  }
  
  # 分布问题
  if(exists("exp_skew") && abs(exp_skew) > 2) {
    issues <- c(issues, "暴露变量严重偏态")
    suggestions <- c(suggestions, "对暴露变量进行对数转换或Box-Cox转换")
  }
  
  if(exists("med_skew") && abs(med_skew) > 2) {
    issues <- c(issues, "中介变量严重偏态")
    suggestions <- c(suggestions, "对中介变量进行对数转换")
  }
  
  # 关联强度问题
  if(abs(cor_exp_out) < 0.1) {
    issues <- c(issues, "暴露与结果关联很弱")
    suggestions <- c(suggestions, "检查暴露变量的测量和分类方法")
  }
  
  if(abs(cor_exp_med) < 0.1) {
    issues <- c(issues, "暴露与中介关联很弱")
    suggestions <- c(suggestions, "考虑其他中介变量或改进暴露测量")
  }
  
  if(abs(cor_med_out) < 0.1) {
    issues <- c(issues, "中介与结果关联很弱")
    suggestions <- c(suggestions, "检查中介变量的选择和测量")
  }
  
  # 统计效能问题
  if(exists("power_exp_out") && power_exp_out < 0.8) {
    issues <- c(issues, "统计效能不足")
    suggestions <- c(suggestions, "增加样本量或寻找更强的效应")
  }
  
  if(length(issues) > 0) {
    cat("   发现的问题:\n")
    for(i in 1:length(issues)) {
      cat("   -", issues[i], "\n")
    }
    
    cat("\n   改进建议:\n")
    for(i in 1:length(suggestions)) {
      cat("   -", suggestions[i], "\n")
    }
  } else {
    cat("   未发现明显问题，数据质量良好\n")
  }
  
  cat("\n=== 诊断完成 ===\n")
  
  # 返回诊断结果
  return(list(
    sample_size = complete_cases,
    correlations = list(
      exposure_outcome = cor_exp_out,
      exposure_mediator = cor_exp_med,
      mediator_outcome = cor_med_out
    ),
    issues = issues,
    suggestions = suggestions
  ))
}

# =============================================================================
# 快速修复函数
# =============================================================================

quick_fix_data <- function(data, exposure_var, mediator_var, outcome_var) {
  
  cat("正在应用快速修复...\n")
  
  fixed_data <- data
  
  # 1. 处理极端异常值（使用Winsorizing）
  winsorize <- function(x, probs = c(0.01, 0.99)) {
    quantiles <- quantile(x, probs, na.rm = TRUE)
    x[x < quantiles[1]] <- quantiles[1]
    x[x > quantiles[2]] <- quantiles[2]
    return(x)
  }
  
  fixed_data[[exposure_var]] <- winsorize(fixed_data[[exposure_var]])
  fixed_data[[mediator_var]] <- winsorize(fixed_data[[mediator_var]])
  fixed_data[[outcome_var]] <- winsorize(fixed_data[[outcome_var]])
  
  # 2. 对数转换（如果变量为正且偏态）
  if(all(fixed_data[[exposure_var]] > 0, na.rm = TRUE)) {
    skew_exp <- moments::skewness(fixed_data[[exposure_var]], na.rm = TRUE)
    if(abs(skew_exp) > 1) {
      fixed_data[[paste0(exposure_var, "_log")]] <- log(fixed_data[[exposure_var]])
      cat("已创建", exposure_var, "的对数转换变量\n")
    }
  }
  
  if(all(fixed_data[[mediator_var]] > 0, na.rm = TRUE)) {
    skew_med <- moments::skewness(fixed_data[[mediator_var]], na.rm = TRUE)
    if(abs(skew_med) > 1) {
      fixed_data[[paste0(mediator_var, "_log")]] <- log(fixed_data[[mediator_var]])
      cat("已创建", mediator_var, "的对数转换变量\n")
    }
  }
  
  # 3. 标准化连续变量
  numeric_vars <- sapply(fixed_data, is.numeric)
  fixed_data[numeric_vars] <- lapply(fixed_data[numeric_vars], function(x) {
    if(length(unique(x[!is.na(x)])) > 10) {  # 只标准化连续变量
      scale(x)[,1]
    } else {
      x
    }
  })
  
  cat("数据修复完成\n")
  return(fixed_data)
}

# =============================================================================
# 使用示例
# =============================================================================

# 示例用法：
# 
# # 1. 诊断分析
# diagnosis <- diagnose_mediation_analysis(
#   data = your_nhanes_data,
#   exposure_var = "log_arsenic",
#   mediator_var = "log_crp", 
#   outcome_var = "anxiety_score",
#   covariates = c("age", "gender", "race", "education", "bmi"),
#   weights = "WTMEC2YR"
# )
# 
# # 2. 快速修复（如果需要）
# fixed_data <- quick_fix_data(
#   data = your_nhanes_data,
#   exposure_var = "arsenic_total",
#   mediator_var = "crp",
#   outcome_var = "anxiety_score"
# )
# 
# # 3. 重新诊断
# diagnosis_after_fix <- diagnose_mediation_analysis(
#   data = fixed_data,
#   exposure_var = "arsenic_total_log",
#   mediator_var = "crp_log",
#   outcome_var = "anxiety_score"
# )

cat("诊断工具加载完成。使用 diagnose_mediation_analysis() 函数开始诊断。\n")