# NHANES砷暴露与焦虑关系的中介效应分析
# 作者: 研究助手
# 日期: 2025-10-02

# 加载必要的包
library(tidyverse)
library(survey)
library(mediation)
library(lavaan)
library(NHANES)
library(psych)
library(VIM)
library(mice)
library(corrplot)

# =============================================================================
# 1. 数据准备和清洗
# =============================================================================

# 加载NHANES数据（这里以2015-2016和2017-2018为例）
# 您需要根据实际情况调整年份和变量

load_nhanes_data <- function() {
  # 人口统计学数据
  demo <- nhanes("DEMO_I") # 2015-2016
  demo2 <- nhanes("DEMO_J") # 2017-2018
  
  # 砷暴露数据（尿液重金属）
  arsenic <- nhanes("UHM_I") # 2015-2016
  arsenic2 <- nhanes("UHM_J") # 2017-2018
  
  # 心理健康数据（焦虑相关）
  mental <- nhanes("DPQ_I") # 抑郁筛查问卷，包含焦虑相关项目
  mental2 <- nhanes("DPQ_J")
  
  # 炎症标志物
  inflam <- nhanes("LAB25_I") # C反应蛋白等
  inflam2 <- nhanes("LAB25_J")
  
  # 氧化应激标志物（如果有的话）
  # 注意：NHANES可能没有直接的氧化应激标志物，可能需要用其他指标代替
  
  # 合并数据
  data_2015 <- demo %>%
    left_join(arsenic, by = "SEQN") %>%
    left_join(mental, by = "SEQN") %>%
    left_join(inflam, by = "SEQN") %>%
    mutate(cycle = "2015-2016")
  
  data_2017 <- demo2 %>%
    left_join(arsenic2, by = "SEQN") %>%
    left_join(mental2, by = "SEQN") %>%
    left_join(inflam2, by = "SEQN") %>%
    mutate(cycle = "2017-2018")
  
  # 合并两个周期的数据
  combined_data <- bind_rows(data_2015, data_2017)
  
  return(combined_data)
}

# 数据清洗函数
clean_nhanes_data <- function(data) {
  cleaned_data <- data %>%
    # 选择成年人（18岁以上）
    filter(RIDAGEYR >= 18) %>%
    
    # 创建砷暴露变量（总砷或无机砷）
    mutate(
      # 尿砷浓度（μg/L）
      arsenic_total = URXUAS,  # 总砷
      arsenic_inorganic = URXUAB + URXUAC,  # 无机砷（砷酸+亚砷酸）
      
      # 肌酐校正的砷浓度
      arsenic_creatinine_adj = arsenic_total / URXUCR * 100,
      
      # 对数转换（处理右偏分布）
      log_arsenic = log(arsenic_total + 1),
      log_arsenic_adj = log(arsenic_creatinine_adj + 1)
    ) %>%
    
    # 创建焦虑相关变量
    mutate(
      # 从抑郁问卷中提取焦虑相关项目
      anxiety_score = rowSums(select(., starts_with("DPQ")), na.rm = TRUE),
      
      # 二分类焦虑变量
      anxiety_binary = ifelse(anxiety_score >= 10, 1, 0)  # 阈值可调整
    ) %>%
    
    # 创建炎症变量
    mutate(
      # C反应蛋白（mg/dL）
      crp = LBXCRP,
      log_crp = log(crp + 1),
      
      # 白细胞计数
      wbc = LBXWBCSI,
      
      # 高炎症状态（CRP > 3 mg/L）
      high_inflammation = ifelse(crp > 0.3, 1, 0)  # 转换为mg/dL
    ) %>%
    
    # 创建协变量
    mutate(
      # 年龄组
      age_group = case_when(
        RIDAGEYR < 30 ~ "18-29",
        RIDAGEYR < 50 ~ "30-49",
        RIDAGEYR < 65 ~ "50-64",
        TRUE ~ "65+"
      ),
      
      # 性别
      gender = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
      
      # 种族/民族
      race = factor(RIDRETH3, levels = c(1, 2, 3, 4, 6, 7),
                   labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
                             "Non-Hispanic Black", "Non-Hispanic Asian", "Other")),
      
      # 教育水平
      education = factor(DMDEDUC2, levels = c(1, 2, 3, 4, 5),
                        labels = c("Less than 9th grade", "9-11th grade", 
                                  "High school", "Some college", "College graduate")),
      
      # 家庭收入
      income = INDHHIN2,
      
      # BMI
      bmi = BMXBMI,
      bmi_category = case_when(
        bmi < 18.5 ~ "Underweight",
        bmi < 25 ~ "Normal",
        bmi < 30 ~ "Overweight",
        TRUE ~ "Obese"
      ),
      
      # 吸烟状态
      smoking = factor(SMQ020, levels = c(1, 2), labels = c("Yes", "No"))
    ) %>%
    
    # 移除缺失值过多的观测
    filter(
      !is.na(arsenic_total),
      !is.na(anxiety_score),
      !is.na(crp),
      !is.na(RIDAGEYR),
      !is.na(gender)
    )
  
  return(cleaned_data)
}

# =============================================================================
# 2. 探索性数据分析
# =============================================================================

exploratory_analysis <- function(data) {
  cat("=== 探索性数据分析 ===\n")
  
  # 基本描述统计
  cat("\n1. 样本基本信息:\n")
  cat("总样本量:", nrow(data), "\n")
  
  # 主要变量的描述统计
  key_vars <- c("arsenic_total", "anxiety_score", "crp", "RIDAGEYR")
  desc_stats <- data %>%
    select(all_of(key_vars)) %>%
    describe()
  print(desc_stats)
  
  # 检查变量分布
  cat("\n2. 变量分布检查:\n")
  
  # 砷暴露分布
  p1 <- ggplot(data, aes(x = arsenic_total)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7) +
    scale_x_log10() +
    labs(title = "砷暴露分布（对数尺度）", x = "总砷浓度 (μg/L)", y = "频数")
  
  # 焦虑评分分布
  p2 <- ggplot(data, aes(x = anxiety_score)) +
    geom_histogram(bins = 30, fill = "lightcoral", alpha = 0.7) +
    labs(title = "焦虑评分分布", x = "焦虑评分", y = "频数")
  
  # CRP分布
  p3 <- ggplot(data, aes(x = crp)) +
    geom_histogram(bins = 50, fill = "lightgreen", alpha = 0.7) +
    scale_x_log10() +
    labs(title = "C反应蛋白分布（对数尺度）", x = "CRP (mg/dL)", y = "频数")
  
  # 保存图形
  ggsave("/workspace/arsenic_distribution.png", p1, width = 8, height = 6)
  ggsave("/workspace/anxiety_distribution.png", p2, width = 8, height = 6)
  ggsave("/workspace/crp_distribution.png", p3, width = 8, height = 6)
  
  # 相关性分析
  cat("\n3. 主要变量相关性:\n")
  cor_vars <- data %>%
    select(log_arsenic, anxiety_score, log_crp, RIDAGEYR, bmi) %>%
    na.omit()
  
  cor_matrix <- cor(cor_vars)
  print(cor_matrix)
  
  # 保存相关性图
  png("/workspace/correlation_matrix.png", width = 800, height = 600)
  corrplot(cor_matrix, method = "color", type = "upper", 
           addCoef.col = "black", tl.cex = 0.8, number.cex = 0.7)
  dev.off()
  
  return(list(desc_stats = desc_stats, cor_matrix = cor_matrix))
}

# =============================================================================
# 3. 中介效应分析
# =============================================================================

# Baron & Kenny方法
baron_kenny_mediation <- function(data) {
  cat("=== Baron & Kenny中介效应分析 ===\n")
  
  # 步骤1: X -> Y (总效应)
  model1 <- lm(anxiety_score ~ log_arsenic + RIDAGEYR + gender + race + 
               education + bmi + smoking, data = data, weights = WTMEC2YR)
  
  cat("\n步骤1 - 总效应 (砷 -> 焦虑):\n")
  print(summary(model1)$coefficients["log_arsenic", ])
  
  # 步骤2: X -> M (a路径)
  model2 <- lm(log_crp ~ log_arsenic + RIDAGEYR + gender + race + 
               education + bmi + smoking, data = data, weights = WTMEC2YR)
  
  cat("\n步骤2 - a路径 (砷 -> 炎症):\n")
  print(summary(model2)$coefficients["log_arsenic", ])
  
  # 步骤3: X + M -> Y (b路径和直接效应)
  model3 <- lm(anxiety_score ~ log_arsenic + log_crp + RIDAGEYR + gender + race + 
               education + bmi + smoking, data = data, weights = WTMEC2YR)
  
  cat("\n步骤3 - 控制中介变量后:\n")
  cat("b路径 (炎症 -> 焦虑):\n")
  print(summary(model3)$coefficients["log_crp", ])
  cat("直接效应 (砷 -> 焦虑):\n")
  print(summary(model3)$coefficients["log_arsenic", ])
  
  # 计算间接效应
  a_coef <- coef(model2)["log_arsenic"]
  b_coef <- coef(model3)["log_crp"]
  indirect_effect <- a_coef * b_coef
  
  cat("\n间接效应 (a × b):", indirect_effect, "\n")
  
  return(list(
    total_effect = model1,
    a_path = model2,
    direct_effect = model3,
    indirect_effect = indirect_effect
  ))
}

# 使用mediation包进行bootstrap检验
bootstrap_mediation <- function(data) {
  cat("=== Bootstrap中介效应分析 ===\n")
  
  # 准备数据（移除缺失值）
  analysis_data <- data %>%
    select(anxiety_score, log_arsenic, log_crp, RIDAGEYR, gender, race,
           education, bmi, smoking, WTMEC2YR) %>%
    na.omit()
  
  cat("分析样本量:", nrow(analysis_data), "\n")
  
  # 中介变量模型
  mediator_model <- lm(log_crp ~ log_arsenic + RIDAGEYR + gender + race + 
                      education + bmi + smoking, 
                      data = analysis_data, weights = WTMEC2YR)
  
  # 结果变量模型
  outcome_model <- lm(anxiety_score ~ log_arsenic + log_crp + RIDAGEYR + gender + race + 
                     education + bmi + smoking, 
                     data = analysis_data, weights = WTMEC2YR)
  
  # Bootstrap中介效应检验
  set.seed(123)
  mediation_results <- mediate(mediator_model, outcome_model, 
                              treat = "log_arsenic", mediator = "log_crp",
                              boot = TRUE, sims = 1000)
  
  cat("\nBootstrap中介效应结果:\n")
  print(summary(mediation_results))
  
  return(mediation_results)
}

# 结构方程模型方法
sem_mediation <- function(data) {
  cat("=== 结构方程模型中介效应分析 ===\n")
  
  # 准备数据
  analysis_data <- data %>%
    select(anxiety_score, log_arsenic, log_crp, RIDAGEYR, gender, race,
           education, bmi, smoking) %>%
    na.omit() %>%
    mutate(
      gender_num = as.numeric(gender) - 1,
      race_white = ifelse(race == "Non-Hispanic White", 1, 0),
      education_num = as.numeric(education)
    )
  
  # 定义SEM模型
  sem_model <- '
    # 中介路径
    log_crp ~ a * log_arsenic + RIDAGEYR + gender_num + race_white + education_num + bmi
    anxiety_score ~ b * log_crp + c * log_arsenic + RIDAGEYR + gender_num + race_white + education_num + bmi
    
    # 间接效应
    indirect := a * b
    
    # 总效应
    total := c + (a * b)
  '
  
  # 拟合模型
  fit <- sem(sem_model, data = analysis_data, se = "bootstrap", bootstrap = 1000)
  
  cat("\nSEM中介效应结果:\n")
  print(summary(fit, standardized = TRUE))
  
  # 参数估计
  params <- parameterEstimates(fit, boot.ci.type = "bca.simple")
  print(params[params$label %in% c("a", "b", "c", "indirect", "total"), ])
  
  return(fit)
}

# =============================================================================
# 4. 敏感性分析和诊断
# =============================================================================

sensitivity_analysis <- function(data) {
  cat("=== 敏感性分析 ===\n")
  
  # 1. 不同砷暴露指标的比较
  cat("\n1. 不同砷暴露指标比较:\n")
  
  # 使用肌酐校正的砷浓度
  model_adj <- lm(anxiety_score ~ log_arsenic_adj + log_crp + RIDAGEYR + gender + race + 
                  education + bmi + smoking, data = data, weights = WTMEC2YR)
  
  cat("肌酐校正砷浓度结果:\n")
  print(summary(model_adj)$coefficients[c("log_arsenic_adj", "log_crp"), ])
  
  # 2. 不同炎症指标的比较
  if("wbc" %in% names(data)) {
    model_wbc <- lm(anxiety_score ~ log_arsenic + log(wbc + 1) + RIDAGEYR + gender + race + 
                    education + bmi + smoking, data = data, weights = WTMEC2YR)
    
    cat("\n使用白细胞计数作为炎症指标:\n")
    print(summary(model_wbc)$coefficients[c("log_arsenic", "log(wbc + 1)"), ])
  }
  
  # 3. 分层分析
  cat("\n2. 分层分析:\n")
  
  # 按性别分层
  male_data <- filter(data, gender == "Male")
  female_data <- filter(data, gender == "Female")
  
  if(nrow(male_data) > 100 & nrow(female_data) > 100) {
    male_med <- bootstrap_mediation(male_data)
    female_med <- bootstrap_mediation(female_data)
    
    cat("男性中介效应:\n")
    print(summary(male_med))
    cat("女性中介效应:\n")
    print(summary(female_med))
  }
  
  # 4. 异常值检验
  cat("\n3. 异常值检验:\n")
  
  # 识别砷暴露的异常值
  arsenic_q99 <- quantile(data$arsenic_total, 0.99, na.rm = TRUE)
  data_no_outliers <- filter(data, arsenic_total <= arsenic_q99)
  
  cat("移除砷暴露异常值后的样本量:", nrow(data_no_outliers), "\n")
  
  # 重新分析
  if(nrow(data_no_outliers) > 500) {
    outlier_results <- bootstrap_mediation(data_no_outliers)
    cat("移除异常值后的中介效应:\n")
    print(summary(outlier_results))
  }
}

# =============================================================================
# 5. 主函数
# =============================================================================

main_analysis <- function() {
  cat("开始NHANES砷暴露中介效应分析...\n")
  
  # 注意：由于这是示例代码，实际运行时需要确保有NHANES数据访问权限
  # 这里提供一个模拟数据生成函数用于测试
  
  # 如果有真实NHANES数据，取消下面的注释
  # data <- load_nhanes_data()
  # cleaned_data <- clean_nhanes_data(data)
  
  # 模拟数据用于演示
  set.seed(123)
  n <- 2000
  cleaned_data <- data.frame(
    SEQN = 1:n,
    arsenic_total = exp(rnorm(n, 2, 1)),
    anxiety_score = rnorm(n, 8, 4),
    crp = exp(rnorm(n, 0, 1)),
    RIDAGEYR = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("Male", "Female"), n, replace = TRUE)),
    race = factor(sample(c("Non-Hispanic White", "Non-Hispanic Black", "Mexican American"), n, replace = TRUE)),
    education = factor(sample(1:5, n, replace = TRUE)),
    bmi = rnorm(n, 28, 5),
    smoking = factor(sample(c("Yes", "No"), n, replace = TRUE)),
    WTMEC2YR = runif(n, 1000, 50000)
  ) %>%
    mutate(
      log_arsenic = log(arsenic_total),
      log_crp = log(crp),
      # 添加一些真实的关联
      anxiety_score = anxiety_score + 0.1 * log_arsenic + 0.05 * log_crp + rnorm(n, 0, 1)
    )
  
  # 探索性分析
  eda_results <- exploratory_analysis(cleaned_data)
  
  # 中介效应分析
  baron_kenny_results <- baron_kenny_mediation(cleaned_data)
  bootstrap_results <- bootstrap_mediation(cleaned_data)
  sem_results <- sem_mediation(cleaned_data)
  
  # 敏感性分析
  sensitivity_analysis(cleaned_data)
  
  # 保存结果
  save(baron_kenny_results, bootstrap_results, sem_results, eda_results,
       file = "/workspace/mediation_analysis_results.RData")
  
  cat("\n分析完成！结果已保存到 mediation_analysis_results.RData\n")
  
  return(list(
    data = cleaned_data,
    eda = eda_results,
    baron_kenny = baron_kenny_results,
    bootstrap = bootstrap_results,
    sem = sem_results
  ))
}

# 运行分析
# results <- main_analysis()