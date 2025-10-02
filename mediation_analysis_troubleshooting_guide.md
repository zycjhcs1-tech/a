# NHANES砷暴露中介效应分析问题诊断与解决方案

## 问题概述

您在使用NHANES数据库研究砷暴露导致焦虑，炎症或氧化应激作为中介变量时，发现没有显著的中介效应或P值很大。这是一个常见问题，本文档将帮助您系统性地诊断和解决这些问题。

## 常见问题及解决方案

### 1. 样本量和统计效能问题

#### 问题表现：
- P值普遍较大（>0.05）
- 置信区间很宽
- 效应量很小

#### 可能原因：
- 样本量不足
- 效应量本身较小
- 变量测量误差大

#### 解决方案：

```r
# 1.1 计算所需样本量
library(pwr)

# 基于预期效应量计算样本量
# 假设预期相关系数为0.1-0.2（小到中等效应）
pwr.r.test(r = 0.15, sig.level = 0.05, power = 0.8)

# 1.2 合并多个NHANES周期数据
# 使用2013-2018年的数据（3个周期）
cycles <- c("2013-2014", "2015-2016", "2017-2018")
combined_data <- map_dfr(cycles, function(cycle) {
  # 加载每个周期的数据
  # ... 数据加载代码
})

# 1.3 使用适当的抽样权重
# NHANES是复杂抽样设计，必须使用权重
library(survey)

# 创建survey设计对象
nhanes_design <- svydesign(
  ids = ~SDMVPSU,        # 主抽样单元
  strata = ~SDMVSTRA,    # 分层变量
  weights = ~WTMEC2YR,   # 抽样权重
  nest = TRUE,
  data = your_data
)

# 使用survey包进行回归分析
model <- svyglm(anxiety_score ~ log_arsenic + covariates, 
                design = nhanes_design)
```

### 2. 变量选择和测量问题

#### 问题表现：
- 暴露变量与结果变量无关联
- 中介变量与暴露或结果变量无关联

#### 可能原因：
- 砷暴露指标选择不当
- 焦虑测量不准确
- 炎症/氧化应激指标不合适

#### 解决方案：

```r
# 2.1 优化砷暴露变量
# 使用多种砷暴露指标
data <- data %>%
  mutate(
    # 总砷
    total_arsenic = URXUAS,
    
    # 无机砷（更有毒性）
    inorganic_arsenic = URXUAB + URXUAC,  # 砷酸 + 亚砷酸
    
    # 肌酐校正（控制尿液稀释）
    arsenic_creatinine = total_arsenic / URXUCR * 100,
    
    # 比重校正
    arsenic_sg_adj = total_arsenic * (1.024 - 1) / (URXUSG - 1),
    
    # 对数转换（处理偏态分布）
    log_arsenic = log(inorganic_arsenic + 1)
  )

# 2.2 改进焦虑测量
# 使用多个焦虑相关指标
anxiety_vars <- c("DPQ030", "DPQ040", "DPQ070", "DPQ080")  # 焦虑相关题目
data$anxiety_composite <- rowMeans(data[anxiety_vars], na.rm = TRUE)

# 或使用GAD-7量表（如果可用）
# 或使用其他心理健康指标作为替代

# 2.3 选择合适的炎症指标
inflammation_indicators <- c(
  "LBXCRP",    # C反应蛋白（最常用）
  "LBXWBCSI",  # 白细胞计数
  "LBXLYPCT",  # 淋巴细胞百分比
  "LBXNEPCT"   # 中性粒细胞百分比
)

# 创建炎症综合指标
data <- data %>%
  mutate(
    crp_log = log(LBXCRP + 1),
    wbc_log = log(LBXWBCSI + 1),
    # 标准化后合并
    inflammation_composite = scale(crp_log)[,1] + scale(wbc_log)[,1]
  )
```

### 3. 混杂因素控制不充分

#### 问题表现：
- 加入协变量后效应消失
- 不同模型结果差异很大

#### 解决方案：

```r
# 3.1 识别重要的混杂因素
confounders <- c(
  # 人口学特征
  "RIDAGEYR",      # 年龄
  "RIAGENDR",      # 性别
  "RIDRETH3",      # 种族/民族
  
  # 社会经济地位
  "DMDEDUC2",      # 教育水平
  "INDHHIN2",      # 家庭收入
  
  # 生活方式
  "SMQ020",        # 吸烟史
  "ALQ101",        # 饮酒
  "PAQ605",        # 体力活动
  
  # 健康状况
  "BMXBMI",        # BMI
  "BPQ020",        # 高血压史
  "DIQ010",        # 糖尿病史
  
  # 饮食因素
  "DR1TKCAL",      # 总热量摄入
  "DR1TSFAT",      # 饱和脂肪
  
  # 其他暴露
  "LBXBPB",        # 血铅
  "LBXBCD"         # 血镉
)

# 3.2 使用有向无环图(DAG)指导变量选择
library(dagitty)
library(ggdag)

# 构建理论模型
dag <- dagitty("dag {
  Arsenic -> Anxiety
  Arsenic -> Inflammation -> Anxiety
  Age -> {Arsenic Inflammation Anxiety}
  Gender -> {Arsenic Inflammation Anxiety}
  SES -> {Arsenic Inflammation Anxiety}
  Smoking -> {Arsenic Inflammation Anxiety}
  BMI -> {Inflammation Anxiety}
}")

# 识别最小充分调整集
adjustmentSets(dag, exposure = "Arsenic", outcome = "Anxiety")

# 3.3 使用倾向性评分方法
library(MatchIt)

# 创建高/低砷暴露组
data$high_arsenic <- ifelse(data$log_arsenic > median(data$log_arsenic, na.rm = TRUE), 1, 0)

# 倾向性评分匹配
match_obj <- matchit(high_arsenic ~ RIDAGEYR + RIAGENDR + RIDRETH3 + DMDEDUC2 + 
                     BMXBMI + SMQ020, data = data, method = "nearest")

matched_data <- match.data(match_obj)
```

### 4. 统计方法问题

#### 问题表现：
- 不同方法结果不一致
- 置信区间包含0

#### 解决方案：

```r
# 4.1 使用多种中介效应分析方法

# Baron & Kenny方法（传统）
baron_kenny <- function(data) {
  # 步骤1-3的分析
  # ...
}

# Bootstrap方法（推荐）
library(mediation)
bootstrap_mediation <- function(data, n_boot = 5000) {
  # 增加bootstrap次数提高精度
  mediate(mediator_model, outcome_model, 
          treat = "arsenic", mediator = "inflammation",
          boot = TRUE, sims = n_boot)
}

# 结构方程模型
library(lavaan)
sem_model <- '
  # 测量模型（如果有潜变量）
  Inflammation =~ crp + wbc + neutrophil_pct
  
  # 结构模型
  Inflammation ~ a * Arsenic + covariates
  Anxiety ~ b * Inflammation + c * Arsenic + covariates
  
  # 间接效应
  indirect := a * b
  total := c + indirect
'

# 4.2 使用稳健标准误
library(sandwich)
library(lmtest)

# 异方差稳健标准误
coeftest(model, vcov = vcovHC(model, type = "HC3"))

# 聚类稳健标准误（如果有聚类结构）
coeftest(model, vcov = vcovCL(model, cluster = data$cluster_var))

# 4.3 贝叶斯方法
library(brms)

bayesian_mediation <- brm(
  bf(anxiety ~ arsenic + inflammation + covariates) +
  bf(inflammation ~ arsenic + covariates),
  data = data,
  chains = 4,
  iter = 4000,
  cores = 4
)
```

### 5. 数据质量问题

#### 解决方案：

```r
# 5.1 处理缺失数据
library(mice)
library(VIM)

# 检查缺失模式
aggr(data[key_variables], col = c('navyblue','red'), 
     numbers = TRUE, sortVars = TRUE)

# 多重插补
mice_imputed <- mice(data[key_variables], m = 5, method = 'pmm')
completed_data <- complete(mice_imputed)

# 5.2 异常值处理
# 使用稳健的异常值检测
library(robustbase)

# Mahalanobis距离
mahal_dist <- mahalanobis(data[numeric_vars], 
                         colMeans(data[numeric_vars], na.rm = TRUE),
                         cov(data[numeric_vars], use = "complete.obs"))

# 标记异常值
data$outlier <- mahal_dist > qchisq(0.999, df = length(numeric_vars))

# 5.3 数据转换
# 处理偏态分布
data <- data %>%
  mutate(
    # Box-Cox转换
    arsenic_bc = forecast::BoxCox(arsenic_total, lambda = 0.3),
    
    # Yeo-Johnson转换（可处理负值）
    crp_yj = car::yjPower(crp, lambda = 0.2),
    
    # 分位数转换
    arsenic_quantile = qnorm(rank(arsenic_total) / (length(arsenic_total) + 1))
  )
```

## 具体改进建议

### 1. 增强统计效能

```r
# 效应量计算
library(effectsize)

# 计算Cohen's d
cohens_d(anxiety_score ~ high_arsenic_group, data = data)

# 计算eta squared
eta_squared(aov(anxiety_score ~ arsenic_quartile, data = data))

# 功效分析
library(pwr)
pwr.f2.test(u = 5, v = 1000, f2 = 0.02, sig.level = 0.05)  # u=预测变量数-1, v=样本量-预测变量数-1
```

### 2. 模型诊断

```r
# 2.1 线性假设检验
library(car)
crPlots(model)  # 成分-残差图

# 2.2 同方差性检验
bptest(model)   # Breusch-Pagan检验
plot(model, which = 1)  # 残差vs拟合值图

# 2.3 正态性检验
shapiro.test(residuals(model))
qqPlot(model)

# 2.4 多重共线性检验
vif(model)  # 方差膨胀因子

# 2.5 影响点诊断
influenceIndexPlot(model)
```

### 3. 敏感性分析

```r
# 3.1 分层分析
# 按性别分层
male_results <- mediation_analysis(filter(data, gender == "Male"))
female_results <- mediation_analysis(filter(data, gender == "Female"))

# 按年龄分层
young_results <- mediation_analysis(filter(data, age < 50))
old_results <- mediation_analysis(filter(data, age >= 50))

# 3.2 不同阈值的敏感性
thresholds <- c(0.9, 0.95, 0.99)
for(thresh in thresholds) {
  outlier_cutoff <- quantile(data$arsenic_total, thresh, na.rm = TRUE)
  subset_data <- filter(data, arsenic_total <= outlier_cutoff)
  results <- mediation_analysis(subset_data)
  print(paste("Threshold:", thresh, "N:", nrow(subset_data)))
  print(results)
}

# 3.3 不同变量组合
variable_combinations <- list(
  c("total_arsenic", "crp"),
  c("inorganic_arsenic", "crp"),
  c("total_arsenic", "wbc"),
  c("arsenic_creatinine_adj", "inflammation_composite")
)

results_comparison <- map(variable_combinations, function(vars) {
  mediation_analysis_custom(data, exposure = vars[1], mediator = vars[2])
})
```

## 结果解释指南

### 1. 效应量解释

- **小效应**: r = 0.1, d = 0.2
- **中等效应**: r = 0.3, d = 0.5  
- **大效应**: r = 0.5, d = 0.8

### 2. 中介效应类型

- **完全中介**: 直接效应不显著，间接效应显著
- **部分中介**: 直接效应和间接效应都显著
- **无中介**: 间接效应不显著

### 3. 临床意义评估

即使统计学显著性较低，也要考虑：
- 效应量的实际意义
- 剂量-反应关系
- 生物学合理性
- 与既往研究的一致性

## 常见错误避免

1. **不要仅依赖p值**: 关注效应量和置信区间
2. **不要忽略复杂抽样设计**: 必须使用适当的权重和设计效应
3. **不要过度调整**: 避免调整中介路径上的变量
4. **不要忽略非线性关系**: 考虑使用样条函数或多项式项
5. **不要忽略交互作用**: 检验暴露与协变量的交互效应

## 报告建议

### 结果报告模板

```
在调整了年龄、性别、种族、教育水平、BMI和吸烟状态后，
砷暴露每增加1个对数单位，焦虑评分增加X分（95% CI: X-X, p = X）。

C反应蛋白在砷暴露与焦虑关系中起部分中介作用，
间接效应为X（95% CI: X-X），占总效应的X%。

敏感性分析显示，在排除异常值、使用不同砷暴露指标和
分层分析中，结果基本一致。
```

## 进一步研究建议

1. **考虑纵向设计**: 如果可能，使用多个时间点的数据
2. **基因-环境交互**: 考虑遗传易感性的调节作用
3. **多重中介模型**: 同时考虑炎症和氧化应激的中介作用
4. **非线性关系**: 使用限制性三次样条等方法
5. **机器学习方法**: 使用随机森林、LASSO等方法进行变量选择

记住：环境流行病学研究中的效应量通常较小，需要大样本量和精确的测量才能检测到显著效应。不显著的结果也可能反映真实的生物学关系。