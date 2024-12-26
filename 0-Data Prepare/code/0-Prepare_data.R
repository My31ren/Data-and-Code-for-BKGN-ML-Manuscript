rm(list=ls())
#install packages
pacman::p_load(
  "bestNormalize",
  "broom",
  "car",
  "caret",
  "cluster",
  "DALEX",
  "DALEXtra",
  "data.table",
  "ddpcr",
  "factoextra",
  "fs",
  "future",
  "furrr",
  "gee",
  "geepack",
  "ggplot2",
  "ggfortify",
  "gridExtra",
  "GGally",
  "gtsummary",
  "ggpattern",
  "Hmisc",
  "lmtest",
  "lubridate",
  "magrittr",
  "mice",
  "mlr3verse",
  "mlr3viz",
  "mlr3extralearners",
  "naniar",
  "patchwork",
  "precrec",
  "rlist",
  "rstatix",
  "R6",
  "rstatix",
  "RCy3",
  "tidyverse",
  "tictoc",
  "vip",
  "vroom",
  "zip",
  "ggsci",
  "RColorBrewer",
  "circlize",
  "C50", #below machine learning
  "cluster",
  "coin",
  "dbarts",
  "e1071",
  "earth",
  "FNN",
  "gbm",
  "glmnet",
  "kernlab",
  "kknn",
  "LiblineaR",
  "lightgbm",
  "MASS",
  "mboost",
  "mgcv",
  "nnet",
  "partykit",
  "ranger",
  "rpart",
  "sandwich",
  "stats",
  "xgboost"
)

# 0-Extract study population -----------------------------------------------
#note: 金属测定时，Zn检出率为100%，因此以其作为“是否测样”的筛选变量，得到的数据集代表人群为“所有参与一血、二血、卵泡液以及头发样本金属元素测定的样品”。N = 121
load(str_c(getwd(), "/input/Article-1.Rdata"))
data %>%
  dplyr::filter(is.na(Zn_1)==F & is.na(Zn_2)==F & is.na(Zn_F)==F & is.na(Zn_H1)==F) -> df


# 1-Define Y X C ----------------------------------------------------------
# Define Y: 未临床妊娠(早期流产)
VarY = "X1临床妊娠_Y0_N1"

# Define C  基本信息
C1 = c("年龄_2Cat", "BMI_3Cat", "民族_2Cat","文化程度_4Cat","职业_2Cat", #基本信息
          "主动吸烟","被动吸烟",# 生活习惯
          "病因_3Cat")# 病因

# Define C  临床诊疗信息
C2 = c("治疗方案_3Cat", "受精方式","X1移植胚胎数量", "X1新鲜冷冻","Gn量", "Gn天数", # 治疗方案相关
  "X1受精率", "X1优质胚胎率", "X1卵裂囊胚", "IVF卵子数", "M2成熟卵子", "受精卵数2PN",# 胚胎质量相关
  "基础E2", "基础LH", "基础FSH", "基础T睾酮","甲功FT4","甲功FT3","甲功TSH", # 激素相关
  "卵泡数", "HCG日_内膜","HCG日_E2", "HCG日_LH"  # 取卵日相关
  )
# Define X  第一次血清金属
df %>%
  dplyr::select(As_1:Li_1) %>%
  colnames() -> X1
# Define X  第二次血清金属
df %>%
  dplyr::select(As_2:Li_2) %>%
  colnames() -> X2
# Define X  卵泡液金属
df %>%
  dplyr::select(As_F:Li_F) %>%
  colnames() -> X3
# Define X  第一次血清PFAS
df %>%
  dplyr::select(PFBA_1:P_62Cl_PFESA_1) %>%
  colnames() -> X4
# Define X  第二次血清PFAS
df %>%
  dplyr::select(PFBA_2:P_62Cl_PFESA_2) %>%
  colnames() -> X5
# Define X  头发
df %>%
  dplyr::select(contains("_H1")) %>%
  colnames() -> X6


# 2-Bind columns -----------------------------
df %>%
  dplyr::select(any_of(VarY),
                地区,# group variable to distinguish train and test set
                any_of(C1), # baseline vars
                any_of(C2), # clinical associated vars
                any_of(X1), # first serum metals
                any_of(X2), # second serum metals
                any_of(X3), # FF metals
                any_of(X4),
                any_of(X5),
                any_of(X6)
                ) -> tmp
tmp %>%
  dplyr::filter(!is.na(PFOA_1)) -> tmp

