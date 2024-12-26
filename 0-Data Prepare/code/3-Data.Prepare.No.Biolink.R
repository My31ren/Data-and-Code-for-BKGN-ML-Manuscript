rm(list=ls())
# install packages
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

# create folders
if(!file.exists(str_c(getwd(),"/input"))){
  dir.create(str_c(getwd(),"/input"))
}
if(!file.exists(str_c(getwd(),"/output"))){
  dir.create(str_c(getwd(),"/output"))
}
# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")
# read data
tmp = readxl::read_xlsx(str_c(inpath, "/1-IVF-n116.xlsx"))
colnames(tmp)

# arrange vars
C1 = c("年龄_2Cat", "BMI_3Cat", "民族_2Cat","文化程度_4Cat","职业_2Cat", #基本信息
       "主动吸烟","被动吸烟",# 生活习惯
       "病因_3Cat")# 病因
# Define C  临床诊疗信息
C2 = c("治疗方案_3Cat", "受精方式","X1移植胚胎数量", "X1新鲜冷冻","Gn量", "Gn天数", # 治疗方案相关
       "X1受精率", "X1优质胚胎率", "X1卵裂囊胚", "IVF卵子数", "M2成熟卵子", "受精卵数2PN",# 胚胎质量相关
       "基础E2", "基础LH", "基础FSH", "基础T睾酮","甲功FT4","甲功FT3","甲功TSH", # 激素相关
       "卵泡数", "HCG日_内膜","HCG日_E2", "HCG日_LH") # 取卵日相关

tmp %>% 
    dplyr::select(matches("^\\w\\w?_1$")) %>% 
    colnames() %>% str_sort() -> metal_1
tmp %>% 
    dplyr::select(matches("^\\w\\w\\w+_1$")) %>% 
    colnames() %>% str_sort() -> pfas_1
tmp %>% 
    dplyr::select(matches("^\\w\\w?_2$")) %>% 
    colnames() %>% str_sort() -> metal_2
tmp %>% 
    dplyr::select(matches("^\\w\\w\\w+_2$")) %>% 
    colnames() %>% str_sort() -> pfas_2
tmp %>% 
    dplyr::select(matches("^\\w\\w?_F$")) %>% 
    colnames() %>% str_sort() -> metal_f
tmp %>% 
    dplyr::select(matches("^\\w\\w?_H1$")) %>% 
    colnames() %>% str_sort() -> metal_h

# arrange
tmp %<>%
  dplyr::select(any_of("X1临床妊娠_Y0_N1"),
                地区,
                any_of(C1),
                any_of(C2),
                any_of(c(metal_1, pfas_1)), 
                any_of(c(metal_2, pfas_2)),
                any_of(metal_f), 
                any_of(metal_h),
                everything())
colnames(tmp)
# 6-export------
tmp %>%
  writexl::write_xlsx(str_c(outpath,"/1-No.Biolink/eSetdata.xlsx"))



