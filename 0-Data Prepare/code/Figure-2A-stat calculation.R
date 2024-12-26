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
getwd()
#load data
node = readxl::read_xlsx(str_c(getwd(), "/output/Figure-2A.Nodes.for.Plotting.xlsx"))
corr = readxl::read_xlsx(str_c(getwd(), "/output/Figure-2A.Correlation.xlsx"))

#correlation within same biometric
corr %>% 
    dplyr::filter(p<0.05) %>% 
    dplyr::filter(
        (str_detect(.$group1, ".1$")&str_detect(.$group2, ".1$"))|
        (str_detect(.$group1, ".2$")&str_detect(.$group2, ".2$"))|
         group1 == group2)

# correlation summary
f1= function(grp1, grp2){
corr %>% 
    dplyr::filter(p<0.05) %>% 
    dplyr::filter(
        (group1 == grp1 & group2 == grp2)|
            (group1 == grp2 & group2 == grp1)) %>% 
    dplyr::mutate(n = nrow(.)) %>% 
    dplyr::select(group1, group2, n) %>% 
    unique()
}
combn(unique(corr$group1),2) %>% t() %>%
    as_tibble() %>% 
    set_names(c("grp1", "grp2")) -> lst

purrr::map2_dfr(lst$grp1, lst$grp2, f1) %>% 
    arrange(-n) -> stat

# correlation with EPL
node %>% 
    dplyr::filter(MinusLogP > -log10(0.05)) %>% 
    arrange(-MinusLogP)
