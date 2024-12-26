
if(!file.exists(str_c(getwd(),"/input"))){
  dir.create(str_c(getwd(),"/input"))
}
if(!file.exists(str_c(getwd(),"/output"))){
  dir.create(str_c(getwd(),"/output"))
}
#define path
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")


# 1-Impute demograpy and clinical vars with NA----------
tmp %>%
  dplyr::select(年龄_2Cat:HCG日_LH) %>%
  mice::mice(m = 5,
             seed = 1,
             method = "cart",
             print = F) %>%
  complete() %>%
  as_tibble() -> tmp2
tmp %<>%
  dplyr::select(-all_of(colnames(tmp2))) %>%
  cbind(tmp2) %>%
  as_tibble()
rm(tmp2)
# 2- Replace exposure vars of zero with LOD/sqrt(2)-------------

# change exposure values from 0 to NA to judge missing proportions.
tmp %>%
  naniar::miss_var_summary()
colnames(tmp)
# install.packages("labelled")
# library(labelled)
tmp %>% 
  mutate(across(.cols = Ag_H1,
                .fns = ~labelled::remove_var_label(.x))) %>%
  mutate(across(.cols = As_1:Zn_H1,
                .fns = ~ na_if(.x, 0))) -> tmp
tmp %>%
  naniar::miss_var_summary()
# metal and PFAS variables to be imputed
tmp %>%
  dplyr::select(As_1:Zn_H1) %>%
  colnames()%>%
  as.matrix() %>% as.vector() -> vars
# check: delete missing that >25%
tmp %>%
  naniar::miss_var_summary() %>%
  dplyr::filter(pct_miss>25) %>%
  dplyr::select(variable) %>%
  as.matrix() %>% as.vector() -> del.lst
# delete
tmp %<>%
  dplyr::select(-all_of(del.lst))
# update list
vars = setdiff(vars, del.lst)

# replace NA as 0 for metals and PFASs
tmp %>%
  as.matrix() %>% as_tibble() %>%
  mutate(
    across(.cols = As_1:Zn_H1,
           .fns = ~ replace_na(.x, 0))) -> tmp

# extract LOD for metals and PFAS
readxl::read_xlsx(str_c(getwd(), "/input/LOD_Total.xlsx"))  -> LODs
LODs %>%
  dplyr::filter(Varname  %in% vars) %>%
  arrange(Varname)  -> lst1

# Function: when values is lower than LOD，replace by LOD/sqrt(2)
LODreplace = function(var){
  df %>%
    dplyr::select(any_of(var)) %>%
    dplyr::rename(v1 = var) %>%
    dplyr::mutate(lod = as.vector(as.matrix(lst[which(lst$Varname==var),2]/sqrt(2)))) %>%
    dplyr::mutate(v1 = if_else(v1 < lod, #if
                                lod/sqrt(2), #yes
                                v1)) %>%
    dplyr::select(v1) %>%
    set_names(var) -> temp
  return(temp)
}
#execute
lst = lst1
df = tmp
purrr::map_dfc(lst1$Varname, LODreplace) -> temp

# complete imputation
tmp %<>%
  dplyr::select(-any_of(colnames(temp))) %>%
  cbind(temp) %>%
  as_tibble()

colnames(tmp)


# 3-Transer var type to factor-------------------
tmp %<>%
  mutate(across(.cols = any_of(c("X1临床妊娠_Y0_N1")),
                .fns = as.factor))
colnames(tmp)

# 1-Demograpical var summary (TableOne)----
install.packages("tableone")
library(tableone)
# Total vars need to be summarized
ttvars = c(C1, C2)
# Categorical var list
catvars <- c("年龄_2Cat", "BMI_3Cat", "民族_2Cat",
             "文化程度_4Cat", "职业_2Cat", #基本信息
             "主动吸烟","被动吸烟", #习惯
             "病因_3Cat", #不孕原因相关
             "治疗方案_3Cat", "受精方式",#治疗方案相关
             "X1移植胚胎数量", "X1新鲜冷冻"#移植相关
)
setdiff(ttvars, catvars)
tab1 <-  print(CreateTableOne(vars = ttvars, strata = "X1临床妊娠_Y0_N1",
                              data = tmp, factorVars = catvars,
                              includeNA = FALSE),
               exact = c("民族_2Cat",
                         "职业_2Cat",
                         "主动吸烟",
                         "病因_3Cat",
                         "X1移植胚胎数量",
                         "治疗方案_3Cat",
                         "文化程度_4Cat"),
               nonnormal = setdiff(ttvars, catvars),
               showAllLevels = TRUE)

rownames(tab1)

tab1 %>%
  as_tibble() %>%
  dplyr::mutate(var = rownames(tab1)) %>%
  set_names("level", "OnPregnancy", "EarlyPregnancyLoss", "P", "TestType", "Var") %>%
  dplyr::select(Var, level, OnPregnancy, EarlyPregnancyLoss, P, TestType) -> tab1
writexl::write_xlsx(tab1, str_c(getwd(), "/output/0-TableOne-results.xlsx"))


colnames(tmp)
# 4-Log transformation
tmp %<>%
  mutate(across(.cols = Ag_1:Zn_H1,
                .fns = log))
tmp %<>%
  mutate(across(.cols = 基础E2:甲功TSH,
                .fns = log))
tmp %<>%
  mutate(across(.cols = HCG日_E2:HCG日_LH,
                .fns = log))

# 5-Scale
# define scale parameters
RangeUpper = 1
RangeLower = 0
Direct = "positive" #direction, either "positive" or "negative"
# Function for scale
FuncScaleRange = function(Var, df){
  df %>%
    dplyr::select(all_of(Var)) %>%
    .[[Var]] %>%
    as_tibble() -> temp

  range(temp$value, na.rm = TRUE) -> rng

  switch(Direct,
         "positive" = {
           temp  %>%
             dplyr::mutate(V2 = (RangeUpper-RangeLower)*(value-rng[1])/
                             (rng[2]-rng[1])+RangeLower) %>%
             dplyr::select(V2) %>%
             purrr::set_names(Var) -> temp1
         },

         "negative" = {
           temp  %>%
             dplyr::mutate(value = (RangeUpper-RangeLower)*(rng[2]-value)/
                             (rng[2]-rng[1])+RangeLower) %>%
             purrr::set_names(Var) -> temp1
         }
  )
  return(temp1)
}
# extract var names that to be scaled
tmp %>%
  dplyr::select(Ag_1:Zn_H1) %>%
  names() -> Vars
# scale execution: 不区分训练集/测试集
purrr::map_dfc(Vars, FuncScaleRange, df = tmp) %>%
  cbind(tmp %>%
          dplyr::select(-all_of(Vars))) %>%
  as_tibble() %>%
  dplyr::select(names(tmp)) -> tmp
tmp

# check
tmp %>%
  naniar::miss_var_summary()
colnames(tmp)
# arrange vars
tmp %<>%
  dplyr::select(any_of(VarY),
                "地区",
                any_of(C1),
                any_of(C2),
                everything())

# #export
tmp %>%
  writexl::write_xlsx(str_c(getwd(),"/input","/1-IVF-n161.xlsx"))

# save
save.image(str_c(getwd(),"/input/eSet.rdata"))





