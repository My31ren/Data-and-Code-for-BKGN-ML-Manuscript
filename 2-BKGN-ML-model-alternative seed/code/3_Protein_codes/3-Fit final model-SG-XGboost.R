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

#-----1. Prepare procedures-----------
#load data
work.dir = str_c(getwd(), "/@-Knowledge.Driven.model-SG model")
inpath = str_c(str_c(work.dir,"/input/3-Protein.Addition"))
outpath = str_c(str_c(work.dir,"/output/3-Protein.Addition"))
load(str_c(outpath, "/eSet_Feature_Select.EXPs.RData"))
eSet$FinalFeatures

work.dir = str_c(getwd(), "/@-Knowledge.Driven.model-SG model")
# define path: No-Biolink
if(!file.exists(str_c(work.dir,"/input/3-Protein.Addition"))){
  dir.create(str_c(work.dir,"/input/3-Protein.Addition"))
}
if(!file.exists(str_c(work.dir,"/output/3-Protein.Addition"))){
  dir.create(str_c(work.dir,"/output/3-Protein.Addition"))
}
inpath = str_c(str_c(work.dir,"/input/3-Protein.Addition"))
outpath = str_c(str_c(work.dir,"/output/3-Protein.Addition"))

# creat folder
if(!file.exists(str_c(outpath,"/0_","Final_Model" ))){
  dir.create(str_c(outpath,"/0_","Final_Model"  ))
}
if(!file.exists(str_c(outpath,"/0_","Final_Model","/Features_Explain"))){
  dir.create(str_c(outpath,"/0_","Final_Model","/Features_Explain"))
}
if(!file.exists(str_c(outpath,"/0_","Final_Model","/Features_Selection"))){
  dir.create(str_c(outpath,"/0_","Final_Model" ,"/Features_Selection"))
}
if(!file.exists(str_c(outpath,"/0_","Final_Model" ,"/Model_Summary"))){
  dir.create(str_c(outpath,"/0_","Final_Model" ,"/Model_Summary"))
}

# define data: Y C X
# Outcome
VarY = "Y"
# n folds for cv resampling
FoldRsmp = 5
AutoTuneN = 50
VarsImpThr = 0.85

##2-SG model ##########################################################

# 2.1 Fit base models ----------------------------
Base.model = function(i, Data, Voca){
  tictoc::tic()

  ##########
  # i = 4
  # eSet$Data = eSet$S2$Data
  # eSet$Voca = eSet$S2$Voca
  #########

  eSet$Data = Data
  eSet$Voca = Voca
  # 2.1.1 Extract exposures-----------------
  eSet$FinalFeatures[[i]] %>%
    dplyr::select(Features) %>%
    as.matrix() %>% as.vector -> VarsX
  eSet$Data %>%
    dplyr::select(VarY, all_of(VarsX)) -> df.train
  # 2.1.2 set train task-----------------
  df.train %>%
    as_task_classif(target = VarY) -> task.base
  # prepare resampling--------------------------
  set.seed(123)#keep rsmp results as the same
  rsmp = rsmp("cv", folds = FoldRsmp)
  rsmp$instantiate(task.base)
  # 2.1.3 set learners----------------------
  # -xgboost(needs to be tuned)-
  learner_xg = lrn(str_c("classif",".xgboost"),
                   predict_sets = c("train", "test"),
                   predict_type = "prob")
  learner_xg$param_set$values$eta = 0.1
  # define autotune
  at = AutoTuner$new(
    learner = learner_xg,
    resampling = rsmp("cv",folds = FoldRsmp),
    measure = msr("classif.acc"),
    search_space = ps(
      subsample = p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = p_dbl(lower = 0.5, upper = 1),
      max_depth = p_int(lower = 3, upper = 10),
      nrounds = p_int(lower = 2, upper = 100)
    ),
    terminator = trm("evals", n_evals = AutoTuneN),
    tuner = tnr("random_search")
  )
  #autotune
  ddpcr::quiet(
    at$train(task.base)
  )
  # use tuned parameters
  learner_xg$param_set$values <- at$model$tuning_instance$result_learner_param_vals
  # 2.1.4 benchmarking: to obtain predict values of Y in each cv test datasets---------------
  # design
  design = benchmark_grid(
    tasks = task.base,
    learners = learner_xg,
    resamplings = rsmp)
  # execute benchmarking
  ddpcr::quiet(
    benchmark(design) -> bmr
  )
  # 2.1.5 extract base model prediction results------------------------
  map(1:5, 
    ~ as.data.table(bmr, 
      reassemble_learners = T, 
      convert_predictions = T,
      predict_sets = "test")$prediction[[.x]]
      ) -> tmp1
  ddpcr::quiet(map_dfr(1:5, 
    ~ bind_cols(tmp1[[.x]]$row_ids, tmp1[[.x]]$truth, tmp1[[.x]]$response, tmp1[[.x]]$prob) %>% 
      set_names(c("row_ids", "truth", "response", "prob.0", "prob.1"))
      ) %>% 
    arrange(row_ids) -> eSet$BaseModel[[list.names(eSet$FinalFeatures)[[i]]]])
  # 2.1.6 extract base model performance and feature importance--------
  measures_bmr = list(msr("classif.acc", predict_sets = "train", id = str_c("acc","_train")),
                        msr("classif.acc", id =str_c("acc","_test")),
                        msr("classif.precision", predict_sets = "train", id = str_c("precision","_train")),
                        msr("classif.precision", id =str_c("precision","_test")),
                        msr("classif.sensitivity", predict_sets = "train", id = str_c("sensitivity","_train")),
                        msr("classif.sensitivity", id =str_c("sensitivity","_test")),
                        msr("classif.specificity", predict_sets = "train", id = str_c("specificity","_train")),
                        msr("classif.specificity", id =str_c("specificity","_test")),
                        msr("classif.npv", predict_sets = "train", id = str_c("npv","_train")),
                        msr("classif.npv", id =str_c("npv","_test"))
                        )
  bmr$score(measures_bmr) %>% 
      as_tibble() %>% 
      group_by(nr)  %>% # different variable combination
      dplyr::mutate(rsmp_method = "cv-5") %>% 
      dplyr::select(nr, learner_id, rsmp_method, acc_train:npv_test) %>% 
      # for some cases, some measures may get NA values (not error, dont worry) since the results is not suitable (such as npv and precision)
      dplyr::mutate(across(.cols =acc_train:npv_test,  ~ mean(.x, na.rm = T)))  %>% 
      unique()  %>% 
      writexl::write_xlsx(str_c(outpath,"/0_Final_Model/Model_Summary/Base.Model-",list.names(eSet$FinalFeatures)[[i]], "-Feature.Importance.xlsx"))

  learner_xg$train(task.base)$importance() %>% 
    sort(decreasing = T)  -> imp_xg
  imp_xg  %>% 
    names() %>% tibble(Features = .) %>%
    mutate(importance = if_else(imp_xg<=0,0,imp_xg)) %>%
    mutate(importance_perct = map_dbl(importance,
                                      ~(./sum(importance))),
           importance_perct_cum = cumsum(importance_perct)) %>%
    left_join(eSet$Voca %>%
                dplyr::select(VarName, VarLabel, OmicGroup) %>%
                rename(Features = "VarName"),
              by = "Features") %>%
    dplyr::select(Features, VarLabel, OmicGroup, everything())  %>% 
    writexl::write_xlsx(str_c(outpath,"/0_Final_Model/Features_Selection/Base.Model-",list.names(eSet$FinalFeatures)[[i]], "Feature.Importance.xlsx"))
  #----------------------------------------------------------------------
  print(str_c("Base model for ", list.names(eSet$FinalFeatures)[[i]] , " group has been established !"))
  tictoc::toc()

  return(eSet)

}
# Execute
purrr::pmap(
       list(as.list(1:6), 
            list(eSet$S1$Data, eSet$S1$Data, eSet$S1$Data, eSet$S2$Data, eSet$FF$Data, eSet$H1$Data), 
            list(eSet$S1$Voca, eSet$S1$Voca, eSet$S1$Voca, eSet$S2$Voca, eSet$FF$Voca, eSet$H1$Voca)
       ),Base.model)


# 2.2----Build SG model-------------------------------------
# 2.2.1 Feature selection for SG model----------

# Extract exposures from base models------------
eSet$BaseModel %>%
  map_dfc( ~ .x  %>% dplyr::select(prob.1)) %>% 
  set_names(str_c("pred.1.", list.names(eSet$BaseModel))) -> base.features

# Extract exposures: for covariates--------------
eSet$FinalFeatures[["Baseline"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX.1
eSet$FinalFeatures[["Clinical.Factors"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX.2
VarsX = c(VarsX.1, VarsX.2) # var names
eSet$S1$Data %>%
  dplyr::select(all_of(VarsX)) -> tmp1 # vars
eSet$S1$Voca %>% 
  dplyr::filter(VarName %in% VarsX) %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> names
eSet$S1$Voca %>% 
  dplyr::filter(VarName %in% colnames(tmp1)) %>% 
  dplyr::select(VarName) %>% 
  as.matrix() %>% as.vector() -> vars
tmp1 %>% 
  dplyr::select(all_of(vars)) -> tmp1 # arrange
colnames(tmp1) = names
# Extract exposures: for serum1-----------------------
eSet$FinalFeatures[["Serum1"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX # var names
eSet$S1$Data %>%
  dplyr::select(VarsX) -> tmp2 # vars
eSet$S1$Voca %>% 
  dplyr::filter(VarName %in% VarsX) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXP1" ~ str_c(VarLabel, ".s1"),
                                     OmicGroup != "EXP1" ~ VarLabel)) %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> names
eSet$S1$Voca %>% 
  dplyr::filter(VarName %in% colnames(tmp2)) %>% 
  dplyr::select(VarName) %>% 
  as.matrix() %>% as.vector() -> vars
tmp2 %>% 
  dplyr::select(any_of(vars)) -> tmp2 #arrange
colnames(tmp2) = names
# Extract exposures: for serum2------------------------------
eSet$FinalFeatures[["Serum2"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX
eSet$S2$Data %>%
  dplyr::select(VarsX) -> tmp3
eSet$S2$Voca %>% 
  dplyr::filter(VarName %in% VarsX) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXP2" ~ str_c(VarLabel, ".s2"),
                                     OmicGroup != "EXP2" ~ VarLabel)) %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> names
eSet$S2$Voca %>% 
  dplyr::filter(VarName %in% colnames(tmp3)) %>% 
  dplyr::select(VarName) %>% 
  as.matrix() %>% as.vector() -> vars
tmp3 %>% 
  dplyr::select(any_of(vars)) -> tmp3#arrange
colnames(tmp3) = names
# Extract exposures: for Hair-----------------------------
eSet$FinalFeatures[["Hair"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX
eSet$H1$Data %>%
  dplyr::select(VarsX) -> tmp4
eSet$H1$Voca %>% 
  dplyr::filter(VarName %in% VarsX) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXPhr" ~ str_c(VarLabel, ".hr"),
                                     OmicGroup != "EXPhr" ~ VarLabel)) %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> names
eSet$H1$Voca %>% 
  dplyr::filter(VarName %in% colnames(tmp4)) %>% 
  dplyr::select(VarName) %>% 
  as.matrix() %>% as.vector() -> vars
tmp4 %>% 
  dplyr::select(any_of(vars)) -> tmp4#Arrange
colnames(tmp4) = names
# Extract exposures: for FF-------------------------
eSet$FinalFeatures[["FF"]] %>%
  as_tibble() %>% 
  dplyr::select(Features) %>%
  as.matrix() %>% as.vector() -> VarsX
eSet$FF$Data %>%
  dplyr::select(VarsX) -> tmp5
eSet$FF$Voca %>% 
  dplyr::filter(VarName %in% VarsX) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXPff" ~ str_c(VarLabel, ".ff"),
                                     OmicGroup != "EXPff" ~ VarLabel)) %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> names
eSet$FF$Voca %>% 
  dplyr::filter(VarName %in% colnames(tmp5)) %>% 
  dplyr::select(VarName) %>% 
  as.matrix() %>% as.vector() -> vars
tmp5 %>% 
  dplyr::select(any_of(vars))  -> tmp5 #ARRANGE
colnames(tmp5) = names
# Combine exposures and Y------------
cbind(eSet$S1$Data %>%
  dplyr::select(Y), tmp1, tmp2, tmp3, tmp4, tmp5) %>% 
  as_tibble() -> eSet$FN$xg$Data
colnames(eSet$FN$xg$Data)

# Create Vocabulary dataset for this new dataset------------
# create voca database
eSet$FN$xg$Voca = 
  tibble(VarName = c("Y", str_c("X", (1:(ncol(eSet$FN$xg$Data)-1))))
          ,VarLabel = colnames(eSet$FN$xg$Data))
# combine omicgroup information
eSet$FN$xg$Voca %>% 
  left_join(eSet$S1$Voca %>% # left join covariates and serum1 information
  dplyr::select(VarLabel, OmicGroup)  %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXP1" ~ str_c(VarLabel, ".s1"),
                                     OmicGroup != "EXP1" ~ VarLabel))
  , by = "VarLabel")  %>% 
  dplyr::mutate(OmicGroup = case_when(
                  VarLabel == "Y" ~ "Outcome",
                  .default = OmicGroup)) %>% 
  left_join(eSet$S2$Voca %>% 
  dplyr::select(VarLabel, OmicGroup) %>%  
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXP2" ~ str_c(VarLabel, ".s2"),
                                     OmicGroup != "EXP2" ~ VarLabel))
  , by = "VarLabel") %>%
  left_join(eSet$H1$Voca %>% 
  dplyr::select(VarLabel, OmicGroup) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXPhr" ~ str_c(VarLabel, ".hr"),
                                     OmicGroup != "EXPhr" ~ VarLabel))
  , by = "VarLabel") %>% 
  left_join(eSet$FF$Voca %>% 
  dplyr::select(VarLabel, OmicGroup) %>% 
  dplyr::mutate(VarLabel = case_when(OmicGroup == "EXPff" ~ str_c(VarLabel, ".ff"),
                                     OmicGroup != "EXPff" ~ VarLabel))
  , by = "VarLabel") -> test
test  %>% 
  print(n=nrow(.))
# tidy
test  %>% 
  dplyr::mutate(across(.cols = everything(), ~str_replace_na(.x, ""))) %>% 
  dplyr::mutate(OmicGroup = case_when(
    OmicGroup.x %in% c("Outcome", "Baseline","Clinical.Factors") ~ OmicGroup.x,
    !OmicGroup.x %in% c("Outcome", "Baseline","Clinical.Factors") ~ str_c(OmicGroup.x,  OmicGroup.y,  OmicGroup.x.x,  OmicGroup.y.y)
   )) %>% 
  dplyr::select(VarName, VarLabel, OmicGroup)  -> eSet$FN$xg$Voca
eSet$FN$xg$Voca %>% 
  print(n=nrow(.))

# Rename variable names as in vocabulary file------------
colnames(eSet$FN$xg$Data)  = eSet$FN$xg$Voca$VarName

# Define train dataset---------------
eSet$FN$xg$Data %>%
  bind_cols(base.features) -> df.train
# Set train task----------------------------------------- 
df.train %>%
  as_task_classif(target = VarY) -> task0.train
# prepare resampling--------------------------
set.seed(123)#keep rsmp results as the same
rsmp = rsmp("cv", folds = FoldRsmp)
rsmp$instantiate(task0.train)

# 2.1.6 Set learners--------------------------------------------
# -xgboost(needs to be tuned)-
learner_xg = lrn(str_c("classif",".xgboost"),
                 predict_sets = c("train", "test"),
                 predict_type = "prob")
learner_xg$param_set$values$eta = 0.1
# define autotune
at = AutoTuner$new(
  learner = learner_xg,
  resampling = rsmp("cv",folds = FoldRsmp),
  measure = msr("classif.acc"),
  search_space = ps(
    subsample = p_dbl(lower = 0.5, upper = 1),
    colsample_bytree = p_dbl(lower = 0.5, upper = 1),
    max_depth = p_int(lower = 3, upper = 10),
    nrounds = p_int(lower = 2, upper = 100)
  ),
  terminator = trm("evals", n_evals = AutoTuneN),
  tuner = tnr("random_search")
)
#autotune
ddpcr::quiet(
  at$train(task0.train)
)
# use tuned parameters
learner_xg$param_set$values <- at$model$tuning_instance$result_learner_param_vals

# importance rank-------------------------------------------
learner_xg$train(task0.train)
# rank importance
learner_xg$importance() %>% sort(decreasing=T) %>%
  as.matrix() %>%
  rownames() %>%
  as_tibble() %>%
  cbind(learner_xg$importance() %>% sort(decreasing=T) %>%
          as.matrix() %>%
          as_tibble()) %>%
  set_names(c("Features", "Importance")) -> importance_xg
# extract rank and generate tasks list----------------------------
tasklist  = map(c(1:nrow(importance_xg)),
             ~ df.train %>%
               dplyr::select(all_of(VarY),
                             all_of(importance_xg$Features[1:.x])) %>%
               as_task_classif(target = VarY)
               )
names(tasklist) = str_c("n_", 1:nrow(importance_xg))

# benchmarking-------------------------------------
# design
design = benchmark_grid(
  tasks = tasklist,
  learners = learner_xg,
  resamplings = rsmp)
# execute benchmarking
ddpcr::quiet(
  benchmark(design) -> bmr
)
# extract benchmarking results-------------------------
measures_bmr = list(msr("classif.acc", predict_sets = "train", id = str_c("acc","_train")),
                    msr("classif.acc", id =str_c("acc","_test")),
                    msr("classif.precision", predict_sets = "train", id = str_c("precision","_train")),
                    msr("classif.precision", id =str_c("precision","_test")),
                    msr("classif.sensitivity", predict_sets = "train", id = str_c("sensitivity","_train")),
                    msr("classif.sensitivity", id =str_c("sensitivity","_test")),
                    msr("classif.specificity", predict_sets = "train", id = str_c("specificity","_train")),
                    msr("classif.specificity", id =str_c("specificity","_test")),
                    msr("classif.npv", predict_sets = "train", id = str_c("npv","_train")),
                    msr("classif.npv", id =str_c("npv","_test"))
                    )
bmr$score(measures_bmr) %>% 
   as_tibble() %>% 
   group_by(nr)  %>% # different variable combination
   dplyr::mutate(rsmp_method = "cv-5") %>% 
   dplyr::select(nr, learner_id, rsmp_method, acc_train:npv_test) %>% 
   # for some cases, some measures may get NA values (not error, dont worry) since the results is not suitable (such as npv and precision)
   dplyr::mutate(across(.cols =acc_train:npv_test,  ~ mean(.x, na.rm = T)))  %>% 
   unique() %>% 
   dplyr::mutate(nr_list = str_c("n_", nr)) %>%
   cbind(unlist(map(1:nrow(importance_xg), ~paste0(importance_xg$Features[1:.x], collapse = ",")))) %>% 
   dplyr::rename(Feature_list = "...15") -> eSet$FeatureSelection[["FinalModel"]]
print(as.data.frame(eSet$FeatureSelection[["FinalModel"]]))

# export benchmarking results------------------------
eSet$FeatureSelection[["FinalModel"]] %>%
  as_tibble() %>% 
  writexl::write_xlsx(str_c(outpath,"/0_","Final_Model","/Model_Summary/Model.Stat-Features.Combinations-SG.Model.xlsx"))
  
# plot-----------------------------------------------
eSet$FeatureSelection[["FinalModel"]] %>%
    as_tibble() %>%
    dplyr::filter(row.names(eSet$FeatureSelection[["FinalModel"]])=="1") %>%
    dplyr::mutate(nr_list = "n_0", 
                  nr = as.integer(0),
                  across(.cols = acc_train:npv_test, ~.x*0)
                 ) %>%
    rbind(eSet$FeatureSelection[["FinalModel"]]) %>%
    dplyr::select(nr, nr_list, contains("_test"))  %>% 
    dplyr::mutate(sum = (acc_test + precision_test + sensitivity_test + specificity_test + npv_test)/5) %>% 
    set_names(c("nr", "nr_list", "ACC", "PPV", "TPR", "TNR", "NPV", "Avg")) %>% 
    tidyr::pivot_longer(cols = ACC:Avg,
                        names_to = "msr", 
                        values_to = "value") %>% 
    dplyr::mutate(msr = factor(msr, level = c("ACC", "PPV", "TPR", "TNR", "NPV", "Avg")), 
                  group = factor(if_else(msr == "Avg", "Avg", "Others"), level = c("Others", "Avg")))  -> plotdata
# extract n with maximum of Avg value
plotdata %>% 
    dplyr::filter(msr == "Avg") %>% 
    dplyr::select(value) %>%
    as.matrix() %>% 
    which.max()  -> row.index
plotdata %>% 
    dplyr::filter(msr == "Avg") %>% 
    .[row.index,"nr"]  %>% as.matrix() -> x.index
plotdata %>% 
    dplyr::filter(msr == "Avg") %>% 
    .[row.index,"value"]  %>% as.matrix() -> y.index

plotdata  %>% 
    ggplot(aes(x = nr,
             y = value,
             group = msr, 
             color = msr)) +
    geom_point()+
    geom_line()+
    scale_x_continuous(
      expand = c(0.01, 0.01), 
      breaks = seq(0, nrow(eSet$FeatureSelection[["FinalModel"]]), by = 1),
      labels = str_c("n_", seq(0, nrow(eSet$FeatureSelection[["FinalModel"]]), by = 1))) +
    scale_y_continuous(
      expand = c(0.01, 0.01), 
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
    ) + 
    geom_segment(aes(x = x.index, y = 0, xend = x.index, yend = y.index),
                 color = "red", 
                 linewidth = 1,
                 linetype="dashed") +
    geom_segment(aes(x = 0, y = y.index, xend = x.index, yend = y.index),
                 color = "red", 
                 linewidth = 1,
                 linetype="dashed") + 
    theme_bw()
ggsave(str_c(outpath,"/0_","Final_Model","/Features_Selection/Feature.Selection.png"), 
         height = 8, 
         width = 9)

# extract final selected features--------------------------------
importance_xg %>%
  dplyr::slice_head(n = x.index) %>%
  left_join(eSet$FN$xg$Voca %>%
          dplyr::select(VarName, VarLabel, OmicGroup) %>%
          rename(Features = "VarName"),
          by = "Features") %>% 
  dplyr::mutate(Importance = if_else(Importance < 0, 0, Importance)) -> eSet$FinalFeatures[["FinalModel"]]
eSet$FinalFeatures[["FinalModel"]] %>%
  writexl::write_xlsx(str_c(outpath,"/0_","Final_Model","/Features_Selection/Selected.Features-SG.Model.xlsx"))
print(str_c("Final features of ", "Final Model", " have been selected !"))
print(eSet$FinalFeatures[["FinalModel"]])

# 2.2.2 Fit final SG model------------------------------
# Extract exposures-----------------
df.train
colnames(df.train)
# 2.1.2 set train task-----------------
names(base.features)
df.train  %>% 
  dplyr::select(Y, any_of(eSet$FinalFeatures[["FinalModel"]]$Features)) %>% 
  as_task_classif(target = VarY) -> task.final_train

# # prepare resampling------------------------------
# set.seed(123)#keep rsmp results as the same
# rsmp = rsmp("cv", folds = FoldRsmp)
# rsmp$instantiate(task.final_train)
# # set learners----------------------
# # -xgboost(needs to be tuned)-
# learner_xg = lrn(str_c("classif",".xgboost"),
#                  predict_sets = c("train", "test"),
#                  predict_type = "prob")
# learner_xg$param_set$values$eta = 0.1
# # define autotune
# at = AutoTuner$new(
#   learner = learner_xg,
#   resampling = rsmp("cv",folds = FoldRsmp),
#   measure = msr("classif.acc"),
#   search_space = ps(
#     subsample = p_dbl(lower = 0.5, upper = 1),
#     colsample_bytree = p_dbl(lower = 0.5, upper = 1),
#     max_depth = p_int(lower = 3, upper = 10),
#     nrounds = p_int(lower = 2, upper = 100)
#   ),
#   terminator = trm("evals", n_evals = AutoTuneN),
#   tuner = tnr("random_search")
# )
# #autotune
# ddpcr::quiet(
#   at$train(task.final_train)
# )
# # use tuned parameters
# learner_xg$param_set$values <- at$model$tuning_instance$result_learner_param_vals

# Benchmarking for establish final model-------------------
# design benchmarking------------------------
design = benchmark_grid(
  tasks = c(task.final_train),
  learners = c(learner_xg),
  resamplings = rsmp)
# execute benchmarking-----------------------
ddpcr::quiet(
  benchmark(design, store_models = T) -> bmr
)
# extract benchmarking results-------------------------
measures_bmr = list(msr("classif.acc", predict_sets = "train", id = str_c("acc","_train")),
                    msr("classif.acc", id =str_c("acc","_test")),
                    msr("classif.precision", predict_sets = "train", id = str_c("precision","_train")),
                    msr("classif.precision", id =str_c("precision","_test")),
                    msr("classif.sensitivity", predict_sets = "train", id = str_c("sensitivity","_train")),
                    msr("classif.sensitivity", id =str_c("sensitivity","_test")),
                    msr("classif.specificity", predict_sets = "train", id = str_c("specificity","_train")),
                    msr("classif.specificity", id =str_c("specificity","_test")),
                    msr("classif.npv", predict_sets = "train", id = str_c("npv","_train")),
                    msr("classif.npv", id =str_c("npv","_test"))
                    )
bmr$score(measures_bmr) %>% 
    as_tibble() %>% 
    group_by(learner_id)  %>% 
    dplyr::mutate(rsmp_method = "cv-5") %>% 
    dplyr::select(learner_id, rsmp_method, acc_train:npv_test) %>% 
    # for some cases, some measures may get NA values (not error, dont worry) since the results is not suitable (such as npv and precision)
    dplyr::mutate(across(.cols =acc_train:npv_test,  ~ mean(.x, na.rm = T)))  %>% 
    unique() -> eSet$ModelStat[["FinalModel.SG"]]

eSet$ModelStat[["FinalModel.SG"]] %>% as.matrix()

# export final model accuracy-----------------------
eSet$ModelStat[["FinalModel.SG"]] %>%
    writexl::write_xlsx(str_c(outpath,"/0_Final_Model/Model_Summary/Measures.Stat-SG.Model.xlsx"))

# 2.2.3 Train final model------------------------------
learner_xg$train(task.final_train)

learner_xg$importance()%>%
  sort(decreasing=T) -> imp_xg
imp_xg %>%
  names() %>% tibble(Features = .) %>%
  mutate(importance = if_else(imp_xg<=0,0,imp_xg)) %>%

  mutate(importance_perct = map_dbl(importance,
                                    ~(./sum(importance))),
         importance_perct_cum = cumsum(importance_perct)) %>%
  left_join(eSet$FN$xg$Voca %>%
              dplyr::select(VarName, VarLabel, OmicGroup) %>%
              rename(Features = "VarName"),
            by = "Features") %>%
  dplyr::select(Features, VarLabel, OmicGroup, everything()) %>%
  writexl::write_xlsx(str_c(outpath,"/0_Final_Model/Features_Selection/importance.xg-SG.Model.xlsx"))

# # 2.2 Explain using DALEX---------------------------------------------
# eSet$FN$xg$Data$Y = as.numeric(eSet$FN$xg$Data$Y)
# set.seed(222)
# ranger_exp = explain_mlr3(learner_xg,
#                           data = eSet$FN$xg$Data[,-1],
#                           y = eSet$FN$xg$Data$Y,
#                           label = "XGBoost",
#                           colorize = F)
# # partial dependence1: threshold 0.95--------------------------------
# importance_xg = readxl::read_xlsx(str_c(outpath,"/0_Final_Model/Features_Selection/importance.xg.xlsx"))
# importance_xg %>% 
#   dplyr::filter(importance_perct_cum <=0.95) %>% 
#   dplyr::select(Features) %>% 
#   as.matrix() %>% as.vector() -> varlist
# eSet$FN$xg$Voca %>% 
#   dplyr::filter(VarName %in% varlist)%>% 
#   writexl::write_xlsx(str_c(outpath,"/0_","Final_Model","/Features_Explain/Selected_features.XGBoost(Threshold0.95).xlsx"))
# # plot-----------------
# rf_pd = model_profile(ranger_exp,
#                       variables = varlist)$agr_profiles
# rf_pd  %>% 
#   as_tibble() %>% 
#   left_join(eSet$FN$xg$Voca %>% dplyr::select(VarLabel, VarName), by = c("_vname_"="VarName"))  %>% 
#   set_names(c("VarName", "_label_", "x", "yhat", "_id_", "VarLabel")) %>% 
#   dplyr::select(VarName, VarLabel, x, yhat) -> plotdata
# #reorder
# importance_xg %>% 
#   dplyr::filter(importance_perct_cum <=0.95)  %>% 
#   dplyr::select(VarLabel)  %>% as.matrix() %>% as.vector -> a1
# plotdata$VarLabel = factor(plotdata$VarLabel, levels =a1)
# #plot
# plotdata %>% 
#   ggplot(aes(x = x, y = yhat)) +
#   geom_smooth(se = T) +
#   facet_wrap(~VarLabel, ncol = 5, scales = "free") + 
#   theme_bw()+
#   theme(strip.text.x = element_text(size = 8),
#   axis.text.x = element_text(size = 9),
#   axis.text.y = element_text(size = 9))
# #save
# ggsave(str_c(outpath,"/0_","Final_Model","/Features_Explain/PD.XGBoost.jpg"),
#       width = 8, height = 16)

##4-Save eSet-------------------------------------------------------------
save.image(str_c(outpath, "/eSet_FitFinalModel-SG.Model.RData"))


load(str_c(outpath, "/eSet_FitFinalModel-SG.Model.RData"))
eSet$ModelStat[["FinalModel.SG"]]%>% 
  as.matrix()
