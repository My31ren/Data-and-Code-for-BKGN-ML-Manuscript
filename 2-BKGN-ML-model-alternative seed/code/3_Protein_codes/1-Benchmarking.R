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

#-----1. Prepare procedures------------------------------
work.dir = str_c(getwd(), "/@-Knowledge.Driven.model-Add.Algorithms")
# define path: No-Biolink
if(!file.exists(str_c(work.dir,"/input/3-Protein.Addition"))){
  dir.create(str_c(work.dir,"/input/3-Protein.Addition"))
}
if(!file.exists(str_c(work.dir,"/output/3-Protein.Addition"))){
  dir.create(str_c(work.dir,"/output/3-Protein.Addition"))
}
inpath = str_c(str_c(work.dir,"/input/3-Protein.Addition"))
outpath = str_c(str_c(work.dir,"/output/3-Protein.Addition"))

# initialize R6 object
eSet <- R6::R6Class(
  "eSet",
  lock_class = FALSE,
  lock_objects = FALSE
)
# import data---------------------------
eSet$S1$Data <- readxl::read_xlsx(str_c(inpath,"/eSetdata.Serum#1.xlsx"))
eSet$S1$Voca <- readxl::read_xlsx(str_c(inpath,"/eSetvoca.Serum#1.xlsx"))

eSet$S2$Data <- readxl::read_xlsx(str_c(inpath,"/eSetdata.Serum#2.xlsx"))
eSet$S2$Voca <- readxl::read_xlsx(str_c(inpath,"/eSetvoca.Serum#2.xlsx"))

eSet$H1$Data <- readxl::read_xlsx(str_c(inpath,"/eSetdata.Hair#1.xlsx"))
eSet$H1$Voca <- readxl::read_xlsx(str_c(inpath,"/eSetvoca.Hair#1.xlsx"))

eSet$FF$Data <- readxl::read_xlsx(str_c(inpath,"/eSetdata.FF.xlsx"))
eSet$FF$Voca <- readxl::read_xlsx(str_c(inpath,"/eSetvoca.FF.xlsx"))



# Check missings features(<66.67% detection rate should be deleted)
eSet$S1$Data %>%
  naniar::miss_var_summary()
eSet$S2$Data %>%
  naniar::miss_var_summary()
eSet$H1$Data %>%
  naniar::miss_var_summary()
eSet$FF$Data %>%
  naniar::miss_var_summary()
# Change var type as defined in voca file
F1 = function(df){
  eSet[[df]]$Data %<>% 
    dplyr::mutate(across(.cols = all_of(eSet[[df]]$Voca %>% 
                                          dplyr::filter(VarType=="Numeric") %>% 
                                          dplyr::select(VarName) %>% as.matrix() %>% as.vector),
                        .fns = as.numeric))
  eSet[[df]]$Data %<>% 
    dplyr::mutate(across(.cols = all_of(eSet[[df]]$Voca %>% 
                                          dplyr::filter(VarType=="Factor") %>% 
                                          dplyr::select(VarName) %>% as.matrix() %>% as.vector),
                        .fns = as.factor))
  return(eSet)
}
F1("S1")
F1("S2")
F1("H1")
F1("FF")

# Scale for covariates
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
# For eSet$S1$Data, extract var names that to be scaled
c(eSet$S1$Data %>%
    dplyr::select(X1:X32) %>% #demographic variables
    names(),
  eSet$S1$Voca %>% 
    dplyr::filter(OmicGroup=="EXP1") %>% 
    dplyr::select(VarName) %>% as.matrix() %>% as.vector()
  ) -> Vars 
# scale execution
purrr::map_dfc(Vars, FuncScaleRange, df = eSet$S1$Data) %>%
  cbind(eSet$S1$Data %>%
          dplyr::select(-all_of(Vars))) %>%
  as_tibble() %>%
  dplyr::select(names(eSet$S1$Data)) -> eSet$S1$Data
# For eSet$S2$Data, extract var names that to be scaled
c(eSet$S2$Data %>%
    dplyr::select(X1:X32) %>% #demographic variables
    names(),
  eSet$S2$Voca %>% 
    dplyr::filter(OmicGroup=="EXP2") %>% 
    dplyr::select(VarName) %>% as.matrix() %>% as.vector()
  ) -> Vars 
# scale execution
purrr::map_dfc(Vars, FuncScaleRange, df = eSet$S2$Data) %>%
  cbind(eSet$S2$Data %>%
          dplyr::select(-all_of(Vars))) %>%
  as_tibble() %>%
  dplyr::select(names(eSet$S2$Data)) -> eSet$S2$Data
# For eSet$FF$Data, extract var names that to be scaled
c(eSet$FF$Data %>%
    dplyr::select(X1:X32) %>% #demographic variables
    names(),
  eSet$FF$Voca %>% 
    dplyr::filter(OmicGroup=="EXPff") %>% 
    dplyr::select(VarName) %>% as.matrix() %>% as.vector()
  ) -> Vars 
# scale execution
purrr::map_dfc(Vars, FuncScaleRange, df = eSet$FF$Data) %>%
  cbind(eSet$FF$Data %>%
          dplyr::select(-all_of(Vars))) %>%
  as_tibble() %>%
  dplyr::select(names(eSet$FF$Data)) -> eSet$FF$Data
# For eSet$H1$Data, extract var names that to be scaled
c(eSet$H1$Data %>%
    dplyr::select(X1:X32) %>% #demographic variables
    names(),
  eSet$H1$Voca %>% 
    dplyr::filter(OmicGroup=="EXPhr") %>% 
    dplyr::select(VarName) %>% as.matrix() %>% as.vector()
  ) -> Vars 
# scale execution
purrr::map_dfc(Vars, FuncScaleRange, df = eSet$H1$Data) %>%
  cbind(eSet$H1$Data %>%
          dplyr::select(-all_of(Vars))) %>%
  as_tibble() %>%
  dplyr::select(names(eSet$H1$Data)) -> eSet$H1$Data


#--2. Function for execute Benchmarking for each Omic Group-------
BchMark = function(OmicGrp, VarsX, Data, Voca){
  # OmicGrp: A chr. that indicate the group that the variables belong to;
  # VarsX: A chr. vector, which indicates variables to be used for benchmarking;
  # Data: data that contains the variables;
  # Voca: vocabulary dataframe that contains the detailed information on variables in Data.

  # creat folder
  if(!file.exists(str_c(outpath,"/1_",OmicGrp ))){
    dir.create(str_c(outpath,"/1_",OmicGrp ))
  }
  if(!file.exists(str_c(outpath,"/1_",OmicGrp ,"/Features_Explain"))){
    dir.create(str_c(outpath,"/1_",OmicGrp ,"/Features_Explain"))
  }
  if(!file.exists(str_c(outpath,"/1_",OmicGrp ,"/Features_Selection"))){
    dir.create(str_c(outpath,"/1_",OmicGrp ,"/Features_Selection"))
  }
  if(!file.exists(str_c(outpath,"/1_",OmicGrp ,"/Model_Summary"))){
    dir.create(str_c(outpath,"/1_",OmicGrp ,"/Model_Summary"))
  }

  # Outcome
  VarY = "Y"
  # n folds for cv resampling
  FoldRsmp = 5
  AutoTuneN = 50
  VarsImpThr = 0.85

  eSet$Data = Data
  eSet$Voca = Voca
  # divide data into train data and test data
  eSet$Data %>%
    dplyr::select(VarY, VarsX) -> df.train

  # set task
  df.train %>%
    as_task_classif(target = VarY) -> task0.train

  # 2.1-Select learner with best performance----
  # 2.1.1 set uniform resampling results----
  set.seed(123)#keep rsmp results as the same
  rsmp = rsmp("cv", folds = FoldRsmp)
  rsmp$instantiate(task0.train)

  # 2.1.2 set learners----
  # # -lasso-
  # learner_lasso = lrn(str_c("classif",".cv_glmnet"),
  #                     predict_type = "prob",
  #                     predict_sets = c("train", "test"),
  #                     id = "lasso")
  # learner_lasso$param_set$values = list(alpha = 1)

  # # -elastic net learner(needs to be tuned)-
  # learner_elasticnet = lrn(str_c("classif",".glmnet"),
  #                          predict_type = "prob",
  #                          predict_sets = c("train", "test"),
  #                          id = "elastic net")
  # # define autotune rules
  # at = AutoTuner$new(
  #   learner = learner_elasticnet,
  #   resampling = rsmp("cv",folds = FoldRsmp),
  #   measure = msr("classif.acc"),
  #   search_space = ps(
  #     alpha = p_dbl(lower = 0, upper = 1),
  #     lambda = p_dbl(lower = 0, upper = 1)
  #   ),
  #   terminator = trm("evals", n_evals = AutoTuneN),
  #   tuner = tnr("random_search")
  # )
  # # autotune
  # ddpcr::quiet(
  #   at$train(task0.train)
  # )
  # # use the tuned parameters
  # learner_elasticnet$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  # -random forest(needs to be tuned)----------------
  learner_rf = lrn(str_c("classif",".ranger"),
                   predict_sets = c("train", "test"),
                   predict_type = "prob")
  # define autotune
  at = AutoTuner$new(
    learner = learner_rf,
    resampling = rsmp("cv",folds = FoldRsmp),
    measure = msr("classif.acc"),
    search_space = ps(
      max.depth = p_int(lower = 6, upper = 30),
      mtry.ratio = p_dbl(lower = 0.3, upper = 0.7),
      num.trees = p_int(lower = 10, upper = 100)
    ),
    terminator = trm("evals", n_evals = AutoTuneN),
    tuner = tnr("random_search")
  )
  # autotune
  ddpcr::quiet(
    at$train(task0.train)
  )
  # set tuned values for parameter
  learner_rf$param_set$values <- at$model$tuning_instance$result_learner_param_vals
  learner_rf$param_set$values$importance <- "permutation"

  # -xgboost(needs to be tuned)----------------------------
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
  # set tuned values for parameters
  learner_xg$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  # -Naive bayes----------------------------------------
  learner_naive.bayes = lrn(str_c("classif",".naive_bayes"),
                        predict_type = "prob",
                        predict_sets = c("train", "test"),
                        id = "naive_bayes")

  # -rpart(needs to be tuned)-----------------------------
  learner_rpart = lrn(str_c("classif",".rpart"),
                      predict_type = "prob",
                      predict_sets = c("train", "test"),
                      id = "rpart")
  # define autotune rules
  at = AutoTuner$new(
    learner = learner_rpart,
    resampling = rsmp("cv",folds = FoldRsmp),
    measure = msr("classif.acc"),
    search_space = ps(
      cp = p_dbl(lower = 0.001, upper = 0.1),
      maxdepth = p_int(lower = 6, upper = 30),
      minsplit = p_int(lower = 1, upper = 10)
    ),
    terminator = trm("evals", n_evals = AutoTuneN),
    tuner = tnr("random_search")
  )
  # autotune
  ddpcr::quiet(
   at$train(task0.train)
  )
  # use tuned parameters
  learner_rpart$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  # -lightgbm(needs to be tuned)-----------------------------
  learner_lightgbm = lrn(str_c("classif",".lightgbm"),
                     predict_type = "prob",
                     predict_sets = c("train", "test"),
                     id = "lightgbm")
  learner_lightgbm$param_set$values$learning_rate = 0.1 # similar as eta in xgboost
  # define autotune rules
  at = AutoTuner$new(
    learner = learner_lightgbm,
    resampling = rsmp("cv",folds = FoldRsmp),
    measure = msr("classif.acc"),
    search_space = ps(
      num_leaves = p_int(lower = 5, upper = 50),
      max_depth = p_int(lower = 3, upper = 10),
      feature_fraction = p_dbl(lower = 0, upper = 1)
    ),
    terminator = trm("evals", n_evals = AutoTuneN),
    tuner = tnr("random_search")
  )
  # autotune
  ddpcr::quiet(
   at$train(task0.train)
  )
  # use tuned parameters
  learner_lightgbm$param_set$values = at$model$tuning_instance$result_learner_param_vals

  # -svm(needs to be tuned)-------------------------------
  learner_svm <- lrn(str_c("classif",".svm"),
                     predict_type = "prob",
                     predict_sets = c("train", "test"),
                     id = "svm")
  learner_svm$param_set$values = list(type = "C-classification")
  # auto-tune svm
  at = AutoTuner$new(
    learner = learner_svm,
    resampling = rsmp("cv",folds = FoldRsmp),
    measure = msr("classif.acc"),
    search_space = ps(
      cost = p_dbl(lower = 0.1, upper = 20),
      # gamma = p_dbl(lower = 0.1, upper = 10),
      kernel = p_fct(c("linear", "radial", "sigmoid"))
    ),
    terminator = trm("evals", n_evals = AutoTuneN),
    tuner = tnr("random_search")
  )
  # autotune
  ddpcr::quiet(
   at$train(task0.train)
  )
  # use tuned parameters
  learner_svm$param_set$values <- at$model$tuning_instance$result_learner_param_vals

  # 2.1.3-benchmarking -----------------
  # design benchmarking
  design = benchmark_grid(
    tasks = task0.train,
    learners = c(
                #  learner_lasso,
                #  learner_elasticnet,
                 learner_rf,
                 learner_xg,
                 learner_naive.bayes,
                 learner_rpart,
                 learner_lightgbm,
                 learner_svm),
    resamplings = rsmp)
  # execute benchmarking-----------
  ddpcr::quiet(
    benchmark(design) -> bmr
  )

  ## if u wanna save the learners?
  # eSet$Learners[[OmicGrp]]$lasso = learner_lasso
  # eSet$Learners[[OmicGrp]]$elasticnet = learner_elasticnet
  # eSet$Learners[[OmicGrp]]$rf = learner_rf
  # eSet$Learners[[OmicGrp]]$xg = learner_xg

  # extract accuracy
  measures_bmr = list(msr("classif.acc", predict_sets = "train", id = str_c("acc","_train")),
                      msr("classif.acc", id =str_c("acc","_test")),
                      #msr("classif.auc", predict_sets = "train", id = str_c("auc","_train")),
                      #msr("classif.auc", id =str_c("auc","_test")),
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
    unique()  ->  eSet$ModelStat[[OmicGrp]]

  # export
  eSet$ModelStat[[OmicGrp]]  %>%
    writexl::write_xlsx(str_c(outpath,"/1_",OmicGrp,"/Model_Summary/Benchmark_result.xlsx"))
  print(str_c("Benchmarking has been done for Subgroup:", OmicGrp, " !"))
  print(eSet$ModelStat[[OmicGrp]])
  #------------------------------------------------------

  return(eSet)

}

# 4----define inputs into function BchMark--------
# Omic groups
unique(eSet$S1$Voca$OmicGroup)
Omics = c("Baseline",
          "Clinical.Factors",
          "Serum1",
          "Serum2",
          "FF",
          "Hair")

# Exposures
Xlist = list()
OmicGrp = "Baseline"
eSet$S1$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp) %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[1]]
OmicGrp = c("CliniBase", "CliniIVF")
eSet$S1$Voca %>%
    dplyr::filter(OmicGroup  %in%  OmicGrp) %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[2]]
OmicGrp = "Serum1"
eSet$S1$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp | OmicGroup == "EXP1") %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[3]] 
OmicGrp = "Serum2"
eSet$S2$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp| OmicGroup == "EXP2") %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[4]] 
OmicGrp = "FF"
eSet$FF$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp| OmicGroup == "EXPff") %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[5]] 
OmicGrp = "Hair"
eSet$H1$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp| OmicGroup == "EXPhr") %>%
    dplyr::select(VarName) %>%
    as.matrix() %>% as.vector()-> Xlist[[6]] 
rm(OmicGrp)

# Data and Voca list
DataList = list(eSet$S1$Data, eSet$S1$Data, eSet$S1$Data, eSet$S2$Data, eSet$FF$Data, eSet$H1$Data)
VocaList = list(eSet$S1$Voca, eSet$S1$Voca, eSet$S1$Voca, eSet$S2$Voca, eSet$FF$Voca, eSet$H1$Voca)

# #execute
# BchMark(Omics[[1]], Xlist[[1]], DataList[[1]], VocaList[[1]])
purrr::pmap(list(Omics, Xlist, DataList, VocaList), BchMark)


#save rData
save.image(str_c(outpath, "/eSet_BenchMark.EXPs.RData"))


# # 5---plot------------
load(str_c(outpath, "/eSet_BenchMark.EXPs.RData"))
eSet$ModelStat
map_dfr(1:6, ~ eSet$ModelStat[[.x]] %>% 
            dplyr::mutate(OmicGroup = list.names(eSet$ModelStat)[[.x]]) %>% 
            dplyr::select(-rsmp_method) %>% 
            dplyr::select(OmicGroup, everything())
       ) -> plotdata
# tidy
plotdata %<>% 
  dplyr::mutate(learner_id = case_when(
                  learner_id == "classif.ranger" ~ "Random Forest", 
                  learner_id == "classif.xgboost" ~ "XGBoost", 
                  learner_id == "lightgbm" ~ "Lightgbm", 
                  learner_id == "naive_bayes" ~ "Naive Bayes", 
                  learner_id == "rpart" ~ "Decision Tree", 
                  learner_id == "svm" ~ "SVM"), 
                OmicGroup = if_else(OmicGroup == "FF", "Follicular.Fluid", OmicGroup)
                ) %>% 
  dplyr::rename(
                ACC_train = "acc_train", 
                ACC_test = "acc_test", 
                PPV_train = "precision_train", 
                PPV_test = "precision_test", 
                NPV_train = "npv_train", 
                NPV_test = "npv_test", 
                TPR_train = "sensitivity_train", 
                TPR_test = "sensitivity_test", 
                TNR_train = "specificity_train", 
                TNR_test = "specificity_test")

#-using ggiraphExtra----------------------------------------------
# devtools::install_github("cardiomoon/ggiraphExtra")
# pacman::p_load("ggiraphExtra")
library("ggiraphExtra")
plotdata %>% 
  dplyr::select(OmicGroup, learner_id, contains("_test")) %>%
  set_names(c("OmicGroup", "learner_id", "ACC", "PPV", "TPR", "TNR", "NPV")) %>% 
  ggiraphExtra::ggRadar(aes(color = learner_id, facet = OmicGroup), 
                        rescale = FALSE) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),  # Customize axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    panel.grid.major = element_line(color = "grey80"),  # Customize grid lines
    #legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10), 
    strip.text.x = element_text(size = 14, face = "bold", color = "black")
  ) +  # Customize legend text
  scale_y_continuous(
    limits = c(0, 0.9),  # Ensure the limits go from 0 to 1
    breaks = seq(0, 8, by = 0.2)  # Set breaks to ensure equal spacing
  ) +
  scale_fill_manual(values = c("#2F4F4F","#FFB366","#0099B4","#ED0000","#925E9F","#42b540"), name = "Algorithms") +  # Custom fill colors
  scale_color_manual(values = c("#2F4F4F","#FFB366","#0099B4","#ED0000","#925E9F","#42b540"), name = "Algorithms")  # Custom line colors
ggsave(str_c(outpath, "/BenchMarking.EXP.Addition-Test.dataset.png"), width = 12, height = 7)


plotdata %>% 
  dplyr::select(OmicGroup, learner_id, contains("_train")) %>%
  set_names(c("OmicGroup", "learner_id", "ACC", "PPV", "TPR", "TNR", "NPV")) %>% 
  ggiraphExtra::ggRadar(aes(color = learner_id, facet = OmicGroup), 
                        rescale = FALSE) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),  # Customize axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    panel.grid.major = element_line(color = "grey80"),  # Customize grid lines
    #legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10), 
    strip.text.x = element_text(size = 14, face = "bold", color = "black")
  ) +  # Customize legend text
  scale_y_continuous(
    limits = c(0, 1),  # Ensure the limits go from 0 to 1
    breaks = seq(0, 1, by = 0.25)  # Set breaks to ensure equal spacing
  ) +
  scale_fill_manual(values = c("#2F4F4F","#FFB366","#0099B4","#ED0000","#925E9F","#42b540"), name = "Algorithms") +  # Custom fill colors
  scale_color_manual(values = c("#2F4F4F","#FFB366","#0099B4","#ED0000","#925E9F","#42b540"), name = "Algorithms")  # Custom line colors
ggsave(str_c(outpath, "/BenchMarking.EXP.Addition-Train.dataset.png"), width = 12, height = 7)




# ggsave(str_c(outpath, "/BenchMarking.EXP.Addition.png"), width = 12, height = 5)


