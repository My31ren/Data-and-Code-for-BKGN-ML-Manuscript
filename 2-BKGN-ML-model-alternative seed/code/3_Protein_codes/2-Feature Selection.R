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

#--2. Function for Feature Selection for each Omic Group------------
Feat.Select = function(OmicGrp, VarsX, Data, Voca){
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


  # 2.1-Select learner with best performance-----------------
  # 2.1.1 set uniform resampling results---------------
  set.seed(123)#keep rsmp results as the same
  rsmp = rsmp("cv", folds = FoldRsmp)
  rsmp$instantiate(task0.train)

  # 2.1.2 set learners------------------------------------

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

 # 2.2----Feature selection using chosen algorithms----------------
  # XGboost
  # train and rank-------------------------------------------
  # train
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
  # extract rank and generate tasks list------------------------
  tasklist  = map(c(1:nrow(importance_xg)),
               ~ df.train %>%
                 dplyr::select(any_of(VarY),
                               any_of(importance_xg$Features[1:.x])) %>%
                 as_task_classif(target = VarY)
                 )
  names(tasklist) = str_c("n_", 1:nrow(importance_xg))
  # benchmarking-----------------------
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
    dplyr::rename(Feature_list = "...15") -> eSet$FeatureSelection[[OmicGrp]]
  print(as.data.frame(eSet$FeatureSelection[[OmicGrp]]))

  # export benchmarking results------------------------
  eSet$FeatureSelection[[OmicGrp]] %>%
    as_tibble() %>% 
    writexl::write_xlsx(str_c(outpath,"/1_",OmicGrp,"/Model_Summary/Model.Stat-Features.Combinations.xlsx"))
  
  # plot-----------------------------------------------
  eSet$FeatureSelection[[OmicGrp]] %>%
    as_tibble() %>%
    dplyr::filter(row.names(eSet$FeatureSelection[[OmicGrp]])=="1") %>%
    dplyr::mutate(nr_list = "n_0", 
                  nr = as.integer(0),
                  across(.cols = acc_train:npv_test, ~.x*0)
                 ) %>%
    rbind(eSet$FeatureSelection[[OmicGrp]]) %>%
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
      breaks = seq(0, nrow(eSet$FeatureSelection[[OmicGrp]]), by = 1),
      labels = str_c("n_", seq(0, nrow(eSet$FeatureSelection[[OmicGrp]]), by = 1))) +
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
  ggsave(str_c(outpath,"/1_",OmicGrp,"/Features_Selection/Feature.Selection.png"))
  
  # extract final selected features--------------------------------
  importance_xg %>%
    dplyr::slice_head(n = x.index) %>%
    left_join(eSet$Voca %>%
            dplyr::select(VarName, VarLabel, OmicGroup) %>%
            rename(Features = "VarName"),
            by = "Features") %>% 
    dplyr::mutate(Importance = if_else(Importance < 0, 0, Importance)) -> eSet$FinalFeatures[[OmicGrp]]
  eSet$FinalFeatures[[OmicGrp]] %>%
    writexl::write_xlsx(str_c(outpath,"/1_",OmicGrp,"/Features_Selection/Selected.Features.xlsx"))

  print(str_c("Final features of ", OmicGrp, " have been selected !"))
  print(eSet$FinalFeatures[[OmicGrp]])

  #####################################################
  return(eSet)

}

# ----3. Execute Feture selection FUNCTION -------------------
unique(eSet$S1$Voca$OmicGroup)
eSet$S1$Voca %<>% 
  dplyr::mutate(OmicGroup = if_else(OmicGroup  %in% c("CliniBase", "CliniIVF"),
                                    "Clinical.Factors",
                                    OmicGroup))
# define data OmicGroup, X
# Omic groups
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
OmicGrp = "Clinical.Factors"
eSet$S1$Voca %>%
    dplyr::filter(OmicGroup == OmicGrp) %>%
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
# Feat.Select(Omics[[1]], Xlist[[1]], DataList[[1]], VocaList[[1]])
purrr::pmap(list(Omics, Xlist, DataList, VocaList), Feat.Select)

#save rData
save.image(str_c(outpath, "/eSet_Feature_Select.EXPs.RData"))


