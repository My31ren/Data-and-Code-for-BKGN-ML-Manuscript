rm(list=ls())
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
# load data--------
work.dir = (str_c(getwd(), ""))

# no biolink----
load(str_c(work.dir, "/@-Knowledge.Driven.model-SG model-Alternative.Seed.20240820/output/1-No.Biolink/eSet_FitFinalModel-SG.Model.RData"))

bmr$score(measures_bmr)

map(1:5, 
    ~ as.data.table(bmr, 
      reassemble_learners = T, 
      convert_predictions = T,
      predict_sets = "test")$prediction[[.x]]
      )  -> dt1
map_dfr(1:5, 
    ~ bind_cols(dt1[[.x]]$row_ids, dt1[[.x]]$truth, dt1[[.x]]$response, dt1[[.x]]$prob) %>% 
      set_names(c("row_ids", "truth", "response", "prob.0", "prob.1"))
      ) %>% 
    arrange(row_ids) -> plt1

#go addition---------
load(str_c(work.dir, "/output/2-GO.Addition/eSet_FitFinalModel-SG.Model.RData"))
map(1:5, 
    ~ as.data.table(bmr, 
      reassemble_learners = T, 
      convert_predictions = T,
      predict_sets = "test")$prediction[[.x]]
      )  -> dt2
map_dfr(1:5, 
    ~ bind_cols(dt2[[.x]]$row_ids, dt2[[.x]]$truth, dt2[[.x]]$response, dt2[[.x]]$prob) %>% 
      set_names(c("row_ids", "truth", "response", "prob.0", "prob.1"))
      ) %>% 
    arrange(row_ids) -> plt2

# protein addition-------
load(str_c(work.dir, "/output/3-Protein.Addition/eSet_FitFinalModel-SG.Model.RData"))
map(1:5, 
    ~ as.data.table(bmr, 
      reassemble_learners = T, 
      convert_predictions = T,
      predict_sets = "test")$prediction[[.x]]
      )  -> dt3
map_dfr(1:5, 
    ~ bind_cols(dt3[[.x]]$row_ids, dt3[[.x]]$truth, dt3[[.x]]$response, dt3[[.x]]$prob) %>% 
      set_names(c("row_ids", "truth", "response", "prob.0", "prob.1"))
      ) %>% 
    arrange(row_ids) -> plt3

# plotting---------------
p1 = roc(plt1$truth~plt1$prob.1)
p2 = roc(plt2$truth~plt2$prob.1)
p3 = roc(plt3$truth~plt3$prob.1)
ggroc(list(No.Biolink = p1, GO.Integration = p2, Protein.Integration = p3), 
      legacy.axes = TRUE, 
      linewidth = 2) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed") + 
  theme_bw() + 
  scale_y_continuous(
    limits = c(0, 1),  # Ensure the limits go from 0 to 1
    expand = expansion(add = c(0.02, 0.02))
  ) + 
  scale_x_continuous(
    limits = c(0, 1),  # Ensure the limits go from 0 to 1
    expand = expansion(add = c(0.02, 0.03))
  ) + 
  scale_color_manual(values = c("#99A8B0","#F47F72", "#8DD2C5")) + 
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks.length.x = unit(0.2, "cm"), 
    axis.ticks.length.y = unit(0.2, "cm"), 
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20),
    legend.key.size = unit(55, "pt"),
    legend.key.height = unit(35, "pt"), 
    legend.position = c(0.75, 0.2)
    )
ggsave(str_c(work.dir, "/output/Fig4-B-ROC.png"),
       height = 9,
       width = 8,
       dpi = 600)

#-retieve AUC values for each group--------
auc(p1)#NoBiolink
auc(p2)#GOAddition
auc(p3)#ProteinAddition

