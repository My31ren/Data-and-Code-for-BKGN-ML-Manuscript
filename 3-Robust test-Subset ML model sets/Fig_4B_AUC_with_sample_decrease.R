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

# - extract original results: 90%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.9", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "90%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.9", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "90%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.9", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "90%") -> m3
rbind(m1, m2, m3) -> m_0.9

# - extract original results: 80%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.8", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "80%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.8", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "80%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.8", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "80%") -> m3
rbind(m1, m2, m3) -> m_0.8

# - extract original results: 70%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.7", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "70%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.7", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "70%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.7", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "70%") -> m3
rbind(m1, m2, m3) -> m_0.7

# - extract original results: 60%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.6", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "60%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.6", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "60%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.6", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "60%") -> m3
rbind(m1, m2, m3) -> m_0.6

# - extract original results: 50%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.5", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "50%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.5", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "50%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.5", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "50%") -> m3
rbind(m1, m2, m3) -> m_0.5

# - extract original results: 40%------
work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.4", 
                    "/output/1-No.Biolink", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Reference",
                  SubData = "40%") -> m1

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.4", 
                    "/output/2-GO.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "GO.Addition",
                  SubData = "40%") -> m2

work.dir = getwd()
load(str_c(work.dir, "/@-Knowledge.Driven.model-using subset data-0.4", 
                    "/output/3-Protein.Addition", 
                    "/eSet_FitFinalModel-SG.Model.rdata"))
eSet$SG.AUC %>% 
    unlist()  %>% 
    as_tibble() %>% 
    dplyr::rename(AUC = "value") %>% 
    dplyr::mutate(SG.model = "Protein.Addition",
                  SubData = "40%") -> m3
rbind(m1, m2, m3) -> m_0.4





rbind(m_0.9, m_0.8, m_0.7, m_0.6, m_0.5, m_0.4)  %>% 
    add_row(AUC = 0.819, SG.model = "Reference", SubData = "100%") %>% 
    add_row(AUC = 0.876, SG.model = "GO.Addition", SubData = "100%") %>% 
    add_row(AUC = 0.823, SG.model = "Protein.Addition", SubData = "100%") %>% 
    dplyr::mutate(SubData = factor(SubData, levels = c("100%", "90%", "80%", "70%", "60%", "50%", "40%")))  -> plotdata

plotdata %>% 
    group_by(SG.model, SubData) %>%
    dplyr::mutate(AUC.mean = mean(AUC), 
                  AUC.sd = sd(AUC)) %>% 
    dplyr::select(-AUC) %>% 
    mutate_all(~if_else(is.na(.), 0, .)) %>% 
    unique() %>% 
    print(n=50) -> plt


plt %>% 
    dplyr::mutate(SG.model = factor(SG.model, levels = c("Reference", "GO.Addition", "Protein.Addition"))) %>% 
    ggplot(aes(x = SubData, y = AUC.mean, group = SG.model, color = SG.model, fill = SG.model)) + 
    geom_line(position=position_dodge(0.5)) + 
    geom_point(size = 7, position=position_dodge(0.5)) + 
    geom_errorbar(ymin = plt$AUC.mean - plt$AUC.sd,
                  ymax = plt$AUC.mean + plt$AUC.sd, 
                  width=.5, 
                  position=position_dodge(0.5)) + 
    scale_y_continuous(limits = c(0.71, 0.95), expand = c(0.01, 0.01, 0.01, 0.00)) + 
    scale_color_manual(values = c("#99A8B0","#F47F72", "#8DD2C5")) + 
    ylab("Mean AUC(sd) values among 10 randomly subdata") + 
    xlab("Proportions of subdata retained for modelling") + 
    theme_bw() + 
    theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks.length.x = unit(0.2, "cm"), 
    axis.ticks.length.y = unit(0.2, "cm"), 
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
    )


ggsave(str_c(getwd(), "/Robust_test_of_SG.model.png"), 
       height = 11, 
       width = 13, 
       dpi = 600)



#--statistical test--------
extract.ttest = function(size){
    plotdata %>% 
        dplyr::filter(SubData == size & SG.model == "Reference") %>% 
        dplyr::select(AUC) -> ref
    plotdata %>% 
        dplyr::filter(SubData == size & SG.model == "GO.Addition") %>% 
        dplyr::select(AUC) -> GO
    plotdata %>% 
        dplyr::filter(SubData == size & SG.model == "Protein.Addition") %>% 
        dplyr::select(AUC) -> Protein

    rst = tibble(SammpleSize = size, 
                 grp1 = "Reference", 
                 grp2 = "GO", 
                 T.test.P = t.test(ref, GO)[[3]]) %>% # 3 corresponds to p value
          add_row(SammpleSize = size, 
                  grp1 = "Reference", 
                  grp2 = "Protein", 
                  T.test.P = t.test(ref, Protein)[[3]]) %>% 
          add_row(SammpleSize = size, 
                  grp1 = "GO", 
                  grp2 = "Protein", 
                  T.test.P = t.test(GO, Protein)[[3]])
              
    return(rst)
}
library(writexl)
map_dfr(c("90%", "80%", "70%", "60%", "50%", "40%"), extract.ttest) %>% 
    write_xlsx(str_c(getwd(), "/Ttest_AUC.xlsx"))
