rm(list=ls())
#load packages
pacman::p_load(
    "gmodels",
    "rstatix",
    "foreign",
    "haven",
    "openxlsx",
    "tidyverse",
    "naniar",
    "Hmisc",
    "mice",
    "tableone",
    "caret",
    "magrittr",
    "R6",
    "lubridate",
    "tictoc",
    "vroom",
    "rstatix",
    "ggplot2",
    "ggfortify",
    "GGally",
    "zip",
    "ddpcr",
    "grid",
    "flextable",
    "ggrepel",
    "paletteer",
    "furrr",
    "rlist")


# 1---load reference model data----------------------------
outpath = stringr::str_c(getwd(), "/2-Interpretable.Model.using.AOP/@Data and Codes/@-Knowledge.Driven.model-SG model-Alternative.Seed.20240820/output")
load(str_c(outpath,"/1-No.Biolink/eSet_Feature_Select.Rdata"))

FSplot = function(OmicGrp, sn, abbre){
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
      #breaks = seq(0, nrow(eSet$FeatureSelection[[OmicGrp]]), by = 1),
      #labels = str_c("n_", seq(0, nrow(eSet$FeatureSelection[[OmicGrp]]), by = 1))
      ) +
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
    theme_bw() + 
    theme(axis.ticks.length.x = unit(0.2,'cm'), 
          axis.ticks.length.y = unit(0.2,'cm'))
  ggsave(str_c(getwd(),"/2-Interpretable.Model.using.AOP/3-Figure S3-S6", "/Feature_", sn, "_", abbre, ".png"), 
         height = 8, 
         width = 9)
}

names(eSet$FeatureSelection)

#Figure S3
FSplot("Baseline", "S3", "BL")
FSplot("Clinical.Factors", "S3", "CF")

#Figure S4
FSplot("Serum1", "S4", "SR1")
FSplot("Serum2", "S4", "SR2")
FSplot("FF", "S4", "FF")
FSplot("Hair", "S4", "HR")


#-2 load GO integrated model-------
outpath = stringr::str_c(getwd(), "/2-Interpretable.Model.using.AOP/@Data and Codes/@-Knowledge.Driven.model-SG model-Alternative.Seed.20240820/output")
load(stringr::str_c(outpath,"/2-GO.Addition/eSet_Feature_Select.GOs.Rdata"))

names(eSet$FeatureSelection)
FSplot("Serum1", "S5", "SR1")
FSplot("Serum2", "S5", "SR2")
FSplot("FF", "S5", "FF")
FSplot("Hair", "S5", "HR")

#-3 load Protein integrated model-------
outpath = stringr::str_c(getwd(), "/2-Interpretable.Model.using.AOP/@Data and Codes/@-Knowledge.Driven.model-SG model-Alternative.Seed.20240820/output")
load(stringr::str_c(outpath,"/3-Protein.Addition/eSet_Feature_Select.EXPs.Rdata"))

names(eSet$FeatureSelection)
FSplot("Serum1", "S6", "SR1")
FSplot("Serum2", "S6", "SR2")
FSplot("FF", "S6", "FF")
FSplot("Hair", "S6", "HR")


