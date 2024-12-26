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

work.dir = (str_c(getwd(), "/@-Knowledge.Driven.model-SG model-Alternative.Seed.20240820"))
#--STEP1: Extract weights from SG models of Demographical, Exposures, GO, and Protein----------
# GO importance import-----------------------
df1 <- readxl::read_xlsx(str_c(work.dir,"/output/2-GO.Addition/0_Final_Model/Features_Selection/Selected.Features-SG.Model.xlsx" ))
df1.s1 <- readxl::read_xlsx(str_c(work.dir,"/output/2-GO.Addition/0_Final_Model/Features_Selection/Base.Model-Serum1Feature.Importance.xlsx" ))
df1.s2 <- readxl::read_xlsx(str_c(work.dir,"/output/2-GO.Addition/0_Final_Model/Features_Selection/Base.Model-Serum2Feature.Importance.xlsx" ))
df1.hr <- readxl::read_xlsx(str_c(work.dir,"/output/2-GO.Addition/0_Final_Model/Features_Selection/Base.Model-HairFeature.Importance.xlsx" ))
df1.ff <- readxl::read_xlsx(str_c(work.dir,"/output/2-GO.Addition/0_Final_Model/Features_Selection/Base.Model-FFFeature.Importance.xlsx" ))

# create lists for execute
df1 %>% 
  dplyr::filter(is.na(VarLabel)==T) %>% 
  dplyr::select(Importance) %>% 
  as.matrix()  %>% as.vector()  -> lst1
lst2 = list(df1.s2, df1.hr, df1.s1, df1.ff)
# function
F_mg = function(x1, x2){
  #x1: weight values
  #x2: datasets
  x2 %>% 
    dplyr::select(Features, importance, VarLabel, OmicGroup) %>% 
    dplyr::mutate(importance = importance*x1) -> tmp
  return(tmp)
}
# merge importance of SG model
df1 %>% 
  bind_rows(map2_dfr(lst1, lst2, F_mg) %>%
    dplyr::rename(Importance = "importance")) %>% # merge
  dplyr::filter(is.na(OmicGroup)==F) %>%
  dplyr::mutate(VarLabel = case_when(OmicGroup=="GO2" ~ str_replace(VarLabel, "GO_", "GOs2_"), 
                                     OmicGroup=="GO1" ~ str_replace(VarLabel, "GO_", "GOs1_"), 
                                     OmicGroup=="GOff" ~ str_replace(VarLabel, "GO_", "GOff_"), 
                                     OmicGroup=="GOhr" ~ str_replace(VarLabel, "GO_", "GOhr_"), 
                                     TRUE ~ VarLabel)) %>% 
  group_by(VarLabel) %>% 
  reframe(Importance = sum(Importance))  %>% # sum by individual group
  arrange(-Importance)  -> weight1

weight1 %<>% 
  dplyr::mutate(OmicGroup = case_when(str_detect(VarLabel, "_1$")==T ~ "Exposure.s1", 
                                      str_detect(VarLabel, "_2$")==T ~ "Exposure.s2", 
                                      str_detect(VarLabel, "_F$")==T ~ "Exposure.ff", 
                                      str_detect(VarLabel, "_H1$")==T ~ "Exposure.hr", 
                                      str_detect(VarLabel, "^GOs1")==T ~ "GO.s1", 
                                      str_detect(VarLabel, "^GOs2")==T ~ "GO.s2",
                                      str_detect(VarLabel, "^GOff")==T ~ "GO.ff",
                                      str_detect(VarLabel, "^GOhr")==T ~ "GO.hr"))

weight1 %>% 
  dplyr::filter(str_detect(OmicGroup, "^Exposure"))

# Protein importance import--------------------------------------
df2 <- readxl::read_xlsx(str_c(work.dir,"/output/3-Protein.Addition/0_Final_Model/Features_Selection/Selected.Features-SG.Model.xlsx" ))
df2.s1 <- readxl::read_xlsx(str_c(work.dir,"/output/3-Protein.Addition/0_Final_Model/Features_Selection/Base.Model-Serum1Feature.Importance.xlsx" ))
df2.s2 <- readxl::read_xlsx(str_c(work.dir,"/output/3-Protein.Addition/0_Final_Model/Features_Selection/Base.Model-Serum2Feature.Importance.xlsx" ))
df2.hr <- readxl::read_xlsx(str_c(work.dir,"/output/3-Protein.Addition/0_Final_Model/Features_Selection/Base.Model-HairFeature.Importance.xlsx" ))
df2.ff <- readxl::read_xlsx(str_c(work.dir,"/output/3-Protein.Addition/0_Final_Model/Features_Selection/Base.Model-FFFeature.Importance.xlsx" ))

# create lists for execute
df2 %>% 
  dplyr::filter(is.na(VarLabel)==T) %>% 
  dplyr::select(Importance) %>% 
  as.matrix()  %>% as.vector()  -> lst1
lst2 = list(df2.s2, df2.hr, df2.ff, df2.s1)
# function
F_mg = function(x1, x2){
  #x1: weight values
  #x2: datasets
  x2 %>% 
    dplyr::select(Features, importance, VarLabel, OmicGroup) %>% 
    dplyr::mutate(importance = importance*x1) -> tmp
  return(tmp)
}
# merge importance of SG model
df2 %>% 
  bind_rows(map2_dfr(lst1, lst2, F_mg) %>%
    dplyr::rename(Importance = "importance") %>%
    dplyr::mutate(VarLabel = case_when(OmicGroup=="EXP2" ~ str_c(VarLabel, ".s2"), 
                                     OmicGroup=="EXP1" ~ str_c(VarLabel, ".s1"), 
                                     OmicGroup=="EXPff" ~ str_c(VarLabel, ".ff"), 
                                     OmicGroup=="EXPhr" ~ str_c(VarLabel, ".hr"), 
                                     TRUE ~ VarLabel))) %>% # merge
  dplyr::filter(is.na(OmicGroup)==F) %>% 
  group_by(VarLabel) %>% 
  reframe(Importance = sum(Importance))  %>% # sum by individual group
  arrange(-Importance)  -> weight2

weight2 %<>% 
  dplyr::mutate(OmicGroup = case_when(str_detect(VarLabel, "_1$")==T ~ "Exposure.s1", 
                                      str_detect(VarLabel, "_2$")==T ~ "Exposure.s2", 
                                      str_detect(VarLabel, "_F$")==T ~ "Exposure.ff", 
                                      str_detect(VarLabel, "_H1$")==T ~ "Exposure.hr", 
                                      str_detect(VarLabel, ".s1$")==T ~ "Protein.s1", 
                                      str_detect(VarLabel, ".s2$")==T ~ "Protein.s2",
                                      str_detect(VarLabel, ".ff$")==T ~ "Protein.ff",
                                      str_detect(VarLabel, ".hr$")==T ~ "Protein.hr", 
                                      TRUE ~ "Clinical.Factors"))

#--STEP2: Extract Linked exposure-GO and exposure-protein from input file----------
vroom::vroom(str_c(work.dir, "/input/TableS3_exe_go.csv")) %>%
    dplyr::filter(stringr::str_detect(source,"EX")|stringr::str_detect(target,"GO")) -> DB_Chemical_GO
colnames(DB_Chemical_GO)
DB_Chemical_GO %>%
    dplyr::select(group_name,target) %>%
    setNames(c("source","target"))  %>% 
    unique() -> expo_go
vroom::vroom(str_c(work.dir, "/input/TableS2_exe_protein.csv")) %>%
    dplyr::filter(stringr::str_detect(source,"EX:E")|stringr::str_detect(target,"EX:P")) -> DB_Chemical_EXP
colnames(DB_Chemical_EXP)
DB_Chemical_EXP %>%
    dplyr::select(group_name,preferred.name.x) %>%
    setNames(c("source","target")) %>% 
    unique() -> expo_protein

#--STEP3: Determine which exposures, GO, and Proteins are linked in model context----------
# GO----------------
weight1 %>% 
  dplyr::filter(str_detect(OmicGroup, "^Exposure")) -> w1.expo
weight1 %>% 
  dplyr::filter(str_detect(OmicGroup, "^GO")) -> w1.go
#s1
w1.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.s1") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.s1#exposure(s1) in model
w1.go %>% 
  dplyr::filter(OmicGroup=="GO.s1") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> go.s1#GO(s1) in model
expo_go %>% 
  dplyr::mutate(source = str_c(source, "_1"), 
                target = str_replace(target, "GO:", "GOs1_")) %>% 
  dplyr::filter(source %in% expo.s1 & target %in% go.s1) %>% 
  unique()  -> edge.GOs1
#s2
w1.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.s2") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.s2#exposure(s2) in model
w1.go %>% 
  dplyr::filter(OmicGroup=="GO.s2") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> go.s2#GO(s2) in model
expo_go %>% 
  dplyr::mutate(source = str_c(source, "_2"), 
                target = str_replace(target, "GO:", "GOs2_")) %>% 
  dplyr::filter(source %in% expo.s2 & target %in% go.s2) %>% 
  unique()  -> edge.GOs2
#ff
w1.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.ff") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.ff#exposure(ff) in model
w1.go %>% 
  dplyr::filter(OmicGroup=="GO.ff") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> go.ff#GO(ff) in model
expo_go %>% 
  dplyr::mutate(source = str_c(source, "_F"), 
                target = str_replace(target, "GO:", "GOff_")) %>% 
  dplyr::filter(source %in% expo.ff & target %in% go.ff) %>% 
  unique()  -> edge.GOff
#hr
w1.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.hr") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.hr#exposure(hr) in model
w1.go %>% 
  dplyr::filter(OmicGroup=="GO.hr") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> go.hr#GO(hr) in model
expo_go %>% 
  dplyr::mutate(source = str_c(source, "_H1"), 
                target = str_replace(target, "GO:", "GOhr_")) %>% 
  dplyr::filter(source %in% expo.hr & target %in% go.hr) %>% 
  unique()  -> edge.GOhr
#linked nodes
rbind(edge.GOs1, edge.GOs2, edge.GOff, edge.GOhr) %>% 
  as.matrix() %>% as.vector() %>% 
  unique()  -> a1
#mark linked nodes with "Is.linked" feature
weight1 %<>% 
  dplyr::mutate(Is.linked = case_when(VarLabel %in% a1 ~ 1, 
                                      TRUE ~ 0))

# Protein----------------
weight2  %>% print(n=100)
weight2 %>% 
  dplyr::filter(str_detect(OmicGroup, "^Exposure")) -> w2.expo
weight2 %>% 
  dplyr::filter(str_detect(OmicGroup, "^Protein")) -> w2.protein
#s1
w2.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.s1") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.s1#exposure(s1) in model
w2.protein %>% 
  dplyr::filter(OmicGroup=="Protein.s1") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> protein.s1#GO(s1) in model
expo_protein %>% 
  dplyr::mutate(source = str_c(source, "_1"), 
                target = str_c(target, ".s1")) %>% 
  dplyr::filter(source %in% expo.s1 & target %in% protein.s1) %>% 
  unique()  -> edge.protein.s1
#s2
w2.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.s2") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.s2#exposure(s1) in model
w2.protein %>% 
  dplyr::filter(OmicGroup=="Protein.s2") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> protein.s2#GO(s1) in model
expo_protein %>% 
  dplyr::mutate(source = str_c(source, "_2"), 
                target = str_c(target, ".s2")) %>% 
  dplyr::filter(source %in% expo.s2 & target %in% protein.s2) %>% 
  unique()  -> edge.protein.s2
#ff
w2.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.ff") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.ff#exposure(s1) in model
w2.protein %>% 
  dplyr::filter(OmicGroup=="Protein.ff") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> protein.ff#GO(ff) in model
expo_protein %>% 
  dplyr::mutate(source = str_c(source, "_F"), 
                target = str_c(target, ".ff")) %>% 
  dplyr::filter(source %in% expo.ff & target %in% protein.ff) %>% 
  unique()  -> edge.protein.ff
#hr
w2.expo %>% 
  dplyr::filter(OmicGroup=="Exposure.hr") %>% 
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> expo.hr#exposure(hr) in model
w2.protein %>% 
  dplyr::filter(OmicGroup=="Protein.hr") %>%   
  dplyr::select(VarLabel) %>% 
  as.matrix() %>% as.vector() -> protein.hr#GO(hr) in model
expo_protein %>% 
  dplyr::mutate(source = str_c(source, "_H1"), 
                target = str_c(target, ".hr")) %>% 
  dplyr::filter(source %in% expo.hr & target %in% protein.hr) %>% 
  unique()  -> edge.protein.hr
#linked nodes
rbind(edge.protein.s1, edge.protein.s2, edge.protein.ff, edge.protein.hr) %>% 
  as.matrix() %>% as.vector() %>% 
  unique()  -> a2
#mark linked nodes with "Is.linked" feature
weight2 %<>% 
  dplyr::mutate(Is.linked = case_when(VarLabel %in% a2 ~ 1, 
                                      TRUE ~ 0))

#--STEP4: combine GO and Protein results in STEP3----------
weight1 %<>% 
  dplyr::mutate(pct = Importance/(.$Importance %>% sum()))
weight2 %<>% 
  dplyr::mutate(pct = Importance/(.$Importance %>% sum()))

# rbind nodes
rbind(weight1, weight2) %>% 
  dplyr::mutate(VarLabel = str_replace(VarLabel, "GOs1_|GOs2_|GOff_|GOhr_", "GO:"), 
                VarLabel = str_replace(VarLabel, ".s1|.s2|.ff|.hr", "")) %>%
  dplyr::mutate(OmicGroup = str_replace(OmicGroup, ".s1|.s2|.ff|.hr", "")) %>% 
  group_by(VarLabel) %>% 
  reframe(pct = sum(pct), 
          OmicGroup = OmicGroup, 
          Is.linked = max(Is.linked)) %>% 
  arrange(-pct)  %>% unique() -> pct
# rbing edges
rbind(edge.GOs1, edge.GOs2, edge.GOff, edge.GOhr) %>% 
rbind(edge.protein.s1, edge.protein.s2, edge.protein.ff, edge.protein.hr) %>% 
  dplyr::mutate(target = str_replace(target, "GOs1_|GOs2_|GOff_|GOhr_", "GO:"), 
                target = str_replace(target, ".s1|.s2|.ff|.hr", "")) %>%
  unique() %>% 
  dplyr::mutate(source.class = case_when(
                  str_detect(source, "_1$") ~ "Serum #1", 
                  str_detect(source, "_2$") ~ "Serum #2", 
                  str_detect(source, "_F$") ~ "Follicular Fluid", 
                  str_detect(source, "_H1$") ~ "Hair"), 
                target.class = case_when(
                  str_detect(target, "^GO") ~ "GO", 
                  TRUE ~ "Protein", 
                ), 
                interaction = "association",
                edge.color = "black") -> edges

# STEP5: network plot----------------------------------
library(RCy3)
# build nodes
pct %>% 
  dplyr::mutate(id = VarLabel, 
                label = VarLabel, 
                group = case_when(
                  str_detect(VarLabel, "_1$") ~ "Serum #1", 
                  str_detect(VarLabel, "_2$") ~ "Serum #2", 
                  str_detect(VarLabel, "_F$") ~ "Follicular Fluid", 
                  str_detect(VarLabel, "_H1$") ~ "Hair", 
                  str_detect(VarLabel, "^GO") ~ "GO",
                  TRUE ~ "Protein"), 
                color = group, 
                pct = pct*100) -> nodes
# export-------------------------
edges %>% 
  writexl::write_xlsx(str_c(work.dir, "/output/Figure4-edges.xlsx"))
nodes %>% 
  writexl::write_xlsx(str_c(work.dir, "/output/Figure4-nodes.xlsx"))

# NOTE: you should open Cytoscape software (we used version 3.9.1)------------------------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(node,
                                  edge,
                                  title = "GO-EXP",
                                  collection = "biolinker")

#default style--------------------------------------------------------
style = "Style1"
defaults <- list(NODE_SHAPE="diamond",NODE_SIZE=20,
                 EDGE_TRANSPARENCY=225,NODE_LABEL_POSITION="S,W,c,0.00,0.00")
nodeLabels <- RCy3::mapVisualProperty('node label','label','p')
nodeShape <- RCy3::mapVisualProperty('Node Shape','group',"d",
                                     c("chemical","EXP","go"),
                                     c("RECTANGLE", "ELLIPSE", "ELLIPSE"))
nodeSize <- RCy3::mapVisualProperty('Node Size','n',"d",
                                    c(20,30,40,50,60),
                                    c(20,30,40,50,60))
nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                             c("chemical","disease","go"),
                                             c(30, 30, 25))
nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                     c("chemical","EXP","go"),
                                     c(40, 40, 25))
# arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','type','d',
#                                        c(1,2),
#                                        c("Arrow","Arrow"))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location',"group",'d',
                                         c("chemical","EXP","go"),
                                         c(3,2,1))
RCy3::createVisualStyle(style,
                        defaults,
                        list(nodeLabels,
                             #nodeFills,
                             nodeShape,
                             nodeSize,
                             nodeLableFontSize,
                             #nodeLableColor,
                             nodeWidth,
                             # arrowShapes,
                             nodeZLocation
                             # edgeWidth
                        ))
RCy3::setVisualStyle(style)
#plot---------------------
RCy3::setNodeColorMapping(mapping.type = 'd',
                          table.column = 'group',
                          table.column.values = c("chemical","disease","go"),
                          c("#D64358","#4EB3D3","#005A32"),
                          style.name = style)