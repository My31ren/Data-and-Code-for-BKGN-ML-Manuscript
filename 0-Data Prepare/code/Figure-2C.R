rm(list=ls())
library(RCy3)
library(tidyverse)
library(readxl)
library(writexl)

# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

#read data
readxl::read_xlsx(str_c(inpath, "/Table-S2_exe_exp.xlsx")) %>%
  dplyr::filter(stringr::str_detect(EXE,"EX")|stringr::str_detect(EXP,"EX"))-> DB_Chemical_EXP
readxl::read_xlsx(str_c(inpath, "/Table-S4_exp_exd.xlsx"))  %>%
  dplyr::filter(stringr::str_detect(EXP,"EX")|stringr::str_detect(EXD,"EX")) -> DB_EXP_Disease


DB_Chemical_EXP %>%
  rename('species' = 'chemical') -> DB_Chemical_EXP

# tidy data -----------------------------
colnames(DB_EXP_Disease)
colnames(DB_Chemical_EXP)

# set disease as Pregnancy loss
DB_EXP_Disease %>%
  dplyr::mutate(EXD = "Pregnancy Loss") -> DB_EXP_Disease

# NODES --------------------------------------
DB_Chemical_EXP %>%
  dplyr::mutate(id = species, label = species, group = "chemical",color = "chemical") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_chem

DB_Chemical_EXP %>%
  dplyr::mutate(id = EXP, label = EXP, group = "EXP",color = "EXP") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_chemexp

DB_EXP_Disease %>%
  dplyr::mutate(id = EXP, label = EXP, group = "EXP",color = "EXP") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_diseaexp
DB_EXP_Disease %>%
  dplyr::mutate(id = EXD, label = EXD, group = "disease",color = "disease") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_disease

rbind(df_disease,df_chem,df_diseaexp,df_chemexp) %>%
  dplyr::distinct(id,.keep_all = TRUE)-> Nodes_EXP


# EDGES -------------------------------
DB_EXP_Disease %>%
  dplyr::mutate(interaction = "association") %>%
  dplyr::mutate(source.class = "EXP") %>%
  dplyr::mutate(target.class = "disease") %>%
  dplyr::mutate(edge.color = "black") %>%
  dplyr::select(EXP,EXD,interaction,source.class,target.class,source,edge.color) %>%
  setNames(c("source","target","interaction","source.class","target.class","database","edge.color")) -> df_edge1

DB_Chemical_EXP %>%
  dplyr::mutate(interaction = "association") %>%
  dplyr::mutate(source.class = "chemical") %>%
  dplyr::mutate(target.class = "EXP") %>%
  dplyr::select(species,EXP,interaction,source.class,target.class,source) %>%
  setNames(c("source","target","interaction","source.class","target.class","database")) %>%
  dplyr::mutate(edge.color = "black") -> df_edge2

rbind(df_edge1,df_edge2) %>%
  dplyr::mutate(type = dplyr::case_when(target.class == "EXP" ~ 1,
                                        target.class == "disease" ~ 2)) -> Edges_EXP
Edges_EXP %>%
  dplyr::mutate(dep = stringr::str_c(source,"-",target)) %>%
  dplyr::distinct(dep,.keep_all = TRUE) -> Edges_EXP

# Set Nodes size according to edges num -------------------------
Edges_EXP %>%
  dplyr::filter(!str_detect(source, '^EX:P')) %>%
  dplyr::count(source) %>%
  dplyr::arrange(desc(n)) %>%
  bind_rows(Edges_EXP %>%
              dplyr::filter(str_detect(target, 'EX:P')) %>%
              dplyr::count(target) %>%
              rename(source = target) %>%
              dplyr::arrange(desc(n)))-> temp
temp %>%
  print(n=40)

Nodes_EXP %>%
  dplyr::left_join(temp,
                   by = c("id" = "source")) %>%
  dplyr::mutate(count = dplyr::case_when(group == "disease" ~ 120,
                                         group == "chemical" & n==1 ~ 20,
                                         group == "chemical" & n==2 ~ 50,
                                         group == "chemical" & n==3 ~ 70,
                                         group == "chemical" & n==4 ~ 90,
                                         group == "chemical" & n==5 ~ 100,
                                         group == "chemical" & n>5 ~ 120,#chemicals
                                         group == "EXP" & n==1 ~ 20,
                                         group == "EXP" & n==2 ~ 30,
                                         group == "EXP" & n==3 ~ 50,
                                         group == "EXP" & n==4 ~ 70,
                                         group == "EXP" & n>4 ~ 90)) -> Nodes_EXP
Nodes_EXP %>%
  dplyr::filter(!is.na(id)) -> Nodes_EXP


#export node and edges
Nodes_EXP %>%
  writexl::write_xlsx(str_c(outpath, '/Figure-2C.Nodes.for.Plotting.xlsx'))
Edges_EXP %>%
  writexl::write_xlsx(str_c(outpath, '/Figure-2C.Edges.for.Plotting.xlsx'))


# REMENBER OPENNING Cytoscape software before running codes below!!! -----------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(Nodes_EXP,
                                  Edges_EXP,
                                  title = "EXP",
                                  collection = "biolinker")

#default style--------------------------------------------------------
style = "Style1"
defaults <- list(NODE_SHAPE="diamond",NODE_SIZE=20,
                 EDGE_TRANSPARENCY=225,NODE_LABEL_POSITION="S,W,c,0.00,0.00")
nodeLabels <- RCy3::mapVisualProperty('node label','label','p')
nodeShape <- RCy3::mapVisualProperty('Node Shape','group',"d",
                                     c("chemical","disease","go"),
                                     c("RECTANGLE", "DIAMOND", "ELLIPSE"))
nodeSize <- RCy3::mapVisualProperty('Node Size','count',"d",
                                    c(20,30,40,50,60),
                                    c(20,30,40,50,60))
nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                             c("chemical","disease","go"),
                                             c(30, 30, 25))
nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                     c("chemical","disease","go"),
                                     c(40, 40, 25))
arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','type','d',
                                       c(1,2),
                                       c("Arrow","Arrow"))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location',"group",'d',
                                         c("chemical","disease","go"),
                                         c(3,2,1))
RCy3::createVisualStyle(style,
                        defaults,
                        list(nodeLabels,
                             nodeShape,
                             nodeSize,
                             nodeLableFontSize,
                             nodeWidth,
                             arrowShapes,
                             nodeZLocation
                        ))
RCy3::setVisualStyle(style)
#plot---------------------
RCy3::setNodeColorMapping(mapping.type = 'd',
                          table.column = 'group',
                          table.column.values = c("chemical","disease","EXP"),
                          c("#D64358","#4EB3D3","#005A32"),
                          style.name = style)