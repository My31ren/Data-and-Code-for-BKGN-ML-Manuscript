rm(list=ls())
library(RCy3)
library(tidyverse)
library(readxl)
library(writexl)

# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

read_csv(str_c(inpath, "/TableS3_exe_go.csv")) %>%
  dplyr::filter(stringr::str_detect(EXE,"EX")|stringr::str_detect(go.id,"GO")) -> DB_Chemical_GO
readxl::read_xlsx(str_c(inpath, "/TableS5_go_exd.xlsx"))  %>%
  dplyr::filter(stringr::str_detect(go.id,"GO")|stringr::str_detect(EXD,"EX")) -> DB_GO_Disease

DB_Chemical_GO %>%
  rename('species' = 'chemical') -> DB_Chemical_GO

# tidy data --------------------------
colnames(DB_GO_Disease)
colnames(DB_Chemical_GO)

# set disease as Pregnancy loss
DB_GO_Disease %>%
  dplyr::mutate(EXD = "Pregnancy Loss") -> DB_GO_Disease

# NODES ---------------------------------
DB_Chemical_GO %>%
  dplyr::mutate(id = species, label = species, group = "chemical",color = "chemical") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_chem

DB_Chemical_GO %>%
  dplyr::mutate(id = go.id, label = go.id, group = "go",color = "go") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_chemgo

DB_GO_Disease %>%
  dplyr::mutate(id = go.id, label = go.id, group = "go",color = "go") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_diseago
DB_GO_Disease %>%
  dplyr::mutate(id = EXD, label = EXD, group = "disease",color = "disease") %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> df_disease

rbind(df_disease,df_chem,df_diseago,df_chemgo) %>%
  dplyr::distinct(id,.keep_all = TRUE)-> Nodes_GO


# EDGES -----------------------------
DB_GO_Disease %>%
  dplyr::mutate(interaction = "association") %>%
  dplyr::mutate(source.class = "GO") %>%
  dplyr::mutate(target.class = "disease") %>%
  dplyr::mutate(edge.color = "black") %>%
  dplyr::select(go.id,EXD,interaction,source.class,target.class,source,edge.color) %>%
  setNames(c("source","target","interaction","source.class","target.class","database","edge.color")) -> df_edge1

DB_Chemical_GO %>%
  dplyr::mutate(interaction = "association") %>%
  dplyr::mutate(source.class = "chemical") %>%
  dplyr::mutate(target.class = "GO") %>%
  dplyr::select(species,go.id,interaction,source.class,target.class,source) %>%
  setNames(c("source","target","interaction","source.class","target.class","database")) %>%
  dplyr::mutate(edge.color = "black") -> df_edge2

rbind(df_edge1,df_edge2) %>%
  dplyr::mutate(type = dplyr::case_when(target.class == "GO" ~ 1,
                                        target.class == "disease" ~ 2)) -> Edges_GO
Edges_GO %>%
  dplyr::mutate(dep = stringr::str_c(source,"-",target)) %>%
  dplyr::distinct(dep,.keep_all = TRUE) -> Edges_GO

# Set Nodes size according to edges num -------------------------
Edges_GO %>%
  dplyr::filter(!str_detect(source, '^GO')) %>%
  dplyr::count(source) %>%
  dplyr::arrange(desc(n)) %>%
  bind_rows(Edges_GO %>%
              dplyr::filter(str_detect(target, 'GO')) %>%
              dplyr::count(target) %>%
              rename(source = target) %>%
              dplyr::arrange(desc(n)))-> temp
temp %>%
  print(n = 700)

Nodes_GO %>%
  dplyr::left_join(temp,
                   by = c("id" = "source")) %>%
  dplyr::mutate(count = dplyr::case_when(group == "disease" ~ 120,
                    group == "chemical" & n>=1 & n<=5 ~ 20,
                    group == "chemical" & n>5 & n<=15 ~ 50,
                    group == "chemical" & n>15 & n<=30 ~ 70,
                    group == "chemical" & n>30 & n<=100 ~ 90,
                    group == "chemical" & n>100 & n<=200 ~ 100,
                    group == "chemical" & n>200 ~ 120,#chemicals
                    group == "go" & n>=1 & n<=5 ~ 1,
                    group == "go" & n>5 & n<=7 ~ 30,
                    group == "go" & n>7 & n<=9 ~ 50,
                    group == "go" & n>9 & n<=11 ~ 70,
                    group == "go" & n>11 ~ 90)) -> Nodes_GO
Nodes_GO %>%
  dplyr::filter(!is.na(id)) -> Nodes_GO

#export node and edges
Nodes_GO %>%
  writexl::write_xlsx(str_c(outpath, '/Figure-2B.Nodes.for.Plotting.xlsx'))
Edges_GO %>%
  writexl::write_xlsx(str_c(outpath, '/Figure-2B.Edges.for.Plotting.xlsx'))


# REMENBER OPENNING Cytoscape software before running codes below!!! -----------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(Nodes_GO,
                                  Edges_GO,
                                  title = "GO",
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
                          table.column.values = c("chemical","disease","go"),
                          c("#D64358","#4EB3D3","#005A32"),
                          style.name = style)