rm(list=ls())
library(RCy3)
library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
library(vroom)
# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

vroom::vroom(str_c(inpath, "/TableS3_exe_go.csv")) %>%
    dplyr::filter(stringr::str_detect(source,"EX")|stringr::str_detect(target,"GO")) -> DB_Chemical_GO
colnames(DB_Chemical_GO)

vroom::vroom(str_c(inpath, "/TableS5_go_disease.csv"))  %>%
    dplyr::filter(stringr::str_detect(source,"GO:")|stringr::str_detect(target,"EX:D")) -> DB_GO_Disease
colnames(DB_GO_Disease)


DB_Chemical_GO %>%
    rename('species' = 'group_name') -> DB_Chemical_GO

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
    dplyr::mutate(id = target, label = target, group = "go",color = "go") %>%
    dplyr::distinct(id,.keep_all = TRUE) %>%
    dplyr::select(id,label,group,color) -> df_chemgo

DB_GO_Disease %>%
    dplyr::mutate(id = source, label = source, group = "go",color = "go") %>%
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
    dplyr::select(source,EXD,interaction,source.class,
                  target.class, edge.color, database) %>%
    setNames(c("source","target","interaction","source.class",
               "target.class","edge.color", "database")) -> df_edge1

DB_Chemical_GO %>%
    dplyr::mutate(interaction = "association") %>%
    dplyr::mutate(source.class = "chemical") %>%
    dplyr::mutate(target.class = "GO") %>%
    dplyr::select(species,target,interaction,source.class,target.class, database) %>%
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
    print(n =200)

Nodes_GO %>%
    dplyr::left_join(temp,
                     by = c("id" = "source")) %>%
    dplyr::mutate(count = dplyr::case_when(group == "disease" ~ 100,
                                           group == "chemical" & n>=1 & n<=5 ~ 5,
                                           group == "chemical" & n>5 & n<=15 ~ 15,
                                           group == "chemical" & n>15 & n<=50 ~ 40,
                                           group == "chemical" & n>50 & n<=100 ~ 70,
                                           group == "chemical" & n>100 & n<=250 ~ 90,
                                           group == "chemical" & n>250 ~ 120,#chemicals
                                           group == "go" & n>=1 & n<=6 ~ 1,
                                           group == "go" & n>6 & n<=10 ~ 20,
                                           group == "go" & n>10 & n<=15 ~ 40,
                                           group == "go" & n>15 & n<=20 ~ 60,
                                           group == "go" & n>20  ~ 90)) -> Nodes_GO
Nodes_GO %>%
    dplyr::filter(!is.na(id)) -> Nodes_GO


#export node and edges
Nodes_GO %>%
    dplyr::select(-count) %>% 
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.B-AOP.GO.Nodes.xlsx'))
Edges_GO %>%
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.B-AOP.GO.Edges.xlsx'))

# ~------------------------- ~-----------
# Fig4A--------
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

# ~-----------------~-------

# Fig2-B2-------
temp %>%
    dplyr::filter(!str_detect(source, "^GO:")) %>%
    arrange(-n) %>% 
    slice(1:10)-> plt
plt$source = factor(plt$source, levels = plt$source)
plt %>%
    ggplot(aes(x=n, y = reorder(source,n), fill = source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colorRampPalette(c("#F47F72", "#fcf6e3"))(10)) +
    scale_x_continuous(limits = c(0,305),expand=expansion(add = c(0.2, 0.2))) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.y = element_text(size = 14))
ggsave(str_c(outpath, "/Fig2-B2-Exposures.png"),
       height = 12,
       width = 8,
       dpi = 600)


temp %>%
    dplyr::filter(str_detect(source, "^GO:")) %>%
    head(n = 10) -> b2
b2$source = factor(b2$source, levels = rev(b2$source))
b2$source
b2 %>%
    ggplot(aes(x=n, y = source, fill = source))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = colorRampPalette(c("#e3f0fc","#7FB2D5"))(10)) +
    scale_x_continuous(limits = c(0,30), breaks = c(0, 10, 20, 30),
                       expand=expansion(add = c(0.2, 0.2))) +
    theme_bw()+ 
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.y = element_text(size = 14))
ggsave(str_c(outpath, "/Fig2-B2-GO.links.png"),
       height = 12,
       width = 7,
       dpi = 600)

temp %>%
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.B-AOP.GO.Links.xlsx'))

# 
# # export links of sensitive biomarkers
# DB_Chemical_GO %>%
#     dplyr::filter(str_detect(species, "selenium|iron|zinc|cobalt")) %>%
#     dplyr::select(species, go.id, ontology, label) %>%
#     unique() %>%
#     arrange(species, ontology) %>%
#     writexl::write_xlsx(str_c(outpath, "/5.Explainer/Table S10-GO pathways for sensitive biomarkers.xlsx"))
