rm(list=ls())
library(RCy3)
library(tidyverse)
library(readxl)
library(writexl)
library(vroom)
# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

#read data
vroom::vroom(str_c(inpath, "/TableS2_exe_protein.csv")) %>%
    dplyr::filter(stringr::str_detect(source,"EX:E")|stringr::str_detect(target,"EX:P")) -> DB_Chemical_EXP
colnames(DB_Chemical_EXP)
vroom::vroom(str_c(inpath, "/TableS4_protein_disease.csv"))  %>%
    dplyr::filter(stringr::str_detect(source,"EX:P")|stringr::str_detect(target,"EX:D")) -> DB_EXP_Disease
colnames(DB_EXP_Disease)

unique(DB_Chemical_EXP$species)
DB_Chemical_EXP %>%
    dplyr::mutate(species = group_name) -> DB_Chemical_EXP

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
    dplyr::mutate(id = preferred.name.x, label = preferred.name.x, group = "EXP",color = "EXP") %>%
    dplyr::distinct(id,.keep_all = TRUE) %>%
    dplyr::select(id,label,group,color) -> df_chemexp

#extract EX:P and prefername data
DB_Chemical_EXP %>% 
    dplyr::select(target, preferred.name.x) %>% 
    rename(source = "target") %>% 
    unique() -> names

DB_EXP_Disease %>%
    left_join(names, by = "source") %>% 
    dplyr::mutate(id = preferred.name.x, label = preferred.name.x, group = "EXP",color = "EXP") %>%
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
    left_join(names, by = "source") %>% 
    dplyr::mutate(interaction = "association") %>%
    dplyr::mutate(source.class = "EXP") %>%
    dplyr::mutate(target.class = "disease") %>%
    dplyr::mutate(edge.color = "black") %>%
    dplyr::select(preferred.name.x,EXD,interaction,source.class,
                  target.class,edge.color, database) %>%
    setNames(c("source","target","interaction","source.class",
               "target.class","edge.color", "database")) -> df_edge1

DB_Chemical_EXP %>%
    dplyr::mutate(interaction = "association") %>%
    dplyr::mutate(source.class = "chemical") %>%
    dplyr::mutate(target.class = "EXP") %>%
    dplyr::select(species,preferred.name.x,interaction,source.class,
                  target.class,database) %>%
    setNames(c("source","target","interaction","source.class",
               "target.class","database")) %>%
    dplyr::mutate(edge.color = "black") -> df_edge2

rbind(df_edge1,df_edge2) %>%
    dplyr::mutate(type = dplyr::case_when(target.class == "EXP" ~ 1,
                                          target.class == "disease" ~ 2)) -> Edges_EXP
Edges_EXP %>%
    dplyr::mutate(dep = stringr::str_c(source,"-",target)) %>%
    dplyr::distinct(dep,.keep_all = TRUE) -> Edges_EXP

# Set Nodes size according to edges num -------------------------
Edges_EXP %>%
    dplyr::filter(source.class == 'chemical') %>%
    dplyr::count(source) %>%
    dplyr::arrange(desc(n)) %>%
    bind_rows(Edges_EXP %>%
                  dplyr::filter(target.class == 'EXP') %>%
                  dplyr::count(target) %>%
                  rename(source = target) %>%
                  dplyr::arrange(desc(n)))-> temp
temp %>%
    print(n=300)

Nodes_EXP %>%
    dplyr::left_join(temp,
                     by = c("id" = "source")) %>%
    dplyr::mutate(count = dplyr::case_when(group == "disease" ~ 100,
                                           group == "chemical" & n<=5 ~ 5,
                                           group == "chemical" & n>5 & n<=10 ~ 15,
                                           group == "chemical" & n>10 & n<=20 ~ 40,
                                           group == "chemical" & n>20 & n<=50 ~ 70,
                                           group == "chemical" & n>50 & n<=80 ~ 90,
                                           group == "chemical" & n>80 ~ 120,#chemicals
                                           group == "EXP" & n>=1 & n<=5 ~ 1,
                                           group == "EXP" & n>5 & n<=12 ~ 20,
                                           group == "EXP" & n>12 & n<=20 ~ 40,
                                           group == "EXP" & n>20 & n<=24 ~ 60,
                                           group == "EXP" & n>24 ~ 90)) -> Nodes_EXP
Nodes_EXP %>%
    dplyr::filter(!is.na(id)) -> Nodes_EXP

#export node and edges
Nodes_EXP %>%
    dplyr::select(-count) %>% 
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.C-AOP.Protein.Nodes.xlsx'))
Edges_EXP %>%
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.C-AOP.Protein.Edges.xlsx'))



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
                                     c("chemical","disease","EXP"),
                                     c("ELLIPSE", "ELLIPSE", "ELLIPSE"))
nodeSize <- RCy3::mapVisualProperty('Node Size','count',"d",
                                    c(20,30,40,50,60),
                                    c(20,30,40,50,60))
nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                             c("chemical","disease","EXP"),
                                             c(30, 30, 25))
nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                     c("chemical","disease","EXP"),
                                     c(40, 40, 25))
arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','type','d',
                                       c(1,2),
                                       c("Arrow","Arrow"))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location',"group",'d',
                                         c("chemical","disease","EXP"),
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



# Fig2-C2-------
temp %>% 
    print(n=200)
temp %>%
    slice(1:35) %>%
    arrange(-n) %>% 
    slice(1:10) -> plt
plt$source = factor(plt$source, levels = plt$source)
plt %>%
    ggplot(aes(x=n, y = reorder(source,n), fill = source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colorRampPalette(c("#F47F72", "#fcf6e3"))(10)) +
    scale_x_continuous(limits = c(0,70),expand=expansion(add = c(0.2, 0.2))) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.y = element_text(size = 14))
ggsave(str_c(outpath, "/Fig2-C2-Exposures.png"),
       height = 12,
       width = 8,
       dpi = 600)


temp %>%
    slice(36:194) %>%
    head(n = 10) -> b2
b2$source = factor(b2$source, levels = rev(b2$source))
b2$source
b2 %>%
    ggplot(aes(x=n, y = source, fill = source))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = colorRampPalette(c("#e5f5f2","#8DD2C5"))(10)) +
    scale_x_continuous(limits = c(0,26), breaks = c(0, 10, 20),
                       expand=expansion(add = c(0.2, 0.2))) +
    theme_bw()+ 
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.y = element_text(size = 14))
ggsave(str_c(outpath, "/Fig2-C2-Protein.links.png"),
       height = 12,
       width = 7,
       dpi = 600)

temp %>%
    writexl::write_xlsx(str_c(outpath, '/SM.Data-Fig2.C-AOP.Protein.Links.xlsx'))



