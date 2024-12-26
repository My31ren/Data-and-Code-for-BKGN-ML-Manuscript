library(RCy3)
library(tidyverse)
library(readxl)
library(rstatix)

# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

# 1-Calculate Figure 2A result----------------
tmp = readxl::read_xlsx(str_c(inpath,"/1-IVF-n116.xlsx"))

#prepare data
colnames(tmp)
tmp %>%
  dplyr::select(As_1:Zn_H1) %>%
  colnames() -> varlist
combn(varlist,2) %>% t() %>%
  as_tibble() -> comb

# function for spearman correlation
F1 = function(v1, v2){
  cor_test(tmp, vars = v1, vars2 = v2, method = "spearman")
}
purrr::map2_dfr(comb$V1, comb$V2, F1) -> cor.result

# Define X  第一次血清金属
tmp %>%
  dplyr::select(As_1:Li_1) %>%
  colnames() -> X1
# Define X  第二次血清金属
tmp %>%
  dplyr::select(As_2:Li_2) %>%
  colnames() -> X2
# Define X  卵泡液金属
tmp %>%
  dplyr::select(As_F:Li_F) %>%
  colnames() -> X3
# Define X  第一次血清PFAS
tmp %>%
  dplyr::select(PFBA_1:P_62Cl_PFESA_1) %>%
  colnames() -> X4
# Define X  第二次血清PFAS
tmp %>%
  dplyr::select(PFBA_2:P_62Cl_PFESA_2) %>%
  colnames() -> X5
# Define X  头发
tmp %>%
  dplyr::select(contains("_H1")) %>%
  colnames() -> X6

cor.result %>%
  dplyr::select(var1, var2, method, cor, p) %>%
  dplyr::mutate(group1 = case_when(
      var1 %in% X1 ~ "Metal.Serum.1",
      var1 %in% X2 ~ "Metal.Serum.2",
      var1 %in% X3 ~ "Metal.FollicularFluid",
      var1 %in% X4 ~ "PFAS.Serum.1",
      var1 %in% X5 ~ "PFAS.Serum.2",
      var1 %in% X6 ~ "Metal.Hair"),
                group2 = case_when(
      var2 %in% X1 ~ "Metal.Serum.1",
      var2 %in% X2 ~ "Metal.Serum.2",
      var2 %in% X3 ~ "Metal.FollicularFluid",
      var2 %in% X4 ~ "PFAS.Serum.1",
      var2 %in% X5 ~ "PFAS.Serum.2",
      var2 %in% X6 ~ "Metal.Hair")) -> out.cor.result
out.cor.result %>%
  writexl::write_xlsx(str_c(outpath, "/Figure-2A.Correlation.xlsx"))

#---------------------


#1. Nodes
varlist %>%
  as_tibble() %>%
  dplyr::mutate(id = value, label = value) %>%
  dplyr::mutate(group = case_when(
    id %in% X1 ~ "Metal.Serum.1",
    id %in% X2 ~ "Metal.Serum.2",
    id %in% X3 ~ "Metal.FollicularFluid",
    id %in% X4 ~ "PFAS.Serum.1",
    id %in% X5 ~ "PFAS.Serum.2",
    id %in% X6 ~ "Metal.Hair"),
    color = group) %>%
  dplyr::distinct(id,.keep_all = TRUE) %>%
  dplyr::select(id,label,group,color) -> nodes

#2. nodes size
tmp$X1临床妊娠_Y0_N1 = as.numeric(tmp$X1临床妊娠_Y0_N1)

purrr::map2_dfr("X1临床妊娠_Y0_N1",
                varlist %>% unique(),
                F1) -> node.size
node.size %>%
  dplyr::mutate(MinusLogP = -log(p+0.00000001),
                id = var2) %>%
  dplyr::select(MinusLogP) %>%
  bind_cols(nodes) %>%
  dplyr::select(id, label, group, color, MinusLogP) -> nodes


#3. Edges
out.cor.result %>%
  dplyr::filter(p < 0.05) %>%
  dplyr::mutate(interaction = "association") %>%
  dplyr::mutate(source.class = group1) %>%
  dplyr::mutate(target.class = group2) %>%
  dplyr::mutate(edge.color = "black") %>%
  dplyr::select(var1,var2,interaction,source.class,target.class,edge.color) %>%
  setNames(c("source","target","interaction","source.class","target.class","edge.color")) %>%
  dplyr::mutate(dep = stringr::str_c(source,"-",target)) %>%
  dplyr::distinct(dep,.keep_all = TRUE) -> edges

#4. export
writexl::write_xlsx(nodes,
                    str_c(outpath, "/Figure-2A.Nodes.for.Plotting.xlsx"))
writexl::write_xlsx(edges,
                    str_c(outpath, "/Figure-2A.Edges.for.Plotting.xlsx"))


# Open sytoscape software----------------------------------------------------------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(nodes,
                                  edges,
                                  title = "GO",
                                  collection = "biolinker")

#default style--------------------------------------------------------
nodes$group %>%
  unique()
nodes$MinusLogP %>%
  summary()
style = "Style1"
defaults <- list(NODE_SHAPE="ELLIPSE",NODE_SIZE=20,
                 EDGE_TRANSPARENCY=225,NODE_LABEL_POSITION="S,W,c,0.00,0.00")
nodeLabels <- RCy3::mapVisualProperty('node label','label','p')
nodeShape <- RCy3::mapVisualProperty('Node Shape','group',"d",
                                     nodes$group %>%
                                       unique(),
                                     c("ELLIPSE", "ELLIPSE", "ELLIPSE",
                                       "ELLIPSE", "ELLIPSE", "ELLIPSE"))
nodeSize <- RCy3::mapVisualProperty('Node Size','MinusLogP',"c",
                                    c(0.01106, 6.83078),
                                    c(100, 600))
nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                             nodes$group %>%
                                               unique(),
                                             c(30, 30, 25, 30, 30, 30))
nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                     nodes$group %>%
                                       unique(),
                                     c(30, 30, 30, 30, 30, 30))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location',"group",'d',
                                         nodes$group %>%
                                           unique(),
                                         c(6, 5, 4, 3, 2, 1))
RCy3::createVisualStyle(style,
                        defaults,
                        list(nodeLabels,
                             nodeShape,
                             nodeSize,
                             nodeLableFontSize,
                             nodeWidth,
                             nodeZLocation))
RCy3::setVisualStyle(style)

#plot---------------------
RCy3::setNodeColorMapping(mapping.type = 'd',
                          table.column = 'group',
                          table.column.values = nodes$group %>%
                            unique(),
                          c("#D64358","#4EB3D3","#005A32",
                            "#D64358","#4EB3D3","#005A32"),
                          style.name = style)
