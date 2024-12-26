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

# 0-Prepare step ----------------------------------
#set work directory
inpath = stringr::str_c(getwd(), "/input")
outpath = stringr::str_c(getwd(), "/output")

if(!file.exists(inpath)){dir.create(inpath)}
if(!file.exists(outpath)){dir.create(outpath)}

#load data
load(str_c(outpath,"/1-No.Biolink/eSet_BenchMark.Rdata"))

eSet$Data
eSet$Voca %>% 
    print(n=188)
# 2-Draw Figure1 ----------------
#prepare data
colnames(tmp)

eSet$Voca %>% 
    dplyr::filter(VarName %in% c(Xlist[3:6] %>% 
                                     unlist())) %>% 
    dplyr::select(VarLabel) %>% as.matrix() %>% as.vector() -> varname
eSet$Data %>% 
    dplyr::select(all_of(Xlist[3:6] %>% 
                             unlist())) %>% 
    set_names(varname) -> tmp

combn(varname,2) %>% t() %>%
    as_tibble() -> comb

F1 = function(v1, v2){
    cor_test(tmp, vars = v1, vars2 = v2, method = "spearman")
}

purrr::map2_dfr(comb$V1, comb$V2, F1) -> cor.result

tmp %>% 
    dplyr::select(matches("^\\w\\w?_1$")) %>% 
    colnames() %>% str_sort() -> metal_1
tmp %>% 
    dplyr::select(matches("^\\w\\w\\w+_1$")) %>% 
    colnames() %>% str_sort() -> pfas_1
tmp %>% 
    dplyr::select(matches("^\\w\\w?_2$")) %>% 
    colnames() %>% str_sort() -> metal_2
tmp %>% 
    dplyr::select(matches("^\\w\\w\\w+_2$")) %>% 
    colnames() %>% str_sort() -> pfas_2
tmp %>% 
    dplyr::select(matches("^\\w\\w?_F$")) %>% 
    colnames() %>% str_sort() -> metal_f
tmp %>% 
    dplyr::select(matches("^\\w\\w?_H1$")) %>% 
    colnames() %>% str_sort() -> metal_h

cor.result %>%
    dplyr::select(var1, var2, method, cor, p) %>%
    dplyr::mutate(group1 = case_when(
        var1 %in% metal_1 ~ "Metal.Serum.1",
        var1 %in% metal_2 ~ "Metal.Serum.2",
        var1 %in% metal_f ~ "Metal.FollicularFluid",
        var1 %in% pfas_1 ~ "PFAS.Serum.1",
        var1 %in% pfas_2 ~ "PFAS.Serum.2",
        var1 %in% metal_h ~ "Metal.Hair"),
        group2 = case_when(
            var2 %in% metal_1 ~ "Metal.Serum.1",
            var2 %in% metal_2 ~ "Metal.Serum.2",
            var2 %in% metal_f ~ "Metal.FollicularFluid",
            var2 %in% pfas_1 ~ "PFAS.Serum.1",
            var2 %in% pfas_2 ~ "PFAS.Serum.2",
            var2 %in% metal_h ~ "Metal.Hair")) -> out.cor.result
out.cor.result %>%
    writexl::write_xlsx(str_c(getwd(), "/output/Fig2-A-Correlation.Results.xlsx"))

#1. Nodes
varname %>%
    as_tibble() %>%
    dplyr::mutate(id = value, label = value) %>%
    dplyr::mutate(group = case_when(
        id %in% metal_1 ~ "Metal.Serum.1",
        id %in% metal_2 ~ "Metal.Serum.2",
        id %in% metal_f ~ "Metal.FollicularFluid",
        id %in% pfas_1 ~ "PFAS.Serum.1",
        id %in% pfas_2 ~ "PFAS.Serum.2",
        id %in% metal_h ~ "Metal.Hair"),
        color = group) %>%
    dplyr::distinct(id,.keep_all = TRUE) %>%
    dplyr::select(id,label,group,color) -> nodes

#2. nodes size
eSet$Data %>% 
    dplyr::select(Y) %>% 
    bind_cols(tmp) -> tmp
tmp$Y = as.numeric(tmp$Y)

purrr::map2_dfr("Y",
                varname %>% unique(),
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

# # delete nodes that excluded from edges
# c(edges$source %>% unique(),
#   edges$target %>% unique()) %>%
#   unique() -> idlist
# nodes %>%


#4. export
writexl::write_xlsx(nodes,
                    str_c(getwd(), "/output/SM.Data-1-Correlation.Plot.Nodes.xlsx"))
writexl::write_xlsx(edges,
                    str_c(getwd(), "/output/SM.Data-1-Correlation.Plot.Edges.xlsx"))



# Open sytoscape software----------------------------------------------------------
RCy3::cytoscapePing()
RCy3::closeSession(F)
RCy3::createNetworkFromDataFrames(nodes,
                                  edges,
                                  title = "Fig2.Correlation",
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
# arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','type','d',
#                                        c(1,2),
#                                        c("Arrow","Arrow"))
nodeZLocation <- RCy3::mapVisualProperty('Node Z Location',"group",'d',
                                         nodes$group %>%
                                             unique(),
                                         c(6, 5, 4, 3, 2, 1))
RCy3::createVisualStyle(style,
                        defaults,
                        list(nodeLabels,
                             #nodeFills,
                             nodeShape,
                             nodeSize,
                             nodeLableFontSize,
                             #nodeLableColor,
                             nodeWidth,
                             #arrowShapes,
                             nodeZLocation
                             # edgeWidth
                        ))
RCy3::setVisualStyle(style)
#plot---------------------
RCy3::setNodeColorMapping(mapping.type = 'd',
                          table.column = 'group',
                          table.column.values = nodes$group %>%
                              unique(),
                          c("#D64358","#4EB3D3","#005A32",
                            "#D64358","#4EB3D3","#005A32"),
                          style.name = style)


# Fig2-A2-Correlation counts----------------------
out.cor.result %>%
    dplyr::mutate(
        group1 = case_when(
            var1 %in% metal_1 ~ "Metal.Serum.1",
            var1 %in% metal_2 ~ "Metal.Serum.2",
            var1 %in% metal_f ~ "Metal.FollicularFluid",
            var1 %in% pfas_1 ~ "PFAS.Serum.1",
            var1 %in% pfas_2 ~ "PFAS.Serum.2",
            var1 %in% metal_h ~ "Metal.Hair"), 
        group2 = case_when(
            var2 %in% metal_1 ~ "Metal.Serum.1",
            var2 %in% metal_2 ~ "Metal.Serum.2",
            var2 %in% metal_f ~ "Metal.FollicularFluid",
            var2 %in% pfas_1 ~ "PFAS.Serum.1",
            var2 %in% pfas_2 ~ "PFAS.Serum.2",
            var2 %in% metal_h ~ "Metal.Hair")) %>% 
    dplyr::filter(p<0.05) -> plotdata

#glance
histogram(plotdata %>%  .$cor %>% abs())
# Fig-2-A2 data
plotdata %>%
    dplyr::filter(group1=="Metal.Serum.1"|group2=="Metal.Serum.1") %>% 
    dplyr::mutate(cor = abs(cor)) %>% 
    dplyr::mutate(corsize = case_when(
        cor <= 0.3 ~ "0-0.3",
        cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
        cor > 0.6 ~ ">0.6"
    )) %>% 
    count(by = corsize) %>% 
    dplyr::mutate(group = "Metal.Serum #1") %>% 
    rbind(
        plotdata %>% 
            dplyr::filter(group1=="Metal.Serum.2"|group2=="Metal.Serum.2") %>% 
            dplyr::mutate(cor = abs(cor)) %>% 
            dplyr::mutate(corsize = case_when(
                cor <= 0.3 ~ "0-0.3",
                cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
                cor > 0.6 ~ ">0.6"
            )) %>% 
            count(by = corsize) %>% 
            dplyr::mutate(group = "Metal.Serum #2") 
    ) %>% 
    rbind(
        plotdata %>% 
            dplyr::filter(group1=="Metal.FollicularFluid"|group2=="Metal.FollicularFluid") %>% 
            dplyr::mutate(cor = abs(cor)) %>% 
            dplyr::mutate(corsize = case_when(
                cor <= 0.3 ~ "0-0.3",
                cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
                cor > 0.6 ~ ">0.6"
            )) %>% 
            count(by = corsize) %>% 
            dplyr::mutate(group = "Metal.FF") 
    )  %>% 
    rbind(
        plotdata %>% 
            dplyr::filter(group1=="Metal.Hair"|group2=="Metal.Hair") %>% 
            dplyr::mutate(cor = abs(cor)) %>% 
            dplyr::mutate(corsize = case_when(
                cor <= 0.3 ~ "0-0.3",
                cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
                cor > 0.6 ~ ">0.6"
            )) %>% 
            count(by = corsize) %>% 
            dplyr::mutate(group = "Metal.Hair") 
    ) %>% 
    rbind(
        plotdata %>% 
            dplyr::filter(group1=="PFAS.Serum.1"|group2=="PFAS.Serum.1") %>% 
            dplyr::mutate(cor = abs(cor)) %>% 
            dplyr::mutate(corsize = case_when(
                cor <= 0.3 ~ "0-0.3",
                cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
                cor > 0.6 ~ ">0.6"
            )) %>% 
            count(by = corsize) %>% 
            dplyr::mutate(group = "PFAS.Serum #1") 
    )  %>% 
    rbind(
        plotdata %>% 
            dplyr::filter(group1=="PFAS.Serum.2"|group2=="PFAS.Serum.2") %>% 
            dplyr::mutate(cor = abs(cor)) %>% 
            dplyr::mutate(corsize = case_when(
                cor <= 0.3 ~ "0-0.3",
                cor > 0.3 & cor <= 0.6 ~ "0.3-0.6",
                cor > 0.6 ~ ">0.6"
            )) %>% 
            count(by = corsize) %>% 
            dplyr::mutate(group = "PFAS.Serum #2") 
    )  -> FigA2
FigA2$by = factor(FigA2$by, levels = c("0-0.3","0.3-0.6",">0.6"))
fills = c("Metal.Serum #1"   = "#FFFFCC",
          "Metal.Serum #2"   = "#CCEBC5",
          "Metal.FF"         = "#C6DBEF",
          "Metal.Hair"       = "#9ECAE1",
          "PFAS.Serum #1"    = "#FC9272",
          "PFAS.Serum #2"    = "#78C679",
          "Follicular Fluid" = "#FEE0D2")
FigA2 %>% 
    arrange(by, -n) %>% 
    dplyr::mutate(order = row_number()) %>% 
    ggplot(aes(x = -order, y = n, fill = group)) +
    geom_bar(stat = "identity") + 
    facet_grid(by~., scale = "free_y") +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(limits = c(0, 850),
                       breaks = c(seq(0, 800, 200)),
                       expand=expansion(add = c(3, 3))) +
    scale_fill_manual(values = fills) + 
    theme_bw() + 
    coord_flip() +
    theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank())


ggsave(str_c(getwd(), "/output/Fig2-A2.png"),
       height = 10,
       width = 12, 
       dpi  =600)




# 提取前十个chemical和Go --------------------------------------------------------
Nodes_GO %>%
    dplyr::filter(!str_detect(id, '^GO|^EX')) %>%
    dplyr::arrange(desc(n)) %>%
    slice(1:10) %>%
    mutate(n_2 = 20:11) %>%
    bind_rows(Nodes_GO %>%
                  dplyr::filter(str_detect(id, '^GO')) %>%
                  dplyr::arrange(desc(n)) %>%
                  slice(1:10) %>%
                  mutate(n_2 = 10:1)) %>%
    mutate(x = 1)-> df_plot

ggplot2::ggplot(data = df_plot,
                aes(x = x, y = n_2,
                    fill = group))+
    ggplot2::geom_point(aes(shape = group),size = 5,color = "grey50",show.legend = F)+
    ggplot2::scale_shape_manual(values = Legend_Shape[Node_Group])+
    ggplot2::scale_fill_manual(values = Legend_Color[Node_Group])+
    ggplot2::geom_text(data = df_plot, aes(label = id),nudge_y = 0,nudge_x = 0.1)+
    ggplot2::geom_text(data = df_plot, aes(label = n),nudge_y = 0,nudge_x = 0.2)+
    xlab(NULL)+
    ylab(NULL)+
    ggplot2::theme_bw()+
    ggplot2::theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()
        
    )  -> count_plot
ggplot2::ggsave("output/go_count_plot.png",
                count_plot,
                #limitsize = FALSE,
                width = 2,
                height = 12,
                dpi = 800
)

# add legend into cytoscape -----------------------------------------------
RCy3::deleteAnnotation(names = "legend1")

RCy3::addAnnotationImage(url = "output/GOLegend.png",
                         network = "GO",
                         x.pos = 900,
                         y.pos = 600,
                         name = "legend1",
                         height = 500,
                         width = 400,
                         canvas = "background",
                         borderColor = "#FFFFFF")

RCy3::deleteAnnotation(names = "legend")
RCy3::addAnnotationImage(url = "output/go_count_plot.png",
                         network = "GO",
                         x.pos = -1500,
                         y.pos = -1000,
                         name = "legend",
                         height = 2000,
                         width = 300,
                         canvas = "background",
                         borderColor = "#FFFFFF")


#RCy3::layoutNetwork(Layout)
RCy3::saveSession(stringr::str_c("output/go_1028"),
                  overwriteFile = T)
RCy3::exportImage(stringr::str_c("output/go_1028"),
                  overwriteFile = T,
                  type = "PNG",
                  zoom=500)
RCy3::exportImage(stringr::str_c("output/go_0830"),
                  overwriteFile = T,
                  type = "SVG",
                  zoom=300)
