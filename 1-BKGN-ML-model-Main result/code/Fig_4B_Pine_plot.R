rm(list=ls())
library(ggplot2)
library(tidyverse)
library(readxl)

# read data--------
work.dir =str_c(getwd())
#import
df1 <- readxl::read_xlsx(str_c(work.dir,"/@-Knowledge-driven-model/output/1-No.Biolink/0_Final_Model/Features_Selection/importance.rf.xlsx" ))
df2 <- readxl::read_xlsx(str_c(work.dir,"/@-Knowledge-driven-model/output/2-GO.Addition/0_Final_Model/Features_Selection/importance.rf.xlsx" ))
df3 <- readxl::read_xlsx(str_c(work.dir,"/@-Knowledge-driven-model/output/3-Protein.Addition/0_Final_Model/Features_Selection/importance.rf.xlsx" ))

#Nobiolink
df1 %>% 
    dplyr::filter(importance > 0) %>% 
    dplyr::arrange(desc(importance)) %>% 
    dplyr::select(VarLabel, OmicGroup, importance) %>% 
    dplyr::rename(Importance_NoBio = "importance") -> tmp1
#GO Addition
df2 %>% 
    dplyr::filter(importance > 0) %>% 
    dplyr::arrange(desc(importance)) %>% 
    dplyr::select(VarLabel, OmicGroup, importance) %>% 
    dplyr::rename(Importance_GoAdd = "importance") -> tmp2
#Protein Addition
df3 %>% 
    dplyr::filter(importance > 0) %>% 
    dplyr::arrange(desc(importance)) %>% 
    dplyr::select(VarLabel, OmicGroup, importance) %>% 
    dplyr::rename(Importance_ProteinAdd = "importance") -> tmp3
tmp1 %>% 
    full_join(tmp2, by = c("VarLabel","OmicGroup" )) %>%
    full_join(tmp3, by = c("VarLabel","OmicGroup" )) -> features_rf
features_rf %>% 
    print(n = 80)

#plot
library(RColorBrewer)
library(ggsci)
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
scale_color_manual(values = getPalette(26)) #要18个离散颜色

# single pine plot
features_rf  %>% 
    dplyr::select(VarLabel, OmicGroup, Importance_NoBio) %>% 
    group_by(OmicGroup) %>% 
    dplyr::arrange(OmicGroup, desc(Importance_NoBio)) %>% 
    dplyr::filter(!is.na(Importance_NoBio)) %>% 
    ggplot(aes(x="", y = Importance_NoBio, fill = reorder(VarLabel, -Importance_NoBio))) +
    geom_bar(width = 0.01, stat = "identity", color = "black") + 
    coord_polar("y", start = 0) + 
    theme_void() +
    scale_fill_futurama()


# ggsave(str_c(getwd(),"/@IVF305_rf/output/Feature_importance_NoBiolink.png" ))

# facet
features_rf  %>% 
    gather(key = "Group", value = "Importance", -c("VarLabel", "OmicGroup"))  %>% 
    group_by(Group) %>% 
    dplyr::arrange(Group, desc(Importance)) %>% 
    dplyr::mutate(Group = str_replace(Group,"Importance_NoBio","NoBiolink"),
                  Group = str_replace(Group,"Importance_GoAdd","GoAddition"),
                  Group = str_replace(Group,"Importance_ProteinAdd","ProteinAddition")) %>% 
    dplyr::filter(!is.na(Importance)) %>% 
    group_by(Group) %>% 
    dplyr::mutate(importance_perct = map_dbl(Importance,
                                    ~(./sum(Importance)))) -> a
a$Group = factor(a$Group, ordered = T, levels = c("NoBiolink",
                                                    "GoAddition",
                                                    "ProteinAddition"))
a %>% 
    ggplot(aes(x="", y = importance_perct, fill = reorder(VarLabel, -importance_perct))) +
    geom_bar(width = 0.1, stat = "identity", color = "black") + 
    coord_polar("y", start = 0) + 
    facet_grid(.~Group)+
    theme_void()+
    scale_fill_d3("category20")
ggsave(str_c(work.dir,"/@-Knowledge-driven-model/output/Feature_importance_total_RF.png"))

a %>% 
    ggplot(aes(x="", y = importance_perct, fill = reorder(VarLabel, -importance_perct))) +
    geom_bar(width = 0.1, stat = "identity", color = "black") + 
    coord_polar("y", start = 0) + 
    facet_grid(.~Group)+
    theme_void()+
    scale_fill_d3("category20") +
    theme(legend.position="none")
ggsave(str_c(work.dir,"/@-Knowledge-driven-model/output/Feature_importance_total_RF_nolegend.png" ))

# export GO and Protein list
features_rf %>% 
    dplyr::filter(OmicGroup %in% c("GO1", "GO2", "GOff", "GOhr",
                                   "EXP1", "EXP2", "EXPff", "EXPhr")) %>% 
    writexl::write_xlsx(str_c(work.dir,"/@-Knowledge-driven-model/output/Selected_GOandProtein_RF.xlsx" ))

