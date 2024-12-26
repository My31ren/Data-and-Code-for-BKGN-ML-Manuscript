rm(list=ls())
pacman::p_load("tidyverse",
               "mlr3verse",
               "rlist",
               "readxl")
# create folders
if(!file.exists(str_c(getwd(),"/input"))){
  dir.create(str_c(getwd(),"/input"))
}
if(!file.exists(str_c(getwd(),"/output"))){
  dir.create(str_c(getwd(),"/output"))
}
# define paths
inpath = str_c(getwd(),"/input")
outpath = str_c(getwd(),"/output")

# 1-Protein: load data ---------------------------------------------------------------
df = readr::read_csv(str_c(inpath, "/TableS2_exe_protein.csv"))

# senario 1: full metal(loid)s forms are considered
df %>%
    dplyr::select(group, group_name, source, target, preferred.name.x) %>% 
    unique() %>%
    dplyr::rename(EXP = preferred.name.x,
                  name = group_name) %>% 
    dplyr::select(name, EXP) %>%
    unique() %>%
    dplyr::mutate(index = 1) %>%
    pivot_wider(names_from = EXP,
                names_sep = "_",
                values_from = index) %>%
    dplyr::select(name, everything())-> index_matrix2
# replace NA as 0
index_matrix2[is.na(index_matrix2)] =0
# save
save(index_matrix2, file = str_c(outpath, "/EXP_matrix_full.rdata"))


# 2-GO: load data ---------------------------------------------------------------
df2 = readr::read_csv(str_c(inpath, "/TableS3_exe_go.csv"))

# senario 1: full metal(loid)s forms are considered 
df2 %>%
  dplyr::select(group_name, target) %>%
  dplyr::mutate(GO =str_replace(.$target, ":", "_"),
                name = group_name) %>%
  dplyr::select(name, GO) %>%
  unique() %>% #1849
  dplyr::mutate(index = 1) %>%
  pivot_wider(names_from = GO,
              names_sep = "_",
              values_from = index) %>%
  dplyr::select(name, everything())-> index_matrix4
# replace NA as 0
index_matrix4[is.na(index_matrix4)] =0
# save
save(index_matrix4, file = str_c(outpath,"/GO_matrix_full.rdata"))

