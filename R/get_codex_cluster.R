.libPaths(c(.libPaths(), "/cluster/home/yliang_jh/sbin/R/library/4.3.0"))
pkgs <- c("fs", "configr", "stringr", 
          "jhtools", "glue", "patchwork", "tidyverse", "dplyr", "Seurat", "magrittr", 
          "readxl", "writexl", "ComplexHeatmap", 
          "data.table", "ggplot2", "ggbeeswarm", "ggdendro", "dendextend", "deldir",
          "sf", "corrplot", "ggpubr", "survival", "survminer", "forestmodel", "BiocParallel", "BiocNeighbors")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
source("./R/func_utils.R")

knitr::opts_chunk$set(message = FALSE)

rds_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/codex"
fig_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/figs/fig2"


# colors setting
config_fn = "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/meta/manuscript/configs/colors.yaml"
config_fn = "/cluster/home/lixiyue_jh/projects/quarto/stomatology_lvjiong/quarto_book/colors.yaml"
config_list <- show_me_the_colors(config_fn, "all")
colors_celltype <- config_list$cell_type

config <- read.config(config_fn)
cell_type_order <- config$cell_type_order

sampleinfo <- readRDS("/cluster/home/jhuang/projects/stomatology/docs/lvjiong/sampleinfo/sampleinfo.rds")
sampleinfo <- readRDS("/cluster/home/lixiyue_jh/projects/quarto/stomatology_lvjiong/quarto_book/doc/sampleinfo.rds")

sampleinfo$codex <- sampleinfo$codex %>% mutate(
    Time = case_when(Time > 60 ~ 60, TRUE ~ Time),
    Status = case_when(Time > 60 ~ 0,TRUE ~ Status)
  ) %>% as.data.frame()


srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
spinfo <- sampleinfo$codex %>% filter(Type == "Tumor")

df <- srt[[c("orig.ident", "celltype")]] %>% count(orig.ident, celltype, name = "n") %>% group_by(orig.ident) %>%
  mutate(n_total = sum(n), freq = n / n_total) %>% ungroup() %>% mutate(sample_id = orig.ident) %>%
  filter(sample_id %in% spinfo$sample_id)

mat <- df %>% mutate(celltype = factor(celltype, levels = intersect(config$cell_type_order, unique(celltype)))) %>%
  select(sample_id, celltype, freq) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()

dist_mat <- dist(mat, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
cl <- cutree(hc, k = 4)
col_dend <- as.dendrogram(hc)
color_use <- config_list$cluster[unique(cl[labels(col_dend)])]

col_dend <- color_branches(col_dend, k = 4, col = color_use)


df$cluster <- NA
df$cluster <- cl[match(df$sample_id, names(cl))]
df <- df %>% mutate(cluster = factor(cluster, levels = unique(cl)))

df_clinical_cluster <- spinfo %>% filter(Type == "Tumor") %>% filter(sample_id %in% rownames(mat)) %>% 
        mutate(`Age level` = case_when(Age>30 & Age<=50 ~ "30-50", 
                               Age>50 & Age<=70 ~ "50-70",
                               Age>70 & Age<=90 ~ "70-90")) %>%
        select(sample_id, Gender, `Age level`, `Tumor site`, `Diff. level`, Metastasis,`Clinical stage`, `T stage`) %>%
        mutate(across(c("Clinical stage", "Metastasis", "T stage"), as.character)) %>%
        mutate(`Diff. level` = factor(.[["Diff. level"]], levels = c("high", "median", "low")))
df_clinical_cluster$cluster <- cl[match(df_clinical_cluster$sample_id, names(cl))]
df_clinical_cluster <- df_clinical_cluster %>% mutate(cluster = factor(cluster, levels = unique(cl)))

df1 <- df_clinical_cluster %>% dplyr::select(sample_id, cluster) %>% rename(cell_type_cluster = cluster)


spinfo <- sampleinfo$codex
df <- srt[[c("orig.ident", "celltype_cn")]] %>%
  na.omit() %>%
  rename(sample_id = orig.ident) %>%
  count(sample_id, celltype_cn, name = "n") %>%
  tidyr::complete(sample_id, celltype_cn, fill = list(n = 0)) %>%
  group_by(sample_id) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  left_join(spinfo)
mat <- tidyr::pivot_wider(df, id_cols = sample_id, names_from = celltype_cn, values_from = freq) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()
mat <- mat * 100

hc <- hclust(dist(mat))
cl <- cutree(hc, k = 3)
dend <- as.dendrogram(hc)

df2 <- data.frame(cluster = as.character(cl), row.names = as.character(names(cl)), sample_id = as.character(names(cl))) %>%
      mutate(cluster = factor(cluster, levels = unique(cl))) %>% rename(celltype_CN_cluster = cluster)


spinfo <- sampleinfo$codex
df <- srt[[c("orig.ident", "subtype_cn")]] %>%
  na.omit() %>%
  rename(sample_id = orig.ident) %>%
  count(sample_id, subtype_cn, name = "n") %>%
  tidyr::complete(sample_id, subtype_cn, fill = list(n = 0)) %>%
  group_by(sample_id) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  left_join(spinfo)
mat <- tidyr::pivot_wider(df, id_cols = sample_id, names_from = subtype_cn, values_from = freq) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix() %>% t()
mat <- mat * 100
mat_scale <- t(scale(t(mat)))
hc <- hclust(dist(t(mat_scale)))
k = 10
cl <- cutree(hc, k = k)
dend <- as.dendrogram(hc)
df3 <- data.frame(cluster = cl,row.names = names(cl), sample_id = as.character(names(cl))) %>% 
  mutate(cluster = factor(cluster, levels = sort(as.numeric(unique(cl)))))  %>% rename(subtype_CN_cluster = cluster)





lst_pci_message <- readRDS(glue("{rds_dir}/pci_message.rds"))
df_res2 <- lst_pci_message$df_res
df_surv2 <- lst_pci_message$df_surv

df4 <- df_surv2 %>% dplyr::select(sample_id)
pair_fig = list()
pair2 <- c("Macrophage_CD44+-Tumor_CD44+", "CD4T_CD44+-Tumor_CD44+", "CD8T_CD44+-Tumor_CD44+")
#pair2 <- df_res2 %>% filter(pval < 0.05) %>% pull(var)
for(pair in pair2){
  cut <- surv_cutpoint(df_surv2, time = "Time", event = "Status", variables = pair)
  cp <- cut$cutpoint$cutpoint
  dt_grp <- cut$data %>% mutate(cutoff = cp) %>% as.data.frame()
  dt_grp$group <- ifelse(dt_grp[, pair] > cp, "High", "Low")
  dt_grp$sample_id <- df_surv2$sample_id

  df4 <- left_join(df4, dt_grp %>% dplyr::select(sample_id, group, cutoff), by = "sample_id")
  names(df4)[names(df4) == "group"] <- paste0(pair, "_group")
  names(df4)[names(df4) == "cutoff"] <- paste0(pair, "_cutoff")
}


df_list <- list(df1, df2, df3, df4)
df <- reduce(df_list, full_join, by = "sample_id")

saveRDS(df, "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/quarto/docs/codex_group.rds")
