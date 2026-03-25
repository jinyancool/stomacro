pkgs <- c("fs", "configr", "stringr", "scider",
          "jhtools", "glue", "patchwork", "tidyverse", "dplyr", "Seurat", "magrittr", "SeuratDisk", "readr",
          "readxl", "writexl", "openxlsx", "ComplexHeatmap", "SpatialExperiment", "SpatialExperimentIO", "imcRtools",
          "data.table", "ggplot2", "ggtext", "ggbeeswarm", "ggdendro", "dendextend", "deldir", "BPCells", "hoodscanR",
          "sf", "corrplot", "ggpubr", "BiocParallel", "clusterProfiler", "KEGG.db", "decoupleR",
          "scater", "scuttle", "CCPlotR", "ggspavis", "reticulate", "zellkonverter")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



set.seed(12387)
rds_dir = "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/xenium"
fig_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/figs/fig4"

# colors setting
config_fn = "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/meta/manuscript/configs/colors.yaml"
config_list <- show_me_the_colors(config_fn, "all")
colors_celltype <- config_list$cell_type

config <- read.config(config_fn)
cell_type_order <- config$cell_type_order

sampleinfo <- readRDS("/cluster/home/jhuang/projects/stomatology/docs/lvjiong/sampleinfo/sampleinfo.rds")


use_python("/cluster/home/lixiyue_jh/.conda/envs/pyr/bin/python")
param <- MulticoreParam(workers = 12, progressbar = TRUE)

run_pathway <- function(sample, rds_dir){
    sce <- readRDS(glue("{rds_dir}/sce_{sample}.rds"))
    sce$ccc <- coalesce(sce$Epithelial.1.subtype, sce$cell_type)
    mu <- aggregateAcrossCells(sce, ids=sce$ccc, statistics="mean", use.assay.type=assayNames(sce))
  
    ks <- expand.grid(ks <- colnames(mu), ks)
    cs <- split(colnames(sce), sce$ccc)
    lr <- grepl("-", rownames(sce))
    lr <- lr & !grepl("total", rownames(sce))
    ss <- strsplit(rownames(sce), "-")
    l <- sapply(ss, .subset, 1)
    r <- sapply(ss, .subset, 2)
    df <- mapply(
        i=ks[, 1],
        j=ks[, 2],
        SIMPLIFY=FALSE,
        FUN = function(i, j) {
            if (!i %in% names(cs) || !j %in% names(cs)) return(NULL)
            source <- sce[lr, cs[[i]]]
            target <- sce[lr, cs[[j]]]
            if (ncol(source) == 0 || ncol(target) == 0) return(NULL)
            sr <- cbind(
                rowMeans(assay(source, "s")),
                rowMeans(assay(target, "r")))
            data.frame(
                source=paste(i), target=paste(j),
                ligand=l[lr], receptor=r[lr],
                score=rowMeans(sr))
        }) |> do.call(what=rbind)
    df <- df %>% as.data.frame() %>% mutate(sample_id = sample) %>% as.data.frame()
    rownames(df) <- NULL
    
    return(list(mu = mu, df = df))
}



run_path_test <- function(df, colUse, groups, commu_groups){
    cc_path <- data.frame(path = character(), cell_type = character(), statistic = numeric(), p_value = numeric(),stringsAsFactors = FALSE)
    for (p in unique(df$path)){
        path_data <- df[df$path == p, ]
        colForce <- setdiff(c("source", "target"), c(colUse))[[1]]
        celltypes <- unique(path_data[[colForce]]) %>% .[. %in% commu_groups]
        for (celltype in celltypes){
            df_p <- path_data %>% filter(.data[[colUse]] %in% groups, .data[[colForce]] %in% celltype)
            mean_df <- aggregate(score ~ ., data = df_p[,c("score", colUse)], FUN = mean, na.rm = TRUE)
            mean_vec <- setNames(mean_df$score, paste0("mean_", mean_df[[colUse]]))
            test_result <- tryCatch({
                kruskal.test(reformulate(colUse, response = "score"), data = df_p)
            }, error = function(e) {
                return(NULL)
            })
            if (!is.null(test_result)){
                cc_path_add <- cbind(data.frame(path = p, cell_type = celltype, statistic = test_result$statistic, p_value = test_result$p.value, stringsAsFactors = FALSE), 
                    as.data.frame(as.list(mean_vec)))
                cc_path <- rbind(cc_path, cc_path_add)
            }
        }
    }
    return(cc_path)
}


run_path_test_other <- function(df, colUse, case_groups, ctrl_groups, commu_groups, case_name, ctrl_name){
    cc_path <- data.frame(path = character(), cell_type = character(), 
        case = character(), mean_case = character(), ctrl = character(), mean_ctrl = character(),
        t_statistic = numeric(), df = numeric(), p_value = numeric(),stringsAsFactors = FALSE)
    for (p in unique(df$path)){
        path_data <- df[df$path == p, ]
        colForce <- setdiff(c("source", "target"), c(colUse))[[1]]
        celltypes <- unique(path_data[[colForce]]) %>% .[. %in% commu_groups]
        for (celltype in celltypes){
            case_data <- path_data[(path_data[[colUse]] %in% case_groups)&(path_data[[colForce]] %in% c(celltype)),][["score"]] %>% na.omit()
            ctrl_data <- path_data[(path_data[[colUse]] %in% ctrl_groups)&(path_data[[colForce]] %in% c(celltype)),][["score"]] %>% na.omit()
            test_result <- tryCatch({
                t.test(case_data, ctrl_data, var.equal = FALSE)
            }, error = function(e) {
                return(NULL)
            })
            if (!is.null(test_result)){
                mean1 <- mean(case_data, na.rm = TRUE)
                mean2 <- mean(ctrl_data, na.rm = TRUE)
                cc_path_add <- data.frame(
                    path = p, cell_type = celltype, case = case_name, mean_case = mean1, ctrl = ctrl_name, mean_ctrl = mean2,
                    t_statistic = test_result$statistic, df = test_result$parameter, p_value = test_result$p.value, stringsAsFactors = FALSE)
                cc_path <- rbind(cc_path, cc_path_add)
            }
        }
    }
    return(cc_path)
}







samples_id <- sampleinfo$xenium %>% pull(sample_id) %>% unique() %>% .[! . %in% c("F4-2", "F14")]

mu_list <- list()
df_list <- list()
rds_h_dir <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/xenium/ccc"
res_list <- bplapply(samples_id, run_pathway, rds_dir = rds_h_dir, BPPARAM = param)
names(res_list) <- samples_id

mu_list <- lapply(res_list, `[[`, "mu")
df_list <- lapply(res_list, `[[`, "df")

saveRDS(mu_list, glue("{rds_dir}/xenium_CCC_Epithelial_basal_mu_list.rds"))
saveRDS(df_list, glue("{rds_dir}/xenium_CCC_Epithelial_basal_df_path_list.rds"))

df <- do.call(rbind, df_list)
rownames(df) <- NULL
df <- df %>% mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-"))
df <- df %>% filter(!sample_id %in% c("F14", "F4-2")) %>% as.data.frame()
saveRDS(df, glue("{rds_dir}/xenium_CCC_Epithelial_basal_df_path.rds"))


df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_basal_df_path.rds"))
meta <- readRDS(glue("{rds_dir}/celltype_anno.rds"))
epi_sub <- intersect(unique(df$source), unique(meta$Epithelial.1.subtype))
epi_sub_exclu <- setdiff(epi_sub, c("basal_cell", "basal_inflam", "basal_invasive"))
celltype_exclu_epi <- setdiff(unique(df$source), epi_sub) %>% .[!. %in% c("Glial", "smoothMC", "skeletalMC")]
celltype_exclu_epi <- c("Macrophage", "T cell", "Fibroblast", "Endothelial")

cc_path_s <- run_path_test(df, "source", c("basal_cell", "basal_inflam", "basal_invasive"), celltype_exclu_epi)
cc_path_r <- run_path_test(df, "target", c("basal_cell", "basal_inflam", "basal_invasive"), celltype_exclu_epi)

cc_path_s_filter <- cc_path_s %>% filter(p_value <= 0.05) %>% arrange(desc(abs(statistic)), p_value) %>% as.data.frame()
cc_path_r_filter <- cc_path_r %>% filter(p_value <= 0.05) %>% arrange(desc(abs(statistic)), p_value) %>% as.data.frame()

results_cc_path <- list("source_filter" = cc_path_s_filter, "target_filter" = cc_path_r_filter)
saveRDS(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_basal_test.rds"))
list2excel(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_basal_test.xlsx"))



cc_path_s_basal_cell <- run_path_test_other(df, "source", c("basal_cell"), c("basal_inflam", "basal_invasive"), celltype_exclu_epi, "basal_cell", "other2")
cc_path_r_basal_cell <- run_path_test_other(df, "target", c("basal_cell"), c("basal_inflam", "basal_invasive"), celltype_exclu_epi, "basal_cell", "other2")
cc_path_s_basal_cell <- cc_path_s_basal_cell %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()
cc_path_r_basal_cell <- cc_path_r_basal_cell %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()

cc_path_s_basal_inflam <- run_path_test_other(df, "source", c("basal_inflam"), c("basal_cell", "basal_invasive"), celltype_exclu_epi, "basal_inflam", "other2")
cc_path_r_basal_inflam <- run_path_test_other(df, "target", c("basal_inflam"), c("basal_cell", "basal_invasive"), celltype_exclu_epi, "basal_inflam", "other2")
cc_path_s_basal_inflam <- cc_path_s_basal_inflam %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()
cc_path_r_basal_inflam <- cc_path_r_basal_inflam %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()


cc_path_s_basal_invasive <- run_path_test_other(df, "source", c("basal_invasive"), c("basal_inflam", "basal_cell"), celltype_exclu_epi, "basal_invasive", "other2")
cc_path_r_basal_invasive <- run_path_test_other(df, "target", c("basal_invasive"), c("basal_inflam", "basal_cell"), celltype_exclu_epi, "basal_invasive", "other2")
cc_path_s_basal_invasive <- cc_path_s_basal_invasive %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()
cc_path_r_basal_invasive <- cc_path_r_basal_invasive %>% filter(p_value <= 0.05) %>% arrange(desc((t_statistic)), p_value) %>% as.data.frame()


results_cc_path <- list("source_filter_basal_cell" = cc_path_s_basal_cell, "target_filter_basal_cell" = cc_path_r_basal_cell, 
                        "source_filter_basal_inflam" = cc_path_s_basal_inflam, "target_filter_basal_inflam" = cc_path_r_basal_inflam,
                        "source_filter_basal_invasive" = cc_path_s_basal_invasive, "target_filter_basal_invasive" = cc_path_r_basal_invasive)
saveRDS(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_basal_test_vs_other2.rds"))
list2excel(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_basal_test_vs_other2.xlsx"))




cc_path <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_basal_test.rds"))
cc_path_s_filter <- cc_path$source_filter
cc_path_r_filter <- cc_path$target_filter

select_path_s <- cc_path_s_filter %>% pull(path) %>% head(200) %>% unique()
path1 <- cc_path_s_filter %>% filter(str_detect(path, "^(JAG|LAM|DLL|FGF|TGFB)")) %>% pull(path) %>% head(200) %>% unique()
path2 <- c("THBS1-SDC4", "TNC-SDC4", "ANGPTL4-SDC4","HGF-MET", "COL4A1-SDC4", "COL4A2-SDC4", "COL4A1-CD44", "COL4A2-CD44")
rm_path <- c("VEGFC-FLT4", "VEGFC-KDR", "VEGFC-FLT4_KDR", "DLL1-NOTCH4", "JAG1-NOTCH4", "CXCL2-CXCR1")
select_path_s <- c(path1,path2) %>% .[!. %in% rm_path] %>% unique()
select_path_r <- cc_path_r_filter %>% pull(path) %>% head(200) %>% unique()
path1 <- cc_path_r_filter %>% filter(str_detect(path, "^(JAG|LAM|DLL|FGF|TGFB)")) %>% pull(path) %>% head(200) %>% unique()
path2 <- c("THBS1-SDC4", "TNC-SDC4", "ANGPTL4-SDC4","HGF-MET")
rm_path <- c("EFNB2-EPHB2", "EFNB2-EPHA4", "PDGFA-PDGFRB", "EDN1-EDNRA", "DLL1-NOTCH4", "DLL4-NOTCH3", "JAG1-NOTCH4")
select_path_r <- c(path1,path2) %>% .[!. %in% rm_path] %>% unique()

df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_basal_df_path.rds"))

cell_cell_order <- c(
  paste("basal_cell", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste("basal_inflam", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste("basal_invasive", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_cell", sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_inflam", sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_invasive", sep="-")
)

df_heat_cp_wide_s <- df %>% 
    filter(source %in% c("basal_cell", "basal_inflam", "basal_invasive"), target %in% celltype_exclu_epi, score > 0) %>%
    group_by(source, target, ligand, receptor) %>%
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_s) %>%
    dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, intersect(cell_cell_order,  colnames(.)), drop = FALSE]
df_heat_cp_wide_r <- df %>% 
    filter(target %in% c("basal_cell", "basal_inflam", "basal_invasive"), source %in% celltype_exclu_epi, score > 0) %>%
    group_by(source, target, ligand, receptor) %>%
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_r) %>%
    dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, intersect(cell_cell_order,  colnames(.)), drop = FALSE]

df_heat_cp_wide <- list(s = df_heat_cp_wide_s, r = df_heat_cp_wide_r)
saveRDS(df_heat_cp_wide, glue("{rds_dir}/xenium_CCC_Epithelial_basal_plot.rds"))



cc_path <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_basal_test_vs_other2.rds"))
cc_path_s_filter_basal_cell <- cc_path$source_filter_basal_cell %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()
cc_path_r_filter_basal_cell <- cc_path$target_filter_basal_cell %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()
cc_path_s_filter_basal_inflam <- cc_path$source_filter_basal_inflam %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()
cc_path_r_filter_basal_inflam <- cc_path$target_filter_basal_inflam %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()
cc_path_s_filter_basal_invasive <- cc_path$source_filter_basal_invasive %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()
cc_path_r_filter_basal_invasive <- cc_path$target_filter_basal_invasive %>% arrange(p_value, desc(t_statistic)) %>% as.data.frame()

path1 <- cc_path_s_filter_basal_cell %>% pull(path) %>% head(100) %>% unique() %>% head(30)
path2 <- cc_path_s_filter_basal_inflam %>% pull(path) %>% head(100) %>% unique() %>% head(30)
path3 <- cc_path_s_filter_basal_invasive %>% pull(path) %>% head(100) %>% unique() %>% head(10)
select_path_s <- c(path1,path2,path3) %>% unique()
select_path_s <- c("CDH5-CDH5","KITLG-KIT","JAG1-NOTCH4","NTF3-NTRK2","FGF1-FGFR3","DLL1-NOTCH4",
                   "WNT5A-FZD10","FGF7-FGFR2","FGF2-FGFR2","THBS2-CD36","THBS1-CD47","CXCL11-ACKR3",
                   "LAMA3-CD44","COL4A1-CD44","LAMC2-CD44","TNFSF10-TNFRSF10A","COL4A1-SDC4","COL4A2-SDC4",
                   "THBS1-SDC4","EREG-EGFR","LAMB3-CD44","LAMC2-CD44","TNC-SDC4","LAMA3-CD44","THBS1-SDC4",
                   "COL4A2-CD44","COL4A1-SDC4","TGFA-EGFR","THBS1-CD47","COL4A1-CD44")
rm_path <- c("CDH5-CDH5", "DLL1-NOTCH4", "JAG1-NOTCH4")
select_path_s <- setdiff(select_path_s, rm_path)

path1 <- cc_path_r_filter_basal_cell %>% pull(path) %>% head(100) %>% unique() %>% head(30)
path2 <- cc_path_r_filter_basal_inflam %>% pull(path) %>% head(100) %>% unique() %>% head(30)
path3 <- cc_path_r_filter_basal_invasive %>% pull(path) %>% head(100) %>% unique() %>% head(10)
select_path_r <- c(path1,path2,path3) %>% unique()
select_path_r <- c(
    "CDH5-CDH5","GAS6-TYRO3","FGF1-FGFR3","NTF3-NTRK2","WNT5A-FZD3",
    "JAG1-NOTCH4","KITLG-KIT","PDGFC-PDGFRA","FGF7-FGFR2","DLL1-NOTCH4",
    "THBS1-CD47","CXCL11-ACKR3","LAMA3-CD44","LAMC2-CD44","ITGB2-ICAM1",
    "COL4A2-SDC4","COL4A1-SDC4","COL4A5-SDC4","THBS1-SDC4","TNFSF12-TNFRSF12A",
    "LAMB3-CD44","TNC-SDC4","COL4A1-CD44","COL4A2-CD44","TGFA-EGFR","THBS2-SDC4","ITGB2-ICAM1")
rm_path <- c("CDH5-CDH5", "DLL1-NOTCH4", "JAG1-NOTCH4", "KITLG-KIT", "COL4A1-CD44","COL4A2-CD44", "PDGFC-PDGFRA", "WNT5A-FZD3", "ITGB2-ICAM1", "TNFSF12-TNFRSF12A")
select_path_r <- setdiff(select_path_r, rm_path)

df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_basal_df_path.rds"))

cell_cell_order <- c(
  paste("basal_cell", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste("basal_inflam", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste("basal_invasive", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_cell", sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_inflam", sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "basal_invasive", sep="-")
)



#cell_cell_order <- unlist(lapply(intersect(config$cell_type_order, unique(c(df$source, df$target))), function(ct) {
#  c(
#    paste("basal_cell", ct, sep="-"),
#    paste("basal_inflam", ct, sep="-"),
#    paste("basal_invasive", ct, sep="-"),
#    paste(ct, "basal_cell", sep="-"),
#    paste(ct, "basal_inflam", sep="-"),
#    paste(ct, "basal_invasive", sep="-")
#  )
#}))

df_heat_cp_wide_s <- df %>% 
    filter(source %in% c("basal_cell", "basal_inflam", "basal_invasive"), target %in% celltype_exclu_epi, score > 0) %>%
    group_by(source, target, ligand, receptor) %>%
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_s) %>%
    dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, intersect(cell_cell_order,  colnames(.)), drop = FALSE]
df_heat_cp_wide_r <- df %>% 
    filter(target %in% c("basal_cell", "basal_inflam", "basal_invasive"), source %in% celltype_exclu_epi, score > 0) %>%
    group_by(source, target, ligand, receptor) %>%
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_r) %>%
    dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, intersect(cell_cell_order,  colnames(.)), drop = FALSE]

df_heat_cp_wide <- list(s = df_heat_cp_wide_s, r = df_heat_cp_wide_r)
saveRDS(df_heat_cp_wide, glue("{rds_dir}/xenium_CCC_Epithelial_basal_plot_2.rds"))
list2excel(df_heat_cp_wide, glue("{rds_dir}/xenium_CCC_Epithelial_basal_plot_2.xlsx"), showRow=T)


