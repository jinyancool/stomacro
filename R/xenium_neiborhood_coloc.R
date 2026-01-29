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

param <- MulticoreParam(workers = 12, progressbar = TRUE)



##  celltype neighborhood
df <- lapply(c("0066253", "0066266"), function(slide){
  f <- glue("{rds_dir}/spe_{slide}.rds")
  spe <- readRDS(f)
  spe_s <- spe[, !is.na(spe$cell_type) & !is.na(spe$`Diff. level`)]
  spe_s <- buildSpatialGraph(spe_s, 
                             coords = spatialCoordsNames(spe_s),
                             img_id = "roi", type = "knn", k = 10)
  spe_s <- aggregateNeighbors(spe_s, 
                              colPairName = "knn_interaction_graph", 
                              aggregate_by = "metadata", count_by = "cell_type")
  df_n <- spe_s$aggregatedNeighbors
  df_n$cell_id <- paste0(colnames(spe_s), "_", spe_s$roi)
  df_n <- df_n %>% as.data.frame(check.names = FALSE)
  return(df_n)
}) %>% rbindlist() %>% as.data.frame(check.names = FALSE)
rownames(df) <- df$cell_id
df <- df %>% dplyr::select(-cell_id)

k_res <- kmeans(df, centers = 10)
dt_cn <- data.table(cell_id = rownames(df),
                    celltype_cn = paste0("CN", k_res$cluster) %>% 
                    factor(levels = paste0("CN", 1:10)))

dt_anno <- readRDS(glue("{rds_dir}/celltyped_meta.rds")) 
dt_anno <- dt_anno %>% filter(!is.na(cell_type) & !is.na(`Diff. level`)) %>% left_join(dt_cn, by = "cell_id")
mat_count <- dt_anno %>%
  dplyr::select(cell_type, celltype_cn) %>%
  as.data.table() %>%
  dcast(celltype_cn ~ cell_type, fun.aggregate = length) %>%
  tibble::column_to_rownames("celltype_cn")
mat <- mat_count/rowSums(mat_count)

neighborhood_list <- list("mat" = mat, "metadata" = dt_anno)
saveRDS(neighborhood_list, glue("{rds_dir}/celltype_neiborhood.rds"))


## celltype neighborhood: Epithelial

df <- lapply(c("0066253", "0066266"), function(slide){
  f <- glue("{rds_dir}/spe_{slide}.rds")
  spe <- readRDS(f)
  spe_s <- spe[, !is.na(spe$Epithelial.sub.supply) & !is.na(spe$`Diff. level`)]
  spe_s <- buildSpatialGraph(spe_s, 
                             coords = spatialCoordsNames(spe_s),
                             img_id = "roi", type = "knn", k = 10)
  spe_s <- aggregateNeighbors(spe_s, 
                              colPairName = "knn_interaction_graph", 
                              aggregate_by = "metadata", count_by = "Epithelial.sub.supply")
  df_n <- spe_s$aggregatedNeighbors
  df_n$cell_id <- paste0(colnames(spe_s), "_", spe_s$roi)
  df_n <- df_n %>% as.data.frame()
  return(df_n)
}) %>% rbindlist() %>% as.data.frame()
rownames(df) <- df$cell_id

dt_anno <- readRDS(glue("{rds_dir}/celltyped_meta.rds"))
epi_cell <- dt_anno %>% filter(!is.na(Epithelial.1.subtype) & !is.na(`Diff. level`)) %>% pull(cell_id)
df <- df[df$cell_id %in% epi_cell, ]

df <- df %>% dplyr::select(-cell_id)
k_res <- kmeans(df, centers = 10)
dt_cn <- data.table(cell_id = rownames(df),
                    subtype_cn = paste0("CN", k_res$cluster) %>% 
                    factor(levels = paste0("CN", 1:10)))
dt_anno <- readRDS(glue("{rds_dir}/celltyped_meta.rds")) 
dt_anno <- dt_anno %>% filter(!is.na(Epithelial.1.subtype) & !is.na(`Diff. level`)) %>% left_join(dt_cn, by = "cell_id")
mat_count <- dt_anno %>%
  dplyr::select(Epithelial.1.subtype, subtype_cn) %>%
  as.data.table() %>%
  dcast(subtype_cn ~ Epithelial.1.subtype, fun.aggregate = length) %>%
  tibble::column_to_rownames("subtype_cn")
mat <- mat_count/rowSums(mat_count)

neighborhood_epi_list <- list("mat" = mat, "metadata" = dt_anno)
saveRDS(neighborhood_epi_list, glue("{rds_dir}/celltype_neiborhood_epithelial.rds"))







## co-localization

spe <- readRDS(glue("{rds_dir}/spe_0066253.rds"))
grp_lst <- bplapply(unique(spe$sample_id), function(x) {
  spe_s <- spe[, spe$roi == x]
  sqe <- readHoodData(spe_s, anno_col = "cell_type")
  nbs <- findNearCells(sqe, k = 5)
  mtx <- scanHoods(nbs$distance)      
  df_grp <- mergeByGroup(mtx, nbs$cells) %>% as.data.frame()
  return(df_grp)
}, BPPARAM = param)
grp <- rbindlist(grp_lst, fill = TRUE) %>% 
  as.matrix()
grp[is.na(grp)] <- 0
rownames(grp) <- colnames(spe)
sqe1 <- readHoodData(spe, anno_col = "cell_type")
sqe1 <- mergeHoodSpe(sqe1, grp)  


spe <- readRDS(glue("{rds_dir}/spe_0066266.rds"))
spe <- spe[, !spe$sample_id %in% c("F4-2", "F14")]
grp_lst <- bplapply(unique(spe$sample_id), function(x) {
  spe_s <- spe[, spe$roi == x]
  sqe <- readHoodData(spe_s, anno_col = "cell_type")
  nbs <- findNearCells(sqe, k = 5)
  mtx <- scanHoods(nbs$distance)      
  df_grp <- mergeByGroup(mtx, nbs$cells) %>%
    as.data.frame()
  return(df_grp)
}, BPPARAM = param)
grp <- rbindlist(grp_lst, fill = TRUE) %>% 
  as.matrix()
grp[is.na(grp)] <- 0
rownames(grp) <- colnames(spe)
sqe2 <- readHoodData(spe, anno_col = "cell_type")
sqe2 <- mergeHoodSpe(sqe2, grp)

sqe <- cbind(sqe1, sqe2)
saveRDS(sqe, glue("{rds_dir}/xenium_sqe_coloc.rds"))


sqe <- readRDS(glue("{rds_dir}/xenium_sqe_coloc.rds"))
dt_celltype <- readRDS(glue("{rds_dir}/celltype_anno.rds"))
celltypes <- dt_celltype %>% filter(!is.na(cell_type)) %>% pull(cell_type) %>% unique()
#colnames(grp)

cor_lst <- list()
cor_lst[["col_matrix_total"]] <- plotColocal(sqe, pm_cols = colnames(grp), return_matrix = TRUE)
for (grp_diff in unique(sqe$`Diff. level`)){
  sqe_s <- sqe[, sqe$`Diff. level` == grp_diff]
  dt_cor_diff <- plotColocal(sqe_s, pm_cols = colnames(grp), return_matrix = TRUE)
  cor_lst[[glue("cor_matrix_diff_{grp_diff}")]] <- dt_cor_diff
}
dt_cor_sample <- lapply(unique(sqe$sample_id), function(sp){
  sqe_s <- sqe[, sqe$sample_id == sp]
  cor <- plotColocal(sqe_s, pm_cols = colnames(grp), return_matrix = TRUE)
  dt_cor_s <- as.data.table(cor, keep.rownames = "from") %>%
    melt(id.vars = "from", variable.name = "to", value.name = "prob")
  dt_cor_s$sample_id <- sp
  return(dt_cor_s)
}) %>% rbindlist()

cor_lst[["cor_prob_sampleid"]] <- dt_cor_sample

dt_cor_sample[, from := as.factor(from)]
dt_cor_sample_simplify <- dt_cor_sample[as.numeric(from) > as.numeric(to), `:=`(from = to, to = from)] %>% unique
cor_lst[["cor_prob_sampleid_simplify"]] <- dt_cor_sample_simplify

saveRDS(cor_lst, glue("{rds_dir}/xenium_sqe_cor.rds"))





