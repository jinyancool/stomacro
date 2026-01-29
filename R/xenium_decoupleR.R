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

param <- MulticoreParam(workers = 4, progressbar = TRUE)


## decouple R analyze tumor erosion zone sample

srt <- readRDS(glue("{rds_dir}/srt_erosion_region.rds"))
srt <- subset(srt, subset = sample_id %notin% c("F4-2", "F14"))
net_d <- readRDS(glue("{rds_dir}/decoupleR.rds"))
net_p <- readRDS(glue("{rds_dir}/ROGENy_pathway_from_decoupleR.rds"))
srt$sam_region <- str_c(srt$sample_id, "@", srt$Region)
srt$sam_region <- factor(srt$sam_region, levels = sort(unique(srt$sam_region)))
srt$stype_region <- str_c(srt$subtype, "@", srt$Region)
srt$stype_region <- factor(srt$stype_region, levels = sort(unique(srt$stype_region)))


select_celltypes <- c("Macrophage", "Epithelial")
results <- bplapply(select_celltypes, function(celltype){
    srt_s <- subset(srt, subset = cell_type == celltype)
    DefaultAssay(srt_s) <- "Xenium"
    mat <- GetAssayData(srt_s, assay = "Xenium", layer = "data")
    mat_dgc <- as(mat, "dgCMatrix") %>% as.matrix()

    #tfs
    acts_tar <- decoupleR::run_ulm(mat = mat_dgc, net = net_d, 
                            .source = 'source', .target = 'target',
                            .mor='mor', minsize = 5)
    srt_s[['tfsulm']] <- acts_tar %>%
        tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
        tibble::column_to_rownames('source') %>%
        Seurat::CreateAssayObject(.)
    DefaultAssay(srt_s) <- "tfsulm"
    srt_s <- Seurat::ScaleData(srt_s)
    srt_s@assays$tfsulm@data <- srt_s@assays$tfsulm@scale.data

    Idents(srt_s) <- "Region"
    df_tfs <- t(as.matrix(srt_s@assays$tfsulm@data)) %>%
        as.data.frame(check.names = FALSE) %>%
        dplyr::mutate(Region = srt_s$Region, subtype = srt_s$subtype, samples = srt_s$sample_id,
                    sam_region = str_c(samples, "@", Region), stype_region = str_c(Region, "@", subtype))
    df_tfs_region <- df_tfs %>% dplyr::select(-sam_region, -subtype, -samples, -stype_region) %>%
        tidyr::pivot_longer(cols = -Region, names_to = "source", values_to = "score") %>%
        dplyr::group_by(Region, source) %>%
        dplyr::summarise(mean = mean(score))
    df_tfs_sample_region <- df_tfs %>% dplyr::select(-Region, -subtype, -samples, -stype_region) %>%
        tidyr::pivot_longer(cols = -sam_region, names_to = "source", values_to = "score") %>%
        dplyr::group_by(sam_region, source) %>%
        dplyr::summarise(mean = mean(score))
    df_tfs_stype_region <- df_tfs %>% dplyr::select(-Region, -subtype, -samples, -sam_region) %>%
        tidyr::pivot_longer(cols = -stype_region, names_to = "source", values_to = "score") %>%
        dplyr::group_by(stype_region, source) %>%
        dplyr::summarise(mean = mean(score))
    n_tfs <- 100
    top_region_tfs <- df_tfs_region %>% dplyr::group_by(source) %>%
        dplyr::summarise(std = stats::sd(mean)) %>%
        dplyr::arrange(-abs(std)) %>%
        head(n_tfs) %>%
        dplyr::pull(source)
    mat_top_acts_tfs_region <- df_tfs_region %>% dplyr::filter(source %in% top_region_tfs) %>%
        tidyr::pivot_wider(id_cols = 'Region', names_from = 'source', values_from = 'mean') %>%
        tibble::column_to_rownames('Region') %>%
        as.matrix()
    
    #pathway
    acts_path <- decoupleR::run_ulm(mat = mat_dgc, net = net_p,
                               .source = 'source', .target = 'target',
                               .mor='weight', minsize = 5)
    srt_s[["pathulm"]] <- acts_path %>% 
        tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
        tibble::column_to_rownames('source') %>%
        Seurat::CreateAssayObject(.)
    DefaultAssay(srt_s) <- "pathulm"
    srt_s <- Seurat::ScaleData(srt_s)
    srt_s@assays$pathulm@data <- srt_s@assays$pathulm@scale.data
    
    df_score <- t(srt_s@assays$pathulm@data) %>%
        as_tibble(rownames = "barcode") %>%
        dplyr::mutate(Region = srt_s$Region, sam_region = srt_s$sam_region, subtype = srt_s$subtype) %>%
        tidyr::pivot_longer(cols = Androgen:p53, names_to = "pathway", values_to = "score")
    df_score_sam <- df_score %>%
        dplyr::group_by(sam_region, Region, pathway) %>%
        dplyr::summarise(mean = mean(score)) %>% ungroup()

    decoupleR_message <- list(path = acts_path, tfs = acts_tar, 
                                path_score = df_score, path_score_sam = df_score_sam,
                                tfs_total = df_tfs, tfs_region = df_tfs_region, 
                                tfs_stype_region = df_tfs_stype_region, tfs_sample_region = df_tfs_sample_region,
                                tfs_top100_mat = mat_top_acts_tfs_region)
    saveRDS(decoupleR_message, glue("{rds_dir}/decoupleR_message_{celltype}.rds"))
    saveRDS(srt_s, glue("{rds_dir}/srt_decoupleR_{celltype}.rds"))

    return(celltype)
}, BPPARAM = param)



