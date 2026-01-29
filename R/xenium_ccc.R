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
run_commot <- function(sample, spe, db, rs, rds_dir){
  sub = spe[, spe$sample_id == sample]
  sub <- sub[unique(c(db$ligand, rs)), ]
  xy <- spatialCoords(sub)
  reducedDim(sub, "spatial") <- xy
  ad <- SCE2AnnData(sub, X_name="logcounts")
  ct <- import("commot")
  ct$tl$spatial_communication(ad,
    database_name="CellChatDB",
    dis_thr=10,  
    df_ligrec=db,
    heteromeric=TRUE,
    pathway_sum=TRUE,
    heteromeric_rule="min",
    heteromeric_delimiter="_")
  ccc <- list(
    s=ad$obsm["commot-CellChatDB-sum-sender"],
    r=ad$obsm["commot-CellChatDB-sum-receiver"])
  
  for (df in ccc) colData(sub)[names(df)] <- df
  as <- lapply(ccc, \(.) {
    names(.) <- gsub("^(s|r)-", "", names(.))
    as(t(as.matrix(.)), "dgCMatrix")
  })
  sce <- SingleCellExperiment(as, colData=colData(sub))
  saveRDS(sce, glue("{rds_dir}/CCC2/sce_{sample}.rds"))
  return(ccc)
}

run_pathway <- function(sample, rds_dir){
    sce <- readRDS(glue("{rds_dir}/sce_{sample}.rds"))
    sce$Epithelial_ccc <- coalesce(sce$group, sce$cell_type)
    mu <- aggregateAcrossCells(sce, ids=sce$Epithelial_ccc, statistics="mean", use.assay.type=assayNames(sce))
  
    ks <- expand.grid(ks <- colnames(mu), ks)
    cs <- split(colnames(sce), sce$Epithelial_ccc)
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
            source <- sce[lr, cs[[i]]]
            target <- sce[lr, cs[[j]]]
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




#spinfo <- sampleinfo$xenium_epithelial_for_cc
#spinfo <- spinfo %>% filter(communication %in% c(1,2)) %>% dplyr::select(Epithelial.1.subtype, group)
#
#
#ct <- import("commot")
#db <- ct$pp$ligand_receptor_database(
#    species="human", 
#    database="CellChat", 
#    signaling_type=NULL)
#names(db) <- c("ligand", "receptor", "pathway", "type")
#
#lr <- c(db$ligand, db$receptor)
#ss <- strsplit(lr, "-")
#gs <- sapply(ss, .subset, 1)
#
#keep <- apply(db, 1, \(.) {
#    rs <- strsplit(.["receptor"], "_")
#    lr <- c(.["ligand"], unlist(rs))
#    all(lr %in% rownames(spe))
#})
#db <- db[keep, ]
#
#ss <- strsplit(db$receptor, "_")
#rs <- sapply(ss, .subset, 1)
#
#
#spe <- readRDS(glue("{rds_dir}/spe_0066253.rds"))
#dt <- left_join(as.data.frame(colData(spe)), spinfo, by = "Epithelial.1.subtype")
#spe$group <- dt$group
#counts_matrix <- assay(spe, "counts")
#keep_genes <- rowSums(counts_matrix>0) > 0
#spe$library_size <- librarySizeFactors(spe)
#min_factor <- min(spe$library_size[spe$library_size > 0]) #0.002
#spe$library_size[spe$library_size <= 0] <- 0.0001
#spe1 <- logNormCounts(spe, size.factors=spe$library_size)
#
#keep1 <- apply(db, 1, \(.) {
#    rs <- strsplit(.["receptor"], "_")
#    lr <- c(.["ligand"], unlist(rs))
#    all(lr %in% rownames(spe1))
#})
#db1 <- db[keep1, ]
#ss1 <- strsplit(db1$receptor, "_")
#rs1 <- sapply(ss1, .subset, 1)
#
#
#spe <- readRDS(glue("{rds_dir}/spe_0066266.rds"))
#spe <- spe[, !spe$sample_id %in% c("F4-2", "F14")]
#dt <- left_join(as.data.frame(colData(spe)), spinfo, by = "Epithelial.1.subtype")
#spe$group <- dt$group
#counts_matrix <- assay(spe, "counts")
#keep_genes <- rowSums(counts_matrix>0) > 0
#spe$library_size <- librarySizeFactors(spe)
#min_factor <- min(spe$library_size[spe$library_size > 0]) #0.002
#spe$library_size[spe$library_size <= 0] <- 0.0001
#spe2 <- logNormCounts(spe, size.factors=spe$library_size)
#
#keep2 <- apply(db, 1, \(.) {
#    rs <- strsplit(.["receptor"], "_")
#    lr <- c(.["ligand"], unlist(rs))
#    all(lr %in% rownames(spe2))
#})
#db2 <- db[keep2, ]
#ss2 <- strsplit(db2$receptor, "_")
#rs2 <- sapply(ss2, .subset, 1)
#
#samples_id <- c(unique(spe1$sample_id), unique(spe2$sample_id))


#ccc_results_1 <- bplapply(unique(spe1$sample_id), run_commot, spe = spe1, db = db1, rs = rs1, rds_dir = rds_dir, BPPARAM = param)
#ccc_results_2 <- bplapply(unique(spe2$sample_id), run_commot, spe = spe2, db = db2, rs = rs2, rds_dir = rds_dir, BPPARAM = param)
#ccc_results <- c(ccc_results_1, ccc_results_2)
#names(ccc_results) <- c(nique(spe1$sample_id), unique(spe2$sample_id))
#saveRDS(ccc_results, glue("{rds_dir}/CCC2/ccc_commot_result.rds"))


spe1 <- readRDS(glue("{rds_dir}/spe_0066253.rds"))
spe2 <- readRDS(glue("{rds_dir}/spe_0066266.rds"))
spe2 <- spe2[, !spe2$sample_id %in% c("F4-2", "F14")]
samples_id <- c(unique(spe1$sample_id), unique(spe2$sample_id))

mu_list <- list()
df_list <- list()
rds_h_dir <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/xenium/ccc"
res_list <- bplapply(samples_id, run_pathway, rds_dir = rds_h_dir, BPPARAM = param)
names(res_list) <- samples_id

mu_list <- lapply(res_list, `[[`, "mu")
df_list <- lapply(res_list, `[[`, "df")

saveRDS(mu_list, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_mu_list.rds"))
saveRDS(df_list, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path_list.rds"))

df <- do.call(rbind, df_list)
rownames(df) <- NULL
df <- df %>% mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-"))
df <- df %>% filter(!sample_id %in% c("F14", "F4-2")) %>% as.data.frame()
saveRDS(df, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path.rds"))


df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path.rds"))
cc_path_s <- data.frame(path = character(),target_cell = character(),mean_stem = character(),mean_diff = character(),
    t_statistic = numeric(),df = numeric(),p_value = numeric(),stringsAsFactors = FALSE)
for (p in unique(df$path)){
  path_data <- df[df$path == p, ]
  targets <- unique(path_data$target) %>% .[! . %in% c("stemness", "differentiation", "Epithelial")]
  for (target in targets){
    source_stem_data <- path_data[(path_data$source %in% c("stemness"))&(path_data$target %in% c(target)),][["score"]] %>% na.omit()
    source_diff_data <- path_data[(path_data$source %in% c("differentiation"))&(path_data$target %in% c(target)),][["score"]] %>% na.omit()
    test_result <- tryCatch({
      t.test(source_stem_data, source_diff_data, var.equal = FALSE)
      }, error = function(e) {
      return(NULL)
      })
    if (!is.null(test_result)){
      mean1 <- mean(source_stem_data, na.rm = TRUE)
      mean2 <- mean(source_diff_data, na.rm = TRUE)
      cc_path_add <- data.frame(
        path = p, target_cell = target, mean_stem = mean1, mean_diff = mean2,
        t_statistic = test_result$statistic, df = test_result$parameter, p_value = test_result$p.value, stringsAsFactors = FALSE)
      cc_path_s <- rbind(cc_path_s, cc_path_add)
    }
  }
}


cc_path_r <- data.frame(path = character(),source_cell = character(),mean_stem = character(),mean_diff = character(),
    t_statistic = numeric(),df = numeric(),p_value = numeric(),stringsAsFactors = FALSE)
for (p in unique(df$path)){
  path_data <- df[df$path == p, ]
  sources <- unique(path_data$source) %>% .[! . %in% c("stemness", "differentiation", "Epithelial")]
  for (source in sources){
    target_stem_data <- path_data[(path_data$source %in% c(source))&(path_data$target %in% c("stemness")),][["score"]] %>% na.omit()
    target_diff_data <- path_data[(path_data$source %in% c(source))&(path_data$target %in% c("differentiation")),][["score"]] %>% na.omit()
    test_result <- tryCatch({
      t.test(target_stem_data, target_diff_data, var.equal = FALSE)
      }, error = function(e) {
      return(NULL)
      })
    if (!is.null(test_result)){
      mean1 <- mean(target_stem_data, na.rm = TRUE)
      mean2 <- mean(target_diff_data, na.rm = TRUE)
      cc_path_add <- data.frame(
        path = p, source_cell = source, mean_stem = mean1, mean_diff = mean2,
        t_statistic = test_result$statistic, df = test_result$parameter, p_value = test_result$p.value, stringsAsFactors = FALSE)
      cc_path_r <- rbind(cc_path_r, cc_path_add)
    }
  }
}


cc_path_s_filter <- cc_path_s %>% filter(p_value <= 0.05) %>% 
    arrange(desc(abs(t_statistic)), p_value) %>% as.data.frame()
cc_path_r_filter <- cc_path_r %>% filter(p_value <= 0.05) %>% 
    arrange(desc(abs(t_statistic)), p_value) %>% as.data.frame()
results_cc_path <- list("source_filter" = cc_path_s_filter, "source_raw" = cc_path_s, "target_filter" = cc_path_r_filter, "target_raw" = cc_path_r)
saveRDS(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_CC_path.rds"))
list2excel(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_CC_path.xlsx"))



cc_path <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_CC_path.rds"))
cc_path_s_filter <- cc_path$source_filter
cc_path_r_filter <- cc_path$target_filter

cc_path_s_filter <- cc_path_s_filter %>% filter(!target_cell %in% c("Epithelial", "stemness", "differentiation", "Glial", "skeletalMC", "smoothMC"))
select_path_s <- cc_path_s_filter %>% head(200) %>% pull(path) %>% unique()
path1 <- c("LAMB3-CD44", "LAMA3-CD44", "LAMC2-CD44", "JAG1-NOTCH4", "JAG1-NOTCH3", "JAG1-NOTCH2", "JAG1-NOTCH1",
          "DLL1-NOTCH4", "DLL1-NOTCH3", "DLL1-NOTCH2", "DLL3-NOTCH4", "DLL4-NOTCH3",
          "FGF1-FGFR2", "FGF1-FGFR3", "FGF9-FGFR2", "FGF1-FGFR4", "FGF2-FGFR3", "TGFA-EGFR")
path2 <- cc_path_s_filter %>% filter(target_cell %in% c("Macrophage","T cell")) %>% pull(path) %>% unique() %>% head(15)
path3 <- cc_path_s_filter %>% filter(!target_cell %in% c("Fibroblast","Macrophage","T cell")) %>% pull(path) %>% unique() %>% head(15)
rm_path <- c("DLL3-NOTCH4","DLL4-NOTCH3","APLN-APLNR","DLL1-NOTCH4","JAG1-NOTCH4","DLL4-NOTCH3","JAG1-NOTCH3","DLL1-NOTCH3")
select_path_s <- c(unique(c(path1,path2,path3))) %>% .[!. %in% rm_path]


cc_path_r_filter <- cc_path_r_filter %>% filter(!source_cell %in% c("Epithelial", "stemness", "differentiation", "Glial", "skeletalMC", "smoothMC"))
path1 <- c("LAMB3-CD44", "LAMA3-CD44", "LAMC1-CD44", "LAMC2-CD44", "JAG1-NOTCH4", "JAG1-NOTCH3", "JAG1-NOTCH2",
          "DLL1-NOTCH4", "DLL1-NOTCH3", "DLL1-NOTCH2", "DLL1-NOTCH1", "DLL4-NOTCH2", "DLL4-NOTCH3", "DLL3-NOTCH1", "DLL3-NOTCH4",
          "FGF1-FGFR2", "FGF1-FGFR3", "FGF1-FGFR4", "FGF9-FGFR2", "FGF2-FGFR2", "FGF2-FGFR3", "FGF9-FGFR2", "FGF7-FGFR2", 
          "EGF-EGFR", "COL4A2-CD44")
path2 <- cc_path_r_filter %>% filter(source_cell %in% c("Macrophage","T cell")) %>% pull(path) %>% unique() %>% head(15)
path3 <- cc_path_r_filter %>% filter(!source_cell %in% c("Fibroblast","Macrophage","T cell")) %>% pull(path) %>% unique() %>% head(15)
rm_path <- c("DLL4-NOTCH2", "DLL4-NOTCH3", "DLL1-NOTCH4", "DLL3-NOTCH1", "DLL3-NOTCH4", "JAG1-NOTCH4")
select_path_r <- c(unique(c(path1,path2,path3))) %>% .[!. %in% rm_path]


df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path.rds"))

df_heat_cp_s <- df %>% 
    filter(source %in% c("stemness", "differentiation"), 
          !target %in% c("Epithelial", "stemness", "differentiation", "Glial", "skeletalMC", "smoothMC"), score > 0) %>%
    group_by(source, target, ligand, receptor) %>%       
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_s)
cell_cell_order <- c(
  paste("stemness", intersect(config$cell_type_order, unique(df$target)), sep="-"),
  paste("differentiation", intersect(config$cell_type_order, unique(df$target)), sep="-")
) %>% .[. %in% unique(df_heat_cp_s$cell_cell)]
df_heat_cp_wide_s <- df_heat_cp_s %>% dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, cell_cell_order, drop = FALSE]


df_heat_cp_r <- df %>% 
    filter(target %in% c("stemness", "differentiation"), 
            !source %in% c("Epithelial", "stemness", "differentiation", "Glial", "skeletalMC", "smoothMC"), score > 0) %>%
    group_by(source, target, ligand, receptor) %>%       
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') %>%
    mutate(path = paste(ligand, receptor, sep = "-"), cell_cell = paste(source, target, sep = "-")) %>%
    filter(path %in% select_path_r)
cell_cell_order <- c(
  paste(intersect(config$cell_type_order, unique(df$source)), "stemness", sep="-"),
  paste(intersect(config$cell_type_order, unique(df$source)), "differentiation", sep="-")
) %>% .[. %in% unique(df_heat_cp_r$cell_cell)]
df_heat_cp_wide_r <- df_heat_cp_r %>% dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, cell_cell_order, drop = FALSE]


df_heat_cp_wide <- list(s = df_heat_cp_wide_s, r = df_heat_cp_wide_r)


saveRDS(df_heat_cp_wide, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_plot.rds"))






