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

param <- BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE)

## find marker gene

# tumor erosion region
srt <- readRDS(glue("{rds_dir}/srt_erosion_region.rds"))
gene_entrez <- bitr(unique(rownames(srt)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

select_celltypes <- c("Macrophage", "Epithelial", "T cell")
names_celltypes <- c("Macrophage", "Epithelial", "T")
results <- bplapply(select_celltypes, function(celltype){
    markers <- FindMarkers(srt, ident.1 = "erosion", group.by = "region", subset.ident = celltype)
    dt <- markers %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
        mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
        as.data.table(keep.rownames = "gene") %>%
        merge(gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")
    return(list(raw = markers, filter = dt))
}, BPPARAM = param)
names(results) <- names_celltypes
results_flat <- setNames(flatten(results), unlist(lapply(names(results), function(x) paste0(x, c("_diff_gene", "_diff_gene_sig")))))

saveRDS(results_flat, glue("{rds_dir}/erosion_region_diff_gene.rds"))
list2excel(results_flat, glue("{rds_dir}/erosion_region_diff_gene.xlsx"), showRow = T)

# tumor erosion sample
srt <- readRDS(glue("{rds_dir}/srt_erosion_sample.rds"))
gene_entrez <- bitr(unique(rownames(srt)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

select_celltypes <- c("Macrophage", "Epithelial", "T cell")
names_celltypes <- c("Macrophage", "Epithelial", "T")
results <- bplapply(select_celltypes, function(celltype){
    markers <- FindMarkers(srt, ident.1 = "granuloma", group.by = "group", subset.ident = celltype)
    dt <- markers %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
        mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
        as.data.table(keep.rownames = "gene") %>%
        merge(gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")
    return(list(raw = markers, filter = dt))
}, BPPARAM = param)
names(results) <- names_celltypes
results_flat <- setNames(flatten(results), unlist(lapply(names(results), function(x) paste0(x, c("_diff_gene", "_diff_gene_sig")))))

saveRDS(results_flat, glue("{rds_dir}/erosion_sample_diff_gene.rds"))
list2excel(results_flat, glue("{rds_dir}/erosion_sample_diff_gene.xlsx"), showRow = T)

#epithelial subtype
srt <- readRDS(glue("{rds_dir}/xenium_sketch_celltyped.rds"))
DefaultAssay(srt) <- "Xenium"
srt_epi <- subset(srt, subset = cell_type == "Epithelial")
srt_epi$subtype <- case_when(str_detect(as.character(srt_epi$Epithelial.1.subtype),"spinous") ~ "spinous", 
                              TRUE ~ as.character(srt_epi$Epithelial.1.subtype))
DefaultAssay(srt_epi) <- "Xenium"
srt_s <- JoinLayers(srt_epi)
Idents(srt_s) <- "subtype"
gene_entrez <- bitr(unique(rownames(srt_s)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


case_groups <- unique(srt_s$subtype) %>% na.omit() %>% .[!. %in% c("spinous", "basal_cell", "mix_epi")]
results <- bplapply(case_groups, function(case){
  markers1 <- FindMarkers(srt_s, ident.1 = case, ident.2 = "spinous", group.by = "subtype", 
                         only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
  markers2 <- FindMarkers(srt_s, ident.1 = case, ident.2 = "basal_cell", group.by = "subtype", 
                         only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

  markers1 <- markers1 %>% mutate(case = case, control = "spinous", gene = rownames(.)) %>% 
      merge(gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL") %>% as.data.frame()
  markers2 <- markers2 %>% mutate(case = case, control = "basal_cell", gene = rownames(.)) %>% 
      merge(gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL") %>% as.data.frame()

  return(list(vs_spinous = markers1, vs_basal_cell = markers2))
}, BPPARAM = param)

names(results) <- case_groups

ls_vs_spinous <- lapply(results, function(x) x$vs_spinous)
ls_vs_basal_cell <- lapply(results, function(x) x$vs_basal_cell)

saveRDS(ls_vs_spinous, glue("{rds_dir}/epithelial_cluster_marker_gene_vs_spinous.rds"))
saveRDS(ls_vs_basal_cell, glue("{rds_dir}/epithelial_cluster_marker_gene_vs_basal_cell.rds"))
list2excel(ls_vs_spinous, glue("{rds_dir}/epithelial_cluster_marker_gene_vs_spinous.xlsx"), showRow = T)
list2excel(ls_vs_basal_cell, glue("{rds_dir}/epithelial_cluster_marker_gene_vs_basal_cell.xlsx"), showRow = T)



### run kegg

erosion_region_diff_gene <- readRDS(glue("{rds_dir}/erosion_region_diff_gene.rds"))
erosion_sample_diff_gene <- readRDS(glue("{rds_dir}/erosion_sample_diff_gene.rds"))
names(erosion_region_diff_gene) <- paste0("region_", names(erosion_region_diff_gene))
names(erosion_sample_diff_gene) <- paste0("sample_", names(erosion_sample_diff_gene))
erosion_diff_gene <- c(erosion_region_diff_gene, erosion_sample_diff_gene)

plot_group <- grep("_sig$", names(erosion_diff_gene), value = TRUE)
results <- BiocParallel::bplapply(plot_group, function(group){
    dt <- erosion_diff_gene[[group]]
    geneLst_up <- dt[dt$type == "up",]$ENTREZID
    geneLst_down <- dt[dt$type == "down",]$ENTREZID
    res_up <- enrichKEGG(geneLst_up, organism = "hsa", keyType = "kegg", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
                       minGSSize = 10,
                       use_internal_data = TRUE)
    res_dn <- enrichKEGG(geneLst_down, organism = "hsa", keyType = "kegg", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
                       minGSSize = 10,
                       use_internal_data = TRUE)
    res_up_readable <- setReadable(res_up, OrgDb = "org.Hs.eg.db", keyType="ENTREZID") %>% as.data.frame()
    res_dn_readable <- setReadable(res_dn, OrgDb = "org.Hs.eg.db", keyType="ENTREZID") %>% as.data.frame()

    return(list(res_up = res_up, res_dn = res_dn, readable_up = res_up_readable, readable_dn = res_dn_readable))
}, BPPARAM = param)
names(results) <- gsub("_diff_gene_sig$", "", plot_group)
results_flat <- setNames(flatten(results), 
                        unlist(lapply(names(results), function(x) paste(x, names(results[[x]]), sep = "_"))))
results_readable <- results_flat[grepl("readable", names(results_flat))]
saveRDS(results_flat, glue("{rds_dir}/kegg_erosion_result.rds"))

results_flat <- readRDS(glue("{rds_dir}/kegg_erosion_result.rds"))
results_readable <- results_flat[grepl("readable", names(results_flat))]
list2excel(results_readable, glue("{rds_dir}/kegg_erosion_results.xlsx"), showRow = T)



results_macro_up <- results_flat[grepl("region_Macrophage", names(results_flat)) & grepl("_up", names(results_flat))]
saveRDS(results_macro_up, glue("{rds_dir}/kegg_erosion_result_region_macrophage_up.rds"))


#### epithelial_subtype_markers
ls_vs_spinous <- readRDS(glue("{rds_dir}/epithelial_cluster_marker_gene_vs_spinous.rds"))
ls_vs_basal_cell <- readRDS(glue("{rds_dir}/epithelial_cluster_marker_gene_vs_basal_cell.rds"))
df_vs_spinous <- do.call(rbind, ls_vs_spinous)
df_vs_basal_cell <- do.call(rbind, ls_vs_basal_cell)

groups <- names(ls_vs_spinous)
results <- BiocParallel::bplapply(groups, function(group){
  dt1 <- ls_vs_spinous[[group]]
  dt2 <- ls_vs_basal_cell[[group]]
  
  res_1 <- enrichKEGG(dt1$ENTREZID, organism = "hsa", keyType = "kegg", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
                       minGSSize = 10,
                       use_internal_data = TRUE)
  res_2 <- enrichKEGG(dt2$ENTREZID, organism = "hsa", keyType = "kegg", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
                      minGSSize = 10,
                      use_internal_data = TRUE)
  readable1 <- setReadable(res_1, OrgDb = "org.Hs.eg.db", keyType="ENTREZID") %>% as.data.frame()
  readable2 <- setReadable(res_2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID") %>% as.data.frame()
  return(list(res_vs_spinous = res_1, res_vs_basal_cell = res_2,
              readable_vs_spinous = readable1, readable_vs_basal_cell = readable2))
}, BPPARAM = param)

names(results) <- groups
saveRDS(results, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene.rds"))

results <- readRDS(glue("{rds_dir}/kegg_epithelial_cluster_marker_gene.rds"))
ls_vs_spinous_res <- lapply(results, function(x) x$res_vs_spinous)
ls_vs_basal_cell_res <- lapply(results, function(x) x$res_vs_basal_cell)
ls_vs_spinous_readable <- lapply(results, function(x) x$readable_vs_spinous)
ls_vs_basal_cell_readable <- lapply(results, function(x) x$readable_vs_basal_cell)


saveRDS(ls_vs_spinous_res, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_spinous.rds"))
saveRDS(ls_vs_basal_cell_res, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_basal_cell.rds"))
list2excel(ls_vs_spinous_readable, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_spinous.xlsx"), showRow = T)
list2excel(ls_vs_basal_cell_readable, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_basal_cell.xlsx"), showRow = T)


ls_vs_basal_cell_res <- readRDS(glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_basal_cell.rds"))
cases <- c("basal_invasive", "basal_KRAS+", "CCND2+SFRP1", "basal_CD36+", "basal_EPCAM+")
control <- "basal_cell"

select_pathway <- list(
    "basal_invasive" = c("hsa04510", "hsa04512", "hsa04151", "hsa04066", "hsa04068", "hsa04370", "hsa04010", "hsa04350"),
    "basal_KRAS+" = c("hsa04015", "hsa04370", "hsa04012", "hsa04014", "hsa04066", "hsa04390", "hsa04151", "hsa04010", "hsa04068"),
    "CCND2+SFRP1" = c("hsa05203", "hsa04151", "hsa05165", "hsa04012", "hsa04330", "hsa04010"),
    "basal_CD36+" = c("hsa05203", "hsa04310", "hsa05165"),
    "basal_EPCAM+" = c("hsa04010", "hsa04340", "hsa04151", "hsa04012", "hsa04390", "hsa04310")
)

results <- bplapply(cases, function(case){
    res <- ls_vs_basal_cell_res[[case]]
    path_case <- select_pathway[[case]]
    result_filter <- res@result %>% filter(ID %in% path_case) %>% mutate(case = case, control = control)
    res_filter <- new("enrichResult", result = result_filter)
    return(list(df = result_filter, res = res_filter))
}, BPPARAM = param)
names(results) <- cases
ls_results_df <- lapply(results, function(x) x$df)
ls_results_res <- lapply(results, function(x) x$res)

results <- list(df = ls_results_df, res = ls_results_res)


saveRDS(results, glue("{rds_dir}/kegg_epithelial_cluster_marker_gene_vs_basal_cell_select.rds"))


