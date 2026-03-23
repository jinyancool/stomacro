.libPaths(c(.libPaths(), "/cluster/home/yliang_jh/sbin/R/library/4.3.0"))
pkgs <- c("fs", "configr", "stringr", 
          "jhtools", "glue", "patchwork", "tidyverse", "dplyr", "Seurat", "magrittr", 
          "readxl", "writexl", "ComplexHeatmap", "openxlsx",
          "data.table", "ggplot2", "ggbeeswarm",
          "sf", "corrplot", "ggpubr", "survival", "survminer", "forestmodel", "BiocParallel", "BiocNeighbors")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



rds_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/codex"



save_list_to_excel <- function(data_list, filename, font_size = 12, font_name = "Arial"){
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = font_size, fontName = font_name)
  for(sheet_name in names(data_list)){
    data <- data_list[[sheet_name]]
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data)
    setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")

    header_style <- createStyle(fontSize = font_size, fontName = font_name, textDecoration = "bold",halign = "center")
    borders_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin", borderColour = "black")

    addStyle(wb, sheet_name, borders_style, rows = 1:(nrow(data) + 1),cols = 1:ncol(data), gridExpand = TRUE)
    addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(data))
    freezePane(wb, sheet_name, firstActiveRow = 2, firstActiveCol = 1)
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}



dtdir <- file.path("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/data/cell_segmentation")
files <- list.files(dtdir, pattern = ".tsv")
data <- lapply(files, function(f){
  dt <- fread(file.path(dtdir, f), colClasses = c(Image = "character"))
  dt <- dt %>% 
    select(Image, `Centroid X µm`, `Centroid Y µm`, contains("Cell: Mean"), -contains("DAPI")) %>%
    rename_with(~ str_remove_all(.x, ": Cell: Mean")) %>%
    mutate(Image = sub(" [(]new[)]", "", Image))
  return(dt)
}) %>% rbindlist()
data <- data %>% group_by(Image) %>% mutate(cell_id = paste0(Image, "_", row_number())) %>% ungroup() %>% relocate(cell_id)
data <- data %>% mutate(sum = rowSums(across(where(is.numeric) & !c(`Centroid X µm`, `Centroid Y µm`)), na.rm = TRUE)) %>%
  relocate(`Centroid X µm`, `Centroid Y µm`, .after = everything())
saveRDS(data, glue("{rds_dir}/codex_cell_segmentation.rds"))








img_cluster_anno <- list()
Images <- list.files("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/doc/anno_split3/") %>% sub(".tsv", "", .)
for(img in Images){
	df_anno <- fread(paste0("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/doc/anno_split3/", img, ".tsv"))
	df_anno$Image <- img
  df_anno <- df_anno[, c("cluster", "celltype", "Image")]
  names(df_anno) <- c("cluster", "celltype", "Image")
	img_cluster_anno[[img]] <- df_anno
}
save_list_to_excel(img_cluster_anno, glue("{rds_dir}/anno_split3_img.xlsx"), font_size = 12, font_name = "Arial")


# cell_anno <- list()
# anno_rds_dir <- "/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/output/anno_split/"
# Images <- list.files("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/doc/anno_split3/") %>% sub(".tsv", "", .)
# for(img in Images){
#   srt <- readRDS(paste0(anno_rds_dir, "srt_", img, "_anno.rds"))
# 	dt <- srt@meta.data[, c("cell_id", "Image", "pg_cluster", "celltype")] %>% as.data.frame()
# 	names(dt) <- c("cell_id", "Image", "pg_cluster", "celltype")
# 	cell_anno[[img]] = dt
# }

# combined_anno <- do.call(rbind, cell_anno)
# rownames(combined_anno) <- combined_anno$cell_id

# saveRDS(combined_anno, glue("{rds_dir}/celltype_anno_split3.rds"))

srat <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
dt_celltype <- srat@meta.data %>% select(cell_id, Image, pg_cluster, celltype)
saveRDS(dt_celltype, glue("{rds_dir}/codex_celltype_final.rds"))



subtype_match <- list()
subtype_anno <- list()
suffix <- c("cd4t", "cd8t", "fib", "tumor")
list_match <- list(cd4t = "CD4 T", cd8t = "CD8 T", mpg = "Macrophage", fib = "Fibroblast", tumor = "Tumor cell")
list_match_tsv <- list(cd4t = "cd4t", cd8t = "cd8t", mpg = "macrophage", fib = "fibroblast", tumor = "tumor")
for(s in suffix){
	srt_s <- readRDS(paste0("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/srt_", s, ".rds"))
  dt <- srt_s@meta.data %>% as.data.frame()
  dt <- dt[, c("cell_id", "Image", "pg_cluster", "celltype", "subtype")]
  celltype <- list_match[[s]]
	subtype_anno[[celltype]] <- dt

  s_tsv = list_match_tsv[[s]]
  df_anno <- read_tsv(paste0("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/doc/anno_subtype/", s_tsv, ".tsv"))
  df_anno$celltype <- celltype
  subtype_match[[celltype]] <- df_anno
}

subtype_anno[["Macrophage"]] <- readRDS(glue("{rds_dir}/subtype_macrophage.rds"))

combined_subtype_anno <- do.call(rbind, subtype_anno)
saveRDS(combined_subtype_anno, glue("{rds_dir}/subtype_anno.rds"))

combined_subtype_match <- do.call(rbind, subtype_match)
saveRDS(combined_subtype_match, glue("{rds_dir}/subtype_match_cluster.rds"))





#f_celltype_cn_cluster <- "/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/output/celltype_cn/celltype_cn_cluster_10.csv"
#celltype_cn_cluster <- read_csv(f_celltype_cn_cluster)
#f_subtype_cn_cluster <- "/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/output/subtype_cn/subtype_cn_cluster_35.csv"
#subtype_cn_cluster <- read_csv(f_subtype_cn_cluster)
#f_celltype_composition_cluster <- "/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/output/celltype_freq/celltype_freq_cluster.csv"
#celltype_composition_cluster <- read_csv(f_celltype_composition_cluster)

#cluster_message <- list(celltype_cn = celltype_cn_cluster, 
#                        subtype_cn = subtype_cn_cluster, 
#                        celltype_composition = celltype_composition_cluster)
#saveRDS(cluster_message, glue("{rds_dir}/codex_cluster_message.rds"))



cluster_message <- readRDS(glue("{rds_dir}/codex_cluster_message_bak.rds"))
subtype_cn <- cluster_message$subtype_cn


rename_list <- list(
  "Macrophage_CD44+" = "Macrophage_CD44+", 
  "Macrophage_CD45RO+_CD44+" = "Macrophage_CD44+_Vimentin+", 
  "Macrophage_CD45RO+PD-1" = "Macrophage_CD14+_Ki67+",
  "Macrophage_FOXP3+" = "Macrophage_IDO1+",
  "Macrophage_GranzymeB+_PD-1+" = "Macrophage_PD-L1+",
  "Macrophage_HLA-DR+" = "Macrophage_HLA-DR+",
  "Macrophage_Ki67+" = "Macrophage_Ki67+",
  "Macrophage_others" = "Macrophage_others"
)

subtype_cn <- subtype_cn %>%rename_with(~ recode(., !!!rename_list, .default = .))
cluster_message$subtype_cn <- subtype_cn
saveRDS(cluster_message, glue("{rds_dir}/codex_cluster_message.rds"))



dt <- fread("/cluster/home/yliang_jh/projects/IMC/codex/oral_lvjiong/doc/erosion_zone_immune_cell_count.csv")
list_sample_match <- list("S1" = "3-12-4", "S2" = "3-4-3", "S3" = "3-6-1")
dt$sample_id <- recode(dt$sample_id, !!!list_sample_match)
saveRDS(dt, glue("{rds_dir}/erosion_zone_immune_cell_count.rds"))




