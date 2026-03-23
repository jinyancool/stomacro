suppressPackageStartupMessages({
  library(magrittr)
  library(Seurat)
  library(sf)
  library(ggplot2)
  library(arrow)
  library(data.table)
  library(readr)
  library(glue)
  library(patchwork)
  library(tidyverse)
  library(dplyr)
  library(BPCells)
  library(SeuratWrappers)
  library(purrr)
  library(readxl)
  library(openxlsx)
  library(ComplexHeatmap)
})

rds_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/codex"


find_group_mks_and_topAnno <- function(srat.obj, assay_name, res_prefix, resoLst, out_dir){

  srat_tmp0 <- srat.obj
  DefaultAssay(srat_tmp0) <- assay_name
  srat_tmp0[[assay_name]] <- JoinLayers(srat_tmp0[[assay_name]])

  mks_lst <- list()
  for(clust_name in paste0(res_prefix, resoLst)){
    srat_tmp <- srat_tmp0
    Idents(srat_tmp) <- srat_tmp@meta.data[[clust_name]]
    mks_lst[[clust_name]] <- Seurat::FindAllMarkers(srat_tmp, assay = assay_name, only.pos = TRUE, min.pct = 0.01,logfc.threshold = 0.2) %>%
      tibble() %>% mutate(cluster_name = clust_name)
    rm(srat_tmp)
    gc()
  }
  writexl::write_xlsx(mks_lst, glue("{out_dir}/cluster_mks_total.xlsx"))

  mks_anno <- list()
  for (clust_name in paste0(res_prefix, resoLst)) {
    mks <- mks_lst[[clust_name]] %>% as.data.frame() %>%
      filter(p_val_adj<=0.01) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) %>% 
      ungroup() %>% dplyr::select(c("cluster", "gene", "avg_log2FC"))
  
    mks_anno[[clust_name]] <- left_join(mks, panel, by = "gene")
  }
  writexl::write_xlsx(mks_anno, glue("{out_dir}/cluster_mks_top50_anno.xlsx"))

  rm(srat_tmp0)
}



Macrophage_7 = c("CD14", "CD44", "HLA-DR", "IDO1", "Ki67", "Vimentin", "PD-L1")



srat <- readRDS(glue("{rds_dir}/srt0_split_anno.rds"))
sub_srat <- subset(srat, subset = celltype == "Macrophage")

sub_srat <- NormalizeData(sub_srat, normalization.method = "LogNormalize", scale.factor = 1000)
sub_srat <- FindVariableFeatures(sub_srat, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
sub_srat <- ScaleData(sub_srat, verbose = FALSE)

meta_old <- sub_srat@meta.data %>% mutate(cell_id = rownames(.)) %>% select(cell_id, subtype)
colnames(meta_old) <- c("cell_id", "subtype_old")
macro_sub <- unique(srat$subtype) %>% na.omit() %>% .[startsWith(., "Macrophage_")]

p <- DoHeatmap(sub_srat, 
          features = Macrophage_7,
          group.by = "subtype",
          slot = "scale.data")
ggsave("/cluster/home/lixiyue_jh/projects/quarto/stomatology_lvjiong/quarto_book/results_old/test_heatmap.pdf", p)


avg_expr <- AggregateExpression(
  sub_srat,
  assays = "CODEX",
  features = Macrophage_7,
  group.by = "subtype",
  slot = "scale.data",
  return.seurat = FALSE
)
avg_matrix <- as.matrix(avg_expr[[1]])
avg_matrix <- scale(avg_matrix)


p <- Heatmap(avg_matrix,
        name = "Expression",
        column_title = "Macrophage Subtypes",
        row_title = "Genes",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        cluster_rows = TRUE,
        cluster_columns = TRUE)

pdf("/cluster/home/lixiyue_jh/projects/quarto/stomatology_lvjiong/quarto_book/results_old/test_heatmap2.pdf",width = 7, height = 7)
draw(p)
dev.off()




srat_new <- readRDS(glue("{rds_dir}/srt_new_split_anno.rds"))
sub_srat_new <- subset(srat_new, subset = celltype == "Macrophage")
macro_sub_new <- unique(srat_new$subtype) %>% na.omit() %>% .[startsWith(., "Macrophage_")]
meta_new <- sub_srat_new@meta.data %>% mutate(cell_id = rownames(.)) %>% select(cell_id, subtype)
colnames(meta_new) <- c("cell_id", "subtype_new")

meta <- left_join(meta_old, meta_new, by="cell_id") %>% select(subtype_old, subtype_new)
meta <- meta[, c("subtype_old", "subtype_new")]

table_2 <- table(meta)
table_2 <- prop.table(table_2, margin = 1)
write.xlsx(table_2,"/cluster/home/lixiyue_jh/projects/quarto/stomatology_lvjiong/quarto_book/results_old/test_macro.xlsx", colNames=TRUE, rowNames=TRUE)


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
meta_rename <- sub_srat@meta.data %>% mutate(subtype = recode(subtype, !!!rename_list)) %>% 
  select(cell_id, Image, pg_cluster, celltype, subtype)
saveRDS(meta_rename, glue("{rds_dir}/subtype_macrophage.rds"))


table(srat$celltype, useNA="ifany")
table(srat_new$celltype, useNA="ifany")

table(srat$subtype, useNA="ifany")
table(srat_new$subtype, useNA="ifany")



srat <- readRDS(glue("{rds_dir}/srt0_split_anno.rds"))
meta_rename <- srat@meta.data %>% mutate(subtype = recode(subtype, !!!rename_list, .default = subtype)) %>% as.data.frame()
meta_rename <- meta_rename[colnames(srat), ]
srat$subtype <- meta_rename$subtype
saveRDS(srat, glue("{rds_dir}/srt_split_anno.rds"))