suppressPackageStartupMessages({
  library(magrittr)
  library(Seurat)
  library(SeuratDisk)
  library(SpatialExperiment)
  library(SpatialExperimentIO)
  library(BiocParallel)
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
})

rds_dir = "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/xenium"


## 
sampleinfo <- readRDS("/cluster/home/jhuang/projects/stomatology/docs/lvjiong/sampleinfo/sampleinfo.rds")
dt_sampleinfo <- sampleinfo$xenium



## preprocess message
dir_sp_old <- "/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/Xenium_sample_section/0066253/"
files <- list.files(dir_sp_old)
dt_sp_old <- lapply(files, function(f){
  f_path <- paste0(dir_sp_old, f)
  dt <- fread(f_path)
  dt$sample_id <- strsplit(f, split = "_") %>% sapply('[[', 1)
  dt$type <- sub(".csv", "", f) %>% strsplit(split = "_") %>% sapply('[[', 2)
  return(dt)
}) %>% rbindlist() %>% arrange(type, .by_group = TRUE) %>% 
  distinct(`Cell ID`, sample_id, .keep_all = TRUE) %>% select(`Cell ID`, sample_id, type) %>% group_by(sample_id) %>% arrange(`Cell ID`, .by_group = TRUE) %>% ungroup()

dir_sp_new <- "/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/New_xenium_sample_section/0066253/"
files <- list.files(dir_sp_new, pattern="*_cells_stats.csv")
dt_sp_new <- lapply(files, function(f){
  f_path <- paste0(dir_sp_new, f)
  dt <- fread(f_path)
  dt$sample_id <- strsplit(f, split = "_") %>% sapply('[[', 1)
  return(dt)
}) %>% rbindlist() %>% group_by(sample_id) %>% arrange(`Cell ID`, .by_group = TRUE) %>% ungroup()

dt_diff <- dt_sampleinfo %>% select(sample_id, diff_level, differentiation, granuloma)
dt_sp <- left_join(dt_sp_new, dt_sp_old, by = c("Cell ID", "sample_id")) %>% left_join(dt_diff, by="sample_id")

srat <- LoadXenium("/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/0066253/outs/", segmentations = "cell", flip.xy = TRUE)
srat <- subset(srat, cells = dt_sp$`Cell ID`)
dt_sp <- dt_sp[match(colnames(srat), dt_sp$`Cell ID`), ] %>% select(-"Cell ID") %>% as.data.frame()
rownames(dt_sp) <- colnames(srat)
dt_sp <- dt_sp %>% as.data.frame()
srat <- AddMetaData(srat, metadata = dt_sp)
srat$slide <- "slide1"
srat$cell_ID <- paste(srat$slide, colnames(srat), sep='-')

saveRDS(srat, glue("{rds_dir}/sv5_xenium_object_0066253.rds"))
dim(srat) # 5001 1397413

#for (subsample in slide1_sample_id){
#  subsrat <- subset(srat, subset=sample_id==subsample)
#  counts_matrix <- GetAssayData(object = subsrat, assay = "Xenium", slot = "counts")
#  write_matrix_dir(mat = counts_matrix, dir = glue("{rds_dir}/matrix/{subsample}_matrix"), overwrite = TRUE)
#}



dir_sp_old <- "/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/Xenium_sample_section/0066266/"
files <- list.files(dir_sp_old)
dt_sp_old <- lapply(files, function(f){
  f_path <- paste0(dir_sp_old, f)
  dt <- fread(f_path)
  dt$sample_id <- strsplit(f, split = "_") %>% sapply('[[', 1)
  dt$type <- sub(".csv", "", f) %>% strsplit(split = "_") %>% sapply('[[', 2)
  return(dt)
}) %>% rbindlist() %>% arrange(type, .by_group = TRUE) %>% 
  distinct(`Cell ID`, sample_id, .keep_all = TRUE) %>% select(`Cell ID`, sample_id, type) %>% group_by(sample_id) %>% arrange(`Cell ID`, .by_group = TRUE) %>% ungroup()

dir_sp_new <- "/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/New_xenium_sample_section/0066266/"
files <- list.files(dir_sp_new, pattern="*_cells_stats.csv")
dt_sp_new <- lapply(files, function(f){
  f_path <- paste0(dir_sp_new, f)
  dt <- fread(f_path)
  dt$sample_id <- strsplit(f, split = "_") %>% sapply('[[', 1)
  return(dt)
}) %>% rbindlist() %>% group_by(sample_id) %>% arrange(`Cell ID`, .by_group = TRUE) %>% ungroup()

dt_diff <- dt_sampleinfo %>% select(sample_id, diff_level, differentiation, granuloma)
dt_sp <- left_join(dt_sp_new, dt_sp_old, by = c("Cell ID", "sample_id")) %>% left_join(dt_diff, by="sample_id")

srat <- LoadXenium("/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/data/0066266/outs/", segmentations = "cell", flip.xy = TRUE)
srat <- subset(srat, cells = dt_sp$`Cell ID`)
dt_sp <- dt_sp[match(colnames(srat), dt_sp$`Cell ID`), ] %>% select(-"Cell ID") %>% as.data.frame()
rownames(dt_sp) <- colnames(srat)
dt_sp <- dt_sp %>% as.data.frame()
srat <- AddMetaData(srat, metadata = dt_sp)
srat$slide <- "slide2"
srat$cell_ID <- paste(srat$slide, colnames(srat), sep='-')

saveRDS(srat, glue("{rds_dir}/sv5_xenium_object_0066266.rds"))
dim(srat) # 5001,1260778

#for (subsample in slide2_sample_id){
#  subsrat <- subset(srat, subset=sample_id==subsample)
#  counts_matrix <- GetAssayData(object = subsrat, assay = "Xenium", slot = "counts")
#  write_matrix_dir(mat = counts_matrix, dir = glue("{rds_dir}/matrix/{subsample}_matrix"), overwrite = TRUE)
#}
srat <- readRDS(glue("{rds_dir}/sv5_xenium_object_0066253.rds"))
meta_slide1 <- srat@meta.data
slide1_sample_id <- unique(srat$sample_id)
srat <- readRDS(glue("{rds_dir}/sv5_xenium_object_0066266.rds"))
meta_slide2 <- srat@meta.data
slide2_sample_id <- unique(srat$sample_id)

combined_meta <- rbind(meta_slide1, meta_slide2)
combined_meta <- combined_meta %>% mutate(area = `Area (µm^2)`) %>% 
    mutate(nCount_per_area = nCount_Xenium/area)
combined_meta$cell_id <- paste(rownames(combined_meta), combined_meta$sample_id, sep = "_")
saveRDS(combined_meta, glue("{rds_dir}/combined_meta.rds"))

dim(meta_slide1) # 1397413,22
dim(meta_slide2) # 1260778,22
dim(combined_meta) # 2658191,25

#list_sample_id <- append(slide1_sample_id,slide2_sample_id)
#countlist = list()
#for (subsample in list_sample_id){
#  mat <- open_matrix_dir(dir = glue("{rds_dir}/matrix/{subsample}_matrix"))
#  colnames(mat) <- paste(colnames(mat),subsample,sep = "_")
#  countlist[[subsample]] = mat
#}
#srat <- CreateSeuratObject(countlist, assay = "Xenium")

srat <- readRDS("/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/xenium/sketch/created_srat_by_sampleid.rds")
srat@meta.data <- srat@meta.data[, c("orig.ident", "nCount_Xenium", "nFeature_Xenium")]
srat$sample_id <- sapply(strsplit(colnames(srat),"_"),function(X){return(X[2])})

combined_meta <- readRDS(glue("{rds_dir}/combined_meta.rds"))
combined_meta_add <- combined_meta[match(colnames(srat), combined_meta$cell_id), ] %>% 
  select(-c(orig.ident, nCount_Xenium, nFeature_Xenium, sample_id)) %>% as.data.frame()
rownames(combined_meta_add) <- colnames(srat)
srat <- AddMetaData(srat, metadata = combined_meta_add)
saveRDS(srat, file = glue("{rds_dir}/created_srat_by_sampleid.rds"))
dim(srat) # 5001,2658191


file_anno <- "/cluster/home/tanghl_jh/projects/stomatology/docs/Total_cell_with_mix_macro_epi_subtypes_v1.csv"
dt_anno <- read_csv(file_anno)
dt_anno <- dt_anno[, c("cell_ID", "cell_ID_slide", "slide", "cell_type_reanno_sub", "mix.1.sub.celltype", "macro.1.sub.celltype", "epi.1.sub.celltype")]
names(dt_anno) <- c("cell_id", "cell_ID", "slide", "cell_type", "mix.1.subtype", "Macrophage.1.subtype", "Epithelial.1.subtype")
dt_anno <- dt_anno %>% 
    mutate(Macrophage.1.subtype = case_when(Macrophage.1.subtype == "TAM" ~ "MGC", 
                                            Macrophage.1.subtype == "mix" ~ NA_character_, 
                                            TRUE ~ Macrophage.1.subtype)) %>%
    mutate(Epithelial.1.subtype = case_when(mix.1.subtype == "mix_corneocyte" ~ "CCND2+SFRP1",
                                            Epithelial.1.subtype == "Invasive" ~ "basal_invasive",
                                            Epithelial.1.subtype == "Invasive_KRAS+" ~ "basal_invasive_KRAS+",
                                            TRUE ~ Epithelial.1.subtype)) %>%
    mutate(cell_type = case_when(cell_type == "Tcell" ~ "T cell", cell_type == "Bcell" ~ "B cell",
                                 cell_type == "Epithelial" & is.na(Epithelial.1.subtype) ~ NA_character_,
                                 cell_type == "Macrophage" & is.na(Macrophage.1.subtype) ~ NA_character_,
                                 TRUE ~ cell_type))
dt_anno$Macrophage.sub.supply <- coalesce(dt_anno$Macrophage.1.subtype, dt_anno$cell_type)
dt_anno$Epithelial.sub.supply <- coalesce(dt_anno$Epithelial.1.subtype, dt_anno$cell_type)
dt_anno <- dt_anno %>% select(cell_id, cell_ID, slide, cell_type, Macrophage.1.subtype, Epithelial.1.subtype, Macrophage.sub.supply, Epithelial.sub.supply)
saveRDS(dt_anno, file = glue("{rds_dir}/celltype_anno.rds"))


## process final srat
srat <- readRDS("/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/xenium/sketch/xenium_celltyped.rds")
srat@meta.data <- srat@meta.data[, c("orig.ident", "nCount_Xenium", "nFeature_Xenium", "sketch_snn_res.0.5", "cluster_0.5")]
dim(srat) # 5001,2436100
dt_celltype <- readRDS(glue("{rds_dir}/celltype_anno.rds")) # 2436100, 8
dt_clinical <- readRDS(glue("{rds_dir}/combined_meta.rds")) # 2658191,25

dt_clinical <- dt_clinical[dt_clinical$cell_id %in% colnames(srat),]
dt_clinical <- dt_clinical[match(colnames(srat), dt_clinical$cell_id), ] %>% select(-c(orig.ident, nCount_Xenium, nFeature_Xenium)) %>% as.data.frame()

dt_meta_add <- left_join(dt_clinical, dt_celltype, by = c("cell_id", "slide", "cell_ID"))
dt_meta_add <- dt_meta_add[match(colnames(srat), dt_meta_add$cell_id), ]
rownames(dt_meta_add) <- colnames(srat)
srat <- AddMetaData(srat, metadata = dt_meta_add) # 5001,2436100
saveRDS(srat, glue("{rds_dir}/xenium_sketch_celltyped.rds"))

dt_meta <- srat@meta.data 
saveRDS(dt_meta, glue("{rds_dir}/celltyped_meta.rds"))









### preprocess spe rds

dir <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/xenium/counts/0066253/outs"
spe <- readXeniumSXE(dir)
dim(spe) # 8313,1398220

dir_sp <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/xenium/counts/sample_section/0066253/"
files <- list.files(dir_sp, pattern = "_cells_stats.csv")
dt_sp <- lapply(files, function(f){
  f_path <- paste0(dir_sp, f)
  dt <- fread(f_path)
  dt$roi <- strsplit(f, split = "_") %>% sapply('[[', 1)
  return(dt)
}) %>% rbindlist()
dim(dt_sp) # 1397413,5

df_sp <- dt_sp[match(colnames(spe), `Cell ID`), ][, -"Cell ID"] %>% DataFrame() # 1398220,4
colData(spe) <- cbind(colData(spe), df_sp)
dt_anno_clinical <- sampleinfo$xenium
df_clinical <- dt_anno_clinical %>% as.data.frame() %>% .[match(spe$roi, .$sample_id),] %>% select(-sample_id) %>% DataFrame() # 1398220,7
colData(spe) <- cbind(colData(spe), df_clinical)

dt_anno_celltype <- readRDS(glue("{rds_dir}/celltype_anno.rds")) # 2436100,6
dt_anno_celltype <- dt_anno_celltype %>% filter(slide == "slide1") # 1286736,6

df_celltype <- dt_anno_celltype %>% 
    mutate(cell_id = sapply(strsplit(cell_id, "_"),function(X){return(X[1])}))
df_celltype <- df_celltype[match(colnames(spe), df_celltype$cell_id), ]
df_celltype <- df_celltype %>% select(-c(cell_id, slide, cell_ID)) %>% DataFrame()
colData(spe) <- cbind(colData(spe), df_celltype)



spe_s <- spe[, !is.na(spe$roi)]
spe_lst <- bplapply(unique(spe_s$roi), function(x) {
  spe_s[, spe_s$roi == x]
}, BPPARAM = SerialParam())
spe_new <- do.call("cbind", spe_lst)
spe_new$sample_id <- spe_new$roi
saveRDS(spe_new, glue("{rds_dir}/spe_0066253.rds"))
dim(spe_new) # 8313,139741


dir <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/xenium/counts/0066266/outs"
spe <- readXeniumSXE(dir)# 8313,1262349

dir_sp <- "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/xenium/counts/sample_section/0066266/"
files <- list.files(dir_sp, pattern = "_cells_stats.csv")
dt_sp <- lapply(files, function(f){
  f_path <- paste0(dir_sp, f)
  dt <- fread(f_path)
  dt$roi <- strsplit(f, split = "_") %>% sapply('[[', 1)
  return(dt)
}) %>% rbindlist()
dim(dt_sp) # 1260778,5

df_sp <- dt_sp[match(colnames(spe), `Cell ID`), ][, -"Cell ID"] %>% DataFrame() # 1262349,4
colData(spe) <- cbind(colData(spe), df_sp)
dt_anno_clinical <- sampleinfo$xenium
df_clinical <- dt_anno_clinical %>% as.data.frame() %>% .[match(spe$roi, .$sample_id),] %>% select(-sample_id) %>% DataFrame() # 1262349,7
colData(spe) <- cbind(colData(spe), df_clinical)

dt_anno_celltype <- readRDS(glue("{rds_dir}/celltype_anno.rds")) # 2436100,8
dt_anno_celltype <- dt_anno_celltype %>% filter(slide == "slide2") # 1149364,8

df_celltype <- dt_anno_celltype %>% 
    mutate(cell_id = sapply(strsplit(cell_id, "_"),function(X){return(X[1])}))
df_celltype <- df_celltype[match(colnames(spe), df_celltype$cell_id), ]
df_celltype <- df_celltype %>% select(-c(cell_id, slide, cell_ID)) %>% DataFrame()
colData(spe) <- cbind(colData(spe), df_celltype)


spe_s <- spe[, !is.na(spe$roi)]
spe_lst <- bplapply(unique(spe_s$roi), function(x) {
  spe_s[, spe_s$roi == x]
}, BPPARAM = SerialParam())
spe_new <- do.call("cbind", spe_lst)
spe_new$sample_id <- spe_new$roi
saveRDS(spe_new, glue("{rds_dir}/spe_0066266.rds"))
dim(spe_new) # 8313,1260778


# select region cell
files <- list.files("/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/doc/erosion_zone")
dt_roi <- lapply(files, function(f){
  sp <- strsplit(f, split = "_") %>% sapply('[[', 1)
  dt <- fread(paste0("/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/doc/erosion_zone/", f))
  dt$sample_id <- sp
  return(dt)
}) %>% rbindlist()
saveRDS(dt_roi, glue("{rds_dir}/select_erosion_region.rds"))



## umap supplyment
dt_celltype <- readRDS(glue("{rds_dir}/celltype_anno.rds"))

dt_anno <- readRDS(glue("{rds_dir}/combined_meta.rds"))
dt_anno <- dt_anno %>% select(cell_id, cell_ID)

srat <- readRDS("/cluster/home/ylxie_jh/projects/stomatology/analysis/lvjiong/human/xenium/macrophage_subset/spe_macrophage_all_anno.rds")
umap_macro <- as.data.frame(Embeddings(srat, reduction = "umap"))
colnames(umap_macro) <- c("UMAP_1", "UMAP_2")
umap_macro <- umap_macro %>% mutate(cell_id = rownames(.)) %>%
    left_join(dt_celltype, by="cell_id") %>% 
    select(UMAP_1, UMAP_2, Macrophage.1.subtype) %>%
    na.omit() %>% as.data.frame()

umap_giotto <- read_csv("/cluster/home/yliang_jh/projects/spatialRNA/xenium/oral_lvjiong/output/giotto/umap_harmony_coordinates.csv")
colnames(umap_giotto) <- c("cell_ID", "UMAP_1", "UMAP_2")
umap_giotto <- umap_giotto %>% left_join(dt_anno, by="cell_ID") %>% as.data.frame()
rownames(umap_giotto) <- umap_giotto$cell_id
umap_giotto <- umap_giotto %>% 
    left_join(dt_celltype, by="cell_id") %>% 
    select(UMAP_1, UMAP_2, cell_type) %>%
    na.omit() %>% as.data.frame()
umap_giotto <- umap_giotto %>% filter(UMAP_1 > -44, UMAP_2 > -50, !(UMAP_1 < -25 & UMAP_2 < -25)) %>% as.data.frame()

umap_epi <- read_csv("/cluster/home/tanghl_jh/projects/stomatology/analysis/oral_lvjiong/human/Spatial_omics/Xenium/Figure/Epithelial/mix_macro_epi/withcorn/Epi_umap_coord.csv")
colnames(umap_epi) <- c("cell_id", "UMAP_1", "UMAP_2", "epi_subtype")
umap_epi <- umap_epi %>% select(cell_id, UMAP_1, UMAP_2) %>% as.data.frame()
rownames(umap_epi) <- umap_epi$cell_id
umap_epi <- umap_epi %>%
    left_join(dt_celltype, by="cell_id") %>% 
    select(UMAP_1, UMAP_2, Epithelial.1.subtype) %>%
    na.omit() %>% as.data.frame()

srat <- readRDS(glue("{rds_dir}/xenium_sketch_celltyped.rds"))
umap_seurat <- as.data.frame(Embeddings(srat, reduction = "full.umap"))
colnames(umap_seurat) <- c("UMAP_1", "UMAP_2")
umap_seurat <- umap_seurat %>% mutate(cell_id = rownames(.)) %>%
    left_join(dt_celltype, by="cell_id") %>% 
    select(UMAP_1, UMAP_2, cell_type) %>% 
    na.omit() %>% as.data.frame()

umap_message <- list("umap_macro" = umap_macro, 
                      "umap_giotto" = umap_giotto,
                      "umap_epi" = umap_epi,
                      "umap_seurat" = umap_seurat)
saveRDS(umap_message, glue("{rds_dir}/umap_message.rds"))



## decouple_R message
net_d <- readr::read_rds("/cluster/home/ztao_jh/projects/be_a_rich_man/decoupleR/decoupleR.rds")
net_p <- readr::read_rds("/cluster/home/ylxie_jh/share/PROGENy_pathway_from_decoupleR.rds")

saveRDS(net_d, glue("{rds_dir}/decoupleR.rds"))
saveRDS(net_p, glue("{rds_dir}/ROGENy_pathway_from_decoupleR.rds"))






