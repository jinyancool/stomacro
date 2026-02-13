.libPaths(c(.libPaths(), "/cluster/home/yliang_jh/sbin/R/library/4.3.0"))
pkgs <- c("fs", "configr", "stringr", 
          "jhtools", "glue", "patchwork", "tidyverse", "dplyr", "Seurat", "magrittr", 
          "readxl", "writexl", "ComplexHeatmap", 
          "data.table", "ggplot2", "ggbeeswarm", "ggdendro", "dendextend", "deldir",
          "sf", "corrplot", "ggpubr", "survival", "survminer", "forestmodel", "BiocParallel", "BiocNeighbors")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
source("/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/quarto/R/func_utils.R")

wkdir <- glue("/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/quarto") %>% setwd()

set.seed(2025)



rds_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/codex"

sampleinfo <- readRDS("/cluster/home/jhuang/projects/stomatology/docs/lvjiong/sampleinfo/sampleinfo.rds")
markers_10 <- c("CD45", "CD20", "CD3e", "Granzyme B", "CD4", "CD8", "CD68", "SMA", "CD31", "Pan-Cytokeratin")


config_fn = "/cluster/home/jhuang/projects/stomatology/analysis/lvjiong/human/meta/manuscript/configs/colors.yaml"
config_list = show_me_the_colors(config_fn)
config <- read.config(config_fn)
cell_type_order <- config$cell_type_order


# process srt_split_anno.rds data
data <- readRDS(glue("{rds_dir}/codex_cell_segmentation.rds"))


data_fil <- data %>% 
  subset(CD3e < 22000 & CD4 < 8000 & CD8 < 7000 & CD14 < 6000 & CD20 < 5000 & 
         CD31 < 40000 & CD34 < 5000 & CD44 < 1600 & CD45 < 10000 & CD45RO < 4000 & 
         CD68 < 15000 & `Collagen IV` < 26000 & `E-cadherin` < 4200 & FOXP3 < 35000 & 
         `Granzyme B` < 4200 & `HLA-DR` < 30000 & IDO1 < 35000 & Ki67 < 7500 & 
         `Pan-Cytokeratin` < 20000 & `PD-1` < 18000 & `PD-L1` < 4200 & SMA < 11000 & 
         sum < 180000 & sum > 2000)


list_sp_rds <- list()
for(img in unique(data_fil$Image)){
  data <- filter(data_fil, Image == img)
  rownames(data) <- data$cell_id

  mat <- data[, markers_10] %>% as.matrix() %>% t
  rownames(mat) <- markers_10
  colnames(mat) <- data$cell_id
  srt <- CreateSeuratObject(mat, assay = "CODEX", project = img) %>% AddMetaData(metadata = select(data, !markers_10))
  
  srt <- NormalizeData(srt, normalization.method = NULL) %>% ScaleData()
  srt <- RunPCA(srt, features = markers_10)
  srt <- RunUMAP(srt, dims = 1:ncol(srt@reductions[["pca"]])) # use all PCs
  
  #mat <- as.matrix(srt[["CODEX"]]$scale.data) %>% t()
  #pg.cluster <- FastPG::fastCluster(mat, k = 25, num_threads = 40)
  #srt$pg_cluster <- pg.cluster$communities
  list_sp_rds[[img]] <- srt
}


list_sp_anno_rds <- list()
dt_anno <- readRDS(glue("{rds_dir}/celltype_anno_split3.rds"))
Images <- unique(dt_anno$Image)
for(img in Images){
  df_anno <- dt_anno[dt_anno$Image == img, ]
  srt <- list_sp_rds[[img]]
  srt$celltype <- NA
  srt$celltype <- df_anno$celltype[match(srt$cell_id, df_anno$cell_id)]
  srt$pg_cluster <- df_anno$pg_cluster[match(srt$cell_id, df_anno$cell_id)]
  Idents(srt) <- srt$celltype
  levels(srt) <- intersect(cell_type_order, unique(srt$celltype))
  list_sp_anno_rds[[img]] <- srt
}

saveRDS(list_sp_anno_rds, glue("{rds_dir}/srt_single_split_anno.rds"))

srt <- merge(list_sp_anno_rds[[1]], list_sp_anno_rds[2:length(list_sp_anno_rds)])

cols <- colnames(srt[[]])
cols <- cols[(which(cols == "Image") + 1):(which(cols == "sum") - 1)]
cols_n <- sub("[.]", "-", cols)
srt_combined <- JoinLayers(srt, assay = "CODEX")
mat1 <- srt_combined[["CODEX"]]$counts
mat2 <- srt_combined[[cols]] %>% t %>% `rownames<-`(cols_n)
mat <- rbind(mat1, mat2)
srt <- CreateSeuratObject(counts = mat, assay = "CODEX", project = "OralCancer", meta.data = select(srt_combined[[]], !cols))
srt <- NormalizeData(srt, normalization.method = NULL) %>% ScaleData(split.by = "Image")
srt <- RunPCA(srt, features = markers_10) 
srt <- RunUMAP(srt, dims = 1:ncol(srt@reductions[["pca"]]))
Idents(srt) <- srt$celltype
levels(srt) <- intersect(cell_type_order, unique(srt$celltype))

# subtype
list_markers <- list("CD4 T" = c("CD45RO", "FOXP3", "CD44", "Granzyme B", "HLA-DR", "Ki67", "PD-1"), 
                     "CD8 T" = c("CD45RO", "CD44", "Granzyme B", "HLA-DR", "Ki67", "PD-1"),
                     "Macrophage" = c("CD14", "CD44", "HLA-DR", "IDO1", "Ki67", "Vimentin", "PD-L1"), 
                     "Tumor cell" = c("CD44", "E-cadherin", "IDO1", "HLA-DR", "Vimentin", "Ki67", "PD-L1"), 
                     "Fibroblast" = c("CD44", "SMA", "IDO1", "HLA-DR", "Collagen-IV", "Ki67", "PD-L1"))


list_sub_celltype <- list()

celltypes <- c("CD4 T", "CD8 T", "Macrophage", "Tumor cell", "Fibroblast")
for (celltype in celltypes){
    marker_list <- list_markers[[celltype]]
    srt_s <- subset(srt, idents = celltype)
    srt_s <- NormalizeData(srt_s, normalization.method = NULL) %>% ScaleData(split.by = "Image")
    srt_s <- RunPCA(srt_s, features = marker_list)
    srt_s <- RunUMAP(srt_s, method = "umap-learn", metric = "correlation", dims = 1:ncol(srt_s@reductions[["pca"]])) # use all PCs

    mat <- as.matrix(srt_s[["CODEX"]]$scale.data[marker_list, ]) %>% t() %>% pmin(2.5) %>% pmax(-2.5)
    pg.cluster <- FastPG::fastCluster(mat, k = 25, num_threads = 40)
    srt_s$pg_cluster <- pg.cluster$communities
    Idents(srt_s) <- srt_s$pg_cluster
    levels(srt_s) <- as.character(0:max(srt_s$pg_cluster))

    df_anno <- readRDS(glue("{rds_dir}/subtype_match_cluster.rds"))
    df_anno <- df_anno[df_anno$celltype == celltype, ]
    srt_s$subtype <- NA
    srt_s$subtype <- df_anno$subtype[match(srt_s$pg_cluster, df_anno$pg_cluster)]
    table(srt_s$subtype, useNA = "ifany")
    Idents(srt_s) <- srt_s$subtype
    if("Ignore" %in% Idents(srt_s)) srt_s <- subset(srt_s, idents = "Ignore", invert = TRUE)
    dt <- srt_s@meta.data[, c("cell_id", "Image", "pg_cluster", "celltype", "subtype")]
    list_sub_celltype[[celltype]] <- dt
}
#combined_subtype_anno <- do.call(rbind, list_sub_celltype)
#saveRDS(combined_subtype_anno, glue("{rds_dir}/subtype_anno.rds"))


combined_subtype_anno <- readRDS(glue("{rds_dir}/subtype_anno.rds"))
srt_s <- subset(srt, idents = "T cell")
cts <- srt_s[["CODEX"]]$scale.data[c("CD4", "CD8"), ]
srt_s$celltype <- ifelse(cts[1, ] > cts[2, ], "CD4 T", "CD8 T")
srt$celltype[match(colnames(srt_s), colnames(srt))] <- srt_s$celltype
srt <- subset(srt, idents = "Ignore", invert = TRUE)

srt$subtype <- NA
srt$subtype[match(combined_subtype_anno$cell_id, colnames(srt))] <- combined_subtype_anno$subtype

idx <- srt$celltype %in% c("B cell", "Blood vessel", "NKT")
srt$subtype[idx] <- srt$celltype[idx]

saveRDS(srt, glue("{rds_dir}/srt_split_anno.rds"))


## calculate celltype CN 
srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
df <- srt[[]] %>% select(c("orig.ident", "Centroid.X.µm", "Centroid.Y.µm", "celltype"))
sps <- unique(df$orig.ident)
lst_dt_celltype <- list()
for(sp in sps){
  df_s <- df %>% filter(orig.ident == sp)
  ct_lst <- bplapply(1:nrow(df_s), nb_freq, df = df_s, distance = 20, BPPARAM = MulticoreParam(40))
  dt_ct <- rbindlist(ct_lst, fill = TRUE)
  dt_ct[, names(dt_ct) := lapply(.SD, nafill, fill = 0)]
  dt_ct$cell_id <- rownames(df_s)
  lst_dt_celltype[[sp]] <- dt_ct
}
dt <- do.call(lst_dt_celltype, rbind)
cols <- colnames(dt[, -"cell_id"])
dt[, (cols) := lapply(.SD, nafill, fill = 0), .SDcols = cols]
ct <- cols[-which(cols == "n")]
dt[, (ct) := lapply(.SD, function(x) x/n), .SDcols = ct]
dt_fil <- dt[n > 5, ] # exclude 8.0% cells

set.seed(2025)
k_res <- kmeans(dt_fil[, ..ct], centers = 10)
dt_fil$cluster <- k_res$cluster
cols_out <- c("cell_id", "cluster", ct)
srt$celltype_cn <- NA
srt$celltype_cn <- paste0("CN", dt_fil$cluster)[match(colnames(srt), dt_fil$cell_id)] %>% factor(., levels = paste0("CN", 1:10))
dt_celltype_cn <- dt_fil[, ..cols_out]

## calculate subtype CN
srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
df <- srt[[]] %>% select(c("orig.ident", "Centroid.X.µm", "Centroid.Y.µm", "subtype")) %>% na.omit()
sps <- unique(df$orig.ident)
lst_dt_subtype <- list()
for(sp in sps){
  df_s <- df %>% filter(orig.ident == sp)
  ct_lst <- bplapply(1:nrow(df_s), nb_freq, df = df_s, distance = 20, ct = "subtype", BPPARAM = MulticoreParam(40))
  dt_ct <- rbindlist(ct_lst, fill = TRUE)
  dt_ct[, names(dt_ct) := lapply(.SD, nafill, fill = 0)]
  dt_ct$cell_id <- rownames(df_s)
  lst_dt_subtype[[sp]] <- dt_ct
}
dt <- do.call(lst_dt_subtype, rbind)
cols <- colnames(dt[, -"cell_id"])
dt[, (cols) := lapply(.SD, nafill, fill = 0), .SDcols = cols]
ct <- cols[-which(cols == "n")]
dt[, (ct) := lapply(.SD, function(x) x/n), .SDcols = ct]
dt_fil <- dt[n > 5, ] # exclude 8.6% cells

set.seed(2025)
k_res <- kmeans(dt_fil[, ..ct], centers = 35)
dt_fil$cluster <- k_res$cluster
cols_out <- c("cell_id", "cluster", ct)
srt$subtype_cn <- NA
srt$subtype_cn <- paste0("CN", dt_fil$cluster)[match(colnames(srt), dt_fil$cell_id)] %>% factor(., levels = paste0("CN", 1:35))
dt_subtype_cn <- dt_fil[, ..cols_out]

saveRDS(srt, glue("{rds_dir}/srt_split_anno.rds"))

cluster_message <- list(celltype_cn = dt_celltype_cn, subtype_cn = dt_subtype_cn)
saveRDS(cluster_message, glue("{rds_dir}/codex_cluster_message.rds"))




## make pci rds


srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
df <- srt[[]] %>% select(c("orig.ident", "Centroid.X.µm", "Centroid.Y.µm")) %>% na.omit() %>% tibble::rownames_to_column("cell_id")
sps <- unique(df$orig.ident)


lst_deldir_edges <- list()
for(sp in sps){
  df_s <- df %>% filter(orig.ident == sp)
  del <- deldir(df_s$Centroid.X.µm, df_s$Centroid.Y.µm)
  edges <- del$delsgs %>% mutate(from = ind1, to = ind2) %>% select(from, to)
  
  edges <- edges %>%
    left_join(df_s %>% mutate(index = row_number()) %>% select(index, cell_id1 = cell_id), 
              by = c("from" = "index")) %>%
    left_join(df_s %>% mutate(index = row_number()) %>% select(index, cell_id2 = cell_id), 
              by = c("to" = "index"))
  lst_deldir_edges[[sp]] <- edges
}

saveRDS(lst_deldir_edges, glue("{rds_dir}/pci_deldir_edgs.rds"))

lst_deldir_edges <- readRDS(glue("{rds_dir}/pci_deldir_edgs.rds"))
df_subtype <- srt[[c("cell_id", "subtype")]]
df_pci2 <- lapply(sps, function(sp){
  edge <- lst_deldir_edges[[sp]]
  edge$cluster1 <- df_subtype$subtype[match(edge$cell_id1, df_subtype$cell_id)]
  edge$cluster2 <- df_subtype$subtype[match(edge$cell_id2, df_subtype$cell_id)]
  df_pci <- pci(edge) #utils.function
  df_pci$sample_id <- sp
  return(df_pci)
}) %>% rbindlist()

spinfo <- sampleinfo$codex

fill_val <- min(df_pci2$log2_lr) - 1
df_surv2 <- df_pci2 %>%
  tidyr::pivot_wider(id_cols = sample_id, names_from = pair, values_from = log2_lr, values_fill = fill_val) %>%
  left_join(spinfo) %>%
  filter(Type == "Tumor")

res_lst2 <- lapply(unique(df_pci2$pair), function(pair){
  formula <- sprintf("Surv(Time, Status) ~ `%s`", pair) %>%
    as.formula()
  cfit <- coxph(formula, data = df_surv2)
  hr <- exp(coef(cfit))
  ci <- exp(confint(cfit))
  pval <- summary(cfit)$coefficients[, "Pr(>|z|)"]
  res <- data.frame(var = pair,
                    hr = hr,
                    ci1 = ci[, 1],
                    ci2 = ci[, 2],
                    pval = pval)
  return(res)
})
df_res2 <- rbindlist(res_lst2)


lst_pci_message <- list(df_pci = df_pci2, 
                        df_surv = df_surv2,
                        df_res = df_res2)

saveRDS(lst_pci_message, glue("{rds_dir}/pci_message.rds"))



## make subtype heatmap of cor and pci
srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
sp_tumor <- sampleinfo$codex %>% filter(Type == "Tumor") %>% pull(sample_id) %>% unique
srt <- subset(srt, subset = !is.na(subtype) & Image %in% sp_tumor)

meta <- srt@meta.data %>% select(Image, subtype)

df <- meta %>% group_by(Image, subtype) %>% summarise(count = n(), .groups = 'drop')
df_total <- df %>% group_by(subtype) %>% summarise(Image = "total",count = n(),.groups = 'drop')
df <- rbind(df, df_total) %>% group_by(Image) %>% 
  mutate(proportion = count / sum(count)) %>% select(-count) %>%
  pivot_wider(names_from = subtype, values_from = proportion, values_fill = 0) %>% as.data.frame()
rownames(df) <- df$Image
df <- df %>% select(-Image)
df_cor <- cor(df, method = "pearson")

#saveRDS(cor_matrix, glue("{rds_dir}/cor_subtype.rds"))

pci_message <- readRDS(glue("{rds_dir}/pci_message.rds"))
df_pci <- pci_message$df_pci
df_pci <- df_pci[, .(
  mean_lr = mean(lr, na.rm = TRUE),
  sd_lr = sd(lr, na.rm = TRUE),
  se_lr = sd(lr, na.rm = TRUE) / sqrt(.N),
  mean_relative_freq = mean(relative_freq, na.rm = TRUE),
  mean_log2_lr = mean(log2_lr, na.rm = TRUE),
  n_samples = .N,
  sample_list = paste(unique(sample_id), collapse = ";")
), by = .(pair, cluster1, cluster2)] %>%
  select(cluster1, cluster2, mean_log2_lr) %>% as.data.frame()
df_pci <- df_pci %>% na.omit()



df_pci_long <- df_pci %>% rbind(df_pci %>% 
          rename(cluster1 = cluster2, cluster2 = cluster1) %>%
          select(cluster1, cluster2, mean_log2_lr)) %>% 
          distinct(cluster1, cluster2, .keep_all = TRUE)
dt_pci_wide <- as.data.table(df_pci_long) %>% dcast(., cluster1 ~ cluster2, value.var = "mean_log2_lr") %>% column_to_rownames("cluster1") %>% as.matrix()

dist_mat <- dist(dt_pci_wide, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
order <- colnames(dt_pci_wide)[hc$order]

df_cor <- df_cor[order, order]
df_cor_long <- df_cor %>% as.data.frame() %>% rownames_to_column("cluster1") %>% 
    pivot_longer(cols = -cluster1,names_to = "cluster2",values_to = "correlation")

df_pci_wide <- as.data.table(df_pci_long) %>% dcast(., cluster1 ~ cluster2, value.var = "mean_log2_lr") %>% column_to_rownames("cluster1") %>% as.data.frame()
df_pci_wide <- df_pci_wide[order, order]


df_cor_long <- df_cor_long %>% mutate(cluster1 = factor(.$cluster1, levels = order),
                                      cluster2 = factor(.$cluster2, levels = order))
df_pci_long <- df_pci_long %>% mutate(cluster1 = factor(.$cluster1, levels = order),
                                      cluster2 = factor(.$cluster2, levels = order))


subtype_cor_pci <- list(df_cor = df_cor, 
                        df_cor_long = df_cor_long, 
                        df_pci = df_pci_wide,
                        df_pci_long = df_pci_long)

saveRDS(subtype_cor_pci, glue("{rds_dir}/subtype_cor_pci_result.rds"))





## make voronoi celltype rds

lst_srt <- readRDS(glue("{rds_dir}/srt_single_split_anno.rds"))
dt_celltype <- readRDS(glue("{rds_dir}/codex_celltype_final.rds"))
dt_celltype <- dt_celltype %>% select(cell_id, celltype)
voronoi_fig <- list()
Images <- c("3-2-4", "3-4-3", "3-5-4")
Images <- names(lst_srt)

lst_sf_invert <- list()
for(img in Images){

  srt <- lst_srt[[img]]

  # voronoi plot
  df_ct <- srt[[]][, c("cell_id", "Image", "pg_cluster", "Centroid X µm", "Centroid Y µm")] %>% 
    left_join(dt_celltype, by = "cell_id") %>%
    mutate(celltype = factor(celltype, levels = intersect(cell_type_order, unique(srt$celltype))))
  
  sf_pt <- df_ct %>% na.omit() %>% st_as_sf(coords = c("Centroid X µm", "Centroid Y µm"), crs = NA)
  bbox <- bb_circle(sf_pt)
  
  # voronoi 
  sf_voronoi <- st_union(sf_pt) %>%
    st_voronoi() %>%
    st_collection_extract("POLYGON")%>%
    st_as_sf() %>%
    st_join(sf_pt) %>%
    st_intersection(bbox)
  
  # reverse y-axis
  scale_mat <- matrix(c(1, 0, 0, -1), nrow = 2)
  sf_invert <- sf_voronoi
  st_geometry(sf_invert) <- st_geometry(sf_invert) * scale_mat
  
  lst_sf_invert[[img]] <- sf_invert
}

saveRDS(lst_sf_invert, glue("{rds_dir}/sf_invert_insitu_celltype_message.rds"))



## make voronoi subtype rds


lst_srt <- readRDS(glue("{rds_dir}/srt_single_split_anno.rds"))
srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))
dt_celltype <- srt@meta.data %>% select(cell_id, celltype, subtype)
dt_CD44subtype <- dt_celltype %>% 
    mutate(CD44_subtype = case_when(subtype %in% c("CD4T_CD44+", "CD8T_CD44+", "Macrophage_CD44+") ~ subtype,
                                    TRUE ~ celltype)) %>%
    select(cell_id, CD44_subtype)


voronoi_fig <- list()
Images <- c("3-4-3", "3-12-4")
Images <- names(lst_srt)

lst_sf_invert <- list()
for(img in Images){

  srt <- lst_srt[[img]]

  # voronoi plot
  df_ct <- srt[[]][, c("cell_id", "Image", "Centroid X µm", "Centroid Y µm")] %>% 
    left_join(dt_CD44subtype, by = "cell_id") %>%
    mutate(CD44_subtype = factor(CD44_subtype, levels = intersect(cell_type_order, unique(dt_CD44subtype$CD44_subtype))))
  
  sf_pt <- df_ct %>% na.omit() %>% st_as_sf(coords = c("Centroid X µm", "Centroid Y µm"), crs = NA)
  bbox <- bb_circle(sf_pt)
  
  # voronoi 
  sf_voronoi <- st_union(sf_pt) %>%
    st_voronoi() %>%
    st_collection_extract("POLYGON")%>%
    st_as_sf() %>%
    st_join(sf_pt) %>%
    st_intersection(bbox)
  
  # reverse y-axis
  scale_mat <- matrix(c(1, 0, 0, -1), nrow = 2)
  sf_invert <- sf_voronoi
  st_geometry(sf_invert) <- st_geometry(sf_invert) * scale_mat
  
  lst_sf_invert[[img]] <- sf_invert
}

saveRDS(lst_sf_invert, glue("{rds_dir}/sf_invert_insitu_CD44_subtype_message.rds"))





## make voronoi subtype rds


srt <- readRDS(glue("{rds_dir}/srt_split_anno.rds"))

lst_voronoi_subtype_cn_freq <- list()

# samples subtype_cn_cluster5
df_cl <- read.csv(glue("{rds_dir}/sample_cluster_heatmap_subtype_cn_k10.csv")) %>% rename(Image = X)
imgs <- df_cl %>% filter(cluster == 5) %>% pull(Image)
srt_s <- subset(srt, subset = Image %in% imgs)

# cn type near cn25 and cn30
df <- srt_s[[]] %>% select(c("orig.ident", "Centroid.X.µm", "Centroid.Y.µm", "subtype_cn")) %>% na.omit()
cns <- paste0("CN", c(25, 30))
for(img in imgs){
  df_s <- df %>% filter(orig.ident == img)
  idx <- which(df_s$subtype_cn %in% cns)
  ct_lst <- bplapply(idx, nb_freq, df = df_s, distance = 15, ct = "subtype_cn", BPPARAM = MulticoreParam(40))
  dt_ct <- rbindlist(ct_lst, fill = TRUE)
  dt_ct[, names(dt_ct) := lapply(.SD, nafill, fill = 0)]
  dt_ct$cell_id <- rownames(df_s[idx, ])
  lst_voronoi_subtype_cn_freq[[img]] <- dt_ct
}



dt_freq <- lapply(imgs, function(img){
  dt <- lst_voronoi_subtype_cn_freq[[img]]
  cols <- colnames(dt[, -c("cell_id")])
  dt[, (cols) := lapply(.SD, nafill, fill = 0), .SDcols = cols]
  ct <- cols[-which(cols == "n")]
  dt[, (ct) := lapply(.SD, function(x) x/n), .SDcols = ct]
}) %>% rbindlist(fill = TRUE)

dt_freq[, Image := tstrsplit(cell_id, split = "_")[[1]]]
dt_freq_sum <- dt_freq[, -c("cell_id", "n")] %>%
  melt(id.vars = "Image", variable.name = "subtype_cn", value.name = "freq") %>%
  .[, .(freq_sum = sum(freq)), by = .(subtype_cn, Image)]
dt_top <- dt_freq_sum[order(-freq_sum), ][, .SD[1:10, ], by = Image] %>%
  .[, .N, by = subtype_cn] %>%
  .[order(-N), ] %>%
  .[1:10, ]
cn_top <- dt_top[, subtype_cn] %>% sort


lst_sf_invert <- list()
for(img in imgs){
  srt_ss <- subset(srt_s, subset = Image == img)
  
  # voronoi plot
  df_ct <- srt_ss[[]][, c("cell_id", "Image", "subtype_cn", "Centroid.X.µm", "Centroid.Y.µm")] %>%
    mutate(subtype_cn = factor(subtype_cn, levels = names(config_list$CN_subtype))) %>%
    mutate(subtype_cn = ifelse(subtype_cn %in% cn_top, as.character(subtype_cn), "other"))
  sf_pt <- df_ct %>% na.omit() %>% st_as_sf(coords = c("Centroid.X.µm", "Centroid.Y.µm"), crs = NA)
  bbox <- bb_circle(sf_pt)
  
  # voronoi 
  sf_voronoi <- st_union(sf_pt) %>%
    st_voronoi() %>%
    st_collection_extract("POLYGON")%>%
    st_as_sf() %>%
    st_join(sf_pt) %>%
    st_intersection(bbox)
  
  # reverse y-axis
  scale_mat <- matrix(c(1, 0, 0, -1), nrow = 2)
  sf_invert <- sf_voronoi
  st_geometry(sf_invert) <- st_geometry(sf_invert) * scale_mat

  lst_sf_invert[[img]] <- sf_invert
  }

saveRDS(lst_sf_invert, glue("{rds_dir}/sf_invert_insitu_subtype_message.rds"))







