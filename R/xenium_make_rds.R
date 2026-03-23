pkgs <- c("fs", "configr", "stringr", "scider",
          "jhtools", "glue", "patchwork", "tidyverse", "dplyr", "Seurat", "magrittr", "SeuratDisk", "readr",
          "readxl", "writexl", "ComplexHeatmap", "SpatialExperiment", "SpatialExperimentIO", "imcRtools",
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




### sketch 
srat <- readRDS(glue("{rds_dir}/created_srat_by_sampleid.rds"))
srat <- subset(srat, subset = nCount_Xenium >= 30)
srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 1000) 
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

srat <- SketchData(srat, assay="Xenium",ncells = 5000,method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(srat), seed=12387)
srat_sketch <- srat
DefaultAssay(srat) <- "sketch"

srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
srat <- ScaleData(srat, verbose = FALSE)
srat <- RunPCA(srat, verbose = FALSE)

srat <- RunUMAP(srat, reduction = "pca",dims = 1:30, 
                return.model = TRUE, verbose = FALSE, reduction.name = 'umap_bfHarmony')

srat <- IntegrateLayers(object = srat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",dims = 1:30, group.by="slide")

srat <- RunUMAP(srat, reduction = "harmony", dims = 1:30, 
                return.model = TRUE, verbose = FALSE, reduction.name = "umap")
srat <- FindNeighbors(srat, reduction = "harmony", dims = 1:30)
srat <- FindClusters(srat, resolution = 0.5, verbose = FALSE)

reso = 0.5
col = glue("sketch_snn_res.{reso}")
srat[[col]] <- factor(srat[[col]][,1],levels = as.character(sort(as.numeric(levels(srat[[col]][,1])))))

srat <- ProjectIntegration(object = srat, sketched.assay = "sketch", assay = "Xenium", 
                            features = VariableFeatures(srat[["sketch"]]), reduction = "harmony")
srat <- ProjectData(object = srat, assay = "Xenium",
  full.reduction = "harmony.full",
  sketched.assay = "sketch",
  sketched.reduction = "harmony.full",
  umap.model = "umap",
  dims = 1:30,
  refdata = list(cluster_0.5 = "sketch_snn_res.0.5"))

col = glue("cluster_{reso}")
srat[[col]] <- factor(srat[[col]][,1],levels = as.character(sort(as.numeric(levels(srat[[col]][,1])))))

srat0 <- srat
DefaultAssay(srat) <- "sketch"
srat[["sketch"]] <- JoinLayers(srat[["sketch"]])

clust_name <- glue("sketch_snn_res.{reso}")
Idents(srat) <- srat@meta.data[[clust_name]]
mks <- Seurat::FindAllMarkers(srat_tmp, assay = 'sketch', only.pos = TRUE, min.pct = 0.01,logfc.threshold = 0.2) %>% tibble() %>% mutate(cluster_name = clust_name)


dt_anno <- readRDS(glue("{rds_dir}/celltype_anno.rds"))
dt_anno <- dt_anno[colnames(srat), ]
srat$cell_type <- dt_anno$cell_type
srat <- subset(srat, subset = (sample_id != "F14") & (sample_id != "F4-2") )
srat_total <- srat


epi_srat <- subset(srat, subset = celltype == "Epithelial")
DefaultAssay(epi_srat) <- "Xenium"
epi_srat <- DietSeurat(epi_srat, assays = c("Xenium"), layers = c("counts"), dimreducs = c("harmony.full", "full.umap"))
epi_srat[["Xenium3"]] <- as(object = epi_srat[["Xenium"]], Class = "Assay")
DefaultAssay(epi_srat) <- "Xenium3"
epi_srat[["Xenium"]] <- NULL
epi_srat <- RenameAssays(object = epi_srat, Xenium3 = "Xenium")

SaveH5Seurat(epi_srat, filename = glue("{rds_dir}/Epithelial_xenium_diet.h5Seurat"))
Convert(glue("{rds_dir}/Epithelial_xenium_diet.h5Seurat"), dest = "h5ad")


spe_mac <- subset(srat, subset = celltype == "Macrophage")
DefaultAssay(spe_mac) <- "Xenium"
spe_mac <- SCTransform(spe_mac, assay = "Xenium", variable.features.n = 1000)
spe_mac <- RunPCA(spe_mac, verbose = FALSE)

spe_mac <- IntegrateLayers(object = spe_mac, method = HarmonyIntegration, normalization.method = "harmony", verbose = F)
spe_mac <- RunUMAP(spe_mac, reduction = "harmony", dims = 1:30,
                   return.model = TRUE, verbose = FALSE, reduction.name = "umap")
spe_mac <- RunTSNE(spe_mac, reduction = "harmony", dims = 1:30,
                   verbose = FALSE, reduction.name = "tsne")

spe_mac <- FindNeighbors(spe_mac, reduction = "harmony", dims = 1:30)
spe_mac <- FindClusters(spe_mac, resolution = 1, verbose = FALSE)
for (i in 1:40) {
  slot(object = spe_mac@assays$SCT@SCTModel.list[[i]], name="umi.assay") <- "Xenium"
}
spe_mac <- PrepSCTFindMarkers(object = spe_mac, assay = "SCT", verbose = TRUE)
all_markers <- FindAllMarkers(object = spe_mac, only.pos = TRUE)



dt_anno <- readRDS(glue("{rds_dir}/celltype_anno.rds"))
dt_anno <- dt_anno[colnames(srat), ]
srat$cell_type <- dt_anno$cell_type
srat$Macrophage.1.subtype <- dt_anno$Macrophage.1.subtype
srat$Epithelial.1.subtype <- dt_anno$Epithelial.1.subtype

saveRDS(srat, glue("{rds_dir}/xenium_sketch_celltyped.rds"))



##  celltype neighborhood
df <- lapply(c("0066253", "0066266"), function(slide){
  f <- glue("{rds_dir}/spe_{slide}.rds")
  spe <- readRDS(f)
  spe_s <- spe[, !is.na(spe$cell_type) & !is.na(spe$diff_level)]
  spe_s <- buildSpatialGraph(spe_s, 
                             coords = spatialCoordsNames(spe_s),
                             img_id = "roi", type = "knn", k = 10)
  spe_s <- aggregateNeighbors(spe_s, 
                              colPairName = "knn_interaction_graph", 
                              aggregate_by = "metadata", count_by = "cell_type")
  df_n <- spe_s$aggregatedNeighbors
  df_n$cell_id <- paste0(colnames(spe_s), "_", spe_s$roi)
  df_n <- df_n %>% as.data.frame()
  return(df_n)
}) %>% rbindlist() %>% as.data.frame()
rownames(df) <- df$cell_id
df <- df %>% dplyr::select(-cell_id)

k_res <- kmeans(df, centers = 10)
dt_cn <- data.table(cell_id = rownames(df),
                    celltype_cn = paste0("CN", k_res$cluster) %>% 
                    factor(levels = paste0("CN", 1:10)))

dt_anno <- readRDS(glue("{rds_dir}/celltyped_meta.rds")) 
dt_anno <- dt_anno %>% filter(!is.na(cell_type) & !is.na(diff_level)) %>% left_join(dt_cn, by = "cell_id")
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
  spe_s <- spe[, !is.na(spe$Epithelial.sub.supply) & !is.na(spe$diff_level)]
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
epi_cell <- dt_anno %>% filter(!is.na(Epithelial.1.subtype) & !is.na(diff_level)) %>% pull(cell_id)
df <- df[df$cell_id %in% epi_cell, ]

df <- df %>% dplyr::select(-cell_id)
k_res <- kmeans(df, centers = 10)
dt_cn <- data.table(cell_id = rownames(df),
                    subtype_cn = paste0("CN", k_res$cluster) %>% 
                    factor(levels = paste0("CN", 1:10)))
dt_anno <- readRDS(glue("{rds_dir}/celltyped_meta.rds")) 
dt_anno <- dt_anno %>% filter(!is.na(Epithelial.1.subtype) & !is.na(diff_level)) %>% left_join(dt_cn, by = "cell_id")
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
  nbs <- findNearCells(sqe, k = 100)
  mtx <- scanHoods(nbs$distance)      
  df_grp <- mergeByGroup(mtx, nbs$cells) %>% as.data.frame()
  return(df_grp)
}, BPPARAM = SerialParam())
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
  nbs <- findNearCells(sqe, k = 100)
  mtx <- scanHoods(nbs$distance)      
  df_grp <- mergeByGroup(mtx, nbs$cells) %>%
    as.data.frame()
  return(df_grp)
}, BPPARAM = SerialParam())
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
for (grp_diff in unique(sqe$diff_level)){
  sqe_s <- sqe[, sqe$diff_level == grp_diff]
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







## immune mediated tumor erosion zone region
dt_roi <- readRDS(glue("{rds_dir}/select_erosion_region.rds"))
dt_roi[, cell_id := paste0(`Cell ID`, "_", sample_id)]

srt <- readRDS(glue("{rds_dir}/xenium_sketch_celltyped.rds"))
srt <- subset(srt, subset = sample_id %notin% c("F4-2", "F14"))
srt$subtype <- coalesce(srt$Epithelial.1.subtype, srt$Macrophage.1.subtype)

sps <- unique(dt_roi$sample_id)
srt_s <- subset(srt, subset = sample_id %in% sps)
DefaultAssay(srt_s) <- "Xenium"
srt_s <- JoinLayers(srt_s)
srt_s$region <- "non-erosion"
srt_s$region[which(colnames(srt_s) %in% dt_roi$cell_id)] <- "erosion"

celltype <- srt_s$cell_type
celltype[is.na(celltype)] <- "Unknown"
Idents(srt_s) <- celltype

saveRDS(srt_s, glue("{rds_dir}/srt_erosion_region.rds"))
srt_s <- readRDS(glue("{rds_dir}/srt_erosion_region.rds"))

# macrophage diff gene
markers_mpg <- FindMarkers(srt_s, ident.1 = "erosion", group.by = "region", subset.ident = "Macrophage")
dt_mpg <- markers_mpg %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
  mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
  as.data.table(keep.rownames = "gene")

gene_entrez <- bitr(unique(dt_mpg$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dt_mpg <- merge(dt_mpg, gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")

# epithelial diff gene
markers_epi <- FindMarkers(srt_s, ident.1 = "erosion", group.by = "region", subset.ident = "Epithelial")
dt_epi <- markers_epi %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
  mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
  as.data.table(keep.rownames = "gene")

gene_entrez <- bitr(unique(dt_epi$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dt_epi <- merge(dt_epi, gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")

erosion_region_diff_gene <- list(
  "Macrophage_diff_gene" = markers_mpg,
  "Macrophage_diff_gene_sig" = dt_mpg,
  "Epithelial_diff_gene" = markers_epi,
  "Epithelial_diff_gene_sig" = dt_epi
)
saveRDS(erosion_region_diff_gene, glue("{rds_dir}/erosion_region_diff_gene.rds"))


# run clusterprofiler: macrophage & epithelial in erosion
erosion_region_diff_gene <- readRDS(glue("{rds_dir}/erosion_region_diff_gene.rds"))

dt_mpg <- erosion_region_diff_gene$Macrophage_diff_gene_sig
map_res <- compareCluster(ENTREZID ~ type,
                      data = dt_mpg,
                      organism = "hsa",
                      fun = "enrichKEGG",
                      use_internal_data = TRUE)
#geneLst <- dt_mpg[dt_mpg$type == "up",]$ENTREZID
#res <- enrichKEGG(geneLst, organism = "hsa", 
#                  keyType = "kegg", pvalueCutoff = 0.05, use_internal_data = TRUE)
#res <- enrichGO(geneLst, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID",pvalueCutoff = 0.05)
mpg_res_readable <- setReadable(map_res, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")


dt_epi <- erosion_region_diff_gene$Epithelial_diff_gene_sig
epi_res <- compareCluster(ENTREZID ~ type,
                      data = dt_epi,
                      organism = "hsa",
                      fun = "enrichKEGG",
                      use_internal_data = TRUE)
epi_res_readable <- setReadable(epi_res, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")


erosion_region_result <- list(
    "Macrophage_enrich_kegg_res_plot" = map_res,
    "Macrophage_enrich_kegg_res_read" = mpg_res_readable,
    "Epithelial_enrich_kegg_res_plot" = epi_res,
    "Epithelial_enrich_kegg_res_read" = epi_res_readable
    )

saveRDS(erosion_region_result, glue("{rds_dir}/erosion_region_result.rds"))




## immune mediated tumor erosion zone sample
srt <- readRDS(glue("{rds_dir}/xenium_sketch_celltyped.rds"))
srt$subtype <- coalesce(srt$Epithelial.1.subtype, srt$Macrophage.1.subtype)
srt <- subset(srt, subset = sample_id %notin% c("F4-2", "F14"))
df_sp <- sampleinfo$xenium
sps <- df_sp %>% filter(!is.na(granuloma)) %>% pull(sample_id)
srt_s <- subset(srt, subset = sample_id %in% sps)
DefaultAssay(srt_s) <- "Xenium"
srt_s <- JoinLayers(srt_s)

group <- df_sp$granuloma[match(srt_s$sample_id, df_sp$sample_id)]
srt_s$group <- ifelse(group == 1, "granuloma", "no-granuloma")

celltype <- srt_s$cell_type
celltype[is.na(celltype)] <- "Unknown"
Idents(srt_s) <- celltype
saveRDS(srt_s, glue("{rds_dir}/srt_erosion_sample.rds"))
srt_s <- readRDS(glue("{rds_dir}/srt_erosion_sample.rds"))

# macrophage diff gene
markers_mpg <- FindMarkers(srt_s, ident.1 = "granuloma", group.by = "group", subset.ident = "Macrophage")
dt_mpg <- markers_mpg %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
  mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
  as.data.table(keep.rownames = "gene")

gene_entrez <- bitr(unique(dt_mpg$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dt_mpg <- merge(dt_mpg, gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")

# epithelial diff gene
markers_epi <- FindMarkers(srt_s, ident.1 = "granuloma", group.by = "group", subset.ident = "Epithelial")
dt_epi <- markers_epi %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
  mutate(type = ifelse(avg_log2FC > 0, "up", "down")) %>%
  as.data.table(keep.rownames = "gene")

gene_entrez <- bitr(unique(dt_epi$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dt_epi <- merge(dt_epi, gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")

erosion_sample_diff_gene <- list(
  "Macrophage_diff_gene" = markers_mpg,
  "Macrophage_diff_gene_sig" = dt_mpg,
  "Epithelial_diff_gene" = markers_epi,
  "Epithelial_diff_gene_sig" = dt_epi
)
saveRDS(erosion_sample_diff_gene, glue("{rds_dir}/erosion_sample_diff_gene.rds"))


# run clusterprofiler: macrophage & epithelial in erosion sample
erosion_sample_diff_gene <- readRDS(glue("{rds_dir}/erosion_sample_diff_gene.rds"))

dt_mpg <- erosion_sample_diff_gene$Macrophage_diff_gene_sig
map_res <- compareCluster(ENTREZID ~ type,
                      data = dt_mpg,
                      organism = "hsa",
                      fun = "enrichKEGG",
                      use_internal_data = TRUE)
mpg_res_readable <- setReadable(map_res, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")

dt_epi <- erosion_sample_diff_gene$Epithelial_diff_gene_sig
epi_res <- compareCluster(ENTREZID ~ type,
                      data = dt_epi,
                      organism = "hsa",
                      fun = "enrichKEGG",
                      use_internal_data = TRUE)
epi_res_readable <- setReadable(epi_res, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")


erosion_sample_result <- list(
    "Macrophage_enrich_kegg_res_plot" = map_res,
    "Macrophage_enrich_kegg_res_read" = mpg_res_readable,
    "Epithelial_enrich_kegg_res_plot" = epi_res,
    "Epithelial_enrich_kegg_res_read" = epi_res_readable)

saveRDS(erosion_sample_result, glue("{rds_dir}/erosion_sample_result.rds"))





## decouple R analyze tumor erosion zone sample

srt <- readRDS(glue("{rds_dir}/srt_erosion_region.rds"))
srt <- subset(srt, subset = sample_id %notin% c("F4-2", "F14"))
net_d <- readRDS(glue("{rds_dir}/decoupleR.rds"))
net_p <- readRDS(glue("{rds_dir}/ROGENy_pathway_from_decoupleR.rds"))
srt$sam_region <- str_c(srt$sample_id, "@", srt$region)
srt$sam_region <- factor(srt$sam_region, levels = sort(unique(srt$sam_region)))
srt$stype_region <- str_c(srt$subtype, "@", srt$region)
srt$stype_region <- factor(srt$stype_region, levels = sort(unique(srt$stype_region)))

## macrophage
srt_mac <- subset(srt, subset = cell_type == "Macrophage")
DefaultAssay(srt_mac) <- "Xenium"
#srt_mac <- JoinLayers(srt_mac, assay = "Xenium")
mat <- GetAssayData(srt_mac, assay = "Xenium", layer = "data")
dim(mat)
mat_dgc <- as(mat, "dgCMatrix") %>% as.matrix()

# tfs
acts_tar <- decoupleR::run_ulm(mat = mat_dgc,
                           net = net_d,
                           .source = 'source',
                           .target = 'target',
                           .mor='mor',
                           minsize = 5)
srt_mac[['tfsulm']] <- acts_tar %>%
  tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(srt_mac) <- "tfsulm"
srt_mac <- Seurat::ScaleData(srt_mac)
srt_mac@assays$tfsulm@data <- srt_mac@assays$tfsulm@scale.data

Idents(srt_mac) <- "region"
df_tfs <- t(as.matrix(srt_mac@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(region = srt_mac$region, subtype = srt_mac$subtype, samples = srt_mac$sample_id,
                sam_region = str_c(samples, "@", region), stype_region = str_c(region, "@", subtype))
df_tfs_region <- df_tfs %>% dplyr::select(-sam_region, -subtype, -samples, -stype_region) %>%
  tidyr::pivot_longer(cols = -region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(region, source) %>%
  dplyr::summarise(mean = mean(score))
df_tfs_sample_region <- df_tfs %>% dplyr::select(-region, -subtype, -samples, -stype_region) %>%
  tidyr::pivot_longer(cols = -sam_region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(sam_region, source) %>%
  dplyr::summarise(mean = mean(score))
df_tfs_stype_region <- df_tfs %>% dplyr::select(-region, -subtype, -samples, -sam_region) %>%
  tidyr::pivot_longer(cols = -stype_region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(stype_region, source) %>%
  dplyr::summarise(mean = mean(score))

n_tfs <- 100
top_region_tfs <- df_tfs_region %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)
mat_top_acts_tfs_region <- df_tfs_region %>% dplyr::filter(source %in% top_region_tfs) %>%
  tidyr::pivot_wider(id_cols = 'region', names_from = 'source', values_from = 'mean') %>%
  tibble::column_to_rownames('region') %>%
  as.matrix()


# pathway
acts_path <- decoupleR::run_ulm(mat = mat_dgc,
                               net = net_p,
                               .source = 'source',
                               .target = 'target',
                               .mor='weight',
                               minsize = 5)


srt_mac[["pathulm"]] <- acts_path %>% 
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
DefaultAssay(srt_mac) <- "pathulm"
srt_mac <- Seurat::ScaleData(srt_mac)
srt_mac@assays$pathulm@data <- srt_mac@assays$pathulm@scale.data

df_score <- t(srt_mac@assays$pathulm@data) %>%
  as_tibble(rownames = "barcode") %>%
  dplyr::mutate(region = srt_mac$region, sam_region = srt_mac$sam_region, subtype = srt_mac$subtype) %>%
  tidyr::pivot_longer(cols = Androgen:p53, names_to = "pathway", values_to = "score")
df_score_sam <- df_score %>%
  dplyr::group_by(sam_region, region, pathway) %>%
  dplyr::summarise(mean = mean(score)) %>% ungroup()



decoupleR_message_macro <- list(path = acts_path, tfs = acts_tar, 
                                path_score = df_score, path_score_sam = df_score_sam,
                                tfs_total = df_tfs, tfs_region = df_tfs_region, 
                                tfs_stype_region = df_tfs_stype_region, tfs_sample_region = df_tfs_sample_region,
                                tfs_top100_mat = mat_top_acts_tfs_region)
saveRDS(decoupleR_message_macro, glue("{rds_dir}/decoupleR_message_Macrophage.rds"))
saveRDS(srt_mac, glue("{rds_dir}/srt_decoupleR_Macrophage.rds"))




## epithelial
srt_epi <- subset(srt, subset = cell_type == "Epithelial")
DefaultAssay(srt_epi) <- "Xenium"
mat <- GetAssayData(srt_epi, assay = "Xenium", layer = "data")
dim(mat)
mat_dgc <- as(mat, "dgCMatrix") %>% as.matrix()

# tfs
acts_tar <- decoupleR::run_ulm(mat = mat_dgc,
                           net = net_d,
                           .source = 'source',
                           .target = 'target',
                           .mor='mor',
                           minsize = 5)
srt_epi[['tfsulm']] <- acts_tar %>%
  tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(srt_epi) <- "tfsulm"
srt_epi <- Seurat::ScaleData(srt_epi)
srt_epi@assays$tfsulm@data <- srt_epi@assays$tfsulm@scale.data

Idents(srt_epi) <- "region"
df_tfs <- t(as.matrix(srt_epi@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(region = srt_epi$region, subtype = srt_epi$subtype, samples = srt_epi$sample_id,
                sam_region = str_c(samples, "@", region), stype_region = str_c(region, "@", subtype))
df_tfs_region <- df_tfs %>% dplyr::select(-sam_region, -subtype, -samples, -stype_region) %>%
  tidyr::pivot_longer(cols = -region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(region, source) %>%
  dplyr::summarise(mean = mean(score))
df_tfs_sample_region <- df_tfs %>% dplyr::select(-region, -subtype, -samples, -stype_region) %>%
  tidyr::pivot_longer(cols = -sam_region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(sam_region, source) %>%
  dplyr::summarise(mean = mean(score))
df_tfs_stype_region <- df_tfs %>% dplyr::select(-region, -subtype, -samples, -sam_region) %>%
  tidyr::pivot_longer(cols = -stype_region, names_to = "source", values_to = "score") %>%
  dplyr::group_by(stype_region, source) %>%
  dplyr::summarise(mean = mean(score))

n_tfs <- 100
top_region_tfs <- df_tfs_region %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)
mat_top_acts_tfs_region <- df_tfs_region %>% dplyr::filter(source %in% top_region_tfs) %>%
  tidyr::pivot_wider(id_cols = 'region', names_from = 'source', values_from = 'mean') %>%
  tibble::column_to_rownames('region') %>%
  as.matrix()


# pathway
acts_path <- decoupleR::run_ulm(mat = mat_dgc,
                               net = net_p,
                               .source = 'source',
                               .target = 'target',
                               .mor='weight',
                               minsize = 5)


srt_epi[["pathulm"]] <- acts_path %>% 
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
DefaultAssay(srt_epi) <- "pathulm"
srt_epi <- Seurat::ScaleData(srt_epi)
srt_epi@assays$pathulm@data <- srt_epi@assays$pathulm@scale.data

df_score <- t(srt_epi@assays$pathulm@data) %>%
  as_tibble(rownames = "barcode") %>%
  dplyr::mutate(region = srt_epi$region, sam_region = srt_epi$sam_region, subtype = srt_epi$subtype) %>%
  tidyr::pivot_longer(cols = Androgen:p53, names_to = "pathway", values_to = "score")
df_score_sam <- df_score %>%
  dplyr::group_by(sam_region, region, pathway) %>%
  dplyr::summarise(mean = mean(score)) %>% ungroup()



decoupleR_message_epi <- list(path = acts_path, tfs = acts_tar, 
                                path_score = df_score, path_score_sam = df_score_sam,
                                tfs_total = df_tfs, tfs_region = df_tfs_region, 
                                tfs_stype_region = df_tfs_stype_region, tfs_sample_region = df_tfs_sample_region,
                                tfs_top100_mat = mat_top_acts_tfs_region)
saveRDS(decoupleR_message_epi, glue("{rds_dir}/decoupleR_message_Epithelial.rds"))
saveRDS(srt_epi, glue("{rds_dir}/srt_decoupleR_Epithelial.rds"))











## Epithelial subtype activated pathway
srt <- readRDS(glue("{rds_dir}/xenium_sketch_celltyped.rds"))
net_d <- readRDS(glue("{rds_dir}/decoupleR.rds"))
net_p <- readRDS(glue("{rds_dir}/ROGENy_pathway_from_decoupleR.rds"))

srt_epi <- subset(srt, subset = cell_type == "Epithelial")
DefaultAssay(srt_epi) <- "Xenium"

srt_s <- JoinLayers(srt_epi)
Idents(srt_s) <- "Epithelial.1.subtype"

markers_epi <- FindAllMarkers(srt_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
dt_epi <- markers_epi %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

gene_entrez <- bitr(unique(dt_epi$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
dt_epi <- merge(dt_epi, gene_entrez, all.x = TRUE, by.x = "gene", by.y = "SYMBOL")

epi_markers <- list(raw = markers_epi, filter = dt_epi)
saveRDS(epi_markers, glue("{rds_dir}/epithelial_cluster_marker_gene.rds"))




mat <- LayerData(srt_s, assay = "Xenium", layer = "data")
dim(mat)
mat_dgc <- as(mat, "dgCMatrix") %>% as.matrix()

# tfs
acts_tar <- decoupleR::run_ulm(mat = mat_dgc,
                           net = net_d,
                           .source = 'source',
                           .target = 'target',
                           .mor='mor',
                           minsize = 5)
srt_s[['tfsulm']] <- acts_tar %>%
  tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(srt_s) <- "tfsulm"
srt_s <- Seurat::ScaleData(srt_s)
srt_s@assays$tfsulm@data <- srt_s@assays$tfsulm@scale.data

Idents(srt_s) <- "Epithelial.1.subtype"
df_tfs <- t(as.matrix(srt_s@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(subtype = srt_s$Epithelial.1.subtype, samples = srt_s$sample_id)
df_tfs_subtype <- df_tfs %>% dplyr::select(-samples) %>% tidyr::pivot_longer(cols = -subtype, names_to = "source", values_to = "score") %>%
  dplyr::group_by(subtype, source) %>% dplyr::summarise(mean = mean(score))

n_tfs <- 100
top_tfs <- df_tfs_subtype %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)
mat_top_acts_tfs <- df_tfs_subtype %>% dplyr::filter(source %in% top_tfs) %>%
  tidyr::pivot_wider(id_cols = 'subtype', names_from = 'source', values_from = 'mean') %>%
  tibble::column_to_rownames('subtype') %>%
  as.matrix()


# pathway
acts_path <- decoupleR::run_ulm(mat = mat_dgc,
                               net = net_p,
                               .source = 'source',
                               .target = 'target',
                               .mor='weight',
                               minsize = 5)


srt_s[["pathulm"]] <- acts_path %>% 
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
DefaultAssay(srt_s) <- "pathulm"
srt_s <- Seurat::ScaleData(srt_s)
srt_s@assays$pathulm@data <- srt_s@assays$pathulm@scale.data

df_score <- t(srt_s@assays$pathulm@data) %>%
  as_tibble(rownames = "barcode") %>%
  dplyr::mutate(subtype = srt_s$Epithelial.1.subtype, samples = srt_s$sample_id) %>%
  tidyr::pivot_longer(cols = Androgen:p53, names_to = "pathway", values_to = "score")
df_score_sam <- df_score %>%
  dplyr::group_by(subtype, pathway) %>%
  dplyr::summarise(mean = mean(score)) %>% ungroup()



decoupleR_message_epi_cluster <- list(path = acts_path, tfs = acts_tar, 
                                path_score = df_score, path_score_sam = df_score_sam,
                                tfs_total = df_tfs, tfs_subtype = df_tfs_subtype,
                                tfs_top100_mat = mat_top_acts_tfs)
saveRDS(decoupleR_message_epi_cluster, glue("{rds_dir}/decoupleR_message_Epithelial_totalSubtype.rds"))
saveRDS(srt_s, glue("{rds_dir}/srt_decoupleR_Epithelial_totalSubtype.rds"))







# spe scider the celltype arounding diff level density of Epithelial 

spe <- readRDS(glue("{rds_dir}/spe_0066253.rds"))

#spe$Epithelial_supply <- coalesce(spe$Epithelial.1.subtype, spe$cell_type)
spe$Epithelial_supply <- case_when(
  is.na(spe$cell_type) ~ "Epithelial",
  spe$Epithelial.1.subtype == "mix_epi" ~ "Unkown",
  !is.na(spe$Epithelial.1.subtype) ~ "Epithelial",
  TRUE ~ spe$cell_type
)
spe$Epithelial_subsupply <- case_when(
  is.na(spe$cell_type) ~ "epi_na",
  !is.na(spe$Epithelial.1.subtype) ~ spe$Epithelial.1.subtype,
  TRUE ~ spe$cell_type
)

coi <- c("Epithelial")
coi <- unique(na.omit(spe$Epithelial_subsupply)) %>% .[. != "mix_epi"]
coi <- unique(na.omit(spe$cell_type)) %>% .[. != "Epithelial"]

samples <- unique(spe$sample_id)
for (sample in samples){
  spe_s <- spe[, spe$sample_id == sample]
  coi <- c("Epithelial")
  coi <- unique(na.omit(spe_s$cell_type)) %>% .[. != "Epithelial"]
  coi <- unique(na.omit(spe_s$Epithelial_subsupply)) %>% .[. != "mix_epi"]
  spe_s <- gridDensity(spe_s, id = "Epithelial_subsupply", grid.length.x = 50, bandwidth = 25)
  spe_s <- findROI(spe_s, coi = coi, min.density = 1)
  spe_cor <- corDensity(spe_s)
  p1 <- plotDensity(spe_s)
  p2 <- plotDensity(spe_s, coi = coi)
  p3 <- scider::plotROI(spe_s, roi = coi)
  pdf(glue("{fig_dir}/scider/xenium_density_subtype_{sample}.pdf"))
  print(p1)
  print(p2)
  print(p3)
  dev.off()

}

results <- corDensity(spe)
results$ROI
results$overall
p <- scider::plotCorHeatmap(results$ROI)
p <- plotCorHeatmap(results$overall)
spe_cor <- corDensity(spe, roi = "Fibroblasts")
p <- scider::plotCorHeatmap(spe_cor)


spe <- getContour(spe, coi = coi, bins = 10)
scider::plotContour(spe, coi = coi, line.width = 0.5)

spe <- getContour(spe, coi = coi, equal.cell = TRUE)
plotContour(spe, coi = coi)

spe <- allocateCells(spe)
plotSpatial(spe, group.by = "fibroblasts_contour", pt.alpha = 0.5)

spe <- allocateCells(spe, contour = coi, to.roi = TRUE)
scider::plotCellCompo(spe, contour = coi, self.included = FALSE)
scider::plotCellCompo(spe, contour = coi, self.included = FALSE, roi = coi)

plotCellCompo(spe, contour = "Fibroblasts", roi = "Fibroblasts")


## cell communication of stemness and differentiation Epithelial

use_python("/cluster/home/lixiyue_jh/.conda/envs/pyr/bin/python")
spinfo <- sampleinfo$xenium_epithelial_for_cc
spinfo <- spinfo %>% filter(communication %in% c(1,2)) %>% dplyr::select(Epithelial.1.subtype, group)


ccc_results <- list()


spe <- readRDS(glue("{rds_dir}/spe_0066253.rds"))
dt <- left_join(as.data.frame(colData(spe)), spinfo, by = "Epithelial.1.subtype")
spe$group <- dt$group

counts_matrix <- assay(spe, "counts")
keep_genes <- rowSums(counts_matrix>0) > 0

spe$library_size <- librarySizeFactors(spe)
min_factor <- min(spe$library_size[spe$library_size > 0]) #0.002
spe$library_size[spe$library_size <= 0] <- 0.0001
spe <- logNormCounts(spe, size.factors=spe$library_size)

ct <- import("commot")
db <- ct$pp$ligand_receptor_database(
    species="human", 
    database="CellChat", 
    signaling_type=NULL)
names(db) <- c("ligand", "receptor", "pathway", "type")

lr <- c(db$ligand, db$receptor)
ss <- strsplit(lr, "-")
gs <- sapply(ss, .subset, 1)

keep <- apply(db, 1, \(.) {
    rs <- strsplit(.["receptor"], "_")
    lr <- c(.["ligand"], unlist(rs))
    all(lr %in% rownames(spe))
})
db <- db[keep, ]

ss <- strsplit(db$receptor, "_")
rs <- sapply(ss, .subset, 1)

for (sample in unique(spe$sample_id)){
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
  ccc_results[sample] <- ccc
  for (df in ccc) colData(sub)[names(df)] <- df
  as <- lapply(ccc, \(.) {
    names(.) <- gsub("^(s|r)-", "", names(.))
    as(t(as.matrix(.)), "dgCMatrix")
  })
  sce <- SingleCellExperiment(as, colData=colData(sub))
  saveRDS(sce, glue("{rds_dir}/CCC/sce_{sample}.rds"))
}

slide1_sample_id <- unique(spe$sample_id)


spe <- readRDS(glue("{rds_dir}/spe_0066266.rds"))
dt <- left_join(as.data.frame(colData(spe)), spinfo, by = "Epithelial.1.subtype")
spe$group <- dt$group

counts_matrix <- assay(spe, "counts")
keep_genes <- rowSums(counts_matrix>0) > 0

spe$library_size <- librarySizeFactors(spe)
min_factor <- min(spe$library_size[spe$library_size > 0]) #0.002
spe$library_size[spe$library_size <= 0] <- 0.00001
spe <- logNormCounts(spe, size.factors=spe$library_size)

ct <- import("commot")
db <- ct$pp$ligand_receptor_database(
    species="human", 
    database="CellChat", 
    signaling_type=NULL)
names(db) <- c("ligand", "receptor", "pathway", "type")

lr <- c(db$ligand, db$receptor)
ss <- strsplit(lr, "-")
gs <- sapply(ss, .subset, 1)

keep <- apply(db, 1, \(.) {
    rs <- strsplit(.["receptor"], "_")
    lr <- c(.["ligand"], unlist(rs))
    all(lr %in% rownames(spe))
})
db <- db[keep, ]

ss <- strsplit(db$receptor, "_")
rs <- sapply(ss, .subset, 1)

for (sample in unique(spe$sample_id)){
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
  ccc_results[sample] <- ccc
  for (df in ccc) colData(sub)[names(df)] <- df
  as <- lapply(ccc, \(.) {
    names(.) <- gsub("^(s|r)-", "", names(.))
    as(t(as.matrix(.)), "dgCMatrix")
  })
  sce <- SingleCellExperiment(as, colData=colData(sub))
  saveRDS(sce, glue("{rds_dir}/CCC/sce_{sample}.rds"))
}

slide2_sample_id <- unique(spe$sample_id)

saveRDS(ccc_results, glue("{rds_dir}/CCC/ccc_commot_result.rds"))


results <- data.frame(
    sample_id = character(),
    column = character(),
    group1 = character(),
    group2 = character(),
    mean_group1 = numeric(),
    mean_group2 = numeric(),
    t_statistic = numeric(),
    df = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )

samples_id <- c(slide1_sample_id, slide2_sample_id) %>% .[! . %in% c("F14", "F4-2")]
for (sample in samples_id){
  sce <- readRDS(glue("{rds_dir}/CCC/sce_{sample}.rds"))
  if (all(is.na(sce$group))) { next }
  sce <- sce[,!is.na(sce$group)]
  df <-  as.data.frame(colData(sce))
  sr_cols <- grep("^(s-|r-|r[.]|s[.])", colnames(df), value = TRUE)
  df1 <- df[, sr_cols]
  sr_cols <- colnames(df1[, !sapply(df1, function(x) all(is.na(x) | x == 0))])
  group1 <- "stemness"
  group2 <- "differentiation"
  for (col_name in sr_cols){
    val_group1 <- df[df$group == group1,][[col_name]]
    val_group2 <- df[df$group == group2,][[col_name]]
    test_result <- tryCatch({
      t.test(val_group1, val_group2, var.equal = FALSE)
      #formula <- as.formula(paste(col_name, "~ group"))
      #t.test(formula, data = df, var.equal = FALSE)
      }, error = function(e) {
      return(NULL)
      })
    if (!is.null(test_result)) {
      mean1 <- mean(val_group1, na.rm = TRUE)
      mean2 <- mean(val_group2, na.rm = TRUE)
      results_add <- data.frame(
        sample_id = sample,
        column = col_name,
        group1 = group1,
        group2 = group2,
        mean_group1 = mean1,
        mean_group2 = mean2,
        t_statistic = test_result$statistic,
        df = test_result$parameter,
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, results_add)
    }
  }
}

results$p_adj <- p.adjust(results$p_value, method = "BH")
results$significance <- cut(
      results$p_adj,
      breaks = c(0, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", "NS"),
      include.lowest = TRUE
    )
results_filter <- results[results$p_adj<=0.05,]

results <- results %>% group_by(sample_id) %>% arrange(desc(abs(t_statistic)), p_adj, .by_group = TRUE) %>% ungroup() %>% as.data.frame()
results_filter <- results_filter %>% group_by(sample_id) %>% arrange(desc(abs(t_statistic)), p_adj, .by_group = TRUE) %>% ungroup() %>% as.data.frame()
results_top10 <- results_filter %>% group_by(sample_id) %>% arrange(desc(abs(t_statistic)), p_adj, .by_group = TRUE) %>% 
  slice_head(n = 10) %>% ungroup() %>% as.data.frame()


ccc_test <- list("filter" = results_filter, "raw" = results, "top10" = results_top10)
saveRDS(ccc_test, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff.rds"))



mu_list <- list()
df_list <- list()

for (sample in samples_id){
  print(sample)
  sce <- readRDS(glue("{rds_dir}/CCC/sce_{sample}.rds"))
  sce$Epithelial_ccc <- coalesce(sce$group, sce$cell_type)
  mu <- aggregateAcrossCells(sce, 
    ids=sce$Epithelial_ccc, 
    statistics="mean",
    use.assay.type=assayNames(sce))
  mu_list[[sample]] <- mu
  
  ks <- expand.grid(ks <- colnames(mu), ks)
  cs <- split(colnames(sce), sce$Epithelial_ccc)
  lr <- grepl("-", rownames(sce))
  lr <- lr & !grepl("total", rownames(sce))
  ss <- strsplit(rownames(sce), "-")
  l <- sapply(ss, .subset, 1)
  r <- sapply(ss, .subset, 2)
  df <- mapply(
    i=ks[, 1], j=ks[, 2],
    SIMPLIFY=FALSE, \(i, j) {
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
  df_list[[sample]] <- df
}

saveRDS(mu_list, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_mu_list.rds"))
saveRDS(df_list, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path_list.rds"))

df <- do.call(rbind, df_list)
rownames(df) <- NULL
df$cell_cell <- paste(df$source, "-", df$target)
df$path <- paste(df$ligand,"-",df$receptor)
df <- df %>% filter(!sample_id %in% c("F14", "F4-2")) %>% as.data.frame()
saveRDS(df, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path.rds"))


cc_path <- data.frame(path = character(),target_cell = character(),mean_stem = character(),mean_diff = character(),
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
        path = p,
        target_cell = target,
        mean_stem = mean1,
        mean_diff = mean2,
        t_statistic = test_result$statistic,
        df = test_result$parameter,
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
      cc_path <- rbind(cc_path, cc_path_add)
    }
  }
}

cc_path_filter <- cc_path %>% filter(p_value <= 0.05) %>% 
    arrange(desc(abs(t_statistic)), p_value) %>% as.data.frame()
results_cc_path <- list("filter" = cc_path_filter, "raw" = cc_path)
saveRDS(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_CC_path.rds"))
save_list_to_excel(results_cc_path, glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_CC_path.xlsx"))


select_path <- cc_path_filter %>% head(200) %>% pull(path) %>% unique()
df <- readRDS(glue("{rds_dir}/xenium_CCC_Epithelial_stem_diff_df_path.rds"))


df_mean <- df %>% group_by(source, target, ligand, receptor) %>% 
    summarise(score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), 
              n_observations = n(), n_samples = n_distinct(sample_id), .groups = 'drop') 

df_mean_select <- df_mean %>% 
      filter(source %in% c("stemness", "differentiation")) %>% 
      filter(!target %in% c("Epithelial", "stemness", "differentiation")) %>%
      filter(score > 0)
df_heat_cc <- df_mean_select %>% group_by(source, target) %>% summarise(score = mean(score), .groups = "drop")
df_heat_cp <- df_mean_select %>% group_by(source, target, ligand, receptor) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>% mutate(path = paste(ligand,"-",receptor), cell_cell = paste(source,"-",target)) %>%
    filter(path %in% select_path) %>% as.data.frame()
cell_cell_order <- c(
  paste("stemness", "-", intersect(config$cell_type_order, unique(df$target))),
  paste("differentiation", "-", intersect(config$cell_type_order, unique(df$target)))
) %>% .[. %in% unique(df_heat_cp$cell_cell)]

df_heat_cp_wide <- df_heat_cp %>% dplyr::select(cell_cell, path, score) %>%
    pivot_wider(names_from = cell_cell, values_from = score, values_fill = NA) %>%
    column_to_rownames("path") %>% as.matrix() %>% .[, cell_cell_order, drop = FALSE]


p11 <- cc_dotplot(df_mean_select, option = 'B', n_top_ints = 200) +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1))+
  scale_x_discrete(labels = function(x) gsub("&rarr;", " → ", x)) +
  scale_color_gradientn(colors = config_list$scale_7)
p1 <- cc_heatmap(df_mean_select, option = 'B', n_top_ints = 100) + 
  scale_fill_gradientn(colors = config_list$scale_7, oob = scales::squish)

p2 <- cc_heatmap(df_mean_select) + 
  scale_fill_gradientn(colors = config_list$scale_7, oob = scales::squish)

p3 <- ggplot(df_heat_cc, aes(target, source, fill = score)) + geom_tile(color = "black") +
  scale_fill_gradientn(colours = config_list$scale_7) +
  theme_classic()
p4 <- cc_circos(df_mean)
p5 <- cc_network(df_mean %>% filter(
    source %in% c("stemness", "differentiation") |
    target %in% c("stemness", "differentiation")
  ), colours = colors_16)
p6 <- ggplot(df_heat_cp, aes(cell_cell, path, fill = score)) + geom_tile(color = "black") +
  scale_fill_gradientn(colours = config_list$scale_7) +
  theme_classic()
p7 <- pheatmap(df_heat_cp_wide, scale = "row", 
  color = config_list$scale_7, border_color = "white",
  cluster_rows = TRUE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 8, fontsize_col = 8, angle_col = "90", main = "CCC"
)

pdf(glue("{fig_dir}/CCC/xenium_ccc_heatmaps_cc_mean.pdf"), width = 8, height = 7)
print(p11)
print(p7)
dev.off()


samples_id <- samples_id %>% .[! . %in% c("F14", "F4-2")]
for (sample in samples_id){
  df_s <- df_list[[sample]]
  df_s_select <- df_s %>% filter(source %in% c("stemness", "differentiation")) %>%
      filter(!target %in% c("Epithelial", "stemness", "differentiation")) %>% filter(score > 0)
  df_s_heat <- df_s_select %>% group_by(source, target) %>% summarise(score = mean(score), .groups = "drop")
  pdf(glue("{fig_dir}/CCC/xenium_ccc_heatmaps_cc_excludeEpi_{sample}.pdf"), width = 12, height = 8)
  p1 <- cc_heatmap(df_s_select, option = 'B', n_top_ints = 100) + 
    scale_fill_gradientn(colors = config_list$scale_7, oob = scales::squish)
  p2 <- cc_heatmap(df_s_select, option = 'A') + 
    scale_fill_gradientn(colors = config_list$scale_7, oob = scales::squish)
  p3 <- ggplot(df_s_heat, aes(target, source, fill = score)) + geom_tile(color = "white") +
    scale_fill_gradientn(colours = config_list$scale_7) +
    theme_classic()
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}





  feature_use <- c()
  for (assay_name in assayNames(mu)){
    assay_data <- assay(mu, assay_name) %>% as.data.frame()
    fea1 <- assay_data %>% dplyr::filter(differentiation > 0.1 | stemness > 0.1) %>% 
        rownames_to_column(var = "feature") %>% pull(feature)
    feature_use <- c(feature_use, fea1)
  }
  feature_use <- unique(feature_use)


  pdf(glue("{fig_dir}/CCC/xenium_ccc_heatmaps_{sample}.pdf"), width = 5, height = 8)
  for (. in assayNames(mu)) {
    tryCatch({
      p <- plotHeatmap(mu, exprs_values=.,
        features=feature_use,
        main=switch(., s="sender", r="receiver"),
        show_colnames=TRUE, center=TRUE, scale=TRUE)
      }, error = function(e) {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
      })
  }
  dev.off()









