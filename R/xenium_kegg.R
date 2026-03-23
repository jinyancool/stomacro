suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(glue)
  library(writexl)
  library(openxlsx)
})

set.seed(12387)
rds_dir = "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/rds/xenium"
fig_dir <- "/cluster/home/lixiyue_jh/projects/stomatology/analysis/lvjiong/human/meta/manuscript/figs/fig4"

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





erosion_region_diff_gene <- readRDS(glue("{rds_dir}/erosion_region_diff_gene.rds"))
dt_mpg_region <- erosion_region_diff_gene$Macrophage_diff_gene_sig
dt_epi_region <- erosion_region_diff_gene$Epithelial_diff_gene_sig
erosion_sample_diff_gene <- readRDS(glue("{rds_dir}/erosion_sample_diff_gene.rds"))
dt_mpg_sample <- erosion_sample_diff_gene$Macrophage_diff_gene_sig
dt_epi_sample <- erosion_sample_diff_gene$Epithelial_diff_gene_sig


map_dt <- list("region_mpg" = dt_mpg_region,
               "region_epi" = dt_epi_region,
               "sample_mpg" = dt_mpg_sample,
               "sample_epi" = dt_epi_sample)

res_total_list <- list()
res_read_list <- list()
for (dt_name in names(map_dt)){
  dt <- map_dt[[dt_name]]
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
  
  name_up_plot <- paste0(dt_name, "_up_plot")
  name_dn_plot <- paste0(dt_name, "_down_plot")
  name_up_read <- paste0(dt_name, "_up_result")
  name_dn_read <- paste0(dt_name, "_down_result")
  
  res_total_list[[name_up_plot]] <- res_up
  res_total_list[[name_dn_plot]] <- res_dn
  res_total_list[[name_up_read]] <- res_up_readable
  res_total_list[[name_dn_read]] <- res_dn_readable
  
  res_read_list[[name_up_read]] <- res_up_readable
  res_read_list[[name_dn_read]] <- res_dn_readable
  
}

saveRDS(res_total_list, glue("{rds_dir}/kegg_result.rds"))

save_list_to_excel(res_read_list, glue("{rds_dir}/kegg_erosion.xlsx"), font_size = 12, font_name = "Arial")




epi_cluster_diff_gene <- readRDS(glue("{rds_dir}/epithelial_cluster_marker_gene.rds"))
select_cluster <- c("basal_CD36+", "basal_EPCAM+", "basal_KRAS+", "basal_invasive", "basal_cell", "CCND2+SFRP1")

dt_diff <- epi_cluster_diff_gene$filter %>% as.data.frame() %>% filter(cluster %in% select_cluster) 
map_res <- compareCluster(ENTREZID ~ cluster,
                          data = dt_diff,
                          organism = "hsa",
                          fun = "enrichKEGG",
                          use_internal_data = TRUE)
mpg_res_readable <- setReadable(map_res, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
saveRDS(map_res, glue("{rds_dir}/kegg_epicluster_selectcluster_result.rds"))
save_list_to_excel(mpg_res_readable, glue("{rds_dir}/kegg_epicluster_selectcluster_result_readable.xlsx"), font_size = 12, font_name = "Arial")

p <- dotplot(map_res, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(glue("{fig_dir}/xenium_kegg_epi_cluster.pdf"), p, width = 7, height = 7)
