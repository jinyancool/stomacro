# function to get polygon from boundary box
bbox_polygon <- function(x) {
  bb <- sf::st_bbox(x)
  
  p <- matrix(
    c(bb["xmin"], bb["ymin"], 
      bb["xmin"], bb["ymax"],
      bb["xmax"], bb["ymax"], 
      bb["xmax"], bb["ymin"], 
      bb["xmin"], bb["ymin"]),
    ncol = 2, byrow = T
  )
  
  sf::st_polygon(list(p))
}

# get circle edge for TMA core
bb_circle <- function(x){
  bb <- sf::st_bbox(x)
  width <- bb["xmax"] - bb["xmin"]
  height <- bb["ymax"] - bb["ymin"]
  radius <- min(width, height) / 2
  
  center <- sf::st_centroid(bbox_polygon(x))
  incircle <- sf::st_buffer(center, dist = radius)
  
  return(incircle)
}

# calculate window-based regional cell constitution
nb_freq <- function(i, df, x = "Centroid.X.µm", y = "Centroid.Y.µm", ct = "celltype",
                    distance = 50){
  # data
  query <- df[i, c(x, y)] %>% as.matrix()
  data <- df[, c(x, y)] %>% as.matrix()
  
  # neighborhood
  nb <- BiocNeighbors::queryNeighbors(data, query, threshold = distance)  
  idx <- nb$index[[1]]
  
  # celltype table
  tbl_ct <- df[idx, ] %>% pull(ct) %>% table
  df_ct <- as.data.frame(tbl_ct) %>%
    rename(celltype = ".") %>%
    tidyr::pivot_wider(names_from = celltype, values_from = Freq) %>%
    mutate(n = length(idx))
  return(df_ct)
}

# calculate pairwise cell-cell interaction
pci <- function(df_edge){
  # count edge number
  n_edge <- df_edge %>%
    mutate(cluster_min = pmin(cluster1, cluster2),
           cluster_max = pmax(cluster1, cluster2)) %>%
    select(-cluster1, -cluster2) %>%
    rename(cluster1 = cluster_min, cluster2 = cluster_max) %>%
    count(cluster1, cluster2, name = "Nij") %>%
    mutate(pair = paste(cluster1, cluster2, sep = "-"))
  
  # calculate log likelihood ratios and relative frequencies
  n_edge <- n_edge %>%
    group_by(cluster1) %>%
    mutate(Ni = sum(Nij)) %>%
    group_by(cluster2) %>%
    mutate(Nj = sum(Nij)) %>%
    ungroup() %>%
    mutate(Nt = sum(Nij)) %>%
    mutate(lr = Nij / Ni * Nt  / Nj,
           relative_freq = Nij / Ni,
           log2_lr = log2(lr)) %>%
    select(pair, cluster1, cluster2, lr, relative_freq, log2_lr)
  
  return(n_edge)
}



