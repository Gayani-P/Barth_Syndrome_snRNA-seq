################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0' 
REPO="http://cran.us.r-project.org"

.libPaths(c("",LIB)) 
.libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0",'/cluster/home/gperer01/R/x86_64-pc-linux-gnu-library/4.0', .libPaths()))

library(reticulate)
library(Seurat)
library(sctransform)
library(ggplot2)
library(plotly)
library(tidyr)
library(plyr)
library(dplyr)
library(cowplot)
library(BiocManager)
library(Matrix)
library(stringr)
library(igraph)
library(tidyverse)
library(httr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(Matrix)
library(devtools)
library(pryr)
library(clustree)


setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/4_ChooseR")
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4
#https://github.com/rbpatt2019/chooseR
################################################################################
#Load Helper Functions: 

#Find_Clusters: 
find_clusters <- function(
  obj,
  reduction = "pca",
  npcs = 16,                                            
  assay = "RNA",
  features = NULL,
  resolution = 0.8,
  verbose = FALSE) {
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = reduction,
    dims = 1:npcs,
    assay = integtrated,
    features = features,
    verbose = verbose,
    graph.name = paste(reduction, assay, sep = ".")
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = resolution,
    graph.name = paste(reduction, assay, sep = "."),
    verbose = verbose
  )
  return(obj)
}

#Multiple_Cluster
multiple_cluster <- function(
  obj,
  n = 100,
  size = 0.8,
  npcs = 16,
  res = 1.2,
  reduction = "pca",
  assay = "RNA") {
  
  # Initialise tibble for data
  clusters <- dplyr::as_tibble(Seurat::Cells(obj))
  clusters <- dplyr::rename(clusters, "cell" = value)
  
  # Get samples
  samples <- n_samples(n, Seurat::Cells(obj), size = size)
  
  # Repeated clusters
  j <- 1
  for (idx in samples) {
    message(paste0("\tClustering ", j, "..."))
    small_obj <- obj[, idx]
    small_obj <- find_clusters(
      small_obj,
      reduction = reduction,
      npcs = npcs,
      resolution = res,
      assay = assay
    )
    clusters <- dplyr::left_join(
      clusters,
      dplyr::as_tibble(Seurat::Idents(small_obj), rownames = "cell"),
      by = "cell"
    )
    j <- j + 1
  }
  return(clusters)
}

#N_Samples:
n_samples <- function(
  n,
  input,
  size = 0.8,
  replace = FALSE,
  simplify = FALSE) {
  splits <- replicate(
    n,
    sample(
      input,
      as.integer(length(input) * size),
      replace = replace
    ),
    simplify = simplify
  )
}

#Find_Matches:
find_matches <- function(col, df) {
  mtchs <- outer(df[[col]], df[[col]], "==")
  # Records drops as imaginary, mtchs as 1, not mtchs as 0
  mtchs[is.na(mtchs)] <- 1i
  return(mtchs)
}

#Percent_Match:
percent_match <- function(x, n = 100) {
  return(Re(x) / (n - Im(x)))
}

#Group_Scores:
group_scores <- function(tbl, clusters) {
  colnames(tbl) <- clusters
  data <- tbl %>%
    tibble::add_column("cell_1" = clusters) %>%
    tidyr::pivot_longer(-cell_1, names_to = "cell_2", values_to = "percent") %>%
    dplyr::group_by(cell_1, cell_2) %>%
    dplyr::summarise("avg_percent" = mean(percent)) %>%
    dplyr::ungroup()
  return(data)
}

#Group_Sil:
group_sil <- function(sil, res) {
  sil <- tibble::as_tibble(sil[, ]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise("avg_sil" = mean(sil_width)) %>%
    tibble::add_column("res" = res)
  return(sil)
}

#Boot_Median
boot_median <- function(x, interval = 0.95, R = 25000, type = "bca") {
  # Define median to take data and indices for use with boot::
  med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }
  
  # Calculate intervals
  boot_data <- boot::boot(data = x, statistic = med, R = R)
  boot_ci <- boot::boot.ci(boot_data, conf = interval, type = type)
  
  # Extract desired statistics
  ci <- list(
    low_med = boot_ci$bca[4],
    med = boot_ci$t0,
    high_med = boot_ci$bca[5]
  )
  return(ci)
}
################################################################################

################################################################################
#Load Data
obj<-readRDS(file = "/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/3_Harmony/out/V2_3_Harmony_18PCs.rds")
dim(obj)
obj<- obj[, sample(colnames(obj), size =40000, replace=F)]
################################################################################
npcs <- 18
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
assay <- "RNA"
reduction <- "harmony"
results_path <- "results/"
#Make sure to Make a Folder in Working Directory called "results"
################################################################################
#Run Pipeline:
for (res in resolutions) {
  #message(paste0("Clustering ", res, "..."))
  #message("\tFinding ground truth...")
  
  #"Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  obj <- find_clusters(obj,
                       reduction = reduction,
                       assay = assay,
                       resolution = res)
  
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]
  
  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
    obj,
    n = 100,
    size = 0.8,
    npcs = npcs,
    res = res,
    reduction = reduction,
    assay = assay
  )
  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  # This is the more time efficient vectorisation
  # However, it exhausts vector memory for (nearly) all datasets
  # matches <- purrr::map(columns, find_matches, df = results)
  # matches <- purrr::reduce(matches, `+`)
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }
  
  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )
  
  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "..."))
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
  
  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
}

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))
################################################################################
# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()
meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
    dplyr::filter(med >= threshold) %>%
    dplyr::arrange(n_clusters) %>%
    tail(n = 1) %>%
    dplyr::pull(res)
)

################################################################################
saveRDS(obj, "V2_4_ChooseR.rds")
################################################################################