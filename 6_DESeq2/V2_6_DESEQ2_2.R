################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0' 
REPO="http://cran.us.r-project.org"

.libPaths(c("",LIB)) 
.libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0",
            '/cluster/home/gperer01/R/x86_64-pc-linux-gnu-library/4.0', .libPaths()))

library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(Seurat)


setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/6_DESEQ2")
#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
################################################################################
#Read Single Cell Object: 
sce<-readRDS("V2_6_DESEQ2_1.rds")

#Check the assays present
assays(sce)

#Check the Counts Matrix
dim(counts(sce))
counts(sce[1:6, 1:6])

#Explore Cell Metadata
dim(colData(sce))
head(colData(sce))
################################################################################
#Determine the number of clusters and cluster names present in data: 
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names
length(cluster_names) #Total number of clusters

#Determine the number of Samples present in data: 
sce$mouse<-as.factor(sce$mouse)
sample_names<-levels(sce$mouse) #sample_id = mouse
sample_names

#Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "mouse")]
head(groups)

#Aggregate across cluster-sample groups
#transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

#Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

#Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

# Understanding tstrsplit()
## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
################################################################################
# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names


# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

saveRDS(counts_ls, "V2_6_DESEQ2_2_counts.rds")
################################################################################
# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(genotype, batch, mouse)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

# Rename rows
rownames(metadata) <- metadata$mouse
head(metadata)

# Number of cells per sample and cluster
t <- table(colData(sce)$mouse,
           colData(sce)$cluster_id)
t[1:6, 1:6]


# Creating metadata list

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

saveRDS(metadata_ls, "V2_6_DESEQ2_2_metadata.rds")
################################################################################