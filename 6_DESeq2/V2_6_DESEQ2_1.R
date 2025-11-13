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
#Load Data:
obj<-readRDS(file = "/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/5_CellAssignment/V2_5_CellAssignment.rds")
Idents(object = obj)<- "final_assignments"
################################################################################
#Create Single Cell Experiment Object: 
counts<-obj@assays$RNA@counts                                                   #extract counts
metadata<-obj@meta.data                                                         #extract metadata
################################################################################
metadata$cluster_id<-factor(obj@active.ident)                                   #create cluster_id column

sce<-SingleCellExperiment(assays = list(counts = counts), 
                          colData = metadata)

saveRDS(sce, "V2_6_DESEQ2_1.rds")
################################################################################