################################################################################
#Make sure gsl/2.6 module is loaded 

LIB='/cluster/tufts/hpc/tools/R/4.0.0' 
REPO="http://cran.us.r-project.org"
.libPaths(c("",LIB)) 

library(scLR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(igraph)
library(CrossTalkeR)
library(factoextra)
library(ggrepel)
library(stringr)
require(devtools)
library(OmnipathR)
library(ggraph)
library(magrittr)


setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/7_scLR")
#https://www.bioconductor.org/packages/release/bioc/vignettes/OmnipathR/inst/doc/paths.html#protein-protein-interactions
################################################################################
#Load Data
obj<-readRDS(file = "/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/5_CellAssignment/V2_5_CellAssignment.rds")
#obj<-subset(obj, downsample = 1000)

#Set Idents
Idents(object = obj)<- "final_assignments"

#Determine the number of cells per cluster per genotype
n_cells<-table(obj@active.ident, obj@meta.data$mouse)
n_cells

#Add cell type assignments to metadata:
obj@meta.data$assignments<-obj@active.ident

#Extract count matrix
countmatrix<-obj@assays$RNA@counts
countmatrix[1:6, 1:10]

#Extract Cell Info
cellinfo<-obj@meta.data
cellinfo<-cellinfo%>%dplyr::select(mouse, genotype, assignments)
head(cellinfo, 10)
################################################################################
#OmnipathR
ppi<-import_intercell_network(transmitter_param = list(categories = c("ligand")),receiver_param = list(categories = c("receptor")))
ppi<-ppi_new%>%filter(entity_type_intercell_source == "protein" & entity_type_intercell_target == "protein")%>%dplyr::select(source_genesymbol, target_genesymbol)

lrpairs0<-ppi%>%dplyr::select(source_genesymbol, target_genesymbol)
output<-scLR(countmatrix, 
               cellinfo, 
               lrpairs.sample = lrpairs0,
               normalization = "ByCT",
               bntest = TRUE,
               sig.bntest = 0.05,
               rhos = "est",
               parallel.use = TRUE, 
               impute.miss.celltype = 0) 

saveRDS(output, "V2_7_scLR_1.rds")
write.csv(output$Rs, "V2_7_scLR_0.csv")
################################################################################
res<-output$Rs
res.sig<-res%>%dplyr::filter(adj.p < 0.05)
res.sig<-res.sig%>%dplyr::filter(bntest.adj.p > 0.05)
write.csv(res.sig, "V2_7_scLR_1_sig.csv")
################################################################################
