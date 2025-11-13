#Load Libraries
################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0'
LIBHOME='/cluster/home/gperer01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIBHOME,LIB,""))

library(DoubletFinder)
library(cowplot)

setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/1_DoubletFinder")
#https://github.com/chris-mcginnis-ucsf/DoubletFinder
#https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
################################################################################
#Load Data
m4685.matrix <- readRDS("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/0_SoupX/out/m4685_SoupXObject_Corrected_Matrix.rds")
dim(m4685.matrix)

#Add Genotype and Mouse Labels to Cell Barcodes
colnames(m4685.matrix)<-paste("HCM-m4685", colnames(m4685.matrix), sep = "-")

#Create Seurat Object
m4685<-CreateSeuratObject(counts = m4685.matrix, 
                          min.cells = 3,          #include features detected in at least 3 cells
                          min.features = 200)     #include cells with at least 200 features(done now to reduce size of dataset--but I'm not using this as final empty droplet threshold) 
rm(m4685.matrix)
################################################################################
#Metadata: 

#Add mito.reads to metadata(Capital "MT" for human/cat data; lowercase "mt" for mice)
m4685[["percent.MT"]] <- PercentageFeatureSet(m4685, pattern = "^mt-")

#Remove Mitochondrial Reads:
dim(m4685)
m4685.rownames <- rownames(m4685)
m4685_mito.rows <- grep(pattern = "^mt-", x = m4685.rownames) 
m4685 <- m4685[-m4685_mito.rows,]
dim(m4685)

#Extract Metadata from Seurat Object
metadata<-m4685@meta.data    

#Make Column for Cell Barcodes(Labels)
metadata$cells<-rownames(metadata)

#Add Batch to Metadata
metadata$batch <-NA
metadata$batch[which(str_detect(metadata$cells, "m4685"))]<-"B"

metadata_to_add<-str_split_fixed(metadata$cells, pattern = "-", n = 3)
colnames(metadata_to_add)<-c("genotype", "human", "barcode")

metadata<-cbind(metadata, metadata_to_add)
rm(metadata_to_add) 

#Add Meta Data Frame back into Seurat Object 
m4685@meta.data<-metadata  
################################################################################
#Pre-QC Plots:

VlnPlot(m4685, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
ggsave(filename = "images/m4685_PreQC_nFeature_nCount_percent.MT.png")

#Plot NUMIs per Cell            
metadata %>% ggplot(aes(color=genotype, x=nCount_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800) +
  ggtitle("m4685: Pre-QC NUMIs per Cell")
ggsave(filename = "images/m4685_PreQC.NUMIs.per.Cell.png")         

#Plot NGenes per Cell
metadata %>% ggplot(aes(color=genotype, x=nFeature_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)+
  ggtitle("m4685: Pre-QC NGenes per Cell")
ggsave(filename = "images/m4685_PreQC.NGenes.per.Cell.png")

#Plot percent.mt per Cell            
metadata %>% ggplot(aes(color=genotype, x=percent.MT, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 20) +
  ggtitle("m4685: Pre-QC Percent Mito per Cell")
ggsave(filename = "images/m4685_PreQC.PercentMito.per.Cell.png")    

################################################################################      
#Quality Control Filters:

#Cell-Level Filtering
filtered_m4685 <- subset(x = m4685, subset = (nCount_RNA >= 500) & 
                           (nFeature_RNA >= 250))
dim(filtered_m4685)

#Gene-Level Filtering: 
#Extract counts
counts <- GetAssayData(object = filtered_m4685, slot = "counts")

#Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

#Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

#Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

#Reassign to filtered Seurat object
filtered_m4685 <- CreateSeuratObject(filtered_counts, meta.data = filtered_m4685@meta.data)
rm(filtered_counts, metadata, nonzero,counts,keep_genes)

################################################################################ 
#Standard Seurat Workflow: 

m4685<-NormalizeData(filtered_m4685, normalization.method = "LogNormalize", scale.factor = 10000)
rm(filtered_m4685)
m4685<-FindVariableFeatures(m4685, selection.method = "vst", nfeatures = 2000)
m4685<-ScaleData(m4685, vars.to.regress = "percent.MT")
m4685<-RunPCA(m4685)

#####################################
#Choose #PCs for Dimension Reduction: 
ElbowPlot(object = m4685, ndims = 40)

#Determine percent of variation associated with each PC
pct <- m4685[["pca"]]@stdev / sum(m4685[["pca"]]@stdev) * 100

#Calculate cumulative percents for each PC
cumu <- cumsum(pct)

#Determine which PC exhibits cumulative percent greater than 90% and %variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

#Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

#Last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

#Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

#Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave("images/m4685_ElbowPlot.png")

rm(plot_df, co1, co2, cumu, pct)
#####################################
m4685<-RunUMAP(m4685, dims = 1:pcs)
m4685<-FindNeighbors(m4685, reduction = "pca")
m4685<-FindClusters(m4685, resolution = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.3, 2.0))
################################################################################
#pK Identification(no ground-truth)
sweep.res.list.m4685 <- paramSweep_v3(m4685, PCs = 1:pcs, sct = FALSE)
sweep.stats.m4685 <- summarizeSweep(sweep.res.list.m4685, GT = FALSE)
bcmvn.m4685 <- find.pK(sweep.stats.m4685)

bcmvn.m4685 #Choose appropriate pK value from this table
pK=as.numeric(as.character(bcmvn.m4685$pK))
BCmetric=bcmvn.m4685$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b", col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
ggsave("images/m4685_Choose.pK.Value.png")

##Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(m4685@meta.data$seurat_clusters)
nExp_poi <- round(0.08*nrow(m4685@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
m4685 <- doubletFinder_v3(m4685, PCs = 1:pcs, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN.name = colnames(m4685@meta.data)[grepl("pANN", colnames(m4685@meta.data))]

m4685 <- doubletFinder_v3(m4685, PCs = 1:pcs, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = pANN.name, sct = FALSE)
DF.name = colnames(m4685@meta.data)[grepl("DF.classification", colnames(m4685@meta.data))]

#Plot Doublets(no Homotypic Doublet Estimate adjustment?)
UMAP_atlas<-DimPlot(m4685, reduction = "umap", label = TRUE, pt.size=.1, shuffle = FALSE, raster = FALSE) + ggtitle("UMAP m4685")
UMAP_by_doublet<-DimPlot(m4685, reduction = "umap", group.by = DF.name[2], pt.size = .1) + ggtitle("UMAP m4685 Doublets") 
plot_grid(UMAP_atlas, UMAP_by_doublet)
ggsave("images/m4685_UMAP.with.Doublets.png")

#Plot UMIs for doublet vs. singlet
VlnPlot(m4685, features = "nFeature_RNA", group.by = DF.name[2], pt.size = 0.1)
ggsave("images/m4685_nCounts.Singlets.vs.Doublets.png")

#Remove Doublets
m4685 = m4685[, m4685@meta.data[, DF.name[2]] == "Singlet"]
UMAP_atlas<-DimPlot(m4685, reduction = "umap", label = TRUE, pt.size=.1, shuffle = FALSE, raster = FALSE)+ggtitle("UMAP m4685 Doublets Removed")
UMAP_atlas
ggsave("images/m4685_UMAP.Doublets.Removed.png")

dim(m4685)
saveRDS(m4685, file = "out/m4685_Doublets.Removed.rds")
