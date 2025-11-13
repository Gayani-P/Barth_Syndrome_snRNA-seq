#Load Libraries
################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0'
LIBHOME='/cluster/home/gperer01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIBHOME,LIB,""))

library(ggplot2)
library(Seurat)
library(harmony)
library(cowplot)

setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/3_Harmony")
################################################################################
#Load Data: 
################################################################################
obj<-readRDS("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/2_QualityControl/V2_2_QualityControl_1.rds")
################################################################################
#Standard Seurat Workflow:
################################################################################
obj<-NormalizeData(obj)
obj<-FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj<-ScaleData(obj, vars.to.regress = "nCount_RNA")
obj<-RunPCA(obj)
 ################################################################################
#PRE-HARMONY Plots:
################################################################################
p1<-DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "genotype")
p2<-VlnPlot(object = obj, features = "PC_1", group.by = "genotype", pt.size = .1)
plot_grid(p1,p2)
ggsave("images/Pre-Harmony_PC1.png")
################################################################################
#Run Harmony:
################################################################################
obj<-RunHarmony(obj, "genotype", plot_convergence = TRUE, epsilon.cluster = -Inf, epsilon.harmony = -Inf)
ggsave("images/Harmony_Convergence_Plot.png")
################################################################################
#POST-HARMONY Plots:
################################################################################
p1 <- DimPlot(object = obj, reduction = "harmony", pt.size = .1, group.by = "genotype")
p2 <- VlnPlot(object = obj, features = "harmony_1", group.by = "genotype", pt.size = .1)
plot_grid(p1,p2)
ggsave("images/Post-Harmony.png")
################################################################################
#Choose #PCs for Dimension Reduction: 
################################################################################
#Determine percent of variation associated with each PC
pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100

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
ggsave("images/Elbow_Plot.png")

rm(plot_df, co1, co2, cumu, pct)
pcs
################################################################################
#Non-Linear Dimension Reduction & Clustering: 
################################################################################
obj<-RunUMAP(obj, dims = 1:pcs, reduction = "harmony") 
obj<-FindNeighbors(obj, dims = 1:pcs, reduction = "harmony")
obj<-FindClusters(obj, resolution = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 5))
################################################################################
#Save Seurat Object:
################################################################################
saveRDS(obj, file = "out/V2_3_Harmony_18PCs.rds")
################################################################################