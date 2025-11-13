#Load Libraries
################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0'
LIBHOME='/cluster/home/gperer01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIBHOME,LIB,""))

library(ggplot2)
library(Seurat)

################################################################################
#Load Data
setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/1_DoubletFinder/out")

m4685<-readRDS(file = "m4685_Doublets.Removed.rds")
m4686<-readRDS(file = "m4686_Doublets.Removed.rds")
m5300<-readRDS(file = "m5300_Doublets.Removed.rds")
m5301<-readRDS(file = "m5301_Doublets.Removed.rds")
m5355<-readRDS(file = "m5355_Doublets.Removed.rds")
m5356<-readRDS(file = "m5356_Doublets.Removed.rds")
m4909<-readRDS(file = "m4909_Doublets.Removed.rds")
m4910<-readRDS(file = "m4910_Doublets.Removed.rds")


setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/2_QualityControl")

#Merge Seurat Objects: 
obj<-merge(m4685, y = c(m4910, m5300, m5301, 
                           m4686, m4909, m5355, m5356))
rm(m4685, m4910, m5300, m5301, 
   m4686, m4909, m5355, m5356)

#Trim Metadata/Add Complexity Metric
obj@meta.data$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
metadata<-obj@meta.data
metadata<-metadata%>%dplyr::select(orig.ident, 
                                   nCount_RNA, 
                                   nFeature_RNA, 
                                   percent.MT, 
                                   cells, 
                                   batch, 
                                   human,
                                   genotype, 
                                   barcode, 
                                   log10GenesPerUMI)
################################################################################
#Making some corrections to Metadata:
names(metadata)[names(metadata) == 'human'] <- 'mouse'
obj@meta.data<-metadata

metadata<-obj@meta.data

for (i in 1:dim(metadata)[1]) {
  if (metadata$mouse[i] == "m4685"){metadata$genotype[i] = "KO"} else if
  (metadata$mouse[i] == "m4910"){metadata$genotype[i] = "KO"} else if
  (metadata$mouse[i] == "m5300"){metadata$genotype[i] = "KO"} else if
  (metadata$mouse[i] == "m5301"){metadata$genotype[i] = "KO"} else if
  (metadata$mouse[i] == "m4686"){metadata$genotype[i] = "WT"} else if
  (metadata$mouse[i] == "m4909"){metadata$genotype[i] = "WT"} else if
  (metadata$mouse[i] == "m5355"){metadata$genotype[i] = "WT"} else if
  (metadata$mouse[i] == "m5356"){metadata$genotype[i] = "WT"}}

obj@meta.data<-metadata
################################################################################
#PRE-QC Plots: 
################################################################################
#Plot NCells per Genotype
metadata %>% ggplot(aes(x=genotype, fill=genotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("TAZ Mice: Pre-QC NCells per Genotype")            
ggsave(filename = "images/PreQC.NCells.per.Genotype.png")

#Plot NUMIs per Cell            
metadata %>% ggplot(aes(color=genotype, x=nCount_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800) +
  ggtitle("TAZ Mice: Pre-QC NUMIs per Cell")
ggsave(filename = "images/PreQC.NUMIs.per.Cell.png")         


#Plot NGenes per Cell
metadata %>% ggplot(aes(color=genotype, x=nFeature_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)+
  ggtitle("TAZ Mice: Pre-QC NGenes per Cell")
ggsave(filename = "images/PreQC.NGenes.per.Cell.png")

#Plot NUMIs vs. NGenes
metadata %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 800) +
  geom_hline(yintercept = 300) +
  ggtitle("TAZ Mice: Pre-QC NUMIs vs NGenes")
ggsave(filename = "images/PreQC.NUMIs.vs.NGenes.png")

#Plot Complexity of Gene Expression(Genes per UMI)
metadata %>% ggplot(aes(x=log10GenesPerUMI, color = genotype, fill=genotype)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("TAZ Mice: Pre-QC NGenes per UMI (Library Complexity)")
ggsave(filename = "images/PreQC.NGenes.per.UMI.(Library.Complexity).png")
################################################################################
#QC Filters: 
################################################################################
#Cell Level Filtering
filtered_human_matrix <- subset(x = obj, subset = (nCount_RNA >= 500) & 
                                  (nFeature_RNA >= 250) & 
                                  (log10GenesPerUMI > 0.80))

#Gene-Level Filtering:
dim(obj)
dim(filtered_human_matrix)
#Extract counts
counts <- GetAssayData(object = filtered_human_matrix, slot = "counts")
dim(counts)

#Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
dim(nonzero)

#Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
length(keep_genes)

#Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
dim(filtered_counts)

#Extract Metadata
metadata<-filtered_human_matrix@meta.data
rm(filtered_human_matrix)

#Reassign to filtered Seurat object
filtered_human_matrix <- CreateSeuratObject(filtered_counts, meta.data = metadata)
dim(filtered_human_matrix)
rm(filtered_counts, metadata, nonzero,counts, keep_genes, obj)
################################################################################
#PostQC Plots:
################################################################################
metadata<-filtered_human_matrix@meta.data 

#Plot NCells per Genotype
metadata %>% ggplot(aes(x=genotype, fill=genotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("TAZ Mice: Post-QC NCells per Genotype")            
ggsave(filename = "images/PostQC.NCells.per.Genotype.png")

#Plot NUMIs per Cell            
metadata %>% ggplot(aes(color=genotype, x=nCount_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800) +
  ggtitle("TAZ Mice: Post-QC NUMIs per Cell")
ggsave(filename = "images/PostQC.NUMIs.per.Cell.png")         


#Plot NGenes per Cell
metadata %>% ggplot(aes(color=genotype, x=nFeature_RNA, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)+
  ggtitle("TAZ Mice: Post-QC NGenes per Cell")
ggsave(filename = "images/PostQC.NGenes.per.Cell.png")

#Plot NUMIs vs. NGenes
metadata %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 800) +
  geom_hline(yintercept = 300) +
  ggtitle("TAZ Mice: Post-QC NUMIs vs NGenes")
ggsave(filename = "images/PostQC.NUMIs.vs.NGenes.png")

#Plot Complexity of Gene Expression(Genes per UMI)
metadata %>% ggplot(aes(x=log10GenesPerUMI, color = genotype, fill=genotype)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("TAZ Mice: Post-QC NGenes per UMI (Library Complexity)")
ggsave(filename = "images/PostQC.NGenes.per.UMI.(Library.Complexity).png")
################################################################################
#Decide Which Variables to Regress in Addition to nUMI: 
################################################################################
cell_cycle_genes <- read.csv(file = "/cluster/tufts/chinlab/single_nuclei_cryab_mice/1_preprocessing/2_Harmony/cell_cycle_genes.csv")
g2m.features<-cell_cycle_genes$phase == "G2/M"
g2m.features<-cell_cycle_genes[g2m.features, ]
g2m.features<-g2m.features$name

s.features<-cell_cycle_genes$phase == "S"
s.features<-cell_cycle_genes[s.features, ]
s.features<-s.features$name

#Normalize
seurat_phase<-NormalizeData(filtered_human_matrix, verbose = TRUE)

#Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.features, 
                                 s.features = s.features)

#Identify top 2000 variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

#Scale the counts
seurat_phase <- ScaleData(seurat_phase)

#Perform PCA
seurat_phase <- RunPCA(seurat_phase)

#Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
ggsave(filename = "images/Cell_Cycle_in_PC1.png")
################################################################################
#Save Seurat Object:
################################################################################
saveRDS(filtered_human_matrix, "V2_2_QualityControl_1.rds")
################################################################################