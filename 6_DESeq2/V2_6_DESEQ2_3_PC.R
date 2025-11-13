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
library(DEGreport)
library(scran)
library(clusterProfiler)
library(GO.db)
library(org.Mm.eg.db) #Modify all instance os Org.Mm.eg.db to Mouse Reference Org.Mm.eg.db
library(xml2)
library(rvest)


setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/6_DESEQ2")
#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
#Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
#sequence count data: removing the noise and preserving large differences.
#Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
#https://github.com/hbc/NGS_Data_Analysis_Course/blob/master/sessionIII/lessons/04_DEG_deseq2.md
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1406-4
##Assigning Celltypes--GO Stats
#https://pubmed.ncbi.nlm.nih.gov/17098774/
#https://www.nature.com/articles/ng0500_25
################################################################################
#Load Data: 
obj<-readRDS(file = "/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/5_CellAssignment/V2_5_CellAssignment.rds")
Idents(object = obj)<- "final_assignments"

sce<-readRDS("V2_6_DESEQ2_1.rds")
metadata_ls<-readRDS("V2_6_DESEQ2_2_metadata.rds")
counts_ls<-readRDS("V2_6_DESEQ2_2_counts.rds")
cluster_names <- levels(colData(sce)$cluster_id)

#Select cell type of interest
cluster_names

#Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))
################################################################################
#CM/FB/EC/PC/SM/MAC/LEC/NEURO/TLYMPH/DC
################################################################################
idx <- which(names(counts_ls) == "PC") #Change for each cell type
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

#Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)

#Genotype Info
for (i in 1:dim(cluster_metadata)[1]) {
  if (cluster_metadata$sample_id[i] == "m5356"){cluster_metadata$genotype[i] <- "WT"}
  if (cluster_metadata$sample_id[i] == "m4909"){cluster_metadata$genotype[i] <- "WT"}
  if (cluster_metadata$sample_id[i] == "m4686"){cluster_metadata$genotype[i] <- "WT"}
  if (cluster_metadata$sample_id[i] == "m5355"){cluster_metadata$genotype[i] <- "WT"}
  
  
  if (cluster_metadata$sample_id[i] == "m5301"){cluster_metadata$genotype[i] <- "KO"}
  if (cluster_metadata$sample_id[i] == "m4685"){cluster_metadata$genotype[i] <- "KO"}
  if (cluster_metadata$sample_id[i] == "m5300"){cluster_metadata$genotype[i] <- "KO"}
  if (cluster_metadata$sample_id[i] == "m4910"){cluster_metadata$genotype[i] <- "KO"}
  
}

#Batch Info
for (i in 1:dim(cluster_metadata)[1]) {
  if (cluster_metadata$sample_id[i] == "m5356"){cluster_metadata$batch[i] <- "A"}
  if (cluster_metadata$sample_id[i] == "m5301"){cluster_metadata$batch[i] <- "A"}
  if (cluster_metadata$sample_id[i] == "m4909"){cluster_metadata$batch[i] <- "B"}
  if (cluster_metadata$sample_id[i] == "m4686"){cluster_metadata$batch[i] <- "B"}
  if (cluster_metadata$sample_id[i] == "m5355"){cluster_metadata$batch[i] <- "B"}
  if (cluster_metadata$sample_id[i] == "m4685"){cluster_metadata$batch[i] <- "B"}
  if (cluster_metadata$sample_id[i] == "m5300"){cluster_metadata$batch[i] <- "C"}
  if (cluster_metadata$sample_id[i] == "m4910"){cluster_metadata$batch[i] <- "C"}
  
}

#Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ genotype)
sizeFactors(dds)
dds <- estimateSizeFactors(dds, type = "poscounts")
sizeFactors(dds)
################################################################################
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "genotype")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "batch") 
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
################################################################################
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ genotype + batch)
sizeFactors(dds)
dds <- estimateSizeFactors(dds, type = "poscounts")
sizeFactors(dds)

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ batch,
                 minmu = 1e-6, 
                 minReplicatesForReplace = Inf)
plotDispEsts(dds_lrt) #Check that dispersion decreases with increasing mean

#Extract results
res_LRT <- results(dds_lrt)

#Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>% 
  as_tibble()

#DON'T FILTER LRT RESULTS BY LFC!!!
padj.cutoff<-0.001
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

write.csv(sig_res_LRT,
          paste0("results/", unique(cluster_metadata$cluster_id), "_", 
                 levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1], "LRT_signif_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
################################################################################
#Plotting DEGs
################################################################################
#Subset results for faster cluster finding (for classroom demo purposes)
clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=1000)


#Obtain rlog values for those significant genes
cluster_rlog <- rld[clustering_sig_genes$gene, ]
cluster_rlog <- assay(cluster_rlog)

cluster_metadata$genotype<-as.factor(cluster_metadata$genotype)


#Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = cluster_metadata, time = "genotype", col = NULL, minc = 0)
ggsave("images/PC_LRT_direction.of.expression.png", units = "px", width = 4500, height =2970)


#Extract the Group genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%filter(cluster == 1)
group2<-clusters$df %>% filter(cluster == 2)
################################################################################
#GENE ONTOLOGY ANALYSIS OF GROUPS:
################################################################################
obj_PC<-subset(x = obj, ident = "PC", downsample = 20000) #Modify here
#Gene Universe
GeneUniverse<-rownames(cluster_counts)
GeneUniverseIDs<-select(org.Mm.eg.db, GeneUniverse, "ENTREZID", "SYMBOL")

#GROUP 1: GO
#DEG list
list<-as.character(group1$gene)  
list.geneIDs<-select(org.Mm.eg.db, list, "ENTREZID", "SYMBOL")%>%filter(!is.na(ENTREZID))

group1.ego<-enrichGO(gene = list.geneIDs$ENTREZID, 
                     universe = GeneUniverseIDs$ENTREZID, 
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP",             # Do "MF" for Molecular Function Gene Ontology                                       
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
head(group1.ego)

result<-group1.ego@result%>%filter(p.adjust < 0.05)
write.csv(result, "results/GO_Terms/PC_group1_BP_sigGOTerms.csv")

#GROUP 2: GO
#DEG list
list2<-as.character(group2$gene)  
list2.geneIDs<-select(org.Mm.eg.db, list2, "ENTREZID", "SYMBOL")%>%filter(!is.na(ENTREZID))

group2.ego<-enrichGO(gene = list2.geneIDs$ENTREZID, 
                     universe = GeneUniverseIDs$ENTREZID, 
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP",             # Do "MF" for Molecular Function Gene Ontology                                       
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
head(group2.ego)

result2<-group2.ego@result%>%filter(p.adjust < 0.05)
write.csv(result, "results/GO_Terms/PC_group2_BP_sigGOTerms.csv")
################################################################################
#Picking features to generate heatmap of
#Features picked by those contributing to GOterms

#Group1 Features:
features.g1 = list()
for (i in 1:dim(result)[1]) {
  list = strsplit(result$geneID[i], "/", fixed = TRUE)
  features.g1 = c(list[[1]], features.g1)
}
features.g1<-unique(features.g1)%>%unlist()

#Group2 Features:
features.g2 = list()
for (i in 1:dim(result2)[1]) {
  list = strsplit(result2$geneID[i], "/", fixed = TRUE)
  features.g2 = c(list[[1]], features.g2)
}
features.g2<-unique(features.g2)%>%unlist()

features<-c(features.g1, features.g2)

#Heatmap
DoHeatmap(
  obj_PC,                                                                   #Modify Here to features
  features = features,
  cells = NULL,
  group.by = "genotype",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE) + ggtitle("PC Group.1 Heatmap (20,000 cells)")
ggsave("images/PC_group.1_heatmap.png", units = "px", width = 4500, height =2970)
###############################################################################

