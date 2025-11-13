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
output<-readRDS("V2_7_scLR_1.rds")
res<-output$Rs
res.sig<-res%>%dplyr::filter(adj.p < 0.05)
res.sig<-res.sig%>%dplyr::filter(bntest.adj.p > 0.05)
################################################################################

################################################################################
#obs.xy.diff: WT- KO
lr.pair<-as.data.frame(str_split_fixed(res.sig$lr.gene.name, pattern = "-", n = 2))
colnames(lr.pair)<-c("Ligand", "Receptor")
cell.pair<-as.data.frame(str_split_fixed(res.sig$lr.cell.name, pattern = "-", n = 2))
colnames(cell.pair)<-c("Ligand.Cluster", "Receptor.Cluster")
res.sig.table<-as.data.frame(cbind(lr.pair$Ligand, 
                                   cell.pair$Ligand.Cluster, 
                                   lr.pair$Receptor, 
                                   cell.pair$Receptor.Cluster, 
                                   as.numeric(res.sig$obs.xy.diff)))
colnames(res.sig.table)<-c("Ligand", "Ligand_Cluster", "Receptor", "Receptor_Cluster", "MeanLR")
rm(lr.pair, cell.pair)

#Filter by which condition the LR pairs are increased in(WT-WT): 
WT<-dplyr::filter(res.sig, obs.xy.diff > 0)
KO<-dplyr::filter(res.sig, obs.xy.diff < 0)
################################################################################
#Store all potential cell-cell interactions increased in WT
################################################################################
celltype.pairs<-c("CM-CM", 
                  "CM-FB",
                  "CM-EC",
                  "CM-PC",
                  "CM-SM",
                  "CM-MAC",
                  "CM-MES",
                  "CM-LEC", 
                  "CM-BLYMPH",
                  "CM-TLYMPH",
                  "CM-ADIPO",
                  "CM-DC",
                  
                  "FB-CM", 
                  "FB-FB",
                  "FB-EC",
                  "FB-PC",
                  "FB-SM",
                  "FB-MAC",
                  "FB-MES",
                  "FB-LEC", 
                  "FB-BLYMPH",
                  "FB-TLYMPH",
                  "FB-ADIPO",
                  "FB-DC",
                  
                  "EC-CM", 
                  "EC-FB",
                  "EC-EC",
                  "EC-PC",
                  "EC-SM",
                  "EC-MAC",
                  "EC-MES",
                  "EC-LEC", 
                  "EC-BLYMPH",
                  "EC-TLYMPH",
                  "EC-ADIPO",
                  "EC-DC",
                  
                  "PC-CM", 
                  "PC-FB",
                  "PC-EC",
                  "PC-PC",
                  "PC-SM",
                  "PC-MAC",
                  "PC-MES",
                  "PC-LEC", 
                  "PC-BLYMPH",
                  "PC-TLYMPH",
                  "PC-ADIPO",
                  "PC-DC",
                  
                  "SM-CM", 
                  "SM-FB",
                  "SM-EC",
                  "SM-PC",
                  "SM-SM",
                  "SM-MAC",
                  "SM-MES",
                  "SM-LEC", 
                  "SM-BLYMPH",
                  "SM-TLYMPH",
                  "SM-ADIPO",
                  "SM-DC",
                  
                  "MAC-CM", 
                  "MAC-FB",
                  "MAC-EC",
                  "MAC-PC",
                  "MAC-SM",
                  "MAC-MAC",
                  "MAC-MES",
                  "MAC-LEC", 
                  "MAC-BLYMPH",
                  "MAC-TLYMPH",
                  "MAC-ADIPO",
                  "MAC-DC",
                  
                  "MES-CM", 
                  "MES-FB",
                  "MES-EC",
                  "MES-PC",
                  "MES-SM",
                  "MES-MAC",
                  "MES-MES",
                  "MES-LEC", 
                  "MES-BLYMPH",
                  "MES-TLYMPH",
                  "MES-ADIPO",
                  "MES-DC",
                  
                  "LEC-CM", 
                  "LEC-FB",
                  "LEC-EC",
                  "LEC-PC",
                  "LEC-SM",
                  "LEC-MAC",
                  "LEC-MES",
                  "LEC-LEC", 
                  "LEC-BLYMPH",
                  "LEC-TLYMPH",
                  "LEC-ADIPO",
                  "LEC-DC",
                  
                  "BLYMPH-CM", 
                  "BLYMPH-FB",
                  "BLYMPH-EC",
                  "BLYMPH-PC",
                  "BLYMPH-SM",
                  "BLYMPH-MAC",
                  "BLYMPH-MES",
                  "BLYMPH-LEC", 
                  "BLYMPH-BLYMPH",
                  "BLYMPH-TLYMPH",
                  "BLYMPH-ADIPO",
                  "BLYMPH-DC",
                  
                  "TLYMPH-CM", 
                  "TLYMPH-FB",
                  "TLYMPH-EC",
                  "TLYMPH-PC",
                  "TLYMPH-SM",
                  "TLYMPH-MAC",
                  "TLYMPH-MES",
                  "TLYMPH-LEC", 
                  "TLYMPH-BLYMPH",
                  "TLYMPH-TLYMPH",
                  "TLYMPH-ADIPO",
                  "TLYMPH-DC",
                  
                  "ADIPO-CM", 
                  "ADIPO-FB",
                  "ADIPO-EC",
                  "ADIPO-PC",
                  "ADIPO-SM",
                  "ADIPO-MAC",
                  "ADIPO-MES",
                  "ADIPO-LEC", 
                  "ADIPO-BLYMPH",
                  "ADIPO-TLYMPH",
                  "ADIPO-ADIPO",
                  "ADIPO-DC",
                  
                  "DC-CM", 
                  "DC-FB",
                  "DC-EC",
                  "DC-PC",
                  "DC-SM",
                  "DC-MAC",
                  "DC-MES",
                  "DC-LEC", 
                  "DC-BLYMPH",
                  "DC-TLYMPH",
                  "DC-ADIPO",
                  "DC-DC")

#create counts column
celltype.pairs<-as.data.frame(celltype.pairs)
colnames(celltype.pairs)<-c("pairs")
celltype.pairs$counts<-0

#count instances of cell-cell interactions
count<-WT%>%dplyr::count(lr.cell.name)
colnames(count)<-c("pairs", "counts")
#Add counts to dataframe with all potential cell-cell interactions: 
count<-left_join(celltype.pairs, count, by = "pairs" )%>%
  dplyr::select(pairs, counts.y)
count$counts.y[is.na(count$counts.y)] = 0.001 #Make 0's 0.001 for graphing purposes downstream
colnames(count)<-c("pairs", "counts")
count$counts<-count$counts * 0.1 #scaling edges

#Ligand/Receptor Counts
LR_celltypes<-as.data.frame(str_split_fixed(celltype.pairs$pairs, pattern = "-", n = 2))
colnames(LR_celltypes)<-c("Ligand", "Receptor")

ligand.count<-LR_celltypes%>%dplyr::count(Ligand)
receptor.count<-LR_celltypes%>%dplyr::count(Receptor)
colnames(ligand.count)<-c("Cell_Type", "L.Count")
ligand.count$L.Count<-0
colnames(receptor.count)<-c("Cell_Type", "R.Count")
receptor.count$R.Count<-0

WT.count<-as.data.frame(str_split_fixed(WT$lr.cell.name, pattern = "-", n = 2))
colnames(WT.count)<-c("Ligand", "Receptor")
L.Count<-WT.count%>%dplyr::count(Ligand)
colnames(L.Count)<-c("Cell_Type", "L.Count")
ligand.count<-left_join(ligand.count, L.Count, by= "Cell_Type")
ligand.count<- ligand.count%>%dplyr::select(-L.Count.x)
colnames(ligand.count)<-c("Cell_Type", "Count")
write.csv(ligand.count, "WT.ligand.count.csv")
#ligand.count<- replace(ligand.count, is.na(ligand.count$Count), 0)
ligand.count$unit<-"ligand"

R.Count<-WT.count%>%dplyr::count(Receptor)
colnames(R.Count)<-c("Cell_Type", "R.Count")
receptor.count<-left_join(receptor.count, R.Count, by= "Cell_Type")
receptor.count<- receptor.count%>%dplyr::select(-R.Count.x)
colnames(receptor.count)<-c("Cell_Type", "Count")
write.csv(receptor.count, "WT.receptor.count.csv")
#receptor.count<- replace(receptor.count, is.na(ligand.count$Count), 0)
receptor.count$unit<-"receptor"

WT.count.final<-rbind(ligand.count, receptor.count)

ggplot(WT.count.final, aes(x = Cell_Type, y =Count, fill = unit)) +
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("WT Ligand & Receptor Counts")

#Graphing Network: 
WT_igraph<-graph(c("CM","CM", 
                   "CM","FB",
                   "CM","EC",
                   "CM","PC",
                   "CM","SM",
                   "CM","MAC",
                   "CM","MES",
                   "CM","LEC", 
                   "CM","BLYMPH",
                   "CM","TLYMPH",
                   "CM","ADIPO",
                   "CM","DC",
                   
                   "FB","CM", 
                   "FB","FB",
                   "FB","EC",
                   "FB","PC",
                   "FB","SM",
                   "FB","MAC",
                   "FB","MES",
                   "FB","LEC", 
                   "FB","BLYMPH",
                   "FB","TLYMPH",
                   "FB","ADIPO",
                   "FB","DC",
                   
                   "EC","CM", 
                   "EC","FB",
                   "EC","EC",
                   "EC","PC",
                   "EC","SM",
                   "EC","MAC",
                   "EC","MES",
                   "EC","LEC", 
                   "EC","BLYMPH",
                   "EC","TLYMPH",
                   "EC","ADIPO",
                   "EC","DC",
                   
                   "PC","CM", 
                   "PC","FB",
                   "PC","EC",
                   "PC","PC",
                   "PC","SM",
                   "PC","MAC",
                   "PC","MES",
                   "PC","LEC", 
                   "PC","BLYMPH",
                   "PC","TLYMPH",
                   "PC","ADIPO",
                   "PC","DC",
                   
                   "SM","CM", 
                   "SM","FB",
                   "SM","EC",
                   "SM","PC",
                   "SM","SM",
                   "SM","MAC",
                   "SM","MES",
                   "SM","LEC", 
                   "SM","BLYMPH",
                   "SM","TLYMPH",
                   "SM","ADIPO",
                   "SM","DC",
                   
                   "MAC","CM", 
                   "MAC","FB",
                   "MAC","EC",
                   "MAC","PC",
                   "MAC","SM",
                   "MAC","MAC",
                   "MAC","MES",
                   "MAC","LEC", 
                   "MAC","BLYMPH",
                   "MAC","TLYMPH",
                   "MAC","ADIPO",
                   "MAC","DC",
                   
                   "MES","CM", 
                   "MES","FB",
                   "MES","EC",
                   "MES","PC",
                   "MES","SM",
                   "MES","MAC",
                   "MES","MES",
                   "MES","LEC", 
                   "MES","BLYMPH",
                   "MES","TLYMPH",
                   "MES","ADIPO",
                   "MES","DC",
                   
                   "LEC","CM", 
                   "LEC","FB",
                   "LEC","EC",
                   "LEC","PC",
                   "LEC","SM",
                   "LEC","MAC",
                   "LEC","MES",
                   "LEC","LEC", 
                   "LEC","BLYMPH",
                   "LEC","TLYMPH",
                   "LEC","ADIPO",
                   "LEC","DC",
                   
                   "BLYMPH","CM", 
                   "BLYMPH","FB",
                   "BLYMPH","EC",
                   "BLYMPH","PC",
                   "BLYMPH","SM",
                   "BLYMPH","MAC",
                   "BLYMPH","MES",
                   "BLYMPH","LEC", 
                   "BLYMPH","BLYMPH",
                   "BLYMPH","TLYMPH",
                   "BLYMPH","ADIPO",
                   "BLYMPH","DC",
                   
                   "TLYMPH","CM", 
                   "TLYMPH","FB",
                   "TLYMPH","EC",
                   "TLYMPH","PC",
                   "TLYMPH","SM",
                   "TLYMPH","MAC",
                   "TLYMPH","MES",
                   "TLYMPH","LEC", 
                   "TLYMPH","BLYMPH",
                   "TLYMPH","TLYMPH",
                   "TLYMPH","ADIPO",
                   "TLYMPH","DC",
                   
                   "ADIPO","CM", 
                   "ADIPO","FB",
                   "ADIPO","EC",
                   "ADIPO","PC",
                   "ADIPO","SM",
                   "ADIPO","MAC",
                   "ADIPO","MES",
                   "ADIPO","LEC", 
                   "ADIPO","BLYMPH",
                   "ADIPO","TLYMPH",
                   "ADIPO","ADIPO",
                   "ADIPO","DC",
                   
                   "DC","CM", 
                   "DC","FB",
                   "DC","EC",
                   "DC","PC",
                   "DC","SM",
                   "DC","MAC",
                   "DC","MES",
                   "DC","LEC", 
                   "DC","BLYMPH",
                   "DC","TLYMPH",
                   "DC","ADIPO",
                   "DC","DC"))

WT_coord <- layout_in_circle(WT_igraph, order=V(WT_igraph))
E(WT_igraph)$width <- c(count$counts) 
V(WT_igraph)$color <- c("azure4", "red", "green", 
                        "lightseagreen", "orange", "violet", 
                        "cyan", "lightcoral", "blue", 
                        "purple", "pink", "burlywood4") #colors for # of cell types

edge.start <- ends(WT_igraph, es=E(WT_igraph), names=F)[,1]
edge.col <- V(WT_igraph)$color[edge.start]
#360 degrees/# cell types --> convert to radians
#Number of 0's = # cell types
edge_ANGLE <- c(6.283,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0, 
                5.759,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                5.235, 	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                4.712,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                4.188,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                3.665,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                3.141,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,	
                2.617,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,	
                2.094,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0, 	
                1.570,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                1.047,  0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                0.523)


WT.plot<-plot(WT_igraph, 
              layout=WT_coord, 
              vertex.size=30, vertex.frame.color="black", vertex.label.color="black", vertex.label.cex=0.8, edge.curved=-0.1, edge.color=edge.col, edge.arrow.size=0, dge.arrow.width=0, edge.loop.angle=edge_ANGLE, main="WT")
################################################################################
#LR interactions significantly increased in KO
################################################################################
#Store all potential cell-cell interactions increased in KO
celltype.pairs<-c("CM-CM", 
                  "CM-FB",
                  "CM-EC",
                  "CM-PC",
                  "CM-SM",
                  "CM-MAC",
                  "CM-MES",
                  "CM-LEC", 
                  "CM-BLYMPH",
                  "CM-TLYMPH",
                  "CM-ADIPO",
                  "CM-DC",
                  
                  "FB-CM", 
                  "FB-FB",
                  "FB-EC",
                  "FB-PC",
                  "FB-SM",
                  "FB-MAC",
                  "FB-MES",
                  "FB-LEC", 
                  "FB-BLYMPH",
                  "FB-TLYMPH",
                  "FB-ADIPO",
                  "FB-DC",
                  
                  "EC-CM", 
                  "EC-FB",
                  "EC-EC",
                  "EC-PC",
                  "EC-SM",
                  "EC-MAC",
                  "EC-MES",
                  "EC-LEC", 
                  "EC-BLYMPH",
                  "EC-TLYMPH",
                  "EC-ADIPO",
                  "EC-DC",
                  
                  "PC-CM", 
                  "PC-FB",
                  "PC-EC",
                  "PC-PC",
                  "PC-SM",
                  "PC-MAC",
                  "PC-MES",
                  "PC-LEC", 
                  "PC-BLYMPH",
                  "PC-TLYMPH",
                  "PC-ADIPO",
                  "PC-DC",
                  
                  "SM-CM", 
                  "SM-FB",
                  "SM-EC",
                  "SM-PC",
                  "SM-SM",
                  "SM-MAC",
                  "SM-MES",
                  "SM-LEC", 
                  "SM-BLYMPH",
                  "SM-TLYMPH",
                  "SM-ADIPO",
                  "SM-DC",
                  
                  "MAC-CM", 
                  "MAC-FB",
                  "MAC-EC",
                  "MAC-PC",
                  "MAC-SM",
                  "MAC-MAC",
                  "MAC-MES",
                  "MAC-LEC", 
                  "MAC-BLYMPH",
                  "MAC-TLYMPH",
                  "MAC-ADIPO",
                  "MAC-DC",
                  
                  "MES-CM", 
                  "MES-FB",
                  "MES-EC",
                  "MES-PC",
                  "MES-SM",
                  "MES-MAC",
                  "MES-MES",
                  "MES-LEC", 
                  "MES-BLYMPH",
                  "MES-TLYMPH",
                  "MES-ADIPO",
                  "MES-DC",
                  
                  "LEC-CM", 
                  "LEC-FB",
                  "LEC-EC",
                  "LEC-PC",
                  "LEC-SM",
                  "LEC-MAC",
                  "LEC-MES",
                  "LEC-LEC", 
                  "LEC-BLYMPH",
                  "LEC-TLYMPH",
                  "LEC-ADIPO",
                  "LEC-DC",
                  
                  "BLYMPH-CM", 
                  "BLYMPH-FB",
                  "BLYMPH-EC",
                  "BLYMPH-PC",
                  "BLYMPH-SM",
                  "BLYMPH-MAC",
                  "BLYMPH-MES",
                  "BLYMPH-LEC", 
                  "BLYMPH-BLYMPH",
                  "BLYMPH-TLYMPH",
                  "BLYMPH-ADIPO",
                  "BLYMPH-DC",
                  
                  "TLYMPH-CM", 
                  "TLYMPH-FB",
                  "TLYMPH-EC",
                  "TLYMPH-PC",
                  "TLYMPH-SM",
                  "TLYMPH-MAC",
                  "TLYMPH-MES",
                  "TLYMPH-LEC", 
                  "TLYMPH-BLYMPH",
                  "TLYMPH-TLYMPH",
                  "TLYMPH-ADIPO",
                  "TLYMPH-DC",
                  
                  "ADIPO-CM", 
                  "ADIPO-FB",
                  "ADIPO-EC",
                  "ADIPO-PC",
                  "ADIPO-SM",
                  "ADIPO-MAC",
                  "ADIPO-MES",
                  "ADIPO-LEC", 
                  "ADIPO-BLYMPH",
                  "ADIPO-TLYMPH",
                  "ADIPO-ADIPO",
                  "ADIPO-DC",
                  
                  "DC-CM", 
                  "DC-FB",
                  "DC-EC",
                  "DC-PC",
                  "DC-SM",
                  "DC-MAC",
                  "DC-MES",
                  "DC-LEC", 
                  "DC-BLYMPH",
                  "DC-TLYMPH",
                  "DC-ADIPO",
                  "DC-DC")

#create counts column
celltype.pairs<-as.data.frame(celltype.pairs)
colnames(celltype.pairs)<-c("pairs")
celltype.pairs$counts<-0

#count instances of cell-cell interactions
count<-KO%>%dplyr::count(lr.cell.name)
colnames(count)<-c("pairs", "counts")
#Add counts to dataframe with all potential cell-cell interactions: 
count<-left_join(celltype.pairs, count, by = "pairs" )%>%
  dplyr::select(pairs, counts.y)
count$counts.y[is.na(count$counts.y)] = 0.001 #Make 0's 0.001 for graphing purKOes downstream
colnames(count)<-c("pairs", "counts")
count$counts<-count$counts * 0.1 #scaling edges

#Ligand/Receptor Counts
LR_celltypes<-as.data.frame(str_split_fixed(celltype.pairs$pairs, pattern = "-", n = 2))
colnames(LR_celltypes)<-c("Ligand", "Receptor")

ligand.count<-LR_celltypes%>%dplyr::count(Ligand)
receptor.count<-LR_celltypes%>%dplyr::count(Receptor)
colnames(ligand.count)<-c("Cell_Type", "L.Count")
ligand.count$L.Count<-0
colnames(receptor.count)<-c("Cell_Type", "R.Count")
receptor.count$R.Count<-0

KO.count<-as.data.frame(str_split_fixed(KO$lr.cell.name, pattern = "-", n = 2))
colnames(KO.count)<-c("Ligand", "Receptor")
L.Count<-KO.count%>%dplyr::count(Ligand)
colnames(L.Count)<-c("Cell_Type", "L.Count")
ligand.count<-left_join(ligand.count, L.Count, by= "Cell_Type")
ligand.count<- ligand.count%>%dplyr::select(-L.Count.x)
colnames(ligand.count)<-c("Cell_Type", "Count")
write.csv(ligand.count, "KO.ligand.count.csv")
#ligand.count<- replace(ligand.count, is.na(ligand.count$Count), 0)
ligand.count$unit<-"ligand"

R.Count<-KO.count%>%dplyr::count(Receptor)
colnames(R.Count)<-c("Cell_Type", "R.Count")
receptor.count<-left_join(receptor.count, R.Count, by= "Cell_Type")
receptor.count<- receptor.count%>%dplyr::select(-R.Count.x)
colnames(receptor.count)<-c("Cell_Type", "Count")
write.csv(receptor.count, "KO.receptor.count.csv")
#receptor.count<- replace(receptor.count, is.na(ligand.count$Count), 0)
receptor.count$unit<-"receptor"

KO.count.final<-rbind(ligand.count, receptor.count)

ggplot(KO.count.final, aes(x = Cell_Type, y =Count, fill = unit)) +
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("KO Ligand & Receptor Counts")


#Graphing Network: 
KO_igraph<-graph(c("CM","CM", 
                   "CM","FB",
                   "CM","EC",
                   "CM","PC",
                   "CM","SM",
                   "CM","MAC",
                   "CM","MES",
                   "CM","LEC", 
                   "CM","BLYMPH",
                   "CM","TLYMPH",
                   "CM","ADIPO",
                   "CM","DC",
                   
                   "FB","CM", 
                   "FB","FB",
                   "FB","EC",
                   "FB","PC",
                   "FB","SM",
                   "FB","MAC",
                   "FB","MES",
                   "FB","LEC", 
                   "FB","BLYMPH",
                   "FB","TLYMPH",
                   "FB","ADIPO",
                   "FB","DC",
                   
                   "EC","CM", 
                   "EC","FB",
                   "EC","EC",
                   "EC","PC",
                   "EC","SM",
                   "EC","MAC",
                   "EC","MES",
                   "EC","LEC", 
                   "EC","BLYMPH",
                   "EC","TLYMPH",
                   "EC","ADIPO",
                   "EC","DC",
                   
                   "PC","CM", 
                   "PC","FB",
                   "PC","EC",
                   "PC","PC",
                   "PC","SM",
                   "PC","MAC",
                   "PC","MES",
                   "PC","LEC", 
                   "PC","BLYMPH",
                   "PC","TLYMPH",
                   "PC","ADIPO",
                   "PC","DC",
                   
                   "SM","CM", 
                   "SM","FB",
                   "SM","EC",
                   "SM","PC",
                   "SM","SM",
                   "SM","MAC",
                   "SM","MES",
                   "SM","LEC", 
                   "SM","BLYMPH",
                   "SM","TLYMPH",
                   "SM","ADIPO",
                   "SM","DC",
                   
                   "MAC","CM", 
                   "MAC","FB",
                   "MAC","EC",
                   "MAC","PC",
                   "MAC","SM",
                   "MAC","MAC",
                   "MAC","MES",
                   "MAC","LEC", 
                   "MAC","BLYMPH",
                   "MAC","TLYMPH",
                   "MAC","ADIPO",
                   "MAC","DC",
                   
                   "MES","CM", 
                   "MES","FB",
                   "MES","EC",
                   "MES","PC",
                   "MES","SM",
                   "MES","MAC",
                   "MES","MES",
                   "MES","LEC", 
                   "MES","BLYMPH",
                   "MES","TLYMPH",
                   "MES","ADIPO",
                   "MES","DC",
                   
                   "LEC","CM", 
                   "LEC","FB",
                   "LEC","EC",
                   "LEC","PC",
                   "LEC","SM",
                   "LEC","MAC",
                   "LEC","MES",
                   "LEC","LEC", 
                   "LEC","BLYMPH",
                   "LEC","TLYMPH",
                   "LEC","ADIPO",
                   "LEC","DC",
                   
                   "BLYMPH","CM", 
                   "BLYMPH","FB",
                   "BLYMPH","EC",
                   "BLYMPH","PC",
                   "BLYMPH","SM",
                   "BLYMPH","MAC",
                   "BLYMPH","MES",
                   "BLYMPH","LEC", 
                   "BLYMPH","BLYMPH",
                   "BLYMPH","TLYMPH",
                   "BLYMPH","ADIPO",
                   "BLYMPH","DC",
                   
                   "TLYMPH","CM", 
                   "TLYMPH","FB",
                   "TLYMPH","EC",
                   "TLYMPH","PC",
                   "TLYMPH","SM",
                   "TLYMPH","MAC",
                   "TLYMPH","MES",
                   "TLYMPH","LEC", 
                   "TLYMPH","BLYMPH",
                   "TLYMPH","TLYMPH",
                   "TLYMPH","ADIPO",
                   "TLYMPH","DC",
                   
                   "ADIPO","CM", 
                   "ADIPO","FB",
                   "ADIPO","EC",
                   "ADIPO","PC",
                   "ADIPO","SM",
                   "ADIPO","MAC",
                   "ADIPO","MES",
                   "ADIPO","LEC", 
                   "ADIPO","BLYMPH",
                   "ADIPO","TLYMPH",
                   "ADIPO","ADIPO",
                   "ADIPO","DC",
                   
                   "DC","CM", 
                   "DC","FB",
                   "DC","EC",
                   "DC","PC",
                   "DC","SM",
                   "DC","MAC",
                   "DC","MES",
                   "DC","LEC", 
                   "DC","BLYMPH",
                   "DC","TLYMPH",
                   "DC","ADIPO",
                   "DC","DC"))

KO_coord <- layout_in_circle(KO_igraph, order=V(KO_igraph))
E(KO_igraph)$width <- c(count$counts) 
V(KO_igraph)$color <- c("azure4", "red", "green", 
                        "lightseagreen", "orange", "violet", 
                        "cyan", "lightcoral", "blue", 
                        "purple", "pink", "burlywood4") #colors for # of cell types

edge.start <- ends(KO_igraph, es=E(KO_igraph), names=F)[,1]
edge.col <- V(KO_igraph)$color[edge.start]
#360 degrees/# cell types --> convert to radians
#Number of 0's = # cell types
edge_ANGLE <- c(6.283,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0, 
                5.759,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                5.235, 	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                4.712,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                4.188,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                3.665,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                3.141,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,	
                2.617,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,	
                2.094,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0, 	
                1.570,	0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                1.047,  0,	0,	0,	0,	0,	0,	0,	0, 0, 0, 0, 0,
                0.523)

KO.plot<-plot(KO_igraph, 
              layout=KO_coord, 
              vertex.size=30, vertex.frame.color="black", vertex.label.color="black", vertex.label.cex=0.8, edge.curved=-0.1, edge.color=edge.col, edge.arrow.size=0, dge.arrow.width=0, edge.loop.angle=edge_ANGLE, main="KO")

