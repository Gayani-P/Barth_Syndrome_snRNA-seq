################################################################################
LIB='/cluster/tufts/hpc/tools/R/4.0.0' 
REPO="http://cran.us.r-project.org"

.libPaths(c("",LIB)) 
.libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0", .libPaths()))

library(Seurat)
library(cowplot)
library(Matrix)

setwd("/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/5_CellAssignment")
#https://github.com/satijalab/seurat/issues/1883
################################################################################
#Load Data(Choose which dimension reduction to generate umap from): 
obj<-readRDS(file = "/cluster/tufts/chinlab/gperer01/to_send/Liam/V2/3_Harmony/out/V2_3_Harmony_18PCs.rds")

#Choose Resolution to use(Picked by ChooseR): 
Idents(object = obj)<- "RNA_snn_res.0.9"

#Determine the number of cells per cluster per genotype
n_cells <- FetchData(obj, 
                     vars = c("ident", "mouse")) %>%
  dplyr::count(ident, mouse) %>%
  tidyr::spread(ident, n)
n_cells

#UMAP
umap_original<-DimPlot(obj, label = TRUE, raster = FALSE)
umap_original
ggsave("images/UMAP_18PCs_0.9res.png")

#UMAP by Genotype 
DimPlot(obj, label = TRUE, split.by = "genotype")  + NoLegend()
ggsave("images/UMAP_18PCs_0.9res_by_genotype.png")
################################################################################  
#Exploratory Plots: Save As 2000 x 2000 png 
################################################################################
CM.markers = c("Mybpc3", "Myh7b", "Myh7", "Tnnt2",
               "Ldb3", "Fhod3", "Pln", "Actc1", 
               "Myocd", "Mhrt", "Nppa", "Myh6", 
               "Atp2a2", "Ryr2", "Myl2")

CM.dotplot<-DotPlot(obj,  features = c(CM.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
CM.dotplot
ggsave("images/CM.markers.png")
################################################################################             
FB.markers = c("Fmo2", "Postn", "Col6a3", "Col3a1",
               "Fbn1", "Col1a1", "Ddr2", "Plekhh2", 
               "Igf1", "Scn7a", "Col1a2", "Pdgfra")

FB.dotplot<-DotPlot(obj,  features = c(FB.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
FB.dotplot
ggsave("images/FB.markers.png")
#################################################################################
EC.markers = c("Pecam1", "Vwf", "Aqp1", "Cdh5", "Tek", "Kdr", "Sox7", "Flt4")

EC.dotplot<-DotPlot(obj,  features = c(EC.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
EC.dotplot
ggsave("images/EC.markers.png")
#################################################################################
PC.markers = c("Pdgfrb", "Abcc9", "Kcnj8", "Vtn", "Steap4", "Notch3")

PC.dotplot<-DotPlot(obj,  features = c(PC.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
PC.dotplot
ggsave("images/PC.markers.png")
################################################################################
SM.markers = c("Myh11", "Acta2", "Tagln", "Cnn1", "Lmod1")

SM.dotplot<-DotPlot(obj,  features = c(SM.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
SM.dotplot
ggsave("images/SM.markers.png")
################################################################################
MAC.markers = c("Ccr5", "Cd163", "Cd14", "Fcgr3","Csf1r")

MAC.dotplot<-DotPlot(obj,  features = c(MAC.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
MAC.dotplot
ggsave("images/MAC.markers.png")
################################################################################
MES.markers = c( "Wt1", "Espn", "Lrrn4", "Upk3b")

MES.dotplot<-DotPlot(obj,  features = c(MES.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
MES.dotplot
ggsave("images/MES.markers.png")
################################################################################
LEC.markers = c("Lyve1", "Mmrn1", "Ccl21a", "Cldn5", "Flt4")

LEC.dotplot<-DotPlot(obj,  features = c(LEC.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
LEC.dotplot
ggsave("images/LEC.markers.png")
################################################################################
ADIPO.markers = c("Plin1", "Acaca", "Fasn", "Scd1")

ADIPO.dotplot<-DotPlot(obj,  features = c(ADIPO.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
ADIPO.dotplot
ggsave("images/ADIPO.markers.png")
################################################################################
BLYMPH.markers = c("Cd19", "Cd22", "Ighd", "Cd79a")

BLYMPH.dotplot<-DotPlot(obj,  features = c(BLYMPH.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
BLYMPH.dotplot
ggsave("images/BLYMPH.markers.png")
################################################################################
TLYMPH.markers = c("Skap1", "Itk", "Cd3e", "Themis")

TLYMPH.dotplot<-DotPlot(obj,  features = c(TLYMPH.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
TLYMPH.dotplot
ggsave("images/TLYMPH.markers.png")
################################################################################
DC.markers = c("Itgax", "Itgam", "Cd24a", "Kit", "Flt3")

DC.dotplot<-DotPlot(obj,  features = c(DC.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
DC.dotplot
ggsave("images/DC.markers.png")
################################################################################
NK.markers = c("Klrb1c", "Ncr1")

NK.dotplot<-DotPlot(obj,  features = c(NK.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
NK.dotplot
ggsave("images/NK.markers.png")
################################################################################
NEURO.markers = c("Nrxn1", "Nrxn3", "Plp1", "Cadm2", "Ntm")

NEURO.dotplot<-DotPlot(obj,  features = c(NEURO.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
NEURO.dotplot
ggsave("images/NEURO.markers.png")
################################################################################
CM.markers = c("Mybpc3", "Fhod3", "Myh7", "Tnnt2")
FB.markers = c("Ddr2", "Pdgfra", "Vim", "Col3a1")
EC.markers = c("Pecam1", "Vwf", "Aqp1", "Cdh5")
PC.markers = c("Pdgfrb", "Abcc9", "Kcnj8", "Vtn")
SM.markers = c("Myh11", "Acta2", "Tagln", "Cnn1")
MAC.markers = c("Ccr5", "Cd163", "Fcgr3","Csf1r")
MES.markers = c( "Wt1", "Espn", "Lrrn4", "Upk3b")
LEC.markers = c("Lyve1", "Mmrn1", "Ccl21a", "Flt4")
ADIPO.markers = c("Plin1", "Acaca", "Fasn", "Scd1")
BLYMPH.markers = c("Cd19", "Cd22", "Ighd", "Cd79a")
TLYMPH.markers = c("Skap1", "Itk", "Cd3e", "Themis")
DC.markers = c("Itgax", "Itgam", "Cd24a", "Flt3")
NK.markers = c("Klrb1c", "Ncr1")
NEURO.markers = c("Nrxn1", "Nrxn3", "Plp1", "Ntm")

#Dotplot by Cluster: 
dotplot<-DotPlot(obj,  features = c(CM.markers, 
                                    EC.markers, 
                                    FB.markers,
                                    SM.markers,
                                    PC.markers,
                                    NEURO.markers,
                                    LEC.markers,
                                    BLYMPH.markers,
                                    TLYMPH.markers,
                                    NK.markers,
                                    ADIPO.markers,
                                    MAC.markers,
                                    DC.markers,
                                    MES.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
plot_grid(dotplot, umap_original)
################################################################################
Idents(object = obj)<- "RNA_snn_res.0.9"
obj<- RenameIdents(object = obj, 
                   "0"="CM", 
                   "1"="FB", 
                   "2"="EC", 
                   "3"="FB", 
                   "4"="CM", 
                   "5"="CM", 
                   "6"="EC",
                   "7"="FB", 
                   "8"="MYL", 
                   "9"="PC", 
                   "10"="MYL", 
                   "11"="CM", 
                   "12"="FB",
                   "13"="CM", 
                   "14"="MES", 
                   "15"="FB", 
                   "16"="LYM", 
                   "17"="MYL",
                   "18"="CM",
                   "19"="PC",
                   "20"="LEC",
                   "21"="MYL", 
                   "22"="MYL", 
                   "23"="PC", 
                   "24"="SM", 
                   "25"="MYL", 
                   "26"= "CM",
                   "27"= "FB",
                   "28"= "MYL", 
                   "29"= "ADIPO",
                   "30"= "CM",
                   "31"= "CM",
                   "32"= "PC",
                   "33"= "CM",
                   "34"= "FB")   

#Dotplot by Cell Type: 
dotplot<-DotPlot(obj,  features = c(CM.markers, 
                                    EC.markers, 
                                    FB.markers,
                                    SM.markers,
                                    PC.markers,
                                    NEURO.markers,
                                    LEC.markers,
                                    BLYMPH.markers,
                                    TLYMPH.markers,
                                    NK.markers,
                                    ADIPO.markers,
                                    MAC.markers,
                                    DC.markers,
                                    MES.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()

umap<-DimPlot(obj, raster = FALSE, label= TRUE)

plot_grid(dotplot, umap)
################################################################################
#LYM Subcluster
obj<- FindSubCluster(obj, 
                     cluster = "LYM", 
                     graph.name = "RNA_snn", 
                     subcluster.name = "LYM_subcluster",  
                     resolution = c(0.9), 
                     algorithm = 1)


sub.umap<-DimPlot(obj, reduction = "umap", group.by = "LYM_subcluster", label = FALSE, label.size = 6)
sub.umap

dotplot<-DotPlot(obj,  features = c(CM.markers, 
                                    EC.markers, 
                                    FB.markers,
                                    SM.markers,
                                    PC.markers,
                                    NEURO.markers,
                                    LEC.markers,
                                    BLYMPH.markers,
                                    TLYMPH.markers,
                                    NK.markers,
                                    ADIPO.markers,
                                    MAC.markers,
                                    DC.markers,
                                    MES.markers), group.by = "LYM_subcluster")  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()

plot_grid(sub.umap, dotplot)

#Add new metadata column for reassigned subclusters: 
meta<-obj@meta.data
meta$LYM_subcluster
dim<-dim(meta)

for (cell in 1:dim[1]) { 
  if (meta$LYM_subcluster[cell] == "LYM_0"){meta$LYM_subcluster[cell] <-"MYL"} else if
  (meta$LYM_subcluster[cell] == "LYM_1"){meta$LYM_subcluster[cell] <-"MYL"} else if
  (meta$LYM_subcluster[cell] == "LYM_2"){meta$LYM_subcluster[cell] <-"BLYMPH"} else if
  (meta$LYM_subcluster[cell] == "LYM_3"){meta$LYM_subcluster[cell] <-"FB"} else if
  (meta$LYM_subcluster[cell] == "LYM_4"){meta$LYM_subcluster[cell] <-"EC"} else if
  (meta$LYM_subcluster[cell] == "LYM_5"){meta$LYM_subcluster[cell] <-"TLYMPH"} else if
  (meta$LYM_subcluster[cell] == "LYM_6"){meta$LYM_subcluster[cell] <-"BLYMPH"} else if
  (meta$LYM_subcluster[cell] == "LYM_7"){meta$LYM_subcluster[cell] <-"MYL"} else if
  (meta$LYM_subcluster[cell] == "LYM_8"){meta$LYM_subcluster[cell] <-"TLYMPH"} else if
  (meta$LYM_subcluster[cell] == "LYM_9"){meta$LYM_subcluster[cell] <-"EC"} else if
  (meta$LYM_subcluster[cell] == "LYM_10"){meta$LYM_subcluster[cell] <-"TLYMPH"}}


meta$LYM_subcluster

obj@meta.data<-meta

#Set LYM_subcluster as active ident:
obj$LYM_subcluster<-as.factor(obj$LYM_subcluster)
obj@active.ident<-obj@meta.data$LYM_subcluster
obj@active.ident
DimPlot(obj, label = TRUE, raster = FALSE)
################################################################################
#MYL Subcluster
Idents(object = obj)<- "LYM_subcluster"
obj<- FindSubCluster(obj, 
                     cluster = "MYL", 
                     graph.name = "RNA_snn", 
                     subcluster.name = "MYL_subcluster",  
                     resolution = c(0.9), 
                     algorithm = 1)


sub.umap<-DimPlot(obj, reduction = "umap", group.by = "MYL_subcluster", label = FALSE, label.size = 6)
sub.umap

dotplot<-DotPlot(obj,  features = c(CM.markers, 
                                    EC.markers, 
                                    FB.markers,
                                    SM.markers,
                                    PC.markers,
                                    NEURO.markers,
                                    LEC.markers,
                                    BLYMPH.markers,
                                    TLYMPH.markers,
                                    NK.markers,
                                    ADIPO.markers,
                                    MAC.markers,
                                    DC.markers,
                                    MES.markers), group.by = "MYL_subcluster")  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()

plot_grid(sub.umap, dotplot)

#Add new metadata column for reassigned subclusters: 
meta<-obj@meta.data
meta$MYL_subcluster
dim<-dim(meta)

for (cell in 1:dim[1]) { 
  if (meta$MYL_subcluster[cell] == "MYL_0"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_1"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_2"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_3"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_4"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_5"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_6"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_7"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_8"){meta$MYL_subcluster[cell] <-"DC"} else if
  (meta$MYL_subcluster[cell] == "MYL_9"){meta$MYL_subcluster[cell] <-"CM"} else if
  (meta$MYL_subcluster[cell] == "MYL_10"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_11"){meta$MYL_subcluster[cell] <-"MAC"} else if
  (meta$MYL_subcluster[cell] == "MYL_12"){meta$MYL_subcluster[cell] <-"EC"} else if
  (meta$MYL_subcluster[cell] == "MYL_13"){meta$MYL_subcluster[cell] <-"EC"} else if
  (meta$MYL_subcluster[cell] == "MYL_14"){meta$MYL_subcluster[cell] <-"LEC"}}


meta$MYL_subcluster

obj@meta.data<-meta

#Set MYL_subcluster as active ident:
obj$MYL_subcluster<-as.factor(obj$MYL_subcluster)
obj@active.ident<-obj@meta.data$MYL_subcluster
obj@active.ident
DimPlot(obj, label = TRUE, raster = FALSE)
################################################################################
#Final Assignments:
obj$final_assignments<-as.factor(obj$MYL_subcluster)
#Set Idents to Final Assignments
Idents(object = obj)<- "final_assignments"
obj@active.ident

final.umap<-DimPlot(obj, label = TRUE, raster = FALSE)
final.umap

dotplot<-DotPlot(obj,  features = c(CM.markers, 
                                    EC.markers, 
                                    FB.markers,
                                    SM.markers,
                                    PC.markers,
                                    NEURO.markers,
                                    LEC.markers,
                                    BLYMPH.markers,
                                    TLYMPH.markers,
                                    NK.markers,
                                    ADIPO.markers,
                                    MAC.markers,
                                    DC.markers,
                                    MES.markers))  + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, size = 10)) +  scale_color_gradient2(low = "blue4", high = "orangered") + coord_flip()
dotplot

#Plot Proportion of Total Number of Cells for each Cell Type in each Genotype: 
pt <- table(Idents(obj), obj$genotype)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

g<-ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") 
g

n_cells <- FetchData(obj, 
                     vars = c("ident", "mouse")) %>%
  dplyr::count(ident, mouse) %>%
  tidyr::spread(ident, n)
n_cells

saveRDS(obj, "V2_5_CellAssignment.rds")
################################################################################