library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

data <- Read10X_h5("Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5")
data1 <- CreateSeuratObject(counts = data,project = "data1",assay = "Spatial")

image <- Read10X_Image(image.dir = "./spatial")
DefaultAssay(object = image) <- 'Spatial'
image <- image[colnames(x=data1)]
data1[['image']] <- image

data1 <- SCTransform(data1,assay = "Spatial", return.only.var.genes = FALSE,verbose = T)

data1 <- RunPCA(data1, assay = "SCT", verbose = FALSE)
ElbowPlot(data1, ndims = 40)
data1 <- FindNeighbors(data1, reduction = "pca", dims = 1:30)
data1 <- FindClusters(data1, verbose = FALSE,resolution = 1)
data1 <- RunUMAP(data1, reduction = "pca", dims = 1:30)

p1 <- DimPlot(data1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(data1, label = TRUE, label.size = 3)
plot_grid(p1, p2)

##########UMAP###############
pdf(file = "dimplot_raw.pdf",width = 6,height = 6)
DimPlot(data1, reduction = "umap",split.by = NULL,label=T, pt.size=0.1, label.size=4)
dev.off()

pdf(file = "Spatial_dimplot.pdf",width = 6,height = 6)
SpatialDimPlot(data1,label=T,label.size=3,pt.size.factor = 2)
dev.off()

DefaultAssay(data1)<-"SCT"
total.markers <- FindAllMarkers(data1, min.pct = 0.5, min.diff.pct = 0.1, logfc.threshold =0.25, only.pos = TRUE, test.use = "wilcox")
write.csv(total.markers,"total.markers.csv")
top10 <- total.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file = "Top10geneHeatmap.pdf",width = 12,height = 14)
DoHeatmap(data1, features = top10$gene,slot="data") + 
  scale_fill_gradientn(colors = c("grey98","gold","orange","red","red3")) + 
  NoLegend()
dev.off()

FeaturePlot(data1,features = c("COL1A1"),label = T,ncol = 1)
VlnPlot(data1,features = c("PTPRC"),pt.size = 0.1,ncol=1)

################Cell annotation###################
cluster_name <- plyr::mapvalues(x = data1@active.ident, from = 0:16, to = c("clu0__Cancer cell","clu1__Cancer cell","clu2__Cancer cell","clu3__Neuron","clu4__Macrophage","clu5__Astrocyte",
                                                                            "clu6__Cancer cell","clu7__Endothelial cell","clu8__Macrophage","clu9__Microglial and Oligodendrocyte","clu10__Cancer cell","Clu11__Fibroblast",
                                                                            "Clu12__Glial cell","Clu13__Cancer cell","clu14__Cancer cell","Clu15__Astrocyte","Clu16__Neuron"))
data1 <- AddMetaData( object = data1, metadata = cluster_name, col.name = "cluster_name")
cluster_name_slim <- factor(unlist(lapply(strsplit(as.character(cluster_name), "__"),function(input){input[2]})),
                            levels = c("Macrophage","Cancer cell","Fibroblast","Endothelial cell","Microglial and Oligodendrocyte","Astrocyte","Glial cell","Neuron"))
data1 <- AddMetaData( object = data1, metadata = cluster_name_slim, col.name = "cluster_name_slim")

Idents(data1) <- "cluster_name_slim"
pdf(file = "dimplot_celltype_slim.pdf",width = 8,height = 6)
DimPlot(data1, reduction = "umap",split.by = NULL,label=F, pt.size=0.1, label.size=4)
dev.off()

pdf(file = "Spatial_dimplot.pdf",width = 6,height = 6)
SpatialDimPlot(data1,label=F,label.size=3,pt.size.factor = 1.5)
dev.off()

markers <- c("SIGLEC1","CD22","CD33","SIGLEC5","SIGLEC7","SIGLEC8","SIGLEC9","SIGLEC10","SIGLEC11")

pdf(file = "Spatialplot.pdf",width = 20,height = 12)
SpatialPlot(data1, features = markers,ncol = 5,pt.size.factor=2)
dev.off()

pdf(file = "VlnPlot.pdf",width = 12,height = 5)
VlnPlot(data1,features = markers,pt.size = 0.1,group.by="cluster_name_slim",stack=T,fill.by="ident")+theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="bold"))+ RotatedAxis()
dev.off()

pdf(file = "DotPlot.pdf",width = 10,height = 4)
DotPlot(data1, features = markers, cols = c("navy", "#F6E921"), dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

saveRDS(data1,"data.rds")
data1 <- readRDS("data.rds")
markers <- c("CD24","CD47","CD274","B2M","SIGLEC10","SIRPA")
pdf(file = "DotPlot.pdf",width = 10,height = 4)
DotPlot(data1, features = markers, cols = c("navy", "#F6E921"), dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()