#Script for scRNA analysis using seurat fr NSCLC 

#Loading libraries

install.packages('Seurat')
install.packages('tidyverse')
install.packages("ggplot2")
install.packages("hdf5r")
install.packages("harmony")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(hdf5r)
library(harmony)

#Loading the dataset 

lc<-Read10X_h5(file = "/Users/sahanabaskar/Desktop/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")

#Selecting the gene expression from the modalities

str(lc)
lc_counts<-lc$`Gene Expression`

#Creating seurat object

lc.obj<-CreateSeuratObject(counts = lc_counts,
                           project = "NSCLC",
                           min.cells = 8,
                           min.features = 200 )

View(lc.obj@meta.data)

#Quality Control--- 

#mitochondrial reads
lc.obj[["percent.mt"]]<-PercentageFeatureSet(lc.obj, pattern = '^MT-')

VlnPlot(lc.obj, features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol=3)

plot1 <- FeatureScatter(lc.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  geom_smooth(method = 'lm')
plot2 <- FeatureScatter(lc.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')

CombinePlots(plots = list(plot1, plot2))

#Filtering

lc.obj<-subset(lc.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization of data

lc.obj<-NormalizeData(lc.obj)

#Highly variable features
lc.obj<-FindVariableFeatures(lc.obj,selection.method = 'vst', nfeatures = 2000)

#top 10 higly variable genes
mostv<-head(VariableFeatures(lc.obj),10)

plot3 <- VariableFeaturePlot(lc.obj)
plot4 <- LabelPoints(plot = plot3, points = mostv, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))

#Scaling
all.genes<-rownames(lc.obj)
lc.obj<-ScaleData(lc.obj, features = all.genes)

#Linear dimensionality reduction
lc.obj<-RunPCA(lc.obj, features = VariableFeatures(object = lc.obj))

print(lc.obj[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(lc.obj, dims=1, cells = 500, balanced = TRUE)

#dimensionality of data
ElbowPlot(lc.obj)

#Clustering
lc.obj<-FindNeighbors(lc.obj, dims = 1:15)
lc.obj<-FindClusters(lc.obj, resolution = c(0.1,0.3,0.5,0.7,1))
View(lc.obj@meta.data)

DimPlot(lc.obj, group.by = 'RNA_snn_res.0.1', label = TRUE)

Idents(lc.obj)<-'RNA_snn_res.0.1'
Idents(lc.obj)

#UMAP
lc.obj<-RunUMAP(lc.obj, dims = 1:15)

DimPlot(lc.obj,reduction = 'umap')

#tSNE
lc.obj <- RunTSNE(lc.obj, dims = 1:15)
DimPlot(lc.obj, reduction = "tsne")

# differentially expressed features for each cluster
lc.markers <- FindMarkers(lc.obj, ident.1 = 0, ident.2 = 1:5, test.use = 'wilcox')

lc.markers <- FindMarkers(lc.obj, ident.1 = 2, min.pct = 0.25)
head(lc.markers, n = 5)

top.markers <- as.data.frame(lc.obj@markers)
write.table(top.markers, file = "markers.txt", sep = "\t")

saveRDS(lc.obj, "/Users/sahanabaskar/Desktop/Rshiny/seuratObj.rds",version = 3)

file.info("/Users/sahanabaskar/Desktop/Rshiny/seurat_object.rds")

file.exists("/Users/sahanabaskar/Desktop/Rshiny/seurat_object.rds")
