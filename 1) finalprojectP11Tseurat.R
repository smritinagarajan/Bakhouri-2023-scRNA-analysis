rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) 
# avoid truncated output in R console and scientific notation
set.seed(42)

dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

library(devtools)
library(png)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)
library(celldex)
library(SingleR)
library(viridis)
library(pheatmap)
library(matrixStats)
library(SeuratData)
library(scRNAseq)
library(scran)
library(dbplyr)

# Load the PBMC dataset
P11T.data <- Read10X(data.dir = "/xdisk/darrenc/cmm_523/nagarajan/cmm523hw8/P11T/")
# Initialize the Seurat object with the raw (non-normalized data).

# P11T <- CreateSeuratObject(counts = P11T.data, project = "CD3_T_I", min.cells = 10, min.features = 200)
# Error received: in validObject(.Object) :  invalid class “LogMap” object: Duplicate rownames not allowed

sum(duplicated(colnames(P11T.data)))
# = 32 duplicated column names

colnames(P11T.data) <- make.unique(colnames(P11T.data))
#fixed the problem

duplicate <- names(as.data.frame(P11T.data[,duplicated(colnames(P11T.data))]))
#asked to make a list of duplicated colnames, and received a null dataframe.

P11T <- CreateSeuratObject(counts = P11T.data, project = "CD3_T_I", min.cells = 50, min.features = 200)
#worked!!

P11T
# An object of class Seurat 
# 14035 features across 7915 samples within 1 assay 
# Active assay: RNA (14035 features, 0 variable features)
# 1 layer present: counts

head(P11T@meta.data, 5)
# no percent.mt?

VlnPlot(P11T, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(P11T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#chose not to subset further. has already been done.

P11T <- NormalizeData(P11T, normalization.method = "LogNormalize", scale.factor = 10000)

P11T <- FindVariableFeatures(P11T, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(P11T), 10)

#plotting variable features
plot1 <- VariableFeaturePlot(P11T)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(P11T)
P11T <- ScaleData(P11T, features = all.genes)

P11T <- RunPCA(P11T, features = VariableFeatures(object = P11T))

#visualize defining cells and features
VizDimLoadings(P11T, dims = 1:2, reduction = "pca")
DimPlot(P11T, dims = c(1,2), reduction = "pca") + NoLegend()
DimHeatmap(P11T, dims = 1:9, cells = 500, balanced = TRUE)

ElbowPlot(P11T, ndims = 50)
#to select how many PCs are needed. I'm saying 30.

#cluster cells
P11T <- FindNeighbors(P11T, dims = 1:20)
P11T <- FindClusters(P11T, resolution = 0.5)

P11T <- RunUMAP(P11T, dims = 1:20)
DimPlot(P11T, reduction = "umap")

saveRDS(P11T, file = "../P11T.rds")

#took a break

P11T <- readRDS("P11T.rds")
#to reload the Seurat object

P11T.markers <- FindAllMarkers(P11T, only.pos = TRUE)

P11T.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_min(p_val, n = 3) %>%
  print(n=21)

# print n = 21 was added to see all the data in the console

VlnPlot(P11T, features = c("CD4", "CD8A", "CD8B", "CD45RA", "CCR7", "PDCD1", "IL2RA", "FOXP3", "GZMA"))

FeaturePlot(P11T, features = c("CD4", "CD8A", "CD8B", "CD45RA", "CCR7", "PDCD1", "IL2RA", "FOXP3", "GZMA"))

P11T.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(P11T, features = top10$gene) + NoLegend()

saveRDS(P11T, file = "../P11T_finalSeurat.rds")

#SingleR for labelling clusters, using celldex's DatabaseImmuneCellExpressionData as reference

raw_counts <- LayerData(P11T, assay = "RNA", layer = 'counts')
norm_counts <- LayerData(P11T, assay = "RNA", layer = 'data')
raw_counts[c('CD4', 'CD8A', 'IL2RA', 'CCR7'),1:5]
norm_counts[c('CD4', 'CD8A', 'IL2RA', 'CCR7'),1:5]

ref <- celldex::DatabaseImmuneCellExpressionData()
unique(ref$label.fine)
# unique(ref$label.main)
#want fine labels, not broader main labels

# ref <- ref[,grepl('Tcells, CD4+ | Tcells, CD8+', ref$label.main)]
# unique(ref$label.main)
# not working

# try instead =
ref <- ref[,ref$label.main %in% c("T cells, CD4+", "T cells, CD8+")]

ct_ann <- SingleR(test = norm_counts,
                  ref = ref, 
                  labels = ref$label.fine,
                  de.method = 'wilcox')
ct_ann %>% head()

#SingleR prediction QC
plotScoreHeatmap(ct_ann)
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

# Add singleR labels to seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs
P11T <- AddMetaData(P11T, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

# Visualise them on the UMAP
P11T <- SetIdent(P11T, value = "SingleR_HCA")
DimPlot(P11T, label = T , repel = T, label.size = 5) + NoLegend()
