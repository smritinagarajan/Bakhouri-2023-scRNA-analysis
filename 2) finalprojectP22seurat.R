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
P22T.data <- Read10X(data.dir = "/xdisk/darrenc/cmm_523/nagarajan/cmm523hw8/P22T/")

# P11T <- CreateSeuratObject(counts = P11T.data, project = "CD3_T_I", min.cells = 10, min.features = 200)
# Error received: in validObject(.Object) :  invalid class “LogMap” object: Duplicate rownames not allowed

sum(duplicated(colnames(P22T.data)))
# = 32 duplicated column names

colnames(P22T.data) <- make.unique(colnames(P22T.data))
#fixed the problem

P22T <- CreateSeuratObject(counts = P22T.data, project = "CD3_T_I", min.cells = 50, min.features = 200)
#worked!!

P22T
# An object of class Seurat 
# 7897 features across 2399 samples within 1 assay 
# Active assay: RNA (7897 features, 0 variable features)
# 1 layer present: counts

head(P22T@meta.data, 5)
# no percent.mt?

#not replicating QC measures. the original authors already filtered out cells.

VlnPlot(P22T, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(P22T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#chose not to subset further. has already been done.

P22T <- NormalizeData(P22T, normalization.method = "LogNormalize", scale.factor = 10000)

sub_P22T <- FindVariableFeatures(sub_P22T, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sub_P22T), 10)

#plotting variable features
plot1 <- VariableFeaturePlot(sub_P22T)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(sub_P22T)
sub_P22T <- ScaleData(sub_P22T, features = all.genes)

sub_P22T <- RunPCA(sub_P22T, features = VariableFeatures(object = sub_P22T))

#visualize defining cells and features
VizDimLoadings(sub_P22T, dims = 1:2, reduction = "pca")
DimPlot(sub_P22T, dims = c(1,2), reduction = "pca") + NoLegend()
DimHeatmap(sub_P22T, dims = 1:9, cells = 500, balanced = TRUE)

ElbowPlot(sub_P22T, ndims = 50)
#to select how many PCs are needed. I'm saying 30.


#cluster cells
sub_P22T <- FindNeighbors(sub_P22T, dims = 1:10)
sub_P22T <- FindClusters(sub_P22T, resolution = 0.5)

sub_P22T <- RunUMAP(sub_P22T, dims = 1:10)
DimPlot(sub_P22T, reduction = "umap")

sub_P22T <- subset(P22T, idents = c("0","5","6","7"), invert = TRUE)
DimPlot(sub_P22T, reduction = "umap")


sub_P22T.markers <- FindAllMarkers(sub_P22T, only.pos = TRUE)

sub_P22T.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_min(p_val, n = 3) %>%
  print(n=24)

# print n = 24 was added to see all the data in the console

VlnPlot(sub_P22T, features = c("CD4", "CD8A", "CD8B", "CD45RA", "CCR7", "PDCD1", "IL2RA", "FOXP3", "GZMA"))

FeaturePlot(sub_P22T, features = c("CD4", "CD8A", "CD8B", "CD45RA", "CCR7", "PDCD1", "IL2RA", "FOXP3", "GZMA"))

sub_P22T.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(sub_P22T, features = top10$gene) + NoLegend()

#SingleR for labelling clusters, using celldex's DatabaseImmuneCellExpressionData as reference

raw_counts <- LayerData(sub_P22T, assay = "RNA", layer = 'counts')
norm_counts <- LayerData(sub_P22T, assay = "RNA", layer = 'data')
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
sub_P22T <- AddMetaData(sub_P22T, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

# Visualise them on the UMAP
sub_P22T <- SetIdent(sub_P22T, value = "SingleR_HCA")
DimPlot(sub_P22T, label = T , repel = T, label.size = 5) + NoLegend()

saveRDS(P22T, file = "../sub_P22T_finalSeurat.rds")
