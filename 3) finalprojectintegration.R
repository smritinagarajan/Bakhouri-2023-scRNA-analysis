# integrate P11T and P22T using https://satijalab.org/seurat/archive/v3.0/immune_alignment.html

P11T <- AddMetaData(object=P11T, metadata="irAE", col.name='outcome')
sub_P22T <- AddMetaData(object=sub_P22T, metadata="none", col.name='outcome')

immune.anchors <- FindIntegrationAnchors(object.list = list(P11T, sub_P22T), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

#graph
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "outcome")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "outcome")

#SingleR for annotation again

raw_counts <- LayerData(immune.combined, assay = "RNA", layer = 'counts')
norm_counts <- LayerData(immune.combined, assay = "RNA", layer = 'data')
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
immune.combined <- AddMetaData(immune.combined, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

# Visualise them on the UMAP
immune.combined <- SetIdent(immune.combined, value = "SingleR_HCA")
DimPlot(immune.combined, split.by = "outcome", label = T , repel = T, label.size = 5) + NoLegend()

saveRDS(immune.combined, file = "../immune.combined.rds")

