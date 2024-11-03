final project edgeR + monocle3

# ran on the shell
# wget https://raw.githubusercontent.com/cusanovich/cmm523_code/refs/heads/main/new_edgeR_patch.R %>% 
#  -O /xdisk/darrenc/cmm_523/nagarajan/cmm523hw8/new_edgeR_patch.R 

library(edgeR)
source("/xdisk/darrenc/cmm_523/nagarajan/cmm523hw8/new_edgeR_patch.R")

#merge P22T and P11T RDS

byoutcome <- merge(P11T, y = sub_P22T, project = "CD3_T")
byoutcome <- JoinLayers(byoutcome)

byoutcome <- Seurat2PB(byoutcome, sample="outcome", cluster="seurat_clusters")
dim(byoutcome)

summary(byoutcome$samples$lib.size)
keep.samples <- byoutcome$samples$lib.size > 5e5
table(keep.samples)
byoutcome <- byoutcome[, keep.samples]

keep.genes <- filterByExpr(byoutcome, group=byoutcome$samples$cluster)
table(keep.genes)
byoutcome <- byoutcome[keep.genes, , keep=FALSE]

byoutcome <- normLibSizes(byoutcome)
head(byoutcome$samples, n=10L)
summary(byoutcome$samples$norm.factors)

edgeRgraph <- as.factor(byoutcome$samples$cluster)
plotMDS(byoutcome, pch=16, col=c(2:8), main="MDS")
legend("bottomleft", legend=paste0("cluster",levels(cluster)), pch=16, col=2:8, cex=0.8)