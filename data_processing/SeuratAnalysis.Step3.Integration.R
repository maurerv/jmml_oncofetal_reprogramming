# Data integration
library(Seurat)

# Finding anchors
DefaultAssay(tempSeurat) <- 'RNA'
seulist <- SplitObject(tempSeurat, split.by = "Dataset")
for (i in 1:length(seulist)) {
  seulist[[i]] <- NormalizeData(seulist[[i]])
  seulist[[i]] <- FindVariableFeatures(seulist[[i]])
}
anchors <- FindIntegrationAnchors(seulist, dims = 1:30)
intSeu <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(intSeu) <- "integrated"
tempSeurat <- intSeu

# Data processing
tempSeurat <-
  ScaleData(tempSeurat, verbose = FALSE, vars.to.regress = 'nCount_RNA')
tempSeurat <- RunPCA(tempSeurat, npcs = 100, verbose = TRUE)

findElbow <- function(x, y, threshold) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > threshold)
  return(indices)
}
x = tempSeurat@reductions$pca@stdev
y = 1:100
indices <- findElbow(x, y, 1e4)
pc = indices[1]
if (pc < 5)
{
  pc = indices[2]
}
tempSeurat <-
  FindNeighbors(tempSeurat, reduction = "pca", dims = 1:pc)
resolution <- c(0.1, 0.2, 0.4, 0.6, 0.8)
tempSeurat <- FindClusters(tempSeurat, resolution = resolution)
tempSeurat <- RunUMAP(tempSeurat, dims = 1:pc)