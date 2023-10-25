# Single cell processing
library(Seurat)
tempSeurat <- NormalizeData(tempSeurat)
tempSeurat <- FindVariableFeatures(tempSeurat)
tempSeurat <-
  ScaleData(object = tempSeurat,
            vars.to.regress = 'nCount_RNA',
            verbose = TRUE)
tempSeurat <- RunPCA(object = tempSeurat, npcs = 100)

# Selecting the number of PC
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

# Clustering
tempSeurat <-
  FindNeighbors(tempSeurat, reduction = "pca", dims = 1:pc)
resolution <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 2)
tempSeurat <- FindClusters(tempSeurat, resolution = resolution)
tempSeurat <- RunUMAP(tempSeurat, dims = 1:pc)
