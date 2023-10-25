library(Seurat)

#dir: path of Cellranger output
matrix <- Read10X(data.dir = dir, gene.column = 1)

nUMI <- Matrix::colSums(matrix)
matrix = matrix[, head(order(nUMI, decreasing = T), 20000)]
tempSeurat <- CreateSeuratObject(matrix, project = Alias)

# Percentage of Mitochandrial Gene
mito.genes <- read.csv('MT.csv', stringsAsFactor = F)
mito.genes <- mito.genes[, 1]
mito.genes <- intersect(mito.genes, rownames(matrix))
percent.mito <-
  Matrix::colSums(matrix[mito.genes,]) / Matrix::colSums(matrix)
tempSeurat <-
  AddMetaData(object = tempSeurat,
              metadata = percent.mito,
              col.name = "percent.mito")

# Quality Control
bUMI = 200
pUMI = Inf
bGene = 100
pGene = Inf
pMT = 0.05
tempSeurat <-
  subset(
    tempSeurat,
    subset = nFeature_RNA > bGene &
      nFeature_RNA < pGene &
      percent.mito < pMT & nCount_RNA > bUMI & nCount_RNA < pUMI
  )

# QC by MAD
library(scater)
matrix <- GetAssayData(object = tempSeurat, slot = "counts")
sce <- SingleCellExperiment(list(counts = matrix))
mito.genes <- read.csv('MT.csv', stringsAsFactor = F)
mito.genes <- mito.genes[, 1]
is.mito <- intersect(mito.genes, rownames(matrix))
mito = rep(FALSE, nrow(matrix))
names(mito) <- rownames(matrix)
mito[is.mito] = TRUE
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = mito))
libsize.drop <-
  isOutlier(sce$total_counts,
            nmads = 3,
            type = "both",
            log = TRUE)
feature.drop <-
  isOutlier(
    sce$total_features_by_counts,
    nmads = 3,
    type = "both",
    log = TRUE
  )
sce <- sce[, !(libsize.drop | feature.drop)]
matrix = counts(sce)

tempSeurat <- CreateSeuratObject(matrix, project = Alias)
mito.genes <- read.csv('MT.csv', stringsAsFactor = F)
mito.genes <- mito.genes[, 1]
mito.genes <- intersect(mito.genes, rownames(matrix))
percent.mito <-
  Matrix::colSums(matrix[mito.genes,]) / Matrix::colSums(matrix)
tempSeurat <-
  AddMetaData(object = tempSeurat,
              metadata = percent.mito,
              col.name = "percent.mito")
