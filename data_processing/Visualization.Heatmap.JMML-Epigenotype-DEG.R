library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd(
        '/icgc/dkfzlsdf/analysis/OE0219_projects/JMMLC/scRNA_Result/Thesis/6-JMML_HSC-Epigenotype/'
)

df.sc <- readRDS('JMML.Methy.df.sc.rds')
gene.ann <- readRDS('./JMML.Methy.gene.ann.rds')
jmml.ann <- readRDS('./JMML.Methy.jmml.ann.rds')
jmml.color <- readRDS('./JMML.Methy.jmml.color.rds')
mg <- read.table('Mark.HMnonHM.txt')
GeneIDMap <-
        read.table('genes.tsv', header = F, stringsAsFactors = F)


col_fun <- colorRamp2(c(-3, 0, 3), c('#2166ac', '#ffffff', '#b2182b'))
genecolor <- c('#b2182b', '#2166ac')
names(genecolor) <- levels(gene.ann$Gene)
genecolor <- list(Gene = genecolor)

rownames(df.sc) <-
        GeneIDMap[match(rownames(df.sc), GeneIDMap$V1), 'V2']
ind <- match(mg$V1, rownames(df.sc))
jmml.ann <- jmml.ann[colnames(df.sc), ]

ha = rowAnnotation(foo = anno_mark(at = ind, labels = as.character(mg$V1)))
Heatmap(
        df.sc,
        col = col_fun,
        right_annotation = ha,
        left_annotation = rowAnnotation(df = gene.ann, col = genecolor),
        top_annotation = HeatmapAnnotation(df = jmml.ann,
                                           col = jmml.color),
        show_column_names = T,
        show_row_names = F,
        cluster_rows = F,
        cluster_columns = F
)
