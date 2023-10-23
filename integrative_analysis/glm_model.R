#!Rscript
""" GLM model of DEG expression by HM vs nonHM DMRs.

    Author: Valentin Maurer <valentin.maurer@dkfz-heidelberg>
"""

library(bsseq)
library(Seurat)
library(data.table)
library(GenomicRanges)
library(ComplexHeatmap)

normRow = function(mat){
  t(apply(mat, MARGIN = 1, function(x){(x - min(x))/(max(x) - min(x))}))
}
reference = readRDS("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/external/haematopoiesis_10x/JMML-REF_Sampled.SeuratV3_CPL5_PC17.rds")
reference = NormalizeData(reference)
reference = reference[, reference@meta.data$predicted.id == "CD34+ HSC"]
reference@meta.data$PID = translate_pids(reference@meta.data$Patient)
reference@meta.data$PID[is.na(reference@meta.data$PID)] = reference@meta.data$Patient[is.na(reference@meta.data$PID)]
reference@meta.data$Epigenotype = factor(
  reference@meta.data$PID,
  levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "FCA", "HCA_BM", "HCA_CB", "JUV"),
  labels = c("LM", "LM", "IM", "IM", "HM", "HM", "HM", "HM", "FET", "ADU", "NEO", "JUV")
)
reference@meta.data$Epigenotype2 = factor(
  reference@meta.data$PID,
  levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "FCA", "HCA_BM", "HCA_CB", "JUV"),
  labels = c("nonHM", "nonHM", "nonHM", "nonHM", "HM", "HM", "HM", "HM", "FET", "ADU", "NEO", "JUV")
)
reference = reference[, !is.na(reference@meta.data$Epigenotype)]
Idents(reference) = "Epigenotype2"
level = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")
label = c("LM", "LM", "IM", "IM", "HM", "HM", "HM", "HM")
egt_level = c("LM", "IM", "HM")

dmrs = readRDS("/omics/groups/OE0219/internal/jmmlc_pbat/data/odcf_md/analysis/200612_DMR_model_sub_repMerged/dmrs_gr_sub_MethDiff.rds")[[1]]
# remove hem dmrs from HM vs non HM DMRs
dmrs_model = readRDS(
  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_TumorNormal/sig_dmrs_5inHalf_sub_anno.rds")
hem_dmrs = reduce(unlist(GRangesList(
  dmrs_model[["perinatal_vs_young_adult"]],
  dmrs_model[["perinatal_vs_juvenile"]],
  dmrs_model[["young_adult_vs_juvenile"]])))
keep = setdiff(1:length(dmrs), subjectHits(findOverlaps(dmrs, hem_dmrs)))
drop = intersect(1:length(dmrs), subjectHits(findOverlaps(dmrs, hem_dmrs)))
dmrs = dmrs[keep, ]

expressionInput <- c(`HM vs non-HM` = length(keep),
                     `Postnatal hematopoiesis` = length(hem_dmrs),
                     `HM vs non-HM&Postnatal hematopoiesis` = length(drop))
pdf("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/venn_HMnonHM_minusHem.pdf")
plot(euler(expressionInput), quantities = TRUE)
dev.off()

temp_new = tree_cols(bsdata)
temp_new[grepl(pattern = "P\\d", donor), class := ""]
temp_new[, class := as.vector(factor(donor,
                                     levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
                                     labels = c("LM", "LM", "IM", "IM", "HM", "HM", "HM", "HM")))
]
temp_new[donor == class, classColor := color]
temp_new[class == "LM", classColor := "#0058b4"]
temp_new[class == "IM", classColor := "#fbbb25"]
temp_new[class == "HM", classColor := "#c33126"]
temp_new[is.na(class), class := donor]
patient_methylation = bsseq::getMeth(collapseBSseq(bsdata, group = temp_new$donor),
                                     type = "raw", what = "perRegion", regions = dmrs)
subgroup_methylation = bsseq::getMeth(collapseBSseq(bsdata, group = temp_new$class),
                                      type = "raw", what = "perRegion", regions = dmrs)[,c("LM", "IM", "HM")]
dmrs_anno = as.data.frame(mcols(dmrs))[1:3]
mcols(dmrs) <- cbind(dmrs_anno, patient_methylation, subgroup_methylation)


markers = fread("/omics/odcf/analysis/OE0219_projects/JMMLC/scRNA_Result/Thesis/6-JMML_HSC-Epigenotype/JMML.HSC.HMnonHM.DEG.csv")
markers = markers[p_val_adj < 0.01]

gtf_gr = rtracklayer::import.gff("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/external/haematopoiesis_10x/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf")
gtf_gr = gtf_gr[gtf_gr$type == "gene"]


allGenesL = fread("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/external/haematopoiesis_10x/degLiftover_hg19.bed")
colnames(allGenesL) = colnames(allGenes)
setdiff(allGenes$gene_name, allGenesL$gene_name)

Idents(reference) = "PID"
avgExpression = AverageExpression(reference, features = allGenesL$gene_id)$RNA
avgExpression = avgExpression[, level]
Idents(reference) = "Epigenotype"
avgExpressionEGT = AverageExpression(reference, features = allGenesL$gene_id)$RNA
avgExpressionEGT = avgExpressionEGT[, egt_level]

REGION = 100 * 10**3
allGenes = allGenesL[gene_id %in% rownames(avgExpression)]
allGenesGR = makeGRangesFromDataFrame(allGenes, keep.extra.columns = T) + REGION

overlaps = findOverlaps(dmrs, allGenesGR)
subHits = subjectHits(overlaps)
queHits = queryHits(overlaps)
hits = unique(subHits)

models = lapply(hits, function(q){
  associatedDMRS = queHits[subHits == q]
  methData = as.data.table(dmrs[associatedDMRS])
  methData[, ID := paste0(seqnames, ":", start, "-", end)]
  dmrDistances = distanceToNearest(
    dmrs[associatedDMRS], allGenesGR[q,]-REGION)@elementMetadata$distance
  names(dmrDistances) = methData$ID

  methSubgroup = methData[, c("ID", "LM", "IM", "HM")]
  methSubgroup = as.data.frame(t(methSubgroup))
  colnames(methSubgroup) = methSubgroup["ID", ,drop=F]
  methSubgroup = methSubgroup[egt_level, ,drop=F]
  rnames = rownames(methSubgroup)
  methSubgroup = apply(methSubgroup, MARGIN = 2, as.numeric)
  methSubgroup = methSubgroup[, !colAnyNAs(methSubgroup), drop = F]
  rownames(methSubgroup) = rnames

  methData = methData[, c("ID", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")]
  methData = as.data.frame(t(methData))
  colnames(methData) = methData["ID", ,drop=F]
  methData = methData[level, ,drop=F]
  rnames = rownames(methData)
  methData = apply(methData, MARGIN = 2, as.numeric)
  methData = methData[, !colAnyNAs(methData), drop = F]
  rownames(methData) = rnames

  modelData = cbind(methData, as.data.frame(avgExpression[allGenes[q, gene_id], ]))
  colnames(modelData)[ncol(modelData)] = "avgExpression"

  modelGLM = glm(avgExpression ~., data = modelData)
  modelGLM = step(modelGLM)
  relDMRs = gsub(pattern = "`", replacement = "", names(modelGLM$coefficients))

  methReturn = methData[, intersect(colnames(methData), relDMRs), drop = F]
  methSubgroup = methSubgroup[, intersect(colnames(methData), relDMRs), drop = F]

  correlations = apply(methReturn, MARGIN = 2, function(x){
    cor(x, avgExpression[allGenes[q, gene_id], ], method = "spearman")})
  correlationsSubgroup = apply(methSubgroup, MARGIN = 2, function(x){
    cor(x, avgExpressionEGT[allGenes[q, gene_id], ], method = "spearman")})

  removedDMRs = setdiff(colnames(methData), colnames(methReturn))
  list(model               = modelGLM,
       methylation         = methReturn,
       methylationSubgroup = methSubgroup,
       expression          = as.data.frame(avgExpression[allGenes[q, gene_id], ]),
       expressionSubgroup  = as.data.frame(avgExpressionEGT[allGenes[q, gene_id], ]),
       correlation         = correlations,
       correlationSubgroup = correlationsSubgroup,
       distanceToGene      = dmrDistances[colnames(methReturn)],
       removedDist         = dmrDistances[removedDMRs],
       ResidNullDeviance   = with(summary(modelGLM), 1 - deviance/null.deviance))
})
names(models) = allGenes[hits, gene_name]

distances = as.vector(unlist(sapply(models, function(mod){mod$distanceToGene})))
distancesRM = as.vector(unlist(sapply(models, function(mod){mod$removedDist})))
wilcox.test(distances, distancesRM)
distances = data.table(
  V1 = c(distances, distancesRM), V2 = c(rep("kept", length(distances)), rep("removed", length(distancesRM)))
)
ggplot(distances, aes(x = V1, color = V2))+
  geom_density()+
  geom_histogram()+
  theme_bw(base_size = 14)+
  xlab("Distance to gene body [bp]")
summary(distances)


predictions = rbindlist(lapply(names(models), function(name){
  mod = models[[name]]
  prediction = stats::predict(mod$model)
  observation = t(mod$expression)
  observation = observation[,names(prediction)]
  correlation = cor(prediction, as.vector(unlist(mod$expression)), method = "pearson")
  rmse = sqrt(mean((prediction - observation)**2))
  mae = mean(abs(prediction - observation))
  rmae = mean(abs(rank(prediction) - rank(observation)))
  data.table(correlation = correlation, rmse = rmse, name = name, mae = mae, rmae = rmae)
}))
predictions[name == "CD52", label := name]
predictions = predictions[!is.na(correlation)]

cellSurfaceMarkers = fread(
  "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv")
cellSurfaceMarkers = cellSurfaceMarkers[grepl(pattern = "HUMAN", protein_name)]
cellSurfaceMarkers[, name := gsub(pattern = "(.*)_(.*)", replacement = "\\1", protein_name)]

top_features = data.table(ensembldb::select(ensdb_103, keys = predictions$name,
                                 keytype = "SYMBOL", columns = c("SYMBOL","UNIPROTID")))
colnames(top_features) = c("name", "uniprot")
top_features = top_features[complete.cases(top_features)]
top_features[, uniprot := gsub(uniprot, pattern = "(.*)\\..*", replacement = "\\1")]
top_features[, uniprot := gsub(uniprot, pattern = "(.*)-.*", replacement = "\\1")]
top_features = top_features[!duplicated(top_features)]
table(top_features$uniprot %in% unique(cellSurfaceMarkers$uniprot))
top_features = top_features[uniprot %in% unique(cellSurfaceMarkers[transmembrane == TRUE, uniprot])]

predictions[name %in% top_features$name, group := "cell-surface"]
predictions[name %in% c("LST1", "LTB", "IGLL1"), group := "cell-surface"]

predictions[is.na(group) & name %in% unique(cellSurfaceMarkers[transmembrane == TRUE, name]), group := "cell-surface"]
predictions[is.na(group), group := as.vector(factor(
  name, levels = gtf_gr$gene_name, labels = gtf_gr$gene_biotype))]
predictions[group == "protein_coding", group := "protein-coding"]
table(predictions$group)/nrow(predictions)

fisher.test(matrix(data = c(sum(temp_markers$group == "cell-surface"), sum(predictions$group == "cell-surface"),
                sum(temp_markers$group != "cell-surface"), sum(predictions$group != "cell-surface")),
       nrow = 2, ncol = 2))

fwrite(predictions, "/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/glmGaussian.csv", sep = ";")


translation_genes = data.table(ensembldb::select(ensdb_103, keys = gtf_gr$gene_name,
                                                 keytype = "SYMBOL", columns = c("SYMBOL","UNIPROTID")))
colnames(translation_genes) = c("name", "uniprot")
translation_genes = translation_genes[complete.cases(translation_genes)]
translation_genes[, uniprot := gsub(uniprot, pattern = "(.*)\\..*", replacement = "\\1")]
translation_genes[, uniprot := gsub(uniprot, pattern = "(.*)-.*", replacement = "\\1")]
translation_genes = translation_genes[!duplicated(translation_genes)]
table(translation_genes$uniprot %in% unique(cellSurfaceMarkers$uniprot))
translation_genes = translation_genes[uniprot %in% unique(cellSurfaceMarkers[transmembrane == TRUE, uniprot])]

temp_markers = markers
temp_markers$group = ""
temp_markers[gene %in% translation_genes$name, group := "cell-surface"]
temp_markers[gene %in% c("LST1", "LTB", "IGLL1"), group := "cell-surface"]
temp_markers[group == "" & gene %in% unique(cellSurfaceMarkers[transmembrane == TRUE, name]), group := "cell-surface"]
temp_markers[group == "", group := as.vector(factor(
  gene, levels = gtf_gr$gene_name, labels = gtf_gr$gene_biotype))]
temp_markers[group == "protein_coding", group := "protein-coding"]
table(temp_markers$group)/nrow(temp_markers)

translation_all = as.data.table(gtf_gr)
translation_all$group = ""
translation_all[gene_name %in% translation_genes$name, group := "cell-surface"]
translation_all[gene_name %in% c("LST1", "LTB", "IGLL1"), group := "cell-surface"]
translation_all[group == "" & gene_name %in% unique(cellSurfaceMarkers[transmembrane == TRUE, name]), group := "cell-surface"]
translation_all[group == "", group := as.vector(factor(
  gene_name, levels = gtf_gr$gene_name, labels = gtf_gr$gene_biotype))]
translation_all[group == "protein_coding", group := "protein-coding"]
table(translation_all$group)/nrow(translation_all)

total = data.table(Group = c("GLM", "DEG", "Total"),
                   CellSurface = c(sum(predictions$group == "cell-surface"), sum(temp_markers$group == "cell-surface"), sum(translation_all$group == "cell-surface")),
                   NotCellSurface = c(sum(predictions$group != "cell-surface"), sum(temp_markers$group != "cell-surface"), sum(translation_all$group != "cell-surface"))
)
total[, Total := CellSurface + NotCellSurface]
total[, Ratio := CellSurface/(CellSurface + NotCellSurface)]
total = total[order(Ratio)]
total[, Source := factor(Group, levels = unique(Group), labels = unique(Group))]
grid = expand.grid(unique(total$Group), unique(total$Group))

stat.test = data.table(group1 = c("Total", "Total", "DEG"),
                       group2 = c("DEG", "GLM", "GLM"),
                       p = c(
                         fisher.test(rbind(total[Group == "Total"][, c("CellSurface", "NotCellSurface")], total[Group == "DEG"][, c("CellSurface", "NotCellSurface")]))$p.value,
                         fisher.test(rbind(total[Group == "Total"][, c("CellSurface", "NotCellSurface")], total[Group == "GLM"][, c("CellSurface", "NotCellSurface")]))$p.value,
                         fisher.test(rbind(total[Group == "DEG"][, c("CellSurface", "NotCellSurface")], total[Group == "GLM"][, c("CellSurface", "NotCellSurface")]))$p.value
                       ))
stat.test[, p.adj := p.adjust(p, method = "BH")]
fwrite(total, "/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/proportion_cellsurface_comparison.csv")

ggplot(total, aes(x = Source, y = Ratio))+
  geom_col(fill = "#BABABA", col = "#000000")+
  ylab("Cell-surface markers in set [%]")+
  theme_bw(base_size = 14)


predictions = fread("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/glmGaussian.csv")
corCutoff = .9
predPlot = ggplot(predictions, aes(x = correlation, y = rmae, color = group, size = group))+
  geom_point()+
  geom_label_repel(predictions[correlation > corCutoff & group == "cell-surface"],
                   mapping = aes(x = correlation, y = rmae, label = name), show.legend = FALSE)+
  ylab("Ranked MAE") + xlab("Pearson correlation prediction vs observation")+
  theme_bw(base_size = 14)+
  scale_color_manual(name = "Type",
                     values = c(`cell-surface`   = "#a41220",
                                `protein-coding` = "#0058b4",
                                antisense        = "#A7D9CB",
                                lincRNA          ="#999999"))+
  scale_size_manual(name = "Type",
                    values = c(`cell-surface` = 3, `protein-coding` = 3,
                               antisense      = 1.5, lincRNA        =1.5))+
  geom_vline(xintercept = corCutoff, linetype = "dashed")+
  theme(legend.position = "bottom", aspect.ratio = 1)
predPlot
ggsave("/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/GLM_ModelEvaluation.pdf",
       width = 8, height = 8, predPlot)


# Per patient satisfying cutoff
candidates = predictions[correlation > corCutoff]$name
methTotal  = rbindlist(lapply(candidates, function(name){
  mod = models[[name]]
  if(ncol(mod$methylation) == 0){
    return(data.table())
  }
  return(data.table(t(mod$methylation), keep.rownames = T))}))
corTotal = as.vector(unlist(sapply(candidates, function(name){
  mod = models[[name]]
  return(mod[["correlation"]])})))
geneTotal = rbindlist(lapply(candidates, function(name){
  mod = models[[name]]
  if(ncol(mod$methylation) == 0){
    return(data.table())
  }
  rbindlist(lapply(1:ncol(mod$methylation), function(i){
    res = data.table(t(mod$expression))
    res$rn = name
    res}))}))

methHtmp = Heatmap(
  data.frame(methTotal[,-"rn"], row.names = make.names(methTotal$rn, unique = T)),
  clustering_distance_columns = cluster_na,
  clustering_distance_rows = cluster_na,
  col = ANNOTATION_COLORSCALE_METH,
  name = "Average DMR methylation",
  row_split = corTotal > 0,
  top_annotation = HeatmapAnnotation(
   border = T,
   Patient = colnames(methTotal[,-"rn"]),
   col = list(
     Patient = COLORSCHEME[colnames(methTotal[,-"rn"])]
   )),
  column_split = label,
  left_annotation = rowAnnotation(
   Correlation = corTotal,
   CellSurfaceMarker = ifelse(geneTotal$rn %in% predictions[group == "cell-surface"]$name, "Yes", "No"),
   Type = factor(geneTotal$rn, levels = predictions$name, labels = predictions$group),
   col = list(Correlation = colorRamp2(c(-1, 0, 1),c("#2166AC", "#FFFFFF", "#B2182B")),
              CellSurfaceMarker = c(Yes = "#000000", No = "#BABABA"),
              Type = c(`cell-surface`   = "#a41220", `protein-coding` = "#0058b4",
                       antisense        = "#A7D9CB", lincRNA          ="#999999"))
  ),
  border = "#000000", show_column_names = F,
  cluster_columns = F,
  cluster_row_slices = F, cluster_column_slices = F,
  show_row_names = F
)
rnaHtmp = Heatmap(
  normRow(data.frame(geneTotal[,-"rn"], row.names = make.names(geneTotal$rn, unique = T))),
  clustering_distance_columns = cluster_na,
  clustering_distance_rows = cluster_na,
  col = colorRampPalette(c("#FFFFFF", "#CA3C1B"))(100),
  name = "Normalized gene expression",
  top_annotation = HeatmapAnnotation(
    border = T,
    Patient = colnames(geneTotal[,-"rn"]),
    col = list(
      Patient = COLORSCHEME[colnames(geneTotal[,-"rn"])]
    )),
  right_annotation = rowAnnotation(
    foo = anno_mark(at = which(geneTotal$rn %in% top_features$name),
    labels = geneTotal$rn[geneTotal$rn %in% top_features$name])),
  column_split = label,
  border = "#000000", show_column_names = F,
  cluster_columns = F, show_row_names = F,
  cluster_row_slices = F, cluster_column_slices = F
)
pdf(file.path(
  "/omics/odcf/analysis/OE0565_projects/tce/jmml_coo/JMMLT/analysis/meth_rna/initPaper/",
  paste0("GLM_ExpressionMethylation_Patient_Candidate_perCorrelation_#CA3C1B.pdf")),
    width = 14, height = 30)
draw(methHtmp + rnaHtmp)
dev.off()
