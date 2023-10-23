#!Rscript
""" Define lineage-specific pseudotime through principle component analysis.

    Author: Valentin Maurer <valentin.maurer@dkfz-heidelberg>
"""

library(Seurat)
library(data.table)
library(ggplot2)
library(parallel)
library(ordinal)

BASEPATH = "/omics/groups/OE0219/internal/Valentin/JMMLT/analysis/meth_rna/initPaper/lineages"
#SEURAT_OBJECT = "/omics/odcf/analysis/OE0219_projects/JMMLC/scRNA_Data/4.All_Dataset/JMML_MNC-Postnatal_REF_Sampled/StandardAnalysis/JMML_MNC-Postnatal_REF_Sampled.Seurat.Standard.rds"
SEURAT_OBJECT = "/omics/groups/OE0219/internal/Valentin/JMMLT/external/haematopoiesis_10x/JMML_MNC-Postnatal_REF_Sampled-Anchor.Seurat.PC1-18.rds"
# 1:num_pcs are tested for association with metadata
num_pcs = 15
REL_FEATURES = c("outcome", "Type", "Genotype", "Subgroup", "PID", "Phase", "nCount_RNA", "nFeature_RNA")
# LINEAGE needs to be in order
LINEAGES = list(
  T = c("CD34+ HSC"="#252525", "CD34+ LMPP"="#99b9f8", "CD34+ CLP"="#7dacfc",
        "CD34+ pre-T"="#9a8af8", "Naive T-cell"="#6f5df6", "CD8 T-cell"="#403cd2"),
  NK = c("CD34+ HSC"="#252525", "CD34+ LMPP"="#99b9f8", "CD34+ CLP"="#7dacfc",
         "NK cells"="#0875f7"),
  B = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
        "CD34+ LMPP"="#99b9f8", "CD34+ CLP"="#7dacfc", "CD34+ pre-B cycling"="#5a9fff",
        "CD34+ Lymphoid UNK"="#0d82fa", "CD34+ pre-B"="#0079dc", "Pro-B"="#005faf",
        "Follicular B cell"="#254887"),
  Ery = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
          "CD34+ MEP"="#c2938a", "CD34+ ERP-Early"="#ca7c6f", "CD34+ ERP"="#de7160",
          "Early-Erythroblast"="#e05646", "Erythroblast"="#e62628"),
  Mega = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
           "CD34+ MEP"="#c2938a", "CD34+ MKP"="#d9d9d9", "Platelet"="#d9d9d9"),
  DC = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
         "CD34+ MDP-2"="#bbce7a", "CD34+ MDP-1"="#bbce5c", "Pre-Dendritic"="#94c529",
         "Dendritic Cell"="#5ca800"),
  Mono = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
           "CD34+ MDP-2"="#bbce7a", "CD34+ MDP-1"="#bbce5c",
           "CD34+ Gran"="#e8c091", "Granulocytic-UNK"="#e1b072", "Monocyte"="#f6be13"),
  Neutro = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
            "CD34+ Gran"="#e8c091", "Granulocytic-UNK"="#e1b072",
            "Immature-Neutrophil"="#fea361", "Neutrophil"="#f57e12"),
  Eo = c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
         "CD34+ Eo/B/Mast"="#9356a9", "Eosinophil"="#a643d2")
)

# LINEAGE needs to be in order
# LINEAGE = c("CD34+ HSC", "CD34+ HSC-cycle", "CD34+ MultiLin",
#             "CD34+ MDP-2", "CD34+ MDP-1", "Pre-Dendritic", "Dendritic Cell")
# Colors, do not need to be in order
# DC_cols <- c("CD34+ HSC"="#252525", "CD34+ HSC-cycle"="#737373", "CD34+ MultiLin"="#b3b3b3",
#              "CD34+ MDP-2"="#bbce7a", "CD34+ MDP-1"="#bbce5c", "Pre-Dendritic"="#94c529", "Dendritic Cell"="#5ca800")

data_total = readRDS(SEURAT_OBJECT)
data_total@meta.data$Dataset = ifelse(data_total@meta.data$Project == "JMML",
                                yes = "Tumor", no = "Normal")
data_total@meta.data$Genotype[is.na(data_total@meta.data$Genotype)] = "WT"
data_total@meta.data$Epigenotype[is.na(data_total@meta.data$Epigenotype)] = "WT"
DefaultAssay(data_total) = "RNA"
Idents(data_total) = "predicted.id"
for(k in names(LINEAGES)){
  print(k)
  LINEAGE_COLS = LINEAGES[[k]]
  LINEAGE = names(LINEAGE_COLS)
  LINEAGE_PATH = file.path(BASEPATH, k)
  dir.create(LINEAGE_PATH, recursive = T, showWarnings = F)

  data <- subset(data_total, idents = LINEAGE)
  data = NormalizeData(data)
  data = data[rowSums(data@assays$RNA@counts) != 0, ]
  data@meta.data$outcome = as.numeric(
    factor(data@meta.data$predicted.id, levels = LINEAGE,
           labels = 1:length(LINEAGE))
  )

  is_healthy = data@meta.data$Project != "JMML"
  data = ScaleData(data, feature = rownames(data))
  assay = data@assays$RNA@data[, is_healthy]
  project = data@meta.data[is_healthy, "Project"]
  celltype =   factor(data@meta.data$predicted.id, levels = LINEAGE,
                      labels = 1:length(LINEAGE))[is_healthy]

  # Fit model predicting celltype based on expression
  models = mclapply(1:nrow(assay), function(i){
    # model = lm(celltype ~ assay[i,])
    model = clm(celltype ~ assay[i, ])
    model_summary = summary(model)
    # p_val = ifelse(dim(model_summary$coefficients)[1] == 1,
                   # NA, no = model_summary$coefficients[, 4][[2]])
    # coefficient = model$coefficients[[2]]
    p_val = model_summary$coefficients[,4][nrow(model_summary$coefficients)][[1]]
    coefficient = model_summary$coefficients[,1][nrow(model_summary$coefficients)][[1]]
    data.table(gene = rownames(assay)[i],
               coefficient = coefficient,
               p_value = p_val
               # correlation = cor(outcome, assay[i,], method = "spearman")
               )
  }, mc.cores = 12)
  models = rbindlist(models)
  models[, p_adjust := p.adjust(p_value, method = "BH")]

  # Select features based on coefficient p-value and coefficient magnitude
  features = models[p_adjust < 0.01 & abs(coefficient) > 1.5, gene]
  # Run PCA on selected features
  temp = ScaleData(data, features = features)
  temp = RunPCA(temp, features = features, verbose = F)

  # Search for association between metadata and 1:num_pcs PCs
  drop = grepl(colnames(temp@meta.data), pattern = "logreg", ignore.case = T)
  drop = drop | grepl(colnames(temp@meta.data), pattern = "lineage", ignore.case = T)
  drop = drop | grepl(colnames(temp@meta.data), pattern = "RNA_snn", ignore.case = T)
  drop = drop | grepl(colnames(temp@meta.data), pattern = "prediction.score", ignore.case = T)
  drop = drop | grepl(colnames(temp@meta.data), pattern = "integrated_snn", ignore.case = T)
  metaAssociation = rbindlist(lapply(colnames(temp@meta.data)[!drop], function(metaf){
    x = temp@meta.data[, metaf]
    mask = !is.na(x)
    x = x[mask]
    if(length(x) == 0){
      return(data.table())
    }
    x_numeric = x
    if(!is.numeric(x)){
      x_numeric = as.numeric(factor(x, levels = unique(x), labels = 1:length(unique(x))))
    }
    rbindlist(lapply(colnames(temp@reductions$pca@cell.embeddings)[1:num_pcs],
                     function(pc){
      y = temp@reductions$pca@cell.embeddings[mask, pc]
      p_val = ifelse(length(unique(x)) == 1,
                     yes = NA, no = kruskal.test(y, x)$p.value)
      data.table("p_val" = p_val,
                 "cor"   = cor(x_numeric, y, method = "spearman"),
                 Feature = metaf, PC = pc)
      }))
  }))

  # Find PCs corresponding to celltype (outcome) and Dataset (tumor/normal)
  metaAssociation[Feature == "Dataset"][order(abs(cor), decreasing = T)]
  metaAssociation[Feature == "outcome"][order(abs(cor), decreasing = T)]
  fwrite(metaAssociation, file.path(LINEAGE_PATH, "pcaAssociation.csv"),
         sep = ";")
  pc_bubble = ggplot(metaAssociation[Feature %in% REL_FEATURES],
                     aes(x=PC, y=Feature, size = abs(cor))) +
    geom_point(alpha=1, shape=21, color="black", fill="black") +
    ylab(element_blank()) +
    xlab(element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = rel(2), color = "black"),
          axis.text.x = element_text(angle = 90),
          axis.title = element_text(size = rel(2), color = "black"),
          legend.text = element_text(size = rel(1.5), color = "black"),
          legend.title = element_text(size = rel(2), color = "black")) +
    scale_radius(name="Cor.", range = c(0, 10))
  ggsave(file.path(LINEAGE_PATH, "pca_bubbleplot.pdf"), pc_bubble,
         width = 8, height = 14)

  # Select PC with highest absolute correlation with outcome (celltype) for plotting
  celltype_pc = metaAssociation[Feature == "outcome"][order(abs(cor), decreasing = T), PC][1]
  disease_pc = metaAssociation[Feature == "Type"][order(abs(cor), decreasing = T), PC][1]
  top_disease = metaAssociation[Feature == "Type"][order(abs(cor), decreasing = T), PC][1:min(5, num_pcs)]
  dims = c(celltype_pc, disease_pc)
  dims = as.numeric(gsub(dims, pattern = "PC_(\\d)", replacement = "\\1"))

  pca_dataset = DimPlot(temp, reduction = "pca", group.by = c("Dataset"),
                        dims = dims, pt.size = 2, raster = F,
                        order = c("Tumor", "Normal"))
  pca_datasetSplit = DimPlot(temp, reduction = "pca", split.by = c("Dataset"),
                        dims = dims, pt.size = 2, raster = F,
                        group.by = "Dataset")
  pca_ids = DimPlot(temp, reduction = "pca", group.by = c("predicted.id"),
                    dims = dims, pt.size = 2, raster = F)+
    scale_color_manual(values = LINEAGE_COLS)
  ggsave(file.path(LINEAGE_PATH, "pcaDataset.pdf"), pca_dataset, width = 8, height = 6)
  ggsave(file.path(LINEAGE_PATH, "pcaDatasetSplit.pdf"), pca_datasetSplit, width = 16, height = 6)
  ggsave(file.path(LINEAGE_PATH, "pcaCelltype.pdf"), pca_ids, width = 8, height = 6)


  ##### Figures per celltype
  pc_data = data.table(pc = temp@reductions$pca@cell.embeddings[, celltype_pc],
                       color = temp@meta.data$predicted.id,
                       project = temp@meta.data$Project)
  # Flip PC so that beginning of lineage is on the left and end of the lineage on the right hand side
  flip = ifelse(
    mean(pc_data[color == LINEAGE[1]]$pc) > mean(pc_data[color == LINEAGE[length(LINEAGE)]]$pc),
    yes = -1, no = 1)
  pc_data[, pc := pc * flip]
  pc_dens = ggplot(pc_data, aes(x = pc, fill = color))+
    geom_density(alpha = .5)+
    scale_fill_manual(values = LINEAGE_COLS)+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab(celltype_pc)
  ggsave(file.path(LINEAGE_PATH, "pcaCelltype_density.pdf"), pc_dens, width = 10, height = 6)

  pc_data = data.table(pc = temp@reductions$pca@cell.embeddings[, celltype_pc],
                       color = temp@meta.data$predicted.id,
                       project = temp@meta.data$Project)
  pc_data = pc_data[abs(pc) < 20]
  pc_data[, pc := pc * flip]
  pc_data[, breaks := cut(pc, breaks = 30)]
  pc_data = pc_data[, .N, by = .(color, breaks)]
  pc_data[ , proportion := N/sum(N), by = .(breaks)][, N := NULL]
  pattern <- "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"
  pc_data[, xstart := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xend := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xlabel := mean(c(xstart, xend))]
  pc_data[, color := factor(color, levels = LINEAGE, labels = LINEAGE)]
  pc_prop = ggplot(pc_data, aes(x = xstart, y = proportion, fill = color))+
    geom_bar(stat = "identity", col = "#000000")+
    scale_fill_manual(values = LINEAGE_COLS)+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab(celltype_pc)+
    ylab("Proportion")
  ggsave(file.path(LINEAGE_PATH, "pcaCelltype_proportion.pdf"), pc_prop,
         width = 10, height = 6)

  # Same thing but with disease pc
  pc_data = data.table(pc = temp@reductions$pca@cell.embeddings[, disease_pc],
                       color = temp@meta.data$predicted.id,
                       project = temp@meta.data$Project)
  flip = ifelse(
    mean(pc_data[color == LINEAGE[1]]$pc) > mean(pc_data[color == LINEAGE[length(LINEAGE)]]$pc),
    yes = -1, no = 1)
  pc_data[, pc := pc * flip]
  pc_dens = ggplot(pc_data, aes(x = pc, fill = color))+
    geom_density(alpha = .5)+
    scale_fill_manual(values = LINEAGE_COLS)+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab(disease_pc)
  ggsave(file.path(LINEAGE_PATH, "pcaCelltype_densityDisease.pdf"), pc_dens, width = 10, height = 6)

  pc_data = pc_data[abs(pc) < 20]
  pc_data[, breaks := cut(pc, breaks = 30)]
  pc_data = pc_data[, .N, by = .(color, breaks)]
  pc_data[ , proportion := N/sum(N), by = .(breaks)][, N := NULL]
  pattern <- "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"
  pc_data[, xstart := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xend := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xlabel := mean(c(xstart, xend))]
  pc_data[, color := factor(color, levels = LINEAGE, labels = LINEAGE)]
  pc_prop = ggplot(pc_data, aes(x = xstart, y = proportion, fill = color))+
    geom_bar(stat = "identity", col = "#000000")+
    scale_fill_manual(values = LINEAGE_COLS)+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab(disease_pc)+
    ylab("Proportion")
  ggsave(file.path(LINEAGE_PATH, "pcaCelltype_proportionDisease.pdf"), pc_prop,
         width = 10, height = 6)

  ##### Figures per type
  pc_data = cbind(temp@reductions$pca@cell.embeddings[, top_disease],
                  data.table(color = temp@meta.data$Type))
  pc_data = data.table::melt(pc_data, id.vars = "color")

  pc_dens = ggplot(pc_data, aes(x = value, fill = color))+
    geom_density(alpha = .5)+
    scale_fill_manual(name = "Type", values = list(Normal = "#BABABA", Tumor = "#c33126"))+
    facet_wrap(~variable, scales = "free")+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab(celltype_pc)
  ggsave(file.path(LINEAGE_PATH, "pcaType_density.pdf"), pc_dens, width = 12, height = 10)

  pc_data = pc_data[abs(value) < 20]
  pc_data[, breaks := cut(value, breaks = 30)]
  pc_data = pc_data[, .N, by = .(color, breaks, variable)]
  pc_data[ , proportion := N/sum(N), by = .(breaks, variable)][, N := NULL]
  pattern <- "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"
  pc_data[, xstart := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xend := as.numeric(gsub(pattern,"\\2", breaks))]
  pc_data[, xlabel := mean(c(xstart, xend))]
  pc_prop = ggplot(pc_data, aes(x = xstart, y = proportion, fill = color))+
    geom_bar(stat = "identity", col = "#000000")+
    facet_wrap(.~variable, scales = "free_x")+
    scale_fill_manual(name = "Type", values = list(Normal = "#BABABA", Tumor = "#c33126"))+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom")+
    xlab("PC")+
    ylab("Proportion")
  ggsave(file.path(LINEAGE_PATH, "pcaType_proportion.pdf"), pc_prop, width = 12, height = 10)

  saveRDS(temp, file.path(LINEAGE_PATH, "seurat_object.rds"))
}
