#!Rscript
""" Relation of JMML to normal HSC development

    Author: Valentin Maurer <valentin.maurer@dkfz-heidelberg>
"""

required_libs = c("MOFA2", "data.table", "parallel", "ggplot2", "ComplexHeatmap",
  "bsseq", "ape", "RColorBrewer")
required_libs = sapply(required_libs, function(x){
  suppressPackageStartupMessages(library(x, character.only = T))}
)
rm(required_libs)
source("/omics/groups/OE0219/internal/Valentin/JMMLT/scripts/accessory_funcs.R")

tree_cols = function(data){
  pheno = readRDS(
    "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq/stats/bsseq_HSC_comb_snpRemoved_repMerged_pheno.rds"
  )
  cluster_cols = data.table(sample = colnames(data), donor = NA)
  cluster_cols[, donor := factor(sample, levels = rownames(pheno), labels = pheno$Donor)]
  cluster_cols[is.na(donor), donor := sample]
  cluster_cols[, Patient := translate_pids(donor)]
  cluster_cols[!is.na(Patient), donor := Patient][, Patient := NULL]
  cluster_cols[, donor := as.character(donor)]
  cluster_cols[!grepl(pattern = "P\\d", donor) & ! donor %in% c("LML", "HML"), donor := factor(
    donor, levels = c("cordblood", "adult_bonemarrow"), labels = c("NEO", "ADU"))]
  cluster_cols[grepl(pattern = "JU", sample), donor := "JUV"]
  cluster_cols[grepl(pattern = "FS", sample), donor := "FES"]
  cluster_cols[grepl(pattern = "FL", sample), donor := "FEL"]
  cluster_cols[, color := COLORSCHEME[donor]]
  cluster_cols[, donorCtype := donor]
  cluster_cols[donor == "ADU", donorCtype := paste(
    donor, gsub(sample, pattern = "(.*)_(.*)_(.*)\\d_NORMAL", replacement = "\\3"),
    sep = "_")]
  cluster_cols[sample == donor, name := donorCtype]
  cluster_cols[sample != donor, name := paste0(donorCtype, ".", 1:.N), by = donorCtype]
  cluster_cols[sample %in% c("P2", "P3"), donorCtype := paste(donorCtype, "scIGMT", sep = "_")]
  return(cluster_cols)
}
tree_cols2 = function(data_colnames){
  pheno = readRDS(
    "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq/stats/bsseq_HSC_comb_snpRemoved_repMerged_pheno.rds"
  )
  cluster_cols = data.table(sample = data_colnames, donor = NA)
  cluster_cols[, donor := factor(sample, levels = rownames(pheno), labels = pheno$Donor)]
  cluster_cols[is.na(donor), donor := sample]
  cluster_cols[, Patient := translate_pids(donor)]
  cluster_cols[!is.na(Patient), donor := Patient][, Patient := NULL]
  cluster_cols[, donor := as.character(donor)]
  cluster_cols[!grepl(pattern = "P\\d", donor) & ! donor %in% c("LML", "HML"), donor := factor(
    donor, levels = c("cordblood", "adult_bonemarrow"), labels = c("NEO", "ADU"))]
  cluster_cols[grepl(pattern = "JU", sample), donor := "JUV"]
  cluster_cols[grepl(pattern = "FS", sample), donor := "FES"]
  cluster_cols[grepl(pattern = "FL", sample), donor := "FEL"]
  cluster_cols[, color := COLORSCHEME[donor]]
  cluster_cols[, donorCtype := donor]
  cluster_cols[donor == "ADU", donorCtype := paste(
    donor, gsub(sample, pattern = "(.*)_(.*)_(.*)\\d_NORMAL", replacement = "\\3"),
    sep = "_")]
  cluster_cols[sample == donor, name := donorCtype]
  cluster_cols[sample != donor, name := paste0(donorCtype, ".", 1:.N), by = donorCtype]
  cluster_cols[sample %in% c("P2", "P3"), donorCtype := paste(donorCtype, "scIGMT", sep = "_")]
  return(cluster_cols)
}


COLORSCHEME <- c(
  ADU = "#E5E5E5", JUV = "#A7D9CB", NEO = "#978474", FEL = "#6D7A9F", FES = "#595959",
  P1 = "#0058b4",  P2 = "#2188c9",  P3 = "#fbbb25",  P4 = "#fca349",
  P5 = "#ff6b36",  P6 = "#e34e2e",  P7 = "#c33126",  P8 = "#a41220",
  LML = "#23395d", HML = "#6f0000", FET = "#BABABA"
)

###############################################################################
BIN = "1kb"; NREGIONS = 10* 10 ** 3; NREGIONS = 5000
regions  = fread(paste0("/omics/groups/OE0219/internal/genomes/Hsapiens/GRCh37/seq/GRCh37.", BIN, ".bed"))
regions  = regions[V1 %in% 1:22][, V1 := paste0("chr", V1)]
regions  = makeGRangesFromDataFrame(regions,
                                    seqnames.field = "V1", start.field = "V2", end.field = "V3")


# Select regions that correlate with developmental state
bsdata <- readRDS("/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/bsseq/bsseq_HSC_combTotal_snpRemoved.rds")
pheno  <- colData(bsdata)
drop   <- pheno$Name %in% c("FL2_HSC_PBAT_2", "JU3_HSC_PBAT_11")
bsdata <- bsdata[, !drop]
temp_anno = tree_cols2(colnames(bsdata))
keep = !temp_anno$donorCtype %in% c("ADU_CMP", "ADU_GMP", "ADU_MEP", "ADU_MPP")
bsdata      <- bsdata[, keep]


temp_meth = bsseq::getMeth(bsdata, type = "raw", regions = regions, what = "perRegion")
temp_new = tree_cols(temp_meth)
keep_regions = complete.cases(temp_meth) &
  rowSds(temp_meth[,!grepl(pattern = "P\\d", temp_new$donor)],na.rm = T) > 0.01
temp_meth = temp_meth[keep_regions, ]

group = temp_new$donorCtype
group[group %in% c("FEL", "FES")] = "FET"
group[group == "ADU_HSC"] = "ADU"
state_meth = temp_meth[, group %in% c("FET", "NEO", "JUV", "ADU")]
groups = tree_cols(state_meth)$donor
groups[groups %in% c("FEL", "FES")] = "FET"
timings = as.numeric(factor(groups,
                            levels = c("FET", "NEO", "JUV", "ADU"), labels = c(1, 2, 3, 4)))
state_cor = cor(t(state_meth), timings, method = "spearman")
# pearson produces similar results but with only ~400 regions
# state_cor = cor(t(state_meth), timings, method = "pearson")
keep_regions_cor = abs(state_cor) >= .9
keep_regions_cor[is.na(keep_regions_cor)] = FALSE

# All samples
not_progenitor = tree_cols(temp_meth)
not_progenitor = !not_progenitor$donorCtype %in% c("ADU_CMP", "ADU_GMP", "ADU_MEP", "ADU_MPP")
temp = temp_meth[keep_regions_cor, not_progenitor]

# Check if the selected sites correlate with TWGBS or WGBS methodolgy (they dont)
seqMethod = grepl(pattern = "C010", colnames(temp)) * 1
method_cor = cor(t(temp), seqMethod, method = "spearman")
table(abs(method_cor) >.5)
keep = abs(method_cor) < .4

# Plot normal tree with method correction
phyl_dist = as.matrix(stats::dist(t(state_meth[keep_regions_cor, ][keep, ]), method = "man"))
temp_anno = tree_cols(phyl_dist)
colnames(phyl_dist) = temp_anno$name; rownames(phyl_dist) = temp_anno$name
temp_anno = data.frame(temp_anno[, c("color"), drop = F], row.names = temp_anno$name)
plot_tree(tree = ape::nj(X = as.matrix(phyl_dist)), col_df = temp_anno, col_id = 1,
          type = 'unrooted', rotate = 240)

temp_cols = tree_cols(state_meth)
state = temp_cols$donor
state[state %in% c("FEL", "FES")] = "FET"
state = factor(state, levels = c("FET", "NEO", "JUV", "ADU"), labels = c("FET", "NEO", "JUV", "ADU"))
state_order = order(state)
state = state[state_order]
temp_cols = temp_cols[state_order, ]

region_anno = rowMeans(state_meth[
  keep_regions_cor, c("FL1_HSC_PBAT_1", "FL3_HSC_PBAT_3", "FL4_HSC_PBAT_4", "FS1_HSC_PBAT_5",
                      "FS2_HSC_PBAT_6", "FS3_HSC_PBAT_7")]) >  rowMeans(state_meth[
                        keep_regions_cor, c("C010_HSC_HSC1_NORMAL", "C010_HSC_HSC2_NORMAL",
                                          "C010_HSC_HSC3_NORMAL")])
region_anno[region_anno == FALSE] = "Hyper"
region_anno[region_anno == TRUE] = "Hypo"

htmp = Heatmap(state_meth[keep_regions_cor, state_order][keep, ],
               top_annotation = HeatmapAnnotation(
                 border = T,
                 HSC = factor(temp_cols$donor, levels = c("FEL", "FES", "NEO", "JUV", "ADU"),
                                 labels = c("FEL", "FES", "NEO", "JUV", "ADU")),
                 col = list(
                   HSC = COLORSCHEME[c("FEL", "FES", "NEO", "JUV", "ADU")])
               ),
               col = ANNOTATION_COLORSCALE_W, show_column_names = F,
               border = "#000000",
               heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
               name = "Average methylation",
               column_split = state, cluster_column_slices = F,
               row_split = region_anno[keep]
)
draw(htmp, heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")


# Merged samples
sample_meth = bsseq::getMeth(collapseBSseq(bsdata, group = temp_new$donorCtype),
                             type = "raw", what = "perRegion",
                             regions = regions[as.vector(keep_regions), ][as.vector(keep_regions_cor), ])

phyl_dist = as.matrix(stats::dist(t(sample_meth[keep, ]), method = "man"))
temp_cols = temp_new[, c("color", "donorCtype")]
temp_cols = temp_cols[!duplicated(temp_cols), ]
temp_cols = data.frame(temp_cols[, c("color"), drop = F], row.names = temp_cols$donorCtype)
temp_cols = temp_cols[rownames(phyl_dist), ,drop=F]
plot_tree(tree = fastme.bal(X = as.matrix(phyl_dist)),
          col_df = temp_cols, col_id = 1, type = "unrooted", rotate = 120)

# Check region methylation in scIGMT
final_regions = as.data.table(
  regions[as.vector(keep_regions), ][as.vector(keep_regions_cor), ]
)[,1:3]
final_regions[, index := .I]
region_meth = cbind(final_regions, sample_meth)
saveRDS(final_regions, "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/DevelopmentalRegions.RDS")
saveRDS(region_meth, "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/DevelopmentalRegions_Meth.RDS")
# for ucsc hg19 to hg38 liftover
# fwrite(final_regions, "~/regions.tsv", sep = "\t", col.names = F)

bed_files = list.files(
  "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq/filteredMethCallsPE/",
  full.names = T, pattern = "*bed"
)
bed_files = bed_files[!grepl(pattern = "i7-only", ignore.case = T, bed_files)]
devRegions = fread("/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/DevelopmentalRegions_hg38.bed")
colnames(devRegions) = c("Chr", "Start", "End", "Index")
devRegions[, Chr := gsub(pattern = "chr", replacement = "", Chr)]
setkey(devRegions, Chr, Start, End)
avg_region_cgmeth = rbindlist(mclapply(bed_files, function(i) {
  data = fread(i)
  setkey(data, chr, start, end)
  data = foverlaps(data, devRegions,
                   type = "within",
                   by.x = c("chr", "start", "end"),
                   by.y = c("Chr", "Start", "End"),
                   nomatch = 0L
                   )
  data = data[, .(n_total = sum(n_total, na.rm = T),
                  n_meth = sum(n_meth, na.rm = T),
                  avg_beta = mean(beta_value, na.rm = T)),
              by = .(Index)]
  data[, sample := paste(translate_bsnames(i)[c("PID", "PLATE", "WELL")], collapse = "_")]
}, mc.cores = 10))
saveRDS(avg_region_cgmeth, "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/DevelopmentalRegion_methylation_scigmt.RDS")

avg_dev_cgmeth = readRDS("/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/DevelopmentalRegion_methylation_scigmt.RDS")
model = load_model("/omics/groups/OE0219/internal/Valentin/JMMLT/processing/meth_rna/mofa/enhancermeth-facs-rna-tssmeth.hdf5")
model@samples_metadata$pid = gsub(
  pattern = "(.*)_(.*)_(.*)", replacement = "\\1", model@samples_metadata$sample)
model@samples_metadata$plate = gsub(
  pattern = "(.*)_(.*)_(.*)", replacement = "\\2", model@samples_metadata$sample)
set.seed(42)
model = run_tsne(model)
x_dimr = as.data.table(plot_dimred(model, method = "TSNE", color_by = "pid",
                                          dot_size = 3, return_data = T))
x_dimr[(x < -8) & (x > -18) & (y > -1) & (y < 7.5), group := "HML"]
x_dimr[(x < 0) & (x > -12) & (y > -12) & (y < -4) , group := "HML"]
x_dimr[(x < 24) & (x > 15) & (y > -5) & (y < 0),    group := "LML"]
x_dimr[is.na(group) & color_by %in% c("P6", "P7"),  group := "HML"]
x_dimr[is.na(group) & color_by == "P2",  group := "LM"]
x_dimr[is.na(group) & color_by == "P3",  group := "IM"]

avg_dev_cgmeth[, group := factor(sample, levels = x_dimr$sample, labels = x_dimr$group)]
avg_dev_cgmeth = avg_dev_cgmeth[!is.na(group)]
avg_dev_cgmeth = avg_dev_cgmeth[, .(avg_beta = mean(avg_beta)), by = .(Index, group)]
avg_dev_cgmeth = dcast(avg_dev_cgmeth, Index ~ group, value.var = "avg_beta")
colnames(avg_dev_cgmeth)[1] = c("index")


total_dev = merge(region_meth[as.vector(unlist(keep)), ], avg_dev_cgmeth, by = c("index"))
colnames(total_dev)[match("ADU_HSC", colnames(total_dev))] = "ADU"
total_dev_meth = total_dev[,5:ncol(total_dev)]
phyl_dist = stats::dist(t(total_dev_meth), method = "man")
tree_anno = data.frame(color = COLORSCHEME[colnames(total_dev_meth)],
                       row.names = colnames(total_dev_meth))
tree_anno["LM", "color"] = "#2188c9"
tree_anno["IM", "color"] = "#fbbb25"
plot_tree(tree = fastme.bal(X = phyl_dist), col_df = tree_anno, col_id = 1,
          type = 'unrooted', rotate = 155)
###############################################################################
