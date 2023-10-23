#!Rscript
""" Identifying and plotting epigentic scars of healthy HSC development

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
bsdata <- readRDS("/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/bsseq/bsseq_HSC_combTotal_snpRemoved.rds")
phenoTotal  <- colData(bsdata)
drop        <- (phenoTotal$Sample_Type == "normal") & (phenoTotal$Celltype != "HSC")
drop        <- drop | phenoTotal$Name %in% c("FL2_HSC_PBAT_2", "JU3_HSC_PBAT_11")
bsdata      <- bsdata[, !drop]
temp_anno = tree_cols2(colnames(bsdata))
keep = !temp_anno$donorCtype %in% c("ADU_CMP", "ADU_GMP", "ADU_MEP", "ADU_MPP")
bsdata      <- bsdata[, keep]

dmrTest = function(bsdata, dmrs, groups){
  ugroups = unique(groups)
  if(length(ugroups) > 2){
    stop("Only 2 groups are supported.")
  }

  groupTrans = data.table(
    name  = colnames(bsdata),
    label = groups
  )
  coverage_m = getCoverage(bsdata, regions = dmrs, type = "M", what = "perRegionTotal")
  coverage = getCoverage(bsdata, regions = dmrs, type = "Cov", what = "perRegionTotal")
  coverage_u = coverage - coverage_m
  rm(coverage)

  coverage_m = as.data.table(coverage_m)
  coverage_u = as.data.table(coverage_u)
  coverage_m[, index := 1:.N]; coverage_u[, index := 1:.N]
  coverage_m = melt(coverage_m, id.vars = "index", value.name = "Methylated")
  coverage_u = melt(coverage_u, id.vars = "index", value.name = "Unmethylated")
  coverage_m[, group := factor(variable, levels = groupTrans$name, labels = groupTrans$label)]
  coverage_u[, group := factor(variable, levels = groupTrans$name, labels = groupTrans$label)]

  coverage_m = coverage_m[, .(Methylated = sum(Methylated, na.rm = T)), by = .(index, group)]
  coverage_u = coverage_u[, .(Unmethylated = sum(Unmethylated, na.rm = T)), by = .(index, group)]
  coverage = merge(coverage_m, coverage_u, by = c("index", "group"))

  meth = getMeth(bsdata, regions = dmrs, type = "raw", what = "perRegion")
  meth = as.data.table(meth)
  meth[, index := 1:.N]
  meth = melt(meth, id.vars = "index", value.name = "Methylation")
  meth[, group := factor(variable, levels = groupTrans$name, labels = groupTrans$label)]
  message("now not using na.rm =T")
  #meth = meth[, .(SDmeth = sd(Methylation, na.rm = T)), by = .(index, group)]
  meth = meth[, .(SDmeth = sd(Methylation)), by = .(index, group)]

  meth = dcast(meth, index ~ group, value.var = "SDmeth")
  colnames(meth) = c("index", paste("SDmeth", colnames(meth)[2:ncol(meth)],sep ="_"))

  res = coverage[, chisq.test(.SD)$p.value, by = index, .SDcols = c("Methylated", "Unmethylated")]
  colnames(res) = c("index", "p_val")
  idx = res$index
  total_meth = dmrs@elementMetadata[idx, ]$meanMethy1+dmrs@elementMetadata[idx, ]$meanMethy2
  rel_methdiff = abs(dmrs@elementMetadata[idx,]$diff.Methy)/total_meth

  res[, p_adj := p.adjust(p_val, method = "BH")]
  res[, meth_diff := dmrs@elementMetadata[index, ]$diff.Methy]
  res[, meth_diff_rel := rel_methdiff]
  return(merge(res, meth, by = "index"))
}


normals = colnames(bsdata)[bsdata$Sample_Type == "normal"]
normals = normals[!grepl(pattern = "CMP", normals)]
normals = normals[!grepl(pattern = "GMP", normals)]
normals = normals[!grepl(pattern = "MEP", normals)]
normals = normals[!grepl(pattern = "MPP", normals)]

dmrs_natal = readRDS(
  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_PrePostnatal/sig_dmrs_5inHalf_sub_anno.rds"
)$Postnatal_vs_Prenatal
anno = tree_cols2(colnames(dmrs_natal@elementMetadata))
anno = anno[sample %in% colnames(bsdata)]
meth_natal = as.data.table(dmrs_natal@elementMetadata[, anno$sample])
#meth_natal = getMeth(bsdata[, normals], region = dmrs_natal, type = "raw", what = "perRegion")
anno_natal = tree_cols2(colnames(meth_natal))
groups_natal = anno_natal$donor
groups_natal[groups_natal %in% c("FES", "FEL")] = "Prenatal"
groups_natal[groups_natal != "Prenatal"] = "Postnatal"

dmrs_neo = readRDS(
  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_NeoPostneo/sig_dmrs_5inHalf_sub_anno.rds"
)$Neo_vs_Postneo
dmrs_neof = readRDS(
  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_FetalNeo/sig_dmrs_5inHalf_sub_anno.rds"
)$Neo_vs_Fetal
dmrs_postf = readRDS(
  "/omics/groups/OE0219/internal/jmmlc_pbat_normals/data/odcf_md/analysis/220512_DMR_hierachy_HSC_FetPostneo//sig_dmrs_5inHalf_sub_anno.rds"
)$Postneo_vs_Fetal

# To remove "Old" signal coming from the fetal to postfetal transition
# this only removes 4 ppotential epigenetic scars
drop  = c(findOverlaps(dmrs_neo, dmrs_postf)@from, findOverlaps(dmrs_neo, dmrs_neof)@from)
keep = setdiff(1:nrow(dmrs_neo@elementMetadata), unique(drop))
dmrs_neo = dmrs_neo[keep, ]

anno = tree_cols2(colnames(dmrs_neo@elementMetadata))
anno = anno[sample %in% colnames(bsdata)]
meth_neo = as.data.table(dmrs_neo@elementMetadata[, anno$sample])
anno_neo = tree_cols2(colnames(meth_neo))
groups_neo = anno_neo$donor
groups_neo[groups_neo == "NEO"] = "Neo"
groups_neo[groups_neo != "Neo"] = "Postneo"


sig_natal = dmrTest(bsdata[, anno_natal$sample], dmrs = dmrs_natal, groups = groups_natal)
sig_neo   = dmrTest(bsdata[, anno_neo$sample], dmrs = dmrs_neo, groups = groups_neo)
sig_natal = sig_natal[complete.cases(sig_natal)]
sig_neo   = sig_neo[complete.cases(sig_neo)]

sig_natal = sig_natal[p_adj < 0.01]
sig_neo   = sig_neo[p_adj < 0.01]
sig_natal = sig_natal[order(p_adj, decreasing = T)]
sig_neo = sig_neo[order(p_adj, decreasing = T)]
sel_natal = sig_natal[SDmeth_Prenatal < .1 & SDmeth_Postnatal < .1]
sel_neo = sig_neo[SDmeth_Neo < .1 & SDmeth_Postneo < .05]

# natal normal
meth_natal_sel = getMeth(
  collapseBSseq(bsdata[,anno_natal$sample], group = anno_natal$donor),
  regions = dmrs_natal[sel_natal$index,], type = "raw", what = "perRegion")

htmp_natal = Heatmap(meth_natal_sel,
                     top_annotation = HeatmapAnnotation(
                       border = T,
                       HSC = colnames(meth_natal_sel),
                       col = list(
                         HSC = COLORSCHEME[colnames(meth_natal_sel)])
                     ),
                     col = ANNOTATION_COLORSCALE_W, show_column_names = F,
                     border = "#000000", name = "Average methylation", column_title = "Prenatal vs Postnatal",
                     heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
                     #column_split = groups_natal,
                     cluster_column_slices = F, use_raster = F,
                     row_split = str_to_title(dmrs_natal[sel_natal$index, ]@elementMetadata$direction)
)
draw(htmp_natal, heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")


# natal all
meth_natal_sel = getMeth(bsdata,
                         regions = dmrs_natal[sel_natal$index,], type = "raw", what = "perRegion")
anno_temp = tree_cols(meth_natal_sel)
htmp_natal = Heatmap(meth_natal_sel,
                     top_annotation = HeatmapAnnotation(
                       border = T,
                       HSC = anno_temp$donor,
                       col = list(
                         HSC = COLORSCHEME[anno_temp$donor])
                     ),
                     col = ANNOTATION_COLORSCALE_W, show_column_names = F,
                     border = "#000000", name = "Average methylation", column_title = "Prenatal vs Postnatal",
                     heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
                     #column_split = groups_natal,
                     cluster_column_slices = F, use_raster = F,
                     row_split = str_to_title(dmrs_natal[sel_natal$index, ]@elementMetadata$direction)
)
draw(htmp_natal, heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")

# neo normals
anno_temp = tree_cols2(normals)
anno_temp$donor[anno_temp$donor %in% c("FEL", "FES")] = "FET"
meth_neo_sel = getMeth(
  collapseBSseq(bsdata[,anno_temp$sample], group = anno_temp$donor),
  regions = dmrs_neo[sel_neo$index,], type = "raw", what = "perRegion")
htmp_neo = Heatmap(meth_neo_sel,
                   top_annotation = HeatmapAnnotation(
                     border = T,
                     HSC = colnames(meth_neo_sel),
                     col = list(
                       HSC = COLORSCHEME[colnames(meth_neo_sel)])
                   ),
                   col = ANNOTATION_COLORSCALE_W, show_column_names = F,
                   border = "#000000", name = "Average methylation",
                   column_title = "Neonatal vs Postneonatal",
                   heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
                   #column_split = groups_neo,
                   cluster_column_slices = F, use_raster = F,
                   row_split = str_to_title(dmrs_neo[sel_neo$index, ]@elementMetadata$direction)
)
draw(htmp_neo, heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")


anno_temp = tree_cols2(colnames(bsdata))
meth_neo_sel = getMeth(
  collapseBSseq(bsdata, group = anno_temp$donorCtype),
  regions = dmrs_neo[sel_neo$index, ], type = "raw", what = "perRegion")
meth_neo_sel = meth_neo_sel[
  ,!colnames(meth_neo_sel) %in% c("ADU_CMP", "ADU_GMP", "ADU_MEP", "ADU_MPP")]
colnames(meth_neo_sel)[match("ADU_HSC", colnames(meth_neo_sel))] = "ADU"
htmp_neo = Heatmap(meth_neo_sel,
                   top_annotation = HeatmapAnnotation(
                     border = T,
                     HSC = colnames(meth_neo_sel),
                     col = list(
                       HSC = COLORSCHEME[colnames(meth_neo_sel)])
                   ),
                   col = ANNOTATION_COLORSCALE_W, show_column_names = F,
                   border = "#000000", name = "Average methylation",
                   column_title = "Neonatal vs Postneonatal",
                   heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
                   cluster_column_slices = F, use_raster = F,
                   row_split = str_to_title(dmrs_neo[sel_neo$index, ]@elementMetadata$direction)
)
draw(htmp_neo, heatmap_legend_side = "bottom", annotation_legend_side = 'right',
     legend_grouping = "original")


tumor_samples = colnames(bsdata)[grepl(pattern = "tumor", colnames(bsdata))]
tumor_anno = tree_cols2(tumor_samples)
tumor_meth_natal = getMeth(collapseBSseq(bsdata[, tumor_samples], group = tumor_anno$donor),
                         regions = dmrs_natal[sel_natal$index, ], type = "raw", what = "perRegion")
tumor_meth_neo = getMeth(collapseBSseq(bsdata[, tumor_samples], group = tumor_anno$donor),
                         regions = dmrs_neo[sel_neo$index, ], type = "raw", what = "perRegion")

total_dmrs = reduce(unlist(GRangesList(dmrs_natal[sel_natal$index,], dmrs_neo[sel_neo$index, ])))
meth_total = getMeth(bsdata, regions = total_dmrs, type = "raw", what = "perRegion")
total_anno = tree_cols(meth_total)
colnames(meth_total) = total_anno$name
rm(meth_total)

natal_meth = getMeth(collapseBSseq(bsdata, group = total_anno$donor),
                     regions = dmrs_natal[sel_natal$index], type = "raw", what = "perRegion")
natal_scar = cbind(as.data.table(dmrs_natal[sel_natal$index])[, 1:3],
                   natal_meth)
neo_meth   = getMeth(collapseBSseq(bsdata, group = total_anno$donor),
                     regions = dmrs_neo[sel_neo$index], type = "raw", what = "perRegion")
neo_scar   = cbind(as.data.table(dmrs_neo[sel_neo$index])[,1:3],
                   neo_meth)

natal_scar[, group := "PrePostnatal"]
neo_scar[, group := "NeoPostneo"]
scars = rbind(natal_scar, neo_scar)
scars[, index := .I]

prepare_liftover = scars[, c("seqnames", "start", "end", "group", "index")]
# ucsc hg19 to hg38 liftover
#fwrite(prepare_liftover, "~/scars.tsv", sep = "\t", col.names = F)

# bed_files = list.files(
#   "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq/filteredMethCallsPE/",
#   full.names = T, pattern = "*bed"
# )
# bed_files = bed_files[!grepl(pattern = "i7-only", ignore.case = T, bed_files)]
# scarRegions = fread("/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/scars_hg38.tsv")
# colnames(scarRegions) = c("Chr", "Start", "End", "Group", "Index")
# scarRegions[, Chr := gsub(pattern = "chr", replacement = "", Chr)]
# setkey(scarRegions, Chr, Start, End)
# avg_scar_cgmeth = rbindlist(mclapply(bed_files, function(i) {
#   data = fread(i)
#   setkey(data, chr, start, end)
#   data = foverlaps(data, scarRegions,
#                    type = "within",
#                    by.x = c("chr", "start", "end"),
#                    by.y = c("Chr", "Start", "End"),
#                    nomatch = 0L
#                    )
#   data = data[, .(n_total = sum(n_total, na.rm = T),
#                   n_meth = sum(n_meth, na.rm = T),
#                   avg_beta = mean(beta_value, na.rm = T)),
#               by = .(Index, Group)]
#   data[, sample := paste(translate_bsnames(i)[c("PID", "PLATE", "WELL")], collapse = "_")]
# }, mc.cores = 10))
# saveRDS(avg_scar_cgmeth, "/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/scar_methylation_scigmt.RDS")
avg_scar_cgmeth = readRDS("/omics/groups/OE0219/internal/Valentin/JMMLT/processing/bisulfiteseq_bulk/scar_methylation_scigmt.RDS")

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

avg_scar_cgmeth[, group := factor(sample, levels = x_dimr$sample, labels = x_dimr$group)]
avg_scar_cgmeth = avg_scar_cgmeth[!is.na(group)]
avg_scar_cgmeth = avg_scar_cgmeth[, .(avg_beta = mean(avg_beta)), by = .(Index, Group, group)]
avg_scar_cgmeth = dcast(avg_scar_cgmeth, Group + Index ~ group, value.var = "avg_beta")
colnames(avg_scar_cgmeth)[1:2] = c("group", "index")

total_scars = merge(scars, avg_scar_cgmeth, by = c("index", "group"))
total_scar_meth = total_scars[,6:ncol(total_scars)]

tree_anno = data.frame(color = COLORSCHEME[colnames(total_scar_meth)],
                       row.names = colnames(total_scar_meth))
tree_anno["LM", "color"] = "#2188c9"
tree_anno["IM", "color"] = "#fbbb25"
TEMP_COLORSCHEME = tree_anno$color
names(TEMP_COLORSCHEME) = rownames(tree_anno)
names(TEMP_COLORSCHEME)[match("HML", names(TEMP_COLORSCHEME))] = "HM-like"
names(TEMP_COLORSCHEME)[match("LML", names(TEMP_COLORSCHEME))] = "LM-like"

hsc = colnames(total_scar_meth)
data_anno = ifelse(colnames(total_scar_meth) %in% c("LM", "IM", "HML", "LML"),
              yes = "scIGMT", no = "Bulk BS-seq")
data_anno = factor(data_anno, levels = c("scIGMT", "Bulk BS-seq"))
temp_cols = list(HSC = TEMP_COLORSCHEME, Data = c(scIGMT = "#252525", `Bulk BS-seq` = "#F0F0F0"))
upper = total_scars$group == "PrePostnatal"
t(apply(total_scar_meth[upper, ], MARGIN = 1, scale_na))
directions = rowMeans(total_scar_meth[upper, c("FEL", "FES")]) > rowMeans(
  total_scar_meth[upper, c("NEO", "JUV", "ADU")])
directions = c(directions,
               rowMeans(total_scar_meth[!upper, c("NEO")]) > rowMeans(
  total_scar_meth[!upper, c("JUV", "ADU")]))
directions = ifelse(directions > 0, yes = "Hypo", no = "Hyper")
hsc = factor(hsc, levels = c("FEL", "FES", "NEO", "JUV", "ADU",
                             "LM", "LML", "IM", "HML",
                             paste0("P", 1:8)),
             labels = c("FEL", "FES", "NEO", "JUV", "ADU",
                        "LM", "LM-like", "IM", "HM-like",
                        paste0("P", 1:8)))

htmp_scar_natal = Heatmap(total_scar_meth[upper, ],
         top_annotation = HeatmapAnnotation(
           border = T,
           HSC  = hsc, Data = data_anno,
           col = temp_cols),
         column_title = "Prenatal vs Postnatal",
         row_split = directions[upper],
         col = ANNOTATION_COLORSCALE_W, show_column_names = T,
         name = "Average methylation",
         border = "#000000",
         heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
         cluster_column_slices = F, use_raster = F,
         cluster_rows = T, cluster_columns = T)
htmp_scar_neo = Heatmap(total_scar_meth[!upper,],
  top_annotation = HeatmapAnnotation(
    border = T,
    HSC  = hsc, Data = data_anno,
    col = temp_cols),
  #row_split = directions[!upper],
  column_title = "Neonatal vs Postneonatal",
  col = ANNOTATION_COLORSCALE_W, show_column_names = T,
  name = "Average methylation",
  border = "#000000",
  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")),
  cluster_column_slices = F, use_raster = F,
  cluster_rows = T, cluster_columns = T)

draw(htmp_scar_natal,
     heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")
draw(htmp_scar_neo,
     heatmap_legend_side = "bottom", annotation_legend_side = 'right', legend_grouping = "original")
###############################################################################
