source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway"

tmp.dir1 = paste(work.dir, "6pseudotime", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
CellDimPlot(srat.merge.rna, group.by = "celltype2_final", split.by = "orig.ident", reduction = "pca_umap", cells.highlight = T, theme_use = "theme_blank")

#############
# cytotrace #
#############
library(CytoTRACE2)

# 如果要cross sample比较CytoTRACE2_Score, 建议每个sample单独跑
cytotrace2_result_list = list()
for (sample in samples.rna) {
  expression_data = srat.merge.rna[["SCT"]]$data[, colnames(srat.merge.rna)[srat.merge.rna$orig.ident == sample]]
  cytotrace2_result_list[[sample]] = cytotrace2(expression_data, species = "human", is_seurat = F, full_model = T, ncores = 4, seed = 6)
}
cytotrace2_result = purrr::reduce(cytotrace2_result_list, rbind)
qsave(cytotrace2_result, paste(tmp.dir1, "cytotrace2_result.qs", sep = "/"))
cytotrace2_result = qread(paste(tmp.dir1, "cytotrace2_result.qs", sep = "/"))

lvls = c("Totipotent", "Pluripotent", "Multipotent", "Oligopotent", "Unipotent", "Differentiated")
srat.merge.rna@meta.data = cbind(srat.merge.rna@meta.data, cytotrace2_result[colnames(srat.merge.rna), ])
srat.merge.rna$CytoTRACE2_Potency = factor(srat.merge.rna$CytoTRACE2_Potency, levels = lvls)
srat.merge.rna$preKNN_CytoTRACE2_Potency = factor(srat.merge.rna$preKNN_CytoTRACE2_Potency, levels = lvls)
pdf(paste(tmp.dir1, "cytotrace2.plot.pdf", sep = "/"), width = 8*1.5, height = 8)
FeatureStatPlot(srat.merge.rna, stat.by = "preKNN_CytoTRACE2_Score", group.by = "orig.ident", plot_type = "box", add_trend = T, palcolor = colors.organoid.sample,
                bg_palcolor = "white",
                comparisons = list(c("R51d40G30", "R51d40G60"), c("R51d40G60", "R51d40G90")), pairwise_method = "wilcox.test", sig_label = "p.format", aspect.ratio = 1)
CellStatPlot(srat.merge.rna, stat.by = "preKNN_CytoTRACE2_Potency", group.by = "orig.ident", plot_type = "area", palette = "RdYlGn", aspect.ratio = 0.8)
CellStatPlot(srat.merge.rna, stat.by = "preKNN_CytoTRACE2_Potency", group.by = "celltype2_final", palette = "RdYlGn", 
             stat_type = "percent", plot_type = "bar", label = F, aspect.ratio = 0.8)
dev.off()

##########
# TDESeq #
##########
srat.merge.rna$preKNN_CytoTRACE2_Score = srat.merge.rna$preKNN_CytoTRACE2_Score
tde = CreateTDEseqObject(counts = srat.merge.rna[["SCT"]]$counts,
                         data = srat.merge.rna[["SCT"]]$data,
                         meta.data = srat.merge.rna@meta.data)
tde_method = "cell"
tde_param = list(sample.var = "orig.ident",
                 stage.var = "orig.ident",
                 fit.model = "lmm",
                 tde.thr = 0.05,
                 num.core = 6)
tde = tdeseq(object = tde, tde.method = tde_method, tde.param=tde_param)

p = plot_pseudotime_heatmap_genger(cds[heat.genes, ], 
                                   labels_row = unlist(marker.set), repel.text = T, repel.degree = 0,
                                   anno.cols = c("timepoint", "celltype.re"), hmcols = colorRampPalette(c("navy","white","firebrick3"))(1000), 
                                   annos.colors = list("timepoint" = colors.timepoint[levels(cds$timepoint)], "celltype.re" = colors.lung.merge[levels(droplevels(cds$celltype.re))]), 
                                   num_clusters = 9, clustering_method = "ward.D2", 
                                   row_tree_order = c(1, 3, 7, 8, 9, 6, 2, 5, 4), row_tree_rev = c(3, 7, 6, 9, 5), row_tree_sort_auto = F, row_tree_anno = F,
                                   nbin = 1000, return.mat = F) # ward.D2* average** complete*** mcquitty median


