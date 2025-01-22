source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("205-30", "205-60", "205-90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar"

tmp.dir1 = paste(work.dir, "3celltype_proportion", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))

OR.res = calc_or(srat.merge.rna@meta.data, samlpe.vec = srat.merge.rna$orig.ident, cluster.vec = srat.merge.rna$celltype2_final)
mat.a = OR.res$OR.dist.tb %>% tibble::column_to_rownames(var = "rid") %>% as.matrix()
mat.a[is.infinite(mat.a)] = sort(mat.a, decreasing = T)[2]
mat.b = OR.res$count.dist.melt.ext.tb[, c(1, 2, 6)]
mat.b = spread(mat.b, key = "cid", value = "adj.p.value") %>% tibble::column_to_rownames(var = "rid")
mat.b = mat.b[rownames(mat.a), colnames(mat.a)]
mat.b = ifelse(mat.b >= 0.05, "NS",
               ifelse(mat.b < 0.0001,"****",
                      ifelse(mat.b < 0.001,"***",
                             ifelse(mat.b < 0.01,"**",
                                    ifelse(mat.b < 0.05,"*","")))))
library(pheatmap)
annotation_row = data.frame(sample = levels(srat.merge.rna$orig.ident))
rownames(annotation_row) = levels(srat.merge.rna$orig.ident)
ann_colors = list(sample = colors.organoid.alveolar.sample)
p = pheatmap(t(mat.a), main = "celltype abundance across sample", scale = "none", fontsize = 10, cellheight = 30, cellwidth = 30, 
             annotation_row = annotation_row, annotation_colors = ann_colors, annotation_names_row = F,
             cluster_row = F, clustering_distance_rows = "correlation", 
             cluster_col = F, border = NA, legend_breaks = c(0.5, 1.5), 
             display_numbers = t(mat.b), fontsize_number = 8, number_color = "black",
             color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100))

pdf(paste(tmp.dir1, "ctp_sample.prop.pdf", sep = "/"), width = 8*2, height = 8)
CellStatPlot(srat.merge.rna, stat.by = "orig.ident", group.by = "celltype2_final", bg_alpha = 0, palcolor = colors.organoid.alveolar.sample, stat_type = "percent", plot_type = "trend", aspect.ratio = 0.8)
CellStatPlot(srat.merge.rna, stat.by = "celltype2_final", group.by = "orig.ident", palette = "Set1", bg_alpha = 0, stat_type = "percent", plot_type = "trend", aspect.ratio = 0.8)
grid::grid.newpage()
p
dev.off()


library(Nebulosa)
pdf(paste(tmp.dir1, "Tip.light.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", cells.highlight = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final == "Tip ETV5+"], palette = "Set1")
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("ETV5"), assay = "MAGIC_SCT", bg_cutoff = 0.01, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("NKX2-1"), assay = "MAGIC_SCT", bg_cutoff = 0.01, bg_color = "lightgrey", pt.size = 1)
plot_density(srat.merge.rna, features = c("ETV5", "NKX2-1"), slot = "data", reduction = "pca_umap", pal = "inferno", joint = T)[[3]]
dev.off()

pdf(paste(tmp.dir1, "AT1.light.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", 
            cells.highlight = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final %in% c("Proliferating AT1", "AT1")], palette = "Set1")
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("CLIC5"), assay = "MAGIC_SCT", bg_cutoff = 0.01, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("CAV1"), assay = "MAGIC_SCT", bg_cutoff = 0.2, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("AGER"), assay = "MAGIC_SCT", bg_cutoff = 0.15, bg_color = "lightgrey", pt.size = 1)
plot_density(srat.merge.rna, features = c("CLIC5", "CAV1", "AGER"), slot = "data", reduction = "pca_umap", pal = "inferno", joint = T)[[4]]
dev.off()

pdf(paste(tmp.dir1, "AT2.light.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", 
            cells.highlight = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final %in% c("Proliferating AT2", "AT2")], palette = "Set1")
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("ABCA3"), assay = "MAGIC_SCT", bg_cutoff = 0.01, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("LAMP3"), assay = "MAGIC_SCT", bg_cutoff = 0.2, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("SFTPC"), assay = "MAGIC_SCT", bg_cutoff = 0.15, bg_color = "lightgrey", pt.size = 1)
plot_density(srat.merge.rna, features = c("ABCA3", "LAMP3", "SFTPC"), slot = "data", reduction = "pca_umap", pal = "inferno", joint = T)[[4]]
dev.off()

srat.alveolar = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final %in% c("Tip ETV5+", "AT2", "AT1", "Proliferating AT2", "Proliferating AT1")])
srat.alveolar$celltype2_final = droplevels(srat.alveolar$celltype2_final)
srat.alveolar$celltype2_final = factor(srat.alveolar$celltype2_final, levels = c("Tip ETV5+", "Proliferating AT2", "Proliferating AT1", "AT2", "AT1"))
pdf(paste(tmp.dir1, "YAP1-TAZ.marker.compare.pdf", sep = "/"), width = 8*2, height = 8)
FeatureStatPlot(srat.alveolar, assay = "SCT", stat.by = "YAP1", group.by = "celltype2_final", split.by = "orig.ident",
                plot_type = "violin", palcolor = colors.organoid.alveolar.sample, add_trend = T, bg_alpha = 0, 
                comparisons = T, sig_label = "p.format", aspect.ratio = 0.6)
FeatureStatPlot(srat.alveolar, assay = "SCT", stat.by = "TAZ", group.by = "celltype2_final", split.by = "orig.ident",
                plot_type = "violin", palcolor = colors.organoid.alveolar.sample, add_trend = T, bg_alpha = 0, 
                comparisons = T, sig_label = "p.format", aspect.ratio = 1)
FeatureStatPlot(srat.merge.rna, assay = "SCT", stat.by = "YAP1", group.by = "orig.ident", 
                plot_type = "violin", palcolor = colors.organoid.alveolar.sample, add_trend = T, bg_alpha = 0, 
                comparisons = list(c("205-30", "205-60"), c("205-60", "205-90")), sig_label = "p.format", aspect.ratio = 1)
FeatureStatPlot(srat.merge.rna, assay = "SCT", stat.by = "TAZ", group.by = "orig.ident", 
                plot_type = "violin", palcolor = colors.organoid.alveolar.sample, add_trend = T, bg_alpha = 0, 
                comparisons = list(c("205-30", "205-60"), c("205-60", "205-90")), sig_label = "p.format", aspect.ratio = 1)
dev.off()

### scREF human
srat.merge = qread("/media/biogenger/D/Projects/CMY/Analysis/human_lung/ref_human_fetal_lung_wgcna/2cnet_Epithelium/1reannotation/srat.reanno.qs")
pdf(paste(tmp.dir1, "YAP1-TAZ.marker.compare.scREF-human.pdf", sep = "/"), width = 8*2, height = 8)
for (ctp in unique(srat.merge$celltype_re_final)) {
  srat.alveolar = subset(srat.merge, cells = colnames(srat.merge)[srat.merge$celltype_re_final %in% ctp])
  srat.alveolar$timepoint = droplevels(srat.alveolar$timepoint)
  p = FeatureStatPlot(srat.alveolar, assay = "RNA", stat.by = "YAP1", group.by = "timepoint", split.by = "timepoint", title = ctp,
                      plot_type = "violin", palcolor = colors.organoid.alveolar.sample, add_trend = T, bg_alpha = 0, 
                      comparisons = F, sig_label = "p.format", aspect.ratio = 1)
  print(p)
}
dev.off()
