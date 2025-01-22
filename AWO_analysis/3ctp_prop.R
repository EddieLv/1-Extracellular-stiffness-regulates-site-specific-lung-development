source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway"

tmp.dir1 = paste(work.dir, "3celltype_proportion", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))

OR.res = calc_or(srat.merge.rna@meta.data, samlpe.vec = srat.merge.rna$orig.ident, cluster.vec = srat.merge.rna$celltype2_final)
mat.a = OR.res$OR.dist.tb %>% tibble::column_to_rownames(var = "rid") %>% as.matrix()
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
ann_colors = list(sample = colors.organoid.airway.sample)
p = pheatmap(t(mat.a), main = "celltype abundance across sample", scale = "none", fontsize = 10, cellheight = 30, cellwidth = 30, 
             annotation_row = annotation_row, annotation_colors = ann_colors, annotation_names_row = F,
             cluster_row = F, clustering_distance_rows = "correlation", 
             cluster_col = F, border = NA, legend_breaks = c(0.5, 1.5), 
             display_numbers = t(mat.b), fontsize_number = 8, number_color = "black",
             color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100))

pdf(paste(tmp.dir1, "ctp_sample.prop.pdf", sep = "/"), width = 8*2, height = 8)
CellStatPlot(srat.merge.rna, stat.by = "orig.ident", group.by = "celltype2_final", bg_alpha = 0, palcolor = colors.organoid.airway.sample, stat_type = "percent", plot_type = "trend", aspect.ratio = 0.8)
CellStatPlot(srat.merge.rna, stat.by = "celltype2_final", group.by = "orig.ident", bg_alpha = 0, palette = "Set1", stat_type = "percent", plot_type = "trend", aspect.ratio = 0.8)
grid::grid.newpage()
p
dev.off()

library(Nebulosa)
DefaultAssay(srat.merge.rna) = "MAGIC_SCT"
pdf(paste(tmp.dir1, "goblet.light.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", cells.highlight = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final == "Goblet"], palette = "Set1")
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("MUC5AC"), assay = "MAGIC_SCT", bg_cutoff = 0.01, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("TFF1"), assay = "MAGIC_SCT", bg_cutoff = 0.2, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("MUC13"), assay = "MAGIC_SCT", bg_cutoff = 0.15, bg_color = "lightgrey", pt.size = 1)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("TSPAN8"), assay = "MAGIC_SCT", bg_cutoff = 0.3, bg_color = "lightgrey", pt.size = 1)
plot_density(srat.merge.rna, features = c("MUC5AC", "TFF1", "MUC13", "TSPAN8"), slot = "data", reduction = "pca_umap", pal = "inferno", joint = T)[[5]]
dev.off()

##
srat.merge.rna = PrepSCTFindMarkers(srat.merge.rna, assay = "SCT")
###################
##### Denoise #####
###################
library(Rmagic)
data.magic.list = list()
data = GetAssayData(srat.merge.rna, assay = "SCT", slot = "data")
for (sample in unique(srat.merge.rna$orig.ident)) {
  message(sample)
  data.tmp = data[, colnames(srat.merge.rna)[srat.merge.rna$orig.ident == sample]]
  data.magic.list[[sample]] = magic(t(data.tmp), genes = "all_genes", knn.dist.method = "euclidean", n.jobs = 8, verbose = T)$result %>% t()
}
data.magic = reduce(data.magic.list, cbind) %>% as.matrix() %>% as.sparse()
srat.merge.rna[["MAGIC_SCT"]] = CreateAssay5Object(data = data.magic[, colnames(srat.merge.rna)])
##

srat.goblet = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final == "Goblet"])
# srat.goblet = merge(x = srat.goblet0, y = map(1:10, ~srat.goblet0))
# srat.goblet = JoinLayers(srat.goblet, assay = "MAGIC_SCT")
DefaultAssay(srat.goblet) = "MAGIC_SCT"
p = FeatureStatPlot(srat.goblet, stat.by = c("MUC5AC", "TFF1", "MUC13", "TSPAN8"), group.by = "orig.ident", plot_type = "box", palette = "Set1", add_trend = T, bg_alpha = 0, 
                    comparisons = list(c("R51d40G30", "R51d40G60"), c("R51d40G60", "R51d40G90")), pairwise_method = "wilcox.test", sig_label = "p.format", aspect.ratio = 1)
vlnplot_genger = function(df_clean, title, y1, y2) {
  print(title)
  library(tidyverse)
  library(ggbeeswarm)
  library(scales)
  library(ggfun)
  
  df_clean$value[df_clean$orig.ident == "R51d40G30"] <- df_clean$value[df_clean$orig.ident == "R51d40G30"] + 1e-10
  
  df_clean_mean <- df_clean %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(mean = mean(value), median = median(value))
  
  signif_out <- c()
  for (i in list(c("R51d40G30", "R51d40G60"), c("R51d40G60", "R51d40G90"))) {
    out <- wilcox.test(df_clean %>% dplyr::filter(orig.ident == i[1]) %>% dplyr::pull(value),
                       df_clean %>% dplyr::filter(orig.ident == i[2]) %>% dplyr::pull(value))
    signif_out <- c(signif_out, out$p.value)
  }
  
  breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    10^(seq.int(low, high))
  }
  
  p = ggplot(data = df_clean) + 
    geom_quasirandom(aes(x = orig.ident, y = value, shape = orig.ident, fill = orig.ident), method = "pseudorandom", size = 3, alpha = 0.8, color = "white") + 
    scale_shape_manual(values = c(23, 22, 21)) + 
    scale_fill_manual(values = colors.organoid.airway.sample) + 
    scale_y_log10(breaks = breaks_log10, labels = trans_format(log10, math_format(10^.x))) +
    annotation_logticks(sides = "l", outside = T) + 
    coord_cartesian(clip = "off") + 
    # scale_x_discrete(labels = c("Day 0", "Day 0", "Day 21", "Day 21", "Day 35", "Day 35")) +
    labs(x = NULL, y = "log10(normalized expression)") + 
    # geom_hline(yintercept = 100, linetype = "dashed") + 
    # mean
    annotate(geom = "segment", x = 0.6, xend = 1.4, y = df_clean_mean$mean[1], yend = df_clean_mean$mean[1], linetype = "dashed", linewidth = 1) + 
    annotate(geom = "segment", x = 1.6, xend = 2.4, y = df_clean_mean$mean[2], yend = df_clean_mean$mean[2], linetype = "dashed", linewidth = 1) + 
    annotate(geom = "segment", x = 2.6, xend = 3.4, y = df_clean_mean$mean[3], yend = df_clean_mean$mean[3], linetype = "dashed", linewidth = 1) + 
    # Case1
    annotate(geom = "segment", x = 1, xend = 2, y = y1, yend = y1, linewidth = 1) +
    annotate(geom = "segment", x = 1, xend = 1, y = 0.8*y1, yend = y1, linewidth = 1) + 
    annotate(geom = "segment", x = 2, xend = 2, y = 0.8*y1, yend = y1, linewidth = 1) + 
    annotate(geom = "text", x = 1.5, y = 1.5*y1, label = paste0("p = ", scientific(signif_out[1])), size = 5) +
    # Case2
    annotate(geom = "segment", x = 2, xend = 3, y = y2, yend = y2, linewidth = 1) + 
    annotate(geom = "segment", x = 2, xend = 2, y = 0.8*y2, yend = y2, linewidth = 1) + 
    annotate(geom = "segment", x = 3, xend = 3, y = 0.8*y2, yend = y2, linewidth = 1) + 
    annotate(geom = "text", x = 2.5, y = 1.5*y2, label = paste0("p = ", scientific(signif_out[2])), size = 5) +
    ggtitle(label = title) + 
    thm_font() + 
    thm_rect() + 
    theme(axis.text = element_text(size = 15),
          axis.text.y.left = element_text(margin = margin(r = 10)),
          legend.background = element_roundrect(color = "#808080", linetype = 1)) + 
    theme(aspect.ratio = 0.9)
  
  return(p)
}

pdf(paste(tmp.dir1, "goblet.marker.compare.pdf", sep = "/"), width = 8*2, height = 8)
wrap_plots(map(1:4, ~vlnplot_genger(p[[.]]$data, 
                                    title = c("MUC5AC", "TFF1", "MUC13", "TSPAN8")[.], 
                                    y1 = c(0.1, 0.3, 2.8, 0.8)[.], 
                                    y2 = c(0.5, 2.5, 8, 2)[.])), ncol = 2)
dev.off()


