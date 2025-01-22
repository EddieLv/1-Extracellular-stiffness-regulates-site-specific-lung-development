source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("205-30", "205-60", "205-90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar"

tmp.dir1 = paste(work.dir, "6pseudotime", sep = "/")
dir.create(tmp.dir1, recursive = T)

min_max = function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
cols.ctp = palette_scp(levels(srat.merge.rna$celltype2_final), palette = "Set1")
srat.merge.rna = RunUMAP(srat.merge.rna, reduction = "pca", reduction.name = "pca_umap", dims = 1:29, umap.method = "uwot", metric = "cosine", min.dist = 0.5, n.neighbors = 20, return.model = T, seed.use = 666)
srat.merge.rna = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final %in% c("Tip ETV5+", "AT2", "Proliferating AT2", "Proliferating AT1", "AT1")])
srat.merge.rna$celltype2_final = factor(srat.merge.rna$celltype2_final, levels = c("Tip ETV5+", "AT2", "Proliferating AT2", "Proliferating AT1", "AT1"))
noise.cells = c("205-30_ATACTTCCATGACAGG-1", "205-60_GCGTGCAAGCAACAGC-1", "205-90_GGGCTACGTCACTACA-1", "205-90_GTGAGCCCAAGGCTTT-1")
srat.merge.rna = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[!colnames(srat.merge.rna) %in% noise.cells])
CellDimPlot(srat.merge.rna, group.by = "celltype2_final", reduction = "pca_umap", label = T, label_insitu = T, palette = "Paired")

srat.merge.rna[["pca_umap"]]@cell.embeddings[srat.merge.rna$celltype2_final == "Tip ETV5+", 2] = srat.merge.rna[["pca_umap"]]@cell.embeddings[srat.merge.rna$celltype2_final == "Tip ETV5+", 2] + 2
srat.merge.rna[["pca_umap"]]@cell.embeddings[srat.merge.rna$celltype2_final == "Tip ETV5+", 1] = srat.merge.rna[["pca_umap"]]@cell.embeddings[srat.merge.rna$celltype2_final == "Tip ETV5+", 1] - 1
CellDimPlot(srat.merge.rna, group.by = "celltype2_final", reduction = "pca_umap", label = T, label_insitu = T, palette = "Set1")

# n = 42
# pca_contribution(srat.merge.rna, n = n, reduction = "pca")
# srat.merge.rna = FindNeighbors(srat.merge.rna, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 15, graph.name = c("SCT_nn", "SCT_snn")) %>%
#   RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 15, return.model = T, seed.use = 666)
# CellDimPlot(srat.merge.rna, group.by = "celltype2_final", reduction = "pca_umap", label = T, label_insitu = T, palette = "Set1")

db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")
anno.df = db.merge.df.list$db.merge.df.final %>% dplyr::filter(species == "human" & database == "MsigDB" & collection %in% "GO_BP")
pathway.list = split(anno.df$pathway.gene, anno.df$pathway.name)
library(UCell)
srat.merge.rna = AddModuleScore_UCell(srat.merge.rna, features = pathway.list[c("Cell Cell Signaling By Wnt", "Hippo Signaling")], assay = "SCT", slot = "data", 
                                      name = "")
write.csv(srat.merge.rna@meta.data, "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/meta.data.csv", row.names = T)
seurat2scanpy(srat.merge.rna, ann.X = "SCT-data", ann.raw.X = "RNA-counts", h5ad_path = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/srat.h5ad")

### slingshot
DefaultAssay(srat.merge.rna) = "SCT"
srat.merge.rna = RunSlingshot(srt = srat.merge.rna, group.by = "celltype2_final", reduction = "pca_umap", start = "Tip ETV5+", end = c("AT2", "AT1"))
CellDimPlot(srat.merge.rna, group.by = "celltype2_final", reduction = "pca_umap", lineages = paste0("Lineage", 1:2), lineages_span = 0.6)

### stream
stream.df = read.csv("/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/stream/stream.obs.csv", row.names = 1)

srat.merge.rna$branch = stream.df[colnames(srat.merge.rna), "branch_id_alias"]
srat.merge.rna$State = 3
srat.merge.rna$State[srat.merge.rna$branch == "('S2', 'S1')"] = 2
srat.merge.rna$State[srat.merge.rna$branch == "('S0', 'S1')"] = 1
srat.merge.rna$State = factor(srat.merge.rna$State, levels = c(1, 2, 3))
srat.merge.rna$Pseudotime = as.numeric(stream.df[colnames(srat.merge.rna), "S0_pseudotime"])

srat.merge.rna$Pseudotime_minmax = min_max(srat.merge.rna$Pseudotime)
srat.merge.rna$CytoTRACE2__minmax = min_max(srat.merge.rna$preKNN_CytoTRACE2_Score)

srat.merge.rna$Lineage1 = NA
srat.merge.rna$Lineage2 = NA
srat.merge.rna$Lineage1[srat.merge.rna$State %in% c(1, 2)] = min_max(srat.merge.rna$Pseudotime[srat.merge.rna$State %in% c(1, 2)])
srat.merge.rna$Lineage2[srat.merge.rna$State %in% c(1, 3)] = min_max(srat.merge.rna$Pseudotime[srat.merge.rna$State %in% c(1, 3)])

p1 = FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("Pseudotime_minmax"), palette = "Spectral", label = F, theme_use = "theme_blank", pt.size = 1.5)
p2 = FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("CytoTRACE2__minmax"), palette = "Spectral", label = F, theme_use = "theme_blank", pt.size = 1.5)
p3 = FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("Lineage1"), palette = "Spectral", label = F, theme_use = "theme_blank", pt.size = 1.5)
p4 = FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("Lineage2"), palette = "Spectral", label = F, theme_use = "theme_blank", pt.size = 1.5)
p5 = CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "State", label = F, theme_use = "theme_blank", palcolor = cols.ctp, pt.size = 1.5)
p6 = CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", label = F, palcolor = cols.ctp,
                 lineages = paste0("Lineage", 1:2), lineages_span = 0.6, theme_use = "theme_blank")
pdf(paste(tmp.dir1, "alveolar.stream.umap.pdf", sep = "/"), width = 8*2, height = 8*1.5)
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()

pdf(paste(tmp.dir1, "alveolar.stream.traj.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2_final", palcolor = cols.ctp,
            label = F, theme_use = "theme_blank", pt.size = 1.5)
dev.off()

db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")
anno.df = db.merge.df.list$db.merge.df.final %>% dplyr::filter(species == "human" & database == "MsigDB" & collection %in% "GO_BP")
pathway.list = split(anno.df$pathway.gene, anno.df$pathway.name)
srat.merge.rna@assays$SCT@meta.features[, "Positive Regulation Of Hippo Signaling"] = rownames(srat.merge.rna[["SCT"]]) %in% pathway.list[["Positive Regulation Of Hippo Signaling"]]

srat.merge.rna$celltype2_final = factor(srat.merge.rna$celltype2_final, levels = c("Tip ETV5+", "AT2", "Proliferating AT2", "Proliferating AT1", "AT1"))
srat.merge.rna = SetIdent(srat.merge.rna, value = "celltype2_final")
srat.merge.rna = PrepSCTFindMarkers(srat.merge.rna, assay = "SCT")
heat.genes = FindAllMarkers(srat.merge.rna, assay = "SCT", slot = "data", logfc.threshold = 0.5, min.pct = 0.3, only.pos = T) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% dplyr::pull(gene) %>% unique()

marker.set = list(
  "('S3', 'S1')" = c("CLIC5", "CAV1", "AGER"),
  "('S2', 'S1')" = c("ABCA3", "LAMP3", "SFTPC"),
  "('S0', 'S1')" = c("MKI67", "TOP2A")
)

lvls = c("Totipotent", "Pluripotent", "Multipotent", "Oligopotent", "Unipotent", "Differentiated")
srat.merge.rna$preKNN_CytoTRACE2_Potency = droplevels(srat.merge.rna$preKNN_CytoTRACE2_Potency)
ht = DynamicHeatmap(
  srt = srat.merge.rna, assay = "SCT", slot = "data", lineages = c("Lineage2", "Lineage1"), 
  features = heat.genes, features_label = unlist(marker.set),
  use_fitted = T, n_split = 4, split_method = "kmeans", split_order = c(2, 1, 3, 4), decreasing = F, reverse_ht = "Lineage2",
  heatmap_palcolor = colorRampPalette(c("navy","white","firebrick3"))(100), pseudotime_palette = "Spectral",
  separate_annotation = list("preKNN_CytoTRACE2_Potency", "celltype2_final"), 
  separate_annotation_palcolor = list( palette_scp(levels(droplevels(srat.merge.rna$preKNN_CytoTRACE2_Potency)), palcolor = palette_scp(lvls, palette = "Paired")),
                                       cols.ctp),
  # feature_annotation = c("Positive Regulation Of Hippo Signaling"), feature_annotation_palcolor = list(c("red", "white")),
  height = 6, width = 2, border = F,
  show_row_names = F, use_raster = T
)
pdf(paste(tmp.dir1, "alveolar.stream.heat.pdf", sep = "/"), width = 8, height = 8)
print(ht$plot)
dev.off()
qsave(srat.merge.rna, "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/srat.alv.qs")
srat.merge.rna = qread("/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/srat.alv.qs")

p = GroupHeatmap(srat.merge.rna, assay = "SCT", slot = "data", 
                 cell_annotation = c("Pseudotime_minmax"), cell_annotation_palcolor = list(colors.organoid.alveolar.sample[levels(srat.merge.rna$orig.ident)]),
                 group.by = "celltype2_final", group_palcolor = list(cols.ctp[levels(srat.merge.rna$celltype2_final)]),
                 split.by = "orig.ident", cell_split_palcolor = list(colors.organoid.alveolar.sample[levels(srat.merge.rna$orig.ident)]),
                 features = c("ETV5", "SFTPC", "SFTPB", "MKI67", "TOP2A", "YAP1", "AGER", "HOPX"), cluster_rows = T, cluster_columns = F, width = 6, height = 1.2)
pdf("/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/6pseudotime/markers.dynamic.pdf", width = 8*2, height = 8)
p$plot
dev.off()

