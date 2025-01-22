source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/upstream_res"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway"

tmp.dir1 = paste(work.dir, "1RNA_res", sep = "/")

min.nCount = 1000; min.features = 500

library(SCP)
##########################
###### Seurat Object #####
##########################
dir.create(tmp.dir1, recursive = T)
srat.list.rna = list()
for (sample in samples.rna) {
  message(sample)
  # read raw counts matrix
  mat = ReadMtx(mtx = paste(upstream.dir, sample, "filter_matrix", "matrix.mtx.gz", sep = "/"), 
                cells = paste(upstream.dir, sample, "filter_matrix", "barcodes.tsv.gz", sep = "/"), cell.column = 1, 
                features = paste(upstream.dir, sample, "filter_matrix", "features.tsv.gz", sep = "/"), feature.column = 1)
  # create seurat object
  seurat = create_seurat_RNA(count_mat = mat, gene_id_input = F, metadata_path = NULL, 
                             assay = "RNA", project = sample, species = "human", 
                             run_cellQC = T, min.nCount = min.nCount, max.nCount.quantile = 0.99, min.features = min.features, max.nFeature.quantile = 0.99,
                             run_geneQC = T, min.cells = 5, genes_remove = NULL, 
                             cellname_prefix = NULL, cellname_suffix = NULL, 
                             isST = F, run_doublet = T)
  seurat$nCount_RNA_log10 = log10(seurat$nCount_RNA)
  seurat$nFeature_RNA_log10 = log10(seurat$nFeature_RNA)
  srat.list.rna[[sample]] = seurat
}
# remove doublets
plots1 = map(names(srat.list.rna), ~CellStatPlot(srat.list.rna[[.]], stat.by = c("pd_DoubletFinder", "pd_DoubletDetection", "pd_scDblFinder", "pd_Scrublet"), plot_type = "upset", stat_level = "doublet", title = .))
wrap_plots(plots1, ncol = length(srat.list.rna))
for (sample in names(srat.list.rna)) {
  srat.list.rna[[sample]]$doublet = "singlet"
  srat.list.rna[[sample]]$doublet[srat.list.rna[[sample]]$pd_DoubletDetection == "doublet" & 
                                    srat.list.rna[[sample]]$pd_DoubletFinder == "doublet" & 
                                    srat.list.rna[[sample]]$pd_scDblFinder == "doublet"] = "doublet"
}
plots2 = map(names(srat.list.rna), ~FeatureStatPlot(srat.list.rna[[.]], stat.by = c("nCount_RNA", "nFeature_RNA"), group.by = "doublet", comparisons = list(c("singlet", "doublet")), sig_label = "p.format", nrow = 2, title = .))
wrap_plots(plots2, ncol = length(srat.list.rna))
plots3 = map(names(srat.list.rna), ~CellDimPlot(srat.list.rna[[.]], reduction = "pca_umap", group.by = "doublet", title = .))
wrap_plots(plots3, ncol = length(srat.list.rna))

doublet_metadata_list = map(names(srat.list.rna), ~srat.list.rna[[.]]@meta.data[c("cell_id", "DoubletFinder", "pd_DoubletFinder", "cxds", "bcds", "hybrid", "Scrublet", 
                                                                                  "pd_Scrublet", "pd_scDblFinder", "DoubletDetection", "pd_DoubletDetection", "doublet")])
names(doublet_metadata_list) = names(srat.list.rna)
qsave(doublet_metadata_list, paste(tmp.dir1, "doublet_metadata_list.qs", sep = "/"))
doublet_metadata_list = qread(paste(tmp.dir1, "doublet_metadata_list.qs", sep = "/"))
# rerun
srat.list.rna = list()
for (sample in samples.rna) {
  message(sample)
  # read raw counts matrix
  mat = ReadMtx(mtx = paste(upstream.dir, sample, "filter_matrix", "matrix.mtx.gz", sep = "/"), 
                cells = paste(upstream.dir, sample, "filter_matrix", "barcodes.tsv.gz", sep = "/"), cell.column = 1, 
                features = paste(upstream.dir, sample, "filter_matrix", "features.tsv.gz", sep = "/"), feature.column = 1)
  # create seurat object
  seurat = create_seurat_RNA(count_mat = mat[, doublet_metadata_list[[sample]]$cell_id[doublet_metadata_list[[sample]]$doublet == "singlet"]], gene_id_input = F, metadata_path = NULL, 
                             assay = "RNA", project = sample, species = "human", 
                             run_cellQC = T, min.nCount = min.nCount, max.nCount.quantile = 0.99, min.features = min.features, max.nFeature.quantile = 0.99,
                             run_geneQC = T, min.cells = 5, genes_remove = NULL, 
                             cellname_prefix = NULL, cellname_suffix = NULL, 
                             isST = F, run_doublet = F)
  seurat$nCount_RNA_log10 = log10(seurat$nCount_RNA)
  seurat$nFeature_RNA_log10 = log10(seurat$nFeature_RNA)
  srat.list.rna[[sample]] = seurat
}
qsave(srat.list.rna, paste(tmp.dir1, "srat.list.rna.qs", sep = "/"))
srat.list.rna = qread(paste(tmp.dir1, "srat.list.rna.qs", sep = "/"))

#######################################
##### merge and SCTransform - 1st #####
#######################################
srat.merge = merge_samples_with_sct(srat.list.rna, assay = "RNA", 
                                    vars.to.regress = c("nFeature_RNA", "pct.mt", "pct.hb"), 
                                    variable.features.n = 2000, 
                                    add.cell.ids = names(srat.list.rna), 
                                    PrepSCT = F)
# cellcycle
cc.list = calc_cellcycle(srat.merge, sample.col = "orig.ident", assay = "SCT", species = "human", nbin = 20)
for (sample in names(cc.list)) {
  srat.list.rna[[sample]]@meta.data = cbind(srat.list.rna[[sample]]@meta.data, cc.list[[sample]])
  srat.list.rna[[sample]]$CC.Difference = srat.list.rna[[sample]]$G2M.Score - srat.list.rna[[sample]]$S.Score
}
srat.merge = merge_samples_with_sct(srat.list.rna, assay = "RNA", 
                                    vars.to.regress = c("nFeature_RNA", "pct.mt", "pct.hb", "CC.Difference"), 
                                    variable.features.n = 2000, 
                                    add.cell.ids = names(srat.list.rna))

##################
##### Gruffi #####
##################
library(gruffi)
library(ggpubr)
DefaultAssay(srat.merge) = "SCT"
srat.merge = RunPCA(srat.merge, assay = "SCT", npcs = 50, verbose = T) %>%
  RunUMAP(reduction = "pca", reduction.name = "umap", dims = 1:36, umap.method = "uwot", metric = "cosine", min.dist = 0.1, n.neighbors = 20, return.model = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:36, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn"))
combined.obj = Seurat.utils::SetupReductionsNtoKdimensions(obj = srat.merge, nPCs = 50, dimensions = 3:2, reduction = "umap")
DefaultAssay(combined.obj) = "SCT"
rm(srat.merge); gc()
combined.obj = AutoFindGranuleResolution(obj = combined.obj, assay = "SCT")
optimal.granule.res = combined.obj@misc$gruffi$optimal.granule.res
# granules with <30 cells are cell-by-cell re-assigned to a neighboring granule
combined.obj = ReassignSmallClusters(combined.obj, ident = optimal.granule.res) # Will be stored in meta data column as "seurat_clusters.reassigned".
optimal.granule.res = combined.obj@misc$gruffi$optimal.granule.res
# calculate GO score
go1 = "GO:0006096" # Glycolysis糖酵解(缺氧或低氧时高)
go2 = "GO:0034976" # ER-stress内质网应激(蛋白质折叠错误积累)
go3 = "GO:0042063" # Gliogenesis胶质发生(越高越好)
combined.obj = AssignGranuleAverageScoresFromGOterm(combined.obj, assay = "SCT", GO_term = go1, save.UMAP = F, new_GO_term_computation = T, clustering = optimal.granule.res, plot.each.gene = F)
combined.obj = AssignGranuleAverageScoresFromGOterm(combined.obj, assay = "SCT", GO_term = go2, save.UMAP = F, new_GO_term_computation = T, clustering = optimal.granule.res, plot.each.gene = F)
combined.obj = AssignGranuleAverageScoresFromGOterm(combined.obj, assay = "SCT", GO_term = go3, save.UMAP = F, new_GO_term_computation = T, clustering = optimal.granule.res, plot.each.gene = F)
# Create score names:
i1 = ParseGruffiGranuleScoreName(goID = go1)
i2 = ParseGruffiGranuleScoreName(goID = go2)
i3 = ParseGruffiGranuleScoreName(goID = go3)
combined.obj = FindThresholdsShiny(obj = combined.obj, stress.ident1 = i1, stress.ident2 = i2, notstress.ident3 = i3)
CellDimPlot(combined.obj, group.by = "is.Stressed") / CellStatPlot(combined.obj, stat.by = "is.Stressed", group.by = "orig.ident")
# filtering
stress_metadata = combined.obj@meta.data[, c("orig.ident", "cell_id", "Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063", "is.Stressed")]
qsave(stress_metadata, paste(tmp.dir1, "stress_metadata.qs", sep = "/"))

#######################################
##### merge and SCTransform - 2st #####
#######################################
doublet_metadata_list = qread(paste(tmp.dir1, "doublet_metadata_list.qs", sep = "/"))
########################################
# 24.06.28: more strict doublet remove #
for (sample in names(doublet_metadata_list)) {
  doublet_metadata_list[[sample]]$doublet = "singlet"
  doublet_metadata_list[[sample]]$doublet[doublet_metadata_list[[sample]]$pd_DoubletDetection == "doublet" & 
                                            doublet_metadata_list[[sample]]$pd_scDblFinder == "doublet"] = "doublet"
}
########################################
stress_metadata = qread(paste(tmp.dir1, "stress_metadata.qs", sep = "/"))
noise.spots1 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/mt_high.txt")
noise.spots1$sample = str_split(noise.spots1$V1, "_CELL", simplify = T)[, 1]
noise.spots1 = noise.spots1 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))
noise.spots2 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/heat_shocked.txt")
noise.spots2$sample = str_split(noise.spots2$V1, "_CELL", simplify = T)[, 1]
noise.spots2 = noise.spots2 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))
noise.spots3 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/ribo_high.txt")
noise.spots3$sample = str_split(noise.spots3$V1, "_CELL", simplify = T)[, 1]
noise.spots3 = noise.spots3 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))
noise.spots4 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/ribo_low.txt")
noise.spots4$sample = str_split(noise.spots4$V1, "_CELL", simplify = T)[, 1]
noise.spots4 = noise.spots4 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))
noise.spots5 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/mt_high_ribo_low.txt")
noise.spots5$sample = str_split(noise.spots5$V1, "_CELL", simplify = T)[, 1]
noise.spots5 = noise.spots5 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))
noise.spots6 = read.table("/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/SAA1_SAA2_high.txt")
noise.spots6$sample = str_split(noise.spots6$V1, "_CELL", simplify = T)[, 1]
noise.spots6 = noise.spots6 %>% dplyr::mutate(cell_id = str_replace(V1, paste0(sample, "_"), ""))

srat.list.rna = list()
for (sample in samples.rna) {
  message(sample)
  # read raw counts matrix
  mat = ReadMtx(mtx = paste(upstream.dir, sample, "filter_matrix", "matrix.mtx.gz", sep = "/"), 
                cells = paste(upstream.dir, sample, "filter_matrix", "barcodes.tsv.gz", sep = "/"), cell.column = 1, 
                features = paste(upstream.dir, sample, "filter_matrix", "features.tsv.gz", sep = "/"), feature.column = 1)
  # create seurat object
  remove.cells1 = doublet_metadata_list[[sample]]$cell_id[doublet_metadata_list[[sample]]$doublet == "doublet"]
  remove.cells2 = stress_metadata$cell_id[stress_metadata$orig.ident == sample & stress_metadata$is.Stressed == T]
  remove.cells3 = noise.spots1$cell_id[noise.spots1$sample == sample]
  remove.cells4 = noise.spots2$cell_id[noise.spots2$sample == sample]
  remove.cells5 = noise.spots3$cell_id[noise.spots3$sample == sample]
  remove.cells6 = noise.spots4$cell_id[noise.spots4$sample == sample]
  remove.cells7 = noise.spots5$cell_id[noise.spots5$sample == sample]
  remove.cells8 = noise.spots6$cell_id[noise.spots6$sample == sample]
  remove.cells = unique(c(remove.cells1, remove.cells2, remove.cells3, remove.cells4, remove.cells5, remove.cells6, remove.cells7, remove.cells8))
  seurat = create_seurat_RNA(count_mat = mat[, !colnames(mat) %in% remove.cells], gene_id_input = F, metadata_path = NULL, 
                             assay = "RNA", project = sample, species = "human", 
                             run_cellQC = T, min.nCount = min.nCount, max.nCount.quantile = 0.99, min.features = min.features, max.nFeature.quantile = 0.99,
                             run_geneQC = T, min.cells = 5, genes_remove = NULL, 
                             cellname_prefix = NULL, cellname_suffix = NULL, 
                             isST = F, run_doublet = F)
  seurat$nCount_RNA_log10 = log10(seurat$nCount_RNA)
  seurat$nFeature_RNA_log10 = log10(seurat$nFeature_RNA)
  srat.list.rna[[sample]] = seurat
}
qsave(srat.list.rna, paste(tmp.dir1, "srat.list.rna.qs", sep = "/"))
srat.list.rna = qread(paste(tmp.dir1, "srat.list.rna.qs", sep = "/"))

srat.merge = merge_samples_with_sct(srat.list.rna, assay = "RNA", 
                                    vars.to.regress = c("nFeature_RNA", "pct.mt", "pct.hb", "pct.ribo"), 
                                    variable.features.n = 2000, 
                                    add.cell.ids = names(srat.list.rna), 
                                    PrepSCT = F)
# cellcycle
cc.list = calc_cellcycle(srat.merge, sample.col = "orig.ident", assay = "SCT", species = "human", nbin = 20)
for (sample in names(cc.list)) {
  srat.list.rna[[sample]]@meta.data = cbind(srat.list.rna[[sample]]@meta.data, cc.list[[sample]])
  srat.list.rna[[sample]]$CC.Difference = srat.list.rna[[sample]]$G2M.Score - srat.list.rna[[sample]]$S.Score
}
srat.merge = merge_samples_with_sct(srat.list.rna, assay = "RNA", 
                                    vars.to.regress = c("nFeature_RNA", "pct.mt", "pct.hb", "CC.Difference"), 
                                    variable.features.n = 2000, 
                                    add.cell.ids = names(srat.list.rna))

###################
##### Denoise #####
###################
library(Rmagic)
data.magic.list = list()
data = GetAssayData(srat.merge, assay = "SCT", slot = "data")
for (sample in unique(srat.merge$orig.ident)) {
  message(sample)
  data.tmp = data[, colnames(srat.merge)[srat.merge$orig.ident == sample]]
  data.magic.list[[sample]] = magic(t(data.tmp), genes = "all_genes", knn.dist.method = "euclidean", n.jobs = 8, verbose = T)$result %>% t()
}
data.magic = reduce(data.magic.list, cbind) %>% as.matrix() %>% as.sparse()
srat.merge[["MAGIC_SCT"]] = CreateAssay5Object(data = data.magic[, colnames(srat.merge)])

#############################
##### Simple Clustering #####
#############################
n = 36
res.list = c(0.8)

DefaultAssay(srat.merge) = "SCT"
srat.merge = RunPCA(srat.merge, assay = "SCT", npcs = 50, features = VariableFeatures(srat.merge), seed.use = 666)
srat.merge = FindNeighbors(srat.merge, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn", verbose = F)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.merge@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.merge@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.merge@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
# k.param: 20, the biggger, the clustering more mixed
srat.merge = RunUMAP(srat.merge, reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.1, n.neighbors = 20, return.model = T, seed.use = 666)

pdf(paste(tmp.dir1, "basic_lighting.pdf", sep = "/"), width = 8*2, height = 8*1.5)
FeatureDimPlot(srat.merge, features = c("nCount_RNA_log10", "nFeature_RNA_log10", "S.Score", 
                                        "G2M.Score", "MKI67", "TOP2A", 
                                        "pct.mt", "pct.ribo", "pct.hb"), reduction = "pca_umap", show_stat = F, ncol = 3)
dev.off()

pdf(paste(tmp.dir1, "simple_clustering.pdf", sep = "/"), width = 8*2, height = 8*2)
CellDimPlot(srat.merge, group.by = c("orig.ident", "SCT_snn_res.0.8"), reduction = "pca_umap", label = T, label_insitu = T, label_repel = T, palette = "Set1")
dev.off()

pdf(paste(tmp.dir1, "simple_lighting.pdf", sep = "/"), width = 8, height = 8*0.5)
for (gene in c("TP63", "MUC5AC", "MUC5B", "FOXJ1", "SCGB1A1", "SCGB3A2", "KRT5", "ACE2", "CTSB", "CTSL", "TMPRSS2")) {
  plots = list()
  if (gene %in% rownames(srat.merge[["MAGIC_SCT"]]$data)) {
    DefaultAssay(srat.merge) = "SCT"
    p1 = FeatureDimPlot(srat.merge, features = gene, reduction = "pca_umap", title = paste0(sample, "_norm"), show_stat = F)
    DefaultAssay(srat.merge) = "MAGIC_SCT"
    p2 = FeatureDimPlot(srat.merge, features = gene, reduction = "pca_umap", title = paste0(sample, "_magic"), show_stat = F)
    plots[[paste0(sample, 1)]] = p1
    plots[[paste0(sample, 2)]] = p2
  } else {
    plots[[paste0(sample, 1)]] = ggplot()
    plots[[paste0(sample, 2)]] = ggplot()
  }
  print(wrap_plots(plots, ncol = 2))
}
dev.off()

# add info
stress_metadata = qread(paste(tmp.dir1, "stress_metadata.qs", sep = "/"))
srat.merge$orig.ident = factor(srat.merge$orig.ident, levels = c("R51d40G30", "R51d40G60", "R51d40G90"))
srat.merge@meta.data[, c("Glycolysis", "ER.Stress", "Gliogenesis")] = stress_metadata[colnames(srat.merge), c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063")]
library(UCell)
srat.merge = AddModuleScore_UCell(srat.merge, features = pcd.list, assay = "SCT", ncores = 4, name = ".pcd")
srat.merge = AddModuleScore_UCell(srat.merge, features = list("heat.shock.score" = heat_shock.genes), assay = "SCT", ncores = 4, name = "")

qsave(srat.merge, paste(tmp.dir1, "srat.merge.qs", sep = "/"))
srat.merge = qread(paste(tmp.dir1, "srat.merge.qs", sep = "/"))
