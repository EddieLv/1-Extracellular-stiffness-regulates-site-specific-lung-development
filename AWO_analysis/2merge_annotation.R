source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/1RNA_res"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway"

tmp.dir1 = paste(work.dir, "2merge_annotation", sep = "/")
tmp.dir2 = paste(tmp.dir1, "annotation_ref", sep = "/")
dir.create(tmp.dir1, recursive = T)
dir.create(tmp.dir2, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.merge.qs", sep = "/"))

DefaultAssay(srat.merge.rna) = "SCT"
### reduction ###
srat.merge.rna = RunPCA(srat.merge.rna, assay = "SCT", npcs = 50, verbose = T)

for (n in 20:50) {
  srat.merge.rna@meta.data = srat.merge.rna@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
  res.list = c(1)
  srat.merge.rna = FindNeighbors(srat.merge.rna, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
    FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn")
  for (i in res.list) {
    level = sort(as.numeric(levels(srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]])))
    srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
  }
  # k.param: 20, the biggger, the clustering more mixed
  srat.merge.rna = RunUMAP(srat.merge.rna, reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
  # [PCs] and [Resolution] are determined!
  res_use = "SCT_snn_res.1"
  srat.merge.rna[["seurat_clusters"]] = srat.merge.rna[[res_use]]
  srat.merge.rna = SetIdent(srat.merge.rna, value = "seurat_clusters")
  p = CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "seurat_clusters", label = T, label_insitu = T, legend.position = "top", legend.direction = "horizontal")
  ggsave(paste(tmp.dir1, "nPC", paste0(n, ".png"), sep = "/"), plot = p, width = 8, height = 8)
}
# basic QC
# pct.mt > 10% -> cell stress or damage
# pct.ribo < 4% & pct.ribo > 45% -> poor cell health or technical artifacts
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.ribo", "pct.hb", "heat.shock.score"))
FeatureStatPlot(srat.merge.rna, group.by = "seurat_clusters", stat.by = c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.ribo", "pct.hb", "heat.shock.score"), stack = T, add_box = T)
# programmed cell death
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), ".pcd")])
# heat shocked cells
markers.de.sc = FindAllMarkers(srat.merge.rna, logfc.threshold = 0.1, test.use = "wilcox", slot = "data", min.pct = 0.3, only.pos = T)
markers.de.sc = markers.de.sc %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 50, order_by = avg_log2FC)
markers.de.sc$heat.shock = ifelse(markers.de.sc$gene %in% heat_shock.genes, T, F)
pheatmap(log2(table(markers.de.sc$cluster, markers.de.sc$heat.shock) + 1), scale = "none")

n = 36
pca_contribution(srat.merge.rna, n, reduction = "pca")

res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0)
srat.merge.rna = FindNeighbors(srat.merge.rna, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn")
for (i in res.list) {
  level = sort(as.numeric(levels(srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.merge.rna@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
# k.param: 20, the biggger, the clustering more mixed
srat.merge.rna = RunUMAP(srat.merge.rna, reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)

library(clustree)
clustree(srat.merge.rna, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)

# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.5"
srat.merge.rna[["seurat_clusters"]] = srat.merge.rna[[res_use]]
srat.merge.rna = SetIdent(srat.merge.rna, value = "seurat_clusters")
p = CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "seurat_clusters", label = T, label_insitu = T)
p

pdf(paste(tmp.dir1, "1basic_lighting.pdf", sep = "/"), width = 8*2.5, height = 8*2)
FeatureDimPlot(srat.merge.rna, features = c("nCount_RNA_log10", "nFeature_RNA_log10", "S.Score", "G2M.Score", 
                                            "MKI67", "TOP2A", "pct.mt", "pct.ribo", "pct.hb", "heat.shock.score"), reduction = "pca_umap", show_stat = F, ncol = 4)
FeatureDimPlot(srat.merge.rna, features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), ".pcd")], reduction = "pca_umap", show_stat = F, ncol = 5)
dev.off()

### annotation ###
# confirm big type
ctp1.geneset.list = list(
  "epi.markers" = c("EPCAM", "KRT8", "KRT18", "FXYD3", "PERP", "CDH1"),
  "mesen.markers" = c("COL6A2","MFAP4","DCN","FHL1","COL1A2","COL3A1"),
  "imm.markers" = c("CD37","CORO1A","LCP1","PTPRC","CD53","LAPTM5"),
  "endo.markers" = c("FLT1", "CDH5", "CLDN5", "EGFL7", "ESAM","IFI27"),
  "neu.markers" = c("CHGA","ASCL1","NNAT","STMN2","GRP","MPZ")
)
library(UCell)
srat.merge.rna = AddModuleScore_UCell(srat.merge.rna, features = ctp1.geneset.list, assay = "SCT", ncores = 4, name = "")

pdf(paste(tmp.dir1, "2big_type_lighting.pdf", sep = "/"), width = 8*1.5, height = 8)
p + FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = "mesen.markers")
p + FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = "epi.markers")
p + FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = "endo.markers")
p + FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = "neu.markers")
p + FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = "imm.markers")
dev.off()

######################
###### celltype1 #####
######################
srat.merge.rna$celltype1 = as.character(srat.merge.rna$seurat_clusters)
srat.merge.rna$celltype1[srat.merge.rna$seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
                                                               12, 13, 14, 15, 17)] = "Epithelium"
srat.merge.rna$celltype1[srat.merge.rna$seurat_clusters %in% c(11)] = "Mesenchyme"
srat.merge.rna$celltype1[srat.merge.rna$seurat_clusters %in% c(16)] = "Neuron"
srat.merge.rna$celltype1 = factor(as.character(srat.merge.rna$celltype1), levels = c("Mesenchyme", "Epithelium", "Neuron"))
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype1")

#####################
###### MP Score #####
#####################
library(UCell)
for (pkl in names(pkls)) {
  mp.df = read.csv(paste("/media/biogenger/D/Projects/CMY/Analysis/human_lung/ref_human_fetal_lung_wgcna", paste0("2cnet_", pkl), "2metaprogram/MP.hotspot.csv", sep = "/"))
  mp.df = mp.df %>% dplyr::filter(Module != "-1")
  mp.list = split(mp.df$X, mp.df$Module)
  names(mp.list) = paste0(pkl, ".MP", names(mp.list))
  srat.merge.rna = AddModuleScore_UCell(srat.merge.rna, features = mp.list, assay = "SCT", ncores = 4, name = "")
}

pdf(paste(tmp.dir2, "Epithelium.MP.pdf", sep = "/"), width = 8*2.5, height = 8*1.5)
FeatureDimPlot(srat.merge.rna, 
               features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Epithelium.MP")], 
               reduction = "pca_umap", show_stat = F, ncol = 5)
dev.off()

pdf(paste(tmp.dir2, "Stromal.MP.pdf", sep = "/"), width = 8*3, height = 8*1.5)
FeatureDimPlot(srat.merge.rna, 
               features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Stromal.MP")], 
               reduction = "pca_umap", show_stat = F, ncol = 6)
dev.off()

######################
###### LungMapv1 #####
######################
library(UCell)
srat.merge.rna = AddModuleScore_UCell(srat.merge.rna, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")

for (ctp in colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "_lungmap")]) {
  p = FeatureDimPlot(srat.merge.rna, features = ctp, reduction = "pca_umap", show_stat = F)
  ggsave(paste(tmp.dir2, "LunngMapV1.topDE", paste0(gsub("\\/", ".", ctp), ".png"), sep = "/"), p, width = 8, height = 8)
}

###########################
###### Classic Marker #####
###########################
library(UCell)
srat.merge.rna = AddModuleScore_UCell(srat.merge.rna, 
                                      features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                            "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                            "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
for (ctp in colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "_classic")]) {
  p = FeatureDimPlot(srat.merge.rna, features = ctp, reduction = "pca_umap", show_stat = F)
  ggsave(paste(tmp.dir2, "LungEpithelium.Classic", paste0(gsub("\\/", ".", ctp), ".png"), sep = "/"), p, width = 8, height = 8)
}

####################
###### ChatGPT #####
####################
# 首先要确定合适的分群
srat.merge.rna = SetIdent(srat.merge.rna, value = "seurat_clusters")
markers.de.sc = FindAllMarkers(srat.merge.rna, logfc.threshold = 1, test.use = "wilcox", slot = "data", min.pct = 0.5, only.pos = T)
markers.de.sc.valid = markers.de.sc %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)
markers.de.sc.classic = markers.de.sc.valid %>% dplyr::filter(gene %in% unlist(gene.sets.genger$human$LGEA[c("epithelium", "mesenchyme")]))
markers.de.sc.up.top = markers.de.sc.valid %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 15, order_by = avg_log2FC)
markers.de.sc.up.final = rbind(markers.de.sc.classic, markers.de.sc.up.top) %>% dplyr::distinct()
table(markers.de.sc.up.final$cluster)

library(openai)
library(GPTCelltype)
res = gptcelltype(split(markers.de.sc.up.final$gene, markers.de.sc.up.final$cluster), 
                  tissuename = "human lung airway epithelial", 
                  model = "gpt-4", mine_url = "https://api.xiaoai.plus")
srat.merge.rna@meta.data$gpt4_label = res[as.character(Idents(srat.merge.rna))]
srat.merge.rna@meta.data$gpt4_label = factor(srat.merge.rna@meta.data$gpt4_label, levels = str_sort(unique(srat.merge.rna@meta.data$gpt4_label), numeric = T))

pdf(paste(tmp.dir2, "GPTCelltype.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "gpt4_label", label = T, label_repel = T)
dev.off()

#######################
###### celltypist #####
#######################
srat.merge.rna = JoinLayers(srat.merge.rna, assay = "RNA")
srat.merge.rna = NormalizeData(srat.merge.rna, assay = "RNA")
seurat2scanpy(srat.merge.rna, ann.X = "RNA-data", ann.raw.X = NULL, h5ad_path = paste(tmp.dir1, "srat.merge.rna.logCPM.h5ad", sep = "/"))

pkls = list(
  "Stromal" = "/media/biogenger/D/Reference/celltypist/human_fetal_lung_ref/Stromal.pkl", 
  "Epithelium" = "/media/biogenger/D/Reference/celltypist/human_fetal_lung_ref/Epithelium.pkl", 
  "Endothelium" = "/media/biogenger/D/Reference/celltypist/human_fetal_lung_ref/Endothelium.pkl", 
  "Neuron" = "/media/biogenger/D/Reference/celltypist/human_fetal_lung_ref/Neuron.pkl", 
  "Immune" = "/media/biogenger/D/Reference/celltypist/human_fetal_lung_ref/Immune.pkl"
)
celltypist.list = list()
for (pkl in names(pkls)) {
  celltypist.list[[pkl]] = run_celltypist(model_path = pkls[[pkl]], h5ad_path = paste(tmp.dir1, "srat.merge.rna.logCPM.h5ad", sep = "/"), out_dir = paste(tmp.dir1, "celltypist", pkl, sep = "/"))
}
srat.merge.rna$celltypist_epithelium_label = celltypist.list$Epithelium$prob_match$majority_voting

pdf(paste(tmp.dir2, "CellTypist.major_vote.pdf", sep = "/"), width = 8, height = 8)
CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltypist_epithelium_label", label = T, label_repel = T)
dev.off()

for (pkl in names(pkls)) {
  celltypist.list[[pkl]] = list(
    "best_match" = read.csv(paste(tmp.dir1, "celltypist", pkl, "celltypist.best_match.csv", sep = "/"), row.names = 1),
    "prob_match" = read.csv(paste(tmp.dir1, "celltypist", pkl, "celltypist.prob_match.csv", sep = "/"), row.names = 1)
  )
  ps = plot_celltypist(srat.merge.rna, group.by = "seurat_clusters", celltypist.df = celltypist.list[[pkl]]$prob_match)[[4]]
  for (p in names(ps)) {
    ggsave(paste(tmp.dir2, "CellTypist.Abundance", pkl, paste0(gsub("\\/", ".", p), ".png"), sep = "/"), ps[[p]], width = 8, height = 8)
  }
}

##################
###### GO-BP #####
##################
db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")
# hyper-test
srat.merge.rna = SetIdent(srat.merge.rna, value = "seurat_clusters")
markers.de.tmp = FindAllMarkers(srat.merge.rna, logfc.threshold = 1, test.use = "wilcox", slot = "data", min.pct = 0.5, only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
table(markers.de.tmp$cluster)
df.enrich = enricher_hyper_genger(gene.list = split(markers.de.tmp$gene, markers.de.tmp$cluster), minGSSize = 10, maxGSSize = 500, 
                                  db.df = db.merge.df.list$db.merge.df.final, pathway.id = "pathway.id", pathway.name = "pathway.name", pathway.gene = "pathway.gene", 
                                  database = c("MsigDB"), collection = c("GO_BP"), specy = "human")
df.enrich$enrich.score = -log10(df.enrich$p_adj)
df.enrich$generatio = df.enrich$hit.gene.size / df.enrich$cluster.gene.size
df.enrich$bgRatio = df.enrich$pathway.gene.size / df.enrich$all.pathway.gene.size
df.enrich$enrichment.foldchange = df.enrich$generatio / df.enrich$bgRatio
df.plot = enricher_plot_genger(plot_df = df.enrich, x.col1 = "list.name", x.col2 = "collection", y.name.col = "pathway.name", y.name.len = 60, 
                               each.sig.n = 10, top = T,
                               color.by = "enrich.score", color.by.legend.name = "-log10(p.adj）", 
                               size.by = "enrichment.foldchange", 
                               plot.path = paste(tmp.dir1, "GOBP_seurat_clusters", sep = "/"), 
                               plot.compare.width = 25, plot.compare.height = 30, 
                               plot.ind.width = 10, plot.ind.height = 8)
write.csv(df.enrich, paste(tmp.dir1, "GOBP_seurat_clusters", "enrichment.csv", sep = "/"), row.names = F)

######################
###### celltype2 #####
######################
# pre-look
p.whole = CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "seurat_clusters", label = T, label_insitu = T)
p.whole
#############
# 2: SCGB1A1+, SCGB3A2+, SCGB3A1+ secretory (MUC16+) (11 pcw) (Cell, fetal lung) (RAS?)
# 6: SCGB1A1+, SCGB3A2LO/-, SCGB3A1+ secretory (MUC16+) (CYP2F1+, MUC4+, KRT4+) (15 pcw) (Cell, fetal lung)
# 17: goblet (MUC13+, TSPAN8+)
# 1: Tissue Homeostasis (SPP1;ALB;SERPINA3;ITGB1;CCN2;CXADR) 肽, 离子, 激素等调节; 上皮增殖
# 2: Response To Chemokine (CXCL1;EDN1;CXCL3;CCL20;CXCL2;CXCL8;ZC3H12A;CCL2;CXCL6) 免疫招募; 上皮凋亡
# 3: Endodermal Cell Differentiation (FN1;LAMA3) 代谢调节; 上皮增殖 
# 4: Neuron Death (CCL2;TUBA1A;SOD2;NCOA7) 免疫招募
# 8: Antiviral Innate Immune Response (IFIT1;MX1) 干扰素反应; 免疫反应
# 13: Regulation Of Mirna Transcription (ETS1;TGFB2;SOX9;KLF4;HIF1A) 上皮增殖和迁移; TGFB调控
# 9: suprabasal (more differentiated than basal which are under them)
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial$suprabasal$up, assay = "MAGIC_SCT")
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = lungmapv1.de.list$`SMG_Basal/Duct`, assay = "MAGIC_SCT")
srat.merge.rna = SetIdent(srat.merge.rna, value = "seurat_clusters")
tmp = FindMarkers(srat.merge.rna, assay = "SCT", ident.1 = 0, min.pct = 0.5, only.pos = T) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC))
rownames(tmp)[1:12]
FeatureDimPlot(srat.merge.rna, reduction = "pca_umap", features = rownames(tmp)[1:12], assay = "SCT", ncol = 3)

plot.cluster_similarity(srat.merge.rna, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

srat.merge.rna$celltype2 = "genger"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(6)]] = "Secretory1"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(5)]] = "Secretory2"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(17)]] = "Goblet"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(1)]] = "A1"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(2)]] = "A2"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(3)]] = "A3"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(4)]] = "A4"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(8)]] = "A5"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(13)]] = "A6"
srat.merge.rna$celltype2[colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(9)]] = "Suprabasal"
p.whole + CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "celltype2", label = T, label_insitu = T)
#############
# 15: ciliated
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(15)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 42
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.1"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c(paste0(c("Ciliated", "Deuterosomal"), "_lungmap"),
                                                              paste0(c("ciliated", "deuterosomal"), "_classic")), ncol = 2)
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(1)]] = "Deuterosomal"
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 2, 3)]] = "Ciliated"

# 16: pnec
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(16)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 42
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.7"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c(paste0(c("PNEC"), "_lungmap"),
                                                              paste0(c("PNEC"), "_classic")), ncol = 2)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c("Epithelium.MP3", "Epithelium.MP4", "Epithelium.MP6", "Epithelium.MP11", "Epithelium.MP12"), ncol = 3)
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1, 2)]] = "Neuroendocrine"

# 11: airway.FB
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(11)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 42
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.5"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$mesenchyme[c("Lipofibroblast", "Pericyte", "secondary_crest_myofibroblast_cell", 
                                                                                    "proliferative_mesenchymal_progenitor", "Chondrocytes", "Mesothelium", "Adventitial")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c(paste0(c("MEC", "SCMF", "Chondrocyte"), "_lungmap"),
                                                              paste0(c("myoepithelial", "secondary_crest_myofibroblast_cell", "Chondrocytes"), "_classic")), ncol = 3)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Stromal.MP")], ncol = 6)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = gene.sets.genger$human$LGEA$mesenchyme$Chondrocytes$up)
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1, 2, 3, 4)]] = "SCMF"

# 12: myoepith
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(12)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 42
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.5"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c(paste0(c("MEC", "Basal"), "_lungmap"),
                                                              paste0(c("myoepithelial", "basal"), "_classic")), ncol = 2)
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1, 2, 3, 4)]] = "MEC"

# 14: proliferating epith
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(14)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 43
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.3"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "orig.ident", cells.highlight = colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(3)])

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c("SOX9", "SOX2", "MKI67", "TOP2A",
                                                              paste0(c("MEC", "Basal"), "_lungmap"), paste0(c("myoepithelial", "basal"), "_classic")), ncol = 3)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Epithelium.MP")], ncol = 5)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c("DTL", "WDR76", "TESC"), assay = "MAGIC_SCT")
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(3)]] = "Proliferating Stalk"
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1)]] = "Basal"
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(2)]] = "Proliferating DTL+ Stalk" # Cell Research

# 0 7 10 Basal(from above): proliferating epith
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(0, 7, 10) | srat.merge.rna$celltype2 == "Basal"])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 43
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 20, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 20, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.3"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)
plot.cluster_similarity(srat.tmp, assay = "SCT", slot = "data", cluster_col = "seurat_clusters", method = "pearson", size = 10)

CellDimPlot(srat.merge.rna, reduction = "pca_umap", group.by = "orig.ident", cells.highlight = colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(14)])

srat.tmp = AddModuleScore_UCell(srat.tmp, features = lungmapv1.de.list, assay = "SCT", ncores = 4, name = "_lungmap")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$epithelium$airway_epithelial[c("basal", "ciliated", "goblet", "secretory", "ionocyte", "PNEC", "myoepithelial", 
                                                                                                      "tuft", "serous", "mucous", "BASC", "terminal_ciliated_duct_basal_cells", 
                                                                                                      "RAS", "deuterosomal", "suprabasal")], assay = "SCT", ncores = 4, name = "_classic")
srat.tmp = AddModuleScore_UCell(srat.tmp, 
                                features = gene.sets.genger$human$LGEA$mesenchyme[c("Lipofibroblast", "Pericyte", "secondary_crest_myofibroblast_cell", 
                                                                                    "proliferative_mesenchymal_progenitor", "Chondrocytes", "Mesothelium", "Adventitial")], assay = "SCT", ncores = 4, name = "_classic")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c("SOX9", "SOX2", "MKI67", "TOP2A", "DTL", "WDR76", "TESC"), ncol = 3, assay = "MAGIC_SCT")
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c(paste0(c("MEC", "SCMF", "Basal"), "_lungmap"), 
                                                              paste0(c("myoepithelial", "proliferative_mesenchymal_progenitor", "basal"), "_classic")), ncol = 3)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Epithelium.MP")], ncol = 5)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = colnames(srat.merge.rna@meta.data)[str_detect(colnames(srat.merge.rna@meta.data), "Stromal.MP")], ncol = 5)
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = c("nCount_RNA_log10", "nFeature_RNA_log10", "S.Score", "G2M.Score", 
                                                              "pct.mt", "pct.ribo", "pct.hb", "heat.shock.score"))
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = heat_shock.genes, assay = "MAGIC_SCT")

srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
tmp = FindMarkers(srat.tmp, assay = "SCT", ident.1 = 5, min.pct = 0.5, only.pos = T, recorrect_umi = FALSE) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC))
rownames(tmp)[1:12]
FeatureDimPlot(srat.tmp, reduction = "pca_umap", features = rownames(tmp)[1:12], assay = "SCT", ncol = 3)

CellStatPlot(srat.tmp, stat.by = "orig.ident", group.by = "seurat_clusters")

srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(4)]] = "Proliferating DTL+ Stalk" # Cell Research
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1, 2, 3, 6, 7)]] = "Basal"
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(5, 8, 9, 10, 11, 12, 13, 14)]] = "Basal" # low quality

# 1 17: goblet -> 主要为了增加G30的Goblet数量
srat.tmp = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$seurat_clusters %in% c(1, 17)])
srat.tmp@meta.data = srat.tmp@meta.data %>% dplyr::select(!starts_with("SCT_snn_res."))
res.list = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0)
srat.tmp = RunPCA(srat.tmp, assay = "SCT", npcs = 50)
n = 42
pca_contribution(srat.tmp, n, reduction = "pca")
srat.tmp = FindNeighbors(srat.tmp, reduction = "pca", dims = 1:n, nn.method = "annoy", annoy.metric = "cosine", k.param = 10, graph.name = c("SCT_nn", "SCT_snn")) %>% 
  FindClusters(algorithm = 1, resolution = res.list, graph.name = "SCT_snn") %>% 
  RunUMAP(reduction = "pca", reduction.name = "pca_umap", dims = 1:n, umap.method = "uwot", metric = "cosine", min.dist = 0.3, n.neighbors = 10, return.model = T, seed.use = 666)
for (i in res.list) {
  level = sort(as.numeric(levels(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]])))
  srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]] = factor(srat.tmp@meta.data[[paste0("SCT_snn_res.", i)]], levels = level)
}
library(clustree)
clustree(srat.tmp, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.0.3"
srat.tmp[["seurat_clusters"]] = srat.tmp[[res_use]]
srat.tmp = SetIdent(srat.tmp, value = "seurat_clusters")
CellDimPlot(srat.tmp, reduction = "pca_umap", group.by = c("seurat_clusters", "celltypist_epithelium_label", "gpt4_label"), label = T, label_insitu = T, ncol = 1)

table(srat.tmp$orig.ident, srat.tmp$seurat_clusters)

srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(5)]] = "Goblet"
srat.merge.rna$celltype2[colnames(srat.tmp)[srat.tmp$seurat_clusters %in% c(0, 1, 2, 3, 4, 6, 7)]] = "A1"

###################
###### Output #####
###################
srat.merge.rna$celltype2 = factor(as.character(srat.merge.rna$celltype2), levels = str_sort(unique(srat.merge.rna$celltype2), numeric = T))
# final annotation
srat.merge.rna$celltype2_final = as.character(srat.merge.rna$celltype2)
srat.merge.rna$celltype2_final[srat.merge.rna$celltype2 %in% c("A1", "A2", "A3", "A4", "A5")] = "Immature Airway Epithelium1"
srat.merge.rna$celltype2_final[srat.merge.rna$celltype2 %in% c("A6")] = "Immature Airway Epithelium2"
srat.merge.rna$celltype2_final[srat.merge.rna$celltype2 %in% c("SCMF")] = "Myofibroblast"
srat.merge.rna$celltype2_final = factor(srat.merge.rna$celltype2_final, levels = c("Immature Airway Epithelium1", "Immature Airway Epithelium2",
                                                                                   "Proliferating Stalk", "Proliferating DTL+ Stalk", 
                                                                                   "Basal", "Suprabasal", "MEC", 
                                                                                   "Deuterosomal", "Ciliated", 
                                                                                   "Goblet", 
                                                                                   "Secretory1", "Secretory2",
                                                                                   "Neuroendocrine", 
                                                                                   "Myofibroblast"))
write.csv(srat.merge.rna@meta.data, paste(tmp.dir1, "metadata.csv", sep = "/"))
qsave(srat.merge.rna, paste(tmp.dir1, "srat.anno.qs", sep = "/"))
srat.merge.rna = qread(paste(tmp.dir1, "srat.anno.qs", sep = "/"))

pdf(paste(tmp.dir1, "3annotation.umap.pdf", sep = "/"), width = 8*2, height = 8)
CellDimPlot(srat.merge.rna, group.by = "orig.ident", reduction = "pca_umap", palcolor = colors.organoid.airway.sample, label = F, theme_use = "theme_blank") + 
  CellDimPlot(srat.merge.rna, group.by = "celltype1", reduction = "pca_umap", palette = "Set1", label = F, theme_use = "theme_blank")
CellDimPlot(srat.merge.rna, group.by = c("seurat_clusters", "celltype2_final"), reduction = "pca_umap", palette = "Set1", label = T, label_repel = T)
dev.off()

#########################
###### MetaNeighbor #####
#########################
library(MetaNeighbor)
mn_data = SummarizedExperiment(assays = list(counts = srat.merge.rna[["RNA"]]$counts), colData = srat.merge.rna@meta.data)
var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$orig.ident)
celltype_NV = MetaNeighborUS(dat = mn_data, study_id = mn_data$orig.ident, cell_type = mn_data$celltype2_final, var_genes = var_genes, fast_version = T)
top_hits = topHits(cell_NV = celltype_NV, dat = mn_data, study_id = mn_data$orig.ident, cell_type = mn_data$celltype2_final, threshold = 0.9)
#分别自定义列注释和行注释函数：
column_ha = HeatmapAnnotation(celltype = str_split(colnames(celltype_NV), "\\|", simplify = T)[,2], 
                              study = str_split(colnames(celltype_NV), "\\|", simplify = T)[,1], 
                              col = list(study = colors.organoid.airway.sample, 
                                         celltype = palette_scp(levels(mn_data$celltype2_final), palette = "Set1")), which = "col")
row_ha = HeatmapAnnotation(study = str_split(rownames(celltype_NV), "\\|", simplify = T)[,1],
                           celltype = str_split(rownames(celltype_NV), "\\|", simplify = T)[,2],
                           col = list(study = colors.organoid.airway.sample,
                                      celltype = palette_scp(levels(mn_data$celltype2_final), palette = "Set1")), which = "row")
pdf(paste(tmp.dir1, "4annotation.metaneighbor.pdf", sep = "/"), width = 8*2.5, height = 8*1.5)
Heatmap(celltype_NV, name = "celltype similarity across orig.ident", 
        col = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)),
        clustering_method_rows = "complete", clustering_method_columns = "complete",
        top_annotation = column_ha, right_annotation = row_ha, 
        row_dend_side = "right", column_dend_side = "top",
        show_row_names = T, row_names_side = "left", 
        show_column_names = T, column_names_side = "bottom",
        width = unit(18, "cm"), height = unit(18, "cm"))
# spearman cor
plot.cluster_similarity(srat.merge.rna, assay = "SCT", slot = "data", cluster_col = "celltype2_final", method = "pearson", size = 20)
dev.off()

##########################
##### DEG pseudobulk #####
##########################
# pseudobulking DE
srat.merge.rna.bulk = AggregateExpression(srat.merge.rna, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "orig.ident", "celltype2_final"))

srat.merge.rna.bulk = SetIdent(srat.merge.rna.bulk, value = "celltype2_final")
markers.de.bulk = FindAllMarkers(srat.merge.rna.bulk, logfc.threshold = 0.3, test.use = "wilcox", slot = "data", min.pct = 0.5, only.pos = T)
# write.csv(markers.de.bulk, file = paste(tmp.dir1, "markers.de.bulk.csv", sep = "/"))
markers.de.bulk = read.csv(paste(tmp.dir1, "markers.de.bulk.csv", sep = "/"), row.names = 1)
markers.de.bulk$cluster = factor(markers.de.bulk$cluster, levels = levels(srat.merge.rna$celltype2_final))

markers.de.bulk.valid = markers.de.bulk %>% dplyr::filter(avg_log2FC > 0.5 & p_val < 0.05) %>% pureCluster(cluster.col = "cluster", feature.col = "gene", value.col = "avg_log2FC")
markers.de.bulk.up.top = markers.de.bulk.valid %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 3, order_by = avg_log2FC)

##########################
##### DEG singlecell #####
##########################
srat.merge.rna = SetIdent(srat.merge.rna, value = "celltype2_final")
markers.de.sc = FindAllMarkers(srat.merge.rna, logfc.threshold = 0.3, test.use = "wilcox", slot = "data", min.pct = 0.5, only.pos = T)
# write.csv(markers.de.sc, file = paste(tmp.dir1, "markers.de.sc.csv", sep = "/"))
markers.de.sc = read.csv(paste(tmp.dir1, "markers.de.sc.csv", sep = "/"), row.names = 1)
markers.de.sc$cluster = factor(markers.de.sc$cluster, levels = levels(srat.merge.rna$celltype2_final))

markers.de.sc.valid = markers.de.sc %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% pureCluster()
markers.de.sc.up.top = markers.de.sc.valid %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 2, order_by = avg_log2FC)

# 记得修改为classic markers
ht1 = GroupHeatmap(srat.merge.rna, group.by = "celltype2_final", assay = "RNA", slot = "data", 
                   features = unique(c(markers.de.sc.valid$gene, unlist(plot.marker.list))), features_label = unlist(plot.marker.list),
                   cluster_rows = F, cluster_columns = F, group_palette = "Set1",
                   heatmap_palcolor = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)), 
                   height = 7, width = 8)
ht2 = Doheatmap_genger(srat.merge.rna, assay = "SCT", slot = "counts", sort.by = "celltype2_final", scale.split.by = "orig.ident", vars.to.regress = NULL,
                       split.by = "celltype2_final", column_gap = unit(0.5, "mm"), 
                       marker.genes = unique(c(markers.de.sc.valid$gene, unlist(plot.marker.list))), markers.highlight = unlist(plot.marker.list), highlight.line.color = "black", highlight.line.lwd = 0.3, 
                       font.col.row = "black", font.size.row = 8, col.name.rot = 90, col.name.size = 0, 
                       annotation.vars = c("nCount_RNA_log10", "nFeature_RNA_log10", "pct.mt", "pct.ribo", "pct.hb", "orig.ident", "celltype2_final"), 
                       annotation.cols = list("nCount_RNA_log10" =  viridis(3, option = "D"), 
                                              "nFeature_RNA_log10" = viridis(3, option = "D"), 
                                              "pct.mt" =  viridis(3, option = "F"), 
                                              "pct.ribo" = viridis(3, option = "F"), 
                                              "pct.hb" = viridis(3, option = "F"), 
                                              "orig.ident" = palette_scp(x = levels(srat.merge.rna$orig.ident), palette = "Paired", type = "discrete"),
                                              "celltype2_final" = palette_scp(x = levels(srat.merge.rna$celltype2_final), palette = "Set1", type = "discrete")),
                       heat.min = NULL, heat.max = NULL, heat.scale = NULL, heat.cols = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100), 
                       plot.mean = F)
pdf(paste(tmp.dir1, "5top.DEGs.heat.pdf", sep = "/"), width = 8*1.5, height = 8)
ht1$plot
draw(ht2)
dev.off()

plot.marker.list = list(
  "Immature Airway Epithelium1" = c("ITIH2", "FXYD2"), # Response To Steroid Hormone;Cytokine Mediated Signaling Pathway;Positive Regulation Of Epithelial Cell Migration
  "Immature Airway Epithelium2" = c("ANKRD1", "FILIP1L"), # Positive Regulation Of Epithelial Cell Migration;Epithelial Cell Development;Regulation Of Rna Splicing
  "Proliferating Stalk" = c("MKI67", "TOP2A", "KIF20A"), # Nuclear Chromosome Segregation;Mitotic Nuclear Division;Cell Cycle G2 M Phase Transition
  "Proliferating DTL+ Stalk" = c("DTL", "PCLAF"), # Dna Replication;Dna Conformation Change;Deoxyribonucleotide Metabolic Process
  "Basal" = c("TP63", "IGFBP4", "COL17A1"), # Insulin Like Growth Factor Receptor Signaling Pathway;Response To Oxygen Levels;Response To Mechanical Stimulus
  "Suprabasal" = c("KRT5", "NOTCH3", "KRT6C"), # Keratinocyte Differentiation;Intermediate Filament Organization;Positive Regulation Of Cell Growth
  "MEC" = c("KRT14", "CXCL14"), # Epidermis Development;Muscle Cell Migration;Epithelial Cell Proliferation
  "Deuterosomal" = c("DEUP1", "CDC20B"), # Microtubule Based Movement;Microtubule Based Transport;Axoneme Assembly
  "Ciliated" = c("FOXJ1", "TMEM190"), # Cilium Organization;Cilium Movement;Cilium Or Flagellum Dependent Cell Motility
  "Goblet" = c("MUC5AC", "TFF1", "MUC13", "TSPAN8"), # Protein Localization To Plasma Membrane;Negative Regulation Of Cell Adhesion;Response To Oxygen Levels
  "Secretory1" = c("SCGB1A1", "SCGB3A1", "SCGB3A2"), # Sterol Biosynthetic Process;Positive Regulation Of Macrophage Chemotaxis;Negative Regulation Of Peptidase Activity
  "Secretory2" = c("CYP2F1", "KRT4", "ALDH1A3"), # Granulocyte Chemotaxis;Negative Regulation Of Hydrolase Activity;Leukocyte Migration Involved In Inflammatory Response
  "Neuroendocrine" = c("CALCA", "GRP", "ASCL1", "CHGA"), # Oxidative Phosphorylation;Atp Metabolic Process;Neuropeptide Signaling Pathway
  "Myofibroblast" = c("ACTA2", "TGFBI", "MGP") # Collagen Fibril Organization;Extracellular Matrix Assembly;Transforming Growth Factor Beta Production
)

# anno.df = db.merge.df.list$db.merge.df.final %>% dplyr::filter(species == "human" & database == "MsigDB" & collection %in% "curated")
# pathway.list = split(anno.df$pathway.gene, anno.df$pathway.name)
# pathway.list[["Reactome Signaling By Notch"]]

plot.marker.df = data.frame()
for (ctp in names(plot.marker.list)) {
  plot.marker.df = rbind(plot.marker.df, data.frame(cluster = ctp, gene = plot.marker.list[[ctp]]))
}
plot.marker.df$cluster = factor(plot.marker.df$cluster, levels = levels(srat.merge.rna$celltype2_final))
library(scRNAtoolVis)
p = jjDotPlot(srat.merge.rna, assay = "SCT", slot = "data", id = "celltype2_final", markerGene = plot.marker.df, 
              anno = T, scale = T, col.min = -2, col.max = 2, xtree = F, ytree = F, plot.margin = c(5, 1, 1, 1))
pdf(paste(tmp.dir1, "5classic_top.DEGs.dot.pdf", sep = "/"), width = 8*2, height = 8*1.5)
print(p)
dev.off()

######################
##### enrichment #####
######################
db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")

# hyper-test
markers.de.sc = read.csv(paste(tmp.dir1, "markers.de.sc.csv", sep = "/"), row.names = 1)
markers.de.sc$cluster = factor(markers.de.sc$cluster, levels = levels(srat.merge.rna$celltype2_final))
markers.de.sc.valid = markers.de.sc %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)
df.enrich = enricher_hyper_genger(gene.list = split(markers.de.sc.valid$gene, markers.de.sc.valid$cluster), minGSSize = 10, maxGSSize = 500, 
                                  db.df = db.merge.df.list$db.merge.df.final, pathway.id = "pathway.id", pathway.name = "pathway.name", pathway.gene = "pathway.gene", 
                                  database = c("MsigDB"), collection = c("curated", "KEGG", "wiki", "GO_BP", "hallmark"), specy = "human")
df.enrich$enrich.score = -log10(df.enrich$p_adj)
df.enrich$generatio = df.enrich$hit.gene.size / df.enrich$cluster.gene.size
df.enrich$bgRatio = df.enrich$pathway.gene.size / df.enrich$all.pathway.gene.size
df.enrich$enrichment.foldchange = df.enrich$generatio / df.enrich$bgRatio
df.plot = enricher_plot_genger(plot_df = df.enrich, x.col1 = "list.name", x.col2 = "collection", y.name.col = "pathway.name", y.name.len = 60, 
                               each.sig.n = 10, top = T,
                               color.by = "enrich.score", color.by.legend.name = "-log10(p.adj）", 
                               size.by = "enrichment.foldchange", 
                               plot.path = paste(tmp.dir1, "enrichment", "hyper-test", sep = "/"), 
                               plot.compare.width = 25, plot.compare.height = 30, 
                               plot.ind.width = 10, plot.ind.height = 8)
write.csv(df.enrich, paste(tmp.dir1, "enrichment", "hyper-test", "enrichment.csv", sep = "/"), row.names = F)

# GSEA
srat.merge.rna = SetIdent(srat.merge.rna, value = "celltype2_final")
names(srat.merge.rna[["RNA"]]@layers) # data是否合并?
df.gsea = enricher_gsea_genger(expr.mat = GetAssayData(srat.merge.rna, assay = "RNA", slot = "data"), cell2label = Idents(srat.merge.rna),  thresh.min = 0, logged = T, minGSSize = 10, maxGSSize = 500, 
                               db.df = db.merge.df.list$db.merge.df.final, pathway.id = "pathway.id", pathway.name = "pathway.name", pathway.gene = "pathway.gene", 
                               database = c("MsigDB"), collection = c("curated", "KEGG", "wiki", "GO_BP", "hallmark"), specy = "human")$mine.df
df.gsea$enrich.score = -log10(df.gsea$p_adj)
df.plot = enricher_plot_genger(plot_df = df.gsea, x.col1 = "list.name", x.col2 = "collection", y.name.col = "pathway.name", y.name.len = 60, 
                               each.sig.n = 10, top = T,
                               color.by = "NES", color.by.legend.name = "-log10(p.adj）", 
                               size.by = "enrich.score", 
                               plot.path = paste(tmp.dir1, "enrichment", "GSEA", sep = "/"), 
                               plot.compare.width = 25, plot.compare.height = 30, 
                               plot.ind.width = 10, plot.ind.height = 8)
write.csv(df.enrich, paste(tmp.dir1, "enrichment", "GSEA", "enrichment.csv", sep = "/"), row.names = F)

pdf("/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/2merge_annotation/big.ctp.umap.pdf", width = 8*2, height = 8)
FeatureDimPlot(
  srt = srat.merge.rna, features = c("epi.markers", "mesen.markers", "neu.markers"),
  label = F,
  reduction = "pca_umap", theme_use = "theme_blank"
)
dev.off()
