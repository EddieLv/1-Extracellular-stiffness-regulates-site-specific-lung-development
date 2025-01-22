source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/flow_sep/airway"

tmp.dir1 = paste(work.dir, "4time_predict", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge1 = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
DefaultAssay(srat.merge1) = "RNA"
srat.merge1 = DietSeurat(srat.merge1, assays = "RNA")
srat.merge1$condition = srat.merge1$orig.ident
srat.bulk1 = AggregateExpression(srat.merge1, assays = "RNA", return.seurat = T, group.by = c("condition"))
srat.bulk1$dataset = "organoid"

srat.merge2 = qread("/media/biogenger/D/Projects/CMY/Analysis/human_lung/ref_human_fetal_lung_wgcna/2cnet_Epithelium/3cytotrace2/srat.reanno.mp.ctyotrace2.qs")
DefaultAssay(srat.merge2) = "RNA"
srat.merge2 = DietSeurat(srat.merge2, assays = "RNA")
srat.merge2$condition = paste0(srat.merge2$dataset, "-", srat.merge2$timepoint)
srat.bulk2 = AggregateExpression(srat.merge2, assays = "RNA", return.seurat = T, group.by = c("condition"))
srat.bulk2$dataset = str_split(srat.bulk2$condition, "-", simplify = T)[, 2]

# merge
srat.bulk = merge(srat.bulk1, srat.bulk2)
srat.bulk = JoinLayers(srat.bulk, assay = "RNA")
srat.bulk$timepoint = str_split(srat.bulk$condition, "-", simplify = T)[, 4]
srat.bulk$timepoint = factor(srat.bulk$timepoint, levels = str_sort(unique(srat.bulk$timepoint), numeric = T))
srat.bulk$stage = "genger"
srat.bulk$stage[srat.bulk$condition %in% samples.rna] = srat.bulk$condition[srat.bulk$condition %in% samples.rna]
srat.bulk$stage[srat.bulk$timepoint %in% paste0(seq(2, 7, 0.01), "w")] = "Embryonic"
srat.bulk$stage[srat.bulk$timepoint %in% paste0(seq(7, 17, 0.01), "w")] = "Pseudoglandular"
srat.bulk$stage[srat.bulk$timepoint %in% paste0(seq(17, 26, 0.01), "w")] = "Canalicular"
srat.bulk$stage[srat.bulk$timepoint %in% paste0(seq(26, 38, 0.01), "w")] = "Saccular"
srat.bulk$stage[srat.bulk$timepoint %in% paste0(seq(38, 99, 0.01), "w")] = "Alveolar"
srat.bulk$stage = factor(srat.bulk$stage, levels = c(samples.rna, "Embryonic", "Pseudoglandular", "Canalicular"))

srat.bulk = srat.bulk %>% 
  NormalizeData(assay = "RNA", normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(assay = "RNA", nfeatures = 1000) %>% 
  ScaleData(assay = "RNA", features = VariableFeatures(srat.bulk)) %>% 
  # RunPCA(assay = "RNA", features = VariableFeatures(srat.bulk), npcs = 20) %>% 
  RunUMAP(assay = "RNA", features = VariableFeatures(srat.bulk), n.components = 2, n.neighbors = 10, min.dist = 0.1)

CellDimPlot(srat.bulk, reduction = "pca", group.by = "dataset", palette = "Set1", pt.size = 2)

library(harmony)
srat.bulk = RunHarmony(srat.bulk, group.by.vars = "dataset", reduction.use = "umap")

CellDimPlot(srat.bulk, reduction = "harmony", group.by = c("stage", "dataset"), palette = "Set1", pt.size = 2,
            cells.highlight = colnames(srat.bulk)[srat.bulk$orig.ident %in% samples.rna], sizes.highlight = 3)

# similarity
library(stats)
library(ComplexHeatmap)
expr.mat = as.matrix(srat.bulk[["RNA"]]$data)
expr.mat.ave = sapply(split(colnames(srat.bulk), srat.bulk$stage)[c("Embryonic", "Pseudoglandular", "Canalicular")], function(cells) {rowMeans(expr.mat[ , cells])})
expr.mat.final = cbind(expr.mat.ave, expr.mat[, samples.rna])
cormatrix = cor(expr.mat.final, method = "spearman")
cordist = as.dist(1 - cormatrix, diag = T, upper = T)
corclus = hclust(d = cordist, method = "ward.D2")
hmap = Heatmap(matrix = cormatrix, cluster_rows = corclus, cluster_columns = corclus, col = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)))
hmap

plot.cluster_similarity(srat.bulk, assay = "RNA", slot = "data", cluster_col = "stage", method = "pearson", size = 10)


