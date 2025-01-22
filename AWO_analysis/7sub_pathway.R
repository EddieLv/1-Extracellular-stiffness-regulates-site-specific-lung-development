source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway"

tmp.dir1 = paste(work.dir, "7sub_pathway", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
srat.merge.rna$hard.num = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G30"] = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G60"] = 1
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G90"] = 2

srat.sub = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[!srat.merge.rna$celltype2_final %in% c("Immature Airway Epithelium1", 
                                                                                                          "Immature Airway Epithelium2")])
srat.sub$celltype2_final = droplevels(srat.sub$celltype2_final)
srat.sub = SetIdent(srat.sub, value = "celltype2_final")
########################
# decoupleR #
########################
library(decoupleR)
df.tmp = read.table("/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/7sub_pathway/progeny-net.txt", header = T) %>% as.tibble() %>% OmnipathR::pivot_annotations()
top = 500
net = df.tmp %>% dplyr::distinct(pathway, genesymbol, .keep_all = TRUE) %>% 
  dplyr::mutate(weight = as.double(weight), p_value = as.double(p_value)) %>% 
  dplyr::select(genesymbol, p_value, pathway, weight) %>% dplyr::group_by(pathway) %>% 
  dplyr::group_split() %>% 
  purrr::map(function(df) {df %>% dplyr::arrange(p_value) %>% head(top)}) %>% 
  dplyr::bind_rows() %>% 
  dplyr::select(pathway, genesymbol, weight, p_value) %>% 
  rlang::set_names(c("source", "target", "weight", "p_value"))
mat = as.matrix(srat.sub[["SCT"]]$data)
# Run mlm
acts = run_mlm(mat = mat, net = net, .source = "source", .target = "target", .mor = "weight", minsize = 5)
acts
# Extract mlm and store it in pathwaysmlm in data
srat.sub[["pathwaysmlm"]] = acts %>%
  pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") %>%
  column_to_rownames("source") %>%
  Seurat::CreateAssayObject(.)
# Change assay
DefaultAssay(srat.sub) = "pathwaysmlm"
# Scale the data
srat.sub = ScaleData(srat.sub)
srat.sub[["pathwaysmlm"]]$data = srat.sub[["pathwaysmlm"]]$scale.data
FeatureDimPlot(srat.sub, assay = "pathwaysmlm", slot = "scale.data", features = "TGFb", reduction = "pca_umap")
# Extract activities from object as a long dataframe
df = t(as.matrix(srat.sub[["pathwaysmlm"]]$data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Idents(srat.sub)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))
# Transform to wide matrix
top_acts_mat = df %>%
  pivot_wider(id_cols = "cluster", names_from = "source", values_from = "mean") %>%
  column_to_rownames("cluster") %>%
  as.matrix()
# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
my_breaks = c(seq(-2, 0, length.out = ceiling(palette_length/2) + 1), seq(0.05, 2, length.out = floor(palette_length/2)))
# Plot

pdf(paste(tmp.dir1, "progeny.pdf", sep = "/"), width = 8, height = 8)
pheatmap::pheatmap(top_acts_mat, border_color = NA, color = my_color, breaks = my_breaks, cellwidth = 25, cellheight = 25, main = "decoupleR")
CellDimPlot(srat.sub, group.by = "celltype2_final", reduction = "pca_umap", label = T, label_insitu = T, palette = "Set1")
FeatureDimPlot(srat.sub, assay = "pathwaysmlm", slot = "scale.data", features = c("Hypoxia", "TNFa", "TGFb", "Estrogen"), reduction = "pca_umap", ncol = 2)
dev.off()

########################
# pathway #
########################
db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")
anno.df = db.merge.df.list$db.merge.df.final %>% dplyr::filter(species == "human" & database == "MsigDB" & collection %in% "GO_BP")
pathway.list = split(anno.df$pathway.gene, anno.df$pathway.name)

### UCell ###
library(UCell)
pathway.list.tmp = pathway.list
names(pathway.list.tmp) = paste0("GOBP-", names(pathway.list))
srat.sub = AddModuleScore_UCell(srat.sub, 
                                features = pathway.list.tmp[paste0("GOBP-", 
                                                                   c("Hypoxia Inducible Factor 1alpha Signaling Pathway",
                                                                     "Positive Regulation Of Tumor Necrosis Factor Mediated Signaling Pathway",
                                                                     "Positive Regulation Of Transforming Growth Factor Beta Production",
                                                                     "Intracellular Estrogen Receptor Signaling Pathway"))], 
                                assay = "SCT", ncores = 8, name = "", maxRank = 2000)

FeatureDimPlot(srat.sub, assay = "pathwaysmlm", slot = "scale.data", features = paste0("GOBP-", 
                                                                                       c("Hypoxia Inducible Factor 1alpha Signaling Pathway",
                                                                                         "Positive Regulation Of Tumor Necrosis Factor Mediated Signaling Pathway",
                                                                                         "Positive Regulation Of Transforming Growth Factor Beta Production",
                                                                                         "Intracellular Estrogen Receptor Signaling Pathway")), reduction = "pca_umap", ncol = 2)


