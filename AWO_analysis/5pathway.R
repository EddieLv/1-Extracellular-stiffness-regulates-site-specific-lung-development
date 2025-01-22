source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("R51d40G30", "R51d40G60", "R51d40G90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/bronchi_old/airway"

tmp.dir1 = paste(work.dir, "5pathway", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
srat.merge.rna$hard.num = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G30"] = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G60"] = 1
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "R51d40G90"] = 2

# srat.sub = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final == "Goblet"])
srat.sub = srat.merge.rna

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
srat.sub = AddModuleScore_UCell(srat.sub, features = pathway.list.tmp, assay = "SCT", ncores = 8, name = "", maxRank = 2000)
# ucell.mat = srat.sub@meta.data[, names(pathway.list.tmp)] %>% as.matrix() %>% t()
# ucell.res = sapply(split(colnames(srat.sub), srat.sub$orig.ident), function(cells) {rowMeans(ucell.mat[, cells])})
# 
# library(TCseq)
# ucell.res = ucell.res[rowMaxs(ucell.res) > 0, ]
# tcseq.res = timeclust(ucell.res, algo = "cm", k = 6, standardize = T)
# wrap_plots(timeclustplot(tcseq.res, membership.color = viridis(100, option = "H"), cols = 1), ncol = 3)
# 
# tcseq.res.df = clustCluster(tcseq.res) %>% as.data.frame() %>% tibble::rownames_to_column(var = "pathway")
# trend.list = split(tcseq.res.df$pathway, tcseq.res.df$.)
# names(trend.list) = paste0("m", names(trend.list))

### GSVA ###
Idents(srat.sub) = "orig.ident"
expr = AverageExpression(srat.sub, assays = "MAGIC_SCT", slot = "data")[[1]]
expr = expr[rowSums(expr) > 0, ]  #过滤细胞表达量全为零的基因
expr = as.matrix(expr)
library(GSVA)
params = gsvaParam(expr, pathway.list, minSize = 1, maxSize = Inf, kcdf = "Gaussian", tau = 1, maxDiff = T, absRanking = F)
gsva.res = gsva(params, verbose = T, BPPARAM = BiocParallel::SerialParam(progressbar = T))

library(TCseq)
gsva.res = gsva.res[rowMaxs(gsva.res) > 0, ]
tcseq.res = timeclust(gsva.res, algo = "cm", k = 6, standardize = T, dist = "correlation")
ps = timeclustplot(tcseq.res, membership.color = viridis(100, option = "H"), cols = 1, categories = "")

tcseq.res.df = clustCluster(tcseq.res) %>% as.data.frame() %>% tibble::rownames_to_column(var = "pathway")
trend.list = split(tcseq.res.df$pathway, tcseq.res.df$.)
names(trend.list) = paste0("m", names(trend.list))

tcseq.res.df$cluster = factor(paste0("trend", tcseq.res.df$.), levels = paste0("trend", 1:max(tcseq.res.df$.)))
tcseq.res.df = tcseq.res.df %>% dplyr::arrange(cluster)
gsva.res = gsva.res[tcseq.res.df$pathway, ]
library(ComplexHeatmap)
library(circlize)
column_ha = HeatmapAnnotation(sample = colnames(gsva.res),
                              col = list(sample = palette_scp(colnames(gsva.res), palette = "Set1")),
                              show_annotation_name = F,
                              which = "col")
row_ha = HeatmapAnnotation(trend = tcseq.res.df$cluster, 
                           col = list(trend = palette_scp(unique(tcseq.res.df$cluster), palette = "Paired")),
                           show_annotation_name = F,
                           which = "row")
p = Heatmap(gsva.res, raster_by_magick = F,
            col = colorRamp2(seq(quantile(as.vector(gsva.res), 0.05), quantile(as.vector(gsva.res), 0.95), length = 100), rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100))),
            cluster_rows = F, cluster_columns = F,
            left_annotation = row_ha,
            top_annotation = column_ha,
            show_row_names = F, 
            show_column_names = F)

for (i in 1:length(ps)) {
  ps[[i]] = ps[[i]] & thm_rect() & theme(axis.text.y = element_blank(), axis.ticks.y.left = element_blank()) & ylab("") & ggtitle("") & 
    scale_x_discrete(expand = c(0.01, 0.01)) & theme(plot.margin = margin(0, 50, 0, 10))
  if (i < length(ps)) {
    ps[[i]] = ps[[i]] & theme(axis.text.x = element_blank(), axis.ticks.x.bottom = element_blank())
  }
}
p2 = ggarrange(plotlist = ps, ncol = 1, common.legend = T, legend = "left")

pdf(paste(tmp.dir1, "pathway.trend.pdf", sep = "/"), width = 8, height = 8*1.2)
print(p)
print(p2)
dev.off()

df.up = tcseq.res.df %>% dplyr::filter(cluster == paste0("trend", 4))
df.up = df.up %>% dplyr::filter(paste0("GOBP-", pathway) %in% rownames(ucell.res))
df.up = cbind(df.up, ucell.res[paste0("GOBP-", df.up$pathway), ])
df.up = cbind(df.up, data.frame("membership" = tcseq.res@membership[df.up$pathway, 4]))
# df.up$up.within.ucell = df.up$R51d40G90 >= df.up$R51d40G60 & df.up$R51d40G60 >= df.up$R51d40G30
df.up = df.up %>% dplyr::arrange(desc(membership))
df.up$ord = 1:nrow(df.up)
df.up = df.up %>% dplyr::mutate(show.text = ifelse(grepl("Hypoxia", pathway), wrapText(pathway, 35), NA))
ggplot(df.up, aes(x = ord, y = membership)) + geom_point() + geom_label_repel(aes(x = ord, y = membership, label = show.text), nudge_y = 0.1)
write.csv(df.up, paste(tmp.dir1, "trend.up.csv", sep = "/"), row.names = F)
FeatureStatPlot(srat.sub, stat.by = paste0("GOBP-", df.up$pathway[str_detect(df.up$pathway, "Hypoxia")]), group.by = "orig.ident", plot_type = "box", palette = "Set1", add_trend = T, bg_alpha = 0, 
                comparisons = list(c("R51d40G30", "R51d40G60"), c("R51d40G60", "R51d40G90")), pairwise_method = "wilcox.test", sig_label = "p.format", aspect.ratio = 1)

df.down = tcseq.res.df %>% dplyr::filter(cluster == paste0("trend", 1))
df.down = df.down %>% dplyr::filter(paste0("GOBP-", pathway) %in% rownames(ucell.res))
df.down = cbind(df.down, ucell.res[paste0("GOBP-", df.down$pathway), ])
df.down = cbind(df.down, data.frame("membership" = tcseq.res@membership[df.down$pathway, 1]))
# df.down$up.within.ucell = df.down$R51d40G90 <= df.down$R51d40G60 & df.down$R51d40G60 <= df.down$R51d40G30
df.down = df.down %>% dplyr::arrange(desc(membership))
df.down$ord = 1:nrow(df.down)
df.down = df.down %>% dplyr::mutate(show.text = ifelse(grepl("Mechani", pathway), wrapText(pathway, 35), NA))
ggplot(df.down, aes(x = ord, y = membership)) + geom_point() + geom_label_repel(aes(x = ord, y = membership, label = show.text), nudge_y = 0.1)
write.csv(df.down, paste(tmp.dir1, "trend.down.csv", sep = "/"), row.names = F)
FeatureStatPlot(srat.sub, stat.by = paste0("GOBP-", df.down$pathway[str_detect(df.down$pathway, "Mechani")][1]), group.by = "orig.ident", 
                plot_type = "box", palette = "Set1", add_trend = T, bg_alpha = 0, aspect.ratio = 1)
GroupHeatmap(srat.sub, group.by = "orig.ident", assay = "RNA", slot = "data", 
             features = pathway.list[["Sensory Perception Of Mechanical Stimulus"]],
             cluster_rows = F, cluster_columns = F, group_palette = "Set1",
             heatmap_palcolor = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)))

########################
# DEG #
########################
# downregulate
grps = combn(samples.rna, 2, simplify = F)
DefaultAssay(srat.sub) = "SCT"
df1 = map(grps, function(x) {
  df = FindMarkers(srat.sub, ident.1 = x[1], ident.2 = x[2], group.by = "orig.ident", assay = "SCT", only.pos = T, min.pct = 0.35, recorrect_umi = FALSE, min.cells.group = 2) %>% 
    tibble::rownames_to_column(var = "gene") %>% dplyr::mutate(group = paste(x, collapse = "-"), condition = "downregulate")
  return(df)
}) %>% reduce(., rbind)
# upregulate
grps = combn(rev(samples.rna), 2, simplify = F)
DefaultAssay(srat.sub) = "SCT"
df2 = map(grps, function(x) {
  df = FindMarkers(srat.sub, ident.1 = x[1], ident.2 = x[2], group.by = "orig.ident", assay = "SCT", only.pos = T, min.pct = 0.35, recorrect_umi = FALSE, min.cells.group = 2) %>% 
    tibble::rownames_to_column(var = "gene") %>% dplyr::mutate(group = paste(x, collapse = "-"), condition = "upregulate")
  return(df)
}) %>% reduce(., rbind)
# combine
df.deg = rbind(df1, df2)
df.deg = df.deg %>% dplyr::filter(p_val_adj < 0.05)
write.csv(df.deg, paste(tmp.dir1, "DEG.sig.csv", sep = "/"), row.names = F)
df.deg = read.csv(paste(tmp.dir1, "DEG.sig.csv", sep = "/"))
library(ggVennDiagram)
pdf(paste(tmp.dir1, "DEG.sig.venn.pdf", sep = "/"), width = 8, height = 8)
grps = combn(rev(samples.rna), 2, simplify = F)
ggVennDiagram(split(df.deg$gene, df.deg$group)[unlist(map(grps, ~paste(., collapse = "-")))], 
              set_color = palette_scp(unlist(map(grps, ~paste(., collapse = "-"))), palette = "Set2"), set_size = 3, 
              edge_lty = "solid", edge_size = 0.5) +
  labs(title = "upregulate") + 
  scale_fill_distiller(palette = "RdBu")
grps = combn(samples.rna, 2, simplify = F)
ggVennDiagram(split(df.deg$gene, df.deg$group)[unlist(map(grps, ~paste(., collapse = "-")))], 
              set_color = palette_scp(unlist(map(grps, ~paste(., collapse = "-"))), palette = "Set2"), set_size = 3, 
              edge_lty = "solid", edge_size = 0.5) +
  labs(title = "downregulate") + 
  scale_fill_distiller(palette = "RdBu")
dev.off()
# venn genes
venn.df = map(unique(df.deg$condition), function(x) {
  df = df.deg %>% dplyr::filter(condition == x)
  df = split(df$gene, df$group) %>% 
    unlist() %>% table() %>% 
    as.data.frame() %>% 
    dplyr::rename("gene" = ".") %>% 
    dplyr::mutate(condition = x)
  return(df)
}) %>% reduce(., rbind) %>% dplyr::arrange(desc(condition), desc(Freq))
write.csv(venn.df, paste(tmp.dir1, "DEG.venn.csv", sep = "/"), row.names = F)
venn.df = read.csv(paste(tmp.dir1, "DEG.venn.csv", sep = "/"))
# GO enrich
use.genes = venn.df %>% dplyr::filter(Freq >= 2) %>% {split(.$gene, .$condition)}
df.enrich = enricher_hyper_genger(gene.list = use.genes, minGSSize = 6, maxGSSize = 500, 
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
df.enrich = read.csv(paste(tmp.dir1, "enrichment", "hyper-test", "enrichment.csv", sep = "/"))

# upregulate
# GO-BP
up.go = c("Positive Regulation Of Epithelial Cell Migration", "Epithelial To Mesenchymal Transition", 
          "Transforming Growth Factor Beta Production", "Response To Mechanical Stimulus")
# curated
up.curated = c("Elvidge Hypoxia Up", "Mechanoregulation And Pathology Of Yaptaz Via Hippo And Nonhippo Mechanisms")

# downregulate
# GO-BP
down.go = c("Cytoplasmic Translation", "Cellular Oxidant Detoxification", "Oxidative Phosphorylation", "Regulation Of Wnt Signaling Pathway")
# curated
down.curated = c("Reactome Cellular Response To Chemical Stress", "Reactome Cellular Response To Heat Stress")

# GO for heatmap
df.tmp = df.enrich %>% 
  dplyr::filter(list.name == "upregulate" & collection %in% c("GO_BP") & pathway.name %in% up.go | 
                  list.name == "upregulate" & collection %in% c("curated") & pathway.name %in% up.curated)
df.tmp$genes.hit = gsub(";", "/", df.tmp$genes.hit)
df.tmp$text = paste0(df.tmp$pathway.name, paste0(" - (", df.tmp$collection, ")"))
df.tmp$text = wrapText(df.tmp$text, 80)
df.tmp$genes.hit = wrapText(df.tmp$genes.hit, 80)
df.tmp = df.tmp %>% dplyr::arrange(desc(hit.gene.size), p_val)
df.tmp$text = factor(df.tmp$text, levels = rev(df.tmp$text))
p1 = ggplot(data = df.tmp, aes(x = hit.gene.size, y = text)) +
  geom_bar(width = 0.8, stat = 'identity', fill = "#FF69B4") +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0.5)) +
  theme(axis.text.y = element_blank()) + 
  geom_text(data = df.tmp, aes(x = 0.1, y = text, label = text), size = 4, hjust = 0, vjust = -0.5) + 
  geom_text(data = df.tmp, aes(x = 0.1, y = text, label = genes.hit , color = -log10(p_val)), size = 3, fontface = "italic", hjust = 0, vjust = 1.5) +
  scale_colour_gradientn(colors = rev(viridis(n = 100, option = "F"))) + 
  labs(title = 'upregulate') + 
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank()) + 
  thm_rect()

df.tmp = df.enrich %>% 
  dplyr::filter(list.name == "downregulate" & collection %in% c("GO_BP") & pathway.name %in% down.go | 
                  list.name == "downregulate" & collection %in% c("curated") & pathway.name %in% down.curated)
df.tmp$genes.hit = gsub(";", "/", df.tmp$genes.hit)
df.tmp$text = paste0(df.tmp$pathway.name, paste0(" - (", df.tmp$collection, ")"))
df.tmp$text = wrapText(df.tmp$text, 80)
df.tmp$genes.hit = wrapText(df.tmp$genes.hit, 80)
df.tmp = df.tmp %>% dplyr::arrange(desc(hit.gene.size), p_val)
df.tmp$text = factor(df.tmp$text, levels = rev(df.tmp$text))
p2 = ggplot(data = df.tmp, aes(x = hit.gene.size, y = text)) +
  geom_bar(width = 0.8, stat = 'identity', fill = "lightblue") +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0.5)) +
  theme(axis.text.y = element_blank()) + 
  geom_text(data = df.tmp, aes(x = 0.1, y = text, label = text), size = 4, hjust = 0, vjust = -0.5) + 
  geom_text(data = df.tmp, aes(x = 0.1, y = text, label = genes.hit , color = -log10(p_val)), size = 3, fontface = "italic", hjust = 0, vjust = 1.5) +
  scale_colour_gradientn(colors = rev(viridis(n = 100, option = "G"))) + 
  labs(title = 'downregulate') + 
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank()) + 
  thm_rect()

pdf(paste(tmp.dir1, "heat.GO.bar.pdf", sep = "/"), width = 8, height = 8*0.5)
p1
p2
dev.off()

# draw heat
genes.sorted1 = df.deg %>% 
  dplyr::filter(condition == "downregulate" & gene %in% use.genes$downregulate) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>% 
  dplyr::pull(gene) %>% unique()
genes.sorted2 = df.deg %>% 
  dplyr::filter(condition == "upregulate" & gene %in% use.genes$upregulate) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>% 
  dplyr::pull(gene) %>% unique()

ht = GroupHeatmap(srat.sub, group.by = "orig.ident", assay = "SCT", slot = "data", 
                  features = c(genes.sorted1, genes.sorted2), 
                  features_label = c("SCGB3A2", "HSBP1", "HSP90AA1", "HSPA1A", "YAP1", "HIF1A", "VEGFA"),
                  feature_split = c(rep("downregulate", length(genes.sorted1)), rep("upregulate", length(genes.sorted2))), 
                  feature_split_palcolor = c("downregulate" = "lightblue", "upregulate" = "#FF69B4"),
                  cluster_rows = F, cluster_columns = F, group_palcolor = list(colors.organoid.airway.sample),
                  heatmap_palcolor = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)), 
                  width = 3, height = 6)
pdf(paste(tmp.dir1, "DEG.heat.pdf", sep = "/"), width = 8, height = 8)
ht$plot
dev.off()
