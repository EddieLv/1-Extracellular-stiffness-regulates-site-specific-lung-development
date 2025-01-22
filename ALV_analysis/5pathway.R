source("/media/biogenger/D/scripts/Downstream/R_seuratV5/initialize_analysis_env.R")
initialize.genger_env()

python_exec = "/home/biogenger/miniconda3/bin/python"

samples.rna = c("205-30", "205-60", "205-90")

upstream.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar/2merge_annotation/"
work.dir = "/media/biogenger/D/Projects/LZY/Analysis/alveolar_old_1/alveolar"

tmp.dir1 = paste(work.dir, "5pathway", sep = "/")
dir.create(tmp.dir1, recursive = T)

library(SCP)
### individual ###
srat.merge.rna = qread(paste(upstream.dir, "srat.anno.qs", sep = "/"))
srat.merge.rna$hard.num = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "205-30"] = 0
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "205-60"] = 1
srat.merge.rna$hard.num[srat.merge.rna$orig.ident == "205-90"] = 2

# srat.sub = subset(srat.merge.rna, cells = colnames(srat.merge.rna)[srat.merge.rna$celltype2_final == "Goblet"])
srat.sub = srat.merge.rna

########################
# pathway #
########################
db.merge.df.list = readRDS("/media/biogenger/D/enrich_database/db.merge.df.list.rds")
anno.df = db.merge.df.list$db.merge.df.final %>% dplyr::filter(species == "human" & database == "MsigDB" & collection %in% "GO_BP")
pathway.list = split(anno.df$pathway.gene, anno.df$pathway.name)

########################
# DEG #
########################
# downregulate
grps = combn(samples.rna, 2, simplify = F)
DefaultAssay(srat.sub) = "SCT"
df1 = map(grps, function(x) {
  df = FindMarkers(srat.sub, ident.1 = x[1], ident.2 = x[2], group.by = "orig.ident", assay = "SCT", only.pos = T, min.pct = 0.25, recorrect_umi = FALSE, min.cells.group = 2) %>% 
    tibble::rownames_to_column(var = "gene") %>% dplyr::mutate(group = paste(x, collapse = "-"), condition = "downregulate")
  return(df)
}) %>% reduce(., rbind)
# upregulate
grps = combn(rev(samples.rna), 2, simplify = F)
DefaultAssay(srat.sub) = "SCT"
df2 = map(grps, function(x) {
  df = FindMarkers(srat.sub, ident.1 = x[1], ident.2 = x[2], group.by = "orig.ident", assay = "SCT", only.pos = T, min.pct = 0.25, recorrect_umi = FALSE, min.cells.group = 2) %>% 
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
# use.genes = list(
#   "all-up" = venn.df$gene[venn.df$Freq == 2 & venn.df$condition == "upregulate"],
#   "G60-G30" = df.deg$gene[df.deg$group == "205-60-205-30"]
# )
# use.genes$`G60-G30` = use.genes$`G60-G30`[!use.genes$`G60-G30` %in% c(df.deg$gene[df.deg$group == "205-90-205-30"], df.deg$gene[df.deg$group == "205-90-205-60"])]
df.enrich = enricher_hyper_genger(gene.list = use.genes, minGSSize = 6, maxGSSize = 500, 
                                  db.df = db.merge.df.list$db.merge.df.final, pathway.id = "pathway.id", pathway.name = "pathway.name", pathway.gene = "pathway.gene", 
                                  database = c("MsigDB", "KEGG"), collection = c("curated", "KEGG", "wiki", "GO_BP", "hallmark"), specy = "human")
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
up.go = c("Morphogenesis Of A Polarized Epithelium", "Positive Regulation Of Cell Cycle", "Erbb Signaling Pathway")
# curated
up.curated = c("Elvidge Hypoxia Up", "Hippo Signaling Regulation Pathways")

# downregulate
# GO-BP
# down.go = c("Regulation Of Cellular Respiration", "Positive Regulation Of Peptidase Activity", "Cytoplasmic Translation", "Positive Regulation Of Type I Interferon Production", "Oxidative Phosphorylation")
down.go = c("Protein Transport Within Lipid Bilayer", "Cell Cell Signaling By Wnt", "Regulation Of Plasma Membrane Organization", "Positive Regulation Of Stem Cell Population Maintenance")
# curated
down.curated = c("Elvidge Hypoxia Up")

# GO for heatmap
df.tmp = df.enrich %>% 
  dplyr::filter(list.name == "all-up" & collection %in% c("GO_BP") & pathway.name %in% up.go | 
                  list.name == "all-up" & collection %in% c("curated") & pathway.name %in% up.curated)
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
  dplyr::filter(list.name == "G60-G30" & collection %in% c("GO_BP") & pathway.name %in% down.go | 
                  list.name == "G60-G30" & collection %in% c("curated") & pathway.name %in% down.curated)
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
  labs(title = 'G60-G30') + 
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
                  features_label = c("HSPA8", "HSPA1B", "ABCA3", "SFTPB", "HOPX", "AGER", "LMO7", "HIF1A", "YAP1"),
                  feature_split = c(rep("downregulate", length(genes.sorted1)), rep("upregulate", length(genes.sorted2))), 
                  feature_split_palcolor = c("downregulate" = "lightblue", "upregulate" = "#FF69B4"),
                  cluster_rows = F, cluster_columns = F, group_palcolor = list(colors.organoid.airway.sample),
                  heatmap_palcolor = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)), 
                  width = 3, height = 6)
pdf(paste(tmp.dir1, "DEG.heat.pdf", sep = "/"), width = 8, height = 8)
ht$plot
dev.off()
