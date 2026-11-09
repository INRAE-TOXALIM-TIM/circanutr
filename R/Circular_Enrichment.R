devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr")
library(kuipeR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ViSEAGO")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GOfuncR")

library(GOSemSim)
library(org.Mm.eg.db)
library(ggplot2)
library(tidyverse)

source("./R/fun/Get_go_name.R")
source("./R/fun/HyperT_clust.R")
source("./R/fun/Chronoplot_text.R")
load("./CircanutrEnv.RData")

GOgenes_df <- as.data.frame(org.Mm.egGO2EG)

GOgenes_df$ENSEMBL <- mapIds(org.Mm.eg.db,
                             GOgenes_df$gene_id,
                             column = "ENSEMBL",
                             keytype = "ENTREZID")

GOgenes_df <- GOgenes_df[GOgenes_df$Evidence != 'IEA' &
                           GOgenes_df$Ontology == 'BP', ]  # possibility to avoid some specific Evidence code (IEA = Inferred from electronic annotation)


GOgenes <- split(GOgenes_df$ENSEMBL, GOgenes_df$go_id)



# AdLib
GOphase_AL <- lapply(GOgenes, function (x) {
  Liste <- dryList$results[rownames(dryList$results) %in% x, "phase_CTRL"]
  Liste
}) # 12k

GOphase_AL <- GOphase_AL[sapply(GOphase_AL, function(x) {
  length(x) > 10
})] # 1193
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_AL <- lapply(GOphase_AL, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "computed")
})

GOkuip_AL <- do.call(rbind.data.frame, GOkuip_AL)
#### 264 significatifs

GOkuip_AL$go_name <- sapply(rownames(GOkuip_AL), Get_go_name)



# DF
GOphase_DF <- lapply(GOgenes, function (x) {
  Liste <- dryList$results[rownames(dryList$results) %in% x, "phase_RF"]
  Liste
}) # 12k

GOphase_DF <- GOphase_DF[sapply(GOphase_DF, function(x) {
  length(x) > 10
})] # 1048
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_DF <- lapply(GOphase_DF, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "computed")
})

GOkuip_DF <- do.call(rbind.data.frame, GOkuip_DF)
# 331 significatifs

GOkuip_DF$go_name <- sapply(rownames(GOkuip_DF), Get_go_name)



# Rhythmic cat. in each diet
GO_AdLib <- rownames(GOkuip_AL[GOkuip_AL$p.value < 0.05, ])
GO_DayFed <- rownames(GOkuip_DF[GOkuip_DF$p.value < 0.05, ])



# Common significative GO cat
GO_common <- intersect(GO_AdLib, GO_DayFed)



# Rhythmic cat. in one of the diets
GO_rhythm <- unique(c(GO_AdLib, GO_DayFed))



# Hyper-geometric test for clustering GO categories between models
HT_models <- lapply(GO_rhythm, function(x) {
  HyperT_clust(dryList, Background = rownames(dryList$results)[rownames(dryList$results) %in% unique(unlist(GOgenes))], GOgenes[x])
})

HT_models <- do.call(rbind.data.frame, HT_models)


# Histogram of enriched categories in each model
HT_models_df <- HT_models |>
  filter(pval < 0.05) |>
  group_by(go_id) |>
  filter(n() == 1) |>
  filter(cluster != 1) |>
  mutate(go_name = Get_go_name(go_id)) |>
  ungroup() |>
  mutate(adj.pval = p.adjust(pval, method = "BH")) |>
  filter(adj.pval < 0.04) # 100 cat


d <- GOSemSim::godata('org.Mm.eg.db', ont = "BP", computeIC = FALSE)

# My_GOsemSim <- GOSemSim::godata(
#   annoDb = org.Mm.eg.db::org.Mm.eg.db,
#   keytype = "ENTREZID",
#   ont = 'BP',
#   computeIC = TRUE,
#   processTCSS = FALSE,
#   cutoff = NULL
# )


my_go_ids <- HT_models_df$go_id

my.Sim.matrix <- sapply(my_go_ids, function(x) {
  GOSemSim::goSim(x, my_go_ids, d, measure = "Wang")
}) |> data.frame()

my.Sim.matrix <- my.Sim.matrix[, !is.na(my.Sim.matrix[1, ])]

row.names(my.Sim.matrix) <- colnames(my.Sim.matrix)

GOheatMap <- pheatmap::pheatmap(
  my.Sim.matrix,
  clustering_method = "complete",
  show_colnames = F,
  silent = T
)

ggpubr::ggarrange(plotlist = list(GOheatMap[[4]]))

my_clust <- GOheatMap$tree_row

###############################################################################
my_clust_cut <- as.data.frame(cutree(my_clust, k = 21)) |>
  rename("GO_cluster" = `cutree(my_clust, k = 21)`) |>
  rownames_to_column(var = "go_id") |>
  mutate(go_id = gsub("GO.", "GO:", go_id)) |>
  left_join(HT_models_df) |>
  group_by(GO_cluster) |>
  arrange(GO_cluster, -m) |>
  slice(1) |>
  filter(!is.na(go_name))
# Remettre à k = 20 si pas de NA
###############################################################################

HT_model_plot <- ggplot(
  HT_models_df |>
    filter(go_id %in% my_clust_cut$go_id),
  aes(
    x = -log10(pval),
    y = fct_reorder(go_name, -log10(pval)),
    fill = log10(q)
  )
) +
  geom_bar(stat = "identity", alpha = 0.66) +
  facet_grid(scales = "free",
             rows = vars(cluster),
             space = "free") +
  scale_fill_viridis_c(
    limits = c(0, NA),
    begin = 0.5,
    end = 1,
    direction = -1,
    option = "inferno"
  ) +
  geom_text(
    aes(x = 0, label = go_name),
    size = 10,
    size.unit = "pt",
    hjust = 0
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    text = element_text(size = 10)
  )


# ggsave(
#   filename = "GO_cat_per_model.png",
#   plot = HT_model_plot,
#   path = "./output/graphs/",
#   width = 13,
#   height = 13,
#   dpi = 600,
#   units = "cm"
# )


path <- "./output/graphs/Chronoplots/"

# Chronoplot of both diets for one category by counts
lapply(my_clust_cut$go_id, function(x) {
  plot <- chronoplot(
    circ.data = data.frame(
      1:max(length(GOphase_AL[[x]]), length(GOphase_DF[[GO_common[x]]])),
      AdLib = GOphase_AL[[x]],
      DayFed = GOphase_DF[[x]]
    ),
    limits = c(0, 24),
    n.class = 12,
    method = "count",
    threshold = 0.05,
    color = c("cyan", "orange"),
    alternative = "greater",
    log10 = F,
    axis.text.size = 8,
    points.size = 0.5
  )
  plot[[2]] +
    labs(title = insert_line_breaks(Get_go_name(x) , 25)) +
    theme(text = element_text(size = 8), legend.position = "none")
  
  # ggsave(
  #   filename = paste("./Count/Chronoplot_", Get_go_name(x), ".png", sep = ""),
  #   path = path,
  #   width = 6,
  #   height = 6,
  #   dpi = 600,
  #   units = "cm"
  # )
})

# plots of p-values
lapply(my_clust_cut$go_id, function(x) {
  plot <- chronoplot(
    circ.data = data.frame(
      1:max(length(GOphase_AL[[x]]), length(GOphase_DF[[GO_common[x]]])),
      AdLib = GOphase_AL[[x]],
      DayFed = GOphase_DF[[x]]
    ),
    limits = c(0, 24),
    n.class = 12,
    method = "stat",
    threshold = 0.05,
    color = c("cyan", "orange"),
    alternative = "greater",
    log10 = F,
    axis.text.size = 8,
    points.size = 0.5
  )
  plot[[2]] +
    labs(title = insert_line_breaks(Get_go_name(x) , 25)) +
    theme(text = element_text(size = 8), legend.position = "none")
  
  ggsave(
    filename = paste("./Pvalue/Chronoplot_Pvalue_", Get_go_name(x), ".png", sep = ""),
    path = path,
    width = 6,
    height = 6,
    dpi = 600,
    units = "cm"
  )
})
