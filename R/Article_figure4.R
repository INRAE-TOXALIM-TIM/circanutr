# Set up #######################################################################

library(tidyverse)
library(mixOmics)
library(edgeR)
library(org.Mm.eg.db)

devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr")
library(kuipeR)

source("~/circanutr/R/fun/Get_go_name.R")
source("~/circanutr/R/fun/HyperT_clust.R")

load("~/circanutr/FastingLm.RData")
load("~/circanutr/RS2205_20250123_120703.RData")

# A)	Venn des gènes rythmiques en Fasting / gènes non rythmiques en Fasting ###
#     comparés aux gènes rythmiques en AdLib
tpm <- read.table(file = "~/circanutr/data/salmon.merged.gene_tpm.tsv", dec = ".", header = T)
rownames(tpm) <- tpm[, 1]   # Row names are gene_id
tpm <- tpm[, -1]

gn_total <- as.vector(rownames(tpm))
# A vector of all genes names from the ENSEMBL GRCm39 Mus musculus genome
# This is the "Mouse" set

OneCount <- read.table(file = "~/circanutr/data/OneCountTable.csv", sep = ";")
OneCount <- OneCount[-1, ]
# Number of mouse genes detected in "liver"Mouse" with >= 1 count in AdLib

FilterByCount <- data.frame(t(rbind(gn_total, OneCount$V2)))
FilterByCount$V2 <- as.logical(FilterByCount$V2)

Liver <- subset(FilterByCount, V2)
# This is the "Liver" set

liver_detected <- as.matrix(
  read.csv2(
    "~/circanutr/data/FilterCountTable.csv",
    stringsAsFactors = FALSE,
    row.names = 1
  )
)
# Number of mouse genes detected in the liver with >= 10 Counts in AdLib
# This is the "Detected" set

liver_Adlib_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(2, 4, 5), ])
# Number of mouse genes detected in the liver and rhythmics in AdLib (loss (2), synchronic (4), and shift (5))
# This is the "Rhythmic" set

liver_FA_rhythmic <- rownames(dryList$parameters[Fasting.lm$comment %in% c("gain", "sync", "shift"), ])
# Number of mouse genes detected in the liver and rhythmics in FA (gain (3), synchronic (4), and shift (5))

## Plot ########################################################################
Venn <- plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `Rhythmic AdLib` = liver_Adlib_rhythmic,
      `Rhythmic FA` = liver_FA_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c("#828282", "#D8D8D8", "#FFFFFF", "#67987E", "#E0115F")
)
# Mouse = 45295 + 12657 + 8524 + 5189 + 3676 + 2957          (GRCm39)
# Liver = 12657 + 8524 + 5189 + 3676 + 2957                  (>=1 Count in AdLib)
# Detected = 8524 + 5189 + 3676 + 2957                       (>=10 Counts in AdLib)
# Rhythmic AdLib = 5189 + 3676                               (In dryR models rhythmic in AdLib diet)
# Rhythmic FA = 5189 + 2957                                  (In dryR models rhythmic in Fasting diet)


ggsave(
  filename = "4A_FA_Venn.svg",
  plot = Venn,
  path = "~/circanutr/output/figures/",
  width = 9,
  height = 9,
  dpi = 600,
  units = "cm"
)

# B)	ACP AdLib et Fasting #####################################################
CT_FA <- dgef$samples[grep("CT|^F$", dgef$samples$Diet), ]

groups_00_CTFA <- paste0(dgef$samples[grep("CT|^F$", dgef$samples$Diet), 11], "_ZT", dgef$samples[grep("CT|^F$", dgef$samples$Diet), 12]) %>% factor

CTFA_colors <- c(
  "#C7E9C0",
  "#A1D99B",
  "#74C476",
  "#41AB5D",
  "#238B45",
  "#005A32",
  "#F4B3CB",
  "#EE7BA6",
  "#E84481",
  "#B11650",
  "#85123D",
  "#501028"
)

pcaCTFA <- pca(t(cpm(
  dgef,
  log = T,
  prior = 1,
  normalized = T
)[, grep("CTRL|^F$", dgef$samples$Diet)]), ncomp = 2)

plotIndiv(
  pcaCTFA,
  comp = 1:2,
  ind.names = T,
  group = groups_00_CTFA,
  style = "ggplot2",
  legend = F,
  scale = F,
  pch = 16,
  col = CTFA_colors,
  legend.title = "Groups",
  legend.title.pch = "Diet",
  cex = 3,
  title = " ",
  Y.label = " ",
  X.label = " "
)

ggsave(
  filename = "4B_FA_ACP.svg",
  path = "~/circanutr/output/figures/",
  width = 10,
  height = 10,
  dpi = 600,
  units = "cm"
)

# Paramètres rythmiques ########################################################
# C.	Courbe d’amplitude des gènes rythmiques ##################################
dryList.Fasting$results |>
  rownames_to_column("gene_name") |>
  dplyr::select(gene_name, contains("amp_") &
                  !contains("relamp_")) |>
  pivot_longer(!gene_name, names_to = "Diet", values_to = "log2FC") |>
  group_by(Diet) |>
  filter(!is.na(log2FC)) |>
  arrange(log2FC) |>
  mutate(`number of genes` = order(rev(log2FC)),
         Diet = gsub("amp_", "", Diet)) |>
  ggplot(aes(
    x = log2FC,
    y = log10(`number of genes`),
    colour = Diet
  )) +
  xlab("Amplitude(log2FC)") +
  geom_line(linewidth = .7) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_color_manual(values = c("#67987E", "#e0115f")) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = .5),
    axis.ticks = element_line(linewidth = .5)
  )

ggsave(
  filename = "4C.FA_Rhythmic_Amplitude.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 5,
  dpi = 600,
  units = "cm"
)

# D. Histogramme du décalage des phases d’expressions des gènes rythmiques. ####
ggplot(
  dryList.Fasting$parameters |>
    filter((!is.na(amp_CTRL) &
              !is.na(amp_F))) |>
    mutate(phase_shift = ifelse(
      abs((phase_CTRL - phase_F)) < 12, abs((phase_CTRL - phase_F)), 24 - abs(phase_CTRL - phase_F)
    )) |>
    mutate(phase_shift = ifelse(phase_shift == 0, NA, phase_shift), ),
  aes(x = phase_shift)
) +
  geom_histogram(bins = 12,
                 breaks = seq(0, 12, by = 1),
                 fill = "#e0115f") +
  theme_classic() +
  xlim(0, 13) +
  scale_x_continuous(breaks = seq(0, 12, 4)) +
  scale_y_continuous(
    breaks = seq(0, 300, 100),
    expand = c(0, 0),
    limits = c(0, 300)
  ) +
  ylab(label = "Gene count") +
  xlab(label = "Phase shift (ZT)") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    text = element_text(size = 10)
  )

ggsave(
  filename = "4D.FA_Phase_shift.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 5,
  dpi = 600,
  units = "cm"
)

# E.	Histogramme des répartitions des phases des gènes rythmiques sur 24h. ####
ggpubr::ggarrange(
  plotlist = list(
    # Ad Libitum
    ggplot(
      dryList.Fasting$results |>
        rownames_to_column("gene_name") |>
        dplyr::select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_F)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_F),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "F")
        )) |>
        filter(rhythmic_in %in% c("CTRL", "both")),
      aes(x = phase_CTRL, y = after_stat(count), fill = rhythmic_in)
    ) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity"
      ) +
      scale_fill_manual(values = c("#67987E", "grey70")) +
      ylab(label = "Genes count") +
      scale_x_continuous(
        breaks = seq(0, 24, 4),
        expand = c(0, 0),
        limits = c(0, 24)
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 550)) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.title = element_blank()
      )
    ,
    # Fasting
    ggplot(
      dryList.Fasting$results |>
        rownames_to_column("gene_name") |>
        dplyr::select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_F)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_F),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "F")
        )) |>
        filter(rhythmic_in %in% c("F", "both")),
      aes(x = phase_F, y = after_stat(count), fill = rhythmic_in)
    ) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity"
      ) +
      scale_fill_manual(values = c("#e0115f", "grey70")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 550)) +
      scale_x_continuous(
        breaks = seq(0, 24, 4),
        expand = c(0, 0),
        limits = c(0, 24)
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  )
)

ggsave(
  filename = "4E.FA_Rhythmic_cumulative_counts.svg",
  path = "~/circanutr/output/figures/",
  width = 9.5,
  height = 4.5,
  dpi = 600,
  units = "cm"
)

Robust.Fasting <- rownames(dryList$results)[(is.na(dryList$results$a_CTRL) &
                                               is.na(dryList.Fasting$results$a_CTRL)) |
                                              (!is.na(dryList$results$a_CTRL) &
                                                 !is.na(dryList.Fasting$results$a_CTRL))]

# Non-circadian diets integration
Table.Fasting <- dryList$parameters |>
  merge(Fasting.lm, by = 0, all.x = TRUE) |>
  mutate(Amp_FA = dryList.Fasting$parameters$amp_F) |>
  filter(Row.names %in% Robust.Fasting) |>
  group_by(chosen_model) |>
  mutate(count = n()) |>
  filter(
    (best_model == "circadian" & comment == "sync") |
      (best_model == "circadian" & comment == "shift") |
      (best_model == "flat" & comment == "flat") |
      (best_model == "circadian" & comment == "gain") |
      best_model == "linear"
  ) |>
  group_by(chosen_model, best_model, comment)

# Models selection
Models.Fasting <- Table.Fasting |>
  summarize(
    Percent = n() / mean(count) * 100,
    Rsqr = mean(Rsqr, na.rm = T),
    Mean_amp = mean(Amp_FA, na.rm = T)
  )

# F.	Tableau des modèles et leurs nombres de gènes ############################
dryList.Fasting$results |>
  group_by(chosen_model) |>
  summarise(model_size = n())

# chosen_model model_size
# 1            1       7936
# 2            2       1756
# 3            3       3293
# 4            4       5641
# 5            5       1720

Table.Fasting |>
  group_by(chosen_model, best_model) |>
  count(comment) |>
  mutate(Percent = n / 20346 * 100)

# G.	Amplitude moyenne d’expression des gènes rythmiques par modèle et par régime ####
dryList.Fasting$parameters |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, amp_CTRL, amp_F, chosen_model) |>
  pivot_longer(!c(gene_id, chosen_model), names_to = "Diet") |>
  group_by(chosen_model, Diet) |>
  summarise(value = mean(value)) |>
  mutate(`Mean\namplitude` = ifelse(is.na(value), 0, value),
         Diet = gsub("amp_", "", Diet)) |>
  ggplot(aes(
    x = Diet,
    y = -chosen_model,
    fill = `Mean\namplitude`
  )) +
  geom_tile() +
  scale_fill_viridis_c(
    begin = 0.5,
    end = 1,
    
    direction = -1,
    option = "inferno"
  ) +
  coord_fixed() +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = "4G_FA_Mean_amplitude_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 5,
  dpi = 600,
  units = "cm"
)

# H.	Heatmap des modèles par régime ###########################################
Heatmaps_Graph <- function(dryList, model) {
  OrderModL <- dryList$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_F, phase_CTRL)) |>
    row.names()
  
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "CTRL"] # get genes with the same order than ordermodL
  dataHM$FA   <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "F"]
  
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "FA"), function(x) {
    pheatmap::pheatmap(
      dataHM[[x]],
      scale = "row",
      color = colorRampPalette(c("blue", "yellow"))(100),
      cluster_rows = F,
      cluster_cols = F,
      show_rownames = F,
      legend = F,
      show_colnames = F,
      silent = T
    )[[4]]
  }))
}

model_size <- dryList.Fasting$results |>
  group_by(chosen_model) |>
  filter(chosen_model != 1) |>
  summarise(model_size = n())

ggpubr::ggarrange(
  plotlist = lapply(2:5, function(i) {
    Heatmaps_Graph(dryList = dryList.Fasting, model = i)
  }),
  ncol = 1,
  heights = 2 * (model_size$model_size / max(model_size$model_size)),
  widths = 1
)

ggsave(
  filename = "4H_FA_Rhythmic_heatmaps_per_model.png",
  path = "~/circanutr/output/figures/",
  width = 3,
  height = 10,
  dpi = 600,
  units = "cm"
)

# I.	Histogramme de la distribution des phases par modèle #####################
ggplot(
  dryList.Fasting$results |>
    rownames_to_column("gene_name") |>
    dplyr::select(gene_name, contains("phase_"), chosen_model) |>
    pivot_longer(!c(gene_name, chosen_model), names_to = "Condition") |>
    filter(chosen_model != 1),
  aes(x = value, y = after_stat(count), fill = Condition)
) +
  geom_histogram(
    bins = 24,
    breaks = seq(0, 24, by = 1),
    position = "identity",
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("#67987E", "#e0115f")) +
  ggtitle(label = "Cumulative counts\nper model\nof rhythmic genes") +
  ylab(label = "Gene count") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 450),
    breaks = c(0, 400, 200)
  ) +
  scale_x_continuous(
    breaks = seq(0, 24, 4),
    expand = c(0, 0),
    limits = c(0, 24)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 10),
    strip.background = element_blank()
  ) +
  facet_wrap(ncol = 1,
             facets = vars(chosen_model),
             scales = "free_y")

ggsave(
  filename = "4I.FA_Rhythmic_cumulative_count_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 3,
  height = 8,
  dpi = 600,
  units = "cm"
)

# J.	Histogramme des fonctions ontologiques retrouvées parmi les modèles ######
GOgenes_df <- as.data.frame(org.Mm.egGO2EG)

GOgenes_df$ENSEMBL <- mapIds(org.Mm.eg.db,
                             GOgenes_df$gene_id,
                             column = "ENSEMBL",
                             keytype = "ENTREZID")

GOgenes_df <- GOgenes_df[GOgenes_df$Evidence != 'IEA' &
                           GOgenes_df$Ontology == 'BP', ]  # possibility to avoid some specific Evidence code (IEA = Inferred from electronic annotation)


GOgenes <- split(GOgenes_df$ENSEMBL, GOgenes_df$go_id)

## AdLib ####
GOphase_AL <- lapply(GOgenes, function (x) {
  Liste <- dryList$results[rownames(dryList$results) %in% x, "phase_CTRL"]
  Liste
}) # 12k

GOphase_AL <- GOphase_AL[sapply(GOphase_AL, function(x) {
  length(x) > 10
})] # 1193
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_AL <- lapply(GOphase_AL, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "fitted")
})

GOkuip_AL <- do.call(rbind.data.frame, GOkuip_AL)
#### 264 significatifs

GOkuip_AL$go_name <- sapply(rownames(GOkuip_AL), Get_go_name)



## FA ####
GOphase_FA <- lapply(GOgenes, function (x) {
  Liste <- dryList.Fasting$results[rownames(dryList.Fasting$results) %in% x, "phase_F"]
  Liste
}) # 12k

GOphase_FA <- GOphase_FA[sapply(GOphase_FA, function(x) {
  length(x) > 10
})] # 1048
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_FA <- lapply(GOphase_FA, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "fitted")
})

GOkuip_FA <- do.call(rbind.data.frame, GOkuip_FA)
# 331 significatifs

GOkuip_FA$go_name <- sapply(rownames(GOkuip_FA), Get_go_name)

# Rhythmic cat. in each diet
GO_AdLib <- rownames(GOkuip_AL[GOkuip_AL$p.value < 0.05, ])
GO_Fasting <- rownames(GOkuip_FA[GOkuip_FA$p.value < 0.05, ])

# Common significative GO cat
GO_common <- intersect(GO_AdLib, GO_Fasting)

# Rhythmic cat. in one of the diets
GO_rhythm <- unique(c(GO_AdLib, GO_Fasting))

# Hyper-geometric test for clustering GO categories between models
HT_models <- lapply(GO_rhythm, function(x) {
  HyperT_clust(dryList.Fasting,
               Background = rownames(dryList.Fasting$results)[rownames(dryList.Fasting$results) %in% unique(unlist(GOgenes))],
               GOgenes[x])
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
  filter(adj.pval < 0.05) # 100 cat


d <- GOSemSim::godata(annoDb = org.Mm.eg.db::org.Mm.eg.db, ont = "BP", computeIC = FALSE)

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

my_clust_cut <- as.data.frame(cutree(my_clust, k = 21)) |>
  dplyr::rename("GO_cluster" = 1) |>
  rownames_to_column(var = "go_id") |>
  mutate(go_id = gsub("GO.", "GO:", go_id)) |>
  left_join(HT_models_df) |>
  group_by(GO_cluster) |>
  arrange(GO_cluster, -m) |>
  dplyr::slice(1) |>
  filter(!is.na(go_name))

ggplot(
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
  # geom_text(
  #   aes(x = 0, label = go_name),
  #   size = 8,
  #   size.unit = "pt",
  #   hjust = 0
  # ) +
  scale_x_continuous(limits = c(0, NA)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    text = element_text(size = 10),
    legend.position = "top"
  )


ggsave(
  filename = "4J_FA_GO_cat_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 5,
  height = 8,
  dpi = 600,
  units = "cm"
)
