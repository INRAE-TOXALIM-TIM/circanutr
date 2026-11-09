# Set up #######################################################################

library(tidyverse)
library(mixOmics)
library(edgeR)
library(org.Mm.eg.db)

devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr")
library(kuipeR)

load("~/circanutr/RS2205_20250123_120703.RData")
load("~/circanutr/RefedLm.RData")

# A)	Venn des gènes rythmiques en Night feeding / gènes non rythmiques en Night feeding ####
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
# Number of mouse genes detected in the liver and rhythmics in AdLib
# This is the "Rhythmic" set

liver_NF_rhythmic <- rownames(dryList$parameters[Refed.lm$comment %in% c("gain", "sync", "shift"), ])
# Number of mouse genes detected in the liver and rhythmics in NF (gain (3), synchronic (4), and shift (5))

## Plot ########################################################################
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `Rhythmic AdLib` = liver_Adlib_rhythmic,
      `Rhythmic NF` = liver_NF_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c("#828282", "#D8D8D8", "#FFFFFF", "#67987E", "#007fff")
)
# Mouse = 78298 = 45295 + 12657 + 8325 + 6299 + 3156 + 2566  (GRCm39)
# Liver = 33003 = 12657 + 8325 + 6299 + 3156 + 2566          (>=1 Count in AdLib)
# Detected = 20346 = 8325 + 6299 + 3156 + 2566               (>=10 Counts in AdLib)
# Rhythmic AdLib = 6299 + 2566                               (In dryR models rhythmic in AdLib diet)
# Rhythmic NF = 6299 + 3156                                  (In dryR models rhythmic in Night-fed diet)

# Not saved, made in post with Illustrator / Affinity


# B)	ACP AdLib et Night feeding ###############################################
CT_FR <- dgef$samples[grep("CTRL|FAR", dgef$samples$Diet), ]

groups_00_CTFR <- paste0(dgef$samples[grep("CTRL|FAR", dgef$samples$Diet), 11], "_ZT", dgef$samples[grep("CTRL|FAR", dgef$samples$Diet), 12]) %>% factor

CTFR_colors <- c(
  "#C7E9C0",
  "#A1D99B",
  "#74C476",
  "#41AB5D",
  "#238B45",
  "#005A32",
  "#C6DBEF",
  "#9ECAE1",
  "#6BAED6",
  "#4292C6",
  "#2171B5",
  "#084594"
)

pcaCTFR <- pca(t(cpm(
  dgef,
  log = T,
  prior = 1,
  normalized = T
)[, grep("CTRL|FAR", dgef$samples$Diet)]), ncomp = 2)

plotIndiv(
  pcaCTFR,
  comp = 1:2,
  ind.names = T,
  group = groups_00_CTFR,
  style = "ggplot2",
  legend = F,
  scale = F,
  pch = 16,
  col = CTFR_colors,
  legend.title = "Groups",
  legend.title.pch = "Diet",
  cex = 3,
  title = " ",
  Y.label = " ",
  X.label = " "
)

ggsave(
  filename = "3B_NF_ACP.png",
  path = "~/circanutr/output/figures/",
  width = 10,
  height = 10,
  dpi = 600,
  units = "cm"
)

# Paramètres rythmiques ####################################################
# C.	Courbe d’amplitude des gènes rythmiques ##################################
dryList.Refed$results |>
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
  # ylab("log10(rhythmic genes)") +
  geom_line(linewidth = .7) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_color_manual(values = c("#67987E", "#007FFF")) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = .5),
    axis.ticks = element_line(linewidth = .5)
  )

ggsave(
  filename = "3C.NF_Rhythmic_Amplitude.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 5,
  dpi = 600,
  units = "cm"
)


# D. Histogramme du décalage des phases d’expressions des gènes rythmiques. ####
ggplot(
  dryList.Refed$parameters |>
    filter((!is.na(amp_CTRL) &
              !is.na(amp_FAR))) |>
    mutate(phase_shift = ifelse(
      abs((phase_CTRL - phase_FAR)) < 12, abs((phase_CTRL - phase_FAR)), 24 - abs(phase_CTRL - phase_FAR)
    )) |>
    mutate(phase_shift = ifelse(phase_shift == 0, NA, phase_shift), ),
  aes(x = phase_shift)
) +
  geom_histogram(bins = 12,
                 breaks = seq(0, 12, by = 1),
                 fill = "#007FFF") +
  theme_classic() +
  xlim(0, 13) +
  scale_x_continuous(breaks = seq(0, 12, 4)) +
  scale_y_continuous(breaks = seq(0, 300, 50), expand = c(0, 0)) +
  ylab(label = "Gene count") +
  xlab(label = "Phase shift") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    text = element_text(size = 10)
  )

ggsave(
  filename = "3D.NF_Phase_shift.svg",
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
      dryList.Refed$results |>
        rownames_to_column("gene_name") |>
        dplyr::select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_FAR)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_FAR),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "FAR")
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
      # xlab(label = "Acrophase (ZT)") +
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
    # Night Feeding
    ggplot(
      dryList.Refed$results |>
        rownames_to_column("gene_name") |>
        dplyr::select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_FAR)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_FAR),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "FAR")
        )) |>
        filter(rhythmic_in %in% c("FAR", "both")),
      aes(x = phase_FAR, y = after_stat(count), fill = rhythmic_in)
    ) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity"
      ) +
      scale_fill_manual(values = c("#007FFF", "grey70")) +
      # xlab(label = "Acrophase (ZT)") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 550)) +
      # ylab(label = "Genes count") +
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
  filename = "3E.NF_Rhythmic_cumulative_counts.svg",
  path = "~/circanutr/output/figures/",
  width = 9.5,
  height = 4.5,
  dpi = 600,
  units = "cm"
)


# Profils de rythmicité des gènes rythmiques (AdLib ou Night feeding) ##########
source("./R/fun/Compare_lm.R")
Refed.lm <- Compare_lm(dryList.Refed, groupTable, "FAR")

# Selecting genes with same intercept in both AdLib dryList and NF dryList
Robust.Refed <- rownames(dryList$results)[(is.na(dryList$results$a_CTRL) &
                                             is.na(dryList.Refed$results$a_CTRL)) |
                                            (!is.na(dryList$results$a_CTRL) &
                                               !is.na(dryList.Refed$results$a_CTRL))]

# Non-circadian diets integration
Table.Refed <- dryList$parameters |>
  merge(Refed.lm, by = 0, all.x = TRUE) |>
  mutate(Amp_FR = dryList.Refed$parameters$amp_FAR) |>
  filter(Row.names %in% Robust.Refed) |>
  group_by(chosen_model) |>
  mutate(DF_count = n()) |>
  filter(
    (best_model == "circadian" & comment == "sync") |
      (best_model == "circadian" & comment == "shift") |
      (best_model == "flat" & comment == "flat") |
      (best_model == "circadian" & comment == "gain") |
      best_model == "linear"
  ) |>
  group_by(chosen_model, best_model, comment)

# Models selection
Models.Refed <- Table.Refed |>
  summarize(
    Percent = n() / mean(DF_count) * 100,
    Rsqr = mean(Rsqr, na.rm = T),
    Mean_amp = mean(Amp_FR, na.rm = T)
  )

# F.	Tableau des modèles et leurs nombres de gènes ############################
dryList.Refed$results |>
  group_by(chosen_model) |>
  summarise(model_size = n())

# chosen_model model_size
# 1            1       7540
# 2            2       1281
# 3            3       2556
# 4            4       7471
# 5            5       1498

Table.Refed |>
  group_by(chosen_model, best_model) |>
  count(comment) |>
  mutate(Percent = n / 20346 * 100)


# G.	Amplitude moyenne d’expression des gènes rythmiques par modèle et par régime ####
dryList.Refed$parameters |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, amp_CTRL, amp_FAR, chosen_model) |>
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
  # scale_x_discrete(label = c("AL", "NF")) +
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
  filename = "3G_Nf_Mean_amplitude_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 5,
  dpi = 600,
  units = "cm"
)


# H.	Heatmap des modèles par régime ###########################################
Heatmaps_Graph <- function(dryList, model) {
  OrderModL <- dryList.Refed$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_FAR, phase_CTRL)) |>
    row.names()
  
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList.Refed$ncounts[OrderModL, ])[, dryList.Refed$group == "CTRL"] # get genes with the same order than ordermodL
  dataHM$FAR   <- as.data.frame(dryList.Refed$ncounts[OrderModL, ])[, dryList.Refed$group == "FAR"]
  
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "FAR"), function(x) {
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

model_size <- dryList.Refed$results |>
  group_by(chosen_model) |>
  filter(chosen_model != 1) |>
  summarise(model_size = n())

ggpubr::ggarrange(
  plotlist = lapply(2:5, function(i) {
    Heatmaps_Graph(dryList = dryList.Refed, model = i)
  }),
  ncol = 1,
  heights = 2 * (model_size$model_size / max(model_size$model_size)),
  widths = 1
)

ggsave(
  filename = "3H_NF_Rhythmic_heatmaps_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 3,
  height = 10,
  dpi = 600,
  units = "cm"
)


# I.	Histogramme de la distribution des phases par modèle #####################
ggplot(
  dryList.Refed$results |>
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
  scale_fill_manual(values = c("#67987E", "#007FFF")) +
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
  filename = "3I.NF_Rhythmic_cumulative_count_per_model.svg",
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



## FAR ####
GOphase_FAR <- lapply(GOgenes, function (x) {
  Liste <- dryList.Refed$results[rownames(dryList.Refed$results) %in% x, "phase_FAR"]
  Liste
}) # 12k

GOphase_FAR <- GOphase_FAR[sapply(GOphase_FAR, function(x) {
  length(x) > 10
})] # 1048
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_FAR <- lapply(GOphase_FAR, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "fitted")
})

GOkuip_FAR <- do.call(rbind.data.frame, GOkuip_FAR)
# 331 significatifs

GOkuip_FAR$go_name <- sapply(rownames(GOkuip_FAR), Get_go_name)

# Rhythmic cat. in each diet
GO_AdLib <- rownames(GOkuip_AL[GOkuip_AL$p.value < 0.05, ])
GO_NightFed <- rownames(GOkuip_FAR[GOkuip_FAR$p.value < 0.05, ])

# Common significative GO cat
GO_common <- intersect(GO_AdLib, GO_NightFed)

# Rhythmic cat. in one of the diets
GO_rhythm <- unique(c(GO_AdLib, GO_NightFed))

# Hyper-geometric test for clustering GO categories between models
HT_models <- lapply(GO_rhythm, function(x) {
  HyperT_clust(dryList.Refed,
               Background = rownames(dryList.Refed$results)[rownames(dryList.Refed$results) %in% unique(unlist(GOgenes))],
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
  filename = "3J_NF_GO_cat_per_model.svg",
  path = "~/circanutr/output/figures/",
  width = 5,
  height = 8,
  dpi = 600,
  units = "cm"
)
