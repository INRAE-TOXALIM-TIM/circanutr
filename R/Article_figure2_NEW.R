# Set up #######################################################################
devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr/")
library(kuipeR)

library(tidyverse)
library(GOSemSim)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)

source("./R/fun/Get_go_name.R")

load("~/circanutr/DryList_BigThree_BGEnv.RData")



# A. Venn AL rhythmic genes #######################################################
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

## Plot ########################################################################
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `Rhythmic` = liver_Adlib_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c("#828282", "#D8D8D8", "#FFFFFF", "#67987E")
)
# Mouse = 78298 = 45295 + 12657 + 11481 + 8865  (GRCm39)
# Liver = 33003 = 12657 + 11481 + 8865          (>=1 Count in AdLib)
# Detected = 20346 = 11481 + 8865               (>=10 Counts in AdLib)
# Rhythmic = 8865                               (In dryR models rhythmic in AdLib diet)



# B. Rhythmics lncRNA compared to other rhythmics RNA in AdLib #################
RawGenesTable <- read.table(
  file = "./data/RawGenesTable.csv",
  dec = ".",
  sep = ";",
  header = TRUE,
  row.names = 1
)

Biotype.Model <- cbind.data.frame(
  gene_id = rownames(dryList$parameters),
  biotype = RawGenesTable$gene_biotype[RawGenesTable$gene_id %in% rownames(dryList$parameters)],
  model = as.factor(dryList$parameters$chosen_model)
) |>
  mutate(
    biotype = ifelse(biotype %in% c("protein_coding", "lncRNA"), biotype, "other"),
    biotype = as.factor(biotype),
    biotype = forcats::fct_recode(biotype, "mRNA" = "protein_coding", "Others" = "other"),
    model_group = ifelse(model %in% c(3, 4, 5), "Rhythmic", "Arrhythmic"),
    model_group = as.factor(model_group)
  )

Biotype.Model$biotype <- factor(Biotype.Model$biotype, levels = c("Others", "lncRNA", "mRNA"))

Biotype.Model.percent <- Biotype.Model |>
  group_by(model_group) |>
  dplyr::count(biotype) |>
  mutate(Percent = n / sum(n) * 100) |>
  mutate(Log10 = log10(Percent))

## Plot ########################################################################
ggplot(Biotype.Model.percent,
       aes(x = Percent, y = model_group, fill = biotype)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_x_continuous(expand = 0, breaks = seq(0, 100, by = 20)) +
  scale_y_discrete(expand = 0.15, limits = rev) +
  coord_flip() +
  scale_fill_manual(values = c("#d9d9d9", "#3a6a52", "#96c9ad")) +
  labs(fill = "RNA types") +
  theme_classic() +
  theme(
    aspect.ratio = 2 / 1,
    text = element_text(size = 14),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(t = 5, b = 5)
  )

ggsave(
  filename = "B_RNA.svg",
  path = "~/circanutr/output/",
  width = 9,
  height = 4,
  dpi = 600,
  units = "cm"
)



# C. Venn AL rhythmic genes ####################################################
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

liver_NF_rhythmic <- rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model %in% c(3, 5, 6, 9:15), ])
# This is the "NF Rhythmic" set

liver_DF_rhythmic <- rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model %in% c(4, 7:15), ])
# This is the "DF Rhythmic" set

## Plot ########################################################################
plot(
  eulerr::euler(
    list(
      `AL Rhythmic` = liver_Adlib_rhythmic,
      `NF Rhythmic` = liver_NF_rhythmic,
      `DF Rhythmic` = liver_DF_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c("#67987E", "#007FFF", "#FF7F00")
)
# Mouse = 78298 = 45295 + 12657 + 9988 + 10358  (GRCm39)
# Liver = 33003 = 12657 + 9988 + 10358          (>=1 Count in AdLib)
# Detected = 20346 = 9988 + 10358               (>=10 Counts in AdLib)
# Rhythmic = 10358                               (In dryR models rhythmic in AdLib diet)

# Not saved, made in post with Illustrator / Affinity


# D. Courbe d’amplitude des gènes rythmiques ###################################
dryList_BigThree$results |>
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
  scale_color_manual(values = c("#67987E", "#007FFF", "#FF7F00")) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = .5),
    axis.ticks = element_line(linewidth = .5)
  )

ggsave(
  filename = "2D.AL_NF_DF_Rhythmic_Amplitude.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 6,
  dpi = 600,
  units = "cm"
)


# E. Histogramme du décalage des phases d’expressions des gènes rythmiques. ####
ggplot(
  dryList_BigThree$parameters |>
    filter((
      !is.na(amp_CTRL) &
        !is.na(amp_FAR) &
        !is.na(amp_RF)
    )) |>
    mutate(phase_shift_NF = ifelse(
      abs((phase_CTRL - phase_FAR)) < 12, abs((phase_CTRL - phase_FAR)), 24 - abs(phase_CTRL - phase_FAR)
    )) |>
    mutate(phase_shift_DF = ifelse(
      abs((phase_CTRL - phase_RF)) < 12, abs((phase_CTRL - phase_RF)), 24 - abs(phase_CTRL - phase_RF)
    )) |>
    mutate(phase_shift_NF = ifelse(phase_shift_NF == 0, NA, phase_shift_NF), ) |>
    mutate(phase_shift_DF = ifelse(phase_shift_DF == 0, NA, phase_shift_DF), )
) +
  geom_histogram(
    aes(x = phase_shift_NF),
    bins = 12,
    breaks = seq(0, 12, by = 1),
    fill = "#007FFF"
  ) +
  geom_histogram(
    aes(x = phase_shift_DF),
    bins = 12,
    breaks = seq(0, 12, by = 1),
    fill = "#FF7F00"
  ) +
  theme_classic() +
  xlim(0, 13) +
  scale_x_continuous(breaks = seq(0, 12, 4)) +
  scale_y_continuous(breaks = seq(0, 1300, 100), expand = c(0, 0)) +
  ylab(label = "Gene count") +
  xlab(label = "Phase shift") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    text = element_text(size = 10)
  )

ggsave(
  filename = "2E.NF_Phase_shift.svg",
  path = "~/circanutr/output/figures/",
  width = 4,
  height = 6,
  dpi = 600,
  units = "cm"
)



# F. Chronoplot ################################################################
Phase.table <- data.frame(
  gene_id = rownames(dryList_BigThree$parameters),
  AdLib = dryList_BigThree$parameters[, "phase_CTRL"],
  DayFed = dryList_BigThree$parameters[, "phase_RF"],
  NightFed = dryList_BigThree$parameters[, "phase_FAR"]
)


chronoplot(
  circ.data = Phase.table,
  distribution = NULL,
  limits = c(0, 24),
  n.class = 6,
  after.stat = TRUE,
  alternative = "greater",
  colors = c("#67987E", "#FF7F00", "#007FFF"),
  fills = c("#67987E", "#FF7F00", "#009FFF"),
  pch = NULL,
  axis.text.size = 8,
  points.size = 2,
  smooth.method = "four_points",
  smoothing = 5,
  w = 1 / 16,
  show.dist = FALSE,
  r.label.pos = "center",
  title = NULL
)


ggsave(
  filename = "2F.Chronoplot_three_diets.svg",
  path = "~/circanutr/output/figures/",
  width = 8,
  height = 6,
  dpi = 600,
  units = "cm"
)



# E. Heatmaps of interesting models ############################################
n_models <- 1:max(dryList_BigThree$parameters$chosen_model)
n_GO_per_model <- sapply(Model_Enrich.list, function(x) {
  dim(x)[1]
})
Mdl2HM <- n_models[n_GO_per_model > 0][-1]
# without flat model


## Heatmaps function
Heatmaps_Graph <- function(dryList, model) {
  OrderModL <- dryList$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_FAR, phase_CTRL)) |>
    row.names()
  
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "CTRL"] # get genes with the same order than ordermodL
  dataHM$RF   <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "RF"] # get genes with the same order than ordermodL
  dataHM$FAR   <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "FAR"] # get genes with the same order than ordermodL
  
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "FAR", "RF"), function(x)
  {
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
  }), ncol = 3)
}

model_size <- dryList_BigThree$results |>
  group_by(chosen_model) |>
  filter(chosen_model %in% Mdl2HM) |>
  summarise(model_size = n())

Heatmaps <- ggpubr::ggarrange(
  plotlist = lapply(Mdl2HM, function(i)
  {
    Heatmaps_Graph(dryList = dryList_BigThree, model = i)
  }),
  ncol = 1,
  heights = 2 * (model_size$model_size / max(model_size$model_size)),
  widths = 1
)

ggsave(
  filename = "BIG_heatmaps_per_model.png",
  plot = Heatmaps,
  path = "./output/graphs/",
  width = 3,
  height = 10,
  dpi = 600,
  units = "cm"
)




# F. TF Over-representation Analysis (Hypergeometric test) #####################
TTRUST.df <- read.table("~/circanutr/data/TTRUST_TF_mouse.tsv")

V1 <- AnnotationDbi::mapIds(
  x = org.Mm.eg.db,
  keys = TTRUST.df[, "V1"],
  column = "SYMBOL",
  keytype = "ALIAS",
  multiVals = "first"
)

V2 <- AnnotationDbi::mapIds(
  x = org.Mm.eg.db,
  keys = TTRUST.df[, "V2"],
  column = "SYMBOL",
  keytype = "ALIAS",
  multiVals = "first"
)

TTRUST.df$TF <- V1
TTRUST.df$Targets <- V2


TTRUST.list <- split(TTRUST.df$Targets, TTRUST.df$TF)

TTRUST.list.phase <- lapply(TTRUST.list, function(x) {
  gene_id = FilterGenesTable$gene_id[FilterGenesTable$gene_name %in% x]
  dryList_BigThree$parameters[gene_id, "phase_CTRL"]
})

TTRUST.list.phase.HT <- lapply(TTRUST.list.phase, function(x) {
  q = sum(!is.na(x))
  m = sum(!is.na(dryList_BigThree$parameters$phase_CTRL))
  n = sum(is.na(dryList_BigThree$parameters$phase_CTRL))
  k = length(x)
  p = phyper(q, m, n, k, lower.tail = FALSE)
  if (k < 3) {
    p = NA
  }
  return(p)
})

names(TTRUST.list.phase.HT)[TTRUST.list.phase.HT < 0.05 &
                              !is.na(TTRUST.list.phase.HT)]

# TTRUST.list.phase.HT[["Bmal1"]]
# 0.0613314
# TTRUST.list.phase.HT[["Ppara"]]
# 0.01879963

Kuip.phase <- lapply(TTRUST.list.phase, function(phase) {
  kuiper_test(x = phase,
              limits = c(0, 24),
              method = "estimated")
})

do.call(rbind, Kuip.phase) |>
  filter(p.value < 0.05)
# TF responsables des pics d'expression geniques circadiens

# List AL (19) : responsables des pics d'expression
#         kuiper.stat     p.value
# Fli1       2.056868 0.006734496
# Foxa2      2.077943 0.005781222
# Foxo3      1.747907 0.049814552
# Foxo4      2.155517 0.003239573
# Gata6      2.051378 0.007005304
# Hmga2      2.000107 0.010056431
# Hnf4a      1.914302 0.017923920
# Mlxipl     1.902594 0.019343229
# Nfe2       1.985261 0.011141086
# Nr2f2      2.209442 0.002131534
# Pou2af1    1.909717 0.018468351
# Pou2f1     1.993253 0.010544725
# Pou2f2     1.927172 0.016471523
# Ppara      1.886441 0.021465003
# Six1       1.831730 0.030259754
# Srebf1     1.979428 0.011595407
# Trp63      1.876022 0.022940550
# Usf1       1.909130 0.018539019
# Ybx1       1.813135 0.033896221


# Liste DF : responsables des pics d'expression
#        kuiper.stat      p.value
# Cebpb     1.884364 0.0217522073
# Fos       1.833606 0.0299125739
# Foxc2     1.991999 0.0106363429
# Gtf2i     1.922824 0.0169500433
# Hdac3     2.035227 0.0078605722
# Irf1      1.867336 0.0242383448
# Jun       2.159572 0.0031406225
# Mlxipl    1.882722 0.0219816677
# Myc       1.869633 0.0238889735
# Nfe2l2    1.948948 0.0142519695
# Nfkb1     2.388623 0.0004832521
# Nfya      2.046560 0.0072511268
# Nr2f2     2.135375 0.0037751570
# Nr3c1     1.882718 0.0219822145
# Nr3c2     1.786378 0.0397934372
# Pgr       2.000660 0.0100179543
# Pou2f1    1.951090 0.0140488599
# Rela      1.984488 0.0112004095
# Runx2     1.937705 0.0153620993
# Smad2     1.769997 0.0438252375
# Sp1       2.372785 0.0005541567
# Sp3       1.896629 0.0201039783
# Srebf1    1.998829 0.0101458824
# Stat1     1.834970 0.0296623042
# Tcf4      1.748261 0.0497133057
# Trp53     1.802538 0.0361340550






# dryList_BigThree$parameters |>
#   filter((!is.na(amp_CTRL) &
#             !is.na(amp_FAR) &
#             !is.na(amp_RF))) |>
#   filter(amp_CTRL != amp_FAR & amp_CTRL != amp_RF)