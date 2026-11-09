# Find mouse transcription factors in circadian dataset


# Set up #######################################################################
library(org.Mm.eg.db)
library(patchwork)
library(igraph)
library(ggraph)
library(tidyverse)

# devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr/")
library(kuipeR)

load("~/circanutr/LiverTfEnv.RData")



# Datasets #####################################################################
# load("~/circanutr/CircanutrEnv.RData")

# GO term 0032922 : circadian regulation of gene expression
CircaGO_table <- GOfuncR::get_anno_genes(go_ids = "GO:0032922", database = "Mus.musculus")
# 74 TF

FilterGenesTable <- read_delim(
  "data/FilterGenesTable.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)[, -1]

load("~/circanutr/DryList_BigThree_BGEnv.RData")

rhythmic_TF <- data.frame(gene_id = FilterGenesTable$gene_id[!is.na(dryList$parameters$phase_CTRL)], gene = FilterGenesTable$gene_name[!is.na(dryList$parameters$phase_CTRL)]) |>
  base::merge(CircaGO_table, by = "gene")
# 50 TF

core_clock <- CircaGO_table |>
  subset(gene %in% c("Bmal1", "Clock", "Cry1", "Cry2", "Nr1d1", "Per1", "Per2", "Rora")) |>
  mutate(gene_name = gene, gene = NULL) |>
  left_join(dplyr::select(FilterGenesTable, gene_id, gene_name), by = "gene_name")

TTRUST_TF <- data.frame(gene_name = unique(read.table("~/circanutr/data/TTRUST_TF_mouse.tsv")[, "V1"])) |>
  left_join(y = FilterGenesTable[, c("gene_id", "gene_name")]) |>
  drop_na()
# n = 634

## TF identifiers index
TTRUST_TF_table <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = TTRUST_TF$gene_name,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "ALIAS"
) |>
  filter(ENSEMBL %in% rownames(dryList_BigThree$parameters))

# sum(TTRUST_TF_table$ENSEMBL %in% rownames(dryList$parameters))
# 652

# nrow(TTRUST_TF_table)
# 652

TLINK_TF <- clusterProfiler::read.gmt(
  "~/circanutr/data/TFLink_Mus_musculus_interactions_SS_GMT_proteinName_v1.0.gmt"
)

## TF identifiers index
TLINK_TF_table <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = as.character(TLINK_TF$term),
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "ALIAS"
) |>
  filter(ENSEMBL %in% rownames(dryList_BigThree$parameters))


TLINK_TF.phase.table <- data.frame(gene_name = unique(TLINK_TF_table$SYMBOL)) |>
  left_join(
    dryList_BigThree$parameters |>
      rownames_to_column("gene_id") |>
      left_join(FilterGenesTable[, c("gene_id", "gene_name")]) |>
      dplyr::select(gene_id, gene_name, phase_CTRL, phase_FAR, phase_RF)
  ) |>
  filter(!is.na(gene_id))

TLINK_TF.shift <- TLINK_TF.phase.table |>
  rowwise() |>
  mutate(shift = min(c(
    abs(phase_CTRL - phase_RF), abs(phase_RF - phase_CTRL)
  )))

TLINK_TF.shift.table <- sapply(TLINK_TF.phase.table$phase_CTRL, function(x)
{
  sapply(TLINK_TF.phase.table$phase_CTRL, function(y)
  {
    if (is.na(x) | is.na(y))
    {
      res <- NA
    }
    else {
      if (x < y) {
        res <- y - x
      } else {
        res <- x - y
      }
      if (res > 12) {
        res <- 24 - res
      }
    }
    res
  })
})
colnames(TLINK_TF.shift.table) <- TLINK_TF.phase.table$gene_name
rownames(TLINK_TF.shift.table) <- TLINK_TF.phase.table$gene_name


TLINK_TF.keep <- TLINK_TF |>
  filter(gene %in% FilterGenesTable$gene_name)


# Core clock genes #############################################################
Core_clock_dataGraph <- as.data.frame(log2(dryList.Refed$ncounts + 1))
colnames(Core_clock_dataGraph) <- paste(groupTable[colnames(Core_clock_dataGraph), "Diet"], groupTable[colnames(Core_clock_dataGraph), "ZT"], groupTable[colnames(Core_clock_dataGraph), "sample_name"], sep = ".")

## Night Fed ####
NF_TF <- data.frame(ZT = rep(0:24, 2), diet = rep(c("CTRL", "FAR"), each = 25)) |>
  left_join(
    dryList_BigThree$parameters[rownames(dryList_BigThree$parameters) %in% core_clock$gene_id, ] |>
      rownames_to_column("gene_id") |>
      dplyr::select(-contains("model"), -contains("relamp")) |>
      pivot_longer(
        cols = !gene_id,
        names_to = "nom_cols",
        values_to = "valeurs"
      ) |>
      separate(nom_cols, c("parameter", "diet")) |>
      mutate(parameter = gsub("mean", "MESOR", parameter)) |>
      pivot_wider(
        id_cols = c(gene_id, diet),
        names_from = parameter,
        values_from = valeurs
      ),
    relationship = "many-to-many"
  ) |>
  mutate(Mdl = as.numeric((
    MESOR + a * cos(2 * pi * ZT / 24) +
      b * sin(2 * pi * ZT / 24)
  ))) |>
  mutate(Mdl = ifelse(is.na(Mdl), MESOR, Mdl)) |>
  left_join(
    Core_clock_dataGraph[rownames(dryList.Refed$parameters) %in% core_clock$gene_id, ] |>
      rownames_to_column("gene_id") |>
      pivot_longer(
        cols = !gene_id,
        names_to = "nom_cols",
        values_to = "valeurs"
      )
    |>
      separate(nom_cols, into = c("diet", "ZT", "sample")) |>
      mutate(ZT = as.integer(ZT)) |>
      group_by(gene_id, diet, ZT) |>
      summarize(
        MEAN = mean(valeurs),
        SEM = sd(valeurs) / sqrt(n() - 1)
      ),
    relationship = "many-to-many"
  ) |>
  mutate(diet = factor(diet, levels = c("CTRL", "FAR"), ordered = T)) |>
  left_join(core_clock[, c("gene_id", "gene_name")], relationship = "many-to-many") |>
  ggplot(aes(x = ZT, y = Mdl)) +
  geom_pointrange(
    aes(
      x = ZT,
      y = MEAN,
      ymin = MEAN - SEM,
      ymax = MEAN + SEM,
      color = diet
    ),
    size = .5,
    linewidth = .25
  ) +
  geom_line(aes(color = diet)) +
  scale_color_manual(values = c("#67987E", "#007FFF")) +
  facet_wrap(facets = vars(gene_name),
             ncol = 2,
             scale = "free_y") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  scale_size_area(max_size = .5) +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )


## Fasted ####
# FA_TF <- data.frame(ZT = rep(0:24, 2), diet = rep(c("CTRL", "F"), each = 25)) |>
#   left_join(
#     dryList.Fasting$parameters[rownames(dryList.Fasting$parameters) %in% core_clock$gene_id, ] |>
#       rownames_to_column("gene_id") |>
#       dplyr::select(-contains("model"), -contains("relamp")) |>
#       pivot_longer(
#         cols = !gene_id,
#         names_to = "nom_cols",
#         values_to = "valeurs"
#       ) |>
#       separate(nom_cols, c("parameter", "diet")) |>
#       mutate(parameter = gsub("mean", "MESOR", parameter)) |>
#       pivot_wider(
#         id_cols = c(gene_id, diet),
#         names_from = parameter,
#         values_from = valeurs
#       ),
#     relationship = "many-to-many"
#   ) |>
#   mutate(Mdl = as.numeric((
#     MESOR + a * cos(2 * pi * ZT / 24) +
#       b * sin(2 * pi * ZT / 24)
#   ))) |>
#   mutate(Mdl = ifelse(is.na(Mdl), MESOR, Mdl)) |>
#   left_join(
#     Core_clock_dataGraph[rownames(dryList.Fasting$parameters) %in% core_clock$gene_id, ] |>
#       rownames_to_column("gene_id") |>
#       pivot_longer(
#         cols = !gene_id,
#         names_to = "nom_cols",
#         values_to = "valeurs"
#       )
#     |>
#       separate(nom_cols, into = c("diet", "ZT", "sample")) |>
#       mutate(ZT = as.integer(ZT)) |>
#       group_by(gene_id, diet, ZT) |>
#       summarize(
#         MEAN = mean(valeurs),
#         SEM = sd(valeurs) / sqrt(n() - 1)
#       ),
#     relationship = "many-to-many"
#   ) |>
#   mutate(diet = factor(diet, levels = c("CTRL", "F"), ordered = T)) |>
#   left_join(core_clock[, c("gene_id", "gene_name")], relationship = "many-to-many") |>
#   ggplot(aes(x = ZT, y = Mdl)) +
#   geom_pointrange(
#     aes(
#       x = ZT,
#       y = MEAN,
#       ymin = MEAN - SEM,
#       ymax = MEAN + SEM,
#       color = diet
#     ),
#     size = .1,
#     linewidth = .25
#   ) +
#   geom_line(aes(color = diet)) +
#   scale_color_manual(values = c("#67987E", "#E0115F")) +
#   facet_wrap(facets = vars(gene_name),
#              ncol = 2,
#              scale = "fixed") +
#   theme_classic() +
#   scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
#   scale_size_area(max_size = .5) +
#   theme(
#     text = element_text(size = 8),
#     legend.position = "none",
#     strip.text = element_text(size = 8),
#     axis.text = element_text(size = 8),
#     strip.background = element_blank(),
#     plot.margin = margin(
#       t = 0,
#       r = 0,
#       b = 0,
#       l = 0
#     ),
#     axis.title.y = element_blank()
#   )


## Day Fed ####
DF_TF <- data.frame(ZT = rep(0:24, 2), diet = rep(c("CTRL", "RF"), each = 25)) |>
  left_join(
    dryList_BigThree$parameters[rownames(dryList_BigThree$parameters) %in% core_clock$gene_id, ] |>
      rownames_to_column("gene_id") |>
      dplyr::select(-contains("model"), -contains("relamp")) |>
      pivot_longer(
        cols = !gene_id,
        names_to = "nom_cols",
        values_to = "valeurs"
      ) |>
      separate(nom_cols, c("parameter", "diet")) |>
      mutate(parameter = gsub("mean", "MESOR", parameter)) |>
      pivot_wider(
        id_cols = c(gene_id, diet),
        names_from = parameter,
        values_from = valeurs
      ),
    relationship = "many-to-many"
  ) |>
  mutate(Mdl = as.numeric((
    MESOR + a * cos(2 * pi * ZT / 24) +
      b * sin(2 * pi * ZT / 24)
  ))) |>
  mutate(Mdl = ifelse(is.na(Mdl), MESOR, Mdl)) |>
  left_join(
    Core_clock_dataGraph[rownames(dryList$parameters) %in% core_clock$gene_id, ] |>
      rownames_to_column("gene_id") |>
      pivot_longer(
        cols = !gene_id,
        names_to = "nom_cols",
        values_to = "valeurs"
      )
    |>
      separate(nom_cols, into = c("diet", "ZT", "sample")) |>
      mutate(ZT = as.integer(ZT)) |>
      group_by(gene_id, diet, ZT) |>
      summarize(
        MEAN = mean(valeurs),
        SEM = sd(valeurs) / sqrt(n() - 1)
      ),
    relationship = "many-to-many"
  ) |>
  mutate(diet = factor(diet, levels = c("CTRL", "RF"), ordered = T)) |>
  left_join(core_clock[, c("gene_id", "gene_name")], relationship = "many-to-many") |>
  ggplot(aes(x = ZT, y = Mdl)) +
  geom_pointrange(
    aes(
      x = ZT,
      y = MEAN,
      ymin = MEAN - SEM,
      ymax = MEAN + SEM,
      color = diet
    ),
    size = .5,
    linewidth = .5
  ) +
  geom_line(aes(color = diet)) +
  scale_color_manual(values = c("#67987E", "#FF7F00")) +
  facet_wrap(facets = vars(gene_name),
             ncol = 2,
             scale = "free_y") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  scale_size_area(max_size = .5) +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    axis.title.y = element_blank()
  )


## Saving plots ####
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/NightFed_CoreClock.png",
#   NF_TF,
#   height = 8,
#   width = 8,
#   units = "cm",
#   dpi = 600
# )
#
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Fasted_CoreClock.png",
#   FA_TF,
#   height = 8,
#   width = 8,
#   units = "cm",
#   dpi = 600
# )
#
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/DayFed_CoreClock.png",
#   DF_TF,
#   height = 8,
#   width = 8,
#   units = "cm",
#   dpi = 600
# )



# Genes of interest ############################################################
Plot_GOI <- function(dryList, GOI, diet1, diet2, color1, color2) {
  dataGraph <- as.data.frame(log2(dryList$ncounts + 1))
  colnames(dataGraph) <- paste(groupTable[colnames(dataGraph), "Diet"], groupTable[colnames(dataGraph), "ZT"], groupTable[colnames(dataGraph), "sample_name"], sep = ".")
  
  Model_line <- data.frame(ZT = rep(0:24, 2), diet = rep(c(diet1, diet2), each = 25)) |>
    left_join(
      dryList$parameters[FilterGenesTable$gene_name == GOI, ] |>
        rownames_to_column("gene_id") |>
        dplyr::select(-contains("model"), -contains("relamp")) |>
        pivot_longer(
          cols = !gene_id,
          names_to = "nom_cols",
          values_to = "valeurs"
        ) |>
        separate(nom_cols, c("parameter", "diet")) |>
        mutate(parameter = gsub("mean", "MESOR", parameter)) |>
        pivot_wider(
          id_cols = c(gene_id, diet),
          names_from = parameter,
          values_from = valeurs
        ),
      relationship = "many-to-many"
    ) |>
    mutate(Mdl = as.numeric((
      MESOR + a * cos(2 * pi * ZT / 24) +
        b * sin(2 * pi * ZT / 24)
    ))) |>
    mutate(Mdl = ifelse(is.na(Mdl), MESOR, Mdl))
  
  Mean_points <- dataGraph[rownames(dryList$parameters) %in% rownames(dryList$parameters)[FilterGenesTable$gene_name == GOI], ] |>
    rownames_to_column("gene_id") |>
    pivot_longer(cols = !gene_id,
                 names_to = "nom_cols",
                 values_to = "valeurs") |>
    separate(nom_cols,
             into = c("diet", "ZT", "sample"),
             sep = "\\.") |>
    mutate(ZT = as.integer(ZT)) |>
    group_by(gene_id, diet, ZT) |>
    summarize(MEAN = mean(valeurs),
              SEM = sd(valeurs) / sqrt(n() - 1)) |>
    mutate(diet = factor(diet, levels = c(diet1, diet2), ordered = T))
  
  ggplot(Model_line, aes(x = ZT, y = Mdl)) +
    geom_line(aes(color = diet)) +
    geom_pointrange(
      data = Mean_points,
      aes(
        x = ZT,
        y = MEAN,
        ymin = MEAN - SEM,
        ymax = MEAN + SEM,
        color = diet
      ),
      size = .1,
      linewidth = .25
    ) +
    scale_color_manual(values = c(color1, color2)) +
    theme_classic() +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
    scale_size_area(max_size = .5) +
    theme(
      text = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      strip.background = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.title.y = element_blank()
    )
}

correlated_genes_to_ref <- function(gene) {
  table_all = models_profiles_data |>
    # dplyr::filter(gene_id == get.gene_id(gene)) |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(`z.score` = scale(cpm, center = TRUE, scale = TRUE)) |>
    dplyr::group_by(gene_id, time) |>
    dplyr::summarise(
      z.score = mean(z.score),
      best_model = unlist(unique(best_model)),
      direction = unique(direction),
      .groups = "keep"
    ) |>
    tidyr::unite('group', c(best_model, direction))
  
  table_gene = models_profiles_data |>
    dplyr::filter(gene_id == get.gene_id(gene)) |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(`z.score` = scale(cpm, center = TRUE, scale = TRUE)) |>
    dplyr::group_by(gene_id, time) |>
    dplyr::summarise(
      z.score = mean(z.score),
      best_model = unlist(unique(best_model)),
      direction = unique(direction),
      .groups = "keep"
    ) |>
    tidyr::unite('group', c(best_model, direction))
  
  gene_phase_AL = as.double(TLINK_TF.shift[TLINK_TF.shift$gene_name == gene, "phase_CTRL"])
  gene_phase_NF = as.double(TLINK_TF.shift[TLINK_TF.shift$gene_name == gene, "phase_FAR"])
  gene_phase_DF = as.double(TLINK_TF.shift[TLINK_TF.shift$gene_name == gene, "phase_RF"])
  
  
  ## Behaviour-like genes ####
  BigThreeVenn <- list(
    `Same in DayFed` = rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model == dryList_BigThree$parameters[get.gene_id(gene), "chosen_model"], ]),
    `Same in Fasting` = rownames(dryList_BigThree$parameters)[FilterGenesTable$gene_name %in% unlist(genes_by_mdl[unique(table_gene$group)])],
    `Transcriptions Factors` = TLINK_TF_table$ENSEMBL
  )
  
  venn::venn(BigThreeVenn, ilabels = "counts")
  
  # FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect)]
  # 62
  
  
  ### Pepsi plot for Phase distance to Bmal1 ####
  Big3.phase.table <- data.frame(
    gene_name = FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect)],
    gene_id = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect)],
    AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect), "phase_CTRL"],
    NightFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect), "phase_FAR"],
    DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn, intersect), "phase_RF"]
  ) |>
    rowwise() |>
    dplyr::mutate(
      Shift_AL = min(c(
        abs(AdLib - gene_phase_AL), abs(gene_phase_AL - AdLib)
      )),
      Shift_NF = min(c(
        abs(NightFed - gene_phase_NF), abs(gene_phase_NF - NightFed)
      )),
      Shift_DF = min(c(
        abs(DayFed - gene_phase_DF), abs(gene_phase_DF - DayFed)
      ))
    ) |>
    dplyr::arrange(Shift_AL)
  
  
  
  #### AL ####
  Big3.phase.table_Shift_AL <- data.frame(to = Big3.phase.table$gene_name, Shift_AL = Big3.phase.table$Shift_AL) |>
    mutate(from = gene, .before = to) |>
    mutate(shift_shorter = ifelse(Shift_AL > 12, 24 - Shift_AL, Shift_AL))
  
  print(gatc_plot(
    Big3.phase.table_Shift_AL,
    gene,
    nudge = -0.1,
    edge.width = 20
  ))
  
  
  #### NF ####
  Big3.phase.table_Shift_NF <- data.frame(to = Big3.phase.table$gene_name, Shift_NF = Big3.phase.table$Shift_NF) |>
    mutate(from = gene, .before = to) |>
    mutate(shift_shorter = ifelse(Shift_NF > 12, 24 - Shift_NF, Shift_NF))
  
  print(gatc_plot(
    Big3.phase.table_Shift_NF,
    gene,
    nudge = -0.1,
    edge.width = 20
  ))
  
  
  #### DF ####
  Big3.phase.table_Shift_DF <- data.frame(to = Big3.phase.table$gene_name, Shift_DF = Big3.phase.table$Shift_DF) |>
    mutate(from = gene, .before = to) |>
    mutate(shift_shorter = ifelse(Shift_DF > 12, 24 - Shift_DF, Shift_DF))
  
  print(gatc_plot(
    Big3.phase.table_Shift_DF,
    gene,
    nudge = -0.1,
    edge.width = 20
  ))
}

## Ppara ####
Plot_GOI(dryList.Refed, "Ppara", "CTRL", "FAR", "#67987E", "#007FFF")
Plot_GOI(dryList.Fasting, "Ppara", "CTRL", "F", "#67987E", "#E0115F")
Plot_GOI(dryList, "Ppara", "CTRL", "RF", "#67987E", "#FF7F00")

correlated_genes_to_ref("Ppara")

# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Ppara_NightFed.png",
#   Ppara_NF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Ppara_Fasted.png",
#   Ppara_FA,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Ppara_DayFed.png",
#   Ppara_DF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )


## Srebf1 ####
Plot_GOI(dryList.Refed, "Srebf1", "CTRL", "FAR", "#67987E", "#007FFF")
Plot_GOI(dryList.Fasting, "Srebf1", "CTRL", "F", "#67987E", "#E0115F")
Plot_GOI(dryList, "Srebf1", "CTRL", "RF", "#67987E", "#FF7F00")

correlated_genes_to_ref("Srebf1")

# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp1_NightFed.png",
#   Srebp1_NF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp1_Fasted.png",
#   Srebp1_FA,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp1_DayFed.png",
#   Srebp1_DF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )


## Srefp 2 ####
Plot_GOI(dryList.Refed, "Srebf2", "CTRL", "FAR", "#67987E", "#007FFF")
Plot_GOI(dryList.Fasting, "Srebf2", "CTRL", "F", "#67987E", "#E0115F")
Plot_GOI(dryList, "Srebf2", "CTRL", "RF", "#67987E", "#FF7F00")

correlated_genes_to_ref("Srebf2")

# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp2_NightFed.png",
#   Srebp2_NF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp2_Fasted.png",
#   Srebp2_FA,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Srebp2_DayFed.png",
#   Srebp2_DF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )


## Chrebp Alpha ####
Plot_GOI(dryList.Refed, "Mlxipl", "CTRL", "FAR", "#67987E", "#007FFF")
Plot_GOI(dryList.Fasting, "Mlxipl", "CTRL", "F", "#67987E", "#E0115F")
Plot_GOI(dryList, "Mlxipl", "CTRL", "RF", "#67987E", "#FF7F00")

correlated_genes_to_ref("Mlxipl")


# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Chrebpa_NightFed.png",
#   Chrebpa_NF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Chrebpa_Fasted.png",
#   Chrebpa_FA,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )
# ggsave(
#   "~/circanutr/output/figures/Transcription Factors/Chrebpa_DayFed.png",
#   Chrebpa_DF,
#   height = 10,
#   width = 10,
#   units = "cm",
#   dpi = 600
# )



# DryList 3 diets : AL, NF, DF ####
# load("~/circanutr/DryList_BigThree_BGEnv.RData")

## Chronoplot TF ####

# Check for duplicated genes
# TLINK_TF.phase.table |>
#   group_by(gene_name) |>
#   summarize(count = n()) |>
#   filter(count >1)


chronoplot(
  circ.data = TLINK_TF.phase.table[, -1],
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
  axis.pos = "center",
  title = "Transcriptions factors"
)

# Time of peak expression for rhythmics mouse transcription\nfactors in AdLib, Night feeding, and Day feeding diets

# ggsave(
#   './output/Chronoplot_TF.png',
#   width = 10,
#   height = 10,
#   units = 'cm'
# )



### Ppara - P53 ####

get.gene_name <- function(idx) {
  res <- sapply(idx, function(x) {
    FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% x]
  })
  res <- unlist(res)
  res
}
# get.gene_name(c("ENSMUSG00000022385", "ENSMUSG00000022383"))

get.gene_id <- function(name) {
  res <- sapply(name, function(x) {
    FilterGenesTable$gene_id[FilterGenesTable$gene_name == x]
  })
  res <- unlist(res)
  res
}
# get.gene_id(TLINK_TF$gene[TLINK_TF$term == "Ppara"])


#### Ppara
TLINK_TF.phase.table <- data.frame(
  gene_name = TLINK_TF.keep$gene[TLINK_TF.keep$term == "Ppara"],
  AdLib = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Ppara"]), "phase_CTRL"],
  DayFed = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Ppara"]), "phase_RF"],
  NightFed = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Ppara"]), "phase_FAR"]
)

chronoplot(
  circ.data = TLINK_TF.phase.table,
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
  axis.pos = "center",
  title = "Ppara"
)


#### P53
TLINK_TF.phase.table <- data.frame(
  gene_name = TLINK_TF.keep$gene[TLINK_TF.keep$term == "Tp53"],
  AdLib = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Tp53"]), "phase_CTRL"],
  DayFed = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Tp53"]), "phase_RF"],
  NightFed = dryList_BigThree$parameters[get.gene_id(TLINK_TF.keep$gene[TLINK_TF.keep$term == "Tp53"]), "phase_FAR"]
)

chronoplot(
  circ.data = TLINK_TF.phase.table,
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
  axis.pos = "center",
  title = "p53"
)


# Number of models
max(dryList_BigThree$results$chosen_model)
# 15 models

## Expression Models ####
# Number of genes per model
dryList_BigThree$results |>
  group_by(chosen_model) |>
  summarise(model_size = n()) # Count genes by model
#   chosen_model  model_size
# 1            1        6875  => flat                                   OUT
# 2            2         776  => flat in NF and DF                      OUT
# 3            3        1275  => Rhythmic in NF only                    OUT
# 4            4         730  => Rhythmic in DF only                    CHECK GENES LIST
# 5            5        3653  => Flat in DF                             MERGE WITH 6  (Per1-like)
# 6            6         343  => Shift in NF + Flat in DF               MERGE WITH 5  (Rora-like)
# 7            7         484  => Flat in NF                             OUT
# 8            8         171  => Flat in NF + Shift in DF               OUT
# 9            9         136  => Flat in AL                             CHECK GENES LIST
# 10           10        972  => Flat in AL + Shift in DF               CHECK GENES LIST 6 CF 9
# 11           11        501  => Sync for all                           CHECK GENES LIST
# 12           12         44  => Shift NF and DF but sync between them  CHECK GENES LIST
# 13           13        283  => Shift in NF                            CHECK GENES LIST - SHORT CHALLLENGE
# 14           14       3386  => Shift in DF                (Bmal1-like)
# 15           15        717  => Shift in NF + Shift in  DF (Per2-like)


FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Bmal1"]
# ENSMUSG00000055116
dryList_BigThree$parameters[get.gene_id("Srebf2"), "chosen_model"]
# 10

FilterGenesTable$gene_id[dryList_BigThree$parameters$chosen_model == 10]

get.gene_id("Rora")
get.gene_name("ENSMUSG00000032238")

dryR::dry_plot(dryList_BigThree, "ENSMUSG00000055116")
dryR::dry_plot(dryList_BigThree, get.gene_id("Rora"))
Fasting.model.graph("Nr1f1")

hist(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model %in% 15, "chosen_model_BICW"])


FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Rora"]
# ENSMUSG00000055116 : Bmal1              14    chosen_model_BICW = 0.9556691
# ENSMUSG00000029238 : Clock              14
# ENSMUSG00000020038 : Cry1               14
# ENSMUSG00000068742 : Cry2               14
# ENSMUSG00000020889 : Nr1d1, Rev-erba    14
# ENSMUSG00000021775 : Nr1d2, Rev-erbb    14
# ENSMUSG00000020893 : Per1               5
# ENSMUSG00000055866 : Per2               15
# ENSMUSG00000032238 : Rora               6

# ENSMUSG00000020538 : Srebf1             10
# ENSMUSG00000022463 : Srebf2             15

dryList_BigThree$parameters["ENSMUSG00000032238", "chosen_model"]
# 14, 5, 15, 6


# dryList.Fasting$parameters[get.gene_id("Srebf2"), "chosen_model"]
# 4 : synchronic (Bmal1, Per2)
# 5 : shift (Per1, Rora)
dryR::dry_plot(dryList.Fasting, "ENSMUSG00000055116")

dryList_BigThree$parameters[dryList.Fasting$parameters$chosen_model == 4, ]
# dryList.Fasting$parameters[dryList_BigThree$parameters$chosen_model == 14, ]



## Sunny plot ####
gatc_plot <- function(data.graph,
                      ref.gene,
                      nudge = 0,
                      edge.width = 12,
                      text.size = 8) {
  # create circular layout
  text.size = text.size * 0.38
  data.ggraph <- graph_from_data_frame(data.graph, directed = FALSE)
  
  layout <- create_layout(data.ggraph, layout = 'circle')
  
  # rotate Bmal1 to the top position
  bmal1_idx <- which(layout$name == ref.gene)
  bmal1_angle <- atan2(layout$y[bmal1_idx], layout$x[bmal1_idx])
  rotation <- pi / 2 - bmal1_angle
  
  x_rot <- layout$x * cos(rotation) - layout$y * sin(rotation)
  y_rot <- layout$x * sin(rotation) + layout$y * cos(rotation)
  layout$x <- x_rot
  layout$y <- y_rot
  
  # rotate the labels
  layout$angle <- atan2(layout$y, layout$x) * 180 / pi
  layout$hjust <- ifelse(layout$x > 0, 0, 1)
  layout$angle <- ifelse(layout$x < 0, layout$angle + 180, layout$angle)
  
  edge_colors <- setNames(data.graph$shift_shorter, data.graph$to)
  layout$shift_shorter <- edge_colors[layout$name]
  
  layout$shift_shorter[is.na(layout$shift_shorter)] <- NA
  
  ggraph(layout) +
    geom_edge_arc(
      aes(color = shift_shorter),
      edge_width = edge.width,
      strength = 0,
      alpha = 0.8,
      lineend = 'round'
    ) +
    scale_edge_color_viridis(
      alpha = 1,
      begin = 1,
      end = 0,
      limits = c(0, 12),
      name = "Shift (ZT)",
      guide = guide_edge_colorbar()
    ) +
    scale_edge_width(range = c(0.5, 2.5), guide = "none") +
    annotate(
      "path",
      x = cos(seq(0, 2 * pi, length.out = 200)),
      y = sin(seq(0, 2 * pi, length.out = 200)),
      color = "white",
      linewidth = 20
    ) +
    geom_node_text(
      aes(
        label = name,
        angle = angle,
        hjust = hjust,
        color = shift_shorter
      ),
      fontface = "bold",
      size = text.size,
      nudge_x = layout$x * nudge,
      nudge_y = layout$y * nudge
    ) +
    scale_color_viridis_c(
      begin = 1,
      end = 0,
      limits = c(0, 12),
      na.value = "grey",
      guide = "none"
    ) +
    geom_node_text(
      data = function(x)
        x[x$name == ref.gene, ],
      aes(
        label = name,
        angle = angle,
        hjust = hjust
      ),
      fontface = "bold",
      size = text.size,
      color = "red",
      nudge_x = layout$x[layout$name == ref.gene] * nudge,
      nudge_y = layout$y[layout$name == ref.gene] * nudge
    ) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "right") +
    coord_equal()
}



## Bmal1 Behaviour-like genes ####
BigThreeVenn_bmal1 <- list(
  `Shift in DayFed` = rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model == 14, ]),
  `Sync in Fasting` = rownames(dryList_BigThree$parameters[Fasting.lm$comment == "sync", ]),
  `Transcriptions Factors` = TTRUST_TF$gene_id
)

venn::venn(BigThreeVenn_bmal1, ilabels = "counts")

FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect)]
# 62


### Pepsi plot for Phase distance to Bmal1 ####
Big3_bmal1.phase.table <- data.frame(
  gene_name = FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect)],
  gene_id = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect)],
  AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect), "phase_CTRL"],
  NightFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect), "phase_FAR"],
  DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_bmal1, intersect), "phase_RF"]
) |>
  rowwise() |>
  mutate(
    Shift_Bmal1_AL = min(c(abs(AdLib - 21.8), abs(21.8 - AdLib))),
    Shift_Bmal1_NF = min(c(abs(NightFed - 22.1), abs(22.1 - NightFed))),
    Shift_Bmal1_DF = min(c(abs(DayFed - 10.3), abs(10.3 - DayFed)))
  ) |>
  arrange(Shift_Bmal1_AL)

# TTRUST_TF.shift[TTRUST_TF.shift$gene_name == "Bmal1", ]
# # gene_name   AdLib   DayFed  NightFed  Shift
# # Bmal1       21.8    10.3    22.1      11.5



#### AL ####
Big3.phase.table_BmalShift_AL <- data.frame(to = Big3_bmal1.phase.table$gene_name,
                                            Shift_Bmal1_AL = Big3_bmal1.phase.table$Shift_Bmal1_AL) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Bmal1_AL > 12, 24 - Shift_Bmal1_AL, Shift_Bmal1_AL))


Big3.phase.table_BmalShift_AL.graph <- graph_from_data_frame(Big3.phase.table_BmalShift_AL, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_BmalShift_AL.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 8,
        color = "NA",
        fill = "white",
        shape = 10,
        stroke = 1.5
      ) +
      geom_node_text(aes(label = name), fontface = "bold", size = 3) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



#### NF ####
Big3.phase.table_BmalShift_NF <- data.frame(to = Big3.phase.table$gene_name,
                                            Shift_Bmal1_NF = Big3.phase.table$Shift_Bmal1_NF) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Bmal1_NF > 12, 24 - Shift_Bmal1_NF, Shift_Bmal1_NF))

Big3.phase.table_BmalShift_NF.graph <- graph_from_data_frame(Big3.phase.table_BmalShift_NF, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_BmalShift_NF.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 9,
        color = "NA",
        fill = "white",
        shape = 9,
        stroke = 1.5
      ) +
      geom_node_text(
        aes(label = name),
        fontface = "bold",
        size = 3
        # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
        # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
      ) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



#### DF ####
Big3.phase.table_BmalShift_DF <- data.frame(to = Big3.phase.table$gene_name,
                                            Shift_Bmal1_DF = Big3.phase.table$Shift_Bmal1_DF) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Bmal1_DF > 12, 24 - Shift_Bmal1_DF, Shift_Bmal1_DF))

Big3.phase.table_BmalShift_DF.graph <- graph_from_data_frame(Big3.phase.table_BmalShift_DF, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_BmalShift_DF.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 9,
        color = "NA",
        fill = "white",
        shape = 9,
        stroke = 1.5
      ) +
      geom_node_text(
        aes(label = name),
        fontface = "bold",
        size = 3
        # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
        # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
      ) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



## Clock Behaviour-like genes ####
BigThreeVenn_clock <- list(
  `Shift in DayFed` = rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model == 14, ]),
  `Sync in Fasting` = rownames(dryList_BigThree$parameters[Fasting.lm$comment == "sync", ]),
  `Transcriptions Factors` = TTRUST_TF$gene_id
)

venn::venn(BigThreeVenn_clock, ilabels = "counts")

FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_clock, intersect)]
# 62

TTRUST_TF.shift[TTRUST_TF.shift$gene_name == "Clock", ]


### Pepsi plot for Phase distance to Bmal1 ####
Big3_clock.phase.table <- data.frame(
  gene_name = FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_clock, intersect)],
  gene_id = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_clock, intersect)],
  AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_clock, intersect), "phase_CTRL"],
  DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_clock, intersect), "phase_RF"]
) |>
  rowwise() |>
  mutate(Shift_Clock_AL = min(c(abs(AdLib - 23.3), abs(23.3 - AdLib))),
         Shift_Clock_DF = min(c(abs(DayFed - 12.6), abs(12.6 - DayFed)))) |>
  arrange(Shift_Clock_AL)

# Clock
# # gene_name   AdLib     DayFed    NightFed  Shift
# # Bmal1       23.3      12.6      23.3      11.5



#### AL ####
Big3.phase.table_ClockShift_AL <- data.frame(to = Big3_clock.phase.table$gene_name,
                                             Shift_Clock_AL = Big3_clock.phase.table$Shift_Clock_AL) |>
  mutate(from = "Clock", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Clock_AL > 12, 24 - Shift_Clock_AL, Shift_Clock_AL))


Big3.phase.table_ClockShift_AL.graph <- graph_from_data_frame(Big3.phase.table_ClockShift_AL, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_ClockShift_AL.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 8,
        color = "NA",
        fill = "white",
        shape = 10,
        stroke = 1.5
      ) +
      geom_node_text(aes(label = name), fontface = "bold", size = 3) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



#### NF ####
# Big3.phase.table_BmalShift_NF <- data.frame(to = Big3.phase.table$gene_name,
#                                             Shift_Bmal1_NF = Big3.phase.table$Shift_Bmal1_NF) |>
#   mutate(from = "Bmal1", .before = to) |>
#   mutate(Shift_shorter = ifelse(Shift_Bmal1_NF > 12, 24 - Shift_Bmal1_NF, Shift_Bmal1_NF))
#
# Big3.phase.table_BmalShift_NF.graph <- graph_from_data_frame(Big3.phase.table_BmalShift_NF, directed = FALSE)
#
# ggpubr::ggarrange(
#   plotlist = list(
#     ggraph(Big3.phase.table_BmalShift_NF.graph, layout = "circle") +
#       geom_edge_arc(
#         aes(color = Shift_shorter),
#         edge_width = 12,
#         strength = 0,
#         alpha = 0.8
#       ) +
#       scale_edge_color_viridis(
#         alpha = 1,
#         begin = 1,
#         end = 0,
#         limits = c(0, 12)
#         # name = "Shift (ZT)"
#       ) +
#       scale_edge_width(range = c(0.5, 2.5), guide = "none") +
#       geom_node_point(
#         size = 9,
#         color = "NA",
#         fill = "white",
#         shape = 9,
#         stroke = 1.5
#       ) +
#       geom_node_text(
#         aes(label = name),
#         fontface = "bold",
#         size = 3
#         # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
#         # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
#       ) +
#       theme_graph(base_family = "sans") +
#       theme(
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         plot.subtitle = element_text(
#           hjust = 0.5,
#           size = 10,
#           color = "grey50"
#         )
#       ) +
#       coord_flip(),
#     ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
#       scale_color_viridis_c(
#         begin = 1,
#         end = 0,
#         limits = c(0, 12)
#       ) +
#       geom_point(alpha = 0) +
#       theme_void()
#   ),
#   widths = c(1, 0.1)
# )



#### DF ####
# Big3.phase.table_BmalShift_DF <- data.frame(to = Big3.phase.table$gene_name,
#                                             Shift_Bmal1_DF = Big3.phase.table$Shift_Bmal1_DF) |>
#   mutate(from = "Bmal1", .before = to) |>
#   mutate(Shift_shorter = ifelse(Shift_Bmal1_DF > 12, 24 - Shift_Bmal1_DF, Shift_Bmal1_DF))
#
# Big3.phase.table_BmalShift_DF.graph <- graph_from_data_frame(Big3.phase.table_BmalShift_DF, directed = FALSE)
#
# ggpubr::ggarrange(
#   plotlist = list(
#     ggraph(Big3.phase.table_BmalShift_DF.graph, layout = "circle") +
#       geom_edge_arc(
#         aes(color = Shift_shorter),
#         edge_width = 12,
#         strength = 0,
#         alpha = 0.8
#       ) +
#       scale_edge_color_viridis(
#         alpha = 1,
#         begin = 1,
#         end = 0,
#         limits = c(0, 12)
#         # name = "Shift (ZT)"
#       ) +
#       scale_edge_width(range = c(0.5, 2.5), guide = "none") +
#       geom_node_point(
#         size = 9,
#         color = "NA",
#         fill = "white",
#         shape = 9,
#         stroke = 1.5
#       ) +
#       geom_node_text(
#         aes(label = name),
#         fontface = "bold",
#         size = 3
#         # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
#         # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
#       ) +
#       theme_graph(base_family = "sans") +
#       theme(
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         plot.subtitle = element_text(
#           hjust = 0.5,
#           size = 10,
#           color = "grey50"
#         )
#       ) +
#       coord_flip(),
#     ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
#       scale_color_viridis_c(
#         begin = 1,
#         end = 0,
#         limits = c(0, 12)
#       ) +
#       geom_point(alpha = 0) +
#       theme_void()
#   ),
#   widths = c(1, 0.1)
# )



## Per2 Behaviour-like genes ####

dryList_BigThree$parameters["ENSMUSG00000055866", "chosen_model"]
# model 15 in dryList 3 diets => Shift in NF + different shift in DF

dryList.Fasting$parameters["ENSMUSG00000055866", "chosen_model"]
# model 4 in fasting => Sync

BigThreeVenn_per2 <- list(
  `Shift in DayFed` = rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model == 15, ]),
  `Sync in Fasting` = rownames(dryList_BigThree$parameters[Fasting.lm$comment == "sync", ]),
  `Transcriptions Factors` = TTRUST_TF$gene_id
)

venn::venn(BigThreeVenn_per2, ilabels = "counts")

FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect)]
# 12
# Ppara Smarca2 Mxi1 Smad9 Ncor2 Foxp1
# Bhlhe40 Ing4 Smarca4 Hsf4 Nrip1 Per2


### Pepsi plot for Phase distance to Per2 ####
Big3_per2.phase.table <- data.frame(
  gene_name = FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect)],
  gene_id = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect)],
  AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect), "phase_CTRL"],
  NightFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect), "phase_FAR"],
  DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(BigThreeVenn_per2, intersect), "phase_RF"]
) |>
  rowwise() |>
  mutate(
    Shift_Per2_AL = min(c(abs(AdLib - 14.2), abs(14.2 - AdLib))),
    Shift_Per2_NF = min(c(abs(NightFed - 16.3), abs(16.3 - NightFed))),
    Shift_Per2_DF = min(c(abs(DayFed - 4.59), abs(4.59 - DayFed)))
  ) |>
  arrange(Shift_Per2_AL)

TTRUST_TF.shift[TTRUST_TF.shift$gene_name == "Per2", ]
# gene_name   AdLib   DayFed  NightFed  Shift
# Per2        14.2    4.59    16.3      9.61


#### AL ####
Big3.phase.table_Per2Shift_AL <- data.frame(to = Big3_per2.phase.table$gene_name,
                                            Shift_Per2_AL = Big3_per2.phase.table$Shift_Per2_AL) |>
  mutate(from = "Per2", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Per2_AL > 12, 24 - Shift_Per2_AL, Shift_Per2_AL))


Big3.phase.table_Per2Shift_AL.graph <- graph_from_data_frame(Big3.phase.table_Per2Shift_AL, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_Per2Shift_AL.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 8,
        color = "NA",
        fill = "white",
        shape = 10,
        stroke = 1.5
      ) +
      geom_node_text(aes(label = name), fontface = "bold", size = 3) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



#### NF ####
Big3.phase.table_Per2Shift_NF <- data.frame(to = Big3_per2.phase.table$gene_name,
                                            Shift_Per2_NF = Big3_per2.phase.table$Shift_Per2_NF) |>
  mutate(from = "Per2", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Per2_NF > 12, 24 - Shift_Per2_NF, Shift_Per2_NF))

Big3.phase.table_Per2Shift_NF.graph <- graph_from_data_frame(Big3.phase.table_Per2Shift_NF, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_Per2Shift_NF.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 9,
        color = "NA",
        fill = "white",
        shape = 9,
        stroke = 1.5
      ) +
      geom_node_text(
        aes(label = name),
        fontface = "bold",
        size = 3
        # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
        # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
      ) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)



#### DF ####
Big3.phase.table_Per2Shift_DF <- data.frame(to = Big3_per2.phase.table$gene_name,
                                            Shift_Per2_DF = Big3_per2.phase.table$Shift_Per2_DF) |>
  mutate(from = "Per2", .before = to) |>
  mutate(Shift_shorter = ifelse(Shift_Per2_DF > 12, 24 - Shift_Per2_DF, Shift_Per2_DF))

Big3.phase.table_Per2Shift_DF.graph <- graph_from_data_frame(Big3.phase.table_Per2Shift_DF, directed = FALSE)

ggpubr::ggarrange(
  plotlist = list(
    ggraph(Big3.phase.table_Per2Shift_DF.graph, layout = "circle") +
      geom_edge_arc(
        aes(color = Shift_shorter),
        edge_width = 12,
        strength = 0,
        alpha = 0.8
      ) +
      scale_edge_color_viridis(
        alpha = 1,
        begin = 1,
        end = 0,
        limits = c(0, 12)
        # name = "Shift (ZT)"
      ) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      geom_node_point(
        size = 9,
        color = "NA",
        fill = "white",
        shape = 9,
        stroke = 1.5
      ) +
      geom_node_text(
        aes(label = name),
        fontface = "bold",
        size = 3
        # nudge_y = c(rep(0.2, 32), rep(-0.2, 31)),
        # nudge_x = c(rep(-0.2, 32), rep(0.2, 31))
      ) +
      theme_graph(base_family = "sans") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10,
          color = "grey50"
        )
      ) +
      coord_flip(),
    ggplot(data = data.frame(x = 0, y = 0), aes(x, y, color = x)) +
      scale_color_viridis_c(
        begin = 1,
        end = 0,
        limits = c(0, 12)
      ) +
      geom_point(alpha = 0) +
      theme_void()
  ),
  widths = c(1, 0.1)
)








# Enrichments per models ####

dryList_BigThree$results |>
  group_by(chosen_model) |>
  summarise(model_size = n()) # Count genes by model
#   chosen_model  model_size
# 1            1        6875  => flat                                   OUT
# 2            2         776  => flat in NF and DF                      OUT
# 3            3        1275  => Rhythmic in NF only                    OUT
# 4            4         730  => Rhythmic in DF only                    CHECK GENES LIST
# 5            5        3653  => Flat in DF                             MERGE WITH 6  (Per1-like)
# 6            6         343  => Shift in NF + Flat in DF               MERGE WITH 5  (Rora-like)
# 7            7         484  => Flat in NF                             OUT
# 8            8         171  => Flat in NF + Shift in DF               OUT
# 9            9         136  => Flat in AL                             CHECK GENES LIST
# 10           10        972  => Flat in AL + Shift in DF               CHECK GENES LIST 6 CF 9
# 11           11        501  => Sync for all                           CHECK GENES LIST
# 12           12         44  => Shift NF and DF but sync between them  CHECK GENES LIST
# 13           13        283  => Shift in NF                            CHECK GENES LIST - SHORT CHALLLENGE
# 14           14       3386  => Shift in DF                (Bmal1-like)
# 15           15        717  => Shift in NF + Shift in  DF (Per2-like)


Mapped_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = FilterGenesTable$gene_id,
  columns = c('ENTREZID'),
  keytype = 'ENSEMBL'
)


Model_Enrich.list <- lapply(1:15, function(x) {
  Model.list <- rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model %in% x, ])
  
  Model_Enrich <- clusterProfiler::enrichGO(
    gene     = Model.list,
    OrgDb        = org.Mm.eg.db,
    universe = rownames(dryList_BigThree$parameters),
    keyType = "ENSEMBL",
    ont          = "BP",
    qvalueCutoff = 0.05
  )
  
  Model_Enrich.table <- data.frame(
    paths =  paste0(Model_Enrich@result$Description),
    count =  Model_Enrich@result$Count,
    qvalue = Model_Enrich@result$p.adjust
  )
  
  return(Model_Enrich.table)
})

Model_Enrich.list <- lapply(Model_Enrich.list, function(df) {
  df[df$qvalue < 0.05, ]
})

## Model 14 (Bmal1) ####
# library(BiocManager,clusterProfiler,org.Mm.eg.db,enrichplot)

Model14.list <- rownames(dryList_BigThree$parameters[dryList_BigThree$parameters$chosen_model %in% 14, ])

Model14_Enrich <- enrichGO(
  gene     = Model14.list,
  OrgDb        = org.Mm.eg.db,
  universe = rownames(dryList_BigThree$parameters),
  keyType = "ENSEMBL",
  ont          = "BP",
  qvalueCutoff = 0.05
)

Model14_Enrich@result |>
  filter(p.adjust < 0.05) |>
  nrow()
# 50

Model14_Enrich.table <- data.frame(
  paths =  paste0(Model14_Enrich@result$Description),
  count =  Model14_Enrich@result$Count,
  qvalue = Model14_Enrich@result$p.adjust
)



# Bmal1 correlation ####


# TTRUST_TF.shift[TTRUST_TF.shift$gene_name == "Bmal1", ]
# # gene_name   AdLib   DayFed  NightFed  Shift
# # Bmal1       21.8    10.3    22.1      11.5
#
# ### Correlate with Bmal1
# TTRUST_TF.shift |>
#   filter(AdLib > 20 & AdLib < 24 &
#            DayFed > 8 & DayFed < 12)
# # 14
#
# ### Anti-correlate with Bmal1
# TTRUST_TF.shift |>
#   filter(AdLib > 8 & AdLib < 12 &
#            DayFed > 20 & DayFed < 24)
# # 28
#
# TTRUST_TF.shift |>
#   filter(AdLib > 8 & AdLib < 12 &
#            DayFed > 8 & DayFed < 12)
# # 1
#
# ### Perte rythmicité en day fed
# TTRUST_TF.shift |>
#   filter(AdLib > 8 & AdLib < 12 &
#            is.na(DayFed))
# # 26
#
#
#
#
## Nb de gènes en fonction de leur heure de décalalge (AL vs DF et inverse)
## exclus ceux qui ne se décalent pas, donc que les rythmiques
# TTRUST_TF.phase.table |>
#   rowwise() |>
#   mutate(shift = min(c(abs(AdLib - DayFed), abs(DayFed - AdLib)))) |>
#   filter(shift != 0) |>
#   ggplot(aes(x = shift)) +
#   geom_histogram()



Bmal_cor <- TTRUST_TF.phase.table |>
  filter(AdLib > 20 & AdLib < 24 &
           DayFed > 8 & DayFed < 12 | gene_name == "Bmal1") |>
  rowwise() |>
  mutate(Shift = min(c(abs(AdLib - DayFed), abs(DayFed - AdLib)))) |>
  dplyr::select(gene_name, Shift) |>
  pivot_wider(names_from = "gene_name", values_from = Shift)

Bmal_cor.table <- data.frame(unlist(lapply(colnames(Bmal_cor), function(x)
{
  sapply(colnames(Bmal_cor), function(y)
  {
    paste(x, y)
  })
}))) |>
  separate(col = 1,
           into = c("from", "to"),
           sep = " ") |>
  mutate(Shift = as.vector(TTRUST_TF.shift.table[colnames(Bmal_cor), colnames(Bmal_cor)])) |>
  filter(from != to)




Bmal_anticor <- TTRUST_TF.phase.table |>
  filter(AdLib > 8 & AdLib < 12 &
           DayFed > 20 & DayFed < 24 | gene_name == "Bmal1") |>
  rowwise() |>
  mutate(Shift = min(c(abs(AdLib - DayFed), abs(DayFed - AdLib)))) |>
  dplyr::select(gene_name, Shift) |>
  pivot_wider(names_from = "gene_name", values_from = Shift)

# all.neg <- function(x) {
#   -1 * abs(x)
# }

Bmal_anticor.table <- data.frame(unlist(lapply(colnames(Bmal_anticor), function(x)
{
  sapply(colnames(Bmal_anticor), function(y)
  {
    paste(x, y)
  })
}))) |>
  separate(col = 1,
           into = c("from", "to"),
           sep = " ") |>
  # mutate(Shift = as.vector(all.neg(TTRUST_TF.shift.table[colnames(Bmal_anticor), colnames(Bmal_anticor)]))) |>
  mutate(Shift = as.vector(TTRUST_TF.shift.table[colnames(Bmal_anticor), colnames(Bmal_anticor)])) |>
  filter(from != to) |>
  mutate(sel_col = paste(from, to), reverse = paste(to, from))


# save.image("~/circanutr/LiverTfEnv.RData")

# # Bmal.table <- rbind(Bmal_cor.table, Bmal_anticor.table)
# #   sort(from)



# mutate(sel_col = sort(paste(from, to)))

Bmal.graph <- igraph::graph_from_data_frame(Bmal.table, directed = FALSE)
# mode = "lower", weighted = TRUE, diag = FALSE)


ggraph(Bmal.graph, layout = "circle") +
  geom_edge_arc(aes(color = Shift), strength = 0.1, alpha = 0.8) +
  scale_edge_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Correlation",
    na.value = "white"
  ) +
  scale_edge_width(range = c(0.5, 2.5), guide = "none") +
  geom_node_point(
    size = 8,
    color = "grey30",
    fill = "white",
    shape = 21,
    stroke = 1.5
  ) +
  geom_node_text(aes(label = name), fontface = "bold", size = 3.5) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "grey50"
    )
  )



# TF leading peak of expression ####

TLINK_Target_Table <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = as.character(TLINK_TF$gene),
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "ALIAS"
) |>
  dplyr::filter(ENSEMBL %in% rownames(dryList_BigThree$parameters)) |>
  dplyr::filter(SYMBOL %in% TLINK_TF.keep$gene)

# TLINK_TF.keep <- TLINK_TF |>
#   filter(gene %in% FilterGenesTable$gene_name)

TLINK_kuip <- TLINK_TF.keep |>
  rename("TF" = term, "SYMBOL" = gene)

TLINK_kuip <- TLINK_kuip |>
  left_join(TLINK_Target_Table, by = "SYMBOL")

TLINK_kuip <- TLINK_kuip |>
  left_join(dryList_BigThree$parameters |>
              rownames_to_column(var = "ENSEMBL"),
            by = "ENSEMBL")

# Nombre NA phase par diet
sum(is.na(TLINK_kuip$phase_CTRL)) # 25384
sum(is.na(TLINK_kuip$phase_RF))   # 30738
sum(is.na(TLINK_kuip$phase_FAR))  # 21067
dim(TLINK_kuip)                   # 58039

TLINK_kuip.per.group <- TLINK_kuip |>
  dplyr::select(c(TF, SYMBOL, ENSEMBL, phase_CTRL, phase_RF, phase_FAR)) |>
  distinct() |>
  pivot_longer(!c(TF, SYMBOL, ENSEMBL),
               names_to = "group",
               values_to = "phase") |>
  group_by(TF, group) |>
  summarize(kuiper_test(phase, limits = c(0, 24)))

Pval_per_diet <- TLINK_kuip.per.group |>
  filter(p.value < 0.05) |>
  group_by(group) |>
  summarize(n())
# group       n()
# phase_CTRL  23
# phase_FAR   32
# phase_RF    31

venn::venn(split(x = TLINK_kuip.per.group$TF[TLINK_kuip.per.group$p.value < 0.05], TLINK_kuip.per.group$group[TLINK_kuip$p.value < 0.05]),
           ilabels = "counts")
# Combien de FT sont des régulateurs de gènes coordonnés en pic d’expression ?

TF_regulators <- lapply(split(x = TLINK_kuip.per.group$TF[TLINK_kuip.per.group$p.value < 0.05], TLINK_kuip.per.group$group[TLINK_kuip.per.group$p.value < 0.05]), function(x) {
  sort(as.character(x))
})

data.frame(TF = sort(unique(as.character(
  TLINK_kuip.per.group$TF
))))

#TODO TLINK_kuip en TLINK_drylistParam
TLINK_kuip_signif <- TLINK_kuip |>
  dplyr::filter(TF %in% TLINK_kuip.per.group$TF[TLINK_kuip.per.group$p.value < 0.05]) |>
  dplyr::mutate(TF = as.character(TF)) |>
  dplyr::select(TF, SYMBOL, phase_CTRL, phase_RF) |>
  rename(Target = SYMBOL) |>
  distinct() |>
  group_by(TF) |>
  filter(n() > 10)

TLINK_kuip_signif.list <- split(TLINK_kuip_signif[, -1], f = TLINK_kuip_signif$TF)

# TODO chronoplot sur eux seulement
# lapply(test.list, function(x) {
TLINK_kuip_signif.res <- list()
for (i in seq_along(TLINK_kuip_signif.list)) {
  TLINK_kuip_signif.res[[i]] <- chronoplot(
    circ.data = TLINK_kuip_signif.list[[i]],
    distribution = NULL,
    limits = c(0, 24),
    n.class = 12,
    after.stat = TRUE,
    alternative = "greater",
    colors = c("cyan", "orange", "green3"),
    fill = c("cyan", "orange", "green3"),
    pch = NULL,
    axis.text.size = 8,
    points.size = 2,
    smooth.method = "four_points",
    smoothing = 10,
    w = 1 / 16,
    show.dist = FALSE,
    axis.pos = "center",
    title = names(TLINK_kuip_signif.list)[i]
  )
}
names(TLINK_kuip_signif.res) <- names(TLINK_kuip_signif.list)

# ggsave(
#   filename = paste("./Count/Chronoplot_", Get_go_name(x), ".png", sep = ""),
#   path = path,
#   width = 6,
#   height = 6,
#   dpi = 600,
#   units = "cm"
# )
# })


# TODO
# TF controlling peak of expr in chronoplot
# faire du kuiper par TF pour voir s'ils sortent signif sans ajustement pval
#     car tirage avec remise donc precision test depend de l'ensemble
# Faire ensuite chronoplots associés si signif


kuipeR::kuiper_test(TLINK_TF.phase.table, limits = c(0, 24))

# save.image("~/circanutr/LiverTfEnv.RData")
