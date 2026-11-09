#' Generate and Save Plots for Rhythmic Gene Data
#'
#' @description This script generates several plots visualizing different aspects of rhythmic gene data. It uses data from a list called `dryList` and generates plots for AdLib and Day Fed conditions.
#'
#' @param dryList A list or object containing results from some analysis. It should have components called `parameters` and `results`, which are data frames with gene names as row names and various columns for amplitudes, phases, etc.
#'
#' @return The script generates and saves multiple plots as PNG files in the specified directory. Each plot visualizes different aspects of the rhythmic gene data.
#'
#' @export

library(tidyverse)

# Fig. A) Euler AdLib vs Day Fed ----
AL_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(2, 4, 5), ])
DF_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(3, 4, 5), ])
Euler <- ggvenn::ggvenn(
  list(`AL` = AL_rhythmic, `DF` = DF_rhythmic),
  fill_color = c("#007FFF", "#FF7F00"),
  text_size = 1.4,
  stroke_size = 0.2,
  auto_scale = T,
  set_name_size = 1
)
ggsave(
  filename = "1-Rhythmic_Euler.png",
  plot = Euler,
  path = "./output/graphs/",
  width = 4,
  height = 3,
  dpi = 600,
  units = "cm"
)


# Fig. B) Amplitude of rhythmic genes count ----
## The plot without the legend box
LogAmp <- dryList$results |>
  rownames_to_column("gene_name") |>
  select(gene_name, contains("amp_") & !contains("relamp_")) |>
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
  xlab("Amplitude (log2FC)") +
  ylab("log10(rhythmic genes)") +
  geom_line(linewidth = .5) +
  scale_color_manual(values = c("#007FFF", "#FF7F00")) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 8),
    axis.line = element_line(color = "black", linewidth = .7)
  )
ggsave(
  filename = "2_1-Rhythmic_Amplitude.svg",
  plot = LogAmp,
  path = "~/../Downloads/Nouvelles figures/",
  width = 3,
  height = 4,
  dpi = 600,
  units = "cm"
)


# Fig. C) Cumulative counts per hour in AL and Day Fed mice ----
Acrophase <- ggpubr::ggarrange(
  plotlist = list(
    # Ad Libitum
    ggplot(
      dryList$results |>
        rownames_to_column("gene_name") |>
        select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_RF)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_RF),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "RF")
        )) |>
        filter(rhythmic_in %in% c("CTRL", "both")),
      aes(x = phase_CTRL, y = after_stat(count), fill = rhythmic_in)
    ) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity"
      ) +
      scale_fill_manual(values = c("#007FFF", "grey70")) +
      ggtitle(label = "AdLib") +
      xlab(label = "Acrophase (ZT)") +
      ylab(label = "Gene count") +
      scale_x_continuous(
        breaks = seq(4, 24, 4),
        expand = c(0, 0),
        limits = c(0, 24)
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 400)) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 8)
      )
    ,
    # Day Trained Feeding
    ggplot(
      dryList$results |>
        rownames_to_column("gene_name") |>
        select(gene_name, contains("phase_")) |>
        filter(!is.na(phase_CTRL) |
                 !is.na(phase_RF)) |>
        mutate(rhythmic_in = ifelse(
          !is.na(phase_CTRL) & !is.na(phase_RF),
          "both",
          ifelse(!is.na(phase_CTRL), "CTRL", "RF")
        )) |>
        filter(rhythmic_in %in% c("RF", "both")),
      aes(x = phase_RF, y = after_stat(count), fill = rhythmic_in)
    ) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity"
      ) +
      scale_fill_manual(values = c("#FF7F00", "grey70")) +
      xlab(label = "acrophase (ZT)") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 400)) +
      ylab(label = "genes count") +
      scale_x_continuous(
        breaks = seq(4, 24, 4),
        expand = c(0, 0),
        limits = c(0, 24)
      ) +
      theme_classic() +
      ggtitle(label = "DF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 8)
      )
  )
)

ggsave(
  filename = "3-Rhythmic_cumulative_counts.png",
  plot = Acrophase,
  path = "./output/graphs/",
  width = 7,
  height = 4,
  dpi = 600,
  units = "cm"
)


# Fig. D) Phase shift histogram ----
Shift <- ggplot(
  dryList$parameters |>
    filter((!is.na(amp_CTRL) &
              !is.na(amp_RF))) |>
    mutate(phase_shift = ifelse(
      abs((phase_CTRL - phase_RF)) < 12, abs((phase_CTRL - phase_RF)), 24 - abs(phase_CTRL - phase_RF)
    )) |>
    mutate(phase_shift = ifelse(phase_shift == 0, NA, phase_shift), ),
  aes(x = phase_shift)
) +
  geom_histogram(bins = 12,
                 breaks = seq(0, 12, by = 1),
                 fill = "#FF7F00") +
  ggtitle(element_text("Phase shift induced\nby DF")) +
  theme_classic() +
  xlim(0, 13) +
  scale_x_continuous(breaks = seq(0, 12, 4)) +
  ylab(label = "Gene count") +
  xlab(label = "Phase shift") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8, hjust = 1),
    text = element_text(size = 8)
  )

ggsave(
  filename = "4-Phase_shift.png",
  plot = Shift,
  path = "./output/graphs/",
  width = 3,
  height = 4,
  dpi = 600,
  units = "cm"
)


# Fig. E) Gene Rhythms modelisation ----
# TODO centrer - réduire les valeurs d'expression
# Core clocks genes
FilterGenesTable <- read.csv2(file = "./data/FilterGenesTable.csv")
Genes2plot <- FilterGenesTable[FilterGenesTable$gene_name %in% c("Bmal1", "Clock", "Cry1", "Cry2", "Per2", "Nr1d1"), ]
dataGraph <- as.data.frame(log2(dryList$ncounts + 1))
colnames(dataGraph) <- paste(groupTable[colnames(dataGraph), "Diet"], groupTable[colnames(dataGraph), "ZT"], groupTable[colnames(dataGraph), "sample_name"], sep = ".")
Rhythms <- data.frame(ZT = rep(0:24, 2), diet = rep(c("CTRL", "RF"), each = 25)) |>
  left_join(
    dryList$parameters[rownames(dryList$parameters) %in% Genes2plot$gene_id, ] |>
      rownames_to_column("gene_id") |>
      select(-contains("model"), -contains("relamp")) |>
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
    dataGraph[rownames(dryList$parameters) %in% Genes2plot$gene_id, ] |>
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
  left_join(Genes2plot[, c("gene_id", "gene_name")], relationship = "many-to-many") |>
  ggplot(aes(x = ZT, y = Mdl)) +
  geom_pointrange(
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
  geom_line(aes(color = diet)) +
  scale_color_manual(values = c("#007FFF", "#FF7F00")) +
  facet_wrap(facets = vars(gene_name),
             ncol = 2,
             scale = "free_y") +
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
ggsave(
  filename = "5-Genes_rhythms.png",
  plot = Rhythms,
  path = "./output/graphs/",
  width = 7,
  height = 6,
  dpi = 600,
  units = "cm"
)


# Fig. F1) Mean amplitude per model of AdLib and Day Fed mice ----
MeanAmp <- dryList$parameters |>
  rownames_to_column("gene_id") |>
  select(gene_id, amp_CTRL, amp_RF, chosen_model) |>
  filter(chosen_model != 1) |>
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
  ylab(label = "Model") +
  geom_tile() +
  scale_fill_viridis_c(
    begin = 0.5,
    end = 1,
    direction = -1,
    option = "inferno"
  ) +
  coord_fixed() +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "right")
ggsave(
  filename = "6-Mean_amplitude_per_model.png",
  plot = MeanAmp,
  path = "./output/graphs/",
  width = 6,
  height = 7.5,
  dpi = 600,
  units = "cm"
)


# Fig. F2) Gene counts per model ----
model_size <- dryList$results |>
  group_by(chosen_model) |>
  filter(chosen_model != 1) |>
  summarise(model_size = n())
Size <- model_size |>
  filter(chosen_model != 1) |>
  mutate(model = factor(paste0("#", chosen_model), levels = rev(paste0("#", chosen_model)))) |>
  ggplot(aes(x = 1, y =  model, label = model_size)) +
  geom_text() +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank()
  ) +
  geom_hline(yintercept = 1:4 + 0.5) +
  ggtitle("Gene count\nper model")

ggsave(
  filename = "7-Gene_count_per_model.png",
  plot = Size,
  path = "./output/graphs/",
  width = 4.5,
  height = 5,
  dpi = 600,
  units = "cm"
)


# Fig. G) Heatmaps of DryR models ----
## Heatmaps function
Heatmaps_Graph <- function(dryList, model) {
  OrderModL <- dryList$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_RF, phase_CTRL)) |>
    row.names()
  
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "CTRL"] # get genes with the same order than ordermodL
  dataHM$RF   <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "RF"] # get genes with the same order than ordermodL
  
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "RF"), function(x) {
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

model_size <- dryList$results |>
  group_by(chosen_model) |>
  filter(chosen_model != 1) |>
  summarise(model_size = n())

Heatmaps <- ggpubr::ggarrange(
  plotlist = lapply(2:5, function(i) {
    Heatmaps_Graph(dryList = dryList, model = i)
  }),
  ncol = 1,
  heights = 2 * (model_size$model_size / max(model_size$model_size)),
  widths = 1
)

ggsave(
  filename = "8-Rhythmic_heatmaps_per_model.png",
  plot = Heatmaps,
  path = "./output/graphs/",
  width = 3,
  height = 10,
  dpi = 600,
  units = "cm"
)


# Fig. H) Gene regulation in DT fed mice per model ----
ModRhythm <- ggplot(
  dryList$results |>
    rownames_to_column("gene_name") |>
    select(gene_name, contains("phase_"), chosen_model) |>
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
  scale_fill_manual(values = c("#007FFF", "#FF7F00")) +
  ggtitle(label = "Cumulative counts\nper model\nof rhythmic genes") +
  xlab(label = "Acrophase (ZT)") +
  ylab(label = "Gene count") +
  scale_x_continuous(
    breaks = seq(4, 24, 4),
    expand = c(0, 0),
    limits = c(0, 24)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 1),
    text = element_text(size = 8),
    strip.background = element_blank()
  ) +
  facet_wrap(ncol = 1,
             facets = vars(chosen_model),
             scales = "free_y")
ggsave(
  filename = "9-Rhythmic_cumulative_count_per_model.png",
  plot = ModRhythm,
  path = "./output/graphs/",
  width = 3,
  height = 7,
  dpi = 600,
  units = "cm"
)
