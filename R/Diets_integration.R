# Setting up ----
library(tidyverse)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)

source("./R/fun/Fasted_Enrichment.R")
source("./R/fun/Do_network.R")

load("./All_Diets.RData")



# Selecting genes with same intercept in both Ad Lib dryList ----
## Fasting ----
Robust.Fasting <- rownames(dryList$results)[(is.na(dryList$results$a_CTRL) &
                                               is.na(dryList.Fasting$results$a_CTRL)) |
                                              (!is.na(dryList$results$a_CTRL) &
                                                 !is.na(dryList.Fasting$results$a_CTRL))]


## Refed ----
Robust.Refed <- rownames(dryList$results)[(is.na(dryList$results$a_CTRL) &
                                             is.na(dryList.Refed$results$a_CTRL)) |
                                            (!is.na(dryList$results$a_CTRL) &
                                               !is.na(dryList.Refed$results$a_CTRL))]

# Integrating non-circadian diets ----
## Fasting ----
Table.Fasting <- dryList$parameters |>
  merge(Fasting.lm, by = 0, all.x = TRUE) |>
  mutate(Amp_FA = dryList.Fasting$parameters$amp_F) |>
  filter(Row.names %in% Robust.Fasting) |>
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


# Table.Fasting |>
#   group_by(chosen_model, best_model) |>
#   count(comment) |>
#   mutate(Percent = n / 20346 * 100)


## Refed ----
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

# Table.Refed |>
#   group_by(chosen_model, best_model) |>
#   count(comment) |>
#   mutate(Percent = n / 20346 * 100)



# Selecting models ----
## Fasting ----
Models.Fasting <- Table.Fasting |>
  summarize(
    Percent = n() / mean(DF_count) * 100,
    Rsqr = mean(Rsqr, na.rm = T),
    Mean_amp = mean(Amp_FA, na.rm = T)
  )

## Refed ----
Models.Refed <- Table.Refed |>
  summarize(
    Percent = n() / mean(DF_count) * 100,
    Rsqr = mean(Rsqr, na.rm = T),
    Mean_amp = mean(Amp_FR, na.rm = T)
  )


# Plots of models by Day Fed profiles ----
## Fasting ----
ggplot(
  Models.Fasting |>
    filter(Rsqr > 0.35 | is.na(Rsqr)) |>
    group_by(chosen_model, best_model) |> 
    mutate(comment = factor(
      comment, levels = c("down", "up", "shift", "sync", "gain", "loss", "flat"), ordered = F
    )) |>
    complete(comment),
  aes(
    x = 1,
    y = comment,
    fill = Mean_amp,
    size = Percent
  )
) +
  geom_point(shape = 21,
             color = "black",
             stroke = 0.2) +
  scale_fill_viridis_c(
    limits = c(0, 2),
    breaks = c(0, 0.5, 1, 1.5, 2),
    direction = -1,
    option = "inferno",
    name = "Mean\nAmplitude"
  ) +
  scale_size_continuous(limits = c(0, NA), range = c(1, 5)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  facet_wrap(nrow = 5,
             ncol = 1,
             facets = vars(chosen_model)) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 12)
  )

ggsave(
  filename = "./output/graphs/Models_Fasting.png",
  width = 6,
  height = 18,
  units = "cm",
  dpi = 600
)


## Refed ----
ggplot(
  Models.Refed |>
    filter(Rsqr > 0.35 | is.na(Rsqr)) |>
    mutate(comment = factor(
      comment, levels = c("down", "up", "shift", "sync", "gain", "loss", "flat"), ordered = F
    )) |>
    complete(comment),
  aes(
    x = 1,
    y = comment,
    fill = Mean_amp,
    size = Percent
  )
) +
  geom_point(shape = 21,
             color = "black",
             stroke = 0.2) +
  scale_fill_viridis_c(
    limits = c(0, 2),
    breaks = c(0, 0.5, 1, 1.5, 2),
    direction = -1,
    option = "inferno",
    name = "Mean\nAmplitude"
  ) +
  scale_size_continuous(limits = c(0, NA), range = c(1, 5)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  facet_wrap(nrow = 5,
             ncol = 1,
             facets = vars(chosen_model)) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 12)
  )

ggsave(
  filename = "./output/graphs/Models_Refed.png",
  width = 6,
  height = 18,
  units = "cm",
  dpi = 600
)



# Enrichment ----

## Interesting groups ----
## Syntax : Profile_FAFR.Model_DF <- Fasted.Enrichment(chosen_model, "comment")
Gain.1  <- Fasted_Enrichment(1, "gain")
Sync.2  <- Fasted_Enrichment(2, "sync")
Shift.2 <- Fasted_Enrichment(2, "shift")
Gain.3  <- Fasted_Enrichment(3, "gain")
Shift.4 <- Fasted_Enrichment(4, "shift")
Sync.5  <- Fasted_Enrichment(5, "sync")
Shift.5 <- Fasted_Enrichment(5, "shift")
Down.5  <- Fasted_Enrichment(5, "down")

## Venn ----
ggvenn::ggvenn(
  list(`FA` = Down.5$FA, `FR` = Down.5$FR),
  fill_color = c("#E0115F", "#67987E"),
  text_size = 3,
  stroke_size = 0.2,
  auto_scale = T,
  set_name_size = 0
)

ggsave(
  filename = "./output/graphs/Fasted_Enrichment_Venn_Down_5.png",
  width = 8,
  height = 6,
  units = "cm",
  dpi = 600
)


## Network ----
d <- GOSemSim::godata(org.Mm.eg.db::org.Mm.eg.db, ont = "BP")

do_network(
  threshold = 0.4,
  go_terms = head(Down.5[["enrichment"]][["intersect_FAFR"]]@result[["ID"]]),
  GO_data = d,
  go_size = head(Down.5[["enrichment"]][["intersect_FAFR"]]@result[["Count"]]),
  go_color = -log10(head(Down.5[["enrichment"]][["intersect_FAFR"]]@result[["p.adjust"]]))
)



# Rhythmic parameters on Refed ----
## Amplitude of rhythmic genes ----
ggplot(
  dryList.Refed$results |>
    rownames_to_column("gene_name") |>
    filter(gene_name %in% Table.Refed$Row.names[Table.Refed$best_model == "circadian"]) |>
    select(gene_name, contains("amp_") &
             !contains("relamp_")) |>
    pivot_longer(!gene_name, names_to = "Diet", values_to = "log2FC") |>
    group_by(Diet) |>
    filter(!is.na(log2FC)) |>
    arrange(log2FC) |>
    mutate(
      `number of genes` = order(rev(log2FC)),
      Diet = gsub("amp_", "", Diet)
    ),
  aes(
    x = log2FC,
    y = log10(`number of genes`),
    colour = Diet
  )
) +
  xlab("Amplitude (log2FC)") +
  ylab("log10(rhythmic genes)") +
  geom_line(linewidth = 0.2) +
  scale_color_manual(values = c("#007FFF", "#67987E")) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "plain", size = 7))

ggsave(
  filename = "./output/graphs/Refed_Amplitude.png",
  width = 3,
  height = 4,
  units = "cm",
  dpi = 600
)


## Phase shift of rhythmic genes ----
ggplot(
  dryList.Refed$parameters |>
    filter((!is.na(amp_CTRL) &
              !is.na(amp_FAR))) |>
    rownames_to_column("gene_name") |>
    filter(gene_name %in% Table.Refed$Row.names[Table.Refed$best_model == "circadian"]) |>
    mutate(phase_shift = ifelse(
      abs((phase_CTRL - phase_FAR)) < 12, abs((phase_CTRL - phase_FAR)), 24 - abs(phase_CTRL - phase_FAR)
    )) |>
    mutate(phase_shift = ifelse(phase_shift == 0, NA, phase_shift), ),
  aes(x = phase_shift)
) +
  geom_histogram(bins = 12,
                 breaks = seq(0, 12, by = 1),
                 fill = "#67987E") +
  ggtitle(element_text("Phase shift induced\nby Refeeding")) +
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
  filename = "./output/graphs/Refed_PhaseShift.png",
  width = 3,
  height = 4,
  units = "cm",
  dpi = 600
)
