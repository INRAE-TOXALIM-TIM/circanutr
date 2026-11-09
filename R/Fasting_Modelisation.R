# Fasting Modelization

library(BiocGenerics)
library(ggplot2)
library(ggraph)
library(igraph)
library(kuipeR)

load("~/circanutr/Fasting_Modelisation_Env.RData")

# Expression values table from normalized counts
countData.mat <- as.matrix(read.csv2(
  "./data/FilterCountTable.csv",
  stringsAsFactors = FALSE,
  row.names = 1
))

# Convert count data to an integer array
countData <- array(as.integer(countData.mat), dim = dim(countData.mat))
colnames(countData) <- colnames(countData.mat)
row.names(countData) <- row.names(countData.mat)

# Table containing the genes names
FilterGenesTable <- read.table("./data/FilterGenesTable.csv",
                               sep = ";",
                               header = TRUE)[, -1]

load("~/circanutr/DryList_BigThree_BGEnv.RData")

# DESeq2 Analysis ##############################################################

groupTable <- read.csv2("./data/GroupsTable.csv", row.names = 1)
diet <- as.factor(groupTable$Diet)

colData <- data.frame(row.names = colnames(countData),
                      diet,
                      time = as.factor(groupTable$ZT))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                      colData = colData,
                                      design =  ~ diet + time)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds.full <- DESeq2::DESeq(
  dds,
  full = BiocGenerics::design(dds),
  betaPrior = F,
  fitType = "parametric",
  test = "Wald",
  parallel = T,
  quiet = T
)

cpm <- BiocGenerics::counts(dds.full, normalized = TRUE)
cpm <- cpm[, diet == "F" |
             (diet == "CTRL" & as.factor(groupTable$ZT) == "0")]
time <- groupTable$ZT[diet == "F" |
                        (diet == "CTRL" &
                           groupTable$ZT == "0")]


# data.lm <- RNAseq$WS[unique(unlist(maineffects)), sapply(RNAseq$WS, is.numeric)]
# data.lm <- data.lm + min(data.lm) + 1

formulas <- list(
  circ = x ~ cos(2 * pi * time / 24) + sin(2 * pi * time / 24),
  flat = x ~ 1,
  lin  = x ~ time,
  exp  = x ~ exp(time),
  log  = x ~ log(time + 1),
  quad = x ~ time + I(time^2)
)

lm_results <- lapply(formulas, function(f) {
  apply(cpm, 1, function(gene_expr) {
    df <- data.frame(x = gene_expr, time = time)
    lm(f, data = df)
  })
})

BIC_models <- do.call(cbind.data.frame, lapply(lm_results, function(x)
  sapply(x, BIC)))

BIC_models$best_model <- apply(BIC_models, 1, function(x) {
  colnames(BIC_models)[x == min(x)]
})

BIC_models$best_model[sapply(BIC_models$best_model, length) > 1] <- "flat"

models.direction <- do.call(cbind.data.frame, lapply(lm_results, function(x)
  sapply(x, function(mdl)
    ifelse(mdl$coefficients[2] < 0, 'down', 'up'))))

models.direction$quad <- sapply(lm_results$quad, function(mdl)
  ifelse(mdl$coefficients[3] > 0, 'down', 'up'))

models.direction$flat <- "flat"
models.direction$circ <- "circ"

models.R.sqr <- do.call(cbind.data.frame, lapply(lm_results, function(x)
  sapply(x, function(mdl)
    summary(mdl)$adj.r.squared)))

row.names(models.direction) <- names(lm_results$circ)

BIC_models$direction <- unlist(mapply(
  function(x, y) {
    res = models.direction[x, y]
    if (is.null(res)) {
      res = NA
    } else{
      res = res
    }
    return(res)
  },
  row.names(models.direction),
  BIC_models$best_model
))

BIC_models$adj.r.sqr <- unlist(mapply(function(x, y) {
  if (length(y) > 1) {
    res = NA
  } else{
    res = models.R.sqr[x, y]
    if (is.null(res)) {
      res = NA
    } else{
      res = res
    }
  }
  return(res)
}, row.names(models.R.sqr), BIC_models$best_model))

# Keeping genes with R sqr > 0.5
models_profiles_data <- BIC_models |>
  dplyr::filter(adj.r.sqr > 0.5) |>
  tibble::rownames_to_column('gene_id') |>
  dplyr::select(gene_id, best_model, direction) |>
  dplyr::left_join(as.data.frame(cpm) |>
                     tibble::rownames_to_column('gene_id')) |>
  tidyr::pivot_longer(!c(gene_id, best_model, direction),
                      values_to = 'cpm',
                      names_to = 'sample_name') |>
  dplyr::left_join(data.frame(time = groupTable$ZT, sample_name = groupTable$sample_name),
                   by = "sample_name")

# Keeping genes with R sqr > 0.5
genes_by_mdl <- BIC_models |>
  dplyr::filter(adj.r.sqr > 0.5) |>
  tibble::rownames_to_column('gene_id') |>
  dplyr::left_join(FilterGenesTable[, c('gene_id', 'gene_name')]) |>
  dplyr::group_by(best_model, direction) |>
  dplyr::summarize(gene_name = list(gene_name)) |>
  tidyr::unite(best_model, c(best_model, direction)) |>
  tidyr::unnest(gene_name)

genes_by_mdl <- split(genes_by_mdl$gene_name, genes_by_mdl$best_model)

predicted_models <- models_profiles_data |>
  dplyr::filter(best_model != "flat", best_model != "circ") |>
  dplyr::group_by(gene_id) |>
  dplyr::mutate(`z.score` = scale(cpm, center = TRUE, scale = TRUE)) |>
  dplyr::group_by(best_model, direction) |>
  dplyr::reframe({
    df <- data.frame(x = z.score, time = time)
    f <- formulas[[unlist(models_profiles_data$best_model)[1]]]
    mdl <- lm(f, data = df)
    data.frame(time      = time,
               z.score   = z.score,
               predicted = predict(mdl))
  }) |>
  tidyr::unite('group', c(best_model, direction)) |>
  dplyr::select(-z.score) |>
  distinct()

## Non-rhythmic expression models ####
predicted_profile_per_gene <- models_profiles_data |>
  dplyr::filter(best_model != "flat", best_model != "circ") |>
  dplyr::group_by(gene_id) |>
  dplyr::mutate(`z.score` = scale(cpm, center = TRUE, scale = TRUE)) |>
  dplyr::group_by(gene_id, time) |>
  dplyr::summarise(
    z.score = mean(z.score),
    best_model = unlist(unique(best_model)),
    direction = unique(direction),
    .groups = "keep"
  ) |>
  tidyr::unite('group', c(best_model, direction)) |>
  ggplot(aes(x = time)) +
  geom_line(aes(y = z.score, group = gene_id),
            alpha = .2,
            color = 'steelblue') +
  geom_point(
    data = predicted_models,
    aes(y = predicted),
    # fun = mean,
    # geom = 'line',
    linewidth = 1.2,
    color = 'red'
  ) +
  facet_wrap(vars(group), ncol = 2) +
  ylab('Z score') +
  xlab('time (ZT)') +
  scale_x_continuous(breaks = c(0, 12, 24)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

predicted_profile_per_gene
predicted_models |>
  filter(group == "log_up") |>
  ggplot(aes(x = time, y = predicted)) + geom_point()


# With only genes having R sqr > 0.5
# ggsave(
#   predicted_profile_per_gene,
#   filename = "./output/figures/Transcription Factors/Fasting_models_Rsqr.svg",
#   width = 6,
#   height = 8,
#   units = "cm",
#   dpi = 600
# )

## Plot per gene ####
Fasting.model.graph <- function(gene) {
  table = models_profiles_data |>
    dplyr::filter(gene_id == get.gene_id(gene)) |>
    # dplyr::group_by(gene_id) |>
    # dplyr::mutate(`z.score` = scale(cpm, center = TRUE, scale = TRUE)) |>
    # dplyr::group_by(gene_id, time) |>
    # dplyr::summarise(
    #   z.score = mean(z.score),
    #   best_model = unlist(unique(best_model)),
    #   direction = unique(direction),
    #   .groups = "keep"
    # ) |>
    tidyr::unite('group', c(best_model, direction))
  
  model = unique(table$group)
  print(unique(table$group))
  
  if (grepl("circ", model)) {
    choice <- y ~ cos(2 * pi * x / 24) + sin(2 * pi * x / 24)
  }
  if (grepl("flat", model)) {
    choice <- y ~ 1
  }
  if (grepl("lin", model)) {
    choice <- y ~ x
  }
  if (grepl("exp", model)) {
    choice <- y ~ exp(x)
  }
  if (grepl("log", model)) {
    choice <- y ~ log(x + 1)
  }
  if (grepl("quad", model)) {
    choice <- y ~ x + I(x^2)
  }
  
  table |>
    ggplot(aes(x = time)) +
    stat_smooth(
      # data = predicted_models |> filter(group == model),
      aes(y = log2(cpm)),
      se = FALSE,
      method = "lm",
      formula = choice,
      # fun = mean,
      # geom = 'line',
      linewidth = 1,
      color = 'red'
    ) +
    stat_summary(aes(y = log2(cpm), group = gene_id),
                 color = 'black',
                 geom = "pointrange") +
    facet_wrap(vars(group), ncol = 2) +
    # ylab('Log2 normalized counts') +
    xlab('Time (h)') +
    scale_x_continuous(breaks = c(0, 12, 24)) +
    theme_bw(base_size = 10) +
    theme(
      aspect.ratio = 1,
      panel.grid.minor = element_blank(),
      text = element_text(size = 12),
      strip.text = element_blank(),
      plot.margin = margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1
      )
    )
}

Fasting.model.graph("Bmal1")
Fasting.model.graph("Clock")
Fasting.model.graph("Cry1")
Fasting.model.graph("Cry2")
Fasting.model.graph("Nr1d1")
Fasting.model.graph("Per1")
Fasting.model.graph("Per2")
Fasting.model.graph("Rora")


# Acrophase modelisation for rhythmic genes in fasting #########################

calc.acrophase <- function(fit, period = 24) {
  sss <- as.numeric(coef(fit)["sin(2 * pi * time/24)"])
  rrr <- as.numeric(coef(fit)["cos(2 * pi * time/24)"])
  # para a et b retrouvés dans dryR
  
  if (sss > 0  & rrr > 0) {
    p <- atan(abs(sss / rrr))
  } else if (sss > 0  & rrr < 0) {
    p <- pi - (atan(abs(sss / rrr)))
  } else if (sss < 0  & rrr > 0) {
    p <- 2 * pi - (atan(abs(sss / rrr)))
  } else if (sss < 0  & rrr < 0) {
    p <- pi + (atan(abs(sss / rrr)))
  }
  p * period / (2 * pi)
}


# Circadian genes during fasting
fast_circ <- models_profiles_data[models_profiles_data$best_model == "circ", -c(2:4)]

plot(x = fast_circ$time[fast_circ$gene_id == "ENSMUSG00000000056"], y = fast_circ$cpm[fast_circ$gene_id == "ENSMUSG00000000056"])

fit <- lapply(unique(fast_circ$gene_id), function(gene) {
  cpm <- fast_circ$cpm[fast_circ$gene_id == gene]
  time <- fast_circ$time[fast_circ$gene_id == gene]
  fit <- lm(cpm ~ cos(2 * pi * time / 24) + sin(2 * pi * time / 24))
  return(fit)
})
names(fit) <- unique(fast_circ$gene_id)

fast_circ.pred <- lapply(fit, function(gene) {
  data.frame(time = 0:24,
             cpm = predict(gene, newdata = data.frame(time = 0:24)))
})

plot(x = fast_circ$time[fast_circ$gene_id == "ENSMUSG00000000056"], y = fast_circ$cpm[fast_circ$gene_id == "ENSMUSG00000000056"])
lines(x = fast_circ.pred$ENSMUSG00000000056[, "time"],
      y = fast_circ.pred$ENSMUSG00000000056[, "cpm"])

acro_fast <- sapply(fit, calc.acrophase)
acro_fast <- data.frame(phase_CTRL = dryList_BigThree$parameters$phase_CTRL[rownames(dryList_BigThree$parameters) %in% unique(fast_circ$gene_id)]
                        ,
                        phase_F = acro_fast)


## Fasting-induced shift from AL ####
ggplot(
  acro_fast |>
    dplyr::mutate(phase_shift = ifelse(
      abs((phase_CTRL - phase_F)) < 12, abs((phase_CTRL - phase_F)), 24 - abs(phase_CTRL - phase_F)
    )) |>
    dplyr::mutate(phase_shift = ifelse(phase_shift == 0, NA, phase_shift), ),
  aes(x = phase_shift)
) +
  geom_histogram(bins = 12,
                 breaks = seq(0, 12, by = 1),
                 fill = "#E0115F") +
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

# ggsave(
#   filename = "./output/figures/Transcription Factors/Fasting_Shift.svg",
#   width = 6,
#   height = 8,
#   units = "cm",
#   dpi = 600
# )



# Transcription Factors during Fasting #########################################
get.gene_name <- function(idx) {
  res <- sapply(idx, function(x) {
    FilterGenesTable$gene_name[FilterGenesTable$gene_id %in% x]
  })
  res <- unlist(res)
  res
}
# get.gene_name(c("ENSMUSG00000022385", "ENSMUSG00000022383"))

get.gene_name("ENSMUSG00000055116")

get.gene_id <- function(name) {
  res <- sapply(name, function(x) {
    FilterGenesTable$gene_id[FilterGenesTable$gene_name == x]
  })
  res <- unlist(res)
  res
}
# get.gene_id(TLINK_TF$gene[TLINK_TF$term == "Ppara"])


TLINK_TF <- clusterProfiler::read.gmt(
  "~/circanutr/data/TFLink_Mus_musculus_interactions_SS_GMT_proteinName_v1.0.gmt"
)
# Column "term" correspond to the TF
# Column "gene" correspond to the TF targets genes


## TF identifiers index
TLINK_TF_table <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = as.character(TLINK_TF$term),
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "ALIAS"
)

## Genes list
Rtm_fast <- data.frame(gene_id = FilterGenesTable$gene_id, gene_name = FilterGenesTable$gene_name) |>
  unique() |>
  dplyr::filter(gene_id %in% fast_circ$gene_id)
# 720 genes staying rhythmic under 24h fasting

Rtm_fast_mdl14 <- data.frame(gene_id = FilterGenesTable$gene_id, gene_name = FilterGenesTable$gene_name) |>
  unique() |>
  filter(gene_id %in% rownames(dryList_BigThree$parameters)[dryList_BigThree$parameters$chosen_model == 14])
# 3386 genes staying rhythmic under 24h fasting that Shift under DF but not NF


## TF lists
TF_rtm_fast <- data.frame(gene_id = TLINK_TF_table$ENSEMBL, gene_name = TLINK_TF_table$SYMBOL) |>
  unique() |>
  dplyr::filter(gene_id %in% fast_circ$gene_id)
# 46 TF staying rhythmic under 24h fasting

TF_Bmal_like <- TF_rtm_fast |>
  filter(gene_id %in% rownames(dryList_BigThree$parameters)[dryList_BigThree$parameters$chosen_model == 14])
# 29 TF staying rhythmic under 24h fasting that Shift under DF but not NF

# write.csv(x = TF_rtm_fast$gene_name,
#           file = "./data/TF_rtm_fast.csv",
#           quote = FALSE)
# write.csv(x = TF_Bmal_like$gene_name,
#           file = "./data/TF_Bmal_like.csv",
#           quote = FALSE)

## Venn ####
Bmal_like_TF.venn <- list(
  `Shift in DayFed` = TF_Bmal_like$gene_id,
  `Transcriptions Factors` = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% unique(TLINK_TF_table$ENSEMBL)],
  `Rhythmic in Fasting` = Rtm_fast$gene_id
)


venn::venn(Bmal_like_TF.venn, ilabels = "counts")



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


## Bmal1 ####
Bmal_like_TF.phase <- data.frame(
  gene_id = TF_Bmal_like$gene_id,
  gene_name = TF_Bmal_like$gene_name,
  AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Bmal_like_TF.venn, intersect), "phase_CTRL"],
  NightFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Bmal_like_TF.venn, intersect), "phase_FAR"],
  DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Bmal_like_TF.venn, intersect), "phase_RF"]
) |>
  rowwise() |>
  mutate(
    Shift_Bmal1_AL = min(c(abs(AdLib - 21.8), abs(21.8 - AdLib))),
    Shift_Bmal1_NF = min(c(abs(NightFed - 22.1), abs(22.1 - NightFed))),
    Shift_Bmal1_DF = min(c(abs(DayFed - 10.3), abs(10.3 - DayFed)))
  ) |>
  arrange(Shift_Bmal1_AL)


### AL
Bmal_like_TF.Bmal_shift_AL <- data.frame(to = Bmal_like_TF.phase$gene_name,
                                         Shift_Bmal1_AL = Bmal_like_TF.phase$Shift_Bmal1_AL) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Bmal1_AL > 12, 24 - Shift_Bmal1_AL, Shift_Bmal1_AL))


gatc_plot(Bmal_like_TF.Bmal_shift_AL,
          'Bmal1',
          nudge = -0.1,
          edge.width = 20)


### NF
Bmal_like_TF.Bmal_shift_NF <- data.frame(to = Bmal_like_TF.phase$gene_name,
                                         Shift_Bmal1_NF = Bmal_like_TF.phase$Shift_Bmal1_NF) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Bmal1_NF > 12, 24 - Shift_Bmal1_NF, Shift_Bmal1_NF))


gatc_plot(Bmal_like_TF.Bmal_shift_NF,
          'Bmal1',
          nudge = -0.1,
          edge.width = 20)


### DF
Bmal_like_TF.Bmal_shift_DF <- data.frame(to = Bmal_like_TF.phase$gene_name,
                                         Shift_Bmal1_DF = Bmal_like_TF.phase$Shift_Bmal1_DF) |>
  mutate(from = "Bmal1", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Bmal1_DF > 12, 24 - Shift_Bmal1_DF, Shift_Bmal1_DF))


gatc_plot(Bmal_like_TF.Bmal_shift_DF,
          'Bmal1',
          nudge = -0.1,
          edge.width = 20)





## Clock ####
TF_Clock_like <- TF_rtm_fast |>
  filter(gene_id %in% rownames(dryList_BigThree$parameters)[dryList_BigThree$parameters$chosen_model == 14])
# 29 TF staying rhythmic under 24h fasting that Shift under DF but not NF

Clock_like_TF.venn <- list(
  `Shift in DayFed` = TF_Clock_like$gene_id,
  `Transcriptions Factors` = FilterGenesTable$gene_id[FilterGenesTable$gene_id %in% unique(TLINK_TF_table$ENSEMBL)],
  `Rhythmic in Fasting` = Rtm_fast$gene_id
)


Clock_like_TF.phase <- data.frame(
  gene_id = TF_Clock_like$gene_id,
  gene_name = TF_Clock_like$gene_name,
  AdLib = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Clock_like_TF.venn, intersect), "phase_CTRL"],
  NightFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Clock_like_TF.venn, intersect), "phase_FAR"],
  DayFed = dryList_BigThree$parameters[FilterGenesTable$gene_id %in% purrr::reduce(Clock_like_TF.venn, intersect), "phase_RF"]
) |>
  rowwise() |>
  mutate(
    Shift_Clock_AL = min(c(abs(AdLib - 21.8), abs(21.8 - AdLib))),
    Shift_Clock_NF = min(c(abs(NightFed - 22.1), abs(22.1 - NightFed))),
    Shift_Clock_DF = min(c(abs(DayFed - 10.3), abs(10.3 - DayFed)))
  ) |>
  arrange(Shift_Clock_AL)


### AL
Clock_like_TF.Clock_shift_AL <- data.frame(to = Clock_like_TF.phase$gene_name,
                                           Shift_Clock_AL = Clock_like_TF.phase$Shift_Clock_AL) |>
  mutate(from = "Clock", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Clock_AL > 12, 24 - Shift_Clock_AL, Shift_Clock_AL))



gatc_plot(Clock_like_TF.Clock_shift_AL,
          'Clock',
          nudge = -0.1,
          edge.width = 20)


### NF
Clock_like_TF.Clock_shift_NF <- data.frame(to = Clock_like_TF.phase$gene_name,
                                           Shift_Clock_NF = Clock_like_TF.phase$Shift_Clock_NF) |>
  mutate(from = "Clock", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Clock_NF > 12, 24 - Shift_Clock_NF, Shift_Clock_NF))


gatc_plot(Clock_like_TF.Clock_shift_NF,
          'Clock',
          nudge = -0.1,
          edge.width = 20)


### DF
Clock_like_TF.Clock_shift_DF <- data.frame(to = Clock_like_TF.phase$gene_name,
                                           Shift_Clock_DF = Clock_like_TF.phase$Shift_Clock_DF) |>
  mutate(from = "Clock", .before = to) |>
  mutate(shift_shorter = ifelse(Shift_Clock_DF > 12, 24 - Shift_Clock_DF, Shift_Clock_DF))

gatc_plot(Clock_like_TF.Clock_shift_DF,
          'Clock',
          nudge = -0.1,
          edge.width = 20)





# Time of peak expression for rhythmics mouse transcription\nfactors in AdLib, Night feeding, and Day feeding diets

# ggsave(
#   './output/Chronoplot_TF.png',
#   width = 10,
#   height = 10,
#   units = 'cm'
# )








# TODO GO enrichment of TF targets ####

# save.image("~/circanutr/Fasting_Modelisation_Env.RData")
