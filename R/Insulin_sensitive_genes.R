library(tidyverse)
library(readr)

load("./CircanutrEnv.RData")

insulino_sensible <- read_delim(
  "./data/genes_IR_and_insulin_dependent.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)[, -1]

insulino_sensible <- insulino_sensible |>
  left_join(y = FilterGenesTable[, c("gene_name", "gene_id")]) |>
  left_join(y = dryList$parameters |>
              select(chosen_model) |>
              rownames_to_column("gene_id"))

# Gene per model
insulino_sensible |>
  group_by(chosen_model) |>
  summarise(n())

# chosen_model `n()`
# <int> <int>
# 1            1   112
# 2            2   144
# 3            3   105
# 4            4    43
# 5            5   313
# 6           NA    47


# Test hypergeometric
liver_Ron_Kahn <- read_delim(
  "./data/liver Ron Kahn.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

Enrichment_in_ins <- do.call(rbind, lapply(1:5, function(x) {
  Background <- intersect(unlist(as.vector(liver_Ron_Kahn[, 1])), rownames(dryList$parameters))
  genes_in_model <- row.names(dryList$parameters)[dryList$parameters$chosen_model == x]
  
  genes_in_model[genes_in_model %in% Background]
  
  m <- dim(insulino_sensible)[1] # the number of white balls in the urn
  n <- length(Background[!Background %in% insulino_sensible$gene_id]) # the number of black balls in the urn
  
  k <- length(genes_in_model) 	# the number of balls drawn from the urn ie nb of gene from GO cat in the RNAseq
  q <- length(insulino_sensible$gene_id[insulino_sensible$gene_id %in% genes_in_model])
  
  data.frame(
    model = x,
    model_size = k,
    sensitive = q,
    p_value = phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
  )
}))

# model         model_size sensitive      p_value
# 1     1       9052       112            1.000000e+00
# 2     2       3496       144            9.999999e-01
# 3     3       2429       105            9.998934e-01
# 4     4       1265        43            9.999838e-01
# 5     5       4104       313            7.830026e-09

IR_sensi_enrichment <- ggplot(Enrichment_in_ins |>
                                mutate(y.min = rev(c(0, 313, 356, 461, 605)), y.max = rev(c(
                                  313, 356, 461, 605, 718
                                )))) + geom_rect(aes(
                                  xmin = 1,
                                  xmax = 2,
                                  ymin = (y.min),
                                  ymax = (y.max),
                                  fill = factor(model)
                                )) +
  geom_text(aes(
    label = paste0(
      "model ",
      model,
      "\n p-value = ",
      format(p_value, scientific = T, digits = 3)
    ),
    x = 1.5,
    y = (y.max + y.min) / 2
  ))

ggsave(
  filename = "IR_sensi_enrichment.png",
  plot = graph,
  path = "./output/graphs/",
  width = 16,
  height = 15,
  dpi = 320,
  units = "cm"
)


# Pie chart of insulin_sensitive genes in models
png(
  filename = "./output/graphs/IR_sensi_pie.png",
  width = 16,
  height = 15,
  res = 320,
  units = "cm"
)

pie(
  Enrichment_in_ins$sensitive,
  labels = c("Arrhythmic", "Loss", "Gain", "Identical", "Reversed"),
  clockwise = T,
  col = c("darkgray", "brown3", "seagreen3", "royalblue", "orange")
)
dev.off()

# Chronogram
Chronogram <- ggpubr::ggarrange(
  plotlist = list(
    dryList$parameters[insulino_sensible$gene_id, ] |>
      ggplot(aes(x = phase_CTRL)) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity",
        fill = "#007FFF"
      ) +
      xlim(0, 24) +
      ylim(0, 50) +
      coord_polar() +
      scale_x_continuous(
        breaks = seq(-4, 28, by = 4),
        expand = c(0, 0),
        lim = c(0, 24)
      ) +
      ylab("Gene count") +
      xlab(label = "Ad libitum") +
      theme_minimal(),
    
    dryList$parameters[insulino_sensible$gene_id, ] |>
      ggplot(aes(x = phase_RF)) +
      geom_histogram(
        bins = 24,
        breaks = seq(0, 24, by = 1),
        position = "identity",
        fill = "#FF7F00"
      ) +
      xlim(0, 24) +
      ylim(0, 50) +
      coord_polar() +
      scale_x_continuous(
        breaks = seq(-4, 28, by = 4),
        expand = c(0, 0),
        lim = c(0, 24)
      ) +
      xlab(label = "Restricted") +
      ylab("Gene count") +
      theme_minimal()
  )
)

ggsave(
  filename = "IR_sensi_chronogram.png",
  plot = Chronogram,
  path = "./output/graphs/",
  width = 16,
  height = 15,
  dpi = 320,
  units = "cm"
)

# Heatmaps
Heatmaps_Graph <- function(dryList, model) {
  OrderModL <- dryList$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_RF, phase_CTRL)) |>
    row.names()
  
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "CTRL"] # get genes with the same order than ordermodL
  dataHM$RF   <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "RF"] # get genes with the same order than ordermodL
  
  dataHM$CTRL <- t(apply(dataHM$CTRL, 1, function(x) {
    (x - min(x)) / (max(x) - min(x))
  }))
  dataHM$RF <- t(apply(dataHM$RF, 1, function(x) {
    (x - min(x)) / (max(x) - min(x))
  }))
  
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "RF"), function(x) {
    pheatmap::pheatmap(
      dataHM[[x]],
      scale = "none",
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
  filename = "IR_sensi_heatmaps.png",
  plot = Heatmaps,
  path = "./output/graphs/",
  width = 3,
  height = 10,
  dpi = 320,
  units = "cm"
)
