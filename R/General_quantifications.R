# Loading packages ----
library(tidyverse)
library(ggplot2)
load("./CircanutrEnv.RData")



# Genes Numbers ----

# Corresponds to the raw abundance table of transcripts from the salmon alignment
tpm <- read.table(file = "./data/salmon.merged.gene_tpm.tsv", dec = ".", header = T)
rownames(tpm) <- tpm[, 1]   # Row names are gene_id
tpm <- tpm[, -1]

# A vector of all genes names from the ENSEMBL GRCm39 Mus musculus genome
# This is the "Mouse" set
gn_total <- as.vector(rownames(tpm))

# Number of mouse genes detected in "liver"Mouse" with >= 1 count in AdLib
OneCount <- read.table(file = "./data/OneCountTable.csv", sep = ";")
OneCount <- OneCount[-1, ]

FilterByCount <- data.frame(t(rbind(gn_total, OneCount$V2)))
FilterByCount$V2 <- as.logical(FilterByCount$V2)

# This is the "Liver" set
Liver <- subset(FilterByCount, V2)

# Number of mouse genes detected in the liver with >= 10 Counts in AdLib
# This is the "Detected" set
liver_detected <- as.matrix(read.csv2(
  "./data/FilterCountTable.csv",
  stringsAsFactors = FALSE,
  row.names = 1
))

# Number of mouse genes detected in the liver and rhythmics in AdLib
# This is the "Rhythmic" set
liver_Adlib_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(2, 4, 5), ])

# Number of mouse genes detected in the liver and rhythmics in Day Feeding
# This is the "Rhythmic DF" set
liver_DayFed_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(3, 4, 5), ])



# Euler Plots ----
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
  fills = c("lavenderblush", "indianred3", "palegreen", "paleturquoise1")
)
# Mouse = 78298 = 45295+12657+11481+8865  (GRCm39)
# Liver = 33003 = 12657+11481+8865        (>=1 Count in AdLib)
# Detected = 20346 = 11481+8865           (>=10 Counts in AdLib)
# Rhythmic = 8865                         (In dryR models rhythmic in AdLib diet)


plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `Rhythmic AdLib` = liver_Adlib_rhythmic,
      `Rhythmic DayFed` = liver_DayFed_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "lavenderblush",
    "indianred3",
    "palegreen",
    "paleturquoise",
    "salmon1"
  )
)
# Mouse = 78298 = 45295+12657+11481+8865  (GRCm39)
# Liver = 33003 = 12657+11481+8865        (>=1 Count in AdLib)
# Detected = 20346 = 11481+8865           (>=10 Counts in AdLib)
# Rhythmic AdLib = 8865                   (In dryR models rhythmic in AdLib diet)
# Rhythmic DayFed = 2429                  (In dryR models rhythmic in DayFed diet)



# Histogram of biotype ----
FilterGenesTable |>
  ggplot(aes(y = fct_infreq(gene_biotype))) +
  geom_bar(stat = "count") +
  theme_classic() +
  labs(x = "Counts") +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(vjust = 2),
    axis.ticks.y = element_blank(),
    text = element_text(size = 12, color = "black")
  ) +
  scale_x_sqrt(
    breaks = c(0, 1000, 5000, 10000, 15000),
    expand = c(0, 0),
    limits = c(0, 15500)
  )

ggsave(
  "./output/graphs/Genes_biotype_across_samples.png",
  width = 20,
  height = 12,
  units = "cm",
  dpi = 600
)



# Transcript counts by biotype ----
Tr_biotypes <- read.table("./data/multiqc_featurecounts_biotype_plot.txt", header = T)

Tr_biotypes |>
  pivot_longer(!Sample, names_to = "Biotype", values_to = "Counts") |>
  ggplot() +
  geom_boxplot(aes(x = Biotype, y = (Counts + 1)),
               outlier.size = 0.1,
               linewidth = 0.2) +
  geom_hline(yintercept = 1e6, linetype = "dotdash") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 50,
      vjust = 1,
      hjust = 1
    ),
    axis.title.y = element_text(vjust = 2),
    axis.ticks.x = element_blank(),
    text = element_text(size = 11, color = "black")
  ) +
  scale_y_log10(
    breaks = c(10, 1e3, 1e5, 1e6, 1e7),
    minor_breaks = c(
      1:10,
      seq(10, 1000, length.out = 10),
      seq(1000, 1e5, length.out = 10),
      seq(1e5, 1e7, length.out = 10)
    )
  )

ggsave(
  "./output/graphs/Tr_biotypes.png",
  width = 13,
  height = 6,
  units = "cm",
  dpi = 600
)


# RNA Numbers ----

## RNA of Raw count table ----
RawGenesTable <- read.table(file = "./data/RawGenesTable.csv", dec = ".", sep = ";")
colnames(RawGenesTable) <- RawGenesTable[1, ]
RawGenesTable <- RawGenesTable[-1, ]
RawGenesTable <- RawGenesTable[, -1]

lncRNA <- RawGenesTable[RawGenesTable$gene_biotype == "lncRNA", ]
mRNA <- RawGenesTable[RawGenesTable$gene_biotype == "protein_coding", ]

pdf("./output/graphs/EulerPlot_RNA.pdf")
set.seed(502)
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `mRNA` = mRNA$gene_id,
      `Rhythmic` = liver_Adlib_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "lavenderblush",
    "indianred3",
    "palegreen",
    "plum",
    "paleturquoise1"
  )
)

set.seed(20)
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `lncRNA` = lncRNA$gene_id,
      `Rhythmic` = liver_Adlib_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "lavenderblush",
    "indianred3",
    "palegreen",
    "palegoldenrod",
    "paleturquoise1"
  )
)

set.seed(20)
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `mRNA` = mRNA$gene_id,
      `lncRNA` = lncRNA$gene_id,
      `Rhythmic` = liver_Adlib_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "lavenderblush",
    "indianred3",
    "palegreen",
    "plum",
    "palegoldenrod",
    "paleturquoise1"
  )
)

set.seed(5050)
plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected),
      `mRNA` = mRNA$gene_id,
      `lncRNA` = lncRNA$gene_id,
      `Rhythmic AdLib` = liver_Adlib_rhythmic,
      `Rhythmic DayFed` = liver_DayFed_rhythmic
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "lavenderblush",
    "indianred3",
    "palegreen",
    "plum",
    "palegoldenrod",
    "paleturquoise1",
    "salmon1"
  )
)
dev.off()


## RNA by model - Euler ----
model_1 <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% 1, ])
model_2 <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% 2, ])
model_3 <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% 3, ])
model_4 <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% 4, ])
model_5 <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% 5, ])

pdf("./output/graphs/EulerPlot_RNA_per_Model.pdf")

set.seed(1234)
plot(
  eulerr::euler(
    list(
      `mRNA` = mRNA$gene_id,
      `Arrhythmic` = model_1,
      `Loss` = model_2,
      `Gain` = model_3,
      `Identical` = model_4,
      `Reversed` = model_5
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "plum",
    "darkgray",
    "brown3",
    "seagreen3",
    "royalblue",
    "orange"
  )
)

set.seed(66)
plot(
  eulerr::euler(
    list(
      `lncRNA` = lncRNA$gene_id,
      `Arrhythmic` = model_1,
      `Loss` = model_2,
      `Gain` = model_3,
      `Identical` = model_4,
      `Reversed` = model_5
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c(
    "palegoldenrod",
    "darkgray",
    "brown3",
    "seagreen3",
    "royalblue",
    "orange"
  )
)
dev.off()


## RNA by Model - Histogram ----

Biotype.Model <- cbind.data.frame(
  gene_id = rownames(dryList$parameters),
  model = as.factor(dryList$parameters$chosen_model),
  biotype = as.character(RawGenesTable$gene_biotype[RawGenesTable$gene_id %in% rownames(dryList$parameters)])
)

B.M.sorted <- Biotype.Model
B.M.sorted$biotype[B.M.sorted$biotype != "protein_coding" &
                     B.M.sorted$biotype != "lncRNA"] <- "other"

B.M.sorted$biotype <- as.factor(B.M.sorted$biotype)
Biotype.Model$biotype <- as.factor(Biotype.Model$biotype)

Biotype.Model.percent <- Biotype.Model |>
  group_by(model) |>
  dplyr::count(biotype) |>
  mutate(Percent = n / sum(n) * 100) |>
  mutate(Log10 = log10(Percent))

B.M.sorted.percent <- B.M.sorted |>
  group_by(model) |>
  dplyr::count(biotype) |>
  mutate(Percent = n / sum(n) * 100) |>
  mutate(Log10 = log10(Percent))



ggplot(Biotype.Model.percent, aes(x = model, y = Percent, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("All biotypes per model") +
  labs(x = "Model", y = "Percentage", fill = "Biotype")

ggplot(Biotype.Model.percent |>
         filter(Percent > 10),
       aes(x = model, y = Percent, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("10% more abundant biotypes per model") +
  labs(x = "Model", y = "Percentage", fill = "Biotype")

ggplot(Biotype.Model.percent |>
         filter(Percent <= 10),
       aes(x = model, y = Percent, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Less than 10% abundant biotypes per model") +
  labs(x = "Model", y = "Percentage", fill = "Biotype")


# Only mRNA and lncRNA are shown
ggplot(B.M.sorted.percent, aes(x = Percent, y = model, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("All biotypes per model") +
  labs(x = "Percentage", fill = "Biotype") +
  scale_fill_manual(values = c("darkturquoise", "gray", "coral")) +
  scale_y_discrete(limits = rev) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme_classic()


ggplot(B.M.sorted.percent |>
         filter(Percent > 10),
       aes(x = model, y = Percent, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("10% more abundant biotypes per model") +
  labs(x = "Model", y = "Percentage", fill = "Biotype")

ggplot(B.M.sorted.percent |>
         filter(Percent <= 10),
       aes(x = model, y = Percent, fill = biotype)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Les than 10% abundant biotypes per model") +
  labs(x = "Model", y = "Percentage", fill = "Biotype")


ggsave(
  filename = "./output/graphs/Biotype_per_model.png",
  width = 20,
  height = 10,
  dpi = 600,
  units = "cm"
)
