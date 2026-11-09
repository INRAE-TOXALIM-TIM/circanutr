# Set up #######################################################################
devtools::install_gitlab("arnaud.polizzi/kuipeR", host = "https://forge.inrae.fr")
library(kuipeR)

library(tidyverse)
library(GOSemSim)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)

source("./R/fun/Get_go_name.R")

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

# Not saved, made in post with Illustrator / Affinity


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
  height = 9,
  dpi = 600,
  units = "cm"
)



# C. Go categories ranking #####################################################
GOgenes_df <- as.data.frame(org.Mm.egGO2EG)

GOgenes_df <- GOgenes_df[GOgenes_df$Evidence != 'IEA' &
                           GOgenes_df$Ontology == 'BP', ]
# possibility to avoid some specific Evidence code (IEA = Inferred from electronic annotation)

GOgenes_df$ENSEMBL <- mapIds(org.Mm.eg.db,
                             GOgenes_df$gene_id,
                             column = "ENSEMBL",
                             keytype = "ENTREZID")


GOgenes <- split(GOgenes_df$ENSEMBL, GOgenes_df$go_id)

GO_rhythmic_AL <- data.frame(dryList$parameters["chosen_model"]) |>
  filter(chosen_model %in% c(3, 4, 5)) |>
  rownames_to_column(var = "ENSEMBL") |>
  inner_join(y = GOgenes_df, by = "ENSEMBL") |>
  filter(go_id != "GO:0008150") |>
  group_by(go_id)  |>
  summarise(count = n_distinct(ENSEMBL)) |>
  mutate(go_name = sapply(go_id, Get_go_name))
  # Set of rhythmic genes in AdLib

####
####

BiocManager::install("rrvgo")
library(rrvgo)

Rhythmic_Matrix <- calculateSimMatrix(
  GO_rhythmic_AL$go_id,
  orgdb = "org.Mm.eg.db",
  ont = "BP",
  method = "Rel"
)

reduced.uniq.AL <- reduceSimMatrix(Rhythmic_Matrix,
                                   "uniqueness",
                                   threshold = 0.7,
                                   orgdb = "org.Mm.eg.db")

reduced.size.AL <- reduceSimMatrix(Rhythmic_Matrix,
                                   "size",
                                   threshold = 0.7,
                                   orgdb = "org.Mm.eg.db")

GO_rhythmic_AL$size <- reduced.size.AL$size[match(GO_rhythmic_AL$go_id, reduced.size.AL$go)]
GO_rhythmic_AL$parent <- reduced.size.AL$parent[match(GO_rhythmic_AL$go_id, reduced.size.AL$go)]

####
####

GO_rhythmic_AL |>
  arrange(desc(count)) |>
  head(10) |>
  ggplot(aes(x = reorder(go_name, count), y = count)) +
  geom_col(fill = "#67987E") +
  scale_y_continuous(expand = 0, breaks = seq(0, 450, by = 100)) +
  coord_flip() +
  ylab(label = "Gene count") +
  theme_classic() +
  theme(
    # aspect.ratio = 1 / 2,
    text = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 12, margin = margin(r = 5)),
    plot.margin = margin(l = 5, r = 5),
    panel.spacing = unit(1, "lines")
  )

ggsave(
  filename = "C_Functionnal_Categories.svg",
  path = "~/circanutr/output/",
  width = 9,
  height = 15,
  dpi = 600,
  units = "cm"
)
