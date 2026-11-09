# Install and load libraries ----
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("dryR", quietly = TRUE)) {
  devtools::install_github("naef-lab/dryR")
}

library(dryR)
library(tidyverse)

# Prepare arguments ----
# Read and prepare count data
countData.mat <- as.matrix(read.csv2(
  "./data/FilterCountTable.csv",
  stringsAsFactors = FALSE,
  row.names = 1
))

# Convert count data to an integer array
countData <- array(as.integer(countData.mat), dim = dim(countData.mat))
colnames(countData) <- colnames(countData.mat)
row.names(countData) <- row.names(countData.mat)

# Read group table
groupTable <- read.csv2("./data/GroupsTable.csv", row.names = 1)
group <- groupTable$Diet
time <- as.numeric(groupTable$ZT)

# AdLib (CTRL) and Day Fed (RF) ----
# CTRL group is AdLib
# Restricted Feeding (RF) group is Day Fed

# Run the analysis for count data (e.g., RNA-Seq data)
dryList <- dryseq(countData[, grepl("CTRL|RF", group)], group[grepl("CTRL|RF", group)], time[grepl("CTRL|RF", group)])

# Number of models
max(dryList$results$chosen_model)
# 5 expression models

# Number of genes per model
dryList$results |>
  group_by(chosen_model) |>
  summarise(model_size = n()) # Count genes by model
# chosen_model model_size
#            1       9052
#            2       3496
#            3       2429
#            4       1265
#            5       4104

# dryList components:
# dryList[["results"]]      # data frame summarizing results
# dryList[["parameters"]]   # coefficients: phase, amplitude and mean for each group
# dryList[["ncounts"]]      # normalized counts
# dryList[["counts"]]       # raw counts
# dryList[["cook"]]         # cook's distance for outlier detection
# dryList[["BICW_rhythm"]]  # BICW for each rhythmic model
# dryList[["BICW_mean"]]    # BICW for each mean model

# Saving Analysis ----

# Date.Time format
# SessionName <- paste("Circanutr",
# _,
# format(Sys.time(), "%Y%m%d_%H%M%S"),
# ".RData",
# sep = ""
# )

# save.image(SessionName)

# load("CircanutrEnv.RData")

# Writing tables of dryList object in a directory
sapply(names(dryList), function(x) {
  write.csv2(dryList[[x]], file = paste0("./output/dryLists/dryList_", x, ".csv"))
})

# Generate heatmaps of each model as PDF
plot_models_rhythm(dryList, "./output/graphs/")

# Additional results ----
# Get gene_id of a gene, Bmal1 example
FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Bmal1"]

# Plot gene expression, Bmal1 example
dry_plot(dryList, "ENSMUSG00000055116")

# Loading the result table
dryList <- list()
dryList$results <- read.csv2(file = "./output/dryLists/dryList_results.csv")

# Non circadian diets ----
# Source the Compare_lm function
source("./R/fun/Compare_lm.R")

## AdLib (CTRL) and Fasting (F) ----
dryList.Fasting <- dryseq(countData[, grepl("CTRL|^F$", group)], group[grepl("CTRL|^F$", group)], time[grepl("CTRL|^F$", group)])
# F group doesn't have ZT0 in this run.

Fasting.lm <- Compare_lm(dryList.Fasting, groupTable, "F")

# save.image("FastingLm.RData")


## AdLib (CTRL) and Refed (FAR) ----
dryList.Refed <- dryseq(countData[, grepl("CTRL|^FAR$", group)], group[grepl("CTRL|^FAR$", group)], time[grepl("CTRL|^FAR$", group)])
# FAR group doesn't have ZT0 in this run.

Refed.lm <- Compare_lm(dryList.Refed, groupTable, "FAR")

# save.image("All_Diets.RData")



# DryList 3 diets : AL, NF, DF ####

#### Pre-run dryList
load("~/circanutr/DryList_BigThree_BGEnv.RData")
# DryList_BigThree <- dryseq(countData[, grepl("CTRL|RF|^FAR$", group)], group[grepl("CTRL|RF|^FAR$", group)], time[grepl("CTRL|RF|^FAR$", group)])


#### Number of models
max(dryList_BigThree$results$chosen_model)
# 15 models

#### Number of genes per model
dryList_BigThree$results |>
  group_by(chosen_model) |>
  summarise(model_size = n()) # Count genes by model
#   chosen_model  model_size
# 1            1        6875
# 2            2         776
# 3            3        1275
# 4            4         730
# 5            5        3653  => Flat in DayFed (Per1-like)
# 6            6         343  => Rora-like
# 7            7         484
# 8            8         171
# 9            9         136
# 10           10        972
# 11           11        501
# 12           12         44
# 13           13        283
# 14           14       3386  => Reversed in DayFed (Bmal1-like)
# 15           15        717  => Shift in NF and DF (Per2-like)


FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Bmal1"]
# ENSMUSG00000055116
dryList_BigThree$parameters[, dryList_BigThree$parameters$chosen_model == 14]
# 14

#### Heatmaps of models
dryR::plot_models_rhythm(dryList_BigThree, "~/circanutr/output/BigThree/")



# dryList components:
# dryList[["results"]]      # data frame summarizing results
# dryList[["parameters"]]   # coefficients: phase, amplitude and mean for each group
# dryList[["ncounts"]]      # normalized counts
# dryList[["counts"]]       # raw counts
# dryList[["cook"]]         # cook's distance for outlier detection
# dryList[["BICW_rhythm"]]  # BICW for each rhythmic model
# dryList[["BICW_mean"]]    # BICW for each mean model

# Saving Analysis ----

# Date.Time format
# SessionName <- paste("Circanutr",
# _,
# format(Sys.time(), "%Y%m%d_%H%M%S"),
# ".RData",
# sep = ""
# )

# save.image(SessionName)

# load("CircanutrEnv.RData")

# Writing tables of dryList object in a directory
sapply(names(dryList_BigThree), function(x) {
  write.csv2(dryList_BigThree[[x]],
             file = paste0("./output/dryLists/dryList_BigThree_", x, ".csv"))
})

# Generate heatmaps of each model as PDF
plot_models_rhythm(dryList_BigThree, "./output/graphs/dryList_BigThree/")

# Additional results ----
# Get gene_id of a gene, Bmal1 example
FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Bmal1"]

# Plot gene expression, Bmal1 example
dry_plot(dryList_BigThree, "ENSMUSG00000055116")

# Loading the result table
dryList_BigThree <- list()
dryList_BigThree$results <- read.csv2(file = "./output/dryLists/dryList_results.csv")



FilterGenesTable$gene_id[FilterGenesTable$gene_name == "Rora"]
# ENSMUSG00000055116 : Bmal1              14
# ENSMUSG00000029238 : Clock              14
# ENSMUSG00000020038 : Cry1               14
# ENSMUSG00000068742 : Cry2               14
# ENSMUSG00000020889 : Nr1d1, Rev-erba    14
# ENSMUSG00000021775 : Nr1d2, Rev-erbb    14
# ENSMUSG00000020893 : Per1               5
# ENSMUSG00000055866 : Per2               15
# ENSMUSG00000032238 : Rora               6
dryList_BigThree$parameters["ENSMUSG00000032238", "chosen_model"]
# 14, 5, 15, 6
