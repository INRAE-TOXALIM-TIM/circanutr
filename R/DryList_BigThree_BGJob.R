# Install and load libraries ----
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("dryR", quietly = TRUE)) {
  devtools::install_github("naef-lab/dryR")
}

library(dryR)
library(tidyverse)

setwd("~/circanutr")

# Prepare arguments ----
# Read and prepare count data
countData.mat <- as.matrix(
  read.csv2(
    "./data/FilterCountTable.csv",
    stringsAsFactors = FALSE,
    row.names = 1
  )
)

# Convert count data to an integer array
countData <- array(as.integer(countData.mat), dim = dim(countData.mat))
colnames(countData) <- colnames(countData.mat)
row.names(countData) <- row.names(countData.mat)

# Read group table
groupTable <- read.csv2("./data/GroupsTable.csv", row.names = 1)
group <- groupTable$Diet
time <- as.numeric(groupTable$ZT)

# DryList 3 diets : AL, NF, DF ####
dryList_BigThree <- dryseq(countData[, grepl("CTRL|RF|^FAR$", group)], group[grepl("CTRL|RF|^FAR$", group)], time[grepl("CTRL|RF|^FAR$", group)])

save.image("./DryList_BigThree_BGEnv.RData")
