library(GOSemSim)
library(org.Mm.eg.db)
library(ggplot2)
library(tidyverse)
library(kuipeR)

source("./R/fun/Get_go_name.R")
load("./CircanutrEnv.RData")

GOgenes_df <- as.data.frame(org.Mm.egGO2EG)

GOgenes_df <- GOgenes_df[GOgenes_df$Evidence != 'IEA' &
                           GOgenes_df$Ontology == 'BP', ]  # possibility to avoid some specific Evidence code (IEA = Inferred from electronic annotation)

GOgenes_df$ENSEMBL <- mapIds(org.Mm.eg.db,
                             GOgenes_df$gene_id,
                             column = "ENSEMBL",
                             keytype = "ENTREZID")



GOgenes <- split(GOgenes_df$ENSEMBL, GOgenes_df$go_id)



# AdLib
GOphase_AL <- lapply(GOgenes, function (x) {
  Liste <- dryList$results[rownames(dryList$results) %in% x, "phase_CTRL"]
  Liste
}) # 12k

GOphase_AL <- GOphase_AL[sapply(GOphase_AL, function(x) {
  length(x) > 10
})] # 1193
#### élimine les cat du modele arrhythmic, et ceux pas dans liver

GOkuip_AL <- lapply(GOphase_AL, function(x) {
  kuiper_test(x, limits = c(0, 24), method = "simulated")
})

GOkuip_AL <- do.call(rbind.data.frame, GOkuip_AL)
#### 264 significatifs

GOkuip_AL$go_name <- sapply(rownames(GOkuip_AL), Get_go_name)

GO_AdLib <- rownames(GOkuip_AL[GOkuip_AL$p.value < 0.05, ])
