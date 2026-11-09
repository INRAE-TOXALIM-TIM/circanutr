# Environment loading ----

packages <- c("BiocManager", "clusterProfiler", "org.Mm.eg.db", "ggplot2")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

if (!requireNamespace("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")
library(enrichplot)

rm(package, packages)

load("./CircanutrEnv.RData")



# GO Gene Set enrichment analysis ----

# Map gene IDs to ENTREZID
Mapped_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = FilterGenesTable$gene_id,
  columns = c('ENTREZID'),
  keytype = 'ENSEMBL'
)



## Ad Lib vs Day Fed ----
pcaALDF_loadings <- read.csv2(file = "./data/RS2205pcaCTRF_loadings.csv")


### PC1 ----

# Create a data frame with gene IDs and PC1 loadings, then merge with mapped genes
ALDF_Mapped_genes.PC1 <- data.frame(ENSEMBL = FilterGenesTable$gene_id, loadings = pcaALDF_loadings$X.PC1) |>
  merge(Mapped_genes, all.x = T)

# Sort gene list by PC1 loadings
ALDF_geneList.PC1 <- sort(
  setNames(object = ALDF_Mapped_genes.PC1$loadings, nm = ALDF_Mapped_genes.PC1$ENTREZID),
  decreasing = T
)

# GSEA for PC1
ALDF_GSEA.PC1 <- gseGO(
  geneList     = ALDF_geneList.PC1,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  seed         = TRUE
) # There are ties in the preranked stats (0.46% of the list.)


#### Enrichment Map
emapplot(pairwise_termsim(ALDF_GSEA.PC1)) # 6 clusters can be seen

#### Reduce redundancy of enriched GO terms
ALDF_GSEA.PC1_simple <- simplify(ALDF_GSEA.PC1,
                                 cutoff = 0.3,
                                 by = "qvalue",
                                 select_fun = min)

# Create a data frame with simplified GSEA results
ALDF_GSEA.PC1_simple_df <- data.frame(
  paths =  paste0(ALDF_GSEA.PC1_simple@result$Description),
  NES =    ALDF_GSEA.PC1_simple@result$NES,
  count =  ALDF_GSEA.PC1_simple@result$setSize,
  qvalue = ALDF_GSEA.PC1_simple@result$qvalue
)


#### Graph of genes per cateGOry
ALDF_GSEA.PC1_Readable <- setReadable(ALDF_GSEA.PC1_simple, OrgDb = "org.Mm.eg.db", keyType = "auto")
clusterProfiler::cnetplot(ALDF_GSEA.PC1_Readable)

#### Plot
ggplot(
  ALDF_GSEA.PC1_simple_df,
  aes(
    x = NES,
    y = forcats::fct_reorder(paths, NES),
    color = -log10(qvalue),
    size = count
  )
) +
  ggtitle("GSEA AL vs DF - PC1") +
  geom_point(stat = 'identity') +
  xlab("NES") +
  ylab("") +
  scale_color_viridis_c(
    option = "inferno",
    begin = 0.5,
    end = 1,
    direction = -1,
    limits = c(0, max(
      -log10(ALDF_GSEA.PC1_simple_df$qvalue), na.rm = TRUE
    ))
  ) +
  scale_size_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  xlim(-2.5, 2.5) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines"))


### PC2 ----

# Create a data frame with gene IDs and PC2 loadings, then merge with mapped genes
ALDF_Mapped_genes.PC2 <- data.frame(ENSEMBL = FilterGenesTable$gene_id, loadings = pcaALDF_loadings$X.PC2) |>
  merge(Mapped_genes, all.x = T)

# Sort gene list by PC2 loadings
ALDF_geneList.PC2 <- sort(
  setNames(object = ALDF_Mapped_genes.PC2$loadings, nm = ALDF_Mapped_genes.PC2$ENTREZID),
  decreasing = T
)

# GSEA for PC2
ALDF_GSEA.PC2 <- gseGO(
  geneList     = ALDF_geneList.PC2,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  seed         = TRUE
) # There are ties in the preranked stats (0.46% of the list.)


#### Enrichment Map
emapplot(pairwise_termsim(ALDF_GSEA.PC2)) # 6 clusters can be seen

#### Reduce redundancy of enriched GO terms
ALDF_GSEA.PC2_simple <- simplify(ALDF_GSEA.PC2,
                                 cutoff = 0.3,
                                 by = "qvalue",
                                 select_fun = min)

# Create a data frame with simplified GSEA results
ALDF_GSEA.PC2_simple_df <- data.frame(
  paths =  paste0(ALDF_GSEA.PC2_simple@result$Description),
  NES =    ALDF_GSEA.PC2_simple@result$NES,
  count =  ALDF_GSEA.PC2_simple@result$setSize,
  qvalue = ALDF_GSEA.PC2_simple@result$qvalue
)


#### Graph of genes per cateGOry
ALDF_GSEA.PC2_Readable <- setReadable(ALDF_GSEA.PC2_simple, OrgDb = "org.Mm.eg.db", keyType = "auto")
clusterProfiler::cnetplot(ALDF_GSEA.PC2_Readable)


#### Plot
ggplot(
  ALDF_GSEA.PC2_simple_df,
  aes(
    x = NES,
    y = forcats::fct_reorder(paths, NES),
    color = -log10(qvalue),
    size = count
  )
) +
  ggtitle("GSEA AL vs DF - PC2") +
  geom_point(stat = 'identity') +
  xlab("NES") + ylab("") +
  scale_color_viridis_c(
    option = "inferno",
    begin = 0.5,
    end = 1,
    direction = -1,
    limits = c(0, max(
      -log10(ALDF_GSEA.PC2_simple_df$qvalue), na.rm = TRUE
    ))
  ) +
  scale_size_continuous(limits = c(0, NA)) +
  scale_size_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  xlim(-2.5, 2.5) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines"))



## Ad Lib vs Fasting ----
pcaALFA_loadings <- read.csv2(file = "./data/RS2205pcaCTF_loadings.csv")

### PC1 ----

# Create a data frame with gene IDs and PC1 loadings, then merge with mapped genes
ALFA_Mapped_genes.PC1 <- data.frame(ENSEMBL = FilterGenesTable$gene_id, loadings = pcaALFA_loadings$X.PC1) |>
  merge(Mapped_genes, all.x = T)

# Sort gene list by PC1 loadings
ALFA_geneList.PC1 <- sort(
  setNames(object = ALFA_Mapped_genes.PC1$loadings, nm = ALFA_Mapped_genes.PC1$ENTREZID),
  decreasing = T
)

# GSEA for PC1
ALFA_GSEA.PC1 <- gseGO(
  geneList     = ALFA_geneList.PC1,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  seed         = TRUE
) # There are ties in the preranked stats (0.46% of the list).


#### Enrichment Map
emapplot(pairwise_termsim(ALFA_GSEA.PC1)) # 9 clusters can be seen

#### Reduce redundancy of enriched GO terms
ALFA_GSEA.PC1_simple <- simplify(ALFA_GSEA.PC1,
                                 cutoff = 0.3,
                                 by = "qvalue",
                                 select_fun = min)

# Create a data frame with simplified GSEA results
ALFA_GSEA.PC1_simple_df <- data.frame(
  paths =  paste0(ALFA_GSEA.PC1_simple@result$Description),
  NES =    ALFA_GSEA.PC1_simple@result$NES,
  count =  ALFA_GSEA.PC1_simple@result$setSize,
  qvalue = ALFA_GSEA.PC1_simple@result$qvalue
)


#### Graph of genes per cateGOry
ALFA_GSEA.PC1_Readable <- setReadable(ALFA_GSEA.PC1_simple, OrgDb = "org.Mm.eg.db", keyType = "auto")
clusterProfiler::cnetplot(ALFA_GSEA.PC1_Readable)


#### Plot
ggplot(
  ALFA_GSEA.PC1_simple_df,
  aes(
    x = NES,
    y =  forcats::fct_reorder(paths, NES),
    color = -log10(qvalue),
    size = count
  )
) + ggtitle("GSEA AL vs FA - PC1") +
  geom_point(stat = 'identity') +
  xlab("NES") + ylab("") +
  scale_color_continuous(limits = c(0, NA)) +
  scale_size_continuous(limits = c(0, NA)) +
  theme_bw()


## PC2 ----

# Create a data frame with gene IDs and PC2 loadings, then merge with mapped genes
ALFA_Mapped_genes.PC2 <- data.frame(ENSEMBL = FilterGenesTable$gene_id, loadings = pcaALFA_loadings$X.PC2) |>
  merge(Mapped_genes, all.x = T)

# Sort gene list by PC2 loadings
ALFA_geneList.PC2 <- sort(
  setNames(object = ALFA_Mapped_genes.PC2$loadings, nm = ALFA_Mapped_genes.PC2$ENTREZID),
  decreasing = T
)

# GSEA for PC2
ALFA_GSEA.PC2 <- gseGO(
  geneList     = ALFA_geneList.PC2,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  minGSSize    = 100,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  seed         = TRUE
) # There are ties in the preranked stats (0.46% of the list.)


#### Enrichment Map
emapplot(pairwise_termsim(ALFA_GSEA.PC2)) # 9 clusters can be seen

#### Reduce redundancy of enriched GO terms
ALFA_GSEA.PC2_simple <- simplify(ALFA_GSEA.PC2,
                                 cutoff = 0.3,
                                 by = "qvalue",
                                 select_fun = min)

# Create a data frame with simplified GSEA results
ALFA_GSEA.PC2_simple_df <- data.frame(
  paths =  paste0(ALFA_GSEA.PC2_simple@result$Description),
  NES =    ALFA_GSEA.PC2_simple@result$NES,
  count =  ALFA_GSEA.PC2_simple@result$setSize,
  qvalue = ALFA_GSEA.PC2_simple@result$qvalue
)


#### Graph of genes per cateGOry
ALFA_GSEA.PC2_Readable <- setReadable(ALFA_GSEA.PC2_simple, OrgDb = "org.Mm.eg.db", keyType = "auto")
clusterProfiler::cnetplot(ALFA_GSEA.PC2_Readable)

#### Plot
ggplot(
  ALFA_GSEA.PC2_simple_df,
  aes(
    x = NES,
    y =  forcats::fct_reorder(paths, NES),
    color = -log10(qvalue),
    size = count
  )
) + ggtitle("GSEA AL vs FA - PC2") +
  geom_point(stat = 'identity') +
  xlab("NES") + ylab("") +
  scale_color_continuous(limits = c(0, NA)) +
  scale_size_continuous(limits = c(0, NA)) +
  theme_bw()
