#################################
# Config session
#################################

## On genobioinfo
cd /work/project/TIM/bioinfo/RS2205/biostat
module load statistics/R/4.3.0
srun --x11 --pty --mem=24G R

## On R
### Project variables
ProjectName <- "RS2205"
path <- "/work/project/TIM/bioinfo"
proj_path <- file.path(path, ProjectName)
wd_path <- file.path(proj_path, "biostat")
dir.create(wd_path)
setwd(wd_path)
lit_dirs2c <- file.path(
  wd_path,
  c(
    "QC/CAH",
    "QC/Distributions",
    "QC/MDS",
    "QC/PCA",
    "QC/pairsPlots",
    "DEG",
    "GO"
  )
)
lapply(lit_dirs2c, function(x)
  dir.create(x, recursive = T, showWarnings = F))


### Saving session
SessionName <- paste(ProjectName,
                    "_",
                    format(Sys.time(), "%Y%m%d_%H%M%S"),
                    ".RData",
                    sep = ""
)
save.image(SessionName)


### Loading session
load("RS2205_20250123_120703.RData")


### Libraries and functions
library(limma)
library(edgeR)
library(data.table)
library(dplyr)
library(ggrepel)
library(colorspace)
library(mixOmics)

source(file.path(fun_path, "QCplots_edgeR.R"))
source(file.path(fun_path, "extract_stat_results_edger.R"))



#################################
# Loading data
#################################

## Salmon counts
	#### sample_names must start with RS***
countFn <- "../nf_results/star_salmon/salmon.merged.gene_counts.tsv"

counts <- fread(countFn, stringsAsFactors = FALSE)	# load counts table

genes <- counts[, c(1:2)]	# gene_ID and gene_name

counts <- counts[, -c(1:2)]	# counts only (without gene_ID nor gene_name)

targets <- fread("RS2205_targets.txt", dec = ",")	# sample table and caracteristics of extraction

targets2 <- as_tibble(colnames(counts)) %>% dplyr::left_join(targets, c("value" =
                                                                          "SampleName")) %>% dplyr::rename(sample_name = value
)

targets2$group <- factor(targets2$group)
groups <- factor(targets2$group)	# group = Diet_ZT

conditions <- as.data.frame(strsplit2(groups, "_"))

colnames(conditions) <- c("Diet", "ZT")

targets2 <- cbind(targets2, conditions)


## Annotations
### GTF
gtf <- rtracklayer::import(
  "/work/project/TIM/bioinfo/RS2205/nf_results/genome/ensembl_mus_musculus_genome_GRCm39_2024-8-17.filtered.gtf"
)

gtf_df <- as.data.frame(gtf)

gtf_dfgene = gtf_df[gtf_df$type == "transcript", ]

genes <- dplyr::left_join(
  genes[, 1:2],
  dplyr::select(
    gtf_dfgene,
    gene_id,
    gene_biotype,
    seqnames,
    start,
    end,
    width,
    strand
  ),
  multiple = "first"
) %>% dplyr::rename(chromosome = seqnames)	# adding cols to genes

	#### Always have id ensembl in 'genes' objects with colnames as ENSEMBL
genes$ENSEMBL <- genes$gene_id


### Free disk spae
lss()
rm(gtf, gtf_dfgene, gtf_df)
lss()
gc()



#################################
# Dge object
#################################

dge <- DGEList(counts = counts,
               genes = genes,
               samples = targets2
)

dim(dge)
	#### [1] 78298   120
	####	 genes	 samples



#################################
# Plots aesthetics
#################################

## RNA libraries
lib_peack_size <-  findInterval(dge$samples$lib_peack_size,
                                sort(dge$samples$lib_peack_size))
lib_mean_size <-  findInterval(dge$samples$lib_mean_size, sort(dge$samples$lib_mean_size))
lib_molarity <-  findInterval(dge$samples$lib_molarity, sort(dge$samples$lib_molarity))
lib_target_prop <-  findInterval(dge$samples$lib_target_prop,
                                 sort(dge$samples$lib_target_prop))

pchlib_serial <- as.numeric(factor(dge$samples$lib_serial)) - 1
pchlib_serial[pchlib_serial == 0] <- 16

## Conditions
colorsDiet <- c("#007FFF", "#E0115F", "#67987E", "#FF7F00")	# colors (#50C878)
pchDiet <- c(15, 17, 18, 16)
pchlib_Diet <- as.numeric(factor(dge$samples$Diet)) - 1
pchlib_Diet[pchlib_Diet == 0] <- 4

pchZT <- c(1, 0, 2, 3, 15, 17, 13)
pchlib_ZT <- factor(c("00", "04", "08", "12", "16", "20", "24"))

	#### Add 00, 04, 08 in group names
groups_00 <- paste0(dge$samples$Diet, "_ZT", dge$samples$ZT) %>% factor
	
	#### Color fading for WS (ordered as CTRL, F, FAR, RF)
gp_colors <- c("#000000", "#012D76", "#0051C6", "#017EFF", "#2EB0FE", "#78C4FE", "#CCE8FE",
			   "#501028", "#85123D", "#B11650", "#E84481", "#EE7BA6", "#F4B3CB", "#FBE1EC",
			   "#10501E", "#12852A", "#16B237", "#44E768", "#96F1AB", "#BCF5C8", "#EDFDF2",
			   "#64492C", "#A7590E", "#DF6F01", "#FE8103", "#FEA346", "#FFD2A8", "#FFF4E8"
			   )
gp_colors[as.numeric(groups_00)]


#################################
# QC raw
#################################

QCplots_edgeR(
  dge,
  wd_path = wd_path,
  normalized = F,
  prefix = "All",
  suffix = "RawSet",
  mds = F,
  pca = T,
  fun_path = fun_path
)

	#### Color pallet
hcl_palettes(plot = TRUE)

	#### Quick view
plot(
  1:10,
  pch = 16,
  cex = 5,
  col =  sequential_hcl(10, palette = c(
  "#000000","#446F90","#609DCC","#75BEF7","#8FCEFE","#AADAFE")
))

dev.off()

dge$samples$ZT <- formatC(sub("ZT", "", dge$samples$ZT) %>% as.numeric,
                          width = 2,
                          flag = 0
)


## Plot
pcresR <- pca(t(cpm(
  dge,
  log = T,
  prior = 1,
  normalized.lib.sizes = F
)), ncomp = 6)

pdf("QC/PCA/All_plotPCA_log2CPM_RawSet_explore.pdf")
plot(pcresR)


## PCA: 1v2, 3v4, 5v6	--	Raw Set
for (i in seq(1, 5, 2)) {
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    group = as.character(dge$samples$Diet),
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 3,
    pch = 16
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    group = as.character(dge$samples$Diet),
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 1
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    group = as.character(dge$samples$Diet),
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = pchZT,
    pch.levels = pchlib_ZT,
    col.per.group = colorsDiet,
    legend.title = "Diet",
    legend.title.pch = "ZT",
    cex = 5
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_peack_size))[lib_peack_size],
    cex = 3,
    pch = 16
  )
  legend(
    "topleft",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dge$samples$lib_peack_size
    ), 1)),
    title = "lib_peack_size"
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_mean_size))[lib_mean_size],
    cex = 3,
    pch = 16
  )
  legend(
    "topleft",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dge$samples$lib_mean_size
    ), 1)),
    title = "lib_mean_size"
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_molarity))[lib_molarity],
    cex = 3,
    pch = 16
  )
  legend(
    "topleft",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dge$samples$lib_molarity
    ), 1)),
    title = "lib_molarity"
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_target_prop))[lib_target_prop],
    cex = 3,
    pch = 16
  )
  legend(
    "topleft",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dge$samples$lib_target_prop
    ), 1)),
    title = "lib_target_prop"
  )
  
  plotIndiv(
    pcresR,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_target_prop))[lib_target_prop],
    cex = 3,
    pch = 16
  )
  legend(
    "topleft",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dge$samples$lib_target_prop
    ), 1)),
    title = "lib_target_prop"
  )
  
}

dev.off()



#################################
# Dge filtering and pretreatment (Selected)
#################################

## Filter by expr for All genes
keep <- filterByExpr(dge, group=dge$samples$group, min.count=10)
table(keep)
	#### FALSE  TRUE
	#### 57952 20346

dgef <- dge[keep, , keep.lib.sizes=FALSE]
dim(dgef)
	#### 20346   120
	#### genes	 samples

save.image(SessionName)

sort(table(dgef$genes$gene_biotype), dec=T)
	# # # #                     protein_coding                             lncRNA
	# # # #                              14162                               3970
	# # # #                                TEC               processed_pseudogene
	# # # #                               1241                                545
	# # # #   transcribed_processed_pseudogene                             snoRNA
	# # # #                                117                                102
	# # # # transcribed_unprocessed_pseudogene             unprocessed_pseudogene
	# # # #                                 52                                 46
	# # # #                              snRNA                           misc_RNA
	# # # #                                 38                                 28
	# # # #                             scaRNA                          IG_C_gene
	# # # #                                 13                                  9
	# # # #     transcribed_unitary_pseudogene                               rRNA
	# # # #                                  8                                  3
	# # # #                          TR_C_gene                            Mt_rRNA
	# # # #                                  3                                  2
	# # # #                           ribozyme                 unitary_pseudogene
	# # # #                                  2                                  2
	# # # #                              miRNA                            Mt_tRNA
	# # # #                                  1                                  1
	# # # #  translated_unprocessed_pseudogene
	# # # #                                  1


### Generating plots
QCplots_edgeR(
  dgef,
  wd_path = wd_path,
  normalized = F,
  prefix = "All",
  suffix = "SelSet",
  boxplot = T,
  densities = T,
  mds = F,
  hc = T,
  pca = T,
  pairsplot = T,
  fun_path = fun_path
)


### Plot
pcresS = pca(t(cpm(
  dgef,
  log = T,
  prior = 1,
  normalized.lib.sizes = F
)), ncomp = 6)

pdf("QC/PCA/All_plotPCA_log2CPM_SelSet.pdf")


## PCA: 1v2, 3v4, 5v6	--	Selected Set
for (i in seq(1, 5, 2)) {
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    group =  as.character(dgef$samples$Diet),
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 3,
    pch = 16
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    group = as.character(dgef$samples$Diet),
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 1
  )
  
    plotIndiv(
    pcresS,
    comp = 1:2,
    ind.names = T,
    group = as.character(dgef$samples$Diet),
	style = "ggplot2",
    legend = T,
    scale = F,
    pch = pchZT,
    pch.levels = pchlib_ZT,
    col.per.group = colorsDiet,
    legend.title = "Diet",
    legend.title.pch = "ZT",
    cex = 5
  )

  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    group = dgef$samples$Diet,
    ellipse = T,
	style = "ggplot2",
	title = "First two principal components PCA",
    legend = T,
    scale = F,
    pch = pchZT,
    pch.levels = pchlib_ZT,
    col.per.group = colorsDiet,
    legend.title = "Diet",
    legend.title.pch = "ZT",
    cex = 5
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    group = dgef$samples$Diet,
    ellipse = F,
	style = "ggplot2",
    legend = T,
    scale = F,
    pch = pchZT,
    pch.levels = pchlib_ZT,
    col.per.group = colorsDiet,
    legend.title = "Diet",
    legend.title.pch = "ZT",
    cex = 5
  )

  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_peack_size))[lib_peack_size],
    cex = 3,
    pch = 16
  )
  legend(
    "topright",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dgef$samples$lib_peack_size
    ), 1)),
    title = "lib_peack_size"
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_mean_size))[lib_mean_size],
    cex = 3,
    pch = 16
  )
  legend(
    "topright",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dgef$samples$lib_mean_size
    ), 1)),
    title = "lib_mean_size"
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_molarity))[lib_molarity],
    cex = 3,
    pch = 16
  )
  legend(
    "topright",
    col = pal(2),
    pch = 19,
    legend = c(round(range(
      dgef$samples$lib_molarity
    ), 1)),
    title = "lib_molarity"
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_target_prop))[lib_target_prop],
    cex = 3,
    pch = 16
  )
  legend(
    "topright",
    col = pal(2),
    pch = 19,
    legend = c(round(
      range(dgef$samples$lib_target_prop), 1
    )),
    title = "lib_target_prop"
  )
  
  plotIndiv(
    pcresS,
    comp = i:(i + 1),
    ind.names = T,
    style = "graphics",
    legend = F,
    scale = F,
    col = pal(length(lib_target_prop))[lib_target_prop],
    cex = 3,
    pch = 16
  )
  legend(
    "topright",
    col = pal(2),
    pch = 19,
    legend = c(round(
      range(dgef$samples$lib_target_prop), 1
    )),
    title = "lib_target_prop"
  )
  
}

dev.off()



#################################
# Normalization
#################################

groups <- dgef$samples$group

design <- model.matrix(~0 + groups)
colnames(design) <- sub("groups","",colnames(design))

dgef$design <- design


## normFactors
dgef <- calcNormFactors(dgef, method = "TMM")
cat("\nNormalization factors:\n")
print(dgef$samples$norm.factors)

mar0 = par("mar")

png("QC/Distributions/All_norm_factors.png", width = 900)

par(mar = mar0 + c(5, 0, 5, 0))
barplot(
  dgef$samples$norm.factors,
  col = as.numeric(groups) + 1,
  names.arg = colnames(dgef),
  las = 2
)
legend(
  "topright",
  horiz = F,
  levels(groups),
  text.col = 2:(length(levels(groups)) + 1),
  cex = 0.8
)

dev.off()
save.image(SessionName)


QCplots_edgeR(
  dgef,
  wd_path = wd_path,
  normalized = T,
  prefix = "All",
  suffix = "WorkingSet",
  boxplot = T,
  densities = T,
  mds = T,
  hc = T,
  pca = T,
  pairsplot = T,
  fun_path = fun_path
)


## Plot
pcres = pca(t(cpm(
  dgef,
  log = T,
  prior = 1,
  normalized = T
)), ncomp = 6)

pdf("QC/PCA/All_plotPCA_log2CPM_WorkingSet.pdf")

plot(pcres)


## PCA: 1v2, 3v4, 5v6	--	Working Set
for (i in seq(1, 5, 2)) {
  plotIndiv(
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = dgef$samples$Diet,
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 3,
    pch = 16
  )
  
  plotIndiv(
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = dgef$samples$Diet,
    style = "graphics",
    legend = T,
    scale = F,
    col.per.group = colorsDiet,
    cex = 1
  )
  
  plotIndiv(
  # symbols per time
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = dgef$samples$Diet,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = pchZT,
    pch.levels = pchlib_ZT,
    col.per.group = colorsDiet,
    legend.title = "Diet",
    legend.title.pch = "ZT",
    cex = 5
  )
  
  plotIndiv(
    # color shading
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = gp_colors,
    legend.title = "Groups",
    legend.title.pch = "Diet",
    cex = 5
  )
  
  plotIndiv(
    # CT colors
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = c(gp_colors[1:6], paste0("gray", seq(90, 73, -1))) ,
    legend.title = "Groups",
    cex = 5
  )
  
  plotIndiv(
    # F colors
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group =  c(paste0("gray", seq(90, 85, -1)), gp_colors[7:12], paste0("gray", seq(84, 73, -1))),
    legend.title = "Groups",
    cex = 5
  )

  plotIndiv(
    # FAR colors
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = c(paste0("gray",seq(90,79,-1)), gp_colors[13:18], paste0("gray",seq(78,73,-1))) ,
    legend.title = "Groups",
    cex = 5
  )

  plotIndiv(
    # RF colors
    pcres,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = c(paste0("gray",seq(90,73,-1)), gp_colors[19:24]) ,
    legend.title = "Groups",
    cex = 5
  )
}

dev.off()



#################################
# Refining comparaisons
#################################

## CTRL vs RF
CT_RF <- dgef$samples[grep("CTRL|RF", dgef$samples$Diet),]

groups_00_CTRF <- paste0(dgef$samples[grep("CTRL|RF", dgef$samples$Diet),11], "_ZT", dgef$samples[grep("CTRL|RF", dgef$samples$Diet),12]) %>% factor

CTRF_colors <- c(
  "#000000","#012D76","#0051C6","#017EFF","#2EB0FE","#78C4FE",
  "#64492C","#A7590E","#DF6F01","#FE8103","#FEA346","#FFD2A8"
  )
  
pcaCTRF <- pca(t(cpm(
    dgef,
    log = T,
    prior = 1,
    normalized = T
  )[, grep("CTRL|RF", dgef$samples$Diet)]), ncomp = 6)

pdf("QC/PCA/All_plotPCA_log2CPM_WS_CTRLvsRF.pdf")
plot(pcaCTRF)
for (i in seq(1, 5, 2)) {
  plotIndiv(
    # color shading
    pcaCTRF,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00_CTRF,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = CTRF_colors,
    legend.title = "Groups",
    legend.title.pch = "Diet",
    cex = 5
  )
  
    plotIndiv(
    # color shading - no legend
    pcaCTRF,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00_CTRF,
    style = "ggplot2",
    legend = F,
    scale = F,
    pch = 16,
    col.per.group = CTRF_colors,
    legend.title = "Groups",
    legend.title.pch = "Diet",
    cex = 5
  )

}
dev.off()

## F vs FAR
F_FAR <- dgef$samples[grep("^F|FAR", dgef$samples$Diet),]

groups_00_FFAR <- paste0(dgef$samples[grep("^F|FAR", dgef$samples$Diet),11], "_ZT", dgef$samples[grep("^F|FAR", dgef$samples$Diet),12]) %>% factor

FFAR_colors <- c(
  "#85123D","#B11650","#E84481","#EE7BA6","#F4B3CB","#FBE1EC",
  "#12852A","#16B237","#44E768","#96F1AB","#BCF5C8","#EDFDF2"
  )

pcaFFAR <- pca(t(cpm(
    dgef,
    log = T,
    prior = 1,
    normalized = T
  )[, grep("^F|FAR", dgef$samples$Diet)]), ncomp = 6)

pdf("QC/PCA/All_plotPCA_log2CPM_WS_FvsFAR.pdf")
plot(pcaFFAR)
for (i in seq(1, 5, 2)) {  
  plotIndiv(
    # Legend
    pcaFFAR,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00_FFAR,
    style = "ggplot2",
    legend = T,
    scale = F,
    pch = 16,
    col.per.group = FFAR_colors,
    legend.title = "Groups",
    legend.title.pch = "Diet",
    cex = 5
  )

  plotIndiv(
    # No legend
    pcaFFAR,
    comp = i:(i + 1),
    ind.names = T,
    group = groups_00_FFAR,
    style = "ggplot2",
    legend = F,
    scale = F,
    pch = 16,
    col.per.group = FFAR_colors,
    legend.title = "Groups",
    legend.title.pch = "Diet",
    cex = 5
  )
  
}
dev.off()



#################################
# Exporting file
#################################

save.image(SessionName)


export_results_matrix <- function(dge, prefix = "All") {
  #### pData
  write.csv2(dge$sample[, 1, drop = F],
             paste0(ProjectName, "_", prefix, "_pData.csv"),
             row.names = T)
  #### WorkingSet
  log2CPM <- cpm(
    dge,
    normalized.lib.sizes = T,
    log = T,
    prior.count = 1
  )
  write.csv2(
    cbind.data.frame(
      Geneid = dge$genes[, 1],
      log2CPM,
      stringsAsFactors = F
    ),
    paste0(ProjectName, "_", prefix, "_WorkingSet_log2CPM.csv"),
    row.names = F
  )
}

export_results_matrix(dgef, prefix = "All")
