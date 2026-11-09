#' Rhythmic_expression
#'
#' @description Creates graphs of the variation of expressions for selected genes through time.
#' 
#' @param dryList A DryR list containing dataframes with results, rhythmic parameters, normalized counts, raw counts, and the Cook's distance produced by the dryseq() function.
#' @param FilterGenesTable The database for gene names. It must contain a column "gene_id" with ENSEMBL codename for genes.
#' @param genes_list A list containing names of genes of interest.
#' 
#' @return A graph for all genes in `gene_list`.
#' @export

Rhythmic_expression <- function(dryList, FilterGenesTable, genes_list) {
  FilterGenesTable <- read.csv2(file = "~/CIRCANUTR/data/FilterGenesTable.csv")
  
  Genes2plot <- FilterGenesTable[FilterGenesTable$gene_name %in% genes_list, ]
  
  dataGraph <- as.data.frame(log2(dryList$ncounts + 1))
  colnames(dataGraph) <- paste(groupTable[colnames(dataGraph), "Diet"], groupTable[colnames(dataGraph), "ZT"], groupTable[colnames(dataGraph), "sample_name"], sep = ".")
  
  data.frame(ZT = rep(0:24, 2), diet = rep(c("CTRL", "RF"), each = 25)) |>
    left_join(
      dryList$parameters[rownames(dryList$parameters) %in% Genes2plot$gene_id, ] |>
        rownames_to_column("gene_id") |>
        select(-contains("model"), -contains("relamp")) |>
        pivot_longer(
          cols = !gene_id,
          names_to = "nom_cols",
          values_to = "valeurs"
        ) |>
        separate(nom_cols, c("parameter", "diet")) |>
        mutate(parameter = gsub("mean", "MESOR", parameter)) |>
        pivot_wider(
          id_cols = c(gene_id, diet),
          names_from = parameter,
          values_from = valeurs
        ),
      relationship = "many-to-many"
    ) |>
    mutate(Mdl = as.numeric((
      MESOR + a * cos(2 * pi * ZT / 24) +
        b * sin(2 * pi * ZT / 24)
    ))) |>
    mutate(Mdl = ifelse(is.na(Mdl), MESOR, Mdl)) |>
    left_join(
      dataGraph[rownames(dryList$parameters) %in% Genes2plot$gene_id, ] |>
        rownames_to_column("gene_id") |>
        pivot_longer(
          cols = !gene_id,
          names_to = "nom_cols",
          values_to = "valeurs"
        ) |>
        separate(nom_cols, into = c("diet", "ZT", "sample")) |>
        mutate(ZT = as.integer(ZT)) |>
        group_by(gene_id, diet, ZT) |>
        summarize(
          MEAN = mean(valeurs),
          SEM = sd(valeurs) / sqrt(n() - 1)
        ),
      relationship = "many-to-many"
    ) |>
    mutate(diet = factor(diet, levels = c("CTRL", "RF"), ordered = T)) |>
    left_join(Genes2plot[, c("gene_id", "gene_name")], relationship = "many-to-many") |>
    ggplot(aes(x = ZT, y = Mdl)) +
    geom_pointrange(aes(
      x = ZT,
      y = MEAN,
      ymin = MEAN - SEM,
      ymax = MEAN + SEM,
      color = diet
    ),
    linewidth = .25) +
    geom_line(aes(color = diet), size = 1) +
    scale_color_manual(values = c("#007FFF", "#FF7F00")) +
    facet_wrap(facets = vars(gene_name),
               ncol = 2,
               scale = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
    theme(
      text = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      strip.background = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      axis.title.y = element_blank()
    )
}
