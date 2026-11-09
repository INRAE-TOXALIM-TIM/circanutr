#' Cumulative_Amp_Plot
#'
#' @description Generates a cumulative plot of the number of rhythmic genes as a function of their amplitude (log2FC) for different diets.
#'
#' @param dryList A list or object containing results from some analysis. It should have a component called `results` which is a data frame with gene names as row names and columns containing amplitude values for different diets.
#'
#' @return A ggplot object that is a cumulative plot of the number of rhythmic genes against their amplitude (log2FC) for different diets.
#' @export
#' 
#' @examples
#' Example usage :
#' dryList <- list(results = data.frame(
#'   amp_LB = c(1.2, 0.5, 2.3),
#'   amp_DT = c(0.8, 1.5, 0.9),
#'   row.names = c("gene1", "gene2", "gene3")
#' ))
#' Cumulative_Amp_Plot(dryList)

Cumulative_Amp_Plot <- function(dryList) {
  # graph cumulative number of rhythmic genes = function of amplitude
  dryList$results |>
    rownames_to_column("gene_name") |>
    select(gene_name, contains("amp_") & !contains("relamp_")) |>
    pivot_longer(!gene_name, names_to = "Diet", values_to = "log2FC") |>
    group_by(Diet) |>
    filter(!is.na(log2FC)) |>
    arrange(log2FC) |>
    mutate(`number of genes` = order(rev(log2FC)),
           Diet = gsub("amp_", "", Diet)) |>
    ggplot(aes(
      x = log2FC,
      y = log10(`number of genes`),
      colour = Diet
    )) +
    xlab("Amplitude (log2FC)") +
    ylab("log10(rhythmic genes)") +
    geom_line() +
    scale_color_manual(values = c("#007FFF", "#FF7F00"), labels = c("LB", "DT")) +
    theme_classic()
}
