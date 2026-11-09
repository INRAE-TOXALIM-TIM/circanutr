#' Circadian Heatmaps
#'
#' @description Generates heatmaps of the expression models proposed by DryR.
#'
#' @param dryList A DryR list containing dataframes with results, rhythmic parameters, normalized counts, raw counts, and the Cook's distance produced by the dryseq() function.
#' @param model An integer for the number of models created with the dryseq() function.
#'
#' @return Heatmaps of each model.
#' @export
#'
#' @examples
#' Examples usages :
#' Circadian_Heatmaps(dryList, model = 1)
#'
#' ggpubr::ggarrange(
#' plotlist = lapply(1:5, function(i) {
#'   Circadian_Heatmaps(dryList = dryList, model = i)
#'   }),
#'    ncol = 1,
#'    heights = 2 * (model_size$model_size / max(model_size$model_size)
#'    )
#' )

Circadian_Heatmaps <- function(dryList, model) {
  # Class genes by their phase
  OrderModL <- dryList$results |>
    filter(chosen_model == model) |>
    arrange(ifelse(is.na(phase_CTRL), phase_RF, phase_CTRL)) |>
    row.names()
  
  # Organize data for the heatmaps
  dataHM <- list()
  dataHM$CTRL <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "CTRL"] # get genes with the same order as OrderModL
  dataHM$RF <- as.data.frame(dryList$ncounts[OrderModL, ])[, dryList$group == "RF"] # get genes with the same order as OrderModL
  
  # Heatmap generation
  ggpubr::ggarrange(plotlist = lapply(c("CTRL", "RF"), function(x) {
    pheatmap::pheatmap(
      dataHM[[x]],
      scale = "row",
      color = colorRampPalette(c("blue", "yellow"))(1000),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      main = paste("model", model, sep = " "),
      legend = FALSE,
      show_colnames = FALSE,
      silent = TRUE
    )[[4]]
  }))
}
