#' Venn Diagram of Rhythmic Genes
#'
#' @description Creates a Venn diagram (technically an Euler diagram)
#' of rhythmic genes for Ad Libitum and Day Fed diets.
#'
#' @param dryList A DryR object containing the models and related data.
#' 
#' @return A plot of the common and exclusive rhythmic genes for the Adlib and DTF diets.
#' @export
#' 
#' @examples
#' Example usage:
#' Venn_Rhythmic(dryList)

Venn_Rhythmic <- function(dryList) {
  # List of genes
  LB_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(2, 4, 5), ])
  DT_rhythmic <- rownames(dryList$parameters[dryList$parameters$chosen_model %in% c(3, 4, 5), ])
  
  # Plot
  ggvenn::ggvenn(list(`LB rhythmic` = LB_rhythmic, `DT rhythmic` = DT_rhythmic),
                 auto_scale = TRUE)
}
