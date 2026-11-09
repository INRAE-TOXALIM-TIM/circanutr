#' Get_gene_id
#'
#' @description Retrieves the gene ID code from the Ensembl database for a given gene name.
#'
#' @param gene_name A character string representing the name of the gene (e.g., "Bmal1").
#'
#' @return The gene ID code from the Ensembl database corresponding to the given gene name.
#' @export
#' 
#' @examples
#' Example usage :
#' gene_name <- "Bmal1"
#' Get_gene_id(gene_name)

Get_gene_id <- function(gene_name) {
  FilterGenesTable$gene_id[FilterGenesTable$gene_name == gene_name]
}
