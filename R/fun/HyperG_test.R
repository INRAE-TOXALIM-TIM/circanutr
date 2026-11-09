#' HyperG_test
#'
#' @description Performs a hypergeometric test to identify enriched GO categories among rhythmic genes from RNA-seq data.
#'
#' @param dryList A list or object containing results from RNA-seq analysis. It should have a component called `results` which contains a data frame with gene names as row names and a column named `phase_CTRL` where non-NA values indicate rhythmic genes.
#' @param GOgenes A list of gene ontology (GO) categories. Each element in the list is a vector of gene names belonging to that GO category.
#'
#' @return A data frame with columns: GOname (name of the GO category), GOsize (number of genes in the GO category present in the background), n_rhythmic (number of genes in the GO category that are rhythmic), p_value (p-value from the hypergeometric test).
#' @export
#' 
#' @examples
#' Example usage :
#' dryList <- list(results = data.frame(phase_CTRL = c(NA, 1, NA, 2), row.names = c("gene1", "gene2", "gene3", "gene4")))
#' GOgenes <- list(go1 = c("gene1", "gene3"), go2 = c("gene2", "gene4"))
#' HyperG_test(dryList, GOgenes)

HyperG_test <- function(dryList, GOgenes) {
  Background <- row.names(dryList$results)[row.names(dryList$results) %in% unique(unlist(GOgenes))]
  All_Rhythmic_genes <- Background[Background %in% row.names(dryList$results[!is.na(dryList$results$phase_CTRL), ])] # rhythmic genes
  
  m <- length(All_Rhythmic_genes) # the number of white balls in the urn
  n <- length(Background[!Background %in% All_Rhythmic_genes]) # the number of black balls in the urn
  
  Rhythmic_GOcat <- sapply(names(GOgenes), function(x) {
    k <- sum(GOgenes[[x]] %in% Background) 	# the number of balls drawn from the urn ie nb of gene from GO cat in the RNAseq
    q <- sum(GOgenes[[x]] %in% All_Rhythmic_genes)
    data.frame(
      GOname = Get_go_name(x),
      GOsize = k,
      n_rhythmic = q,
      p_value = phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
    )
  })
  Rhythmic_GOcat <- t(Rhythmic_GOcat)
  Rhythmic_GOcat <- as.data.frame(Rhythmic_GOcat)
  
  return(Rhythmic_GOcat)
}
