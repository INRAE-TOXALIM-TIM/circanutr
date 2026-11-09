#' HyperT_clust
#'
#' @description Hyper-geometric test for clustering identified GO categories groups between drylist models.
#'
#' @param dryList A list or object containing results from some analysis. It should have a component called `results` which itself has a componentcalled `chosen_model`.
#' @param Background A vector or list of background genes or terms againstwhich the hypergeometric test is performed.
#' @param GOgenes A list of gene ontology (GO) categories or terms that have been identified as significant or of interest.
#'
#' @return A data frame with columns: cluster (cluster ID), m (number of background genes in GO categories), n (number of background genes not in GO categories), k (number of genes in each cluster that are in the background), q (number of genes in each cluster that are in the GO categories), pval (p-value from the hypergeometric test), go_id (GO category identifiers).
#' @export
#' 
#' @examples
#' Example usage :

HyperT_clust <- function(dryList, Background, GOgenes) {
  HT <- data.frame(
    cluster = 1:max(dryList$results$chosen_model),
    m = sum(unlist(GOgenes) %in% Background)
  )
  HT$n <- length(Background) - HT$m
  HT$k <- sapply(HT$cluster, function(x) {
    sum(dryList$results$chosen_model[rownames(dryList$results) %in% Background] == x)
  })
  HT$q <- sapply(HT$cluster, function(x) {
    sum(rownames(dryList$results)[dryList$results$chosen_model == x] %in% unlist(GOgenes))
  })
  HT$pval <- phyper(HT$q, HT$m, HT$n, HT$k, lower.tail = F)
  HT$go_id <- names(GOgenes)
  return(HT)
}
