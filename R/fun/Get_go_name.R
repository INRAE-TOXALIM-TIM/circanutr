#' Get GO Name
#'
#' @description Retrieves the name of a Gene Ontology (GO) term given its GO ID.
#'
#' @param go_id A character string representing the GO ID of the desired gene, such as "GO:0006631".
#' 
#' @return A character string representing the GO name from the Gene Ontology database, or NA if the GO ID is not found.
#' @export
#' 
#' @examples
#' # Example usage:
#' # Get_go_name("GO:0006631")

Get_go_name <- function(go_id) {
  pre_id <- GOfuncR::get_names(go_id)
  if (is.na(pre_id$go_name)) {
    output <- NA
  } else {
    GOtable <- GOfuncR::get_names(go_id)
    output <- GOtable$go_name[GOtable$go_id == go_id]
  }
  
  return(output)
}
