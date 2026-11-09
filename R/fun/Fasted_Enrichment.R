Fasted_Enrichment <- function(model_DF, model_FAFR) {
  gene.list.FA <- Table.Fasting$Row.names[Table.Fasting$chosen_model == model_DF &
                                            Table.Fasting$comment == model_FAFR]
  
  gene.list.FR <- Table.Refed$Row.names[Table.Refed$chosen_model == model_DF &
                                          Table.Refed$comment == model_FAFR]
  
  graph <- venn::venn(list(FA = gene.list.FA, FR = gene.list.FR), ilabels = "counts")
  
  Enrich.list <- list(
    intersect_FAFR = intersect(gene.list.FA, gene.list.FR),
    FA_exclusive = setdiff(gene.list.FA, gene.list.FR),
    FR_exclusive = setdiff(gene.list.FR, gene.list.FA)
  )
  
  enrichment <- lapply(Enrich.list, function(x)
    clusterProfiler::enrichGO(
      gene = x,
      keyType = "ENSEMBL",
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      ont = "BP",
      pvalueCutoff =  0.05,
      pAdjustMethod = "BH",
      universe = intersect(Table.Fasting$Row.names, Table.Refed$Row.names)
    ))
  
  
  return(list(
    FA = gene.list.FA,
    FR = gene.list.FR,
    graph = graph,
    enrichment = enrichment
  ))
  
}
