do_network <- function(threshold,
                       go_terms,
                       GO_data,
                       plot.title = NULL,
                       go_size = NULL,
                       go_color = NULL) {
  if (length(go_terms) > 2) {
    go_network <- igraph::make_empty_graph(n = length(go_terms), directed = FALSE)
    V(go_network)$name <- go_terms
    
    # Calcul de la matrice de similarité
    my.Sim.matrix <- sapply(go_terms, function(x) {
      GOSemSim::goSim(x, go_terms, GO_data, measure = "Wang")
    }) |> as.data.frame()
    
    my.Sim.matrix <- my.Sim.matrix[, !is.na(my.Sim.matrix[1, ])]
    row.names(my.Sim.matrix) <- colnames(my.Sim.matrix)
    
    for (i in 1:(length(go_terms) - 1)) {
      for (j in (i + 1):length(go_terms)) {
        if (my.Sim.matrix[i, j] > threshold) {
          go_network <- igraph::add_edges(go_network, edges = c(i, j))
        }
      }
    }
    
    term_df <- AnnotationDbi::select(
      GO.db::GO.db,
      keys = go_terms,
      columns = "TERM",
      keytype = "GOID"
    )
    
    V(go_network)$label <- term_df$TERM[match(V(go_network)$name, term_df$GOID)]
    
    
    if (!is.null(go_size)) {
      go_network <- go_network  |>
        set_vertex_attr('size', index = V(go_network), value = go_size)
    } else {
      V(go_network)$size <- 5
    }
    
    if (!is.null(go_color)) {
      go_network <- go_network |>
        set_vertex_attr('color',
                        index = V(go_network),
                        value = go_color)
    } else {
      V(go_network)$color <- "steelblue" # couleur par défaut
    }
    
    if (!is.null(go_color)) {
      graph.out <- ggraph(go_network, layout = "fr") +
        geom_edge_link(alpha = 0.5) +
        geom_node_point(aes(size = size, color = color), alpha = 0.8) +
        geom_node_text(aes(label = label), repel = TRUE) +
        theme_void() +
        labs(title = plot.title)
    } else {
      graph.out <- ggraph(go_network, layout = "fr") +
        geom_edge_link(alpha = 0.5) +
        geom_node_point(aes(size = size),
                        color = "steelblue",
                        alpha = 0.8) +
        geom_node_text(aes(label = label), repel = TRUE) +
        theme_void() +
        labs(title = plot.title)
    }
    
    
    comp <- igraph::components(go_network)
    df_comp <- data.frame(
      GOID = V(go_network)$name,
      TERM = V(go_network)$label,
      node = comp$membership
    )
    
    if (!is.null(go_size)) {
      df_comp$go_size = go_size
    }
    
    if (!is.null(go_color)) {
      df_comp$go_color = go_color
    }
    
    return(list(network = graph.out, results = df_comp))
    
  } else {
    return(NULL)
  }
}
