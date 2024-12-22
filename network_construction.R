run.SCINET.clusters <- function(
  ace,
  specificity.slot.name,
  G = NULL,
  min.edge.weight = 2,
  spec.sample_no = 1000,
  thread_no = 8,
  compute.topo.specificity = TRUE
) {

    ACTIONetExperiment:::.check_and_load_package("SCINET")

    print("Preprocessing the baseline interactome")
    if (is.null(G)) {
        if (!exists("PCNet")) {
            data("PCNet")
        }
        Adj = PCNet
    } else if (is.matrix(G) | ACTIONetExperiment:::is.sparseMatrix(G)) {
        Adj = as(G, "sparseMatrix")
        Adj@x = rep(1, length(Adj@x))
    } else if (igraph::is.igraph(G)) {
        Adj = as(igraph::get.adjacency(G), "sparseMatrix")
    }


    if (!(specificity.slot.name %in% names(rowMaps(ace)))) {
        message(sprintf("%s does not exist in rowMaps(ace)", specificity.slot.name))
    }

    gene.scores = as.matrix(log1p(rowMaps(ace)[[specificity.slot.name]]))

    common.genes = intersect(rownames(gene.scores), rownames(Adj))
    cat(length(common.genes), "genes found between scRNA-seq and given network\n")
    if (length(common.genes) == 0) {
        print("No common genes found. Check rownames (or vertex names) for the input graph")
        return(ace)
    }
    A = gene.scores[common.genes, ]
    G = Adj[common.genes, common.genes]


    print("Constructing networks")
    gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
    cellstate.nets = SCINET::construct_cell_networks(
      net = G,
      gene_activities = gene.activity.scores,
      thread_no = thread_no
    )
    cellstate.nets.list = as.list(cellstate.nets)

    print("Post-processing networks\n")
    cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
        G.Adj@x[G.Adj@x < min.edge.weight] = 0
        filter.mask = Matrix::colSums(G.Adj) == 0

        G = igraph::graph_from_adjacency_matrix(
          adjmatrix = G.Adj[!filter.mask, !filter.mask],
          mode = "undirected",
          weighted = TRUE
        )

        V(G)$name = common.genes[!filter.mask]
        if (compute.topo.specificity == TRUE) {
            z.scores = SCINET::topo.spec(G, spec.sample_no)
            V(G)$specificity = 1/(1 + exp(-z.scores))
        }

        return(G)
    })

    if (is.null(colnames(gene.scores))) {
        names(cellstate.nets.list.igraph) = 1:ncol(gene.scores)
    } else {
        names(cellstate.nets.list.igraph) = colnames(gene.scores)
    }
    return(cellstate.nets.list.igraph)
}

