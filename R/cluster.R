#' @title clusters
#' @include utilities.R

#' Clusters SNN graph using leiden algorithm.
#'
#' @name find_clusters
#' @param bower BOWER class..
#' @param gr igraph object.
#' @param resolution value for leiden clustering.
#' @param ... passed to leiden::leiden.
#' @description
#' Performs clustering of SNN graph of genesets.
#' @return Cluster assignment for each geneset.
#' @examples
#' \donttest{
#' bwr <- find_clusters(bwr)
#' }
#' @import leiden
#' @export

find_clusters.BOWER <- function(bower, resolution = 3, ...){
    if (!is.null(bower@graph)){
        clusters <- leiden::leiden(bower@graph, resolution = 3, ...)
        igraph::V(bower@graph)$cluster <- clusters
        if (length(bower@genesets) > 0){
            ds <- data.frame(geneset_size = unlist(lapply(bower@genesets, length)))
            idx <- match(igraph::V(bower@graph)$name, row.names(ds))
            igraph::V(bower@graph)$geneset_size <- ds$geneset_size[idx]
        }
        bower@clusters <- clusters
        bower@.graph_data <- .graph_to_data(bower@graph)
        return(bower)
    } else {
        stop('Graph slot not found. Please run snn_graph first.')
    } 
}

#' @name find_clusters
#' @export
find_clusters.igraph <- function(gr, resolution = 3, ...){
    clusters <- leiden::leiden(gr, resolution = 3, ...)    
    return(clusters)
}

#' Set cluster to BOWER class
#'
#' @name set_clusters
#' @param bower bower class
#' @param clusters Cluster assignment for each geneset.
#' @description
#' Manually sets externally/alternatively determined cluster assignment for each geneset.
#' @return bower class
#' @examples
#' \donttest{
#' bower <- set_clusters(bower, c(1,1,1,2,2,2,3,3,3))
#' }
#' @export

set_clusters.BOWER <- function(bower, clusters){    
    bower@clusters <- clusters
    if (!is.null(bower@graph)){
        igraph::V(bower@graph)$cluster <- clusters
        if (length(bower@genesets) > 0){
            ds <- data.frame(geneset_size = unlist(lapply(bower@genesets, length)))
            idx <- match(igraph::V(bower@graph)$name, row.names(ds))
            igraph::V(bower@graph)$geneset_size <- ds$geneset_size[idx]
        }
        bower@.graph_data <- .graph_to_data(bower@graph)
    }    
    return(bower)
}
