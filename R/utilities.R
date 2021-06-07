#########################
# Miscellaneous functions
#########################

#' @param x vector.
#' 

.is_blank <- function(x) x == ""


#' @param gs list of vectors.
#' @param min_size minimum size of geneset otherwise filter out.
#' 

.filter_geneset <- function(gs, min_size) {
    gs <- lapply(gs, function(x) x[!.is_blank(x)])
    if (length(which(is.na(gs))) > 0){
        gs <- gs[!is.na(gs)]
    }
    if (length(which(lapply(gs, length) <= min_size)) > 0) {
        gs <- gs[-which(lapply(gs, length) <= min_size)]
    }
    return(gs)
}


#' @param ... does nothing.
#' 

.emptyBOWER <- function(...) {
    out <- new("BOWER",
               genesets=list(),
               graph=list())
    out
}


#' @param snn KNN graph.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See ?igraph::graph_from_adjacency_matrix.
#' @param weighted This argument specifies whether to create a weighted graph from an adjacency matrix. If it is NULL then an unweighted graph is created and the elements of the adjacency matrix gives the number of edges between the vertices. If it is a character constant then for every non-zero matrix entry an edge is created and the value of the entry is added as an edge attribute named by the weighted argument. If it is TRUE then a weighted graph is created and the name of the edge attribute will be weight. See ?igraph::graph_from_adjacency_matrix.
#' @param ... passed to igraph::graph_from_adjacency_matrix.
#' @return Returns a matrix of tf-idf score of tokens.
#'

.create_graph <- function(snn, mode = 'undirected', weighted = TRUE, ...) {
    requireNamespace('igraph')
    gr <- igraph::graph_from_adjacency_matrix(snn, mode= "undirected", weighted = TRUE, ...)
    gr <- igraph::simplify(gr)
    return(gr)
}


#' @import devtools
#' @export

init <- function()
{   
    setwd("~/Documents/GitHub/bowerbird")
    requireNamespace('roxygen2')
    requireNamespace('devtools')
    devtools::document()
    setwd('..')
    devtools::install('bowerbird')
    setwd('bowerbird')
}

