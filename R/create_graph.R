#' @title create_graph
#'

#' Creates a k-nearest-neighbor graph from a tf-idf matrix.
#'
#' @param snn KNN graph.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See ?igraph::graph_from_adjacency_matrix.
#' @param weighted This argument specifies whether to create a weighted graph from an adjacency matrix. If it is NULL then an unweighted graph is created and the elements of the adjacency matrix gives the number of edges between the vertices. If it is a character constant then for every non-zero matrix entry an edge is created and the value of the entry is added as an edge attribute named by the weighted argument. If it is TRUE then a weighted graph is created and the name of the edge attribute will be weight. See ?igraph::graph_from_adjacency_matrix.
#' @param ... passed to igraph::graph_from_adjacency_matrix.
#' @details
#' Given a list of text, it creates a sparse matrix consisting of tf-idf score for tokens from the text. See `https://github.com/saraswatmks/superml/blob/master/R/TfidfVectorizer.R`. A k shortest-nearest neighbor graph is then computed using the overlap of of the terms.
#' @return Returns a matrix of tf-idf score of tokens.
#' @examples
#' gr <- create_graph(mat)
#' @export
#'

create_graph <- function(snn, mode = 'undirected', weighted = TRUE, ...) {
	requireNamespace('igraph')
	gr <- igraph::graph_from_adjacency_matrix(snn, mode= "undirected", weighted = TRUE, ...)
	gr <- igraph::simplify(gr)
	return(gr)
}