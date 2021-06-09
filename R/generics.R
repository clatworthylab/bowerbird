#' @param x .gmx, .gmt, .csv, .tsv, .txt, data.frame or list object containing gene sets.
#' @param ... Arguments passed to other methods.
#
#' @return Returns a list of gene sets.
#'

.read_geneset <- function(x, ...) {
    UseMethod(generic = "read_geneset")
}

#' @param x either a BOWER class or a list.
#' @param ... Arguments passed to other methods.
#
#' @return Returns a SNN geneset graph.
#' @export
#'

snn_graph <- function(x, ...) {
    UseMethod(generic = "snn_graph")
}

#' @param x either a BOWER class or an igraph object.
#' @param ... Arguments passed to other methods.
#
#' @return Summarizes labels for genesets.
#' @export
#'

summarize_clusters <- function(x, ...) {
    UseMethod(generic = "summarize_clusters")
}

#' @param x either a BOWER class or an igraph object.
#' @param ... Arguments passed to other methods.
#
#' @return Return leiden clustering assignment for each geneset.
#' @export
#'

find_clusters <- function(x, ...) {
    UseMethod(generic = "find_clusters")
}

#' @param x BOWER class.
#' @param ... Arguments passed to other methods.
#
#' @return Set cluster labels for each gene set
#' @export
#'

set_clusters <- function(x, ...) {
    UseMethod(generic = "set_clusters")
}