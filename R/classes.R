#' The FAME container class
#'
#' @slot genesets A list of containing vectors of genes.
#' @slot graph An igraph object that represents the kNN graph.
#'
#' @aliases FAME
#' @rdname FAME
#' @export
#'

setClass("FAME",
    slots=c(
        genesets = "list",
        graph = "list" # this should be a list or an igraph object
        ),
    prototype = list(
        genesets = list(),
        graph = list()
        )
)