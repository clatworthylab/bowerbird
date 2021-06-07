#' The BOWER container class
#'
#' @slot genesets A list of containing vectors of genes.
#' @slot graph An igraph object that represents the kNN graph.
#'
setClassUnion("listORcharacterORdataframe", c("list", "character", "data.frame"))
#' @aliases BOWER
#' @rdname BOWER
#' @export
#'

setClass("BOWER",
    slots=c(
        genesets = "listORcharacterORdataframe",
        graph = "list" # this should be a list or an igraph object
        ),
    prototype = list(
        genesets = list(),
        graph = list()
        )
)