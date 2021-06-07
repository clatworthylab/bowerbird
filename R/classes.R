#' The BOWER container class
#'
#' @slot genesets A list of containing vectors of genes.
#' @slot graph An igraph object that represents the kNN graph.
#'
setClassUnion("listORcharacterORdataframe", c("list", "character", "data.frame"))
setClassUnion("characterORfactorORNULL", c("character", "factor", "NULL"))
setClass("igraph")
setClassUnion("listORigraph", c("list", "igraph", "NULL"))
#' @aliases BOWER
#' @rdname BOWER
#' @export
#'

setClass("BOWER",
    slots=c(
        genesets = "listORcharacterORdataframe",
        graph = "listORigraph", # this should be a list or an igraph object
        clusters = 'characterORfactorORNULL'
        ),
    prototype = list(
        genesets = list(),
        graph = list(),
        clusters = c()
        )
)