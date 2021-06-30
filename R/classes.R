#' The BOWER container class
#'

#' @title BOWER
#' @description
#' The BOWER class for Functional annotation + Gene Module Summarization

#' @rdname BOWER
setClassUnion("genesets", c("list", "character", "data.frame"))

#' @rdname BOWER
setClassUnion("clusters", c("character", "factor", "numeric", "NULL"))

setClass("igraph")
#' @rdname BOWER
setClassUnion("graph", c("igraph", "NULL"))

#' @rdname BOWER
setClassUnion("coregenes", c("list", "NULL"))

#' @rdname BOWER
setClassUnion("scores", c("data.frame", "matrix", "list", "NULL"))

setClass("layout_tbl_graph")
setClass("layout_ggraph")
setClassUnion("hidden", c("data.frame", "layout_tbl_graph", "layout_ggraph", "NULL"))

#' @slot genesets A list of containing vectors of genes.
#' @slot graph An igraph object that represents the kNN graph.
#' @slot clusters A vector holding the cluster labels for each geneset.
#' @slot coregenes A list of containing vectors of core genes.
#' @slot scores A dataframe or list containg the geneset testing results.
#'
#' @aliases BOWER
#' @rdname BOWER
#' @export
#'
setClass("BOWER",
    slots=c(
        genesets = "genesets",
        graph = "graph", # this should be a list or an igraph object
        clusters = 'clusters',
        coregenes = "coregenes",
        scores = "scores", # not hidden but use the same class
        .graph_data = "hidden" # hidden slot for graph in dataframe format
        ),
    prototype = list(
        genesets = list(),
        graph = NULL,
        clusters = NULL,
        coregenes = NULL,
        scores = NULL,
        .graph_data = NULL
        )
    )

setMethod("show","BOWER",function(object) {
    cat("BOWER class\n")
    cat("number of genesets: ", length(object@genesets),"\n")
    if (!is.null(object@graph)){
        cat("genesets kNN Graph: \n")
        print(object@graph)
    }
    if (!is.null(object@clusters)){
        cat("number of geneset clusters: ", length(unique(object@clusters)),"\n")
    }
    if (!is.null(object@coregenes)){
        cathead <- function(x){
            return(cat(head(x), "...\n"))
        }
        cat("Core genes:\n")
        cat("\tFirst six genes shown\n")
        for (i in seq_along(object@coregenes)){
            cat("\t",names(object@coregenes)[i], ": ")
            cathead(object@coregenes[[i]])
        }
    }
})
