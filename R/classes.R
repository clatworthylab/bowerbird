#' The BOWER container class
#'

#' @title genesets-class
#' @rdname genesets-class
#' @description list or character or data.frame
setClassUnion("genesets", c("list", "character", "data.frame"))

#' @title clusters-class
#' @rdname clusters-class
#' @description character or factor or NULL
setClassUnion("clusters", c("character", "factor", "numeric", "NULL"))

setClass("igraph")
#' @title graph-class
#' @rdname graph-class
#' @description igraph or NULL
setClassUnion("graph", c("igraph", "NULL"))

#' @title coregenes-class
#' @rdname coregenes-class
#' @description list or NULL
setClassUnion("coregenes", c("list", "NULL"))


setClass("layout_tbl_graph")
setClass("layout_ggraph")
setClassUnion("hidden", c("data.frame", "layout_tbl_graph", "layout_ggraph", "NULL"))

#' @slot genesets A list of containing vectors of genes.
#' @slot graph An igraph object that represents the kNN graph.
#' @slot clusters A vector holding the cluster labels for each geneset.
#' @slot .graph_data Hidden slot for graph in dataframe format.
#'
#' @title BOWER
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
        .graph_data = "hidden" # hidden slot for graph in dataframe format
        ),
    prototype = list(
        genesets = list(),
        graph = NULL,
        clusters = NULL,
        coregenes = NULL,
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
