#' @title bower
#' @name bower
#' @rdname bower
#' @description
#' A flexible function to initiate a BOWER class from a geneset file or a list of genesets.
#' @include utilities.R
#' @examples
#' gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
#' bwr <- bower(gmt_file)
#' bwr
#' @export
#'

bower <- function(genesets, graph=NULL, clusters = NULL, ...){
    requireNamespace('S4Vectors')
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(genesets) == 0 & length(graph) == 0){
        out <- .emptyBOWER()
    } else if (length(graph) > 0){
        if (length(genesets) > 0){
            if (length(clusters) > 0){
                out <- new('BOWER', genesets = .sanitize(read_geneset(genesets, ...)), graph = graph, clusters = clusters, ...)
            } else {
                out <- new('BOWER', genesets = .sanitize(read_geneset(genesets, ...)), graph = graph, clusters = NULL, ...)
            }            
        } else {
            if (length(clusters) > 0){
                out <- new('BOWER', genesets = list(), graph = graph, clusters = clusters, ...)
            } else {
                out <- new('BOWER', genesets = list(), graph = graph, clusters = NULL, ...)
            }
        }
    } else {
        if (length(clusters) > 0){
            out <- new('BOWER', genesets = .sanitize(read_geneset(genesets, ...)), graph = NULL, clusters = clusters, ...)
        } else {
            out <- new('BOWER', genesets = .sanitize(read_geneset(genesets, ...)), graph = NULL, clusters = NULL, ...)
        }
    }

    if (length(out@graph) > 0){
        out@.graph_data <- .graph_to_data(out@graph)
    }
    
    return(out)
}
