#' @include utilities.R
#'
#' @export
#' @rdname bower
#'

bower <- function(genesets = list(), graph=list(), clusters = c(), ...){
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
                out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = graph, clusters = clusters)
            } else {
                out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = graph, clusters = c())
            }            
        } else {
            if (length(clusters) > 0){
                out <- new('BOWER', genesets = list(), graph = graph, clusters = clusters)
            } else {
                out <- new('BOWER', genesets = list(), graph = graph, clusters = c())
            }
        }
    } else {
        if (length(clusters) > 0){
            out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = list(), clusters = clusters)
        } else {
            out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = list(), clusters = c())
        }
    }

    out
}