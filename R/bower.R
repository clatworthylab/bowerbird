#' @include utilities.R
#'
#' @export
#' @rdname bower
#'

bower <- function(genesets = list(), graph=list(), ...){
    requireNamespace('S4Vectors')
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(genesets) == 0 & length(graph) == 0){
        out <- .emptyBOWER()
    } else if (length(graph) > 0){
        if (length(geneset) > 0){
            out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = graph)
        } else {
            out <- new('BOWER', genesets = list(), graph = graph)
        }
    } else {
        out <- new('BOWER', genesets = read_geneset(genesets, ...), graph = list())
    }

    out
}