#' @include utilities.R
#'
#' @export
#' @rdname scfamous
#'

fame <- function(..., geneset = list(), graph=list()){
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(list(...)) == 0){
        fame <- .emptyFAME()
    } else if (length(graph) > 0){
        if (length(geneset) > 0){
            fame <- new('FAME', geneset = read_geneset(geneset, ...), graph = graph)
        } else {
            fame <- new('FAME', geneset = list(), graph = graph)
        }         
    } else {
        fame <- new('FAME', geneset = read_geneset(geneset, ...), graph = list())
    }

    fame
}