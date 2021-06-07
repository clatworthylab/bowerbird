#' @include utilities.R
#'
#' @export
#' @rdname bower
#'

bower <- function(..., geneset = list(), graph=list()){
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(list(...)) == 0){
        bower <- .foolseldom()
    } else if (length(graph) > 0){
        if (length(geneset) > 0){
            bower <- new('BOWER', geneset = read_geneset(geneset, ...), graph = graph)
        } else {
            bower <- new('BOWER', geneset = list(), graph = graph)
        }         
    } else {
        bower <- new('BOWER', geneset = read_geneset(geneset, ...), graph = list())
    }

    bower
}