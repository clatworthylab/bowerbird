#' @export
#' @rdname scFAME
#' 
fame <- function(..., geneset = list(), graph=list()){
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(list(...)) == 0){
        fame <- .emptyFAN()
    }

    fame
}