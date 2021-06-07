#########################
# Miscellaneous functions
#########################

#' @param x vector.
#' 
.is_blank <- function(x) x == ""

#' @param gs list of vectors.
#' @param min_size minimum size of geneset otherwise filter out.
#' 
.filter_geneset <- function(gs, min_size) {
    gs <- lapply(gs, function(x) x[!.is_blank(x)])
    if (length(which(is.na(gs))) > 0){
        gs <- gs[!is.na(gs)]
    }
    if (length(which(lapply(gs, length) <= min_size)) > 0) {
        gs <- gs[-which(lapply(gs, length) <= min_size)]
    }
    return(gs)
}

#' @name misc
#' @import devtools
#' @export

init <- function()
{   
    setwd("~/Documents/GitHub/scFAME")
    requireNamespace('roxygen2')
    requireNamespace('devtools')
    devtools::document()
    setwd('..')
    devtools::install('scFAME')
    setwd('scFAME')
}

