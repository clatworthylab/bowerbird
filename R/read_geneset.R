#' @title read_geneset
#' 
#' @include utilities.R
#' @include generics.R
#' 

#' Reads various geneset file formats and returns a list of vectors containing gene names.
#'
#' @name read_geneset
#' @param x file,
#' @param min_size minimum size of geneset otherwise filter out.
#' @return Returns a list of gene sets
#' @examples
#' \donttest{
#' gs <- read_geneset('geneset.gmt')
#' }
#' @export
#'

read_geneset <- function(x, min_size = 15) {
    if (class(x) == 'character'){
        if (tools::file_ext(x) == 'gmt') {
            class(x) <- c(class(x), 'gmt')
        } else if (tools::file_ext(x) == 'gmx') {
            class(x) <- c(class(x), 'gmx')
        } else if (tools::file_ext(x) == 'csv') {
            class(x) <- c(class(x), 'csv')
        } else if (tools::file_ext(x) == 'tsv') {
            class(x) <- c(class(x), 'tsv')
        } else if (tools::file_ext(x) == 'txt') {
            class(x) <- c(class(x), 'tsv')
        }        
    }
    return(.read_geneset(x, min_size))
}
