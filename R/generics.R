#' @param x .gmx, .gmt, .csv, .tsv, .txt, data.frame or list object containing gene sets.
#' @param ... Arguments passed to other methods
#
#' @return Returns a list of gene sets.
#'

.read_geneset <- function(x, ...) {
    UseMethod(generic = "read_geneset")
}