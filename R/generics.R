#' @rdname misc
.read_geneset <- function(x, ...) {
    UseMethod(generic = "read_geneset")
}


#' @export
#'

snn_graph <- function(...) {
    UseMethod(generic = "snn_graph")
}


#' @export
#'

summarize_clusters <- function(...) {
    UseMethod(generic = "summarize_clusters")
}


#' @export
#'

find_clusters <- function(...) {
    UseMethod(generic = "find_clusters")
}


#' @export
#'

set_clusters <- function(...) {
    UseMethod(generic = "set_clusters")
}


#' @export
#'

extract_core <- function(...) {
    UseMethod(generic = "extract_core")
}


#' @export
#'

enrich_genesets <- function(...) {
    UseMethod(generic = "enrich_genesets")
}

#' @export
#'

plot_graph <- function(...) {
    UseMethod(generic = "plot_graph")
}