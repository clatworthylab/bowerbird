#' documenting scFA shortcut
#'
#' @name scFA
#' @return for internal use only. quickly redocumenting scFA
#' @examples
#' \donttest{
#' init_scFA() # leaves folder, install, and change back to folder
#' }
#' @import devtools
#' @export
init_scFA <- function() {
	requireNamespace('devtools')
	devtools::document()
	setwd('..')
	devtools::install('scFA')
	setwd('scFA')
}
