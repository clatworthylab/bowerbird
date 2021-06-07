#' documenting scFA shortcut
#'
#' @name scFA
#' @return for internal use only. quickly redocumenting scFA
#' @examples
#' \donttest{
#' init() # leaves folder, install, and change back to folder
#' }
#' @import devtools
#' @export
init <- function() {
	requireNamespace('devtools')
	devtools::document()
	setwd('..')
	devtools::install('scFA')
	setwd('scFA')
}
