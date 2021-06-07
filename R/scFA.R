#' shortcut to initiate building of scFA
#'
#' @return for internal use only. quickly redocumenting scFA
#' @examples
#' \donttest{
#'init_scFA()
#' }
#' @import devtools
#' @export

.init_scFA <- function() {
	requireNamespace('devtools')
	devtools::document()
	setwd('..')
	devtools::install('scFA')
	setwd('scFA')
}

init_scFA <- function()
{
	setwd("~/Documents/GitHub/scFA")
	requireNamespace('roxygen2')
	.init_scFA()
}

