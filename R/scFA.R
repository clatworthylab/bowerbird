#' shortcut to initiate building of scFA
#'
#' @name scFA
#' @return for internal use only. runs init scFA
#' @examples
#' \donttest{
#' scFA()
#' }
#' @export
scFA <- function()
{
	setwd("~/Documents/GitHub/scFA")
	requireNamespace('scFA')
	requireNamespace('roxygen2')
	init_scFA()
}