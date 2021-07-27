#' @include utilities.R

#' Summarize the core genes for main terms.
#' 
#' @title extract_core
#' @description Summarize the outersection genes for the. summary terms.
#' @name extract_core
#' @param bower processed bower object.
#' @param inplace whether or not to return an updated BOWER class or return the output as a list.
#' @examples
#' gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
#' bwr <- bower(gmt_file)
#' bwr <- snn_graph(bwr)
#' bwr <- find_clusters(bwr)
#' bwr <- summarize_clusters(bwr)
#' extract_core(bwr)
#' @return Returns a list of genes that outersect geneset clusters in coregenes slot or as a list.
#' @import dplyr

#' @export
extract_core.BOWER <- function(bower, inplace = TRUE){
	gs <- split(bower@.graph_data, bower@.graph_data$cluster)
	gs <- lapply(gs, function(x) x %>% dplyr::select(name) %>% unlist %>% as.character)	
	# find which cluster only has 1 geneset
	separate <- which(lapply(gs, length) < 2)
	gs <- lapply(gs, function(x) bower@genesets[x])
	# remove temporarily
	if (length(separate) > 0){
		gs_2 <- gs[-separate]
		gs_2 <- lapply(gs_2, function(x){
			ll <- combn(x, 2, simplify = FALSE)
			out <- lapply(ll, function(z) intersect(z[[1]], z[[2]]))
			out <- unique(do.call(c, out))
			return(out)	
		})
		for (i in seq_along(gs[separate])){
			gs_2[separate[i]] <- gs[separate][i][[1]]
		}
	} else {
		gs_2 <- lapply(gs, function(x){
			ll <- combn(x, 2, simplify = FALSE)
			out <- lapply(ll, function(z) intersect(z[[1]], z[[2]]))
			out <- unique(do.call(c, out))
			return(out)	
		})
	}
	summarydat <- bower@.graph_data[,c('cluster', 'labels')]
	summarydat <- summarydat[!is.na(summarydat$labels), ]
	summarydat <- summarydat[order(summarydat$cluster), ]
	names(gs_2) <- summarydat$labels
	if (inplace){
		bower@coregenes <- gs_2
		return(bower)
	} else {
		return(gs_2)
	}
}
