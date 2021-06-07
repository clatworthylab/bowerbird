#' @title summarize_clusters
#'

#' Summarize the terms of the cluster using pagerank algorithm
#'
#' @param graph geneset overlap graph.
#' @param cluster vector of cluster labels for each geneset.
#' @param pattern search pattern to remove from the terms. Unless specified, will default to built-in pattern. 
#' @param sep separator used/found in gene set names to be changed to blank spaces. Default value is underscore ('_').
#' @param ... passed to igraph::graph_from_adjacency_matrix.
#' @details
#' Given a list of text, it creates a sparse matrix consisting of tf-idf score for tokens from the text. See `https://github.com/saraswatmks/superml/blob/master/R/TfidfVectorizer.R`. A k shortest-nearest neighbor graph is then computed using the overlap of of the terms.
#' @return Returns a matrix of tf-idf score of tokens.
#' @examples
#' terms <- summarize_clusters(mat)
#' @import textrank udpipe dplyr pbmcapply
#' @export
#'

summarize_clusters <- function(graph, cluster, pattern = NULL, sep = NULL, ncpus = NULL){	
  requireNamespace('igraph')
	if (is.null(pattern)){
		pattern = '^GO_|^KEGG_|^REACTOME_|^HALLMARK_|POSITIVE_|NEGATIVE_|REGULATION_OF|^GOBP_'
	} 
	if (is.null(pattern)){
		sep
	}
	if (!is.null(cluster)){
		igraph::V(graph)$cluster <- cluster
	}
	if (!is.null(ncpus)){
		require(parallel)
		n_cpus = parallel::detectCores()
	}


	df <- data.frame(name = igraph::V(graph)$name, cluster = igraph::V(graph)$cluster)
	df$name <- gsub(pattern, '', df$name)
	df_split <- split(df, df$cluster)
	df_split <- pbmclapply(df_split, function(x) {
		x <- x %>% select(name) %>% unlist %>% as.character
		x <- gsub(sep, ' ', x)
		return(x)
	})

	tagger <- udpipe_download_model("english", model_dir = tempdir())
	tagger <- udpipe_load_model(tagger$file_model)
	res <- pbmclapply(df_split, function(x) {
		annt <- udpipe_annotate(tagger, paste(x, collapse = '.\n'))
		annt <- as.data.frame(annt)
		annt$textrank_id <- unique_identifier(annt, c("doc_id", "paragraph_id", "sentence_id"))
		sentences <- unique(annt[, c("textrank_id", "sentence")])
		if (nrow(sentences) > 1){
      terminology <- annt[, c("textrank_id", "lemma")]  
      tr <- textrank_sentences(data = sentences, terminology = terminology)
    } else {
      sentences <- rbind(sentences, sentences)
      sentences$textrank_id <- c(1,2)
      terminology1 <- annt[, c("textrank_id", "lemma")]
      terminology2 <- annt[, c("textrank_id", "lemma")]
      terminology2$textrank_id <- 2
      terminology <- rbind(terminology1, terminology2)
      tr <- textrank_sentences(data = sentences, terminology = terminology)
    }
		s <- summary(tr, n = 1, keep.sentence.order = TRUE)
		s <- gsub('[.]', '', s)
		return(s)
	}, mc.cores = n_cpus)

	tmp <- do.call(rbind, res)
	result <- tmp[cluster]
	return(result)
}
