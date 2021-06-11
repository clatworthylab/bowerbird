#' @title summarize
#' @include utilities.R
#'

#' Summarize the terms of the cluster using pagerank algorithm
#'
#' @param graph geneset overlap graph.
#' @param cluster vector of cluster labels for each geneset.
#' @param pattern search pattern to remove from the terms. Unless specified, will default to built-in pattern. 
#' @param sep separator used/found in gene set names to be changed to blank spaces. Default value is underscore ('_').
#' @param ncpus number of cores used for parallelizing reconstruction.
#' @param ... passed to igraph::graph_from_adjacency_matrix.
#' @details
#' Given a list of text, it creates a sparse matrix consisting of tf-idf score for tokens from the text. See `https://github.com/saraswatmks/superml/blob/master/R/TfidfVectorizer.R`. A k shortest-nearest neighbor graph is then computed using the overlap of of the terms.
#' @return Returns a matrix of tf-idf score of tokens.
#' @examples
#' \donttest{
#' terms <- summarize_clusters(mat)
#' }
#' @import textrank udpipe dplyr pbmcapply
#' @export
#'

summarize_clusters.BOWER <- function(bower, cluster = NULL, pattern = NULL, sep = NULL, ncpus = NULL){	
  requireNamespace('igraph')
  requireNamespace('parallel')
	if (is.null(pattern)){
		pattern = '^GO_|^KEGG_|^REACTOME_|^HALLMARK_|POSITIVE_|NEGATIVE_|REGULATION_OF|^GOBP_'
	} 
	if (is.null(sep)){
		sep = '_'
	}
	if (!is.null(cluster)){
		cl <- cluster
		igraph::V(bower@graph)$cluster <- cl
	} else {
		cl <- bower@clusters
		igraph::V(bower@graph)$cluster <- cl
	}

	if (is.null(ncpus)){
		n_cpus = parallel::detectCores()
	}

	df <- data.frame(name = igraph::V(bower@graph)$name, cluster = igraph::V(bower@graph)$cluster)
	df$name <- gsub(pattern, '', df$name)
	df_split <- split(df, df$cluster)
	df_split <- pbmclapply(df_split, function(x) {
		x <- x %>% select(name) %>% unlist %>% as.character
		x <- gsub(sep, ' ', x)
		return(x)
	}, mc.cores = n_cpus)

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
	igraph::V(bower@graph)$terms <- tmp[cl]

  # and also only create a label for the centroid node
  data = .graph_to_data(bower@graph)
  data$`_orig_index` = row.names(data)
  data <- split(data, data$terms)
  datax <- lapply(data, function(x) {
    centroid <- .closest_to_centroid(x, 'x', 'y')$`_orig_index`
    return(centroid)})
  idx <- as.numeric(unlist(datax))
  tmp2 <- tmp[cl]
  tmp2[-idx] <- ""
  igraph::V(bower@graph)$labels <- tmp2
  bower@.graph_data <- .graph_to_data(bower@graph)
  
	return(bower)
}

#' @export
summarize_clusters.igraph <- function(graph, cluster = NULL, pattern = NULL, sep = NULL, ncpus = NULL){	
  requireNamespace('igraph')
  requireNamespace('parallel')
	if (is.null(pattern)){
		pattern = '^GO_|^KEGG_|^REACTOME_|^HALLMARK_|POSITIVE_|NEGATIVE_|REGULATION_OF|^GOBP_'
	} 
	if (is.null(sep)){
		sep = '_'
	} 
	if (!is.null(cluster)){
		cl <- cluster
		igraph::V(graph)$cluster <- cluster
	} else {
		cl <- igraph::V(graph)$cluster
		igraph::V(graph)$cluster <- cl
	}

	if (is.null(ncpus)){
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
  result1 <- tmp[cl]
  igraph::V(graph)$terms <- tmp[cl]
  # and also only create a label for the centroid node
  data = .graph_to_data(graph)
  data$`_orig_index` = row.names(data)
  data <- split(data, data$terms)
  datax <- lapply(data, function(x) {
    centroid <- .closest_to_centroid(x, 'x', 'y')$`_orig_index`
    return(centroid)})
  idx <- as.numeric(unlist(datax))
  tmp2 <- tmp[cl]
  tmp2[-idx] <- ""
  labels <- tmp2

	return(list(result, labels))
}
