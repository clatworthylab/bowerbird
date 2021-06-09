#' @title summarize_clusters
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
#' terms <- summarize_clusters(mat)
#' @import textrank udpipe dplyr pbmcapply
#' @export
#'

summarize_clusters.BOWER <- function(bower, cluster = NULL, pattern = NULL, sep = NULL, ncpus = NULL){	
  requireNamespace('igraph')
  requireNamespace('parallel')
	if (is.null(pattern)){
		pattern = '^GO_|^KEGG_|^REACTOME_|^HALLMARK_|POSITIVE_|NEGATIVE_|REGULATION_OF|^GOBP_'
	} 
	if (is.null(pattern)){
		sep
	}
	if (!is.null(cluster)){
		igraph::V(bower@graph)$cluster <- cluster		
	} else {
		igraph::V(bower@graph)$cluster <- bower@clusters
	}

	if (!is.null(ncpus)){
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
	
	igraph::V(bower@graph)$labels <- tmp[cluster]

	return(bower)
}

#' @export
summarize_clusters.igraph <- function(graph, cluster = NULL, pattern = NULL, sep = NULL, ncpus = NULL){	
  requireNamespace('igraph')
  requireNamespace('parallel')
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


#' @export
closest_to_centroid.BOWER <- function(bower, x_column, y_column, n = 1, id_column = NULL) {
  requireNamespace('stats')
  # Initial tests
  if (missing(bower@graph)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' is not defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' is not defined.")
  }
  coln <- colnames(data)
  if (!x_column %in% coln) {
    stop(x_column, " is not one o the columns in 'data'.")
  }
  if (!y_column %in% coln) {
    stop(y_column, " is not one o the columns in 'data'.")
  }  

  if (!is.numeric(data[,x_column])){
    data[,x_column] <- as.numeric(data[,x_column])
  }
  if (!is.numeric(data[,y_column])){
    data[,y_column] <- as.numeric(data[,y_column])
  }
  # Detection of sets of points
  if(!is.null(id_column)) {
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  } else {
    data <- cbind(data, id_column = 1)
    id_column <- "id_column"
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  }

  # Detection of points closest to centroid
  ucent <- lapply(bs, function(x) {
    gblock <- data[bda == x, ]
    ## Process for more than 2 points
    if (nrow(gblock) > 2) {
      cent <- apply(gblock[, c(x_column, y_column)], 2, mean)
      level <- 0.01
      if (nrow(gblock) > 20) {
        ## Process for more than 20 points
        covm <- stats::cov(gblock[, c(x_column, y_column)])
        ndim <- length(cent); sigma_i <- solve(covm) / stats::qchisq(level,
                                                                     df = ndim)
        stds <- 1 / sqrt(eigen(sigma_i)$values)
        hl <- cent + stds; ll <- cent - stds
        c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
          gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
        con <- sum(c1)

        if (con <= 2) {
          ## Loop for measuring distances among distinct amount of points
          while (con == 0) {
            if (level > 0.98) {
                ## Distance in E space
                ds <- stats::mahalanobis(x = gblock[, c(x_column, y_column)], center = cent, cov = covm, tol = 0.0000009)
                break()
            }
            ## Dtecting whether closest point was detected or not
            level <- level + 0.01
            sigma_i <- solve(covm) / stats::qchisq(level, df = ndim)
            stds <- 1 / sqrt(eigen(sigma_i)$values)
            hl <- cent + stds
            ll <- cent - stds
            c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
              gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
            con <- sum(c1)
            if (con > 0) {
              break()
            }
          }
        }
      } else {
        c1 <- rep(TRUE, nrow(gblock))
      }

      # Returning results according to condition
      if (level > 0.98) {
        return(gblock[which(ds == sort(ds)[1:n])[1:n], ])
      }
        ds <- as.matrix(stats::dist(rbind(cent, gblock[c1, c(x_column, y_column)])))[-1, 1]
        return(gblock[which(ds %in% sort(ds)[1:n])[1:n], ])
    } else {
      return(gblock[1, ])
    }
  })
  return(do.call(rbind, ucent))
}

#' @export
closest_to_centroid.BOWER <- function(data, x_column, y_column, n = 1, id_column = NULL) {
  requireNamespace('stats')
  # Initial tests
  if (missing(data)) {
    stop("Argument 'data' is not defined.")
  }
  if (missing(x_column)) {
    stop("Argument 'x_column' is not defined.")
  }
  if (missing(y_column)) {
    stop("Argument 'y_column' is not defined.")
  }
  coln <- colnames(data)
  if (!x_column %in% coln) {
    stop(x_column, " is not one o the columns in 'data'.")
  }
  if (!y_column %in% coln) {
    stop(y_column, " is not one o the columns in 'data'.")
  }  

  if (!is.numeric(data[,x_column])){
    data[,x_column] <- as.numeric(data[,x_column])
  }
  if (!is.numeric(data[,y_column])){
    data[,y_column] <- as.numeric(data[,y_column])
  }
  # Detection of sets of points
  if(!is.null(id_column)) {
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  } else {
    data <- cbind(data, id_column = 1)
    id_column <- "id_column"
    bda <- data[, id_column]
    bs <- sort(unique(bda))
  }

  # Detection of points closest to centroid
  ucent <- lapply(bs, function(x) {
    gblock <- data[bda == x, ]
    ## Process for more than 2 points
    if (nrow(gblock) > 2) {
      cent <- apply(gblock[, c(x_column, y_column)], 2, mean)
      level <- 0.01
      if (nrow(gblock) > 20) {
        ## Process for more than 20 points
        covm <- stats::cov(gblock[, c(x_column, y_column)])
        ndim <- length(cent); sigma_i <- solve(covm) / stats::qchisq(level,
                                                                     df = ndim)
        stds <- 1 / sqrt(eigen(sigma_i)$values)
        hl <- cent + stds; ll <- cent - stds
        c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
          gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
        con <- sum(c1)

        if (con <= 2) {
          ## Loop for measuring distances among distinct amount of points
          while (con == 0) {
            if (level > 0.98) {
                ## Distance in E space
                ds <- stats::mahalanobis(x = gblock[, c(x_column, y_column)], center = cent, cov = covm, tol = 0.0000009)
                break()
            }
            ## Dtecting whether closest point was detected or not
            level <- level + 0.01
            sigma_i <- solve(covm) / stats::qchisq(level, df = ndim)
            stds <- 1 / sqrt(eigen(sigma_i)$values)
            hl <- cent + stds
            ll <- cent - stds
            c1 <- gblock[, x_column] >= ll[1] & gblock[, x_column] <= hl[1] &
              gblock[, y_column] >= ll[2] & gblock[, y_column] <= hl[2]
            con <- sum(c1)
            if (con > 0) {
              break()
            }
          }
        }
      } else {
        c1 <- rep(TRUE, nrow(gblock))
      }

      # Returning results according to condition
      if (level > 0.98) {
        return(gblock[which(ds == sort(ds)[1:n])[1:n], ])
      }
        ds <- as.matrix(stats::dist(rbind(cent, gblock[c1, c(x_column, y_column)])))[-1, 1]
        return(gblock[which(ds %in% sort(ds)[1:n])[1:n], ])
    } else {
      return(gblock[1, ])
    }
  })
  return(do.call(rbind, ucent))
}