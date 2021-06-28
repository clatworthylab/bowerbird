#########################
# Miscellaneous functions
#########################

#' @param x vector.
#' 

.is_blank <- function(x) x == ""


#' @param gs list of vectors.
#' @param min_size minimum size of geneset otherwise filter out.
#' 

.filter_geneset <- function(gs, min_size) {
    gs <- lapply(gs, function(x) x[!.is_blank(x)])
    if (length(which(is.na(gs))) > 0){
        gs <- gs[!is.na(gs)]
    }
    if (length(which(lapply(gs, length) <= min_size)) > 0) {
        gs <- gs[-which(lapply(gs, length) <= min_size)]
    }
    return(gs)
}


#' @param ... does nothing.
#' 

.emptyBOWER <- function(...) {
    out <- new("BOWER",
               genesets=list(),
               graph=NULL,
               clusters=NULL,
               .graph_data=NULL)
    out
}


#' @param snn KNN graph.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See ?igraph::graph_from_adjacency_matrix.
#' @param weighted This argument specifies whether to create a weighted graph from an adjacency matrix. If it is NULL then an unweighted graph is created and the elements of the adjacency matrix gives the number of edges between the vertices. If it is a character constant then for every non-zero matrix entry an edge is created and the value of the entry is added as an edge attribute named by the weighted argument. If it is TRUE then a weighted graph is created and the name of the edge attribute will be weight. See ?igraph::graph_from_adjacency_matrix.
#' @param ... passed to igraph::graph_from_adjacency_matrix.
#'

.create_graph <- function(snn, mode = 'undirected', weighted = TRUE, ...) {
    requireNamespace('igraph')
    gr <- igraph::graph_from_adjacency_matrix(snn, mode= "undirected", weighted = TRUE, ...)
    gr <- igraph::simplify(gr)
    return(gr)
}

#' @param gr igraph
#' @param ... passed to ggraph::ggraph.
#' 
.graph_to_data <- function(gr, ...) {
    requireNamespace('ggraph')
    if (length(list(...)) == 0){
        g <- ggraph::ggraph(gr, 'igraph', algorithm = 'fr')
    } else {
        g <- ggraph::ggraph(gr, ...)
    }

    return(g$data)
}

#' @param gr data matrix or data.frame that contains at least two columns.
#' @param x_column name of x column/axis.
#' @param y_column name of y column/axis.
#' @param n (numeric) number of points that are close to the centroid to be detected. Default = 1.
#' @param id_column (character or numeric) name or numeric index of the column in data containing identifiers of one or distinct sets of points. If, NULL, the default, only one set is assumed.
#' 

.closest_to_centroid <- function(data, x_column, y_column, n = 1, id_column = NULL) {
  # adapted from https://rdrr.io/github/claununez/biosurvey/src/R/selection_helpers.R
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
