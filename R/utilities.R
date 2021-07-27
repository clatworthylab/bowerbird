#' @description
#' Miscellaneous utility functions.
#'
#########################
# Miscellaneous functions
#########################
#' @param x vector.
#'

.is_blank <- function(x) x == ""


#' @param x named list.
#'

.sanitize <- function(x) {
  names(x) <- make.unique(names(x), sep = '_')
  return(x)
}


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
  dat <- g$data
  for (i in seq_along(colnames(dat))){
    if (colnames(dat)[i] %in% names(igraph::vertex_attr(gr)))
      if (colnames(dat)[i] == 'cluster') {
        dat[,i] <- factor(dat[,i])
      } else if (colnames(dat)[i] == 'labels') {
        dat[,i][dat[,i] == ""] <- NA
        dat[,i] <- factor(dat[,i])
      } else if (!is.numeric(dat[,i])){
        dat[,i] <- factor(dat[,i])
      }
  }
  return(dat)
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


#' @param sce Single cell object. Accepts either a Seurat object of SingleCellExperiment object.
#' @param mode mode of conversion. Accepts one of AUCell, SingleCellExperiment, Seurat, scanpy.
#' @param sce_assay name of assay in SingleCellExperiment object.
#' @param seurat_assay name of assay in Seurat object.
#' @return appropriate single cell object for enrichment analysis
#'
.typecheck <- function(sce, mode, sce_assay = 'logcounts', seurat_assay = 'RNA'){
  cls <- class(sce)
  requireNamespace(c('Seurat', 'SummarizedExperiment'))
  if (mode[1] == 'AUCell' | mode[1] == 'SingleCellExperiment'){
    if (cls == 'Seurat'){
      scex <- Seurat::as.SingleCellExperiment(sce)
    } else if (cls == 'SingleCellExperiment') {
      scex <- sce
    }
  } else if (mode[1] == 'Seurat'){
    if (cls == 'SingleCellExperiment'){
      scex <- Seurat::as.Seurat(sce)
    } else if (cls == 'Seurat') {
      scex <- sce
    }
  } else if (mode[1] == 'scanpy'){
    requireNamespace(c('Matrix', 'reticulate'))
    sc <- tryCatch(reticulate::import('scanpy'), error = function(e) stop('Please install scanpy to use this mode.'))
    if (cls == 'SingleCellExperiment'){
      if (sce_assay %in% names(SummarizedExperiment::assays(sce))){
        mat <- SummarizedExperiment::assays(sce)[[sce_assay]]
      } else {
        stop(paste0(sce_assay, ' not found in SingleCellExperiment object.'))
      }
    } else if (cls == 'Seurat'){
      if (seurat_assay %in% names(sce)){
        mat <- sce[[seurat_assay]]@data
      } else {
        stop(paste0(seurat_assay, ' not found in Seurat object.'))
      }
    }
    if (!"dgCMatrix" %in% class(mat)){
      mat <- Matrix::Matrix(mat, sparse = TRUE)
    }
    obs <- data.frame(row.names = colnames(mat), barcode = colnames(mat))
    var <- data.frame(row.names = rownames(mat), gene_ids = rownames(mat))
    scex <- sc$AnnData(X = Matrix::t(mat), obs = obs, var = var)
    scex$raw <- scex
  }
  return(scex)
}


#' @param deg deg table. must at least have a column for gene symbol, log fold changes and pvalues.
#' @param gene_symbol gene_symbol
#' @param logfoldchanges logfoldchanges
#' @param pvals pvals
#' @param remove_mito_ribo boolean. whether or not to remove mitochondial and ribosomal genes from consideration. Default is TRUE.
#' @return sorted ranked gene list.
#'
.makeRankGeneList <- function(deg, gene_symbol, logfoldchanges, pvals, remove_mito_ribo = TRUE){
  if (remove_mito_ribo) {
    y <- grepl('^RPS|^RPL|^MRPL|^MRPS|^MT-|^Rps|^Rpl|^Mrpl|^Mrps|^mt-', deg[,gene_symbol])
    deg <- deg[!y, ]
  }
  deg <- deg[, c(gene_symbol, logfoldchanges, pvals)]
  deg$nedegog10pval <- -log10(deg[, pvals])
  rank <- unlist(deg$nedegog10pval * sign(deg[,logfoldchanges]))
  rank[rank == Inf] = 300
  rank[rank == -Inf] = -300
  names(rank) <- deg[,gene_symbol]
  rank <- rev(sort(rank))
  return(rank)
}


#' @param x numerical vector
#'
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


#' @param res fgsea result data.frame.
#' @param query column name in fgsea results.
#'
.retrieve_gsea <- function(res, query){
  tmp <- lapply(res, function(x) {
    q <- tryCatch(data.frame(pathway = x[,'pathway'], query = x[,query]), error = function(e) data.frame(pathway = x[,'pathway'], query = do.call(rbind, lapply(x[,query], paste0, collapse = '|'))))
    return(q)
  })
  for (i in seq_along(tmp)){
    colnames(tmp[[i]]) <- c('pathway', names(tmp)[i])
  }
  tmp <- Reduce(function(...) merge(..., by='pathway', all.x=TRUE), tmp)
  row.names(tmp) <- tmp$pathway
  tmp <- tmp[, -1]
  tmp$name <- row.names(tmp)
  return(tmp)
}


#' @param filename, string. The filename of the file in the package cache.
#' @param mustWork, logical. Whether an error should be created if the file does not exist. If mustWork=FALSE and the file does not exist, the empty string is returned.
#' @import pkgfilecache
#' @return Access a single file from the package cache by its file name.
#'
#' @export
get_optional_data_filepath <- function(filename, mustWork=TRUE) {
  pkg_info = get_pkg_info("bowerbird");
  return(get_filepath(pkg_info, filename, mustWork=mustWork));
}


#' @import udpipe pkgfilecache rappdirs
#' @return Access a single file from the package cache by its file name.
#' @export
check_udpipemodel <- function(){
  pkg_info = get_pkg_info("bowerbird");
  local_filenames = 'english-ewt-ud-2.5-191206.udpipe' # need to modify this if it gets updated?
  files_exist = are_files_available(pkg_info, local_filenames)
  if (files_exist){
    tagger <- data.frame(
                         language = 'english-ewt',
                         file_model = get_optional_data_filepath(local_filenames),
                         url = 'https://raw.githubusercontent.com/jwijffels/udpipe.models.ud.2.5/master/inst/udpipe-ud-2.5-191206/english-ewt-ud-2.5-191206.udpipe',
                         download_failed = FALSE,
                         download_message = "OK")
  } else {
    tagger <- udpipe_download_model("english", model_dir = file.path(user_data_dir(), 'bowerbird'))
  }
  return(tagger)
}