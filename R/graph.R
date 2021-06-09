#' @title snn_graph
#' @include utilities.R

#' Creates a k-nearest-neighbor graph from a tf-idf matrix.
#'
#' @name snn_graph
#' @param gs genesets in list.
#' @param max_features use top features sorted by count to be used in bag of words matrix. The default value is set to 100.
#' @param remove_stopwords a list of stopwords to use, by default it uses its inbuilt list of standard stopwords. The default value is FALSE.
#' @param k the maximum number of nearest neighbors to search. The default value is set to 5.
#' @param ... passed to superml::TfIdfVectorizer
#' @details
#' Given a list of text, it creates a sparse matrix consisting of tf-idf score for tokens from the text. See `https://github.com/saraswatmks/superml/blob/master/R/TfidfVectorizer.R`. A k shortest-nearest neighbor graph is then computed using the overlap of of the terms.
#' @return Returns a matrix of tf-idf score of tokens.
#' @examples
#' \donttest{
#' gs <- scFA::read_geneset('geneset.gmt')
#' mat <- snn_graph(gs)
#' }
#' @import superml
#' @export

snn_graph.BOWER <- function(bower, max_features = 100, remove_stopwords = FALSE, k = 5, ...){
    tfv <- TfIdfVectorizer$new(max_features = max_features, remove_stopwords = remove_stopwords, ...)
    tfv_mat <- tfv$fit_transform(bower@genesets)
    rownames(tfv_mat) <- names(bower@genesets)

    requireNamespace('FNN')
    # make a SNN graph
    knn_idx <- FNN::get.knn(tfv_mat, k = k)$nn.index
    rownames(knn_idx) <- names(bower@genesets)
    snn <-  matrix(0, ncol = length(names(bower@genesets)), nrow = length(names(bower@genesets)))
    colnames(snn) <- rownames(snn) <- names(bower@genesets)

    for(i in rownames(knn_idx)){
        i_vec <- knn_idx[i, ]
        for(j in i_vec){
            j <- colnames(snn)[j]
            j_vec <- knn_idx[j, ]
            snn[i, j] <- length(intersect(i_vec, j_vec))/length(union(i_vec, j_vec))
        }
    }
    bower@graph <- .create_graph(snn)    
    return(bower)
}

#' @name snn_graph
#' @export
snn_graph.list <- function(gs, max_features = 100, remove_stopwords = FALSE, k = 5, ...){
    tfv <- TfIdfVectorizer$new(max_features = max_features, remove_stopwords = remove_stopwords, ...)
    tfv_mat <- tfv$fit_transform(gs)
    rownames(tfv_mat) <- names(gs)

    requireNamespace('FNN')
    # make a SNN graph
    knn_idx <- FNN::get.knn(tfv_mat, k = k)$nn.index
    rownames(knn_idx) <- names(gs)
    snn <-  matrix(0, ncol = length(names(gs)), nrow = length(names(gs)))
    colnames(snn) <- rownames(snn) <- names(gs)

    for(i in rownames(knn_idx)){
        i_vec <- knn_idx[i, ]
        for(j in i_vec){
            j <- colnames(snn)[j]
            j_vec <- knn_idx[j, ]
            snn[i, j] <- length(intersect(i_vec, j_vec))/length(union(i_vec, j_vec))
        }
    }
    return(.create_graph(snn))
}

