#' @include utilities.R
#' 

#' Reads various geneset file formats and returns a list of vectors containing gene names
#'
#' @import dplyr
#' 
read_geneset.gmt <- function(x, min_size) {

    n <- max(count.fields(x, sep='\t'), na.rm=TRUE)
    x <- readLines(x)

    .splitvar <- function(x, sep, n) {
        var <- unlist(strsplit(x, split='\t'))
        length(var) = n
        return(var)
    }

    x <- do.call(cbind, lapply(x, .splitvar, sep='\t', n=n))
    x <- apply(x, 1, paste, collapse='\t')
    out <- as.list(read.csv(text=x, sep='\t', check.names = FALSE))
    out <- lapply(out, function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- lapply(out, function(x) x[-1])
    out <- .filter_geneset(out, min_size)
    return(out)
}


read_geneset.gmx <- function(x, min_size) {
    out <- as.list(read.csv(x, sep = '\t', check.names = FALSE, na=c("", NA)))
    out <- lapply(out, function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- lapply(out, function(x) x[-1])
    out <- .filter_geneset(out, min_size)
    return(out)
}


read_geneset.csv <- function(x, min_size) {
    out <- as.list(read.csv(x, check.names = FALSE, na=c("", NA)))
    out <- lapply(out, function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- .filter_geneset(out, min_size)
    return(out)
}


read_geneset.tsv <- function(x, min_size) {
    out <- as.list(read.delim(x, sep = '\t', check.names = FALSE, na=c("", NA)))
    out <- lapply(out, function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- .filter_geneset(out, min_size)
    return(out)
}

read_geneset.data.frame <- function(x, min_size) {
    out <- lapply(as.list(x), function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- .filter_geneset(out, min_size)
    return(out)
}

read_geneset.list <- function(x, min_size) {
    if (is.null(names(x))){
        stop('Please provide names to input list')
    }
    out <- lapply(x, function(x) x %>% unlist %>% as.character)
    out <- lapply(out, function(x) x[!is.na(x)])
    out <- .filter_geneset(out, min_size)
    return(out)
}
