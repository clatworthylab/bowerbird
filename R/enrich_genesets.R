#' @include utilities.R

#' Performs gene set testing.
#'
#' @title enrich_genesets
#' @description
#' Runs geneset tests and store in BOWER class.
#' @name enrich_genesets
#' @param bower Processed BOWER object.
#' @param list list containing differentially expressed gene testing results in a data frame.
#' @param sce a single cell object in the format of a Seurat object or SingleCellExperiment object.
#' @param groupby Column name in the meta.data/colData of the single cell objects specifying the group to average the enrichment score. If not specified, not cluster average will be calculated.
#' @param core boolean. Whether or not to use the coregenes of genesets slot. Default is FALSE (use genesets).
#' @param standardize whether or not to standardize mean enrichment values to 0 to 1. Only used if groupby is not NULL.
#' @param mode choice of enrichment test to perform.
#' @param gene_symbol column name for gene_symbol for gsea.
#' @param logfoldchanges column name for logfoldchanges for gsea.
#' @param pvals column name for pvals for gsea.
#' @param remove_mito_ribo boolean. whether or not to remove mitochondial and ribosomal genes from consideration for gsea. Default is TRUE.
#' @param nperm number of permuation iterations for gsea. Default is 10000.
#' @param minSize minimum geneset size for gsea. Default is 0.
#' @param maxSize maximum geneset size for gsea. Default is 1000.
#' @param sce_assay name of assay in SingleCellExperiment object. 
#' @param seurat_assay name of assay in Seurat object.
#' @param ncpus number of cores used for parallelizing geneset testing.
#' @param aucMaxRank_pct percentage to use for aucMaxRank in AUCell::AUCell_calcAUC.
#' @param ... passed to fgsea::fgsea, AUCell::AUCell_buildRankings or Seurat::AddModuleScore
#' @details
#' 
#' @return Returns a dataframe of average gene set scores.
#' @examples
#' \donttest{
#' bwr <- enrich_genesets(list, bwr)
#' bwr <- enrich_genesets(seurat, bwr)
#' bwr <- enrich_genesets(sce, bwr)
#' }
#' @import fgsea dplyr
#' @export
#' 

enrich_genesets.list <- function(list, bower, core = FALSE, gene_symbol = 'X1', logfoldchanges = 'logfoldchanges', pvals = 'pvals', remove_mito_ribo = TRUE, minSize = 0, maxSize = 1000, ...){
  scex <- lapply(list, .makeRankGeneList, gene_symbol = gene_symbol, logfoldchanges = logfoldchanges, pvals = pvals, remove_mito_ribo = remove_mito_ribo)
  if (core) {
    res <- lapply(scex, function(x) fgsea(pathways = bower@coregenes, stats = x, minSize = minSize, maxSize = maxSize, ...))
  } else {
    res <- lapply(scex, function(x) fgsea(pathways = bower@genesets, stats = x, minSize = minSize, maxSize = maxSize, ...))
  }
  res <- lapply(res, as.data.frame)
  
  NES <- .retrieve_gsea(res, 'NES')
  ES <- .retrieve_gsea(res, 'ES')
  pval <- .retrieve_gsea(res, 'pval')
  padj <- .retrieve_gsea(res, 'padj')
  leadingEdge <- .retrieve_gsea(res, 'leadingEdge')
  
  df <- data.frame(row.names = bwr@.graph_data$name, name = bwr@.graph_data$name)

  NES <- df %>% left_join(NES, by = 'name')
  ES <- df %>% left_join(ES, by = 'name')
  pval <- df %>% left_join(pval, by = 'name')
  padj <- df %>% left_join(padj, by = 'name')
  leadingEdge <- df %>% left_join(leadingEdge, by = 'name')

  NES <- NES[,-1]
  ES <- ES[,-1]
  pval <- pval[,-1]
  padj <- padj[,-1]
  leadingEdge <- leadingEdge[,-1]

  bower@scores <- list(NES = NES, ES = ES, pval = pval, padj = padj, leadingEdge = leadingEdge, full_fgsea = res)
  return(bower)
}

#' @rdname enrich_genesets
#' @export
enrich_genesets.Seurat <- function(sce, bower, groupby = NULL, core = FALSE, standardize = TRUE, mode = c('AUCell', 'Seurat', 'scanpy'), sce_assay = 'logcounts', seurat_assay = 'RNA', ncpus = NULL, aucMaxRank_pct = 5, ...) {
  if (mode[1] == 'AUCell') {
    requireNamespace(c('AUCell', 'SingleCellExperiment', 'SummarizedExperiment'))
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (is.null(ncpus)){
      n_cpus = parallel::detectCores()
    }
    cells_rankings <- AUCell::AUCell_buildRankings(SummarizedExperiment::assays(scex)[[sce_assay]], nCores = n_cpus, plotStats = FALSE, verbose = FALSE, ...)
    if (core){
      cells_AUC <- AUCell::AUCell_calcAUC(bower@coregenes, cells_rankings, aucMaxRank = nrow(cells_rankings) * (aucMaxRank_pct/100), verbose = FALSE)
    } else {
      cells_AUC <- AUCell::AUCell_calcAUC(bower@genesets, cells_rankings, aucMaxRank = nrow(cells_rankings) * (aucMaxRank_pct/100), verbose = FALSE)
    }
    for (i in 1:nrow(cells_AUC@assays@data$AUC)) {
      geneSetsAUC <- cells_AUC@assays@data$AUC[i, ]
      SummarizedExperiment::colData(scex)[, row.names(cells_AUC@assays@data$AUC)[i]] <- geneSetsAUC
    }
    AUC_mat <- t(data.frame(cells_AUC@assays@data$AUC))  
    
    if (!is.null(groupby)){
      check <- groupby %in% colnames(SummarizedExperiment::colData(scex))
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(SummarizedExperiment::colData(scex)[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(SummarizedExperiment::colData(scex)[,groupby]), function(x){
          rm <- rowMeans(t(AUC_mat)[, SummarizedExperiment::colData(scex)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(SummarizedExperiment::colData(scex)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(SummarizedExperiment::colData(scex)[,groupby]), function(x){
          rm <- rowMeans(t(AUC_mat)[, SummarizedExperiment::colData(scex)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(SummarizedExperiment::colData(scex)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- AUC_mat
    }
  } else if (mode[1] == 'Seurat') {
    requireNamespace('Seurat')
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (core){
      scex <- Seurat::AddModuleScore(scex, features = bower@coregenes, assay = seurat_assay, name = 'X__bower__', ...)
      gs_names <- names(bower@coregenes)
    } else {
      scex <- Seurat::AddModuleScore(scex, features = bower@genesets, assay = seurat_assay, name = 'X__bower__', ...)
      gs_names <- names(bower@genesets)
    }
    idx <- grep('X__bower__', colnames(scex@meta.data))
    module_mat <- scex@meta.data[, idx]
    colnames(module_mat) <- gs_names  
    if (!is.null(groupby)){
      check <- groupby %in% colnames(scex@meta.data)
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(scex@meta.data[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(scex@meta.data[,groupby]), function(x){
        rm <- rowMeans(t(module_mat)[, scex@meta.data[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(scex@meta.data[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(scex@meta.data[,groupby]), function(x){
          rm <- rowMeans(t(module_mat)[, scex@meta.data[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(scex@meta.data[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- module_mat
    }
  } else if (mode[1] == 'scanpy') {
    requireNamespace('reticulate')
    sc <- tryCatch(reticulate::import('scanpy'), error = function(e) stop('Please install scanpy to use this mode.'))
    sce1 <- .typecheck(sce, mode = 'SingleCellExperiment', sce_assay = sce_assay, seurat_assay = seurat_assay)
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (core) {
      for (i in seq_along(bower@coregenes)) {
        sc$tl$score_genes(scex, unlist(bower@coregenes[i]), score_name=names(bower@coregenes)[i])
      }
      mat <- data.frame(scex$obs[names(bower@coregenes)])
    } else {
      for (i in seq_along(bower@genesets)) {
        sc$tl$score_genes(scex, unlist(bower@genesets[i]), score_name=names(bower@genesets)[i])
      }
      mat <- data.frame(scex$obs[names(bower@genesets)])
    }      
    if (!is.null(groupby)){
      check <- groupby %in% colnames(SummarizedExperiment::colData(sce1))
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(SummarizedExperiment::colData(sce1)[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(SummarizedExperiment::colData(sce1)[,groupby]), function(x){
          rm <- rowMeans(t(mat)[, SummarizedExperiment::colData(sce1)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(SummarizedExperiment::colData(sce1)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(SummarizedExperiment::colData(sce1)[,groupby]), function(x){
          rm <- rowMeans(t(mat)[, SummarizedExperiment::colData(sce1)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(SummarizedExperiment::colData(sce1)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- mat
    }
  }
  if (standardize && !is.null(groupby)){    
    scores <- t(apply(as.data.frame(scores), 1, range01))
  }
  if (!core) {
    if (!is.null(groupby)){
      requireNamespace('igraph')
      mcn <- match(igraph::V(bower@graph)$name, rownames(scores))
      bower@scores <- scores[mcn,]
    } else {
      mcn <- match(igraph::V(bower@graph)$name, colnames(scores))
      bower@scores <- t(scores[,mcn])
    }    
  } else {
    if (!is.null(groupby)){
      bower@scores <- scores
    } else {
      bower@scores <- t(scores)
    }
  }
  return(bower)
}

#' @rdname enrich_genesets
#' @export
enrich_genesets.SingleCellExperiment <- function(sce, bower, groupby = NULL, core = FALSE, standardize = TRUE, mode = c('AUCell', 'Seurat', 'scanpy'), sce_assay = 'logcounts', seurat_assay = 'RNA', ncpus = NULL, aucMaxRank_pct = 5, ...) {
  if (mode[1] == 'AUCell') {
    requireNamespace(c('AUCell', 'SingleCellExperiment', 'SummarizedExperiment'))
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (is.null(ncpus)){
      n_cpus = parallel::detectCores()
    }
    cells_rankings <- AUCell::AUCell_buildRankings(SummarizedExperiment::assays(scex)[[sce_assay]], nCores = n_cpus, plotStats = FALSE, verbose = FALSE, ...)
    if (core){
      cells_AUC <- AUCell::AUCell_calcAUC(bower@coregenes, cells_rankings, aucMaxRank = nrow(cells_rankings) * (aucMaxRank_pct/100), verbose = FALSE)
    } else {
      cells_AUC <- AUCell::AUCell_calcAUC(bower@genesets, cells_rankings, aucMaxRank = nrow(cells_rankings) * (aucMaxRank_pct/100), verbose = FALSE)
    }
    for (i in 1:nrow(cells_AUC@assays@data$AUC)) {
      geneSetsAUC <- cells_AUC@assays@data$AUC[i, ]
      SummarizedExperiment::colData(scex)[, row.names(cells_AUC@assays@data$AUC)[i]] <- geneSetsAUC
    }
    AUC_mat <- t(data.frame(cells_AUC@assays@data$AUC))  
    
    if (!is.null(groupby)){
      check <- groupby %in% colnames(SummarizedExperiment::colData(scex))
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(SummarizedExperiment::colData(scex)[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(SummarizedExperiment::colData(scex)[,groupby]), function(x){
          rm <- rowMeans(t(AUC_mat)[, SummarizedExperiment::colData(scex)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(SummarizedExperiment::colData(scex)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(SummarizedExperiment::colData(scex)[,groupby]), function(x){
          rm <- rowMeans(t(AUC_mat)[, SummarizedExperiment::colData(scex)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(SummarizedExperiment::colData(scex)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- AUC_mat
    }
  } else if (mode[1] == 'Seurat') {
    requireNamespace('Seurat')
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (core){
      scex <- Seurat::AddModuleScore(scex, features = bower@coregenes, assay = seurat_assay, name = 'X__bower__', ...)
      gs_names <- names(bower@coregenes)
    } else {
      scex <- Seurat::AddModuleScore(scex, features = bower@genesets, assay = seurat_assay, name = 'X__bower__', ...)
      gs_names <- names(bower@genesets)
    }
    idx <- grep('X__bower__', colnames(scex@meta.data))
    module_mat <- scex@meta.data[, idx]
    colnames(module_mat) <- gs_names  
    if (!is.null(groupby)){
      check <- groupby %in% colnames(scex@meta.data)
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(scex@meta.data[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(scex@meta.data[,groupby]), function(x){
        rm <- rowMeans(t(module_mat)[, scex@meta.data[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(scex@meta.data[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(scex@meta.data[,groupby]), function(x){
          rm <- rowMeans(t(module_mat)[, scex@meta.data[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(scex@meta.data[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- module_mat
    }
  } else if (mode[1] == 'scanpy') {
    requireNamespace('reticulate')
    sc <- tryCatch(reticulate::import('scanpy'), error = function(e) stop('Please install scanpy to use this mode.'))
    sce1 <- .typecheck(sce, mode = 'SingleCellExperiment', sce_assay = sce_assay, seurat_assay = seurat_assay)
    scex <- .typecheck(sce, mode = mode, sce_assay = sce_assay, seurat_assay = seurat_assay)
    if (core) {
      for (i in seq_along(bower@coregenes)) {
        sc$tl$score_genes(scex, unlist(bower@coregenes[i]), score_name=names(bower@coregenes)[i])
      }
      mat <- data.frame(scex$obs[names(bower@coregenes)])
    } else {
      for (i in seq_along(bower@genesets)) {
        sc$tl$score_genes(scex, unlist(bower@genesets[i]), score_name=names(bower@genesets)[i])
      }
      mat <- data.frame(scex$obs[names(bower@genesets)])
    }      
    if (!is.null(groupby)){
      check <- groupby %in% colnames(SummarizedExperiment::colData(sce1))
      if (!check){
        stop(paste0(groupby, ' not found in single cell data object. Please try again.'))
      }
      checkfactor <- class(SummarizedExperiment::colData(sce1)[,groupby]) == 'factor'
      if (checkfactor){
        scores <- do.call(cbind, lapply(levels(SummarizedExperiment::colData(sce1)[,groupby]), function(x){
          rm <- rowMeans(t(mat)[, SummarizedExperiment::colData(sce1)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- levels(SummarizedExperiment::colData(sce1)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      } else {
        scores <- do.call(cbind, lapply(unique(SummarizedExperiment::colData(sce1)[,groupby]), function(x){
          rm <- rowMeans(t(mat)[, SummarizedExperiment::colData(sce1)[,groupby] %in% x])
          return(rm)
          })
        )
        colnames(scores) <- unique(SummarizedExperiment::colData(sce1)[,groupby])
        if (core) {
          rownames(scores) <- names(bower@coregenes)
        } else {
          rownames(scores) <- names(bower@genesets)
        }
      }
    } else {
      scores <- mat
    }
  }
  if (standardize && !is.null(groupby)){    
    scores <- t(apply(as.data.frame(scores), 1, range01))
  }
  if (!core) {
    if (!is.null(groupby)){
      requireNamespace('igraph')
      mcn <- match(igraph::V(bower@graph)$name, rownames(scores))
      bower@scores <- scores[mcn,]
    } else {
      mcn <- match(igraph::V(bower@graph)$name, colnames(scores))
      bower@scores <- t(scores[,mcn])
    }    
  } else {
    if (!is.null(groupby)){
      bower@scores <- scores
    } else {
      bower@scores <- t(scores)
    }
  }
  return(bower)
}
