library(ktplots)
library(ggplot2)
data(kidneyimmune)
seu <- kidneyimmune
sce <- as.SingleCellExperiment(seu)
celltypes <- levels(seu)
gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
bwr <- bower(gmt_file)
bwr <- snn_graph(bwr)
bwr <- find_clusters(bwr)
bwr <- summarize_clusters(bwr)
degs <- Seurat::FindAllMarkers(seu)
degs <- split(degs, degs$cluster)
seu$celltype2 <- as.character(seu$celltype)
sce$celltype2 <- as.character(sce$celltype)

test_that("enrich genesets deg", {
    # so in practice, there should be one DEG table per comparison in a list.
    bwr <- enrich_genesets(degs, bwr, gene_symbol = 'gene', logfoldchanges = 'avg_log2FC',  pvals = 'p_val')
    expect_visible(bwr@scores)
    expect_type(bwr@scores, 'list')
    expect_length(bwr@scores, 6)
    for (ds in celltypes){
        p <- plot_graph(bwr, colorby = ds, mode = 'gsea', gsea.slot = 'NES')
        expect_true(is.ggplot(p))
        expect_length(p$layers, 3)
    }

    bwr <- enrich_genesets(degs, bwr, core = TRUE, gene_symbol = 'gene', logfoldchanges = 'avg_log2FC',  pvals = 'p_val')
    expect_visible(bwr@scores)
    expect_type(bwr@scores, 'list')
    expect_length(bwr@scores, 6)
})
expect_visible(bwr@scores)
test_that("enrich genesets seurat", {
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'AUCell')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'AUCell', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'Seurat')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'Seurat', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'scanpy')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype', mode = 'scanpy', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'AUCell')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'AUCell', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'Seurat')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'Seurat', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'scanpy')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(seu, bwr, groupby = 'celltype2', mode = 'scanpy', core = TRUE)
    expect_visible(bwr@scores)
    
})

test_that("enrich genesets sce", {
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'AUCell')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'AUCell', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'Seurat')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'Seurat', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'scanpy')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype', mode = 'scanpy', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'AUCell')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'AUCell', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'Seurat')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'Seurat', core = TRUE)
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'scanpy')
    expect_visible(bwr@scores)
    expect_is(bwr@scores, 'matrix')
    
    bwr <- enrich_genesets(sce, bwr, groupby = 'celltype2', mode = 'scanpy', core = TRUE)
    expect_visible(bwr@scores)
    
})
