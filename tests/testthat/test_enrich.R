test_that("enrich genesets deg", {
	library(ktplots)
	data(kidneyimmune)
	celltypes <- levels(kidneyimmune)
	gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
	bwr <- bower(gmt_file)
	bwr <- snn_graph(bwr)
	bwr <- find_clusters(bwr)
	bwr <- summarize_clusters(bwr)

	degs <- Seurat::FindAllMarkers(kidneyimmune)
	degs <- split(degs, degs$cluster) # so in practice, there should be one DEG table per comparison in a list.
	bwr <- tryCatch(enrich_genesets(degs, bwr, gene_symbol = 'gene', logfoldchanges = 'avg_logFC',  pvals = 'p_val'), error = function(e) {
		enrich_genesets(degs, bwr, gene_symbol = 'gene', logfoldchanges = 'avg_log2FC',  pvals = 'p_val')
	})
	
	plot_list <- lapply(celltypes, function(ds){
		g <- plot_graph(bwr, colorby = ds, mode = 'gsea', gsea.slot = 'NES')
	})
	cowplot::plot_grid(plotlist = plot_list, scale = 0.9)
})

test_that("enrich genesets seurat", {
	library(ktplots)
	data(kidneyimmune)
	celltypes <- levels(kidneyimmune)
	gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
	bwr <- bower(gmt_file)
	bwr <- snn_graph(bwr)
	bwr <- find_clusters(bwr)
	bwr <- summarize_clusters(bwr)

	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'AUCell')

	## Seurat::AddModuleScore
	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'Seurat')

	## scanpy.tl.score_genes
	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'scanpy')
})

test_that("enrich genesets sce", {
	library(ktplots)
	data(kidneyimmune)
	celltypes <- levels(kidneyimmune)

	kidneyimmune <- as.SingleCellExperiment(kidneyimmune)

	gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
	bwr <- bower(gmt_file)
	bwr <- snn_graph(bwr)
	bwr <- find_clusters(bwr)
	bwr <- summarize_clusters(bwr)

	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'AUCell')

	## Seurat::AddModuleScore
	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'Seurat')

	## scanpy.tl.score_genes
	bwr <- enrich_genesets(kidneyimmune, bwr, groupby = 'celltype', mode = 'scanpy')
})