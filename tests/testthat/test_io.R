test_that("read_geneset works", {
	file <- system.file("extdata", "h.all.v7.4.symbols.csv", package = "bowerbird")
	bwr <- bower(file)
	expect_s4_class(bwr, "BOWER")

	df = read.csv(file)
	bwr <- bower(file)
	expect_s4_class(bwr, "BOWER")

	file <- system.file("extdata", "h.all.v7.4.symbols.tsv", package = "bowerbird")
	bwr <- bower(file)
	expect_s4_class(bwr, "BOWER")

	file <- system.file("extdata", "h.all.v7.4.symbols.gmx", package = "bowerbird")
	bwr <- bower(file)
	expect_s4_class(bwr, "BOWER")

	file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
	bwr <- bower(file)
	expect_s4_class(bwr, "BOWER")

	bwr <- snn_graph(bwr)
	bwr <- find_clusters(bwr)
	bwr <- summarize_clusters(bwr)
	bwr <- bower(bwr@coregenes)
	expect_s4_class(bwr, "BOWER")
})
