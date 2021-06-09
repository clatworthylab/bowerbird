# bowerbird
A single-cell package for Functional Annotation + Gene Module Summarization. 
 
*Bowerbirds are famous for the elaborate and sometimes whimsical structures that males build to court females.*

I think this is vaguely similar to my attempts to visualise pathway analyses and their associated with celltype identities/states...

<a href="https://blog.nature.org/science/2021/01/04/bowerbirds-meet-the-bird-worlds-kleptomaniac-love-architects/"><img src="https://blog.nature.org/science/files/2020/11/32487196918_8dd537c82a_k.jpg" alt="bower" width="400"/></a> <img src="img/logo.png" alt="bower" width="200"/> 

Picture credits: 
A satin bowerbird’s bower, decorated with blue objects. © Stefan Marks / [Flickr](https://www.flickr.com/photos/stefan_marks/32487196918/)

## Installation
```R
devtools::install_github('clatworthylab/bowerbird', auth_token = "insert_your_personal_github_access_token")
```

## Quick Usage
```R
library(bowerbird)
# download from msigdb website
bwr <- bower('h.all.v7.4.symbols.gmt') # this performs a read_geneset step internally, which accepts .gmt, .gmx, .csv, .tsv, .txt, or R objects as list or data.frame format.
bwr <- snn_graph(bwr)
bwr <- find_clusters(bwr)
bwr <- summarize_clusters(bwr)
bwr
# Number of genesets:  50 
# Geneset KNN Graph: 
# IGRAPH 5c4d1b0 UNW- 50 124 -- 
# + attr: name (v/c), cluster (v/n), terms (v/c), labels (v/c), weight (e/n)
# + edges from 5c4d1b0 (vertex names):
#  [1] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_HYPOXIA                          
#  [2] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_TGF_BETA_SIGNALING               
#  [3] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_IL6_JAK_STAT3_SIGNALING          
#  [4] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_APOPTOSIS                        
#  [5] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_MYOGENESIS                       
#  [6] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_COMPLEMENT                       
#  [7] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#  [8] HALLMARK_TNFA_SIGNALING_VIA_NFKB--HALLMARK_INFLAMMATORY_RESPONSE            
# + ... omitted several edges
# Number of geneset clusters:  10 
```

## Examples of different ways of loading in genesets

Example 1) Doing some manual editing of a predefined geneset.
```R
library(bowerbird)
# download from msigdb website
file <- 'c5.go.bp.v7.4.symbols.gmt'
# manualy read in to do some fine adjustments/filtering
geneset <- read_geneset(file) # reads in gene file manually
# do a bit of manual filtering
geneset <- geneset[grep('B_CELL|T_CELL|NATURAL_KILLER|ANTIBODY|ANTIGEN|LYMPHOCYTE|IMMUNE|INTERFERON|TOLL|INNATE|ADAPTIVE', names(geneset))]
geneset <- geneset[!grepl('TROPHOBLAST_CELL|FAT_CELL|ENT_CELL', names(geneset))]

bwr <- bower(geneset) 
bwr <- snn_graph(bwr)
bwr <- find_clusters(bwr)
bwr <- summarize_clusters(bwr)
```

Example 2) Extracting from `msigdbr`.
```R
library(bowerbird)
library(msigdbr)
GO <- data.frame(msigdbr::msigdbr(category = "C5", subcategory = "GO:BP"))
genesets <- GO[grep('B_CELL|T_CELL|NATURAL_KILLER|ANTIBODY|ANTIGEN|LYMPHOCYTE|IMMUNE|INTERFERON|TOLL|INNATE|ADAPTIVE', GO$gs_name), ]
genesets <- genesets[!grepl('TROPHOBLAST_CELL|FAT_CELL|ENT_CELL', genesets$gs_name), ]

# convert to list
gs_list <- lapply(unique(genesets$gs_name), function(x) genesets[genesets$gs_name %in% x, "gene_symbol"])
names(gs_list) <- unique(genesets$gs_name)

bwr <- bower(gs_list)
bwr <- snn_graph(bwr)
bwr <- find_clusters(bwr)
bwr <- summarize_clusters(bwr)
```

## Downstream analysis e.g. with AUCell
```R
# load a dummy dataset
library(ktplots)
data(kidneyimmune)

# an example workflow with AUCell
library(igraph)
library(ggraph)
library(AUCell)
library(Seurat)
library(SingleCellExperiment)
# kidneyimmune <- subset(kidneyimmune, idents = c('NK cell', 'B cell', 'Mast cell', 'CD8T cell', 'CD4T cell', 'NKT cell'))
sce <- as.SingleCellExperiment(kidneyimmune)
cells_rankings <- AUCell_buildRankings(counts(sce), nCores=3, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(bwr@genesets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

for (i in 1:nrow(cells_AUC@assays@data$AUC)) {
    geneSetsAUC <- cells_AUC@assays@data$AUC[i, ]
    colData(sce)[, row.names(cells_AUC@assays@data$AUC)[i]] <- geneSetsAUC
  }

AUC_mat <- t(data.frame(cells_AUC@assays@data$AUC))

scores <- do.call(cbind, lapply(unique(sce$celltype), function(x){
  rm <- rowMeans(t(AUC_mat)[, sce$celltype %in% x])
  return(rm)
}))

colnames(scores) <- unique(sce$celltype)
rownames(scores) <-names(bwr@genesets)
scores <- t(apply(as.data.frame(scores), 1, ktplots::range01))

# plot
plot_list <- lapply(colnames(scores), function(ds){
  set.seed(100)
  V(bwr@graph)$score <- scores[, ds]
  g <- ggraph(bwr@graph, 'igraph', algorithm = 'fr')
  g + # geom_edge_fan(aes(alpha = ..index..)) + 
    geom_edge_link(aes(width = weight), alpha = .25) + 
    geom_node_point(aes(color = score), size= 2) + 
    geom_node_text(aes(label = labels), size = 1) + 
    theme_bw() + 
    theme_void() + 
    scale_color_viridis() + 
    # scale_color_viridis(limits=c(0, 0.1), oob=scales::squish) + 
    scale_edge_width(range = c(0, 1)) +
    ggtitle(ds)
})
names(plot_list) <- colnames(scores)
cowplot::plot_grid(plotlist = plot_list)
```
