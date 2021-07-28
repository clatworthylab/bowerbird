#' @title plot_graph
#' @description
#' Plot the geneset kNN network as a graph.
#'
#' @name plot_graph
#' @rdname plot_graph
#' @param bower BOWER object
#' @param colorby name of attribute for colors of nodes.
#' @param node.size name of attribute for size of nodes or numerical value/vector.
#' @param node.label name of attribute for colors of nodes.
#' @param node.label.size size of node label.
#' @param node.max_size maximum node size.
#' @param edge.width name of attribute for width of edges.
#' @param edge.width.thickness.range scale for range of widths of edges.
#' @param edge.alpha transparency of edges.
#' @param guides.colour boolean. whether or not to display guide for vertex colour.
#' @details
#' Just a ggraph wrapper to plot the graph.
#' @return Returns a ggraph plot.
#' @examples
#' gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "bowerbird")
#' bwr <- bower(gmt_file)
#' bwr <- snn_graph(bwr)
#' bwr <- find_clusters(bwr)
#' bwr <- summarize_clusters(bwr, ncpus = 1)
#' plot_graph(bwr)
#' @import ggraph ggplot2
#' @export

plot_graph.BOWER <- function(bower, colorby = 'cluster', mode = c('default','gsea'), gsea.slot = c('padj', 'pval', 'NES', 'ES'), gsea.pval.cutoff = .25, node.size = 'geneset_size', node.label = 'labels', node.label.size = 3, title = NULL, title.size = 3, edge.width = 'weight', edge.width.thickness.range = c(0, 1), edge.alpha = .25, dot_max = 8, guides.colour = FALSE){
    if(is.null(bower@scores)){
        if (colorby %in% colnames(bower@.graph_data)){
            g <- ggraph(bower@.graph_data) + geom_edge_link(aes(width = get(edge.width)), alpha = edge.alpha)
            if (is.numeric(node.size)){
                g <- g + geom_node_point(aes(color = get(colorby)), size = node.size)
            } else {
                g <- g + geom_node_point(aes(color = get(colorby), size = get(node.size)))
            }
            g <- g + 
            scale_size_area(max_size = dot_max) +
            geom_node_text(aes(label = get(node.label)), size = node.label.size) +
            theme_bw() +
            theme_void() +
            scale_edge_width(range = edge.width.thickness.range) + 
            guides(colour=guides.colour, size=guide_legend(title=node.size), edge_width=guide_legend(title=edge.width))
            if (!is.null(title)){
                g <- g + ggtitle(title) + theme(plot.title = element_text(size = title.size))
            }
        } else {
            stop(paste0(colorby, " does not exist in the BOWER object."))
        }
    } else {
        if (!colorby %in% colnames(bower@.graph_data)){
            if(colorby %in% colnames(bower@scores) | colorby %in% colnames(bower@scores[[1]])) {
                if (mode[1] == 'default') {
                    bower@.graph_data[,colorby] = bower@scores[,colorby]
                } else if (mode[1] == 'gsea'){
                    tmp <- bower@scores[[gsea.slot[1]]][,colorby]                
                    if (gsea.slot[1] %in% c('pval', 'padj')){
                        tmp[tmp >= gsea.pval.cutoff] <- NA
                    }
                    bower@.graph_data[,colorby] = tmp
                }
            } else {
                stop(paste0(colorby, " does not exist in the BOWER object."))
            }            
        }
        g <- ggraph(bower@.graph_data) + geom_edge_link(aes(width = get(edge.width)), alpha = edge.alpha)
        if (is.numeric(node.size)){
            g <- g + geom_node_point(aes(color = get(colorby)), size = node.size)
        } else {
            g <- g + geom_node_point(aes(color = get(colorby), size = get(node.size)))
        }
        if (node.label %in% colnames(g$data)){
            g <- g + geom_node_text(aes(label = get(node.label)), size = node.label.size)
        } else {
            stop(paste0(node.label, " does not exist in the BOWER object."))
        }
        g <- g + 
        scale_size_area(max_size = dot_max) +
        theme_bw() +
        theme_void() +
        scale_edge_width(range = edge.width.thickness.range) + 
        guides(colour=guides.colour, size=guide_legend(title=node.size), edge_width=guide_legend(title=edge.width))
        if (!is.null(title)){
            g <- g + ggtitle(title) + theme(plot.title = element_text(size = title.size))
        }
    }
    return(g)
}

#' @rdname plot_graph
#' @param ... passed to ggraph::ggraph
#' @export
plot_graph.igraph <- function(gr, colorby = 'cluster', node.size = 'geneset_size', node.label = 'labels', node.label.size = 3, title = NULL, title.size = 3, edge.width = 'weight', edge.width.thickness.range = c(0, 1), edge.alpha = .25, dot_max = 8, guides.colour = FALSE, ...){
    if (colorby %in% names(igraph::vertex_attr(gr))){
        g <- ggraph(gr, ...) + geom_edge_link(aes(width = get(edge.width)), alpha = edge.alpha)
        if (is.numeric(node.size)){
            g <- g + geom_node_point(aes(color = get(colorby)), size = node.size)
        } else {
            if (node.size %in% names(igraph::vertex_attr(gr))){
                g <- g + geom_node_point(aes(color = get(colorby), size = get(node.size)))
            } else {
                stop(paste0(node.size, " does not exist in the igraph vertex attributes."))
            }
        }        

        if (node.label %in% colnames(g$data)){
            g <- g + geom_node_text(aes(label = get(node.label)), size = node.label.size)
        } else {
            stop(paste0(node.label, " does not exist in the igraph vertex attributes."))
        }
        g <- g + 
        scale_size_area(max_size = dot_max) +
        theme_bw() +
        theme_void() +
        scale_edge_width(range = edge.width.thickness.range) + 
        guides(colour=guides.colour, size=guide_legend(title=node.size), edge_width=guide_legend(title=edge.width))
        if (!is.null(title)){
            g <- g + ggtitle(title) + theme(plot.title = element_text(size = title.size))
        }
    } else {
        stop(paste0(colorby, " does not exist in the igraph object."))
    } 
    return(g)
}

# unused code to split the plots
#
# colors = scales::hue_pal()(length(unique(V(bwr@graph)$cluster)))
# names(colors) <- unique(V(bwr@graph)$cluster)
# gl <- list()
# for (i in seq_along(g)){
# gl[[i]] <- ggraph(g[[i]], 'igraph', algorithm = 'fr')
# lbl <- unique(gl[[i]]$data$labels)
# lbl <- lbl[!lbl %in% ""]
# gl[[i]] <- gl[[i]] + geom_edge_link(aes(width = weight), alpha = .25) +
#     geom_node_point(aes(size = geneset_size), color = colors[unique(gl[[i]]$data$cluster)]) + 
#     scale_size_area(max_size = 8) +
#     theme_bw() +
#     theme_void() +
#     scale_edge_width(range = c(0, 1)) + 
#     guides(colour=FALSE) + 
#     ggtitle(lbl) +
#     theme(legend.position = 'none', plot.title = element_text(size = 4)) +
#     coord_cartesian(clip = "off")
# }

# library(cowplot)
# plot_grid(plotlist = gl, scale = .9)
