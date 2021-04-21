library(grid)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tibble)

# for monocle object
monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()) +
    #theme(axis.line.x = element_line(size=0.25, color="black")) +
    #theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

myplot_cell_clusters <- function(cds,
x=1,
y=2,
color_by="Cluster",
markers=NULL,
show_cell_names=FALSE,
cell_size=1.5,
cell_name_size=2,
n.col=3,
legend.position="right",
pch.use=21,
viridis.option = NULL,
hl.low.color="#d3d3d3",
hl.high.color="red",
outline.color="grey30",
outline.size=0.2,
manualcolor=NULL,
...){
    if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
        stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
    }
    
    gene_short_name <- NULL
    sample_name <- NULL
    data_dim_1 <- NULL
    data_dim_2 <- NULL
    
    #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
    lib_info <- pData(cds)
    
    tSNE_dim_coords <- reducedDimA(cds)
    data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- colnames(cds)
    data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")
    
    data_df$CellType <- as.factor(data_df$CellType)
    data_df$CellType <- factor(data_df$CellType, levels=c("SCprogenitor","IFE","placode","Unknown"))
    data_df <- data_df[order(data_df$CellType, decreasing=TRUE),]
    
    
    markers_exprs <- NULL
    if (is.null(markers) == FALSE){
        markers_fData <- subset(fData(cds), gene_short_name %in% markers)
        if (nrow(markers_fData) >= 1){
            cds_subset <- cds[row.names(markers_fData),]
            if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
                integer_expression <- TRUE
            }
            else {
                integer_expression <- FALSE
                
            }
            if (integer_expression) {
                cds_exprs <- exprs(cds_subset)
                
                if (is.null(sizeFactors(cds_subset))) {
                    stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
                }
                cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
                
                cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
            }
            else {
                cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
            }
            markers_exprs <- cds_exprs
            #markers_exprs <- reshape2::melt(as.matrix(cds_exprs))
            colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
            markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
            #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
            markers_exprs$feature_label<-factor(markers_exprs$feature_label,levels=markers)
            sortlist <- order(markers_exprs$value)
            markers_exprs <- markers_exprs[sortlist,]
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
        data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
        sortlist <- order(data_df$value)
        data_df <- data_df[sortlist,]
        
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label, ncol = n.col)
    }else{
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
    }
    
    # FIXME: setting size here overrides the marker expression funtionality.
    # Don't do it!
    
    if (pch.use != "16") {
        
        if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
            g <- g + geom_point(aes(fill=log10(value + 0.1)), size=I(cell_size),
            shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
            if (is.null(viridis.option) == FALSE) {g=g+scale_fill_viridis_c(name = paste0("log10(value + 0.1)"), option=viridis.option)}
            else {g=g+scale_fill_gradient(low = hl.low.color, high = hl.high.color)}

        }else if (is.null(manualcolor)==TRUE) {
            g <- g + geom_point(aes_string(fill = color_by), size=I(cell_size),
            shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
        }else {
            g <- g + geom_point(aes_string(fill = color_by), size=I(cell_size),
            shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE) +
            scale_fill_manual(values=manualcolor)
        }
        
    }
    else{
        if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
            g <- g + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE)
            if (is.null(viridis.option) == FALSE) {g=g+scale_color_viridis_c(name = paste0("log10(value + 0.1)"), option=viridis.option)}
            else {g=g+scale_colour_gradient(low = hl.low.color, high = hl.high.color)}
            
        }else if (is.null(manualcolor)==TRUE) {
            g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
        }else {
            g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE) +
            scale_color_manual(values=manualcolor)
        }
    }

    
    g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste("Component", x)) +
    ylab(paste("Component", y)) +
    theme(legend.position=legend.position, legend.key.height=grid::unit(0.35, "in"), legend.key.width=grid::unit(0.4, "in")) +
    guides(color = guide_colorbar(label.theme = element_text(size = 10, angle = 0)))+
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(text = element_text(size = 15)) +
    theme(plot.title=element_text(face="bold", size=cell_name_size))
    return(g)
}




myplot_cell_trajectory <- function(cds, x=1, y=2, color_by="State", show_tree=TRUE,show_backbone=TRUE, backbone_color="black", markers=NULL, use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names=FALSE, show_state_number = FALSE, cell_size=1.5, cell_link_size=0.75, cell_name_size=2, state_number_size = 2.9, show_branch_points=TRUE, theta = 0, n.col=3, legend.position="right", viridis.option = NULL, pch.use=21, outline.color="grey30", outline.size=0.2, pt.color=brewer.pal(9,"YlGn"), hl.low.color="#d3d3d3", hl.high.color="red",manualcolor=NULL) {
    requireNamespace("igraph")
    gene_short_name <- NA
    sample_name <- NA
    sample_state <- pData(cds)$State
    data_dim_1 <- NA
    data_dim_2 <- NA
    
    #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
    lib_info_with_pseudo <- pData(cds)
    
    if (is.null(cds@dim_reduce_type)){
        stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }
    
    if (cds@dim_reduce_type == "ICA"){
        reduced_dim_coords <- reducedDimS(cds)
    } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
        reduced_dim_coords <- reducedDimK(cds)
    } else {
        stop("Error: unrecognized dimensionality reduction method.")
    }
    
    ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
    
    dp_mst <- minSpanningTree(cds)
    
    if (is.null(dp_mst)){
        stop("You must first call orderCells() before using this function")
    }
    
    edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>% dplyr::select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    dplyr::left_join(ica_space_df %>% dplyr::select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
    
    data_df <- t(monocle::reducedDimS(cds)) %>%
    as.data.frame() %>%
    dplyr::select_(data_dim_1 = x, data_dim_2 = y) %>%
    rownames_to_column("sample_name") %>%
    dplyr::mutate(sample_state) %>%
    dplyr::left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")
    
    data_df$CellType <- as.factor(data_df$CellType)
    data_df$CellType <- factor(data_df$CellType, levels=c("SCprogenitor","IFE","placode","Unknown"))
    data_df <- data_df[order(data_df$CellType, decreasing=TRUE),]
    
    
    return_rotation_mat <- function(theta) {
        theta <- theta / 180 * pi
        matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
    }
    rot_mat <- return_rotation_mat(theta)
    
    cn1 <- c("data_dim_1", "data_dim_2")
    cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
    cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
    data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
    edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
    edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
    
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
            colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
            markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
            #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
            markers_exprs$feature_label<-factor(markers_exprs$feature_label,levels=markers)
        }
    }
    
    if (pch.use != "16") {
        
        if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
            data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
            sortlist <- order(data_df$value)
            data_df <- data_df[sortlist,]
            
            if(use_color_gradient) {
                if(markers_linear){
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                    geom_point(aes(fill= value), size=I(cell_size),shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
                    if (is.null(viridis.option) == FALSE) {g=g+scale_fill_viridis_c(name = paste0("value"),option=viridis.option)}
                    else {g=g+scale_fill_gradient(name = paste0("value"),low = hl.low.color, high = hl.high.color)}
                    #scale_fill_viridis(name = paste0("value"), ...) +
                    #scale_fill_gradient(name = paste0("value"),low="#d3d3d3", high="red") +
                    g=g+facet_wrap(~feature_label, ncol = n.col)
                } else {
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                    geom_point(aes(fill=log10(value + 0.1)), size=I(cell_size),shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
                    if (is.null(viridis.option) == FALSE) {g=g+scale_fill_viridis_c(name = paste0("log10(value + 0.1)"),option=viridis.option)}
                    else {g=g+scale_fill_gradient(name = paste0("log10(value + 0.1)"),low = hl.low.color, high = hl.high.color)}
                    #scale_fill_viridis(name = name = paste0("log10(value + 0.1)"), ...) +
                    #scale_fill_gradient(name = name = paste0("log10(value + 0.1)"),low="#d3d3d3", high="red") +
                    g=g+facet_wrap(~feature_label, ncol = n.col)
                }
            } else {
                if(markers_linear){
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label, ncol = n.col)
                } else {
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label, ncol = n.col)
                }
            }
        } else if (is.null(manualcolor)==TRUE) {
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
        } else {
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))+scale_fill_manual(values=manualcolor)}
        
    } else {
        if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
            data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
            if(use_color_gradient) {
                if(markers_linear){
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                    geom_point(aes(color= value), size=I(cell_size), na.rm = TRUE)
                    
                    if (is.null(viridis.option) == FALSE) {g=g+scale_color_viridis_c(name = paste0("value"),option=viridis.option)}
                    else {g=g+scale_colour_gradient(name = paste0("value"),low = hl.low.color, high = hl.high.color)}
                    #scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~feature_label, ncol = n.col)
                    g=g+facet_wrap(~feature_label, ncol = n.col)
                    
                    
                } else {
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                    geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE)
                    
                    if (is.null(viridis.option) == FALSE) {g=g+scale_color_viridis_c(name = paste0("log10(value + 0.1)"),option=viridis.option)}
                    else {g=g+scale_colour_gradient(name = paste0("log10(value + 0.1)"),low = hl.low.color, high = hl.high.color)}
                    #scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label, ncol = n.col)
                    g=g+facet_wrap(~feature_label, ncol = n.col)
                    
                }
            } else {
                if(markers_linear){
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label, ncol = n.col)
                } else {
                    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label, ncol = n.col)
                }
            }
        } else if (is.null(manualcolor)==TRUE) {
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
        } else {
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))+scale_color_manual(values=manualcolor)
        }
    }
    
    
    
    
    if (show_tree){
        g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
    }
    
    
    
    
    # FIXME: setting size here overrides the marker expression funtionality.
    # Don't do it!
    if (pch.use != "16") {
        
        if (is.null(markers_exprs) == TRUE){
            if(use_color_gradient) {
                g <- g + geom_point(aes_string(fill = color_by), size=I(cell_size),
                shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
                #g <- g + scale_fill_gradient(low="#d3d3d3", high="red")
                g <- g + scale_fill_gradientn(colours=pt.color)
            } else {
                g <- g + geom_point(aes_string(fill = color_by), size=I(cell_size),
                shape=pch.use, color=outline.color, stroke = outline.size, na.rm = TRUE)
            }
        }
        
    } else {
        if (is.null(markers_exprs) == TRUE){
            if(use_color_gradient) {
                g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
                #g <- g + scale_fill_gradient(low="#d3d3d3", high="red")
                g <- g + scale_colour_gradientn(colours=pt.color)
            } else {
                g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
                #g <- g + scale_colour_gradientn(colours=pt.color)
            }
        }
    }
    
    
    
    if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
        mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
        branch_point_df <- ica_space_df %>%
        dplyr::slice(match(mst_branch_nodes, sample_name)) %>%
        dplyr::mutate(branch_point_idx = seq_len(n()))
        
        g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
        size=5, na.rm=TRUE, branch_point_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
        size=4, color="white", na.rm=TRUE, branch_point_df)
    }
    if (show_cell_names){
        g <- g + geom_text(aes(label=sample_name), size=cell_name_size)
    }
    if (show_state_number){
        g <- g + geom_text(aes(label = sample_state), size = state_number_size)
    }
    
    g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste("Component", x)) +
    ylab(paste("Component", y)) +
    theme(legend.position=legend.position, legend.key.height=grid::unit(0.35, "in"), legend.key.width=grid::unit(0.4, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
    return(g)}

