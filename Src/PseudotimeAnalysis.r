#### Functions for Pseudotime analysis

# Run Diffusion map with library("destiny")
# Original code from Stévant I et.al. Cell Rep. 2019.
run_diffMap <- function(data=data, condition=condition, sigma="local", k = 20){
    destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
    destinyObj$condition <- factor(condition)
    dm <- DiffusionMap(destinyObj, sigma, k)
    return(dm)
}

plot_eigenVal <- function(dm=dm){
    linepad <- .5
    plot(
    eigenvalues(dm),
    ylim = 0:1,
    pch = 20,
    xlab ='Diffusion component (DC)',
    ylab ='Eigenvalue'
    )
}

# Lineage prediction with library("slingshot")
# Original code from Stévant I et.al. Cell Rep. 2019.
rankKeepNA <- function(x) {
    return(
    ifelse(
    is.na(x),
    NA,
    rank(
    x,
    na.last = TRUE,
    ties.method="random"
    )
    )
    )
}

get_pseudotime <- function(pseudotime, wthres=wthres){
    pseudoT <- list()
    for(lineage in 1:length(pseudotime@curves))local({
        curve <- pseudotime@curves[[lineage]]
        lambda <- curve$lambda
        weight <- curve$w
        ps <- curve$lambda
        ps[weight < wthres] <- NA
        ps <- rankKeepNA(ps)
        pseudoT[[lineage]] <<- ps
    })
    df <- t(do.call("rbind",pseudoT))
    colnames(df) <- names(pseudotime@curves)
    return(df)
}

#### Visualization
library(grid)
library(tidyr)
library(ggplot2)
library(rgl)
library("plot3D")

# This function was produced by modifying code (Stévant I et.al. Cell Rep. 2019).
# bty = c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours, size=0.5, type = "s", legend = NULL){
    cond <- factor(condition)
    col <- factor(condition)
    levels(col) <- colours
    col <- as.vector(col)
    DCs <- paste("DC",dc, sep="")
    
    data <- data.frame(
    dm@eigenvectors[,DCs[1]],
    dm@eigenvectors[,DCs[2]],
    dm@eigenvectors[,DCs[3]]
    )
    colnames(data) <- DCs
    
    plot3d(
    data,
    bg=col,
    col=col,
    size=size,
    type = type,
    box = FALSE
    )
    
    if (is.null(legend)==FALSE){
        legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
    }
}

# bty = c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
plot_dm_3D_in_2D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours, outline.color = "grey20", size=1, pch = 21, theta = 40, phi = 40, bty ="b"){
    #cond <- factor(condition)
    cols <- factor(condition)
    levels(cols) <- colours
    cols <- as.vector(cols)
    DCs <- paste("DC",dc, sep="")
    
    data <- data.frame(
    dm@eigenvectors[,DCs[1]],
    dm@eigenvectors[,DCs[2]],
    dm@eigenvectors[,DCs[3]]
    )
    colnames(data) <- DCs
    
    if(pch == 21){
        scatter3D(x=as.matrix(data[,DCs[1]]), y=as.matrix(data[,DCs[2]]), z=as.matrix(data[,DCs[3]]), bg = cols, col = outline.color, pch=pch, cex = size, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, bty = bty)
    } else {
        scatter3D(x=as.matrix(data[,DCs[1]]), y=as.matrix(data[,DCs[2]]), z=as.matrix(data[,DCs[3]]), col = condition, col = cols, pch=pch, cex = size, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, bty = bty)
    }
    
}


# myFeaturePlot3D source code
SetQuantile <- function(cutoff, data) {
    if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
        this.quantile <- as.numeric(x = sub(
        pattern = 'q',
        replacement = '',
        x = as.character(x = cutoff)
        )) / 100
        data <- unlist(x = data)
        data <- data[data > 0]
        cutoff <- quantile(x = data, probs = this.quantile)
    }
    return(as.numeric(x = cutoff))
}

#object: Seurat object
#dataframe: diffusion map results
#DCno: DC
#features.use: gene name

myFeaturePlot3D <- function(dataframe, object, DCno, features.use, pt.size = 0.5, min.exp = -Inf, max.exp = Inf,hl.low.color="#d3d3d3", hl.high.color="red", type = "s") {
    data.plot=data.frame(dataframe[,DCno], FetchData(object = object, vars = features.use))
    data.plot$cell <- rownames(x = data.plot)
    data.plot %>% gather(key = "gene", value = "expression", -DCno, -cell) -> data.plot
    data.plot %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled.expression = scale(expression)) -> data.plot
    min.exp <- SetQuantile(cutoff = min.exp, data = data.plot$scaled.expression)
    max.exp <- SetQuantile(cutoff = max.exp, data = data.plot$scaled.expression)
    data.plot$gene <- factor(x = data.plot$gene, levels = features.use)
    data.plot$scaled.expression <- MinMax(data = data.plot$scaled.expression,min = min.exp,max = max.exp)
    sortlist <- order(data.plot$scaled.expression)
    data.plot <- data.plot[sortlist,]

    p=ggplot(data.plot,aes(x=DCno[1],y=DCno[2]))
    p=p+geom_point(aes(color = scaled.expression), shape=16, size=1)
    p=p+scale_colour_gradient(low = hl.low.color, high = hl.high.color)
    data.plot$color <- ggplot_build(p)$data[1][[1]]$colour
    p3 <- plot3d(data.plot[DCno], col = data.plot$color, size = pt.size, box=FALSE, type=type)
    return(p3)
}


#object: Seurat object
#dataframe: diffusion map results
#DCno: DC
#features.use: gene name
#bty = c("b", "b2", "f", "g", "bl", "bl2", "u", "n")

myFeaturePlot3D_in_Plot2D <- function(dataframe, object, DCno, features.use, pt.size = 1, min.exp = -Inf, max.exp = Inf, hl.low.color="#d3d3d3", hl.high.color="red", outline.color = "grey20", pch=21, theta = 40, phi = 40, bty ="b") {
    
    data.plot=data.frame(dataframe[,DCno], FetchData(object = object, vars = features.use))
    data.plot$cell <- rownames(x = data.plot)
    data.plot %>% gather(key = "gene", value = "expression", -DCno, -cell) -> data.plot
    data.plot %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled.expression = scale(expression)) -> data.plot
    min.exp <- SetQuantile(cutoff = min.exp, data = data.plot$scaled.expression)
    max.exp <- SetQuantile(cutoff = max.exp, data = data.plot$scaled.expression)
    data.plot$gene <- factor(x = data.plot$gene, levels = features.use)
    data.plot$scaled.expression <- MinMax(data = data.plot$scaled.expression,min = min.exp,max = max.exp)
    sortlist <- order(data.plot$scaled.expression)
    data.plot <- data.plot[sortlist,]
    
    p=ggplot(data.plot,aes(x=DCno[1],y=DCno[2]))
    p=p+geom_point(aes(color = scaled.expression), shape=16, size=1)
    p=p+scale_colour_gradient(low = hl.low.color, high = hl.high.color)
    data.plot$color <- ggplot_build(p)$data[1][[1]]$colour
    
    if (pch==21){
        p3 <- scatter3D(x=as.matrix(data.plot[,DCno[1]]), y=as.matrix(data.plot[,DCno[2]]), z=as.matrix(data.plot[,DCno[3]]), bg = data.plot$color, size = pt.size, col = outline.color, pch=pch, theta = theta, phi = phi, fov=90, cex = pt.size, xlab = DCno[1], ylab = DCno[2], zlab = DCno[3], bty = bty, main=features.use)
        
    } else {
        p3 <- scatter3D(x=as.matrix(data.plot[,DCno[1]]), y=as.matrix(data.plot[,DCno[2]]), z=as.matrix(data.plot[,DCno[3]]), colvar=data.plot$scaled.expression, col = data.plot$color, size = pt.size, pch=pch, theta = theta, phi = phi, fov=90, cex = pt.size, xlab = DCno[1], ylab = DCno[2], zlab = DCno[3], bty = bty,main=features.use)
    }
    
    return(p3)
}






#### plot_cells_per_lineage in tSNEplot

highlight_cells_per_lineage_tSNEplot <- function(object, highlight.cells, pt.size = 2, pt.shape = 21, outline.color = "grey30", bg.color = "#d3d3d3", hi.color, outline.size = 0.2,title="tSNE map"){
    
    cell.use.s <- highlight.cells
    
    data.plot=data.frame(object@reductions$tsne@cell.embeddings,row.names(object@reductions$pca@cell.embeddings),object@active.ident)
    colnames(data.plot) <- c("tSNE_1","tSNE_2","Cell","ident")
    
    g=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2))
    g=g+geom_point(size=pt.size, shape=pt.shape, color=outline.color, fill=bg.color, stroke = outline.size)
    g=g+geom_point(data=subset(data.plot, data.plot$Cell %in% cell.use.s), aes(fill =ident), color=outline.color, size=pt.size, shape=pt.shape, stroke = outline.size)
    g=g+scale_fill_manual(values = hi.color,name="")+ggtitle(title)
    g=g+theme_grey()+
    theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=16),
    legend.text = element_text(size =16),
    legend.title = element_text(size =16 ,face="bold"),
    plot.title = element_text(size=18, face="bold", hjust = 0.5),
    aspect.ratio=1
    )
    return(g)
}




#### plot_cells_per_lineage in 3D diffusionmap

## ex
# curves <- slingCurves(male_lineage)
# add.curve <- curves$curve1

#bty = c("b", "b2", "f", "g", "bl", "bl2", "u", "n")

highlight_cells_per_lineage_dmplot3D <- function(object, dm=dm, dc=c(1:3), highlight.cells, condition=condition, colours=colours, lineage.curve=NULL, size=1, pch = 21, outline.color = "grey20", bg.color = "#d3d3d3", title="Diffusion map",phi = 40, theta = 40, add.curve=NULL, bty ="b"){
    
    cell.use.s <- highlight.cells
    
    cols <- factor(condition)
    levels(cols) <- colours
    cols <- as.vector(cols)
    DCs <- paste("DC",dc, sep="")
    
    data <- data.frame(
    dm@eigenvectors[,DCs[1]],
    dm@eigenvectors[,DCs[2]],
    dm@eigenvectors[,DCs[3]],
    Cell=row.names(object@reductions$pca@cell.embeddings),
    Ident=condition,
    Color=cols
    )
    colnames(data) <- c(DCs,"Cell","ident","Color")
    
    data2=subset(data, data$Cell %in%　cell.use.s)
    
    scatter3D(x=as.matrix(data[,DCs[1]]), y=as.matrix(data[,DCs[2]]), z=as.matrix(data[,DCs[3]]), bg = bg.color, col = outline.color, pch=pch, cex = size, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, main =title, bty = bty)
    scatter3D(x=as.matrix(data2[,DCs[1]]), y=as.matrix(data2[,DCs[2]]), z=as.matrix(data2[,DCs[3]]), bg = as.matrix(data2[,"Color"]), col = outline.color, pch=pch, cex = size, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, add=TRUE, bty = bty)
    
    if(is.null(add.curve)==FALSE){
        c <- add.curve
        scatter3D(x = c$s[c$ord,DCs[1]], y =c$s[c$ord,DCs[2]], z = c$s[c$ord,DCs[3]], col = "black", type = "l", ticktype = "detailed", lwd = 4, add = TRUE, xlab = DCs[1], ylab = DCs[2], zlab = DCs[3], theta =theta, phi = phi, bty = bty)
    }
    
}
