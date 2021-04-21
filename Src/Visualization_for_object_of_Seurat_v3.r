

library(grid)
library(tidyr)
library(ggplot2)

# myTSNEPlot source code
# default, group.by = "ident"

myTSNEPlot <- function(object, pt.size = 2,pt.shape = 21, pt.color=NULL, outline.color = "grey30",outline.size=0.2,label.size = 8,label.color="black",viridis.option=NULL, viridis.dir=1, do.label=TRUE,group.by = "ident", no.rect = TRUE, theme.grey=TRUE, title=NULL, legend_ncol=1, xlab="tSNE_1", ylab="tSNE_2", no.legend=FALSE){
    data.plot <- data.frame(object@reductions$tsne@cell.embeddings,ident=factor(Idents(object = object)))
    colnames(data.plot) <- c("x","y","ident")
    data.plot$cell <- rownames(x = data.plot)
    
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(
        object = object,
        vars = group.by
        )[data.plot$cell, 1])
        data.plot$ident <- ident.use
    }
    
    p=ggplot(data.plot,aes(x=x,y=y))
    
    if (pt.shape != "16") {
        p=p+geom_point(size=pt.size, pch = pt.shape, color=outline.color, stroke = outline.size, aes(fill=factor(data.plot$ident)))
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_fill_viridis_d(name ="",option=viridis.option, direction = viridis.dir)}
        if (is.null(pt.color) == FALSE) {p=p=p+scale_fill_manual(name ="",values = pt.color)}
    }
    else{
        p=p+geom_point(size=pt.size, pch = pt.shape, aes(color=factor(data.plot$ident)))
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_colour_viridis_d(name ="",option=viridis.option, direction = viridis.dir)}
        if (is.null(pt.color) == FALSE) {p=p=p+scale_colour_manual(name ="",values = pt.color)}
    }
    
    if (no.rect == TRUE) {
        p=p+theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white",color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.grid.major = element_line(color = NA),
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        aspect.ratio=1,
        legend.position = 'none')
    }else if (theme.grey != TRUE){
        p=p+theme_bw()+theme(
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1))+
        #labs(title=title, fill = group.by, color = group.by)+
        labs(title=title, fill = "", color = "")+
        xlab(xlab)+ylab(ylab)
        if(no.legend != TRUE){p=p+guides(fill=guide_legend(ncol=legend_ncol, byrow=TRUE))}
        else{p=p+theme(legend.position = 'none')}
    }else{
        p=p+theme_grey()+
        theme(aspect.ratio=1)+
        #labs(title=title, fill = group.by, color = group.by)+
        labs(title=title, fill = "", color = "")+
        xlab(xlab)+ylab(ylab)
        if(no.legend != TRUE){p=p+guides(fill=guide_legend(ncol=legend_ncol, byrow=TRUE))}
        else{p=p+theme(legend.position = 'none')}
    }
    
    if (do.label==TRUE) {
        data.plot %>% dplyr::group_by(ident) %>% dplyr::summarize(x = median(x = x), y = median(x = y)) -> centers
        p <- p +
        geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
        geom_text(data = centers, mapping = aes(label = ident), size = label.size, color=label.color, fontface='bold')}
    
    return(p)
}


# myFeaturePlot source code

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


myFeaturePlot <- function(object, features.use, pt.size = 1, pch.use = 16, nCol = 1, min.exp = -Inf, max.exp = Inf,
hl.low.color="#d3d3d3", hl.high.color="red", outline.color="grey70", viridis.option = NULL, outline.size=0.2, title.size=30,legend.position=NULL) {
    if (is.null(nCol)) {
        nCol <- 2
        if (length(features.use) == 1) ncol <- 1
        if (length(features.use) > 6) nCol <- 3
        if (length(features.use) > 9) nCol <- 4
    }
    num.row <- floor(length(features.use) / nCol - 1e-5) + 1
    pList <- lapply(features.use,function(x) SingleFeaturePlot(object,x,pt.size,pch.use,min.exp,max.exp,hl.low.color,hl.high.color,outline.color,viridis.option,outline.size,title.size,legend.position))
    return(plot_grid(plotlist = pList, ncol = nCol))
}
SingleFeaturePlot <-  function(object, features.use, pt.size, pch.use, min.exp, max.exp,hl.low.color,hl.high.color,outline.color,viridis.option,outline.size,title.size,legend.position){
    
    data.plot=data.frame(object@reductions$tsne@cell.embeddings,
    FetchData(object = object, vars = features.use))
    data.plot$cell <- rownames(x = data.plot)
    data.plot %>% tidyr::gather(key = "gene", value = "expression", -tSNE_1, -tSNE_2, -cell) -> data.plot
    data.plot %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled.expression = scale(expression)) -> data.plot
    min.exp <- SetQuantile(cutoff = min.exp, data = data.plot$scaled.expression)
    max.exp <- SetQuantile(cutoff = max.exp, data = data.plot$scaled.expression)
    data.plot$gene <- factor(x = data.plot$gene, levels = features.use)
    data.plot$scaled.expression <- MinMax(
    data = data.plot$scaled.expression,
    min = min.exp,
    max = max.exp
    )
    sortlist <- order(data.plot$scaled.expression)
    data.plot <- data.plot[sortlist,]
    
    p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2))
    
    if (pch.use != "16") {
        p=p+geom_point(aes(fill = scaled.expression),
        shape=pch.use, size=pt.size, color=outline.color, stroke = outline.size)
        if (is.null(viridis.option) == FALSE) {p=p+scale_fill_viridis_c(option=viridis.option)}
        else {p=p+scale_fill_gradient(low = hl.low.color, high = hl.high.color)}
    }
    else{
        p=p+geom_point(aes(color = scaled.expression), shape=pch.use, size=pt.size)
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_colour_viridis_c(option=viridis.option)}
        else {p=p+scale_colour_gradient(low = hl.low.color, high = hl.high.color)}
    }
    
    p=p+theme(panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white",color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.grid.major = element_line(color = NA),
    panel.border = element_blank(),
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.title=element_text(face="bold.italic", size=title.size, hjust = 0.5))
    
    if (is.null(legend.position)==TRUE){
        p=p+theme(legend.position = 'none')
    }
    else {p=p+theme(legend.position = legend.position)}
    
    p=p+labs(title = features.use)
    return(p)}



# mySingleFeaturePlot source code

mySingleFeaturePlot <-  function(object, features.use, pt.size = 2, pch.use = 21, min.exp = -Inf, max.exp = Inf, hl.low.color="#d3d3d3", hl.high.color="red", viridis.option = NULL, outline.color="grey30", outline.size=0.2, title.size=20, no.rect = FALSE, theme.grey=TRUE, legend_ncol=1, legend.position=NULL, xlab="tSNE_1", ylab="tSNE_2"){
    data.plot=data.frame(object@reductions$tsne@cell.embeddings,
    FetchData(object = object, vars = features.use))
    data.plot$cell <- rownames(x = data.plot)
    data.plot %>% gather(key = "gene", value = "expression", -tSNE_1, -tSNE_2, -cell) -> data.plot
    data.plot %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled.expression = scale(expression)) -> data.plot
    min.exp <- SetQuantile(cutoff = min.exp, data = data.plot$scaled.expression)
    max.exp <- SetQuantile(cutoff = max.exp, data = data.plot$scaled.expression)
    data.plot$gene <- factor(x = data.plot$gene, levels = features.use)
    data.plot$scaled.expression <- MinMax(
    data = data.plot$scaled.expression,
    min = min.exp,
    max = max.exp
    )
    sortlist <- order(data.plot$scaled.expression)
    data.plot <- data.plot[sortlist,]
    
    p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2))
    
    if (pch.use != "16") {
        p=p+geom_point(aes(fill = scaled.expression),
        shape=pch.use, size=pt.size, color=outline.color, stroke = outline.size)
        if (is.null(viridis.option) == FALSE) {p=p+scale_fill_viridis_c(option=viridis.option)}
        else {p=p+scale_fill_gradient(name = paste0("value"),low = hl.low.color, high = hl.high.color)}
    }
    else{
        p=p+geom_point(aes(color = scaled.expression), shape=pch.use, size=pt.size)
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_colour_viridis_c(option=viridis.option)}
        else {p=p+scale_colour_gradient(name = paste0("value"),low = hl.low.color, high = hl.high.color)}
    }
    
    if (no.rect == TRUE) {
        p=p+theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white",color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.grid.major = element_line(color = NA),
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(face="bold.italic", size=title.size))
    }else if (theme.grey != TRUE){
        p=p+theme_bw()+theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1),
        plot.title=element_text(face="bold.italic", size=title.size))+
        #labs(title=title, fill = group.by, color = group.by)+
        xlab(xlab)+ylab(ylab)
    }else{
        p=p+theme_grey()+
        theme(plot.title=element_text(face="bold.italic", size=title.size))+
        #labs(title=title, fill = group.by, color = group.by)+
        xlab(xlab)+ylab(ylab)
    }
    
    if (is.null(legend.position)==TRUE){
        p=p+theme(legend.position = 'none')
    }
    else {p=p+theme(legend.position = legend.position)}
    
    p=p+labs(title = features.use)
    return(p)}



# StagePlot3Mod source code
# Scaling data in all stage (all data)

library(grid)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# pch.use=16 or 21
StagePlot3Mod <- function(object, stage.use, features.use, pt.size = 1, pch.use = 21, nCol = 1, min.exp = -Inf, max.exp = Inf,
bgp.color="#dcdcdc", hl.low.color="yellow", hl.high.color="red", outline.color="grey80", outline.size=0.2, title.size=30)
{
    if (is.null(nCol)) {
        nCol <- 2
        if (length(stage.use) == 1) ncol <- 1
        if (length(stage.use) > 6) nCol <- 3
        if (length(stage.use) > 9) nCol <- 4
    }
    num.row <- floor(length(stage.use) / nCol - 1e-5) + 1
    pList <- lapply(stage.use,function(x)
    SingleStagePlot4(object,x, features.use, pt.size, pch.use, min.exp, max.exp, bgp.color, hl.low.color, hl.high.color, outline.color, outline.size, title.size))
    return(plot_grid(plotlist = pList, ncol = nCol))
}
SingleStagePlot4 <-  function(object,stage.use.s, features.use, pt.size, pch.use, min.exp, max.exp, bgp.color, hl.low.color, hl.high.color,
outline.color, outline.size, title.size)
{#Prepare dataframe for plot
    #data.plot=data.frame(object@reductions$tsne@cell.embeddings,stage=object$stage_new,
    data.plot=data.frame(object@reductions$tsne@cell.embeddings,stage=object$stage,
    FetchData(object = object, vars = features.use))
    data.plot$cell <- rownames(x = data.plot)
    data.plot %>% gather(key = "gene", value = "expression", -tSNE_1, -tSNE_2, -stage, -cell) -> data.plot
    data.plot %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled.expression = scale(expression)) -> data.plot
    min.exp <- SetQuantile(cutoff = min.exp, data = data.plot$scaled.expression)
    max.exp <- SetQuantile(cutoff = max.exp, data = data.plot$scaled.expression)
    data.plot$gene <- factor(x = data.plot$gene, levels = features.use)
    data.plot$scaled.expression <- MinMax(
    data = data.plot$scaled.expression,
    min = min.exp,
    max = max.exp
    )
    sortlist <- order(data.plot$scaled.expression)
    data.plot <- data.plot[sortlist,]
    data.plot$scaled.expression01 <- range01(data.plot$scaled.expression)
    data.plot$scaled.color <- seq_gradient_pal(hl.low.color, hl.high.color)(data.plot$scaled.expression01)
    
    p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2))
    
    if (pch.use != "16") {
        p=p+geom_point(color=outline.color, fill=bgp.color, size=pt.size, shape=pch.use, stroke = outline.size)
        p=p+geom_point(data=subset(data.plot, data.plot$stage==stage.use.s),
        fill=subset(data.plot$scaled.color, data.plot$stage==stage.use.s),
        shape=pch.use, size=pt.size, color=outline.color, stroke = outline.size)}
    else{
        p=p+geom_point(colour=bgp.color,size=pt.size, shape=pch.use, stroke = outline.size)
        p=p+geom_point(data=subset(data.plot, data.plot$stage==stage.use.s),
        color=subset(data.plot$scaled.color, data.plot$stage==stage.use.s), shape=pch.use, size=pt.size, stroke = outline.size)}
    
    p=p+theme(panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white",color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.grid.major = element_line(color = NA),
    panel.border = element_blank(),
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),axis.title.y=element_blank(),
    legend.position = 'none',
    plot.title=element_text(face="bold", size=title.size))
    
    p=p+labs(title = stage.use.s)
    return(p)
}








# highlightTSNEPlot source code

highlightTSNEPlot <- function(object, pt.size = 2,pt.shape = 21,bg.color="#d3d3d3",outline.color = "grey30",outline.size=0.2,label.size = 8,label.color="black",do.label=FALSE,group.by = "ident", no.rect = FALSE, theme.grey=TRUE, title=NULL, legend_ncol=1, highlight.use=NULL, highlight.cell=NULL, highlight.col=NULL, no.legend = FALSE, xlab="tSNE_1", ylab="tSNE_2"){
    
    data.plot=data.frame(object@reductions$tsne@cell.embeddings,object@active.ident)
    colnames(data.plot) <- c("x","y","ident")
    
    data.plot$cell <- rownames(x = data.plot)
    
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(
        object = object,
        vars = group.by
        )[data.plot$cell, 1])
        data.plot$ident <- ident.use
    }
    
    p=ggplot(data.plot,aes(x=x,y=y))
    
    if (is.null(highlight.cell) == FALSE) {
        
        if (pt.shape != "16") {
            p=p+geom_point(size=pt.size, shape=pt.shape, color=outline.color, fill=bg.color, stroke = outline.size)
            
            if (is.null(highlight.col) == TRUE) {
                p=p+geom_point(data=subset(data.plot, row.names(data.plot) %in% highlight.cell), fill="red", color=outline.color,
                size=pt.size, shape=pt.shape, stroke = outline.size)
            } else {
                p=p+geom_point(data=subset(data.plot, row.names(data.plot) %in% highlight.cell), fill=highlight.col, color=outline.color,
                size=pt.size, shape=pt.shape, stroke = outline.size)
            }
        }
        
        else{
            p=p+geom_point(size=pt.size, shape=pt.shape, color=bg.color)
            
            if (is.null(highlight.col) == TRUE) {
                p=p+geom_point(data=subset(data.plot, row.names(data.plot) %in% highlight.cell), color="red", size=pt.size, shape=pt.shape)
            } else {
                p=p+geom_point(data=subset(data.plot, row.names(data.plot) %in% highlight.cell), fill=highlight.col, color=outline.color,
                size=pt.size, shape=pt.shape, stroke = outline.size)
            }
        }
        
    } else {
        p=ggplot(data.plot,aes(x=x,y=y))
        if (pt.shape != "16") {
            p=p+geom_point(size=pt.size, shape=pt.shape, color=outline.color, fill=bg.color, stroke = outline.size)
            p=p+geom_point(data=subset(data.plot, data.plot$ident %in% highlight.use), aes(fill =ident), color=outline.color,
            size=pt.size, shape=pt.shape, stroke = outline.size)
            
            if (is.null(highlight.col) == FALSE) {
                p=p+scale_fill_manual(values = highlight.col)
            }
        }
        
        else{
            p=p+geom_point(size=pt.size, shape=pt.shape, color=bg.color)
            p=p+geom_point(data=subset(data.plot, data.plot$ident %in% highlight.use), aes(color =ident), size=pt.size, shape=pt.shape)
            
            if (is.null(highlight.col) == FALSE) {
                p=p+scale_colour_manual(values = highlight.col)
            }
        }
    }
    
    if (no.rect == TRUE) {
        p=p+theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white",color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.grid.major = element_line(color = NA),
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = 'none')
    }else if (theme.grey != TRUE){
        p=p+theme_bw()+theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1))+
        labs(title=title, fill = group.by, color = group.by)+
        xlab(xlab)+ylab(ylab)
        if(no.legend != TRUE){p=p+guides(fill=guide_legend(ncol=legend_ncol, byrow=TRUE))}
        else{p=p+theme(legend.position = 'none')}
    }else{
        p=p+theme_grey()+
        labs(title=title, fill = group.by, color = group.by)+
        xlab(xlab)+ylab(ylab)
        if(no.legend != TRUE){p=p+guides(fill=guide_legend(ncol=legend_ncol, byrow=TRUE))}
        else{p=p+theme(legend.position = 'none')}
    }
    
    if (do.label!=FALSE) {
        data.plot %>% dplyr::group_by(ident) %>% dplyr::summarize(x = median(x = x), y = median(x = y)) -> centers
        #p <- p + geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
        #geom_text(data = centers, mapping = aes(label = ident), size = label.size, color=label.color, fontface='bold')}
        p <- p + geom_point(data=subset(centers, centers$ident %in% highlight.use), mapping = aes(x = x, y = y), size = 0, alpha = 0) +
        geom_text(data=subset(centers, centers$ident %in% highlight.use), mapping = aes(label = ident), size = label.size, color=label.color, fontface='bold')}
    
    return(p)
}


# myComplexHeatmap source code
library(ComplexHeatmap)
library(viridis)

SetIfNull <- function(x, default) {
    if(is.null(x = x)){
        return(default)
    } else {
        return(x)
    }
}


myComplexHeatmap <- function(object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, genelist = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", group.order = NULL, heatmap.col = viridis(50), ann.colors, label.set=NULL){
    
    # Prepare data matrix
    if (is.null(x = data.use)) {
        if (use.scaled) {
            data.use <- GetAssayData(object = object, slot = "scale.data")
        } else {
            data.use <- GetAssayData(object = object)
        }
    }
    # note: data.use should have cells as column names, genes as row names
    cells.use <- SetIfNull(x = cells.use, default = colnames(x = object))
    cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
    if (length(x = cells.use) == 0) {
        stop("No cells given to cells.use present in object")
    }
    genes.use <- SetIfNull(x = genelist$gene, default = rownames(x = data.use))
    genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
    if (length(x = genes.use) == 0) {
        stop("No genes given to genes.use present in object")
    }
    
    if (is.null(x = group.by) || group.by == "ident") {
        cells.ident <- Idents(object = object)[cells.use]
    } else {
        cells.ident <- factor(x = FetchData(
        object = object,
        cells = cells.use,
        vars = group.by
        )[, 1])
        names(x = cells.ident) <- cells.use
    }
    cells.ident <- factor(
    x = cells.ident,
    labels = intersect(x = levels(x = cells.ident), y = cells.ident)
    )
    data.use <- data.use[genes.use, cells.use, drop = FALSE]
    if ((!use.scaled)) {
        data.use = as.matrix(x = data.use)
        if (disp.max==2.5) disp.max = 10;
    }
    data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
    data.use <- as.data.frame(x = t(x = data.use))
    data.use$cell <- rownames(x = data.use)
    colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
    data.use$ident <- cells.ident[data.use$cell]
    
    if(!is.null(group.order)) {
        if(length(group.order) == length(levels(data.use$ident)) && all(group.order %in% levels(data.use$ident))) {
            data.use$ident <- factor(data.use$ident, levels = group.order)
        }
        else {
            stop("Invalid group.order")
        }
    }
    
    data.use <- data.use[order(data.use$ident),]
    
    # Make color for annotation
    ann_colors = list(ident=c(ann.colors))
    names(ann_colors$ident) <- levels(data.use$ident)
    
    # Make annotation for column
    ann.col <- ann.colors
    names(ann.col) <- levels(data.use$ident)
    #annotation_col = HeatmapAnnotation(ident=anno_simple(data.use$ident, height = unit(1, "cm")), col = ann_colors, show_legend = FALSE)
    annotation_col = HeatmapAnnotation(ident = data.use$ident, col = ann_colors, show_legend = FALSE, show_annotation_name = FALSE)
    
    # Make annotation for row
    rownames(genelist) <- make.unique(genelist$gene)
    SubsetGeneList <- genelist[colnames(data.use)[1:(length(colnames(data.use))-2)],]
    SubsetGeneList$cluster <- factor(SubsetGeneList$cluster, levels = levels(data.use$ident))
    annotation_row = rowAnnotation(ident=SubsetGeneList$cluster, show_legend = FALSE, col=ann_colors, show_annotation_name = FALSE)
    
    # Make heatmap
    #colfunc<-colorRampPalette(c("black","#440154"))
    #heatmap.col <- c(colfunc(20),viridis(50))
    
    # Make legend
    #lgd = Legend(col_fun = heatmap.col, title = "", legend_height = unit(6, "cm"), at = c(-2.5, -2, -1, 0, 1, 2, 2.5))
    #breaks=bks,
    #gaps_col = cumsum(table(data.use$ident)),
    
    ph <- Heatmap(t(data.use[,1:length(genes.use)]),
    use_raster = T,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "Heatmap",
    column_split = data.use$ident,
    column_gap = unit(0.3, "mm"),
    row_split = SubsetGeneList$cluster,
    row_gap = unit(0.3, "mm"),
    heatmap_legend_param = list(title = "",
    legend_height = unit(6, "cm"), at = c(disp.min, -2, -1, 0, 1, 2, disp.max)),
    top_annotation = annotation_col,
    left_annotation = annotation_row,
    col = heatmap.col)
    
    if (!is.null(label.set)) {
        position <- which(is.element(rownames(t(data.use)),label.set)==TRUE)
        label.set2 <- rownames(t(data.use))[position]
        ph+rowAnnotation(link = anno_mark(at = position, labels = label.set2),
        width = unit(1, "cm")+ max_text_width(label.set2,gp = gpar(fontsize = 20, fontface="italic")))
    } else {
        return(ph)
    }
}
