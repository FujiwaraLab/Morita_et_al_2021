# Load required libraries

library(grid)
library(tidyr)
library(ggplot2)


# myTSNEPlot source code

myTSNEPlot <- function(object, pt.size = 2,pt.shape = 21, pt.color=NULL, outline.color = "grey30",outline.size=0.2,label.size = 8,label.color="black",viridis.option=NULL, viridis.dir=1, do.label=TRUE,group.by = "ident", no.rect = TRUE, theme.grey=TRUE, title=NULL, legend_ncol=1, xlab="tSNE_1", ylab="tSNE_2", no.legend=FALSE){
    data.plot <- data.frame(object@dr$tsne@cell.embeddings,ident=factor(object@ident))
    colnames(data.plot) <- c("x","y","ident")
    data.plot$cell <- rownames(x = data.plot)
    
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(
        object = object,
        vars.all = group.by
        )[data.plot$cell, 1])
        data.plot$ident <- ident.use
    }
    
    p=ggplot(data.plot,aes(x=x,y=y))
    
    if (pt.shape != "16") {
        p=p+geom_point(size=pt.size, pch = pt.shape, color=outline.color, stroke = outline.size, aes(fill=factor(data.plot$ident)))
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_fill_viridis_d(option=viridis.option, direction = viridis.dir)}
        if (is.null(pt.color) == FALSE) {p=p=p+scale_fill_manual(values = pt.color)}
    }
    else{
        p=p+geom_point(size=pt.size, pch = pt.shape, aes(color=factor(data.plot$ident)))
        if (is.null(viridis.option) == FALSE) {p=p=p+scale_colour_viridis_d(option=viridis.option, direction = viridis.dir)}
        if (is.null(pt.color) == FALSE) {p=p=p+scale_colour_manual(values = pt.color)}
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
    
    if (do.label==TRUE) {
        data.plot %>% dplyr::group_by(ident) %>% dplyr::summarize(x = median(x = x), y = median(x = y)) -> centers
        p <- p +
        geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
        geom_text(data = centers, mapping = aes(label = ident), size = label.size, color=label.color, fontface='bold')}
    
    return(p)
}


# MultiPlotList source code

library(grid)

MultiPlotList <- function(plots, file, cols=1, layout=NULL) {
    # Make a list from the ... arguments and plotlist
    #plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col))
        }
    }
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
    #MultiPlotList(pList, layout=matrix(1:20, nrow=num.row, ncol=nCol, byrow=T))
    #rp()
}
SingleFeaturePlot <-  function(object, features.use, pt.size, pch.use, min.exp, max.exp,hl.low.color,hl.high.color,outline.color,viridis.option,outline.size,title.size,legend.position){
    data.plot=data.frame(object@dr$tsne@cell.embeddings,
    FetchData(object = object, vars.all = features.use))
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
    plot.title=element_text(face="bold", size=title.size))
    
    if (is.null(legend.position)==TRUE){
        p=p+theme(legend.position = 'none')
    }
    else {p=p+theme(legend.position = legend.position)}
    
    p=p+labs(title = features.use)
    return(p)}

# highlightTSNEPlot source code

highlightTSNEPlot <- function(object, pt.size = 2,pt.shape = 21,bg.color="#d3d3d3",outline.color = "grey30",outline.size=0.2,label.size = 8,label.color="black",do.label=FALSE,group.by = "ident", no.rect = FALSE, theme.grey=TRUE, title=NULL, legend_ncol=1, highlight.use=NULL, highlight.cell=NULL, highlight.col=NULL, no.legend = FALSE, xlab="tSNE_1", ylab="tSNE_2"){
    data.plot <- data.frame(object@dr$tsne@cell.embeddings,ident=factor(object@ident))
    colnames(data.plot) <- c("x","y","ident")
    data.plot$cell <- rownames(x = data.plot)
    
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(
        object = object,
        vars.all = group.by
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
        p <- p +
        geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
        geom_text(data = centers, mapping = aes(label = ident), size = label.size, color=label.color, fontface='bold')}
    
    return(p)
}
