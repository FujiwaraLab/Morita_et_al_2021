###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("taRifx")
library("matrixStats")
library("ggplot2")
library("Rtsne")
library("fpc")
library("factoextra")
library("monocle")
library("viridis")
library("gplots")
library("RColorBrewer")
library("destiny")
library("slingshot")
library("rgl")
library("scatterplot3d")
library("made4")
library("pheatmap")
library("matrixStats")
library("statmod")
library("FactoMineR")
library("jackstraw")
library("ReactomePA")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")
library("ggpubr")


############
isSparseMatrix <- function(x){
    class(x) %in% c("dgCMatrix", "dgTMatrix")
}

############
smartEsApply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
    parent <- environment(FUN)
    if (is.null(parent))
    parent <- emptyenv()
    e1 <- new.env(parent=parent)
    multiassign(names(pData(X)), pData(X), envir=e1)
    environment(FUN) <- e1
    
    if (isSparseMatrix(exprs(X))){
        res <- sparseApply(exprs(X), MARGIN, FUN, convert_to_dense, ...)
    }else{
        res <- apply(exprs(X), MARGIN, FUN, ...)
    }
    
    if (MARGIN == 1)
    {
        names(res) <- row.names(X)
    }else{
        names(res) <- colnames(X)
    }
    
    res
}

############
calculate_NB_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
{
    expr_hint <- expr_selection_func(f_expression)
    if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
        disp_guess_fit <- disp_func(expr_hint)
        
        # For NB: Var(Y)=mu*(1+mu/k)
        f_expression_var <- var(f_expression)
        f_expression_mean <- mean(f_expression)
        
        disp_guess_meth_moments <- f_expression_var - f_expression_mean
        disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
        
        #return (max(disp_guess_fit, disp_guess_meth_moments))
        return (disp_guess_fit)
    }
    return (NULL)
}

############
sparseApply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
    if (convert_to_dense){
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }
    }else{
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }
    }
    
    return(res)
    
}

############
fit_model_helper <- function(x,
modelFormulaStr,
expressionFamily,
relative_expr,
disp_func=NULL,
verbose=FALSE,
...){
    
    modelFormulaStr <- paste("f_expression", modelFormulaStr,
    sep = "")
    orig_x <- x
    # FIXME: should we be using this here?
    # x <- x + pseudocount
    if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        if (relative_expr) {
            x <- x/Size_Factor
        }
        f_expression <- round(x)
        if (is.null(disp_func) == FALSE) {
            disp_guess <- calculate_NB_dispersion_hint(disp_func,
            round(orig_x))
            if (is.null(disp_guess) == FALSE && disp_guess >
            0 && is.na(disp_guess) == FALSE) {
                size_guess <- 1/disp_guess
                if (expressionFamily@vfamily == "negbinomial")
                expressionFamily <- negbinomial(isize=1/disp_guess, ...)
                else
                expressionFamily <- negbinomial.size(size=1/disp_guess, ...)
            }
        }
    }
    else if (expressionFamily@vfamily %in% c("uninormal", "binomialff")) {
        f_expression <- x
    }
    else {
        f_expression <- log10(x)
    }
    tryCatch({
        if (verbose) {
            FM_fit <- VGAM::vglm(as.formula(modelFormulaStr),
            family = expressionFamily, epsilon=1e-1)
        }
        else {
            FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr),
            family = expressionFamily, epsilon=1e-1))
        }
        FM_fit
    }, error = function(e) {
        print (e);
        # If we threw an exception, re-try with a simpler model.  Which one depends on
        # what the user has specified for expression family
        #print(disp_guess)
        backup_expression_family <- NULL
        if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
            disp_guess <- calculate_NB_dispersion_hint(disp_func, round(orig_x), expr_selection_func = max)
            backup_expression_family <- negbinomial()
        }else if (expressionFamily@vfamily %in% c("uninormal")){
            backup_expression_family <- NULL
        }else if (expressionFamily@vfamily %in% c("binomialff")){
            backup_expression_family <- NULL
        }else{
            backup_expression_family <- NULL
        }
        if (is.null(backup_expression_family) == FALSE){
            
            #FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, trace=T, epsilon=1e-1, checkwz=F)
            test_res <- tryCatch({
                if (verbose){
                    FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, epsilon = 1e-1, checkwz= TRUE)
                }else{
                    FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, epsilon = 1e-1, checkwz = TRUE))
                }
                FM_fit
            },
            #warning = function(w) { FM_fit },
            error = function(e) {
                #print (e);
                NULL
            })
            #print(test_res)
            test_res
        } else {
            #print(e);
            NULL
        }
    })
}


##########
buildBranchCellDataSet.mod <- function(cds,
pseudotime,
lineageA=2,
lineageB=3,
progenitor_method = c('sequential_split', 'duplicate'))
{
    lineage1 <- sort(pseudotime[!is.na(pseudotime[,lineageA]),lineageA])
    lineage2 <- sort(pseudotime[!is.na(pseudotime[,lineageB]),lineageB])
    
    ordered_lineage1 <- lineage1[order(lineage1, decreasing = FALSE)]
    ordered_lineage2 <- lineage2[order(lineage2, decreasing = FALSE)]
    
    lineage1_cells <- names(ordered_lineage1)
    lineage2_cells <- names(ordered_lineage2)
    
    ancestor_cells_for_branch <- Reduce(intersect, list(lineage1_cells,lineage2_cells))
    
    branch_cell <- names(rev(ordered_lineage1[ancestor_cells_for_branch])[1])
    
    descendents1 <- names(which( ordered_lineage1 > max(ordered_lineage1[ancestor_cells_for_branch])))
    descendents2 <- names(which( ordered_lineage2 > max(ordered_lineage2[ancestor_cells_for_branch])))
    
    path_to_root1 <- unique(c(ancestor_cells_for_branch, branch_cell, descendents1))
    path_to_root2 <- unique(c(ancestor_cells_for_branch, branch_cell, descendents2))
    #paths_to_root <- list(ordered_lineage1,ordered_lineage2)
    
    all_cells_in_subset <- unique(c(path_to_root1,path_to_root2))
    
    common_ancestor_cells <- intersect(path_to_root1, path_to_root2)
    # if (length(paths_to_root) > 2){
    #   for (i in seq(3,length(paths_to_root))){
    #     common_ancestor_cells <- intersect(intersect(paths_to_root[i], paths_to_root[i-1]), common_ancestor_cells)
    #   }
    # }
    
    #when n-center used, this creates problems
    cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, branch_cells)
    
    #State <- pData(cds)$State
    Pseudotime <- pData(cds)$Pseudotime
    
    cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))]
    
    pData <- pData(cds)
    
    pData$Pseudotime <- rep(-1, length(rownames(pData)))
    
    ###**********
    
    pseudotime_path1 <- 0:(length(path_to_root1)-1)
    names(pseudotime_path1) <- path_to_root1
    pseudotime_path1 <- 100 * pseudotime_path1 / max(pseudotime_path1)
    
    pseudotime_path2 <- 0:(length(path_to_root2)-1)
    names(pseudotime_path2) <- path_to_root2
    branch_pseudotime <- pseudotime_path2[branch_cell][[1]]
    pseudotime_path2[common_ancestor_cells] <- pseudotime_path1[common_ancestor_cells]
    pseudotime_path2_descendents <- pseudotime_path2[which(!names(pseudotime_path2) %in% common_ancestor_cells)]
    
    scale.factor <- (100-pseudotime_path2[branch_cell][[1]])/length(pseudotime_path2_descendents)
    
    descendent_modified_pseudo <- (pseudotime_path2_descendents-(length(common_ancestor_cells)-1)) * scale.factor + pseudotime_path2[branch_cell][[1]]
    pseudotime_path2[which(!names(pseudotime_path2) %in% common_ancestor_cells)] <- descendent_modified_pseudo
    
    total.pseudotime <- c(pseudotime_path1,pseudotime_path2)
    total.pseudotime <- total.pseudotime[all_cells_in_subset]
    
    pData$Pseudotime <- total.pseudotime[rownames(pData)]
    
    ###**********
    
    pData$original_cell_id <- row.names(pData)
    
    pData[common_ancestor_cells, "State"] <- 1
    pData[descendents1, "State"] <- 2
    pData[descendents2, "State"] <- 3
    pData[, "State"] <- as.factor(pData[, "State"])
    
    progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 'Pseudotime'])
    
    paths_to_root <- list(ordered_lineage1[path_to_root1],ordered_lineage2[path_to_root2])

    
    if (progenitor_method == 'duplicate') {
        ancestor_exprs <- exprs(cds)[,common_ancestor_cells]
        expr_blocks <- list()
        
        # Duplicate the expression data
        for (i in 1:length(paths_to_root)) { #duplicate progenitors for multiple branches
            if (nrow(ancestor_exprs) == 1)
            exprs_data <- t(as.matrix(ancestor_exprs))
            else exprs_data <- ancestor_exprs
            
            colnames(exprs_data) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
            expr_lineage_data <- exprs(cds)[,setdiff(names(paths_to_root[[i]]), common_ancestor_cells)]
            exprs_data <- cbind(exprs_data, expr_lineage_data)
            expr_blocks[[i]] <- exprs_data
        }
        
        # Make a bunch of copies of the pData entries from the common ancestors
        ancestor_pData_block <- pData[common_ancestor_cells,]
        
        pData_blocks <- list()
        
        weight_vec <- c()
        for (i in 1:length(paths_to_root)) {
            weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
            weight_vec_block <- rep(1, length(common_ancestor_cells))
            
            #pData <- rbind(pData, pData[common_ancestor_cells, ])
            new_pData_block <- ancestor_pData_block
            # new_pData_block$Lineage <- lineage_states[i]
            # new_pData_block$State <- lineage_states[i]
            
            row.names(new_pData_block) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
            
            pData_lineage_cells <- pData[setdiff(names(paths_to_root[[i]]), common_ancestor_cells),]
            # pData_lineage_cells$Lineage <- lineage_states[i]
            # pData_lineage_cells$State <- lineage_states[i]
            
            weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
            
            weight_vec <- c(weight_vec, weight_vec_block)
            
            new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
            new_pData_block$Branch <- paste('Branch', i, sep = '_')
            pData_blocks[[i]] <- new_pData_block
        }
        pData <- do.call(rbind, pData_blocks)
        exprs_data <- do.call(cbind, expr_blocks)
    }
    else if(progenitor_method == 'sequential_split') {
        pData$Branch <- names(paths_to_root)[1]
        
        branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
        pData[common_ancestor_cells[branchA], 'Branch'] <- names(paths_to_root)[1]
        branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
        pData[common_ancestor_cells[branchB], 'Branch'] <- names(paths_to_root)[2]
        
        # Duplicate the root cell to make sure both regression start at pseudotime zero:
        zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
        exprs_data <- cbind(exprs(cds), 'duplicate_root' = exprs(cds)[, zero_pseudotime_root_cell])
        pData <- rbind(pData, pData[zero_pseudotime_root_cell, ])
        row.names(pData)[nrow(pData)] <- 'duplicate_root'
        pData[nrow(pData), 'Branch'] <- names(paths_to_root)[2]
        
        weight_vec <- rep(1, nrow(pData))
        
        for (i in 1:length(paths_to_root)){
            path_to_ancestor <- paths_to_root[[i]]
            branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
            pData[branch_cells,]$Branch <- names(paths_to_root)[i]
        }
    }
    
    pData$Branch <- as.factor(pData$Branch)
    
    pData$State <- factor(pData$State)
    Size_Factor <- pData$Size_Factor
    
    fData <- fData(cds)
    
    colnames(exprs_data) <- row.names(pData) #check this
    cds_subset <- newCellDataSet(as.matrix(exprs_data),
    phenoData = new("AnnotatedDataFrame", data = pData),
    featureData = new("AnnotatedDataFrame", data = fData),
    expressionFamily=cds@expressionFamily,
    lowerDetectionLimit=cds@lowerDetectionLimit)
    pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
    pData(cds_subset)$Size_Factor <- Size_Factor
    
    cds_subset@dispFitInfo <- cds@dispFitInfo
    
    return (cds_subset)
}




#########
##cds: object returned by buildBranchCellDataSet.mod
BEAM.mod <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
relative_expr = TRUE,
cores = 1,
verbose = FALSE,
...) {

    cds_subset <- cds
    
    branchTest_res <- differentialGeneTest(cds_subset,
    fullModelFormulaStr = fullModelFormulaStr,
    reducedModelFormulaStr = reducedModelFormulaStr,
    cores = cores,
    relative_expr = relative_expr,
    verbose=verbose)
    
    cmbn_df <- branchTest_res[, 1:4]
    fd <- fData(cds)[row.names(cmbn_df),]
    cmbn_df <- cbind(cmbn_df, fd)
    
    return(cmbn_df)
    
}


############
##cds: object returned by buildBranchCellDataSet.mod
plot_genes_branched_heatmap.mod <- function(cds_subset,
branch_states=NULL,
branch_labels = c("Cell fate 1", "Cell fate 2"),
cluster_rows = TRUE,
hclust_method = "ward.D2",
num_clusters = 6,
hmcols = NULL,
branch_colors = c('#979797', '#F05662', '#7990C8'),
add_annotation_row = NULL,
add_annotation_col = NULL,
show_rownames = FALSE,
use_gene_short_name = TRUE,
scale_max=3,
scale_min=-3,
norm_method = c("log", "vstExprs"),

trend_formula = '~sm.ns(Pseudotime, df=3) * Branch',

return_heatmap=FALSE,
cores = 1, ...) {
    
    new_cds <- cds_subset
    
    if(is.null(branch_states)) {
        progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, 'State']
        branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
    }
    
    col_gap_ind <- 101
    # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
    # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
    
    newdataA <- data.frame(Pseudotime = seq(0, 100,
    length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
    newdataB <- data.frame(Pseudotime = seq(0, 100,
    length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
    
    BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,
    relative_expr = T, new_data = rbind(newdataA, newdataB))
    
    BranchA_exprs <- BranchAB_exprs[, 1:100]
    BranchB_exprs <- BranchAB_exprs[, 101:200]
    
    #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
    common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == setdiff(pData(new_cds)$State, branch_states),])
    BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
    BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
    BranchB_num <- BranchA_num
    
    norm_method <- match.arg(norm_method)
    
    # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
    if(norm_method == 'vstExprs') {
        BranchA_exprs <- vstExprs(new_cds, expr_matrix=BranchA_exprs)
        BranchB_exprs <- vstExprs(new_cds, expr_matrix=BranchB_exprs)
    } else if(norm_method == 'log') {
        BranchA_exprs <- log10(BranchA_exprs + 1)
        BranchB_exprs <- log10(BranchB_exprs + 1)
    }
    
    heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
    
    heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
    heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
    heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
    heatmap_matrix[is.nan(heatmap_matrix)] = 0
    heatmap_matrix[heatmap_matrix>scale_max] = scale_max
    heatmap_matrix[heatmap_matrix<scale_min] = scale_min
    
    heatmap_matrix_ori <- heatmap_matrix
    heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch
    
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    
    exp_rng <- range(heatmap_matrix) #bks is based on the expression range
    bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
    if(is.null(hmcols)) {
        hmcols <- blue2green2red(length(bks) - 1)
    }
    
    # prin  t(hmcols)
    ph <- pheatmap(heatmap_matrix,
    useRaster = T,
    cluster_cols=FALSE,
    cluster_rows=TRUE,
    show_rownames=F,
    show_colnames=F,
    #scale="row",
    clustering_distance_rows=row_dist,
    clustering_method = hclust_method,
    cutree_rows=num_clusters,
    silent=TRUE,
    filename=NA,
    breaks=bks,
    color=hmcols
    #color=hmcols#,
    # filename="expression_pseudotime_pheatmap.pdf",
    )
    #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, branch_labels, BranchA_num, BranchP_num, BranchB_num, file = 'heatmap_matrix')
    
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
    
    if(!is.null(add_annotation_row)) {
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
        # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
    }
    
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(branch_labels[1], BranchA_num),
    rep("Pre-branch",  2 * BranchP_num),
    rep(branch_labels[2], BranchB_num)))
    
    colnames(annotation_col) <- "Cell Type"
    
    if(!is.null(add_annotation_col)) {
        annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])
    }
    
    names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
    
    annotation_colors=list("Cell Type"=branch_colors)
    
    names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
    
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            
            row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        } else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    } else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }
    
    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) <- row_ann_labels
    
    ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
    useRaster = T,
    cluster_cols=FALSE,
    cluster_rows=TRUE,
    show_rownames=show_rownames,
    show_colnames=F,
    #scale="row",
    clustering_distance_rows=row_dist, #row_dist
    clustering_method = hclust_method, #ward.D2
    cutree_rows=num_clusters,
    # cutree_cols = 2,
    annotation_row=annotation_row,
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    gaps_col = col_gap_ind,
    treeheight_row = 20,
    breaks=bks,
    fontsize = 6,
    color=hmcols,
    border_color = NA,
    silent=FALSE)
    
    grid::grid.rect(gp=grid::gpar("fill", col=NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap){
        return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, heatmap_matrix = heatmap_matrix,
        heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols,
        annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col,
        ph_res = ph_res))
    }
    return(ph_res)
}




###############
##cds: object returned by buildBranchCellDataSet.mod
plot_genes_branched_pseudotime.mod <- function (cds,
branch_states = NULL,
branch_point=1,
branch_labels = NULL,
method = "fitting",
min_expr = NULL,
cell_size = 0.75,
nrow = NULL,
ncol = 1,
panel_order = NULL,
color_by = "State",
expression_curve_linetype_by = "Branch",
trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch",
reducedModelFormulaStr = NULL,
label_by_short_name = TRUE,
relative_expr = TRUE,
#gene_pairs = NULL,
...)
{
    Branch <- NA
    if (is.null(reducedModelFormulaStr) == FALSE) {
        #pval_df <- branchTest(cds,
        #branch_states=branch_states,
        #branch_point=branch_point,
        #fullModelFormulaStr = trend_formula,
        #reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)",
        #...)
        
        cds_subset <- cds
        pval_df <- differentialGeneTest(cds_subset,
        fullModelFormulaStr = trend_formula,
        reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)",
        cores = 1,
        relative_expr = TRUE,
        verbose=FALSE)
        fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
    }
    if("Branch" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Branch is in the model formula we will duplicate the "progenitor" cells
        cds_subset <- cds
    }
    else {
        cds_subset <- cds
        pData(cds_subset)$Branch <- pData(cds_subset)$State
    }
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
    }
    if (integer_expression) {
        CM <- exprs(cds_subset)
        if (relative_expr){
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(CM)))
    }
    else {
        cds_exprs <- reshape2::melt(exprs(cds_subset))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- round(cds_exprs$expression)
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
    # trend_formula <- paste("adjusted_expression", trend_formula,
    #     sep = "")
    cds_exprs$Branch <- as.factor(cds_exprs$Branch)
    
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
    
    full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
    relative_expr = T, new_data = new_data)
    colnames(full_model_expectation) <- colnames(cds_subset)
    
    cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
    if(!is.null(reducedModelFormulaStr)){
        reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
        relative_expr = T, new_data = new_data)
        colnames(reduced_model_expectation) <- colnames(cds_subset)
        cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
    }
    
    # FIXME: If you want to show the bifurcation time for each gene, this function
    # should just compute it. Passing it in as a dataframe is just too complicated
    # and will be hard on the user.
    # if(!is.null(bifurcation_time)){
    #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
    # }
    if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
        levels = panel_order)
    }
    cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
    cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
    
    if(!is.null(reducedModelFormulaStr)){
        cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
        cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
    }
    
    cds_exprs$State <- as.factor(cds_exprs$State)
    cds_exprs$Branch <- as.factor(cds_exprs$Branch)
    
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    # if (!is.null(bifurcation_time)) {
    #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
    #                       color = "black", linetype = "longdash")
    # }
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
    }
    if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
    pval, nrow = nrow, ncol = ncol, scales = "free_y")
    else q <- q + scale_y_log10() + facet_wrap(~feature_label,
    nrow = nrow, ncol = ncol, scales = "free_y")
    if (method == "loess")
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
    method = "loess")
    else if (method == "fitting") {
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
        linetype = "Branch"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
    }
    
    if(!is.null(reducedModelFormulaStr)) {
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
        color = 'black', linetype = 2, data =  cds_exprs)
    }
    
    q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
    
    q <- q + monocle_theme_opts()
    q + expand_limits(y = min_expr)
}
