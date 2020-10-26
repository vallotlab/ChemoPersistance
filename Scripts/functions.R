violin_plot_cooccurence <- function(
    mat_init,mat_pers = NULL, sub_regions=rownames(mat_pers),
    sub_cells = list(colnames(mat_pers),colnames(mat_init)),
    levels=c("Untreated","Persister"),
    colors = c("#dfdfdfff","#118675ff")){
    
    if(!is.null(mat_pers)) mat = Matrix::cBind(mat_pers,mat_init)
    else mat = mat_init
    
    # mat = mat[which(rownames(mat) %in% overexpressed_MM468_persister_regions),]
    mat = mat[which(rownames(mat) %in% sub_regions),]
    mat = mat[,which(colnames(mat) %in% unlist(sub_cells))]
    bin_mat = mat
    bin_mat[(bin_mat>0)]=1
    
    coocurrence_persister_genes_score = as.data.frame(Matrix::colSums(bin_mat) / nrow(bin_mat))
    colnames(coocurrence_persister_genes_score) = "coocurrence_score"
    rownames(coocurrence_persister_genes_score) = colnames(bin_mat)
    coocurrence_persister_genes_score$sample = ""
    coocurrence_persister_genes_score[sub_cells[[1]],"sample"] = levels[1]
    if(!is.null(mat_pers)) coocurrence_persister_genes_score[sub_cells[[2]],"sample"] = levels[2]
    
    coocurrence_persister_genes_score$sample = factor(
        coocurrence_persister_genes_score$sample,levels=levels)
    
    p = ggplot(coocurrence_persister_genes_score, aes(x=sample,y=coocurrence_score,
                                                  fill=sample)) + 
        geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=colors) +
        stat_summary(fun=median, geom="point", size=2, color="black") + geom_jitter(size=0.2, alpha=0.2) 
    if(!is.null(mat_pers)) p = p + ggpubr::stat_compare_means(method = "t.test",ref.group = levels[1])
    p
}

plot_reduced_dim_scExp_devel <- function(scExp, color_by = "sample_id", reduced_dim = c("PCA", 
                                                                                  "TSNE", "UMAP"),
                                   select_x = "Component_1",
                                   select_y = "Component_2",
                                   downsample = 5000,
                                   transparency = 0.6,
                                   size = 1)
{
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(color_by), is.character(reduced_dim), 
              is.character(select_x), is.character(select_y), is.numeric(5000),
              is.numeric(transparency))
    
    if (!reduced_dim[1] %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - ", reduced_dim[1], " is not present in object, please run normalize_scExp first."))
    
    if (!color_by %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by must be present in colnames of colData(scExp).")
    
    if (!paste0(color_by, "_color") %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by's color column must be present 
         in colnames of colData(scExp). Please run color_scExp first.")
    
    if (!select_x %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_x must be present in colnames of PCA of scExp.")
    
    if (!select_y %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_y must be present in colnames of PCA of scExp.")
    
    set.seed(2047)
    if(ncol(scExp) > downsample ) scExp = scExp[,sample(ncol(scExp),downsample,replace = F)]
    
    plot_df = as.data.frame(cbind(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]), 
                                  SingleCellExperiment::colData(scExp)))
    
    p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + geom_point(alpha = transparency,
                                                                              size = size, aes(color = SingleCellExperiment::colData(scExp)[, color_by])) +
        labs(color = color_by) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA))
    
    if (color_by == "total_counts")
    {
        p <- p + scale_color_gradientn(colours = matlab.like(100))
    } else
    {
        
        cols = unique(as.character(
            SingleCellExperiment::colData(scExp)[,paste0(color_by, "_color")]))
        names(cols) = unique(as.character(
            SingleCellExperiment::colData(scExp)[,color_by]))
        p <- p + scale_color_manual(values = cols)
    }
    return(p)
}
