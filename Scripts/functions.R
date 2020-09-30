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
