geco.hclustAnnotHeatmapPlot.withColumn <- function (
    x = NULL, hc = NULL, hc_row = NULL, hmColors = NULL, 
    anocol = NULL, anorow = NULL, xpos = c(0.375, 0.9, 0.3745, 0.885, 0.05, 0.25),
    ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
    dendro.cex = 1, xlab.cex = 0.8, 
    hmRowNames = FALSE, hmRowNames.cex = 0.5, 
    hmColNames = FALSE, hmColNames.cex = 0.5,
    hmCategNamesRows = FALSE, hmCategNamesRows.cex = 0.5) 
{
    par(fig = c(xpos[1], xpos[2], ypos[5], ypos[6]), new = FALSE, 
        mar = c(0, 0, 1.5, 0))
    plot(hc, main = "", sub = "", las = 2, cex = dendro.cex, cex.axis = dendro.cex)
    par(fig = c(xpos[3], xpos[4], ypos[3], ypos[4]), new = TRUE, 
        mar = rep(0, 4))
    geco.imageCol(anocol, xlab.cex = xlab.cex, ylab.cex = 0)
    if (hmColNames) {
        axis(side = 3, padj = 0.5, hadj = 0.5,
             lwd = 0, at = seq(0, 1, length.out = ncol(x)), 
             labels = colnames(x), las = 2, cex.axis = hmColNames.cex)
    }
    par(fig = c(xpos[3], xpos[4], ypos[1], ypos[2]), new = TRUE, 
        mar = rep(0, 4))
    image(t(x), axes = FALSE, xlab = "", ylab = "", col = hmColors)
    par(fig = c(xpos[5], xpos[6], ypos[1] - 0.015, ypos[3] + 
                    0.015), new = TRUE, mar = c(0, 0, 0, 0))
    plot(as.phylo(hc_row), main = "", sub = "", las = 2, cex = dendro.cex, cex.axis = dendro.cex)
    par(fig = c(xpos[6] + 0.015, xpos[1], ypos[1], ypos[3]), 
        new = TRUE, mar = rep(0, 4))
    geco.imageCol(t(anorow)[, dim(anorow)[1]:1], xlab.cex = 0, 
                  ylab.cex = 0)
    box()
    if (hmRowNames) {
        axis(2, hadj = 0.1, lwd = 0, at = seq(0, 1, length.out = nrow(x)), 
             labels = rownames(x), las = 1, cex.axis = hmRowNames.cex)
    }
    if (hmCategNamesRows) {
        axis(1, hadj = 0.75, lwd = 0, at = seq(0, 1, length.out = ncol(anorow)), 
             labels = colnames(anorow), las = 2, cex.axis = hmCategNamesRows.cex)
    }
    
}

enrich_for_TF_ChEA3 <- function(genes_of_interest, n_random = 100, all_genes = unique(toTable(org.Hs.eg.db::org.Hs.egSYMBOL)[,2])){
    url = "https://maayanlab.cloud/chea3/api/enrich/"
    encode = "json"
    
    list_TF_enrichment = list()
    for(i in seq_len(n_random+1)){
        if(i != 1){
            genes = sample(all_genes, length(genes_of_interest))
        } else{
            genes = genes_of_interest
        }
        
        #POST to ChEA3 server
       
        httr::set_config(config(ssl_verifypeer = 0L))
        url = "https://maayanlab.cloud/chea3/api/enrich/"
        encode = "json"
        payload = list(query_name = "myQuery", gene_set = genes)
        
        #POST to ChEA3 server
        response = POST(url = url, body = payload, encode = encode)
        json = content(response, "text")
        
        #results as list of R dataframes
        results = fromJSON(json)
        
        #results as list of R dataframes
        list_TF_enrichment[[i]] = results$`Integrated--meanRank`
        
        if(i %% 10 == 0){cat("Done - ", i ," / ", n_random, ".\n")}
    }
    names(list_TF_enrichment) = c("Genes_of_interest", paste0("random_genes_",seq_len(n_random)))
    return(list_TF_enrichment)
}


TF.IDF.custom <- function(data, scale = 10000, log = TRUE, verbose = TRUE) {
    if (class(x = data) == "data.frame") {
        data <- as.matrix(x = data)
    }
    if (class(x = data) != "dgCMatrix") {
        data <- as(object = data, Class = "dgCMatrix")
    }
    if (verbose) {
        message("Performing TF-IDF normalization")
    }
    npeaks <- Matrix::colSums(x = data)
    tf <- t(x = t(x = data) / npeaks)
    # log transformation
    idf <- 1+ ncol(x = data) / Matrix::rowSums(x = data)
    norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
    if(log) norm.data = log1p(norm.data * scale) else norm.data = 
        norm.data * scale
    return(norm.data)
}


dice. = function (x, y) {
    M.11 = sum(x == 1 & y == 1)
    M.10 = sum(x == 1 & y == 0)
    M.01 = sum(x == 0 & y == 1)
    return (2*M.11 / (2*M.11 + M.10 + M.01))
}

dice <- function(m){
    dice_mat = matrix(0,ncol = nrow(m), nrow= nrow(m),
                      dimnames = list(rownames(m),
                                      rownames(m))) 
    for(row in 1:nrow(m)){
        for(col in 1:nrow(m)){
            dice_mat[row,col] = dice.(m[row,],m[col,])
        }
    }
    return(dice_mat)
}

jaccard. = function (x, y) {
    M.11 = sum(x == 1 & y == 1)
    M.10 = sum(x == 1 & y == 0)
    M.01 = sum(x == 0 & y == 1)
    return (M.11 / (M.11 + M.10 + M.01))
}

jaccard <- function(m){
    jac_mat = matrix(0,ncol = nrow(m), nrow= nrow(m),
                     dimnames = list(rownames(m),
                                     rownames(m))) 
    for(row in 1:nrow(m)){
        for(col in 1:nrow(m)){
            jac_mat[row,col] = jaccard.(m[row,],m[col,])
        }
    }
    return(jac_mat)
}

find_clusters_louvain_scExp <- function(scExp, k =100){
    g = scran::buildSNNGraph(scExp, k=k, use.dimred = 'UMAP')
    clust <- igraph::cluster_louvain(g)$membership
    clust
    cell_clusters = paste0("C",clust)
    cell_clusters = as.factor(cell_clusters)
    return(cell_clusters)
}

find_clusters_louvain = function(reduced_dim, k =100){
    ncells <- nrow(reduced_dim)
    u <- matrix(rpois(100*ncells, 5), ncol=ncells)
    sce <- SingleCellExperiment(assays=list(counts=u),
                                reducedDims=SimpleList(UMAP=as.matrix(reduced_dim) ))
    g = buildSNNGraph(sce, k=k, use.dimred = 'UMAP')
    clust <- igraph::cluster_louvain(g)$membership
    clust
    cell_clusters = paste0("C",clust)
    cell_clusters = as.factor(cell_clusters)
    return(cell_clusters)
}


import_scExp_gz <- function(datadir, pattern){
    files = list.files(datadir, pattern = pattern, full.names = TRUE)
    dir.create(file.path(datadir,"tmp"))
    files = sapply(files, function(f){
        tmp = file.path(datadir,"tmp",paste0(gsub(pattern,"",basename(f)),".tsv"))
        gunzip(f, tmp, overwrite = TRUE, remove = FALSE)
    })
    out <- ChromSCape::import_scExp(files)
    unlink(file.path(datadir,"tmp"),recursive = TRUE,force =  TRUE)
    return(out)
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

makeAttr <- function(graph, default, valNodeList) {
    tmp <- nodes(graph)
    x <- rep(default, length(tmp)); names(x) <- tmp
    
    if(!missing(valNodeList)) {
        stopifnot(is.list(valNodeList))
        allnodes <- unlist(valNodeList)
        stopifnot(all(allnodes %in% tmp))
        for(i in seq(valNodeList)) {
            x[valNodeList[[i]]] <- names(valNodeList)[i]
        }
    }
    return(x)
}

show_in_excel <- function(.data){
    tmp <- paste0(tempfile(),".csv")
    write.csv(.data, tmp)
    fs::file_show(path = tmp)
}

VennDiagram_3 <- function(genesSet1, genesSet2,genesSet3,  savePNG=F,png_title,title,groups) {
    genesSet1  = as.vector(genesSet1)
    genesSet2  = as.vector(genesSet2)
    genesSet3  = as.vector(genesSet3)
    
    aera1 = length(genesSet1)
    aera2 = length(genesSet2)
    aera3 = length(genesSet3)
    
    n12 = length(intersect(genesSet1,genesSet2))
    n23 = length(intersect(genesSet2,genesSet3))
    n13 = length(intersect(genesSet1,genesSet3))
    n123 = length(intersect(intersect(genesSet1,genesSet2),genesSet3))
    
    if(savePNG == TRUE){
        png(paste(png_title,title,".png",sep=""))
        grid.newpage()
        draw.triple.venn(area1 = aera1, area2 = aera2, area3 = aera3, n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = groups, 
                         fill = c("light green", "red", "orange"),euler.d=T, scaled=T)
        dev.off()
    }
    else{
        grid.newpage()
        draw.triple.venn(area1 = aera1, area2 = aera2, area3 = aera3, n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = groups, 
                         fill = c("light green", "red", "orange"),euler.d = T, scaled=T)
    }
}

VennDiagram_2 <- function(genesSet1, genesSet2, savePNG=F,png_title="VennDiagram",title="Main",groups) {
    genesSet1  = as.vector(genesSet1)
    genesSet2  = as.vector(genesSet2)
    aera1 = length(genesSet1)
    aera2 = length(genesSet2)
    
    
    cross.area = length(intersect(genesSet1,genesSet2))
    
    if(savePNG == TRUE){
        png(paste(png_title,title,".png",sep=""))
        
        grid.newpage()
        draw.pairwise.venn(area1 = aera1, area2 = aera2, cross.area = cross.area, category = groups,
                           , fill = c("light green", "red"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
        dev.off()
    }
    else{
        grid.newpage()
        draw.pairwise.venn(area1 = aera1, area2 = aera2, cross.area = cross.area, category =groups,
                           , fill = c("light green", "red"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
    }
}

violin_plot_cooccurence <- function(
    mat_init,mat_pers = NULL,
    sub_cells = list(colnames(mat_pers),colnames(mat_init)),
    levels=c("Untreated","Persister"),
    colors =c("#afafafff", "#1AB8AD")){
    
    if(!is.null(mat_pers)) mat = Matrix::cBind(mat_pers,mat_init)
    else mat = mat_init
    
    # mat = mat[which(rownames(mat) %in% overexpressed_MM468_persister_regions),]
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
