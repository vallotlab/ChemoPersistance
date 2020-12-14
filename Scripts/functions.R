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

#' Get an Igraph from a given KEGG pathway
#'
#' @param KEGG_pathway_id 
#' @param plot 
#' @param plot_file 
#'
#' @return
#' @export
#'
#' @examples
GetKEGGigraph <- function(KEGG_pathway_id, plot = FALSE, plot_file=NULL){
    stopifnot(is.character(KEGG_pathway_id))
    
    # retrieve pathway thanks to KEGGgraph
    tmp <- tempfile()
    res = KEGGgraph::retrieveKGML(KEGG_pathway_id, organism="hsa", destfile=tmp, method="wget", quiet=TRUE)
    mapkG <- KEGGgraph::parseKGML2Graph(res,expandGenes=TRUE, genesOnly = TRUE)
    
    outs <- sapply(KEGGgraph::edges(mapkG), length) > 0
    ins <- sapply(KEGGgraph::inEdges(mapkG), length) > 0
    ios <- outs | ins
    ## translate the KEGG IDs into Gene Symbol
    if(require(org.Hs.eg.db)) {
        ioGeneID <- KEGGgraph::translateKEGGID2GeneID(names(ios))
        nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
    } else {
        nodesNames <- names(ios)
    }
    names(nodesNames) <- names(ios)
    
    mapkG_igraph = igraph::igraph.from.graphNEL(mapkG, name = TRUE, weight = TRUE,
                                        unlist.attrs = TRUE)
    mapkG_igraph = igraph::simplify(mapkG_igraph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
    Isolated = which(igraph::degree(mapkG_igraph)==0)
    mapkG_igraph = igraph::delete.vertices(mapkG_igraph, Isolated)
    
    # Minimum spanning tree graph from pathway
    mstree = igraph::mst(mapkG_igraph)
    V(mstree)$id <- seq_len(vcount(mstree))-1
    roots <- sapply(igraph::decompose(mstree), function(x) {
        V(x)$id[ igraph::topo_sort(x)[1]+1 ] })
    
    if(plot & !is.null(plot_file)){
        # Change names to gene name for graph
        mapkG@nodes = nodesNames
        names(mapkG@edgeL) = nodesNames
        
        mapkG_igraph_gene = igraph::igraph.from.graphNEL(mapkG, name = TRUE, weight = TRUE,
                                                 unlist.attrs = TRUE)
        mapkG_igraph_gene = igraph::simplify(mapkG_igraph_gene, remove.multiple = TRUE, remove.loops = TRUE,
                                edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
        Isolated = which(igraph::degree(mapkG_igraph_gene)==0)
        mapkG_igraph_gene = igraph::delete.vertices(mapkG_igraph_gene, Isolated)
        
        # Minimum spanning tree graph from pathway
        mstree_gene = igraph::mst(mapkG_igraph_gene)
        pdf(file.path(plot_file))
        plot(mstree_gene, layout = igraph::layout_nicely(mstree_gene),
             vertex.color= ifelse(,"red","grey"), vertex.size = 3.75,
             vertex.label.cex=0.25, edge.arrow.width=0.25, edge.arrow.size=0.25, edge.width=0.5)
        dev.off()

    }
    
    return(mstree)
}

#' For each gene in a given list, returns the 'upstream score'
#'
#' @param gene_list 
#'
#' @return
#' @export
#'
#' @examples
upstream_score_KEGG <- function(gene_list){
    ##Get the Entrez gene IDs associated with those symbols
    EG_IDs = mget(gene_list, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
    
    ##Then get the KEGG IDs associated with those entrez genes.
    KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA)
    
    results <- foreach::foreach(KEGG_id = names(KEGG_IDs), .combine=rbind,
                                .packages=c('KEGGREST',"KEGGgraph",'igraph')) %dopar% 
        {
        # results=data.frame("upstream_score"=0,"KEGG_id"="")
         # for(KEGG_id in names(KEGG_IDs)) {
            KEGG_pathway_id = KEGG_IDs[[KEGG_id]]
            print(KEGG_id)
            if(is.na(KEGG_pathway_id[1])) {
                ret = data.frame("upstream_score"=NA, "KEGG_id"=KEGG_id)
                return(ret)
                # results = rbind(results,ret)
            } else {
                list_topo_sorted <- lapply(KEGG_pathway_id, function(id){
                    mstree_sorted <- data.frame("upstream_score"=0,"KEGG_id"="")
                     try({
                        mstree <- GetKEGGigraph(id, plot=FALSE)
                        mstree_sorted <- igraph::as_ids(igraph::topo_sort(mstree))
                        }, TRUE)
                    
                    return(mstree_sorted)
                })
                names(list_topo_sorted) <- paste0(KEGG_id,"_",KEGG_pathway_id)
                list_topo_sorted = lapply(list_topo_sorted, function(x){
                    n = 1 - (which(x==paste0("hsa:",KEGG_id))/length(x))
                    if(length(n)==0) return(NA) else return(n[1])
                } )
                df_topo_sorted = as.data.frame(t(as.data.frame(list_topo_sorted,drop=F)))
                df_topo_sorted$KEGG_id = rep(KEGG_id,nrow(df_topo_sorted))
                colnames(df_topo_sorted)[1] = "upstream_score"
                if(!is.null(df_topo_sorted)) return(df_topo_sorted)
                # if(!is.null(df_topo_sorted)) results=rbind(results,df_topo_sorted)
            }

        }
    results$KEGG_pathway = rownames(results)
    results$KEGG_pathway[grep("_",results$KEGG_pathway,invert = T)] = NA
    results$KEGG_pathway = gsub(".*_","",results$KEGG_pathway)
    
    results$Gene = sapply(mget(results$KEGG_id, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
    results$Pathway = ""
    
    l = sapply(paste0("path:hsa",results$KEGG_pathway[which(!is.na(results$KEGG_pathway))]),function(x){
        print(x)
        ret = ""
        try({
            query = keggGet(x)
            ret = query[[1]]$PATHWAY_MAP
        }, TRUE)
        return(ret)
        })
    results$Pathway[which(!is.na(results$KEGG_pathway))] = l
    return(results)
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
