#' Title
#'
#' @param dist_mat 
#' @param sigma 
#' @param K 
#' @param nEV 
#' @param t 
#' @param PLOT 
#'
#' @return
#' @export
#'
#' @examples
#' 
diffusion_map <- function(dist_mat, 
                          sigma = NULL, 
                          K = 10, 
                          nEV = 30,
                          t = 1, 
                          PLOT = FALSE){
  K <- min(K, nrow(dist_mat))
  dists <- as.matrix(dist_mat)
  sigma_list <- c()
  if (is.null(sigma)){
    for (i in 1:nrow(dists)){
      sigma_list[i] <- sort(dists[i,])[K]
    }
  }
  if (!is.null(sigma)){
    sigma_list <- rep(sigma, nrow(dists))
  }
  
  affinity_matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
  for (i in 1:nrow(affinity_matrix)){
    dist_vec <- dists[i, ]
    dist_vec[is.na(dist_vec)] <- 10^6
    affinity_matrix[i, ] <- exp(-dist_vec^2/(sigma_list[i]^2))
  }
  affinity_matrix_2 <- (affinity_matrix + t(affinity_matrix))/2
  normalized_vec <- sqrt(1/apply(affinity_matrix_2, 1, sum))
  affinity_matrix_3 <- t(affinity_matrix_2 * normalized_vec) * normalized_vec
  
  N_EV <- min(nrow(affinity_matrix_3), nEV)
  E_list <- rARPACK::eigs_sym(affinity_matrix_3, k = N_EV)
  
  if (PLOT){
    plot(1:N_EV, E_list$values[1:N_EV])
  }
  
  diffu_emb <- matrix(0, nrow = nrow(E_list$vector), ncol = N_EV)
  for(i in 1:N_EV){
    diffu_emb[,i] <- E_list$vector[,i]*normalized_vec*(E_list$values[i]^t)
  }
  colnames(diffu_emb) <- paste0("EV_", 1:ncol(diffu_emb))
  rownames(diffu_emb) <- colnames(dist_mat)
  
  res <- list()
  res[["diffu_emb"]] <- diffu_emb
  res[["eigenvalues"]] <- E_list$values[1:N_EV]
  
  res
}


RunDM <- function(object, 
                  reduction = "pca",
                  dims = 1:30,
                  sigma = NULL, 
                  K = 10, 
                  n.components = 30,
                  t = 1, 
                  reduction.key = "DM_",
                  PLOT = FALSE){
  data_S <- object
  dist_mat <- as.matrix(dist(data_S@reductions[[reduction]]@cell.embeddings[, dims]))
  dm_emb <- diffusion_map(dist_mat,
                          sigma = sigma, 
                          K = K, 
                          nEV = n.components+1, #remove the first trivial EV
                          t = t, 
                          PLOT = PLOT)[["diffu_emb"]][, 2:(n.components+1)]
  colnames(dm_emb) <- paste0("DM_", 1:ncol(dm_emb))
  
  data_S[["dm"]] <- CreateDimReducObject(embeddings = dm_emb, key = reduction.key, assay = DefaultAssay(data_S))
  
  data_S                      
}

construct_cost_matrix <- function(object,
                                  cell.embedding = "dm",
                                  dims = 1:5,
                                  K = 10){
  data_S <- object
  cell_embedding <- data_S@reductions[[cell.embedding]]@cell.embeddings[, dims]
  
  message("Constructing KNN graph")
  require(FNN)
  knn_result <- get.knn(cell_embedding, k = K)
  KNN_adj_mat <- matrix(0, nrow = nrow(cell_embedding), ncol = nrow(cell_embedding))
  for (i in 1:nrow(KNN_adj_mat)){
    KNN_adj_mat[i, knn_result$nn.index[i,]] <- 1
  }
  
  message("Constructing cost matrix")
  require(igraph)
  knn_graph <- graph_from_adjacency_matrix(KNN_adj_mat, mode = "undirected")
  graph_dist_mat <- distances(knn_graph)
  #graph_dist_mat[1:15, 1:15]
  if (is.infinite(max(graph_dist_mat))) stop("Please increase K.")
  message(max(graph_dist_mat))
  message(sum(graph_dist_mat))
  
  graph_dist_mat
}


coarse_graining <- function(cell_embedding, #The cell embedding to construct KNN graph for coarse graining
                            gene_expression, #The cell-by-gene expression matrix
                            cell_embedding.visualization = NULL,
                            cost_matrix = NULL,
                            N = 1000, #The granularity
                            random.seed = 1){
  
  message("Run k-means")
  set.seed(random.seed)
  km.res <- kmeans(cell_embedding, N)
  
  message("Update matrices")
  
  #####KNN graph
  KNN_membership_mat <- matrix(0, nrow = N, ncol = nrow(cell_embedding))
  for (i in 1:ncol(KNN_membership_mat)){
    KNN_membership_mat[km.res$cluster[i], i] <- 1
  }
  KNN_membership_mat <- KNN_membership_mat/apply(KNN_membership_mat, 1, sum)
  
  #####Coarse-grain the cell embedding for visualization
  if (is.null(cell_embedding.visualization)) cell_embedding.visualization <- cell_embedding
  cell_embedding.visualization_updated <- KNN_membership_mat %*% cell_embedding.visualization
  
  #####Coarse-grain the gene expression matrix
  gene_expression_updated <- biclust::binarize(KNN_membership_mat, 0) %*% gene_expression
  
  #####Coarse-grain the cost matrix for the calculation of EMD
  if(is.null(cost_matrix)) cost_matrix <- as.matrix(dist(cell_embedding))
  cost_matrix_updated <- KNN_membership_mat %*% cost_matrix %*% t(KNN_membership_mat)
  
  #####Collect the output 
  result <- list()
  result[["cell_embedding.visualization"]] <- cell_embedding.visualization_updated
  result[["gene_expression"]] <- gene_expression_updated
  result[["cost_matrix"]] <- cost_matrix_updated
  result[["genes"]] <- colnames(gene_expression_updated)
  
  result
}


coarse_grain <- function(object,
                         cost.matrix,
                         genes,
                         N = 500,
                         assay = "RNA",
                         cell.embedding = "dm",
                         dims = 1:5){
  data_S <- object
  cell_embedding <- data_S@reductions[[cell.embedding]]@cell.embeddings[, dims]
  gene_expression <- t(as.matrix(data_S@assays[[assay]]@data[genes, ]))
  cg_output <- coarse_graining(cell_embedding = cell_embedding,
                               gene_expression = gene_expression,
                               cell_embedding.visualization = cell_embedding,
                               cost_matrix = cost.matrix,
                               N = N)
  
  cg_output
}


read_emd_mat <- function(dir.path, 
                         file_name = "emd.csv"){
  emd_mat <- as.matrix(read.csv(paste0(dir.path, file_name), header = F))
  genes <- read.csv(paste0(dir.path, "gene_names.csv"), header = F)[,1]
  rownames(emd_mat) <- genes
  colnames(emd_mat) <- genes
  
  emd_mat
}

get_gene_embedding <- function(dist_mat, 
                               sigma = NULL, 
                               K = 10, 
                               nEV = 30,
                               t = 1, 
                               PLOT = FALSE){
  dm_res <- diffusion_map(dist_mat, 
                          sigma, 
                          K, 
                          nEV+1,
                          t, 
                          PLOT)
  dm_res[["diffu_emb"]] <- dm_res[["diffu_emb"]][, 2:(nEV+1)]
  colnames(dm_res[["diffu_emb"]]) <- paste0("DM_", 1:ncol(dm_res[["diffu_emb"]]))
  dm_res[["eigenvalues"]] <- dm_res[["eigenvalues"]][2:(nEV+1)]
  dm_res
}


get_diffusion_matrix <- function(dist_mat,
                                 K = 10){
  K <- 10
  dists <- as.matrix(dist_mat)
  sigma_list <- c()
  for (i in 1:nrow(dists)){
    sigma_list[i] <- sort(dists[i,])[K]
  }
  
  affinity_matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
  for (i in 1:nrow(affinity_matrix)){
    dist_vec <- dists[i, ]
    dist_vec[is.na(dist_vec)] <- 10^6
    affinity_matrix[i, ] <- exp(-dist_vec^2/(sigma_list[i]^2))
  }
  affinity_matrix_2 <- (affinity_matrix + t(affinity_matrix))/2
  normalized_vec <- 1/apply(affinity_matrix_2, 1, sum)
  diffusion_mat <- affinity_matrix_2 * normalized_vec
  #print(apply(diffusion_mat, 1, sum)[1:10])
  diffusion_mat
}


construct_pseudoorder <- function(dist_mat, 
                                  subset,
                                  max.id = NULL){
  emd <- dist_mat[subset, subset]
  dm_emb <- RunDM_emd_v2(emd, PLOT = F)
  
  pseudoorder <- rank(dm_emb[,2])
  names(pseudoorder) <- rownames(dm_emb)
  if (!is.null(max.id)){
    #if(names(pseudoorder)[which.max(pseudoorder)] != max.id){
    if (pseudoorder[max.id] <= max(pseudoorder)/2){
      pseudoorder <- rank(-dm_emb[,2])
      names(pseudoorder) <- rownames(dm_emb)
    }
  }
  
  pseudoorder_all <- rep(0, nrow(dist_mat))
  names(pseudoorder_all) <- rownames(dist_mat)
  pseudoorder_all[names(pseudoorder)] <- pseudoorder
  
  pseudoorder_all
}


identify_trajectories <- function(gene_embedding,
                                  dist_mat,
                                  dims = 1:5,
                                  N,
                                  t = NULL,
                                  cutoff.list = NULL){
  dist_to_origin <- sqrt(apply(gene_embedding[, dims]^2, 1, sum))
  df_ggplot <- as.data.frame(gene_embedding[, dims])
  df_ggplot$selected <- "Other"
  nBranch <- N
  genes <- rownames(gene_embedding)
  diffusion_mat <- get_diffusion_matrix(dist_mat)
    
  for (i in 1:nBranch){
    if (length(which(df_ggplot$selected == "Other")) == 0) stop("Wrong!")
    
    if (length(which(df_ggplot$selected != "Other")) != 0) dist_to_origin[which(df_ggplot$selected != "Other")] <- -Inf
    
    message(genes[which.max(dist_to_origin)])
    
    seed <- rep(0, nrow(gene_embedding))
    seed[which.max(dist_to_origin)] <- 1
    
    seed.diffused <- seed
    max_T <- t #ceiling(nrow(dm_emb_all)/nBranch/10/5)
    for (ii in 1:max_T){
      seed.diffused <- diffusion_mat %*% matrix(seed.diffused, ncol = 1)
    }
    
    #print(seed.diffused[1:10])
    cutoff <- cutoff.list[i]
    
    df_ggplot$selected[which(-log(seed.diffused) < cutoff & df_ggplot$selected == "Other")]  <- paste0("Trajectory-", i)
    print(table(df_ggplot$selected))
    print(
      ggplot(df_ggplot, aes(x = DM_1, y = DM_2, color = selected)) + geom_point() + ggtitle("Sequential trajectory identification") #+ scale_color_gradientn(colors = c("lightgray", "red"))
    )
    df_ggplot <- cbind(df_ggplot, 
                       construct_pseudoorder(dist_mat, 
                                             which(df_ggplot$selected  == paste0("Trajectory-", i)),
                                             max.id = genes[which.max(dist_to_origin)]))
    
    colnames(df_ggplot)[ncol(df_ggplot)] <- paste0("Pseudoorder", i)
  }
  
  df_ggplot
}



get_intermediate <- function(orig_list, N){
  total_N <- length(orig_list)
  step <- floor(total_N/N)
  orig_list[c(1:N) * step]
}

add_gene_set_score <- function(object, 
                               gene_trajectory, 
                               N_bin = 5, 
                               trajectories = 1:2,
                               assay = "RNA",
                               prefix = "Trajectory",
                               use.module_score = F,
                               binarize = F,
                               reverse = NULL,
                               echo = F){
  data_S <- object
  prefix0 <- prefix
  for (trajectory in trajectories){
    prefix <- paste0(prefix0, trajectory, "_genes")
    gene_trajectory_reordered <- gene_trajectory[order(gene_trajectory[, paste0("Pseudoorder", trajectory)]),]
    genes <- rownames(gene_trajectory_reordered)[which(gene_trajectory_reordered[, paste0("Pseudoorder", trajectory)] > 0)]
    if (! is.null(reverse)){
      if (isTRUE(reverse[trajectory])) genes <- rev(genes)
    }
    length(genes)
    
    step <- length(genes)/N_bin
    metadata <- data.frame(rep(0, ncol(data_S)))
    for (i in 1:N_bin){
      start <- ceiling((i-1)*step+1)
      end <- min(ceiling(i*step), length(genes))
      #print(start:end)
      genes_subset <- genes[start:end]
      if (echo) message(paste(genes_subset, collapse = ", "))
      if (!use.module_score){
        normalized_gc <- data_S[[assay]]@data[genes_subset,]#/apply(data_S[[assay]]@data[genes_subset,],1,sum)
        #print(apply(normalized_gc,1,sum))
        if (binarize) normalized_gc <- (normalized_gc > 0)
        metadata[,i] <- apply(normalized_gc,2,sum)
        metadata[,i] <- metadata[,i]/length(genes_subset)#sum(metadata[,i])
      }
      if (use.module_score){
        feature_list <- list()
        feature_list[["tmp"]] <- genes_subset
        data_S <- AddModuleScore(data_S, features = feature_list, assay = assay, name = "tmp")
        metadata[,i] <- data_S$tmp1
      }
    }
    rownames(metadata) <- colnames(data_S)
    colnames(metadata) <- paste0(prefix, 1:N_bin)
    #print(head(metadata))
    data_S <- AddMetaData(data_S, metadata)
  }
  
  data_S
}




visualize_genes <- function(gene_embedding,
                            alpha = 0.5,
                            bty = "b2", 
                            colvar = NULL,
                            main = "trajectory", 
                            pch = 19, 
                            cex = 0.5, 
                            cex.text = 0.75,
                            theta = 0, 
                            phi = 0,
                            col = NULL,
                            genes){
  if (F){
    alpha = 0.5
    bty = "b2"
    colvar = as.integer(as.factor(gene_trajectory$selected))-1
    main = "trajectory"
    pch = 19
    cex = 0.5
    theta = 90
    phi = 30
    col = ramp.col(c("lightgray", hue_pal()(4)))
  }
  labels <- rownames(gene_embedding)
  labels[which(!labels %in% genes)] <- NA
  par(mar=c(1,1,1,1), mfrow = c(1,1))
  scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], 
            alpha = alpha,
            bty = bty, 
            colvar = colvar,
            main = main, 
            pch = pch, 
            cex = cex, 
            theta = theta, 
            phi = phi,
            col = col,
            pch = pch)
  if (F){
    scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], 
              alpha = alpha,
              bty = bty, 
              colvar = colvar,
              main = main, 
              pch = pch, 
              theta = theta, 
              phi = phi,
              col = col,
              pch = pch, cex = 0.5, add = TRUE)
    
  }
  text3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], labels = labels,
         add = TRUE, plot = TRUE, colkey = FALSE, cex = cex.text)
  
}

