generate_circular <- function(maxT = 10,
                              N_cells = 1000,
                              N_genes = 500,
                              meanlog = 0, sdlog = 0.25,
                              scale_l = 1,
                              scale=25, seed=1){
  #sigma ~ LogNormal(meanlog, sdlog^2)
  #epsilon ~ Beta(alpha, beta)
  
  N <- N_genes
  set.seed(seed)
  sigma <- scale_l * rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  set.seed(seed+1)
  alpha <- scale * rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  print(alpha)
  
  set.seed(seed+10+1)
  cell_pt <- sort(runif(n = N_cells, min = 0, max = maxT))
  set.seed(seed+10+2)
  gene_pt <- sort(runif(n = N_genes, min = 0, max = maxT))
  
  gc_mat <- matrix(0, nrow = N_genes, ncol = N_cells)
  
  for (i in 1:N){
    mean <- gene_pt[i]
    sd <- sigma[i]
    scale <- alpha[i]
    
    dist <- pmin(pmin(abs(cell_pt-mean), abs(cell_pt+maxT-mean)), abs(cell_pt-maxT-mean))
    
    expectation <- floor(scale * exp(-dist^2/(2*sd^2)))
    set.seed(seed + i)
    gc_mat[i, ] <- rpois(n = length(expectation), lambda = expectation)
    #gene_distribution_list[[i]] <- pmax(0, expectation + epsilon))
  }
  
  #gc_mat <- do.call(rbind, gene_distribution_list)
  
  rownames(gc_mat) <- unlist(gene_pt)
  colnames(gc_mat) <- unlist(cell_pt)
  gc_mat
}


generate_CC_tree <- function(maxT = 10,
                             N_cells = 15*rep(100, 7),
                             N_genes = 2*rep(100, 7),
                             N_gene_CC = 400,
                             meanlog = 0, sdlog = 0.25,
                             scale = 25,
                             seed = 1,
                             N_genes_noise = 500,
                             sigma_noise = 2.5,
                             sort = TRUE){
  cell_pt_list <- list()
  gene_pt_list <- list()
  Np <- length(N_cells)
  for (i in 1:Np){
    set.seed(seed+10*i+1)
    cell_pt_list[[i]] <- sort(runif(n = N_cells[i], min = 0, max = maxT))
    set.seed(seed+10*i+2)
    gene_pt_list[[i]] <- sort(runif(n = N_genes[i], min = 0, max = maxT))
  }
  
  cell_pt_mat <- matrix(0, nrow = Np, ncol = sum(N_cells))
  gene_pt_mat <- matrix(0, nrow = Np, ncol = sum(N_genes))
  
  cell_pt_mat[1,] <- maxT
  cell_pt_mat[1,1:N_cells[1]] <- cell_pt_list[[1]]
  gene_pt_mat[1,] <- maxT
  gene_pt_mat[1,1:N_genes[1]] <- gene_pt_list[[1]]
  
  N_cells_summed <- unlist(lapply(1:Np, FUN = function(i){
    sum(N_cells[1:i])
  }))
  N_genes_summed <- unlist(lapply(1:Np, FUN = function(i){
    sum(N_genes[1:i])
  }))
  
  for (i in 2:Np){
    cell_pt_mat[i,(N_cells_summed[i-1]+1):N_cells_summed[i]] <- cell_pt_list[[i]]
    gene_pt_mat[i,(N_genes_summed[i-1]+1):N_genes_summed[i]] <- gene_pt_list[[i]]
    if (i %in% c(4,5)){
      cell_pt_mat[2,(N_cells_summed[i-1]+1):N_cells_summed[i]] <- maxT
      gene_pt_mat[2,(N_genes_summed[i-1]+1):N_genes_summed[i]] <- maxT
    }
    if (i %in% c(6,7)){
      cell_pt_mat[3,(N_cells_summed[i-1]+1):N_cells_summed[i]] <- maxT
      gene_pt_mat[3,(N_genes_summed[i-1]+1):N_genes_summed[i]] <- maxT
    }
  }
  
  gc_mat <- matrix(0, nrow = sum(N_genes), ncol = sum(N_cells))
  
  N <- sum(N_genes)
  set.seed(seed)
  sigma <- 2*rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  set.seed(seed+1)
  alpha <- scale * rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  
  for (i in 1:N){
    mean <- gene_pt_mat[, i]
    sd <- sigma[i]
    scale <- alpha[i]
    
    dist <- unlist(lapply(1:ncol(cell_pt_mat), function(i){
      sum(abs(cell_pt_mat[, i] - mean))
    }))
    
    expectation <- floor(scale * exp(-dist^2/(2*sd^2)))
    
    set.seed(seed + i)
    gc_mat[i,] <- rpois(n = length(expectation), lambda = expectation)
    
  }
  #set.seed(1000*seed)
  #gc_mat <- rbind(gc_mat, matrix(rpois(n = N_genes_noise * ncol(gc_mat), lambda = lambda_noise), ncol = ncol(gc_mat)))
  #rownames(gc_mat) <- c(unlist(gene_pt_list), rep(-1, N_genes_noise))
  
  #set.seed(1000*seed)
  #gc_mat <- matrix(floor(pmax(0, as.numeric(gc_mat) + rnorm(n = length(as.numeric(gc_mat)), sd = sigma_noise))), nrow = nrow(gc_mat))
  
  rownames(gc_mat) <- unlist(gene_pt_list)
  colnames(gc_mat) <- unlist(cell_pt_list)
  gc_mat_tree <- gc_mat
  
  gc_mat_ring <- generate_circular(maxT = 20,
                                   N_cells = sum(N_cells),
                                   N_genes = N_gene_CC,
                                   scale_l = 2)
  
  set.seed(10000)
  gc_mat_ring <- gc_mat_ring[, sample(1:ncol(gc_mat_ring))]
  
  gc_mat <- rbind(gc_mat_tree, gc_mat_ring)
  colnames(gc_mat) <- paste0(colnames(gc_mat_tree),"|",colnames(gc_mat_ring))
  gc_mat
}


gc_mat_CC_tree <- generate_CC_tree()

dim(gc_mat_CC_tree)
max(gc_mat_CC_tree)

set.seed(1)
pheatmap::pheatmap(log(gc_mat_CC_tree[, sort(sample(1:ncol(gc_mat_CC_tree), 5000))]+1), cluster_cols = F, cluster_rows = F, color = viridis(100), show_rownames = F, show_colnames = F)

#pheatmap::pheatmap(log(gc_mat_CC_tree+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)

data_S <- CreateSeuratObject(gc_mat_CC_tree)
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)#, nfeatures = nrow(data_S))
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 100, verbose = F)
plot(data_S@reductions$pca@stdev)
DimPlot(data_S, reduction = "pca", dims = c(1,2))
data_S$cell_pseudotime_tree <- as.numeric(unlist(strsplit(colnames(data_S), split = "|", fixed = T))[(1:ncol(data_S))*2-1])
data_S$cell_pseudotime_CC <- as.numeric(unlist(strsplit(colnames(data_S), split = "|", fixed = T))[(1:ncol(data_S))*2])

data_S <- RunDM(data_S, K = 10)
DimPlot(data_S, reduction = "dm", dims = c(3,4))
par(mar = c(1,1,1,1))
scatter3D(10*data_S@reductions$dm@cell.embeddings[,21],
          10*data_S@reductions$dm@cell.embeddings[,22],
          10*data_S@reductions$dm@cell.embeddings[,23], alpha = 0.55,
          bty = "g",#colvar = as.integer(as.factor(data_S$condition)),# col = col_pal, 
          main = "trajectory", pch = 19, cex = 0.25, theta = 0+90+90, phi = 45)#, col = hue_pal()(length(unique(data_S$condition))))#ramp.col(c("red", "blue", "gray")) ) #colkey = FALSE, 

data_S <- RunTSNE(
  data_S, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/home/jz437/git/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
DimPlot(data_S, reduction = "tsne")


FeaturePlot(data_S, features = "cell_pseudotime", reduction = "tsne", order = F) & scale_color_viridis(option = "B") & NoAxes()

FeaturePlot(data_S, features = sample(rownames(data_S), size = 5), ncol = 5, pt.size = 1.5, reduction = "tsne", order = T) & scale_color_viridis(option = "D") & NoAxes()

#phate_emb <- phate(data_S@reductions$pca@cell.embeddings)

#DoHeatmap(data_S, features = sample(rownames(data_S), size = 50), size = 0)

data_S <- RunUMAP(
  data_S, dims = 1:30, n.neighbors = 30, min.dist = 0.5#, a = 5, b = 10, repulsion.strength = 5
)
DimPlot(data_S, reduction = "umap", cols = "black") & NoLegend()

selected_features <- as.character(sort(as.numeric(tail(rownames(data_S), 400)))[c(1, 101, 201, 301, 400)])
FeaturePlot(data_S, features = selected_features, ncol = 5, pt.size = 1, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()
FeaturePlot(data_S, features = rownames(data_S)[c(101,201,309,601,701)], ncol = 5, pt.size = 1, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()

FeaturePlot(data_S, features = "cell_pseudotime_tree", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()
FeaturePlot(data_S, features = "cell_pseudotime_CC", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()

data_S_CC_tree -> data_S
data_S_CC_tree_v2 -> data_S
cost_matrix_orig <- construct_cost_matrix(data_S, K = 10, dims = 1:10, cell.embedding = "pca")
cost_matrix_orig <- construct_cost_matrix(data_S, K = 10, dims = 1:30, cell.embedding = "dm")
#####Select genes
assay <- "RNA"
DefaultAssay(data_S) <- assay
#data_S <- FindVariableFeatures(data_S)
#all_genes <- data_S@assays[[assay]]@var.features
#require(biclust)
#expr_percent <- apply(biclust::binarize(as.matrix(data_S[[assay]]@data[all_genes, ]), threshold = 0.001), 1, sum)/ncol(data_S)
#genes <- all_genes#[which(expr_percent > 0.01 & expr_percent < 0.5)]
genes <- rownames(data_S)
length(genes)

N <- 1000
cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:10, cell.embedding = "pca")
cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:30, cell.embedding = "dm")

cell_embedding.visualization <- cg_output[["cell_embedding.visualization"]]     
gene_expression  <- cg_output[["gene_expression"]]  
cost_matrix <- cg_output[["cost_matrix"]] 
genes <- cg_output[["genes"]]

################################################
#####Step3: EMD computation in Python
################################################
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_dm30_k10_N", 1000, "/")
if (T){
  dir.create(dir.path, recursive = T)
  write.table(genes, paste0(dir.path, "gene_names.csv"), row.names = F, col.names = F, sep = ",")
  write.table(cost_matrix, paste0(dir.path, "ot_cost.csv"), row.names = F, col.names = F, sep = ",")
  require(Matrix)
  Matrix::writeMM(Matrix(gene_expression, sparse = T), paste0(dir.path, "gene_expression.mtx"))
}



################################################
#####Step5: Visualizing gene trajectory
################################################
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_CC_tree_dm30_k10_N", 1000, "/")
emd_mat <- read_emd_mat(dir.path, file_name = "emd.csv")
str(emd_mat)
#emd_mat <- emd_mat[genes, genes]
plot(1:nrow(emd_mat), sort(emd_mat[1,]))
genes <- rownames(emd_mat)

gene_embedding <- get_gene_embedding(emd_mat, K = 5)$diffu_emb
gene_embedding <- get_gene_embedding(emd_mat[rownames(data_S)[1401:1800], rownames(data_S)[1401:1800]], K = 5)$diffu_emb 

custom.config <- umap.defaults
set.seed(1)
umap.output <- umap(emd_mat, n_neighbors = 200, min_dist = 0.3, input="dist")
plot(umap.output$layout[,1], umap.output$layout[,2])
df_plot <- as.data.frame(cbind(umap.output$layout, as.numeric(rownames(emd_mat))))
colnames(df_plot) <- c("x", "y", "gene_PT")
ggplot(data = df_plot, aes(x=x, y=y, color=gene_PT)) +
  geom_point() & scale_color_viridis(option = "C") &
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_blank(), 
       axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

selected_features <- as.character(sort(as.numeric(tail(rownames(data_S), 400)))[c(1, 101, 201, 301, 400)])
FeaturePlot(data_S, features = selected_features, ncol = 5, pt.size = 1, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()
FeaturePlot(data_S, features = rownames(data_S)[c(101,201,309,601,701)], ncol = 5, pt.size = 1, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()

genes <- c(rownames(data_S)[c(101,201,309,601,701)],selected_features)
genes %in% rownames(emd_mat)
labels <- rep("", nrow(emd_mat))
for (i in 1:10){
  labels[which(rownames(emd_mat) == genes[i])] <- paste0("g",i)
}

df_plot <- as.data.frame(cbind(cbind(umap.output$layout, as.numeric(gsub(".1$", "", rownames(emd_mat)))), labels))
colnames(df_plot) <- c("x", "y", "gene_PT", "label")
ggplot(data = df_plot, aes(x=as.numeric(x), y=as.numeric(y), color=as.numeric(gene_PT), label = label)) +
  geom_point(size = 2.5) + scale_color_viridis(option = "C") + geom_text(size = 10, color = "black") + #ggrepel::geom_text_repel() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggplot(data = df_plot, aes(x=as.numeric(x), y=as.numeric(y), color=as.numeric(gene_PT), label = label)) +
  geom_point(size = 2.5) + scale_color_viridis(option = "C") + #geom_text(size = 10, color = "black") + #ggrepel::geom_text_repel() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())


scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = 0, phi = 90)
scatter3D(gene_embedding[,3],
          gene_embedding[,5],
          gene_embedding[,6], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = 0, phi = 0)

scatter3D(gene_embedding[,10],
          gene_embedding[,11],
          gene_embedding[,9], #alpha = 0.55,
          bty = "b2", colvar = as.numeric(genes),
          main = "Gene Trajectory", pch = 19, cex = 0.5, theta = 0, phi = 0)#,
#col = ramp.col(c("lightgray", hue_pal()(5))))


gene_trajectory <- identify_trajectories(gene_embedding, emd_mat, N = 3, t = 3, cutoff.list = c(8,10.5,12.5))
table(gene_trajectory$selected)

scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 135, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(5))))


