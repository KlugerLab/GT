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


gc_mat_circular <- generate_circular(maxT = 15)

dim(gc_mat_circular)
max(gc_mat_circular)

set.seed(1)
pheatmap::pheatmap(log(gc_mat_circular+1), cluster_cols = F, cluster_rows = F, color = viridis(100), show_rownames = F, show_colnames = F)

#pheatmap::pheatmap(log(gc_mat_circular+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)

data_S <- CreateSeuratObject(gc_mat_circular)
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)#, nfeatures = nrow(data_S))
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 100, verbose = F)
plot(data_S@reductions$pca@stdev)
DimPlot(data_S, reduction = "pca", dims = c(1,2), pt.size = 2)
#data_S$cell_pseudotime_tree <- as.numeric(unlist(strsplit(colnames(data_S), split = "|", fixed = T))[(1:ncol(data_S))*2-1])
#data_S$cell_pseudotime_CC <- as.numeric(unlist(strsplit(colnames(data_S), split = "|", fixed = T))[(1:ncol(data_S))*2])

data_S <- RunDM(data_S, K = 25)
DimPlot(data_S, reduction = "dm", dims = c(1,2))
DimPlot(data_S, reduction = "dm", dims = c(1,2),pt.size = 2.5, cols = "black") & NoLegend() & NoAxes()
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
data_S <- RunTSNE(
  data_S, seed.use=3, dims = 1:30, perplexity = 100, nthreads=12
)
DimPlot(data_S, reduction = "tsne")
DimPlot(data_S, reduction = "tsne", cols = "black", pt.size = 2.5) & NoLegend() & NoAxes()
FeaturePlot(data_S, features = rownames(data_S)[c(1, 101, 201, 301, 401)], ncol = 5, pt.size = 2.5, reduction = "tsne", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()



FeaturePlot(data_S, features = "cell_pseudotime", reduction = "tsne", order = F) & scale_color_viridis(option = "B") & NoAxes()

FeaturePlot(data_S, features = sample(rownames(data_S), size = 5), ncol = 5, pt.size = 1.5, reduction = "tsne", order = T) & scale_color_viridis(option = "D") & NoAxes()

#phate_emb <- phate(data_S@reductions$pca@cell.embeddings)

#DoHeatmap(data_S, features = sample(rownames(data_S), size = 50), size = 0)

data_S <- RunUMAP(
  data_S, dims = 1:30, n.neighbors = 250, min.dist = 0.5#, a = 5, b = 10, repulsion.strength = 5
)
DimPlot(data_S, reduction = "umap", cols = "black") & NoLegend()
DimPlot(data_S, reduction = "umap", cols = "black", pt.size = 2.5) & NoLegend() & NoAxes()
dim(data_S)

FeaturePlot(data_S, features = rownames(data_S)[c(1, 101, 201, 301, 401)], ncol = 5, pt.size = 1.5, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()

FeaturePlot(data_S, features = "cell_pseudotime_tree", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()
FeaturePlot(data_S, features = "cell_pseudotime_CC", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()

data_S_circular_v4 -> data_S
cost_matrix_orig <- construct_cost_matrix(data_S, K = 25, dims = 1:5, cell.embedding = "dm")

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

N <- 500
#cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:10, cell.embedding = "pca")
cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:5, cell.embedding = "dm")

cell_embedding.visualization <- cg_output[["cell_embedding.visualization"]]     
gene_expression  <- cg_output[["gene_expression"]]  
cost_matrix <- cg_output[["cost_matrix"]] 
genes <- cg_output[["genes"]]

################################################
#####Step3: EMD computation in Python
################################################
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm30_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm5_k10_N", 500, "/")
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
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm30_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_circular_dm5_k10_N", 500, "/")
emd_mat <- read_emd_mat(dir.path, file_name = "emd.csv")
str(emd_mat)
#emd_mat <- emd_mat[genes, genes]
plot(1:nrow(emd_mat), sort(emd_mat[1,]))
genes <- rownames(emd_mat)

gene_embedding <- get_gene_embedding(emd_mat, K = 5)$diffu_emb
gene_embedding <- get_gene_embedding(emd_mat[rownames(data_S)[1401:1800], rownames(data_S)[1401:1800]], K = 5)$diffu_emb 


custom.config <- umap.defaults
umap.output <- umap(emd_mat, n_neighbors = 50, min_dist = 0.3, input="dist")
plot(umap.output$layout[,1], umap.output$layout[,2])
df_plot <- as.data.frame(cbind(umap.output$layout, as.numeric(rownames(emd_mat))))

set.seed(1)
tsne_emb <- Rtsne::Rtsne(emd_mat, is_distance = TRUE, perplexity = 100)
str(tsne_emb)

genes <- rownames(data_S)[c(1, 101, 201, 301, 401)]
genes %in% rownames(emd_mat)
labels <- rep("", nrow(emd_mat))
for (i in 1:5){
  labels[which(rownames(emd_mat) == genes[i])] <- paste0("g",i)
}

df_plot <- as.data.frame(cbind(cbind(tsne_emb$Y, as.numeric(rownames(emd_mat))), labels))
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


