generate_cylinder <- function(maxT = 10,
                          N_cells = 5000,
                          N_genes = rep(500,2),
                          meanlog = 0, sdlog = 0.25,
                          scale = 25,
                          seed = 1,
                          sigma_noise = 2.5,
                          sort = TRUE){
  gc_mat_linear <- generate_linear(maxT = maxT,
                                  N_cells = 5000,
                                   N_genes = N_genes[1])
  gc_mat_circular <- generate_circular(maxT = maxT,
                                       N_cells = 5000,
                                   N_genes = N_genes[2])
  set.seed(1)
  gc_mat_cylinder <- rbind(gc_mat_linear, gc_mat_circular[,sample(1:ncol(gc_mat_circular))])
  colnames(gc_mat_cylinder) <- paste0(Tc1, "|", Tc2)
  
  gc_mat_cylinder
}

gc_mat_cylinder <- generate_cylinder(maxT = 15)

dim(gc_mat_cylinder)
max(gc_mat_cylinder)

set.seed(1)
pheatmap::pheatmap(log(gc_mat_cylinder+1), cluster_cols = F, cluster_rows = F, color = viridis(100), show_rownames = F, show_colnames = F)

#pheatmap::pheatmap(log(gc_mat_cylinder+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)

data_S <- CreateSeuratObject(gc_mat_cylinder)
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)#, nfeatures = nrow(data_S))
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 100, verbose = F)
plot(data_S@reductions$pca@stdev)
DimPlot(data_S, reduction = "pca", dims = c(1,2))

data_S <- RunDM(data_S, K = 10)
DimPlot(data_S, reduction = "dm", dims = c(1,2))
par(mar = c(1,1,1,1))
scatter3D(10*data_S@reductions$dm@cell.embeddings[,21],
          10*data_S@reductions$dm@cell.embeddings[,22],
          10*data_S@reductions$dm@cell.embeddings[,23], alpha = 0.55,
          bty = "g",#colvar = as.integer(as.factor(data_S$condition)),# col = col_pal, 
          main = "trajectory", pch = 19, cex = 0.25, theta = 0+90+90, phi = 45)#, col = hue_pal()(length(unique(data_S$condition))))#ramp.col(c("red", "blue", "gray")) ) #colkey = FALSE, 

data_S <- RunTSNE(
  data_S, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/home/jz437/git/FIt-SNE/bin/fast_tsne", ann_not_vpcylinder=FALSE, nthreads=12
)
DimPlot(data_S, reduction = "tsne")


FeaturePlot(data_S, features = "cell_pseudotime", reduction = "tsne", order = F) & scale_color_viridis(option = "B") & NoAxes()

FeaturePlot(data_S, features = sample(rownames(data_S), size = 5), ncol = 5, pt.size = 1.5, reduction = "tsne", order = T) & scale_color_viridis(option = "D") & NoAxes()

#phate_emb <- phate(data_S@reductions$pca@cell.embeddings)

#DoHeatmap(data_S, features = sample(rownames(data_S), size = 50), size = 0)

data_S <- RunUMAP(
  data_S, dims = 1:30, n.neighbors = 500, min.dist = 1#, a = 5, b = 10, repulsion.strength = 5
)
DimPlot(data_S, reduction = "umap", cols = "black") & NoLegend() & NoAxes()
dim(data_S)
FeaturePlot(data_S, features = rownames(data_S)[c(51, 151, 251, 351, 451)], ncol = 5, pt.size = 0.5,  reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()
FeaturePlot(data_S, features = rownames(data_S)[c(551,651, 751, 851, 951)], ncol = 5, pt.size = 0.5, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()

FeaturePlot(data_S, features = "cell_pseudotime_cylinder", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()
FeaturePlot(data_S, features = "cell_pseudotime_CC", reduction = "umap", order = F) & scale_color_viridis(option = "B") & NoAxes()

data_S_cylinder_v4 <- data_S
data_S_cylinder_v5 -> data_S
cost_matrix_orig <- construct_cost_matrix(data_S, K = 10, dims = 1:5, cell.embedding = "dm")
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
#cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:10, cell.embedding = "pca")
cg_output <- coarse_grain(data_S, cost_matrix_orig, genes, N = N, assay = "RNA", dims = 1:5, cell.embedding = "dm")

cell_embedding.visualization <- cg_output[["cell_embedding.visualization"]]     
gene_expression  <- cg_output[["gene_expression"]]  
cost_matrix <- cg_output[["cost_matrix"]] 
genes <- cg_output[["genes"]]

################################################
#####Step3: EMD computation in Python
################################################
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm30_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k10_N", 500, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k10_N", 1000, "_v2/")
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
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k15_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_pca10_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm30_k10_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k10_N", 500, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/simulation/", "new_cylinder_dm5_k10_N", 1000, "_v2/")
emd_mat <- read_emd_mat(dir.path, file_name = "emd.csv")
str(emd_mat)
#emd_mat <- emd_mat[genes, genes]
plot(1:nrow(emd_mat), sort(emd_mat[1,]))
genes <- rownames(emd_mat)

gene_embedding <- get_gene_embedding(emd_mat, K = 5)$diffu_emb
#gene_embedding <- get_gene_embedding(emd_mat[rownames(data_S)[1401:1800], rownames(data_S)[1401:1800]], K = 5)$diffu_emb 

custom.config <- umap.defaults
set.seed(1)
umap.output <- umap(emd_mat, n_neighbors = 200, min_dist = 0.3, input="dist")
plot(umap.output$layout[,1], umap.output$layout[,2])
df_plot <- as.data.frame(cbind(umap.output$layout, as.numeric(gsub(".1$", "", rownames(emd_mat)))))
colnames(df_plot) <- c("x", "y", "gene_PT")
ggplot(data = df_plot, aes(x=x, y=y, color=gene_PT)) +
  geom_point() & scale_color_viridis(option = "C") &
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

FeaturePlot(data_S, features = rownames(data_S)[c(51, 151, 251, 351, 451)], ncol = 5, pt.size = 0.5,  reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()
FeaturePlot(data_S, features = rownames(data_S)[c(551,651, 751, 851, 951)], ncol = 5, pt.size = 0.5, reduction = "umap", order = F) & scale_color_viridis(option = "D") & NoAxes() #& NoLegend()

genes <- rownames(data_S)[c(551,651, 751, 851, 951, 51, 151, 251, 351, 451)]
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





