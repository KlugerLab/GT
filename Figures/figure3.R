data_S_pbmc_10k_MM_reduced -> data_S
if (F){
  #####
  ##### Clustering
  #####
  assay <- "RNA"
  DefaultAssay(data_S) <- assay
  data_S <- FindNeighbors(data_S, dims = 1:30)
  data_S <- FindClusters(data_S, resolution = 0.3)
  data_S$cluster <- data_S$RNA_snn_res.0.3 #integrated
  DimPlot(data_S, reduction = "umap", label = T, label.size = 5, group.by = "cluster")
  
  DefaultAssay(data_S) <- "RNA"
  Idents(data_S) <- "cluster"
  cluster_markers <- FindAllMarkers(
    data_S, only.pos = T, min.diff.pct = 0.1
  )
  cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
  cluster_markers <- cluster_markers[which(cluster_markers$p_val_adj <= 0.05), ]
  top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
  DefaultAssay(data_S) <- "RNA"
  DoHeatmap(data_S, features = top10$gene, slot = "scale.data") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
  
  DotPlot(data_S, features = unique(top10$gene), group.by = "cluster") & RotatedAxis()
  FeaturePlot(data_S, features = unique(top10$gene), ncol = 8)
  cluster_relabel <- c("0" = "CD14+CD16- monocytes",
                       "1" = "CD14+HLA-DR(high) monocytes",
                       "2" = "CD14-CD16+ monocytes",
                       "3" = "CD14- dendritic cells (type-2)")
  data_S$celltype_new <- cluster_relabel[as.character(data_S$cluster)]
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  DimPlot(data_S, group.by = "celltype_new", cols = cbPalette[2:5]) & NoAxes()
  FeaturePlot(data_S, features = c("CD14", "HLA-DMA", "FCGR3A", "CD1C"), ncol = 4) & NoAxes() & scale_color_viridis()
  
  Idents(data_S) <- "celltype_new"
  DoHeatmap(data_S, features = top10$gene, group.colors = cbPalette[2:5][c(3,4,2,1)], slot = "scale.data", size = 2.5) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
  DotPlot(data_S, features = unique(top10$gene), group.by = "celltype_new") & RotatedAxis()
  
}
DimPlot(data_S)
#####Select genes
assay <- "RNA"
DefaultAssay(data_S) <- assay
data_S <- FindVariableFeatures(data_S, nfeatures = 2000)
all_genes <- data_S@assays[[assay]]@var.features
require(biclust)
expr_percent <- apply(biclust::binarize(as.matrix(data_S[[assay]]@data[all_genes, ]), threshold = 0), 1, sum)/ncol(data_S)

dir.path <- paste0("/data/rihao/project_with_xiuyuan/pbmc_10k_v3/result/data/", "MM_dm5_graph_N", 1000, "_v2/")
emd_mat <- read_emd_mat(dir.path, file_name = "emd.csv")
str(emd_mat)
plot(1:nrow(emd_mat), sort(emd_mat[1,]))
gene_embedding <- get_gene_embedding(emd_mat, K = 5)$diffu_emb
plot(gene_embedding[,1],
     gene_embedding[,2])
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = 180, phi = 30)
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = -90, phi = 30)

gene_labels <- paste("-----", rownames(gene_embedding))
names(gene_labels) <- rownames(gene_embedding)
genes_selected <- c("CCR2", "ICAM2", "FCGR3A",  "SELL", "C1QA", "C1QB",
                    "CD1C", "CLEC10A", "CD2", "CD72", "CCR5", 
                    "PKIB", "RETN", "CLEC5A", "CSF1R")

DefaultAssay(data_S) <- "alra"
FeaturePlot(data_S, features = genes_selected, ncol = 5) & NoAxes() #& scale_color_viridis()

DefaultAssay(data_S) <- "RNA"
genes_selected <- c("CD2")
genes_selected %in% rownames(gene_embedding)


gene_labels[names(gene_labels) %notin% genes_selected] <- ""

gene_trajectory <- identify_trajectories(gene_embedding, emd_mat, N = 3, t = 5, cutoff.list = c(8.5,10.5,11))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(10,25,30), quantile = 0.1)#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(5,14,6), quantile = 0.02)#,7.75))
table(gene_trajectory$selected)
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(11,21,8), quantile = 0.02, K=5)#,7.75))
table(gene_trajectory$selected)
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(11,21,8), quantile = 0.05, K=5)#,7.75))
table(gene_trajectory$selected)

gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(12,22,10), quantile = 0.025, K = 5)#,7.75))
table(gene_trajectory$selected)


scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3)))) #&   #"lightgray", 
text3D(gene_embedding[,1],
       gene_embedding[,2],
       gene_embedding[,3],  labels = gene_labels,
       add = T, colkey = FALSE, cex = 0.5)
scatter3D(gene_embedding[,2],
          gene_embedding[,3],
          gene_embedding[,4], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1.5, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3)))) #&   #"lightgray", 
scatter3D(gene_embedding[,2],
          gene_embedding[,3],
          gene_embedding[,4], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = .75, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3)))) #&   #"lightgray", 

DefaultAssay(data_S) <- "alra"
FeaturePlot(data_S, features = c("SLC4A3", "FIBCD1", "SERPINB10"), ncol = 3) & NoAxes() #& scale_color_viridis()
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 180, phi = 45,
          col = ramp.col(c("lightgray",hue_pal()(3)))) #&   # 
sort(apply(gene_embedding[which(gene_trajectory$selected == "Trajectory-3"),]^2, 2, sum))
par(mar = c(1,1,1,1)+2.25)
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder1,
          main = "", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))
scatter3D(gene_embedding[,2],
          gene_embedding[,6],
          gene_embedding[,3], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder2,
          main = "", pch = 19, cex = 1, theta = 0, phi = 90,
          col = ramp.col(viridis(8)))

genes_selected <- c("CCR2", "ICAM2", "FCGR3A",  "SELL", "C1QA", "C1QB",
                    "CD1C", "CLEC10A", "CD2", "CD72", "CCR5", 
                    "PKIB", "RETN", "CLEC5A", "CSF1R")

gene_list <- list()
for (i in 1:3){
  message(paste0("Trajectory-", i))
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  message(paste(rev(genes), collapse = ", "))
  gene_list[[i]] <- genes
}
for (i in 1:3){
  print(genes_selected[which(genes_selected %in% gene_list[[i]])])
}



N_traj <- 3
N_bin <- 5
data_S <- add_gene_set_score(data_S, gene_trajectory, N_bin = N_bin, trajectories = 1:N_traj, binarize = T, assay = "alra", echo = T, reverse = c(F, F, T))
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:N_bin), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()

genes <- all_genes[which(expr_percent > 0.005 & expr_percent < 0.25)]
gene_trajectory_sub <- gene_trajectory[genes, ]
data_S <- add_gene_set_score(data_S, gene_trajectory_sub, N_bin = N_bin, trajectories = 1:N_traj, binarize = T, assay = "alra", echo = T, reverse = c(F, F, T))
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",3,"_genes", 1:N_bin), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()

