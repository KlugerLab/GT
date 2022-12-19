data_S_WLS_combined_E14.5 -> data_S
DimPlot(data_S, reduction = "umap", group.by = "orig.ident", shuffle = T, cols = c("darkorange2", "chartreuse3")) & NoAxes()
DimPlot(data_S, group.by = "cell_type", cols = c( "darkred", "darkblue", "bisque4")) & NoAxes()

DefaultAssay(data_S) <- "alra"
FeaturePlot(data_S, features = c("Dkk2", "Dkk1", "Sox2", "Lef1", "Ptch1", "Ccnb1"), order = T, ncol = 3) & NoAxes()
DimPlot(data_S)
DefaultAssay(data_S) <- "RNA"
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
  data_S$cell_type <- "UD"
  data_S$cell_type[which(data_S$cluster %in% c(0,4))] <- "LD"
  data_S$cell_type[which(data_S$cluster %in% c(7))] <- "DC"
  DimPlot(data_S, reduction = "umap", label = T, label.size = 5, group.by = "cell_type")
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  DimPlot(data_S, group.by = "cell_type", cols = cbPalette[2:4], split.by = "orig.ident") & NoAxes()
  FeaturePlot(data_S, features = c("CD14", "HLA-DMA", "FCGR3A", "CD1C"), ncol = 4) & NoAxes() & scale_color_viridis()
  DimPlot(data_S, group.by = "orig.ident") & NoAxes()
  
}

dir.path <- paste0("/data/rihao/project_with_xiuyuan/dermal_WLS/", "ALL_E14_dm10_graph_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/dermal_WLS/", "ALL_E14_dm5_k10_graph_N", 1000, "/")
dir.path <- paste0("/data/rihao/project_with_xiuyuan/dermal_WLS/", "ALL_E14_dm10_k10_graph_N", 1000, "/")
emd_mat <- read_emd_mat(dir.path, file_name = "emd.csv")
str(emd_mat)
plot(1:nrow(emd_mat), sort(emd_mat[1,]))








#gene_embedding <- get_gene_embedding(emd_mat, K = 25)$diffu_emb
gene_embedding <- get_gene_embedding(emd_mat, K = 5)$diffu_emb #E14.5 combined

scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0)
scatter3D(gene_embedding[,1],
          gene_embedding[,3],
          gene_embedding[,4], alpha = 0.55,
          bty = "b2", 
          main = "trajectory", pch = 19, cex = 0.5, theta = 180, phi = 90)

#E14 combined new
gene_trajectory <- identify_trajectories(gene_embedding, emd_mat, N = 3, t = 3, cutoff.list = c(7.5,12,7.5))#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(3,16,2), quantile = 0.05)#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(2,12,2), quantile = 0.02)#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(8,16,2), quantile = 0.02, K=5)#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(11,22,3), quantile = 0.05, K=5)#,7.75))
gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(16,22,3), quantile = 0.2, K=5)#,7.75))

gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(9,16,4), quantile = 0.02, K=5)#,7.75))

#gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(19,12,6), quantile = 0.02, K=5)#,7.75))
#gene_trajectory <- identify_trajectories_v2(gene_embedding, emd_mat, N = 3, t = c(1,12,2), quantile = 0.01)#,7.75))
table(gene_trajectory$selected)

plot(gene_embedding[,1],
     gene_embedding[,2])

zoomin_view_cells <- which(gene_embedding[,1] > -0.01)
scatter3D(gene_embedding[,2],
          gene_embedding[,3],
          gene_embedding[,5], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
scatter3D(gene_embedding[zoomin_view_cells,1],
          gene_embedding[zoomin_view_cells,2],
          gene_embedding[zoomin_view_cells,5], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory[zoomin_view_cells,]$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 0, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))

sort(apply(gene_embedding[which(gene_trajectory$selected == "Trajectory-3"),]^2, 2, sum))
sort(apply(gene_embedding[which(gene_trajectory$selected != "Other"),]^2, 2, sum))
scatter3D(gene_embedding[,6],
          gene_embedding[,1],
          gene_embedding[,12], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
scatter3D(gene_embedding[,6],
          gene_embedding[,1],
          gene_embedding[,12], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder3,
          main = "", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
scatter3D(gene_embedding[,2],
          gene_embedding[,3],
          gene_embedding[,5], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))


scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
scatter3D(gene_embedding[,2],
          gene_embedding[,1],
          gene_embedding[,5], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
scatter3D(gene_embedding[,5],
          gene_embedding[,3],
          gene_embedding[,6], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(5))))


scatter3D(gene_embedding[,5],
          gene_embedding[,3],
          gene_embedding[,6], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder3,
          main = "", pch = 19, cex = 0.5, theta = 90, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))

gene_labels <- paste("-----", rownames(gene_embedding))
names(gene_labels) <- rownames(gene_embedding)
genes_selected <- c("Sox2", "Top2a", "Ccnb1", "Rrm2" ,    "Rad51"   , "Rad51ap1" ,"Brip1"   , "E2f8" ,    "Blm" ,     "Rrm1",
                    gene_list[[3]][which(toupper(gene_list[[3]]) %in% g2m.genes)]
)
gene_labels[names(gene_labels) %notin% genes_selected] <- ""

scatter3D(gene_embedding[,5],
          gene_embedding[,3],
          gene_embedding[,6], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 0, phi = 0,
          col = ramp.col(c(hue_pal()(3)))) #&   #"lightgray", 
text3D(gene_embedding[,5],
       gene_embedding[,3],
       gene_embedding[,6],  labels = gene_labels,
       add = T, colkey = FALSE, cex = 0.5)

scatter3D(gene_embedding[,6],
          gene_embedding[,1],
          gene_embedding[,12], #alpha = 0.55,
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 0.5, theta = 45, phi = 0,
          col = ramp.col(c("lightgray", hue_pal()(3))))
text3D(gene_embedding[,6],
       gene_embedding[,1],
       gene_embedding[,12],  labels = gene_labels,
       add = T, colkey = FALSE, cex = 0.5)

N_traj <- 3
gene_list <- list()
for (i in 1:N_traj){
  message(paste0("Trajectory-", i))
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  gene_list[[paste0("Trajectory-", i)]] <- genes
  message(paste(genes, collapse = ", "))
}
for (i in 1:N_traj) print("Sox2" %in% gene_list[[i]])
for (i in 1:N_traj) print("Top2a" %in% gene_list[[i]])

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes
g2m.genes
gene_list[[3]][which(toupper(gene_list[[3]]) %in% s.genes)]
#i <- 3
#data_S_WLS_combined_E14.5$T3_proportion <- apply(as.matrix(data_S_WLS_combined_E14.5[["RNA"]]@data[gene_list[[i]],])>0, 2, sum)/length(gene_list[[i]])
#FeaturePlot(data_S_WLS_combined_E14.5, features = c("T1_proportion", "T2_proportion", "T3_proportion"))


N_bin <- 7
data_S_WLS_combined_E14.5_v3 -> data_S
data_S <- add_gene_set_score(data_S, gene_trajectory, N_bin = N_bin, trajectories = 1:N_traj, binarize = T, assay = "alra", echo = T)


FeaturePlot(data_S, pt.size = 0.01, features = paste0("Trajectory", 3,"_genes", 1:N_bin), ncol = N_bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoAxes() & theme(plot.title = element_text(size = 10)) & NoLegend() 

FeaturePlot(subset(data_S, orig.ident != "E14WLS"), pt.size = 0.1, features = paste0("Trajectory",3,"_genes", 1:N_bin), ncol = N_bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoAxes() & theme(plot.title = element_text(size = 10)) & NoLegend() 









data_S_WLS_combined_E14.5 <- AddModuleScore(data_S_WLS_combined_E14.5, features = gene_list)
FeaturePlot(data_S_WLS_combined_E14.5, features = c("Cluster1", "Cluster2", "Cluster3"))

data_S <- data_S_WLS_combined_E14.5
VlnPlot(data_S, features = c("T1_proportion", "T2_proportion", "T3_proportion"), pt.size = 0.0, group.by = "orig.ident")
VlnPlot(data_S, features = c("Cluster1", "Cluster2", "Cluster3"), group.by = "orig.ident")

#data_S <- subset(data_S_WLS_combined_E14.5, T1_proportion > 0.1)
data_S <- subset(data_S_WLS_combined_E14.5, Cluster3 > 0)
RidgePlot(data_S, features = "Cluster3", group.by = "orig.ident")
data_S <- data_S[gene_list[[3]],]

DefaultAssay(data_S) <- "RNA"
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S, nfeatures = 2000) #blood cells
#data_S <- FindVariableFeatures(data_S, nfeatures = 5000) #neurons
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 30, verbose = F)
plot(data_S@reductions$pca@stdev)
data_S <- RunUMAP( data_S, dims = 1:30, min.dist = 1)#, min.dist = 1#, spread = 2#, n.neighbors = 30)
DimPlot(data_S, reduction = "umap", group.by = "orig.ident", shuffle = T) #& NoAxes()

#data_S$PC3 <- data_S@reductions$pca@cell.embeddings[,3]
#RidgePlot(data_S, features = "PC3", group.by = "orig.ident")



