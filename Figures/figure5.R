data_S_WLS_combined_E14.5 -> data_S

N_bin <- 7
data_S <- add_gene_set_score(data_S, gene_trajectory, N_bin = N_bin, trajectories = 1:N_traj, binarize = T, assay = "alra", echo = T)


FeaturePlot(data_S, pt.size = 0.01, features = paste0("Trajectory",3,"_genes", 1:N_bin), ncol = N_bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoAxes() & theme(plot.title = element_text(size = 10)) & NoLegend() 

FeaturePlot(subset(data_S, orig.ident == "E14WLS"), pt.size = 0.1, features = paste0("Trajectory",2,"_genes", 1:N_bin), ncol = N_bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoAxes() & theme(plot.title = element_text(size = 10)) & NoLegend() 

RidgePlot(data_S, features = paste0("Trajectory",1,"_genes", 1:N_bin), group.by = "orig.ident")
#pheatmap::pheatmap(cor(t(as.matrix(data_S[["RNA"]]@data[gene_list[[2]],]))), show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F)



data_S <- subset(data_S_WLS_combined_E14.5, orig.ident == "E14WLS" & cell_type %in% c("UD", "DC"))
dim(data_S)
summary_vec <- apply(as.matrix(data_S@assays$RNA@data[gene_list[[1]],])>0, 1, sum) > 0.01*ncol(data_S)
names(summary_vec[which(!summary_vec )])




dists <- as.matrix(dist(data_S@reductions$pca@cell.embeddings[,1:30]))
rownames(dists) <- colnames(data_S)
colnames(dists) <- colnames(data_S)

smooth_KNN <- function(orig.vec, dists, K = 5){
  orig.cells <- rownames(dists)[which(orig.vec == 1)]
  neighbors <- c()
  for (i in orig.cells){
    neighbors <- c(neighbors, names(sort(dists[i,])[1:(K+1)]))
  }
  neighbors <- unique(neighbors)
  new.vec <- orig.vec
  names(new.vec) <- rownames(dists)
  new.vec[neighbors] <- 1
  new.vec
}

#data_S <- RunALRA(data_S)
cutoff <- 0.5
data_S$tmp1 <- binarize(data_S$Trajectory1_genes1, cutoff)
data_S$tmp2 <- binarize(data_S$Trajectory1_genes2, cutoff)
data_S$tmp3 <- binarize(data_S$Trajectory1_genes3, cutoff)
data_S$tmp4 <- binarize(data_S$Trajectory1_genes4, cutoff)
data_S$tmp5 <- binarize(data_S$Trajectory1_genes5, cutoff)
data_S$tmp6 <- binarize(data_S$Trajectory1_genes6, cutoff)
data_S$tmp7 <- binarize(data_S$Trajectory1_genes7, cutoff)

if (F){
  data_S$tmp1 <- smooth_KNN(data_S$tmp1, dists, K = 5)
  data_S$tmp2 <- smooth_KNN(data_S$tmp2, dists, K = 5)
  data_S$tmp3 <- smooth_KNN(data_S$tmp3, dists, K = 5)
  data_S$tmp4 <- smooth_KNN(data_S$tmp4, dists, K = 5)
  data_S$tmp5 <- smooth_KNN(data_S$tmp5, dists, K = 5)
  data_S$tmp6 <- smooth_KNN(data_S$tmp6, dists, K = 5)
  data_S$tmp7 <- smooth_KNN(data_S$tmp7, dists, K = 5)
}

#data_S$tmp6 <- binarize(data_S$Trajectory1_genes6)
data_S$all <- data_S$tmp7
data_S$tmp6 <- pmax(data_S$tmp6 - data_S$all, 0)
data_S$all <- pmax(data_S$tmp6, data_S$all)
data_S$tmp5 <- pmax(data_S$tmp5 - data_S$all, 0)
data_S$all <- pmax(data_S$tmp5, data_S$all)
data_S$tmp4 <- pmax(data_S$tmp4 - data_S$all, 0)
data_S$all <- pmax(data_S$tmp4, data_S$all)
data_S$tmp3 <- pmax(data_S$tmp3 - data_S$all, 0)
data_S$all <- pmax(data_S$tmp3, data_S$all)
data_S$tmp2 <- pmax(data_S$tmp2 - data_S$all, 0)
data_S$all <- pmax(data_S$tmp2, data_S$all)
data_S$tmp1 <- pmax(data_S$tmp1 - data_S$all, 0)
FeaturePlot(subset(data_S, orig.ident != "E14WLS"), pt.size = 0.1, features = paste0("tmp",1:N_bin), ncol = N_bin, order = T, cols = c("lightgray", "darkred"))  & NoLegend()& NoAxes() & theme(plot.title = element_text(size = 10))
DimPlot(subset(data_S, orig.ident != "E13WLS"), pt.size = 0.1, group.by = paste0("tmp",1:N_bin), ncol = N_bin, order = T, cols = c("lightgray", "darkred"))  & NoLegend()& NoAxes() & theme(plot.title = element_text(size = 10))
data_S$stage <- "Other"
data_S$stage[which(data_S$tmp1 == 1)] <- "Stage1"
data_S$stage[which(data_S$tmp2 == 1)] <- "Stage2"
data_S$stage[which(data_S$tmp3 == 1)] <- "Stage3"
data_S$stage[which(data_S$tmp4 == 1)] <- "Stage4"
data_S$stage[which(data_S$tmp5 == 1)] <- "Stage5"
data_S$stage[which(data_S$tmp6 == 1)] <- "Stage6"
data_S$stage[which(data_S$tmp7 == 1)] <- "Stage7"
DimPlot(data_S, group.by = "stage", split.by = "orig.ident", cols = c("lightgray",hue_pal()(7))) & NoAxes()
#DimPlot(data_S, group.by = c("tmp1", "stage"), cols = c("lightgray",hue_pal()(7)))

table(data_S$orig.ident, data_S$stage)
data_S_sub <- subset(data_S, orig.ident == "E14WLS")
table(data_S_sub$Phase, data_S_sub$stage)
summary_table <- t(table(data_S_sub$Phase, data_S_sub$stage))/apply(table(data_S_sub$Phase, data_S_sub$stage), 2, sum)
plot(1:7, summary_table[-1,1])

df_plot <- data.frame(proportion = summary_table[-1,1], 
                      cell_number = as.numeric(table(data_S_sub$stage)[-1])
)
df_plot$sd <- sqrt(df_plot$proportion*(1-df_plot$proportion)/df_plot$cell_number)

df_plot_CTL <- df_plot
df_plot_MUT <- df_plot  

df_plot <- rbind(df_plot_CTL, df_plot_MUT)
df_plot$sample <- c(rep("CTL", 7), rep("MUT", 7))
df_plot$linetype <- c(rep("orig", 11), rep("dashed", 3))
df_plot$sample_new <- paste0(df_plot$sample, "-", df_plot$linetype)
df_plot$stage <- c(1:7, 1:7)
df_plot <- rbind(df_plot, df_plot[11,])
df_plot[15, "linetype"] <- "dashed"
df_plot[15, "sample_new"] <- "MUT-dashed"
df_plot$cell_number_tmp <- c(df_plot$cell_number[-1],0)

df_plot[12:15, "cell_number_tmp"] <- 0

p<-ggplot(df_plot[1:11,], aes(x=stage, y=proportion, group=sample_new)) +
  geom_line(aes(color=sample, size = cell_number_tmp, linetype=linetype))+
  geom_point(aes(color=sample, size = cell_number)) + 
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.0))
p +
  geom_line(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample), linetype = "dashed") +
  geom_point(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample, size = cell_number)) + 
  geom_errorbar(data = df_plot[12:15,], aes(x=stage, y=proportion, ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.0)) +
  theme_classic()


p<-ggplot(df_plot[1:11,], aes(x=stage, y=proportion, group=sample_new)) +
  geom_line(aes(color=sample, linetype=linetype))+
  geom_point(aes(color=sample, size = cell_number)) + 
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.0))
p +
  theme_classic()
p +
  geom_line(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample), linetype = "dashed") +
  geom_point(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample, size = cell_number)) + 
  geom_errorbar(data = df_plot[12:15,], aes(x=stage, y=proportion, ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.0)) +
  theme_classic()

if (F) {
  p<-ggplot(df_plot[1:11,], aes(x=stage, y=proportion, group=sample_new)) +
    geom_line(aes(color=sample, size = cell_number_tmp, linetype=linetype))+
    geom_point(aes(color=sample, size = cell_number)) + 
    geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                  position=position_dodge(0.0))
  p +
    geom_line(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample), linetype = "dashed") +
    geom_point(data = df_plot[12:15,], aes(x=stage, y=proportion, color=sample, size = cell_number)) + 
    geom_errorbar(data = df_plot[12:15,], aes(x=stage, y=proportion, ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                  position=position_dodge(0.0)) +
    theme_classic()
}


df_plot <- data.frame(proportion = summary_table[-1,1], 
                      cell_number = as.numeric(table(data_S_sub$stage)[-1])
)
df_plot$sd <- sqrt(df_plot$proportion*(1-df_plot$proportion)/df_plot$cell_number)
df_plot$sample <- rep("CTL", 7)
df_plot$stage <- 1:7
df_plot$cell_number_tmp <- c(df_plot$cell_number[-1],0)
p<-ggplot(df_plot, aes(x=stage, y=proportion, group=sample)) +
  geom_line(aes(color=sample))+ #, size = cell_number_tmp
  geom_point(aes(color=sample, size = cell_number)) + 
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.1, 
                position=position_dodge(0.0))+
  theme_classic()
p




data_S_WLS_combined_E14.5_v2 <- data_S
data_S_WLS_combined_E14.5_v3 <- data_S
#########Lef1 level change across stages
data_S <- subset(data_S_WLS_combined_E14.5_v3, orig.ident == "E14WLS")
data_S$lef1 <- data_S@assays$alra@data["Lef1",]
lef1_df <- data_S@meta.data[, c("lef1", "stage", "orig.ident")]
lef1_df <- lef1_df %>% filter(stage != "Other")
#lef1_df <- summarize(lef1_df, avarage = mean(lef1))
lef1_df_mean <- lef1_df %>% group_by(stage) %>% dplyr::summarise(mean = mean(lef1),
                                                          sd = sd(lef1))
lef1_df <- merge(lef1_df, lef1_df_mean)

lef1_df 
lef1_df %>% 
  ggplot(mapping = aes(x = stage, y = lef1)) + 
  geom_boxplot(aes(fill = mean), width = 0.5) +
  scale_fill_gradient2(low = "white", high = "forestgreen", limits=c(0,1.8)) +
  geom_point(data = lef1_df_mean, 
             mapping = aes(x = stage, y = mean),
             color="black") +
  geom_line(data = lef1_df_mean, 
            mapping = aes(x = stage, y = mean, group=1))+ ylim(c(0,3))+
  theme_classic()




