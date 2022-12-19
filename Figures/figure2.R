set.seed(1)
Tc2 <- sort(runif(n = 50, min = 0, max = 5.5))
set.seed(2)
Tg2 <- sort(runif(n = 15, min = 0, max = 5.5))

gene_cell_mat <- generate_linear(Tc2, Tg2, seed = 1000)
pheatmap::pheatmap(log(gene_cell_mat+1), cluster_cols = F, cluster_rows = F, color = viridis(100))

