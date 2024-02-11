########################## TRAJ ANALYSIS MARCH 30, 2023 - Apr 10
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/")
combined.all <- readRDS("./combinedALL.rds")

## DE TRAJECTORY 
library(slingshot); library(SingleCellExperiment)
library(RColorBrewer); library(scales)
library(viridis); library(UpSetR)
library(pheatmap); library(msigdbr)
library(fgsea); library(knitr)
library(ggplot2); library(gridExtra)
library(tradeSeq); library(condiments)

# combined.all <- SetIdent(combined.all, value = "sample")
# infected <- WhichCells(combined.all, idents = "Infected")
# infected <- subset(combined.all, cells = infected)
# mock <- WhichCells(combined.all, idents = "Mock")
# mock <- subset(combined.all, cells = mock)
# combined.all <- SetIdent(combined.all, value = "main")
# infected <- SetIdent(infected, value = "main")
# mock <- SetIdent(mock, value = "main")
# 
# counts <- GetAssayData(combined.all, slot="counts", assay="RNA") 
# genes.percent.expressed <- rowMeans(counts>0 )*100  
# genes.filter <- names(genes.percent.expressed[genes.percent.expressed>5])  #select genes expressed in at least 1% of cells
# counts.sub <- counts[genes.filter,]
DefaultAssay(combined.all) <- "RNA"
combined.all <- FindVariableFeatures(combined.all, nfeatures = 500)
var.features <- combined.all@assays$RNA@var.features
# combined.all <- subset(combined.all, features = combined.all@assays$RNA@var.features)

sce <- as.SingleCellExperiment(combined.all, assay = "RNA")
sce <- sce[var.features,]
# sce_inf <- as.SingleCellExperiment(infected, assay = "RNA")
# sce_mock <- as.SingleCellExperiment(mock, assay = "RNA")

shuffle <- sample(ncol(sce))
layout(matrix(1, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce)$sample)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n", 
       legend = levels(factor(colData(sce)$sample)))

# Regions with a high score indicate that the local cell distribution according to treatment label 
# is unbalanced compared the overall distribution. Here, we see that, while there are some small regions of imbalance, 
# the global path along the development axis is well-balanced. 
# This means that we can fit a global trajectory to the full dataset. 

scores <- condiments::imbalance_score(
  reducedDims(sce)$UMAP, 
  condition = colData(sce)$sample,
  k = 20, smooth = 40)

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topright", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)

sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$integrated_snn_res.2, 
                 start.clus = '12', approx_points = 300)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$traj.clusters, 
                 start.clus = '15', approx_points = 300)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$main, dist.method = "mnn",
                 start.clus = 'HSC', approx_points = 300, extend = 'n')
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$trajectory_clusters,
                 start.clus = '3', approx_points = 300,  dist.method = "mnn")

sce_inf <- slingshot(sce_inf, reducedDim = 'UMAP', clusterLabels = colData(sce_inf)$integrated_snn_res.2,
                     start.clus = 'HSC', approx_points = 150)
sce_mock <- slingshot(sce_mock, reducedDim = 'UMAP', clusterLabels = colData(sce_mock)$integrated_snn_res.2,
                      start.clus = 'HSC', approx_points = 300)



# sce <- SlingshotDataSet(sce)
# sce_inf <- SlingshotDataSet(sce_inf)
# sce_mock <- SlingshotDataSet(sce_mock)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(combined.all$main, hue_pal())
# inf_cell_colors_clust <- cell_pal(infected$main, hue_pal())
# mock_cell_colors_clust <- cell_pal(mock$main, hue_pal())
cols <- unique(cell_colors_clust)
cols <- c("#F8766D", "#EC8239", "#DB8E00", "#C79800", "#AEA200", "#8FAA00", "#64B200", "#00B81B",
          "#00BD5C", "#00C085", "#00C1A7", "#00BFC4", "#00BADE", "#00B2F3", "#00A6FF", "#7C96FF", 
          "#B385FF", "#D874FD", "#EF67EB", "#FD61D3", "#FF63B6", "#FF6B94")
vars <- as.data.frame(table(combined.all$main))
levels(vars$Var1)

# integrated
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 16, bty = "n", cex = 2 / 3)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = T)

layout(matrix(1:2, nrow = 1))
par(mar = c(3,4,3,6))
windows(width = 5, height = 6)
# infected
plot(reducedDim(SlingshotDataSet(sce_inf)), col = inf_cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 16, bty = "n", cex = 2 / 3, inset = c(-0.5, 0), xpd = T)
lines(SlingshotDataSet(sce_inf), lwd = 2, type = 'lineages', col = 'black', show.constraints = T )
#mock
plot(reducedDim(SlingshotDataSet(sce_mock)), col = mock_cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 16, bty = "n", cex = 2 / 3, inset = c(-0.5, 0), xpd = T)
lines(SlingshotDataSet(sce_mock), lwd = 2, type = 'lineages', col = 'black', show.constraints = T)


crv1 <- getCurves(sce)
crv1
crv <- SlingshotDataSet(crv1)
crv@lineages <- crv@lineages[c(1,2,3,4,5,6,7,9,10,11,12,14)]
crv@curves <- crv@curves[c(1,2,3,4,5,6,7,9,10,11,12,14)]
par(mar = c(3,4,3,10))
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 20, bty = "n", cex = 1,  inset = c(-0.23, 0), xpd = T)
lines(crv, lwd = 2,  col = 'black')

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 1:5],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 1:5])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 6:7],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 6:7])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", c(12,14)],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", c(12,14)])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 10:11],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 10:11])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 6],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 6])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 7],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 7])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 8],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 8])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 10],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 10])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 11],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 11])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 12],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 12])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 13],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 13])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 14],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 14])

combined.all$pt1 <- sce$slingPseudotime_1
combined.all$pt2 <- sce$slingPseudotime_2
combined.all$pt3 <- sce$slingPseudotime_3
combined.all$pt4 <- sce$slingPseudotime_4
combined.all$pt5 <- sce$slingPseudotime_5
combined.all$pt6 <- sce$slingPseudotime_6
combined.all$pt7 <- sce$slingPseudotime_7
combined.all$pt8 <- sce$slingPseudotime_8
combined.all$pt9 <- sce$slingPseudotime_9
combined.all$pt10 <- sce$slingPseudotime_10
combined.all$pt11 <- sce$slingPseudotime_11
combined.all$pt12 <- sce$slingPseudotime_12
combined.all$pt13 <- sce$slingPseudotime_13
combined.all$pt14 <- sce$slingPseudotime_14

combined.all <- SetIdent(combined.all, value = "main")
FeaturePlot(combined.all, c("pt1", "pt2", "pt3", "pt4", "pt5", "pt6"), label = T, label.size = 2.5)
FeaturePlot(combined.all, c("pt7", "pt8", "pt9", "pt10", "pt11", "pt12", "pt13", "pt14"), label = T, label.size = 2.5)

###########3 Plot Density
layout(matrix(c(1, 1, 2, 3), 2))
# integrated
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 16, bty = "n", cex = 2.7 / 3)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = T)

par(mar = c(4.5, 4, 1, 1))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_6, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_7, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_12, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_14, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_9, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))

ds <- list(Mock = density(na.omit(slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])),
           HIV_Infected = density(na.omit(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9])))
xlim <- range(c(ds$Mock$x, ds$HIV_Infected$x))
ylim <- range(c(ds$Mock$y, ds$HIV_Infected$y))
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")
polygon(c(min(ds$Mock$x),ds$Mock$x,max(ds$Mock$x)),
        c(0,ds$Mock$y,0), col = rgb(0,0,0,.5))
polygon(c(min(ds$HIV_Infected$x),ds$Mock$x,max(ds$Mock$x)),
        c(0,ds$HIV_Infected$y,0), col = alpha(brewer.pal(4,'Set1')[1], alpha = .5))
legend("topright", legend = c("Mock", "Infected"), 
       fill = alpha(c(1, brewer.pal(3, "Set1")[1]), alpha = .5), bty = "n")

par(mar = c(5, 4, 4, 2) + .1)

# plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
#      col = alpha(c(3, 4)[factor(colData(sce)$sample)][shuffle], alpha = .5))
# lines(SlingshotDataSet(sce), type = 'lineages', show.constraints = TRUE)
# legend("topright", pch = 16, col = 3:4, bty = "n", legend = levels(factor(colData(sce)$sample)))
# layout(1)

### Plot all Curves
sce1 <- SlingshotDataSet(sce)
nc <- 3
pt <- slingPseudotime(sce1)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sce1), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sce1, lwd = 2, col = 'black', type = 'lineages') 
}

df <- as.data.frame(sce@assays@data$counts)
write.csv(df, "/storage1/fs1/liang.shan/Active/scRNA/expr_matrix.csv")
combined.all@meta.data$pseutotime1 <- sce@colData$slingPseudotime_1
combined.all@meta.data$pseutotime2 <- sce@colData$slingPseudotime_2
combined.all@meta.data$pseutotime3 <- sce@colData$slingPseudotime_3
combined.all@meta.data$pseutotime4 <- sce@colData$slingPseudotime_4
combined.all@meta.data$pseutotime5 <- sce@colData$slingPseudotime_5
combined.all@meta.data$pseutotimec6 <- sce@colData$slingPseudotime_c6
combined.all@meta.data$pseutotime7 <- sce@colData$slingPseudotime_7
combined.all@meta.data$pseutotime8 <- sce@colData$slingPseudotime_8
combined.all@meta.data$pseutotime9 <- sce@colData$slingPseudotime_9
combined.all@meta.data$pseutotime10 <- sce@colData$slingPseudotime_10
combined.all@meta.data$pseutotime11 <- sce@colData$slingPseudotime_11
combined.all@meta.data$pseutotime12 <- sce@colData$slingPseudotime_12
combined.all@meta.data$pseutotime13 <- sce@colData$slingPseudotime_13
combined.all@meta.data$pseutotime14 <- sce@colData$slingPseudotime_14
md <- as.data.frame(combined.all@meta.data)
write.csv(md, "/storage1/fs1/liang.shan/Active/scRNA/metaData.csv")
genes <- as.data.frame(row.names(sce))
write.csv(genes, "/storage1/fs1/liang.shan/Active/scRNA/gene_names.csv")
cell_groups <- t(read.csv("/storage1/fs1/liang.shan/Active/scRNA/cell_groupings.csv"))
cell_groups <- cell_groups[-c(1),]
write.table(cell_groups, "/storage1/fs1/liang.shan/Active/scRNA/transpose_cell_groups.csv", sep = ",", col.names = F)
pseudotimes <- t(read.csv("/storage1/fs1/liang.shan/Active/scRNA/pseudotimes.csv"))
pseudotimes <- pseudotimes[-c(1),]
pseudotimes[is.na(pseudotimes)] <- 0
write.table(pseudotimes, "/storage1/fs1/liang.shan/Active/scRNA/transpose_pseudotimes.csv", sep = ",", col.names = F)
write.table(combined.all@reductions$umap@cell.embeddings, "/storage1/fs1/liang.shan/Active/scRNA/cell_embeddings.csv", sep = ",", col.names = F)
# # Find number of knots
# icMat <- evaluateK(counts = sce@assays@data$counts,
#                    pseudotime = sce$slingshot@assays@data$pseudotime,
#                    cellWeights = sce$slingshot@assays@data$weights,
#                    conditions = factor(colData(sce)$sample),
#                    nGenes = 200,
#                    k = 5:13)

# fit the NB-GAMs (Negative Binomial General Additive Models) using ? knots, 
# based on the pseudotime and cell-level weights estimated by Slingshot. 
# We use the conditions argument to fit separate smoothers for each condition.
set.seed(12)
sce <- fitGAM(sce, conditions = factor(colData(sce)$sample))
saveRDS(sce, "./fitg_sce.rds")
mean(rowData(sce)$tradeSeq$converged)
sce <- readRDS("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/fitg_sce.rds")

# test the null hypothesis that gene expression is not a function of pseudotime,
# i.e., whether the estimated smoothers are significantly varying as a function of pseudotime 
# within each lineage. The lineages=TRUE argument specifies that we would like the results for each lineage separately,
# asides from the default global test, which tests for significant associations across all lineages in the trajectory
# simultaneously. Further, we specify a log2 fold change cut-off to test against using the l2fc argument.

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))

library(tradeSeq)
assocRes <- rowData(sce)$assocRes
mockGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionMock, "fdr") <= 0.05)
]
infectedGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionInfected, "fdr") <= 0.05)
]

length(mockGenes)
length(infectedGenes)

UpSetR::upset(fromList(list(mock = mockGenes, Infected = infectedGenes)))

### based on mean smoother
yhatSmooth <- predictSmooth(sce, gene = mockGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)

yhatSmooth <- predictSmooth(sce, gene = infectedGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)

## the hierarchical trees constructed here, can also be used for 
## clustering of the genes according to their average expression pattern.
cl <- sort(cutree(heatSmooth$tree_row, k = 8))
table(cl)

conditions <- colData(sce)$sample
pt1 <- colData(sce)$slingPseudotime_1

### based on fitted values (plotting takes a while to run)
yhatCell <- predictCells(sce, gene=infectedGenes)
yhatCellMock <- yhatCell[,conditions == "Infected"]
# order according to pseudotime
ooMock <- order(pt1[conditions == "Infected"], decreasing=FALSE)
yhatCellMock <- yhatCellMock[,ooMock]
pheatmap(t(scale(t(yhatCellMock))), cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames=FALSE)

# Differential Analysis across conditions
condRes <- conditionTest(sce, l2fc = log2(1))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

# plot genes
oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
plotSmoothers(sce, assays(sce)$counts,
              gene = rownames(assays(sce)$counts)[oo[1]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce, assays(sce)$counts,
              gene = rownames(assays(sce)$counts)[oo[nrow(sce)]],
              alpha = 1, border = TRUE)

### based on mean smoother
yhatSmooth <- predictSmooth(sce, gene = conditionGenes, nPoints = 50, tidy = FALSE)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_TGF <- pheatmap(yhatSmoothScaled[, 51:100],
                           cluster_cols = FALSE,
                           show_rownames = FALSE, show_colnames = FALSE, main = "TGF-Beta", legend = FALSE,
                           silent = TRUE
)

matchingHeatmap_mock <- pheatmap(yhatSmoothScaled[heatSmooth_TGF$tree_row$order, 1:50],
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = FALSE, show_colnames = FALSE, main = "Mock",
                                 legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_TGF[[4]], matchingHeatmap_mock[[4]], ncol = 2)

######################### 
# https://bustools.github.io/BUS_notebooks_R/slingshot.html#differential_expression
library(tidymodels)
combined.all <- FindVariableFeatures(combined.all)
pal <- viridis(100, end = 0.95)
top_hvg <- HVFInfo(combined.all) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(variance)) %>% 
  top_n(300, variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(combined.all, slot = "data")[top_hvg,])
dat_use <- t(df[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sce)[,1], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)

model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)

val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)

summary(dat_use_df$pseudotime)

var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:4]
top_genes <- names(var_imp)[5:8]
top_genes <- names(var_imp)[9:12]
top_genes <- names(var_imp)[13:16]

par(mfrow = c(2, 2))
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(reducedDim(SlingshotDataSet(sce)), col = colors, 
       pch = 16, cex = 0.5, main = top_genes[i])
  lines(SlingshotDataSet(sce), lwd = 2, col = 'black', type = 'lineages')
}

##################### April 12 - slingshot
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#finding_differentially_expressed_genes
BiocParallel::register(BiocParallel::SerialParam())

# combined.all <- FindVariableFeatures(combined.all)
# Save the objects as separate matrices for input in slingshot
dimred <- combined.all@reductions$umap@cell.embeddings
clustering <- combined.all$main
counts <- as.matrix(combined.all@assays$RNA@counts[combined.all@assays$RNA@var.features, ])

# The fitting of GAMs can take quite a while, so for demonstration purposes we first do a very stringent filtering of the genes. 
# In an ideal experiment, you would use all the genes, or at least those defined as being variable.
lineages <- getLineages(data = dimred, clusterLabels = clustering)

# Plot the lineages            
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd = 3, col = "black")

curves <- getCurves(lineages, approx_points = 500, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

# library(tradeSeq)

# Removing some genes to speed up the computations for this tutorial
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
# 904 9700

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)

# Define function to plot
library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

# Genes that change with pseudotime

pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)


########################### Trajectory April c6
sce <- as.SingleCellExperiment(combined.all, assay = "RNA")
library(scater)
colLabels(sce) <- combined.all$main
by.cluster <- aggregateAcrossCells(sce, ids=colLabels(sce))
centroids <- reducedDim(by.cluster, "UMAP")

# Set clusters=NULL as we have already aggregated above.
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")

plotUMAP(sce, colour_by="main") + 
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="UMAP")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

common.pseudo <- averagePseudotime(tscan.pseudo) 
plotUMAP(sce, colour_by=I(common.pseudo), 
         text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

pseudo.all <- quickPseudotime(sce, use.dimred="TSNE")
head(pseudo.all$ordering)

pseudo.og <- quickPseudotime(sce, use.dimred="UMAP", outgroup=TRUE)
set.seed(10101)
plot(pseudo.og$mst)

pseudo.mnn <- quickPseudotime(sce, use.dimred="UMAP", with.mnn=TRUE)
mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)
plotUMAP(sce, colour_by=I(mnn.pseudo), text_by="label", text_colour="red") +
  geom_line(data=pseudo.mnn$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

library(slingshot)
sce.sling <- slingshot(sce, reducedDim='UMAP', clusterLabels = colData(sce)$integrated_snn_res.2,
                       start.clus = '12', approx_points = 300)
head(sce.sling$slingPseudotime_1)

embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

plotUMAP(sce.sling, colour_by="slingPseudotime_1") +
  geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2)

sce.sling2 <- slingshot(sce, cluster=colLabels(sce), reducedDim='UMAP')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)

sce <- runUMAP(sce, dimred="PCA")
reducedDim(sce.sling2, "UMAP") <- reducedDim(sce, "UMAP")

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=Dim.1, y=Dim.2), size=1.2)
}

gg

