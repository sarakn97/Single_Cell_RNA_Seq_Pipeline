#### scRNA RSV & Lactobacillus Study
#### FEBRUARY 2024
#### 4 Samples : Lactobacillus, RSV, Co-Infection, Mock
#### Leyao Wang Lab

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
#library(DoubletFinder)

# set working directory
setwd("/storage1/fs1/leyao.wang/Active/scRNA_Feb2024/")
combined.all <- readRDS("./combined_all.rds")

######################### RSV
# read in sample
RSV = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Feb2024/SR003735_RSV_AO_scRNA_RSV_72/outs/filtered_feature_bc_matrix/")
dim(RSV)
# 36601 5427

# Create Sample1 Seurat Object (Infected)
RSV <- CreateSeuratObject(counts = RSV, project = "RSV", min.cells = 3, min.features = 200)
dim(RSV)
# 27094 5425
head(RSV)

######################### CoInf
# read in sample
CoInf = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Feb2024/SR003735_Lact_RSV_AO_scRNA_LacRSV_72/outs/filtered_feature_bc_matrix/")
dim(CoInf)
# 36601 5427

# Create Sample1 Seurat Object (Infected)
CoInf <- CreateSeuratObject(counts = CoInf, project = "CoInf", min.cells = 3, min.features = 200)
dim(CoInf)
# 27680 5097
head(CoInf)

######################### LACT
# read in sample
LACT = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Feb2024/SR003735_Lact_AO_scRNA_Lac_72/outs/filtered_feature_bc_matrix/")
dim(LACT)
# 36601 5427

# Create Sample1 Seurat Object (Infected)
LACT <- CreateSeuratObject(counts = LACT, project = "LACT", min.cells = 3, min.features = 200)
dim(LACT)
# 27680 5097
head(LACT)

######################### Mock
# read in sample
Mock = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Feb2024/SR003735_Mock_AO_scRNA_Mock_72/outs/filtered_feature_bc_matrix/")
dim(Mock)
# 36601 12031

# Create Sample1 Seurat Object (Infected)
Mock <- CreateSeuratObject(counts = Mock, project = "Mock", min.cells = 3, min.features = 200)
dim(Mock)
head(Mock)

######################################################################################
# Identify Condition
Mock[["sample"]] <- "Mock"
RSV[["sample"]] <- "RSV"
CoInf[["sample"]] <- "CoInfection"
LACT[["sample"]] <- "Lactobacillus"

# Merge Samples
merged_seurat <- merge(x=Mock, y= c(RSV, CoInf, LACT))

sum(grepl("MT",rownames(merged_seurat)))
# 330
rownames(merged_seurat)[grepl("MT-",rownames(merged_seurat))]

# calculate the proportion of transcripts mapping to mitochondrial genes
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "MT-")

max(merged_seurat@meta.data$percent.mt)
# 90.1974

# Plot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue", "green", "purple"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue", "green", "purple"))

# Filter based on quality plots 
Filtered <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt > -Inf & percent.mt <20)
VlnPlot(Filtered, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filtered, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue", "green", "orange"))
FeatureScatter(object=Filtered, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue", "green", "orange"))

# Sums post Filtering
sum(Filtered$sample=="CoInfection")
# 10116
sum(Filtered$sample=="RSV")
# 5018
sum(Filtered$sample=="Mock")
# 11490
sum(Filtered$sample=="Lactobacillus")
# 7782

Filtered@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",Filtered@assays$RNA@counts@Dimnames[[1]])
Filtered@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filtered@assays$RNA@counts@Dimnames[[1]])

split_seurat <- SplitObject(Filtered, split.by = "sample")
split_seurat <- split_seurat[c("Mock", "RSV", "CoInfection", "Lactobacillus")]

######################### DOUBLET FINDER #######################################
# Pre-Process Seurat Object
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
  split_seurat[[i]] <- RunPCA(split_seurat[[i]])
  split_seurat[[i]] <- RunUMAP(split_seurat[[i]], dims=1:20)
}

# Doublet Finder Infected
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- split_seurat[[i]] %>%
    FindNeighbors(dims = 1:20)%>%
    FindClusters(resolution=0.6)
}

saveRDS(split_seurat, "./split_seurat_march2024.rds")
################### Find Doublets in Mock
# (pk identification) 
sweep.res.list_1 <- paramSweep_v3(split_seurat[[1]], PCs= 1:10, sct=T)
sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
bcmvn_1 <- find.pK(sweep.stats_1)
#pK = 0.05

###### Add Labels to pK plot
pK=as.numeric(as.character(bcmvn_1$pK))
BCmetric=bcmvn_1$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[1]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.086*nrow(split_seurat[[1]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.28, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
colnames(split_seurat[[1]]@meta.data)

# Plot Doublets & Singlets
DimPlot(split_seurat[[1]], reduction = "umap", group.by = c("DF.classifications_0.25_0.04_185"),label = TRUE,repel = TRUE)

saveRDS(split_seurat, "./split_seurat_march2024.rds")
################### Find Doublets in RSV
# (pk identification) 
sweep.res.list_2 <- paramSweep_v3(split_seurat[[2]], PCs= 1:10, sct=T)
sweep.stats_2 <- summarizeSweep(sweep.res.list_2, GT=F)
bcmvn_2 <- find.pK(sweep.stats_2)
#pK = 0.05

###### Add Labels to pK plot
pK=as.numeric(as.character(bcmvn_2$pK))
BCmetric=bcmvn_2$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
######

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[2]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*nrow(split_seurat[[2]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
colnames(split_seurat[[2]]@meta.data)

# Plot Doublets & Singlets
DimPlot(split_seurat[[2]], reduction = "umap", group.by = c("DF.classifications_0.25_0.03_196"),label = TRUE,repel = TRUE)

saveRDS(split_seurat, "./split_seurat_march2024.rds")
################### Find Doublets in Co-Infection
# (pk identification) 
sweep.res.list_3 <- paramSweep_v3(split_seurat[[3]], PCs= 1:10, sct=T)
sweep.stats_3 <- summarizeSweep(sweep.res.list_3, GT=F)
bcmvn_3 <- find.pK(sweep.stats_3)
#pK = 0.05

###### Add Labels to pK plot
pK=as.numeric(as.character(bcmvn_3$pK))
BCmetric=bcmvn_3$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
######

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[3]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.077*nrow(split_seurat[[3]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[3]] <- doubletFinder_v3(split_seurat[[3]], PCs=1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
colnames(split_seurat[[3]]@meta.data)

# Plot Doublets & Singlets
DimPlot(split_seurat[[3]], reduction = "umap", group.by = c("DF.classifications_0.25_0.06_779"),label = TRUE,repel = TRUE)

saveRDS(split_seurat, "./split_seurat_march2024.rds")
################### Find Doublets in Lactobacillus
# (pk identification) 
sweep.res.list_4 <- paramSweep_v3(split_seurat[[4]], PCs= 1:10, sct=T)
sweep.stats_4 <- summarizeSweep(sweep.res.list_4, GT=F)
bcmvn_4 <- find.pK(sweep.stats_4)
#pK = 0.05

###### Add Labels to pK plot
pK=as.numeric(as.character(bcmvn_4$pK))
BCmetric=bcmvn_4$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
######

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[4]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0595*nrow(split_seurat[[4]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[4]] <- doubletFinder_v3(split_seurat[[4]], PCs=1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
colnames(split_seurat[[4]]@meta.data)

# Plot Doublets & Singlets
DimPlot(split_seurat[[4]], reduction = "umap", group.by = c("DF.classifications_0.25_0.08_463"),label = TRUE,repel = TRUE)
saveRDS(split_seurat, "./split_seurat_march2024.rds")
####################################### REMOVE DOUBLETS #######################
split_seurat <- readRDS("./split_seurat_march2024.rds")

split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.28_988!="Doublet") 
dim(split_seurat[[1]]) # 26415 10502

split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset = DF.classifications_0.25_0.03_196!="Doublet") 
dim(split_seurat[[2]]) # 25498 4822

split_seurat[[3]] <- subset(x = split_seurat[[3]], 
                             subset =DF.classifications_0.25_0.06_779!="Doublet") 
dim(split_seurat[[3]]) # 26466 9337

split_seurat[[4]] <- subset(x = split_seurat[[4]], 
                            subset = DF.classifications_0.25_0.08_463!="Doublet") 
dim(split_seurat[[4]]) # 26351 7319

saveRDS(split_seurat, "./split_seurat_march2024.rds") # pre-down-sample

########################################################## Downsample ##########
# set.seed(111)
# split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 3700, replace=F)]
# split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 3700, replace=F)]
# # split_seurat[[3]] <- split_seurat[[3]][, sample(colnames(split_seurat[[3]]), size = 3700, replace=F)]
# 
# for (i in 1:length(split_seurat)){
#   split_seurat[[i]] <- RunTSNE(split_seurat[[i]], dims = 1:20)
# }
############################ Integration ##########################################
# select most variable features
split_seurat <- readRDS("./split_seurat_march2024.rds")
integ_features <- SelectIntegrationFeatures(object.list=split_seurat, nfeatures=5000)

# prepare SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list=split_seurat, anchor.features=integ_features)

# Find Anchors
integ_anchors <- FindIntegrationAnchors(object.list=split_seurat, normalization.method = "SCT",
                                        anchor.features=integ_features)

# Integrate
combined.all <- IntegrateData(anchorset=integ_anchors, normalization.method="SCT")

saveRDS(combined.all, "./Integrated_combinedall_noDownSample.rds") 
combined.all <- readRDS("./Integrated_combinedall_noDownSample.rds")

DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) # reduction = "pca",
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = c(0.8))
combined.all <- FindVariableFeatures(combined.all, selection.method = "vst", nfeatures = 2000)
combined.all[["percent.mt"]] <- PercentageFeatureSet(combined.all, pattern = "MT-",assay = "RNA")

combined.all <- NormalizeData(combined.all, assay = "RNA")
combined.all <- ScaleData(combined.all, assay = "RNA")

DefaultAssay(combined.all) <- "RNA"
DimPlot(combined.all, reduction = "umap", label = T,repel = TRUE)
########## UMAP plots
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
combined.all <- SetIdent(combined.all, value = combined.all$integrated_snn_res.0.6)
DimPlot(combined.all, reduction = "umap", label = T,repel = TRUE)

combined.all <- SetIdent(combined.all, value = "sample")
DimPlot(combined.all, reduction = "umap", label = F,repel = TRUE)
######## rename genes so Clusters can be recognized
### combined.all
combined.all@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$RNA@counts@Dimnames[[1]])
combined.all@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$RNA@counts@Dimnames[[1]])

combined.all@assays$SCT@counts@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$SCT@counts@Dimnames[[1]])
combined.all@assays$SCT@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$SCT@data@Dimnames[[1]])
rownames(combined.all@assays$SCT@scale.data) <- gsub("GRCh38-","",rownames(combined.all@assays$SCT@scale.data))

combined.all@assays$integrated@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$integrated@data@Dimnames[[1]])
rownames(combined.all@assays$integrated@scale.data) <- gsub("GRCh38-","",rownames(combined.all@assays$integrated@scale.data))

saveRDS(combined.all, "./combined_all.rds") 
combined.all <- readRDS("./combined_all.rds")
################ BEGIN ANNOTATIONS + EXPLORATION ################################
cluster_marks <- FindAllMarkers(combined.all)
write.csv(cluster_marks, "./cluster_markers_res0.3.csv")
cluster_marks <- read.csv("./cluster_markers_res0.3.csv")

cilia <- WhichCells(combined.all, idents = c(9,10))
secretory <- WhichCells(combined.all, idents = c(6,7,11))
basal <- WhichCells(combined.all, idents = c(0,4))
prolif_basal <- WhichCells(combined.all, idents = c(12))
suprabasal <- WhichCells(combined.all, idents = c(1,2,5,8))
idk <- WhichCells(combined.all, idents = c(3))

combined.all@meta.data[cilia, "annotated"] <- "Ciliated"
combined.all@meta.data[secretory, "annotated"] <- "Secretory"
combined.all@meta.data[basal, "annotated"] <- "Basal"
combined.all@meta.data[suprabasal, "annotated"] <- "Suprabasal"
combined.all@meta.data[prolif_basal, "annotated"] <- "Proliferating-Basal"
combined.all@meta.data[idk, "annotated"] <- "Unknown"

combined.all <- SetIdent(combined.all, value = "annotated")
DimPlot(combined.all, label = T)

################################### SingleR Annotations ########################
library(SingleR)
library(celldex)

# establish reference
ref <- BlueprintEncodeData()
hpca.se <- HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(combined.all)
combined.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)
combined.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)
table(combined.main$pruned.labels)
table(combined.fine$pruned.labels)
combined.all@meta.data$combined.main <- combined.main$pruned.labels
combined.all@meta.data$combined.fine <- combined.fine$pruned.labels
combined.all <- SetIdent(combined.all, value = "combined.main")
DimPlot(combined.all, label = T , repel = T, label.size = 3) + NoLegend()
combined.all <- SetIdent(combined.all, value = "combined.fine")
DimPlot(combined.all, label = T , repel = T, label.size = 3) + NoLegend()
combined.all <-  subset(combined.all, combined.fine != "")

saveRDS(combined.all, "./combined_all.rds") 

###################################### Manual Annotation #######################
FeaturePlot(combined.all, features = c("SCGB3A1", "TUBA1A", "S100A2", "TP63", "KRT5", "HELLS", "CFTR", "MKI67", "MUC1"), label = T)
FeaturePlot(combined.all, features = c("RNASE1", "MUC1", "SCGB1A1", "TFF3", "POU2F3", "TRPM5", "TUBB4B", "PIFO", "CXCL14"), label = T)
FeaturePlot(combined.all, features = c("TUBB4B", "PIFO", "CXCL14"), label = T)
FeaturePlot(combined.all, features = c("HES1", "SERPINB3", "TPPP3", "FOXJ1", "HES6"), label = T)
FeaturePlot(combined.all, features = c("ASCL3", "RGS13", "GNAT3", "LRMP", "MYB", "ASCL1", "FOXI1"), label = T)
FeaturePlot(combined.all, features = c("KRT15", "KRT17", "FXYD2", "TMPRSS11E"), label = T)

genes <- c("KRT5", "KRT14", "TP63", "DLK2", "SERPINB4", "KRT19", "NOTCH3", "FABP5", "TFCP2L1", "S100A9", 
           "SCCB1A1", "TSPAN8", "MUC5AC", "MUC1", "S100P", "PSCA", "MES6", "DEUP1", "FOXJ1","PIFO", "TUBB4B",
           "RYR3", "TRPXL", "PALLD", "CYP24A1", "LTF", "LYZ", "PIP", "AZGP1", "MUC5B", "AQP1", "GNG11",
           "BPIFB2", "ACKR1", "FBLN1", "DCN", "C1R", "DES", "CNN1", "ACTA2", "MYL9", "RERGL", "PDGFRB", 
           "TYROBP", "APOC1", "C1QA", "SDS", "FNIP2", "CST3", "MS4A6A", "HLADPB1", "TPSAB1", "CPA3", "HPGDS",
           "IL32", "CD3D", "CCL5", "LTB", "MS4A1", "CD79A", "IGJ", "SSR4", "MZB1", "TOP2A", "MKI67")

combined.all <- SetIdent(combined.all, value = combined.all$integrated_snn_res.0.6)
DoHeatmap(combined.all, features = genes, assay = "RNA", angle = 90, slot="data") + NoLegend()

p = VlnPlot(combined.all,genes,
            pt.size=0,stack=T,flip=T,fill.by='ident',assay='RNA')+
  NoLegend()+
  scale_x_discrete(position='top')+
  theme(panel.background=element_rect(color='black',size=0.8),
        panel.spacing=unit(0,'lines'),
        strip.text.y.right=element_text(size=12,hjust=0),
        axis.text.x=element_text(size=12,angle=45,hjust=0, face= "bold"),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=element_line(size=0.8),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #plot.margin=margin(2,2,2,2,'pt'),
        plot.title=element_blank())+
  coord_cartesian(clip='off')

####################################################### Differential Expression
combined.all <- SetIdent(combined.all, value = "sample")
rsvVmock <- FindMarkers(combined.all, ident.1 = "RSV", ident.2 = "Mock")
CoinfVmock <- FindMarkers(combined.all, ident.1 = "CoInfection", ident.2 = "Mock")
lacVmock <- FindMarkers(combined.all, ident.1 = "Lactobacillus", ident.2 = "Mock")
coVrsv <- FindMarkers(combined.all, ident.1 = "CoInfection", ident.2 = "RSV")
rsvVlac <- FindMarkers(combined.all, ident.1 = "RSV", ident.2 = "Lactobacillus")

write.csv(rsvVmock, "./rsv_vs_mock.csv")
write.csv(CoinfVmock, "./co_vs_mock.csv")
write.csv(lacVmock, "./lac_vs_mock.csv")
write.csv(coVrsv, "./co_vs_rsv.csv")
write.csv(rsvVlac, "./rsv_vs_lac.csv")

rsvVlac_marks <- read.csv("./rsv_vs_lac.csv")
rsvVmock_marks <- read.csv("./rsv_vs_mock.csv")
coVmock_marks <- read.csv("./co_vs_mock.csv")
lacVmock_marks <- read.csv("./lac_vs_mock.csv")
coVrsv_marks <- read.csv("./co_vs_rsv.csv")

row.names(rsvVlac_marks) <- rsvVlac_marks$X
row.names(rsvVmock_marks) <- rsvVmock_marks$X
row.names(coVmock_marks) <- coVmock_marks$X
row.names(lacVmock_marks) <- lacVmock_marks$X
row.names(coVrsv_marks) <- coVrsv_marks$X

# Significant Genes
sum(rsvVlac_marks$p_val_adj  < 0.05) # 2930   # 260
sum(coVrsv_marks$p_val_adj  < 0.05) # 3472    # 297
sum(rsvVmock_marks$p_val_adj  < 0.05) # 1203  # 264
sum(coVmock_marks$p_val_adj  < 0.05) # 3327   # 759
sum(lacVmock_marks$p_val_adj  < 0.05) # 3275  # 216

sum(rsvVlac_marks$p_val_adj  < 0.05 & rsvVlac_marks$avg_log2FC > 0) # 170
sum(rsvVlac_marks$p_val_adj  < 0.05 & rsvVlac_marks$avg_log2FC < 0) # 90

sum(coVrsv_marks$p_val_adj  < 0.05 & coVrsv_marks$avg_log2FC > 0) # 173
sum(coVrsv_marks$p_val_adj  < 0.05 & coVrsv_marks$avg_log2FC < 0) # 124

sum(rsvVmock_marks$p_val_adj  < 0.05 & rsvVmock_marks$avg_log2FC > 0) # 172
sum(rsvVmock_marks$p_val_adj  < 0.05 & rsvVmock_marks$avg_log2FC < 0) # 92

sum(coVmock_marks$p_val_adj  < 0.05 & coVmock_marks$avg_log2FC > 0) # 487
sum(coVmock_marks$p_val_adj  < 0.05 & coVmock_marks$avg_log2FC < 0) # 272

sum(lacVmock_marks$p_val_adj  < 0.05 & lacVmock_marks$avg_log2FC > 0) # 137
sum(lacVmock_marks$p_val_adj  < 0.05 & lacVmock_marks$avg_log2FC < 0) # 79


##################### VOLCANO PLOTS ############################################
combined.all <- SetIdent(combined.all, value = "sample")
comock <- subset(combined.all, idents = c("CoInfection", "Mock"))
library(SCpubr)
SCpubr::do_VolcanoPlot(sample = comock,
                       de_genes = coVmock_marks,
                       n_genes = 20,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)
rsvmock <- subset(combined.all, idents = c("RSV", "Mock"))
SCpubr::do_VolcanoPlot(sample = rsvmock,
                       de_genes = rsvVmock_marks,
                       n_genes = 20,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)
lacmock <- subset(combined.all, idents = c("Lactobacillus", "Mock"))
SCpubr::do_VolcanoPlot(sample = lacmock,
                       de_genes = lacVmock_marks,
                       n_genes = 20,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)

corsv <- subset(combined.all, idents = c("CoInfection", "RSV"))
SCpubr::do_VolcanoPlot(sample = corsv,
                       de_genes = coVrsv_marks,
                       n_genes = 20,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)
rsvlac <- subset(combined.all, idents = c("Lactobacillus", "RSV"))
SCpubr::do_VolcanoPlot(sample = rsvlac,
                       de_genes = rsvVlac_marks,
                       n_genes = 20,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)

############################ GSEA ##############################################
library(dplyr)
library(tibble)
library(ggplot2)
library(pathview)
library(gage)
library(gageData)
library(annotate)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

#### GAGE
coVmock_marks$entrez = mapIds(org.Hs.eg.db, keys=coVmock_marks$X, column="ENTREZID", keytype="SYMBOL", multiVals="first")
rsvVmock_marks$entrez = mapIds(org.Hs.eg.db, keys=rsvVmock_marks$X, column="ENTREZID", keytype="SYMBOL", multiVals="first")
lacVmock_marks$entrez = mapIds(org.Hs.eg.db, keys=lacVmock_marks$X, column="ENTREZID", keytype="SYMBOL", multiVals="first")
rsvVlac_marks$entrez = mapIds(org.Hs.eg.db, keys=rsvVlac_marks$X, column="ENTREZID", keytype="SYMBOL", multiVals="first")
coVrsv_marks$entrez = mapIds(org.Hs.eg.db, keys=coVrsv_marks$X, column="ENTREZID", keytype="SYMBOL", multiVals="first")

# Prepare Input for GAGE - all log fold change values and entrez IDs
fc_comock <- coVmock_marks$avg_log2FC
names(fc_comock) <- coVmock_marks$entrez
fc_rsvmock <- rsvVmock_marks$avg_log2FC
names(fc_rsvmock) <- rsvVmock_marks$entrez
fc_lacmock <- lacVmock_marks$avg_log2FC
names(fc_lacmock) <- lacVmock_marks$entrez
fc_corsv <- coVrsv_marks$avg_log2FC
names(fc_corsv) <- coVrsv_marks$entrez
fc_rsvlac <- rsvVlac_marks$avg_log2FC
names(fc_rsvlac) <- rsvVlac_marks$entrez

# GO - Gene Ontology Database
data(go.sets.hs)
data(go.subs.hs)
# GO - Biological Process
gobpsets = go.sets.hs[go.subs.hs$BP]
gobp_comock = gage(fc_comock, gsets=gobpsets, same.dir=TRUE)
gobp_rsvmock = gage(fc_rsvmock, gsets=gobpsets, same.dir=TRUE)
gobp_lacmock = gage(fc_lacmock, gsets=gobpsets, same.dir=TRUE)
gobp_rsvlac = gage(fc_rsvlac, gsets=gobpsets, same.dir=TRUE)
gobp_corsv = gage(fc_corsv, gsets=gobpsets, same.dir=TRUE)

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msig_comock = gage(fc_comock, gsets=msigdb, same.dir=TRUE)
msig_rsvmock = gage(fc_rsvmock, gsets=msigdb, same.dir=TRUE)
msig_lacmock = gage(fc_lacmock, gsets=msigdb, same.dir=TRUE)
msig_rsvlac = gage(fc_rsvlac, gsets=msigdb, same.dir=TRUE)
msig_corsv = gage(fc_corsv, gsets=msigdb, same.dir=TRUE)

##### final bar plot up-regulated
################################ CO V MOCK ######################################
# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_comock$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
which(msig_comock$greater[,4] < 0.05)
par(mar=c(3,24,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:2]), col="red", horiz=TRUE, cex.names=1, las =2, main = "MSigDB Up-Regulated Biological Processes (Co-Inf v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_comock$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_comock$greater[,4] < 0.05))
par(mar=c(3,24,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:10]), col="red", horiz=TRUE, cex.names=1, las =2, main = "GO Up-Regulated Biological Processes (Co-Inf v. Mock)", xlab = "-log(adj.p.val)")

# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_comock$less[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
which(msig_comock$less[,4] < 0.05)
par(mar=c(3,24,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:2]), col="blue", horiz=TRUE, cex.names=1, las =2, main = "MSigDB Down-Regulated Biological Processes (Co-Inf v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_comock$less[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_comock$less[,4] < 0.05))
par(mar=c(3,35,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:10]), col="blue", horiz=TRUE, cex.names=1, las =2, main = "GO Down-Regulated Biological Processes (Co-Inf v. Mock)", xlab = "-log(adj.p.val)")

################################ RSV V MOCK #####################################
# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_rsvmock$greater[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
sum(which(msig_rsvmock$greater[,4] < 0.05))
par(mar=c(5,22,4,7)+.1)
barplot(neg.p.sort_spl[1:3], col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Up-Regulated Biological Processes (RSV v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_rsvmock$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_rsvmock$greater[,4] < 0.05))
par(mar=c(3,23,3,10)+.1)
barplot(neg.p.sort_lung[1:10], col="red", horiz=TRUE, cex.names=1, las =2, main = "GO Up-Regulated Biological Processes (RSV v. Mock)", xlab = "-log(adj.p.val)")

# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_rsvmock$less[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
sum(which(msig_rsvmock$less[,4] < 0.05)) 
par(mar=c(5,22,4,7)+.1)
barplot(neg.p.sort_spl[1:7], col="blue", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Down-Regulated Biological Processes (RSV v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_rsvmock$less[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
sum(which(gobp_rsvmock$less[,4] < 0.05))
par(mar=c(5,32,4,7)+.1)
barplot(neg.p.sort_spl[1:15], col="blue", horiz=TRUE, cex.names=0.9, las =2, main = "GO Down-Regulated Biological Processes (RSV v. Mock)", xlab = "-log(adj.p.val)")

################################ LAC V MOCK #####################################
# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_lacmock$greater[,4])
neg.p.sort_blood <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_blood)
which(msig_lacmock$greater[,4] < 0.05)
par(mar=c(5,19,4,7)+.1)
barplot(neg.p.sort_blood[1], col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Up-Regulated Biological Processes (Lact v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_lacmock$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_lacmock$greater[,4] < 0.05))
par(mar=c(3,27,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:15]), col="red", horiz=TRUE, cex.names=1, las =2, main = "GO Up-Regulated Biological Processes (Lact v. Mock)", xlab = "-log(adj.p.val)")

# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_lacmock$less[,4])
neg.p.sort_blood <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_blood)
which(msig_lacmock$less[,4] < 0.05)
par(mar=c(5,18,4,7)+.1)
barplot(neg.p.sort_blood[1:10], col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Down-Regulated Biological Processes (Lact v. Mock)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_lacmock$less[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_lacmock$less[,4] < 0.05))
par(mar=c(3,33,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:15]), col="blue", horiz=TRUE, cex.names=1, las =2, main = "GO Down-Regulated Biological Processes (Lact v. Mock)", xlab = "-log(adj.p.val)")


################################ CO V RSV ######################################
# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_corsv$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(msig_corsv$greater[,4] < 0.05))
par(mar=c(3,24,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:3]), col="red", horiz=TRUE, cex.names=1, las =2, main = "MSigDB Up-Regulated Biological Processes (Co-Inf v. RSV)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_corsv$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
which(gobp_corsv$greater[,4] < 0.05)
par(mar=c(3,24,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:15]), col="red", horiz=TRUE, cex.names=1, las =2, main = "GO Up-Regulated Biological Processes (Co-Inf v. RSV)", xlab = "-log(adj.p.val)")

# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_corsv$less[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
which(msig_corsv$less[,4] < 0.05)
par(mar=c(3,24,3,10)+.1)
barplot(neg.p.sort_lung[1:1], col="blue", horiz=TRUE, cex.names=1, las =2, main = "MSigDB Down-Regulated Biological Processes (Co-Inf v. RSV)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_corsv$less[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_corsv$less[,4] < 0.05))
par(mar=c(3,33,3,10)+.1)
barplot(neg.p.sort_lung[1:15], col="blue", horiz=TRUE, cex.names=1, las =2, main = "GO Down-Regulated Biological Processes (Co-Inf v. RSV)", xlab = "-log(adj.p.val)")

################################ RSV V LAC #####################################
# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_rsvlac$greater[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
which(msig_rsvlac$greater[,4] < 0.05) 
par(mar=c(5,22,4,7)+.1)
barplot(neg.p.sort_spl[1:2], col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Up-Regulated Biological Processes (RSV v. Lactobacillus)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_rsvlac$greater[,4])
neg.p.sort_lung <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_lung)
sum(which(gobp_rsvlac$greater[,4] < 0.05))
par(mar=c(3,27,3,10)+.1)
barplot(rev(neg.p.sort_lung[1:10]), col="red", horiz=TRUE, cex.names=1, las =2, main = "GO Up-Regulated Biological Processes (RSV v. Lactobacillus)", xlab = "-log(adj.p.val)")

# MSigDB Pathway Analysis Pathways
neg.pval <-  -log(msig_rsvlac$less[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
which(msig_rsvlac$less[,4] < 0.05) 
par(mar=c(5,27,4,7)+.1)
barplot(neg.p.sort_spl[1:1], col="blue", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Down-Regulated Biological Processes (RSV v. Lactobacillus)", xlab = "-log(adj.p.val)")

# GO Pathway Analysis Pathways
neg.pval <-  -log(gobp_rsvlac$less[,4])
neg.p.sort_spl <- sort(neg.pval, decreasing =  T)
head(neg.p.sort_spl)
sum(which(gobp_rsvlac$less[,4] < 0.05))
par(mar=c(5,26,4,7)+.1)
barplot(neg.p.sort_spl[1:15], col="blue", horiz=TRUE, cex.names=0.9, las =2, main = "GO Down-Regulated Biological Processes (RSV v. Lactobacillus)", xlab = "-log(adj.p.val)")


########################## COMBINED PLOT #######################################
df_l <- data.frame(neg.p.sort_lung)
df_s <- data.frame(neg.p.sort_spl)
df_b <- data.frame(neg.p.sort_blood)
df_l$gs <- row.names(df_l)
df_s$gs <- row.names(df_s)
df_b$gs <- row.names(df_b)
df_l$cell <- "Lung"
df_s$cell <- "Spleen"
df_b$cell <- "Blood"
df_l$c <- "blue"
df_s$c <- "green"
df_b$c <- "red"
colnames(df_l) <- c("neg.p", "gs","cell", "col"); colnames(df_s) <- c("neg.p", "gs", "cell", "col"); colnames(df_b) <- c("neg.p", "gs", "cell", "col")
df_l <- df_l[1:2,]; df_s <- df_s[1:2,]; df_b <- df_b[1:2,]
df <- rbind(df_l, df_s, df_b)

par(mar=c(5,22,5,1)+.1)
barplot(df$neg.p, col=df$col, horiz=TRUE, cex.names=1, cex.axis = 1.5, las =2, font.lab =2, cex.lab = 1.5,main = "MSigDB Up-Regulated Biological Processes", xlab = "-log(adj.p.val)", names.arg = df$gs, font = 2, font.axis = 2) + legend("topright",
                                                                                                                                                                                                                                            legend = c(as.expression(bquote(bold("Blood"))),as.expression(bquote(bold("Spleen"))), as.expression(bquote(bold("Lung")))),
                                                                                                                                                                                                                                            fill = c("red", 'green', "blue"),    
                                                                                                                                                                                                                                            border = "black",
                                                                                                                                                                                                                                            cex = 1.1) 

################################################################################
#################### PROPORTIONS PLOT ##########################################
DefaultAssay(rsvmock) <- "integrated"
rsvmock <- RunUMAP(rsvmock, dims = 1:20) 
rsvmock <- FindNeighbors(rsvmock, reduction = "pca", dims = 1:20)
rsvmock <- FindClusters(rsvmock, resolution = c(0.5, 0.8))
DefaultAssay(rsvmock) <- "RNA"
DimPlot(rsvmock, label = T, pt.size =0.7)
rsvmock <- NormalizeData(rsvmock, assay = "RNA")
rsvmock <- ScaleData(rsvmock, assay = "RNA")

rsvmock <- SetIdent(rsvmock, value = "integrated_snn_res.0.8")
prop = prop.table(table(Idents(rsvmock),rsvmock$sample),1) %>% data.frame()
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Mock',],-prop)$cluster)
prop$cluster = factor(prop$cluster, levels = ranked)
counts = table(Idents(rsvmock)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster, levels = ranked)
prop$group <- factor(prop$group, levels=c('RSV', 'Mock'))

prop_plot = ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='RSV',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='Mock',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,1.3))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('red3','royalblue4'),labels=c('RSV','Mock'), breaks = c("RSV", "Mock"))+
  labs(y='Proportion (%)', x= "Cluster")+
  theme(panel.background=element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=16,hjust=0.5,color='black'),
        axis.text.x=element_text(size=16,hjust=1,color='black',angle=45,vjust=1, face = "bold"),
        legend.position='right',
        legend.text=element_text(size=18,color='black', face = "bold"),
        legend.title=element_blank(),
        axis.title.y=element_text(size=18,hjust=0.2,color='black',margin=margin(r=10), face = "bold"),
        axis.title.x=element_text(size=18,color='black', face = "bold"), plot.title = element_text(size=20, face = "bold"))
prop_plot

prop = as.data.frame(table(Idents(rsvmock),rsvmock$sample))
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Control',],-prop)$cluster)
prop$cluster = factor(prop$cluster,levels=ranked)
counts = table(Idents(rsvmock)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster,levels=ranked)

# Cell Count instead of Proportion
prop_plot = ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='Infected',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='Control',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,1))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('royalblue4','red3'),labels=c('Control','Infected'))+
  labs(y='Cell Count')+
  theme(panel.background=element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=15,hjust=0.5,color='black'),
        axis.text.x=element_text(size=18,hjust=1,color='black',angle=45,vjust=1, face = "bold"),
        legend.position='top',
        legend.text=element_text(size=18,color='black', face = "bold"),
        legend.title=element_blank(),
        axis.title.y=element_text(size=18,hjust=0.2,color='black',margin=margin(r=10), face = "bold"),
        axis.title.x=element_blank())
prop_plot

######## Proportions Plot
library(ggsci)
library(ggbreak)
df <- as.data.frame(table(combined.all$annotated, combined.all$sample))
df$Var1 <- factor(df$Var1, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
p <- ggplot(df) +
  aes(x = Var1, fill = Var2, weight = Freq, by = Var2) +
  geom_bar(position = "fill") +xlab("Cluster") + ylab("Proportion") 
p +scale_fill_lancet() + theme(  axis.text.y=element_text(size=16,color='black', face = "bold"),
            axis.text.x=element_text(size=14,color='black', face = "bold"),
            legend.position='right',
            legend.text=element_text(size=14,color='black', face = "bold"),
            legend.title=element_blank(),
            axis.title.y=element_text(size=20,color='black',face = "bold",angle = 90), axis.title.x=element_text(size=20,color='black',face = "bold"))

df$Var1 <- factor(df$Var1, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
ggplot(df) +aes(x = Var1, weight = Freq, fill=Var2) + geom_bar() +
  theme(strip.text = element_text(size = 10, face = "bold")) + scale_y_break(c(500,1000), scales=2) +
  xlab("Cluster") + ylab("proportion") + ylab("Cell Count") +scale_fill_lancet() + theme(  axis.text.y=element_text(size=16,color='black', face = "bold"),
                                                                                       axis.text.x=element_text(size=14,color='black', face = "bold"),
                                                                                       legend.position='right',
                                                                                       legend.text=element_text(size=14,color='black', face = "bold"),
                                                                                       legend.title=element_blank(),
                                                                                       axis.title.y=element_text(size=20,color='black',face = "bold",angle = 90), axis.title.x=element_text(size=20,color='black',face = "bold"))
#  angle =35, hjust =1, vjust = 1,

df <- as.data.frame(table(rsvmock$annotated, rsvmock$sample))
p <- ggplot(df) +
  aes(x = Var2, fill = Var1, weight = Freq, by = Var1) +
  geom_bar(position = "fill") +xlab("Cluster") + ylab("Proportion") 
p +scale_fill_lancet() + theme(  axis.text.y=element_text(size=16,color='black', face = "bold"),
                                 axis.text.x=element_text(size=14,color='black', face = "bold", angle = 35, hjust = 1),
                                 legend.position='right',
                                 legend.text=element_text(size=14,color='black', face = "bold"),
                                 legend.title=element_blank(),
                                 axis.title.y=element_text(size=20,color='black',face = "bold",angle = 90), axis.title.x=element_text(size=20,color='black',face = "bold"))


################################################################################
##### Differential Expression by Cell-Type #####################################
ciliated <- subset(combined.all, idents = "Ciliated")
DefaultAssay(ciliated) <- "integrated"
ciliated <- RunUMAP(ciliated, dims = 1:20) 
ciliated <- FindNeighbors(ciliated, reduction = "pca", dims = 1:20)
ciliated <- FindClusters(ciliated, resolution = c(0.5, 0.8))
DefaultAssay(ciliated) <- "RNA"
DimPlot(ciliated, label = F, pt.size =0.7)
ciliated <- NormalizeData(ciliated, assay = "RNA")
ciliated <- ScaleData(ciliated, assay = "RNA")
saveRDS(ciliated, "ciliated_cells.rds")


ciliated <- SetIdent(ciliated, value = "sample")
cilia_rsvVmock <- FindMarkers(ciliated, ident.1 = "RSV", ident.2 = "Mock")
cilia_CoinfVmock <- FindMarkers(ciliated, ident.1 = "CoInfection", ident.2 = "Mock")
cilia_lacVmock <- FindMarkers(ciliated, ident.1 = "Lactobacillus", ident.2 = "Mock")
cilia_coVrsv <- FindMarkers(ciliated, ident.1 = "CoInfection", ident.2 = "RSV")
cilia_rsvVlac <- FindMarkers(ciliated, ident.1 = "RSV", ident.2 = "Lactobacillus")

write.csv(cilia_rsvVmock, "./cilia_rsv_vs_mock.csv")
write.csv(cilia_CoinfVmock, "./cilia_co_vs_mock.csv")
write.csv(cilia_lacVmock, "./cilia_lac_vs_mock.csv")
write.csv(cilia_coVrsv, "./cilia_co_vs_rsv.csv")
write.csv(cilia_rsvVlac, "./cilia_rsv_vs_lac.csv")

rsvVlac_marks <- read.csv("./cilia_rsv_vs_lac.csv")
rsvVmock_marks <- read.csv("./cilia_rsv_vs_mock.csv")
coVmock_marks <- read.csv("./cilia_co_vs_mock.csv")
lacVmock_marks <- read.csv("./cilia_lac_vs_mock.csv")
coVrsv_marks <- read.csv("./cilia_co_vs_rsv.csv")

row.names(rsvVlac_marks) <- rsvVlac_marks$X
row.names(rsvVmock_marks) <- rsvVmock_marks$X
row.names(coVmock_marks) <- coVmock_marks$X
row.names(lacVmock_marks) <- lacVmock_marks$X
row.names(coVrsv_marks) <- coVrsv_marks$X

# Significant Genes
sum(rsvVlac_marks$p_val_adj  < 0.05) # 96
sum(coVrsv_marks$p_val_adj  < 0.05) # 101
sum(rsvVmock_marks$p_val_adj  < 0.05) # 142
sum(coVmock_marks$p_val_adj  < 0.05) # 185
sum(lacVmock_marks$p_val_adj  < 0.05) # 174

sum(rsvVlac_marks$p_val_adj  < 0.05 & rsvVlac_marks$avg_log2FC > 0) # 29
sum(rsvVlac_marks$p_val_adj  < 0.05 & rsvVlac_marks$avg_log2FC < 0) # 67

sum(coVrsv_marks$p_val_adj  < 0.05 & coVrsv_marks$avg_log2FC > 0) # 101
sum(coVrsv_marks$p_val_adj  < 0.05 & coVrsv_marks$avg_log2FC < 0) # 0

sum(rsvVmock_marks$p_val_adj  < 0.05 & rsvVmock_marks$avg_log2FC > 0) # 113
sum(rsvVmock_marks$p_val_adj  < 0.05 & rsvVmock_marks$avg_log2FC < 0) # 29

sum(coVmock_marks$p_val_adj  < 0.05 & coVmock_marks$avg_log2FC > 0) # 165
sum(coVmock_marks$p_val_adj  < 0.05 & coVmock_marks$avg_log2FC < 0) # 20

sum(lacVmock_marks$p_val_adj  < 0.05 & lacVmock_marks$avg_log2FC > 0) # 161
sum(lacVmock_marks$p_val_adj  < 0.05 & lacVmock_marks$avg_log2FC < 0) # 13


# DESeq
table(Idents(combined.all),combined.all$sample)
library(edgeR)
hsc = subset(combined.all,idents='HSC')
Idents(hsc) = hsc$sample
DefaultAssay(hsc) = 'RNA'
hsc = NormalizeData(hsc) %>% ScaleData()
mtx = as.matrix(GetAssayData(hsc,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=hsc$sample)
hsc$sample = factor(hsc$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + hsc$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

write.csv(pval, "./hsc_pval.csv")

