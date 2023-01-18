# scRNA Analysis Blood

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Blood")
# saveRDS(combined.all, "./combined_all.rds")
combined.all <- readRDS("./combined_all.rds")

# read in Infected sample
Blood_HIV = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Blood/INF/filtered_feature_bc_matrix")
dim(Blood_HIV)
#36601 4646

# Create Sample1 Seurat Object (Infected)
infected_Blood <- CreateSeuratObject(counts = Blood_HIV, project = "Blood_infected", min.cells = 3, min.features = 200)
dim(infected_Blood)
#17692 4642
head(infected_Blood)

# read sample 2
Blood_NI = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Blood/NI/filtered_feature_bc_matrix")
dim(Blood_NI)
#36601 7117


# create sample2 suerat object
NotInf_Blood <- CreateSeuratObject(counts = Blood_NI, project = "NotInf_Blood", min.cells = 3, min.features = 200)
dim(NotInf_Blood)
# 18918 7106
head(NotInf_Blood)

######################################################################################
# Identify Condition
infected_Blood[["sample"]] <- "Infected"
NotInf_Blood[["sample"]] <- "Not"

# Merge Samples
merged_seurat <- merge(x=infected_Blood, y=NotInf_Blood, add.cell.id = c("Infected", "Not"))
View(merged_seurat@meta.data)

sum(grepl("MT",rownames(merged_seurat)))
# 226
rownames(merged_seurat)[grepl("MT-",rownames(merged_seurat))]

# calculate the proportion of transcripts mapping to mitochondrial genes
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "MT-")

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

#save(merged_seurat, file="./merged_Blood_seurat.RData")
#load("./merged_Blood_seurat.RData")
max(merged_seurat@meta.data$percent.mt)
# 58.39616

# PLot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Filter based on quality plots 
Filter_Blood <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt > -Inf & percent.mt <10)
VlnPlot(Filter_Blood, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filter_Blood, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Sums post Filtering
sum(Filter_Blood$sample=="Infected")
# 4361
sum(Filter_Blood$sample=="Not")
# 6708

Filter_Blood@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38filtered_feature_bc_matrix-","",Filter_Blood@assays$RNA@counts@Dimnames[[1]])
Filter_Blood@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filter_Blood@assays$RNA@counts@Dimnames[[1]])

######### Visualize Mitochondria
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
############

split_seurat <- SplitObject(Filter_Blood, split.by = "sample")
split_seurat <- split_seurat[c("Infected", "Not")]
#options(future.globals.maxSize = 4000 * 1024^2)

############################################################## DOUBLET FINDER
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


################### Find Doublets in Infected Blood
# (pk identification) 
sweep.res.list_1 <- paramSweep_v3(split_seurat[[1]], PCs= 1:10, sct=T)
sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
bcmvn_1 <- find.pK(sweep.stats_1)
#pK = 0.3

###### Add Labels to pK plot
pK=as.numeric(as.character(bcmvn_1$pK))
BCmetric=bcmvn_1$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
######

load("./firstsweep_Blood.RData")
# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[1]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0338*nrow(split_seurat[[1]]@meta.data)) # Adjust according to Doublet Rate  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.3, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_147", sct=T)


# PLot Doublets & Singlets
DimPlot(split_seurat[[1]], reduction = "umap", group.by = c("DF.classifications_0.25_0.3_147",
                                                            "DF.classifications_0.25_0.3_128"),label = TRUE,repel = TRUE)
sum(split_seurat[[1]]$DF.classifications_0.25_0.3_147=="Doublet") #147

################### Doublet Finder for Non-Infected Blood
# (pk identification) 
sweep.res.list_2 <- paramSweep_v3(split_seurat[[2]], PCs= 1:10, sct=T)
sweep.stats_2 <- summarizeSweep(sweep.res.list_2, GT=F)
bcmvn_2 <- find.pK(sweep.stats_2)
# pK = 0.07
######## label Pk chart
pK=as.numeric(as.character(bcmvn_2$pK))
BCmetric=bcmvn_2$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
#######scRNA_sara/scRNA_Blood_sara/Blood_NI_mixREF

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[2]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0516*nrow(split_seurat[[2]]@meta.data)) # Adjust according to Doublet Rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder ---- edit pK
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.07, nExp = nExp_poi.adj, reuse.pANN ="pANN_0.25_0.07_346", sct=T)

# Plot Doublets vs. Singlets
DimPlot(split_seurat[[2]], reduction = "umap", group.by = c("DF.classifications_0.25_0.07_346",
                                                            "DF.classifications_0.25_0.07_299"),label = TRUE,repel = TRUE)

# REMOVE DOUBLETS #
split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.3_147!="Doublet") 
dim(split_seurat[[1]])
#4214
split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset =DF.classifications_0.25_0.07_346 !="Doublet") 
dim(split_seurat[[2]])
#6362

########################################################## Downsample (4200 cells/sample)
set.seed(111)
split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 4200, replace=F)]
split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 4200, replace=F)]

for (i in 1:length(split_seurat)){
  split_seurat[[i]] <- RunTSNE(split_seurat[[i]], dims = 1:20)
}


########################################################### Integration
# select most variable features
integ_features <- SelectIntegrationFeatures(object.list=split_seurat, nfeatures=10000)

# prepare SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list=split_seurat, anchor.features=integ_features)

# Find Anchors
integ_anchors <- FindIntegrationAnchors(object.list=split_seurat, normalization.method = "SCT",
                                        anchor.features=integ_features)

# Integrate
combined.all <- IntegrateData(anchorset=integ_anchors, normalization.method="SCT")

# Save integrated seurat object & workspace
# saveRDS(combined.all, "./combined_all.rds") # 10/19
combined.all <- readRDS("./combined_all.rds")

combined.all@meta.data

DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) # reduction = "pca",
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = 0.7)
combined.all <- RunTSNE(object = combined.all)
combined.all <- FindVariableFeatures(combined.all, selection.method = "vst", nfeatures = 2000)
combined.all[["percent.mt"]] <- PercentageFeatureSet(combined.all, pattern = "MT-",assay = "RNA")


VlnPlot(object = combined.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(combined.all@meta.data, 5)
rownames(combined.all)


# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(combined.all, assay='RNA', slot='counts')
writeMM(counts_matrix, file='/storage1/fs1/leyao.wang/Active/scRNA_Blood_sara/counts_Blood_organoid.mtx')
rownames(counts_matrix)


######## rename genes so Clusters can be recognized
### combined.all
combined.all@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$RNA@counts@Dimnames[[1]])
combined.all@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$RNA@counts@Dimnames[[1]])

combined.all@assays$SCT@counts@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$SCT@counts@Dimnames[[1]])
combined.all@assays$SCT@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$SCT@data@Dimnames[[1]])
rownames(combined.all@assays$SCT@scale.data) <- gsub("GRCh38-","",rownames(combined.all@assays$SCT@scale.data))

combined.all@assays$integrated@data@Dimnames[[1]] <- gsub("GRCh38-","",combined.all@assays$integrated@data@Dimnames[[1]])
rownames(combined.all@assays$integrated@scale.data) <- gsub("GRCh38-","",rownames(combined.all@assays$integrated@scale.data))

### split_seurat rename
for (i in 1:length(split_seurat)) {
  split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]])
  split_seurat[[i]]@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]])
  
  split_seurat[[i]]@assays$SCT@counts@Dimnames[[1]] <- gsub("GRCh38-","",split_seurat[[i]]@assays$SCT@counts@Dimnames[[1]])
  split_seurat[[i]]@assays$SCT@data@Dimnames[[1]] <- gsub("GRCh38-","",split_seurat[[i]]@assays$SCT@data@Dimnames[[1]])
  rownames(split_seurat[[i]]@assays$SCT@scale.data) <- gsub("GRCh38-","",rownames(split_seurat[[i]]@assays$SCT@scale.data))
  
  split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]] <- gsub("mm10---","",split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]])
  split_seurat[[i]]@assays$RNA@data@Dimnames[[1]] <- gsub("mm10---","",split_seurat[[i]]@assays$RNA@counts@Dimnames[[1]])
  
  split_seurat[[i]]@assays$SCT@counts@Dimnames[[1]] <- gsub("mm10---","",split_seurat[[i]]@assays$SCT@counts@Dimnames[[1]])
  split_seurat[[i]]@assays$SCT@data@Dimnames[[1]] <- gsub("mm10---","",split_seurat[[i]]@assays$SCT@data@Dimnames[[1]])
  rownames(split_seurat[[i]]@assays$SCT@scale.data) <- gsub("mm10---","",rownames(split_seurat[[i]]@assays$SCT@scale.data))
}

# Annotations
# ORIG RESOLUTION = combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.5" )
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.7")
DimPlot(combined.all, label=TRUE)
head(combined.all@meta.data)
DefaultAssay(combined.all) <- "RNA"

# ORIGINAL 0.5 RESOLUTION
NKcells <- WhichCells(combined.all, idents= "4")
Myeloid <- WhichCells(combined.all, idents= c("6"))
Bcells <- WhichCells(combined.all, idents= c("2", "1", "8"))
Tcells <- WhichCells(combined.all, idents= c("3", "0", "5", "7", "10", "9"))
Tcell <- subset(combined.all, cells=Tcells)
Bcell <- subset(combined.all, cells=Bcells)
myeloid <- subset(combined.all, cells= Myeloid)
NKcells <- subset(combined.all, cells=NKcells)


# 0.7 RESOLUTION
NKcells <- WhichCells(combined.all, idents= "4")
Myeloid <- WhichCells(combined.all, idents= c("7"))
Bcells <- WhichCells(combined.all, idents= c("2", "0", "9"))
Tcells <- WhichCells(combined.all, idents= c("3", "10", "1", "6", "8", "12", "5", "11"))
Tcell <- subset(combined.all, cells=Tcells)
Bcell <- subset(combined.all, cells=Bcells)
myeloid <- subset(combined.all, cells= Myeloid)
NKcell <- subset(combined.all, cells=NKcells)
combined.all@meta.data$annotations2 <- "NA"
combined.all@meta.data[NKcells, "annotations2"] <- "NK"
combined.all@meta.data[Tcells, "annotations2"] <- "T Cell"
combined.all@meta.data[Myeloid, "annotations2"] <- "Myeloid"
combined.all@meta.data[Bcells, "annotations2"] <- "B Cell "
combined.all <- SetIdent(combined.all, value = "annotations")
DimPlot(combined.all, label=TRUE)


# T cells
FeaturePlot(combined.all, features = c("CD8A", "CD4", "CD8B", "CD40LG", "CD3D", "CD3E"), cols= c("yellow", "blue"), label = TRUE)
# nk cells
FeaturePlot(combined.all, features = c("GNLY", "NKG7"), cols= c("yellow", "blue"), label = TRUE)
combined.all@meta.data[NKcells, "annotations"] <- "NK"
# b cell
FeaturePlot(combined.all, features = c("CD24", "CD79B", "CD19", "CD27", "CD38", "MS4A1", "VPREB3", "RAG1"), cols= c("yellow", "blue"), label = TRUE)
# myeloid
FeaturePlot(combined.all, features = c( "FCGR3A", "CD14", "FCN1",  "CX3CR1", "S100A12"), cols= c("yellow", "blue"), label = TRUE)

#Tcells
# if using annotations2
Tcells <- WhichCells(combined.all, idents= c("3", "10", "1", "6", "8", "12", "5", "11"))
Tcell <- subset(combined.all, cells=Tcells)

DefaultAssay(Tcell) <- "integrated"
DefaultAssay(Tcell) <- "RNA"
Tcell <- RunPCA(Tcell, npcs = 30, verbose = FALSE)
Tcell <- RunUMAP(Tcell, dims = 1:20) 
Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:20)
Tcell <- FindClusters(Tcell, resolution = 0.6)
Tcell <- SetIdent(Tcell, value = "integrated_snn_res.0.6")
DimPlot(Tcell, label=TRUE)

DoHeatmap(Tcell, features = c('CD40LG', "CD4","CD8A", "CD8B", "CCR7", "IL7R", "SELL", "DTHD1", "LYST", "LEF1", "CCR4", "GZMH", "GZMK",
                              "GNLY", "NKG7", 'TNIP3', "CCL5", "FOXP3", "TNFRSF18", "CTLA4", "CXCR6", "CXCR3", 
                              "CCL4", "CCL3", "CMC1", "KLRC2", "FGFBP2", "PDCD1", "TIGIT", "HAVCR2", "LAG3", 
                              "PRDM1", "CD28", "IL2RB", "IFNG", "CD27", "GZMA", "GZMB", "IFNG", "MKI67", "TOP2A", "KLRB1",
                              "USP10", "LGALS1", 'TIMP1', "GZMA", "GZMK", "CD27", "LYAR", "CST7", "IL2RA"), assay = "RNA", slot = "data", angle = 90) + NoLegend()

FeaturePlot(Tcell, features = c("CCR7", "IL7R", "SELL", "DTHD1", "LYST", "LEF1", "GZMH", "GZMK", "KLRG1"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("GNLY", "NKG7", 'CCR7', "CCL5", "FOXP3", "LEF1", "CTLA4", "CXCR6", "CXCR3", "GZMH", "GZMK"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("CD14", "ITGAM", "CCR2", "S100A9", "S100A8"), cols= c("yellow", "blue"), label = TRUE)

cd4t <- WhichCells(Tcell, idents = c("2", "4")) # 293 Cells
cd8t <- WhichCells(Tcell, idents = c("0", "1", "6", "3", "5", "7", "8")) # 1198 Cells
cd4 <- subset(Tcell, cells = cd4t)
cd8 <- subset(Tcell, cells = cd8t)

# CD4
DefaultAssay(cd4) <- "integrated"
cd4 <- RunPCA(cd4, npcs = 30, verbose = FALSE)
cd4 <- RunUMAP(cd4, dims = 1:20) # reduction = "pca",
cd4 <- FindNeighbors(cd4, reduction = "pca", dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 0.5)
cd4<- SetIdent(cd4, value = "integrated_snn_res.1")
DimPlot(cd4, label=T, label.size = 6)
DoHeatmap(cd4, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B", "CCR5", "S100A4", 
                            "IL32", "ISG15", "KLRB1", "IL7R"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd4_mas <- WhichCells(cd4, idents = c("0"))
cd4MAS <- subset(cd4, cells= cd4_mas)
DefaultAssay(cd4MAS) <- "integrated"
cd4MAS <- RunPCA(cd4MAS, npcs = 30, verbose = FALSE)
cd4MAS <- RunUMAP(cd4MAS, dims = 1:20) # reduction = "pca",
cd4MAS <- FindNeighbors(cd4MAS, reduction = "pca", dims = 1:20)
cd4MAS <- FindClusters(cd4MAS, resolution = 0.5)
cd4MAS<- SetIdent(cd4MAS, value = "integrated_snn_res.0.5")
DimPlot(cd4MAS, label=T, label.size = 6)
DoHeatmap(cd4MAS, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B", "CCR5", "S100A4", 
                            "IL32", "ISG15", "KLRB1", "IL7R"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd4naivecm <- WhichCells(cd4, idents = c("0", "1", "4", "6"))
tm_treg_temra <- WhichCells(cd4, idents = c("2", "3"))
ilc <- WhichCells(cd4, idents = c("5"))
combined.all@meta.data[cd4naivecm, "annotations"] <- "CD4+ Naive/CM"
combined.all@meta.data[tm_treg_temra, "annotations"] <- "CD4+"
combined.all@meta.data[ilc, "annotations"] <- "ILC"

cd4$celltype.cd4 <- paste(Idents(cd4), cd4$sample, sep = "_")
Idents(cd4) <- "celltype.cd4"
table(cd4@meta.data$celltype.cd4)
DimPlot(cd4, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

# CD8
DefaultAssay(cd8) <- "integrated"
cd8 <- RunPCA(cd8, npcs = 30, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:20) # reduction = "pca",
cd8 <- FindNeighbors(cd8, reduction = "pca", dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 0.6)
cd8<- SetIdent(cd8, value = "integrated_snn_res.0.6")
DimPlot(cd8, label=T, label.size = 6)

DoHeatmap(cd8, features = c("CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                            "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY",
                            "FGFBP2", "NKG7", "CCL4", "KLRC2", "TRDC", "ZNF683", "TRDV1", "CD3D", "CD3E", "MKI67"), assay = "RNA", slot="data", angle = 90) + NoLegend()

cd8$celltype.cd8 <- paste(Idents(cd8), cd8$sample, sep = "_")
Idents(cd8) <- "celltype.cd8"
table(cd8@meta.data$celltype.cd8)
DimPlot(cd8, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

naivecd8 <- WhichCells(cd8, idents = c("0", "1", "3", "4", "7"))
combined.all@meta.data[naivecd8, "annotations"] <- "CD8+ Naive/CM"
ccl5 <- WhichCells(cd8, idents = c("2", "5", '8'))
combined.all@meta.data[ccl5, "annotations"] <- "CD8+ CCL5+"
ccl5cyc <- WhichCells(cd8, idents = c("6"))
combined.all@meta.data[ccl5cyc, "annotations"] <- "CD8+ CCL5+ cycle"


cd8naive <- WhichCells(Tcell, idents = c("0", "1", "5"))
cd4treg <- WhichCells(Tcell, idents = c("4"))
cd8em <- WhichCells(Tcell, idents = c("3", "7", "8"))
cd4naive <- WhichCells(Tcell, idents = c("2"))
c6 <- WhichCells(Tcell, idents = c("6"))
#c8 <- WhichCells(Tcell, idents = c("8"))
c9 <- WhichCells(Tcell, idents = c("9"))
DimPlot(combined.all, cells.highlight = c9)
# combined.all$annotations <- "NA"
combined.all@meta.data[cd8naive, "annotations"] <- "CD8+ Naive"
combined.all@meta.data[cd8em, "annotations"] <- "CD8+ Eff/Mem"
combined.all@meta.data[cd4naive, "annotations"] <- "CD4+ Naive"
combined.all@meta.data[cd4treg, "annotations"] <- "CD4+ Treg"
combined.all@meta.data[c8, "annotations"] <- "C8"
combined.all@meta.data[c9, "annotations"] <- "C9"
combined.all@meta.data[c6, "annotations"] <- "CD8+ EM/EMRA"

# DENDRITIC?
FeaturePlot(combined.all, features = c("IRF8", "HLA-DRA"), cols= c("yellow", "blue"), label = TRUE)

# B cells
NormalizeData(Bcell)
ScaleData(Bcell)
DefaultAssay(Bcell) <- "integrated"
DefaultAssay(Bcell) <- "RNA"
Bcell <- RunPCA(Bcell, npcs = 30, verbose = FALSE)
Bcell <- RunUMAP(Bcell, dims = 1:20) 
Bcell <- FindNeighbors(Bcell, reduction = "pca", dims = 1:20)
Bcell <- FindClusters(Bcell, resolution = 0.5)
Bcell <- SetIdent(Bcell, value = "integrated_snn_res.0.5")
DimPlot(Bcell, label=TRUE)
FeaturePlot(Bcell, features = c("CD24", "CD79B", "CD19", "CD27", "CD38", "MS4A1"), cols= c("yellow", "blue"), label = TRUE)
DoHeatmap(Bcell, features = c("CD24", "CD79B", "CD19", "CD27", "CD38", "MS4A1", 'RAG1', "VPREB3", "SELL"), assay = "RNA", slot="data", angle = 90) + NoLegend()

FeaturePlot(Bcell, features = c("MS4A1"), cols= c("yellow", "blue"), label = TRUE)
combined.all@meta.data[Bcells, "annotations"] <- "B Cell"
bcell_markers <- FindAllMarkers(Bcell)

# NK CELLS
DefaultAssay(NKcells) <- "integrated"
DefaultAssay(NKcell) <- "RNA"
NKcells <- RunPCA(NKcells, npcs = 30, verbose = FALSE)
NKcells <- RunUMAP(NKcells, dims = 1:20) 
NKcells <- FindNeighbors(NKcells, reduction = "pca", dims = 1:20)
NKcells <- FindClusters(NKcells, resolution = 0.5)
NKcells <- SetIdent(NKcells, value = "integrated_snn_res.0.5")
DimPlot(NKcells, label=TRUE)
FeaturePlot(NKcell, features = c("CCL3", "FGFBP2", "KIR2DL4"), cols= c("yellow", "blue"), label = TRUE)
DoHeatmap(NKcell, features = c("NCAM1", "GNLY", "NKG7", "EOMES", "GZMB",
                               "CD3D", "CD3E", "CCL3", "FGFBP2", "KIR2DL4", "IL7R", "XCL1", "KLRB1", "KLRF1"), assay = "RNA", slot="data", angle = 90) + NoLegend()
ccl3 <- WhichCells(NKcell, idents = c("1"))
combined.all@meta.data[ccl3, "annotations"] <- "NK_CCL3+"
fgfbp2 <- WhichCells(NKcell, idents = c("0"))
combined.all@meta.data[fgfbp2, "annotations"] <- "NK_FGFBP2+"
kir2dl4 <- WhichCells(NKcell, idents = c("2"))
combined.all@meta.data[kir2dl4, "annotations"] <- "NK_KIR2DL4+"

NKcells <- SetIdent(NKcells, value = "annotations")
NKcells$celltype <- paste(Idents(NKcells), NKcells$sample, sep = "_")
Idents(NKcells) <- "celltype"
table(NKcells@meta.data$celltype)
DimPlot(NKcells, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

#myeloid
DefaultAssay(myeloid) <- "integrated"
DefaultAssay(myeloid) <- "RNA"
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:20) 
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.7)
myeloid <- SetIdent(myeloid, value = "annotations")
DimPlot(myeloid, label=TRUE)
FeaturePlot(myeloid, features = c( "FCGR3A", "CD14", "FCN1",  "CX3CR1", "S100A12", "FCGR3B", "GATA2", "SELL"), cols= c("yellow", "blue"), label = TRUE)
cd14 <- WhichCells(myeloid, idents = c("0", "2"))
cd16 <- WhichCells(myeloid, idents = c("1,3"))
combined.all@meta.data[cd14, "annotations"] <- "CD14 Mono"
combined.all@meta.data[cd16, "annotations"] <- "CD16 Mono"

DoHeatmap(myeloid, features = c( "FCGR3A", "CD14", "FCN1",  "CX3CR1", "S100A12", "FCGR3B", "GATA2", "SELL"), assay = "RNA", slot="data", angle = 90) + NoLegend()

myeloid$celltype.mye <- paste(Idents(myeloid), myeloid$sample, sep = "_")
Idents(myeloid) <- "celltype.mye"
table(myeloid@meta.data$celltype.mye)
DimPlot(myeloid, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

#NA annotations
na <- WhichCells(combined.all, idents= "NA")
combined.all@meta.data[na, "annotations"] <- "B Cell"
Na <- subset(combined.all, cells= na)
Na <- RunPCA(Na, npcs = 30, verbose = FALSE)
Na <- RunUMAP(Na, dims = 1:20) 
Na <- FindNeighbors(Na, reduction = "pca", dims = 1:20)
Na <- FindClusters(Na, resolution = 0.6)
Na <- SetIdent(Na, value = "integrated_snn_res.0.6")
DimPlot(Na, label=T)

# Dim PLot
combined.all <- SetIdent(combined.all, value = "annotations" )
DimPlot(combined.all, label=TRUE)

# Plot UMAP split by infection status
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

# add split annotations
combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"

# table cell counts by annotation + infection
table(combined.all@meta.data$celltype.stim)
table(combined.all@meta.data$annotations)

# Find DE genes between groups by infection 
cd8_EM_DEgenes <- FindMarkers(combined.all, ident.1 = "CD8+ Eff/Mem_Infected", ident.2 = "CD8+ Eff/Mem_Not", verbose = FALSE)
head(cd8_EM_DEgenes, n = 15)

cd4_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ _Infected", ident.2 = "CD4+ _Not", verbose = FALSE)
saveRDS(cd4_DEgenes, "./cd4DEgenes.rds")
cd16_DEgenes <- FindMarkers(combined.all, ident.1 = "CD16 Mono_Infected", ident.2 = "CD16 Mono_Not", verbose = FALSE)
saveRDS(cd16_DEgenes, "./cd16DEgenes.rds")
cd4_naive_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ Naive_Infected", ident.2 = "CD4+ Naive_Not", verbose = FALSE)
saveRDS(cd4_naive_DEgenes, "./cd4naiveDEgenes.rds")
cd8_undet_DEgenes <- FindMarkers(combined.all, ident.1 = "CD8+ ?_Infected", ident.2 = "CD8+ ?_Not", verbose = FALSE)
saveRDS(cd8_undet_DEgenes, "./cd8unDEgenes.rds")
cd14_DEgenes <- FindMarkers(combined.all, ident.1 = "CD14 Mono_Infected", ident.2 = "CD14 Mono_Not", verbose = FALSE)
cd4treg_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ Treg_Infected", ident.2 = "CD4+ Treg_Not", verbose = FALSE)
bcell_DEgenes <- FindMarkers(combined.all, ident.1 = "B Cell_Infected", ident.2 = "B Cell_Not", verbose = FALSE)

# Volcano PLot
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
cd8EM <- subset(combined.all, idents = "CD8+ Eff/Mem")
Idents(cd8EM) <- "sample"
avg.cd8EM <- as.data.frame(log10(AverageExpression(cd8EM, verbose = FALSE)$RNA))
avg.cd8EM$gene <- rownames(avg.cd8EM)
genes.to.label = c("GNLY", "IFI27", "B2M", "GZMB", "CCL4", "NKG7", "TMSB4X")
p1 <- ggplot(avg.cd8EM, aes(Infected, Not)) + geom_point() + ggtitle("CD8+ Effector/Memory")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd4N <- subset(combined.all, idents = "CD4+ Naive")
Idents(cd4N) <- "sample"
avg.cd4N <- as.data.frame(log10(AverageExpression(cd4N, verbose = FALSE)$RNA))
avg.cd4N$gene <- rownames(avg.cd4N)
genes.to.label = c("ISG15", "CRIP1", "IFI27", "S100A6", "S100A10", "MALAT1")
p2 <- ggplot(avg.cd4N, aes(Infected, Not)) + geom_point() + ggtitle("CD4+ Naive")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd16M <- subset(combined.all, idents = "CD16 Mono")
Idents(cd16M) <- "sample"
avg.cd16M <- as.data.frame(log10(AverageExpression(cd16M, verbose = FALSE)$RNA))
avg.cd16M$gene <- rownames(avg.cd16M)
genes.to.label = c("IGLC3", "CRIP1", "IFI27", "FOS", "HLA-B")
p3 <- ggplot(avg.cd16M, aes(Infected, Not)) + geom_point() + ggtitle("CD16 MONOCYTE")
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd14M <- subset(combined.all, idents = "CD14 Mono")
Idents(cd14M) <- "sample"
avg.cd14M <- as.data.frame(log10(AverageExpression(cd14M, verbose = FALSE)$RNA))
avg.cd14M$gene <- rownames(avg.cd14M)
genes.to.label = c("GNLY", "NEAT1", "IFI27", "CCL5", "NKG7", "MT2A")
p4 <- ggplot(avg.cd14M, aes(Infected, Not)) + geom_point() + ggtitle("CD14 MONOCYTE")
p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd4treg <- subset(combined.all, idents = "CD4+ Treg")
Idents(cd4treg) <- "sample"
avg.cd4treg <- as.data.frame(log10(AverageExpression(cd4treg, verbose = FALSE)$RNA))
avg.cd4treg$gene <- rownames(avg.cd4treg)
genes.to.label = c("CXCL10", "NALT1", "AURKC", "MIF-AS1")
p5 <- ggplot(avg.cd4treg, aes(Infected, Not)) + geom_point() + ggtitle("CD4 TREG")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

bcell <- subset(combined.all, idents = "B Cell")
Idents(bcell) <- "sample"
avg.bcell <- as.data.frame(log10(AverageExpression(bcell, verbose = FALSE)$RNA))
avg.bcell$gene <- rownames(avg.bcell)
genes.to.label = c("CCDC136", "TSHR", "CPAMD8")
p4 <- ggplot(avg.bcell, aes(Infected, Not)) + geom_point() + ggtitle("B Cell")
p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)



#CD11C
FeaturePlot(combined.all, label = T, cols=c("yellow","blue"), features = c("ITGAX"))
AverageExpression(combined.all, features = "ITGAX", slot="data", assays = "RNA")
DoHeatmap(combined.all, features = c("ITGAX"), assay = "RNA", slot="data", angle = 90) + NoLegend()

# HONGBO
FeaturePlot(Tcell, cols= c("grey", "red"), label=T, features = c("CCR5", "TCF7","CX3CR1", "BACH2"))
DoHeatmap(Tcell, features = c("CCR5", "TCF7", "CX3CR1", "BACH2"), assay = "RNA", slot="data", angle = 90) + NoLegend()

FeaturePlot(NKcells, cols=c("grey", "red"), label=T, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                  "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"))
DoHeatmap(NKcells, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"), assay = "RNA", slot="data", angle = 90) + NoLegend()


cd4 <- subset(combined.all, idents = c("CD4+ Naive/CM", "CD4+"))
cd4 <- ScaleData(cd4, assay = "RNA")
cd4 <- SetIdent(cd4, value = "sample")
FeaturePlot(cd4, cols= c("grey", "red"), label=T, features = c("CX3CR1", "BACH2"))
DoHeatmap(cd4, features = c("CCR5", "TCF7", "CX3CR1", "BACH2"), assay = "RNA", slot="scale.data", angle = 90)



CCR5 = GetAssayData(object = cd4, 
                    assay = "RNA", slot = "data")["CCR5",]
TCF7 = GetAssayData(object = cd4, 
                    assay = "RNA", slot = "data")["TCF7",]
CD4 = GetAssayData(object = cd4, 
                   assay = "RNA", slot = "data")["CD4",]

names(which(TCF7>0 & CCR5<=0)) # 760
ccr5 <- names(which(CCR5>0 & TCF7<=0)) #9
names(which(TCF7<=0 & CCR5<=0)) # 208
dp <- names(which(CCR5>0 & TCF7>0)) # 6
DimPlot(cd4, cells.highlight = dp, label=T)

CCR5 <- subset(combined.all, cells = c(ccr5, dp))
CCR5@meta.data$group <- "NA"
CCR5@meta.data[ccr5, "group"] <- "CCR5+ TCF7-"
CCR5@meta.data[dp, "group"] <- "CCR5+ TCF7+"
CCR5 <- NormalizeData(CCR5, assay = "RNA")
CCR5 <- ScaleData(CCR5, assay = "RNA")

DefaultAssay(CCR5) <- "integrated"
CCR5 <- RunPCA(CCR5, npcs = 30, verbose = FALSE)
CCR5 <- RunUMAP(CCR5, dims = 1:20) 
CCR5 <- FindNeighbors(CCR5, reduction = "pca", dims = 1:20)
CCR5 <- FindClusters(CCR5, resolution = 0.3)
CCR5 <- SetIdent(CCR5, value = "group")
DimPlot(CCR5)
ccr5_df <- as.data.frame(AverageExpression(CCR5, verbose = FALSE, slot = "scale.data", features = c("BACH2", "CX3CR1", "PDCD1", "TOX", "CXCL13", "TIGIT", 
                                                                                                    "CTLA4", "TNFRSF9", "HAVCR2", "LAG3", "CST7", "GZMK", 
                                                                                                    "GZMA", "NKG7", "IFNG", "PRF1", "GZMB", "GNLY", "PRDM1", 
                                                                                                    "BATF", "CD28", "GZMH", "KLRK1", "KLRB1", "CTSW", "BTLA", 
                                                                                                    "FOXP3", "TNFRSF4", "EEF1A1", "SELL", "CCR7", "IL6R", 
                                                                                                    "IGFBP4", "IGFL2I", "IL7R") )$RNA)
write.csv(ccr5_df, "/storage1/fs1/liang.shan/Active/scRNA_sara/CCR5_expr_blood.csv")

DoHeatmap(CCR5, features = c("BACH2", "CX3CR1", "PDCD1", "TOX", "CXCL13", "TIGIT", 
                             "CTLA4", "TNFRSF9", "HAVCR2", "LAG3", "CST7", "GZMK", 
                             "GZMA", "NKG7", "IFNG", "PRF1", "GZMB", "GNLY", "PRDM1", 
                             "BATF", "CD28", "GZMH", "KLRK1", "KLRB1", "CTSW", "BTLA", 
                             "FOXP3", "TNFRSF4", "EEF1A1", "SELL", "CCR7", "IL6R", 
                             "IGFBP4", "IGFL2I", "IL7R"), assay = "RNA", slot="data", angle = 90) 

# NK Cells
nkcells <- WhichCells(combined.all, idents = c("NK_CCL3+", "NK_FGFBP2+", "NK_KIR2DL4+"))
NKcells <- subset(combined.all, cells = nkcells)
NKcells<- ScaleData(NKcells, assay = "RNA")
Idents(NKcells) <- "sample"
table(NKcells$sample)
RidgePlot(NKcells, features = c("KIR2DL1",  "KIR2DL3"))

FeaturePlot(NKcells, cols=c("grey", "red"), label=T, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                  "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"))
DoHeatmap(NKcells, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"), assay = "RNA", slot="scale.data", angle = 90) + NoLegend()

cd158_df <- as.data.frame(AverageExpression(NKcells, verbose = FALSE, slot = "scale.data", features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                                                    "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3") )$RNA)



# CD166
combined.all <- NormalizeData(combined.all, assay = "RNA")
combined.all <- ScaleData(combined.all, assay = "RNA")

cd166 <- as.data.frame(AverageExpression(combined.all, features = c("ALCAM"), verbose = FALSE, slot = "data", assays = "RNA"))
write.csv(cd166, "./CD166_expression.csv")

combined.all <- SetIdent(combined.all, value = "sample")
not <- WhichCells(combined.all, idents = "Not")
combined.all@meta.data[not, "sample"] <- "Mock"

Idents(combined.all) <- "celltype.stim"
interest <- WhichCells(combined.all, idents = c("CD14 Mono_Infected", "CD14 Mono_Mock", "CD16 Mono_Infected", "CD16 Mono_Mock", "CD4+_Infected", "CD4+_Mock"))
INTEREST <- subset(combined.all, cells = interest)
inf14 <- WhichCells(combined.all, idents = c("CD14 Mono_Infected"))
mock14 <- WhichCells(combined.all, idents = c("CD14 Mono_Mock"))
inf16 <- WhichCells(combined.all, idents = c("CD16 Mono_Infected"))
mock16 <- WhichCells(combined.all, idents = c("CD16 Mono_Mock"))

FetchData(combined.all, cells = inf16,  vars =c("ALCAM"), slot="data")
inf14<- FetchData(combined.all, cells = inf14,  vars =c("ALCAM"), slot="counts")
mock16 <- as.data.frame(FetchData(combined.all, cells = mock16,  vars =c("ALCAM"), slot="data"))
mock14 <- as.data.frame(FetchData(combined.all, cells = mock14,  vars =c("ALCAM"), slot="counts"))

FeaturePlot(combined.all, features = "ALCAM", cols =c("grey", "red"), split.by = "sample", label = TRUE, label.size = 3)
my_levels <- c("CD14 Mono_Infected", "CD14 Mono_Mock", "CD16 Mono_Infected", "CD16 Mono_Mock", "CD4+_Infected", "CD4+_Mock")
# Re-level object@ident
INTEREST@active.ident <- factor(x = INTEREST@active.ident, levels = my_levels)

RidgePlot(INTEREST, features = "ALCAM", slot = "data")

my_levels <- c("B Cell_Infected", "B Cell_Mock", "C9_Infected", "C9_Mock", "CD14 Mono_Infected", "CD14 Mono_Mock", "CD16 Mono_Infected", "CD16 Mono_Mock",
               "CD4+ Naive/CM_Infected", "CD4+ Naive/CM_Mock", "CD4+_Infected", "CD4+_Mock", "CD8+ CCL5+ cycle_Infected", "CD8+ CCL5+ cycle_Mock", "CD8+ CCL5+_Infected", "CD8+ CCL5+_Mock",
               "CD8+ Naive/CM_Infected", "CD8+ Naive/CM_Mock", "ILC_Infected", "ILC_Mock", "NK_CCL3+_Infected", "NK_CCL3+_Mock", "NK_FGFBP2+_Infected", "NK_FGFBP2+_Mock", "NK_KIR2DL4+_Infected", "NK_KIR2DL4+_Mock")
combined.all@active.ident <- factor(x = combined.all@active.ident, levels = my_levels)
FeaturePlot(combined.all, features = "ALCAM", cols = c("grey", "red"), label = TRUE)

AverageExpression(INTEREST, verbose = FALSE, slot = "counts", assays = "RNA", features = "ALCAM")

Idents(INTEREST) <- "annotations"
VlnPlot(INTEREST, features = c("ALCAM"), split.by = "sample", 
        pt.size = 0, combine = FALSE, assay = "RNA", slot = "counts")

combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all$celltype.stim)

combined.all <- SetIdent(combined.all, value = "annotations" )
DimPlot(combined.all, label=TRUE)




########################################################### RITU NK CELL ANALYSIS
combined.all <- SetIdent(combined.all, value = "annotations")
DimPlot(combined.all, label = T)

FeaturePlot(combined.all, features = c("GNLY", "PRF1", "TYROBP", "TRDC"), label =T)

NK <- subset(combined.all, idents = c("NK_CCL3+", "NK_KIR2DL4+", "NK_FGFBP2+"))
FeaturePlot(NK, features = c("GNLY", "PRF1", "TYROBP", "TRDC"), label = T)

DefaultAssay(NK) <- "integrated"
NK <- RunPCA(NK, npcs = 30, verbose = FALSE)
NK <- RunUMAP(NK, dims = 1:20) 
NK <- FindNeighbors(NK, reduction = "pca", dims = 1:20)
NK <- FindClusters(NK, resolution = 0.6)
NK <- SetIdent(NK, value = "integrated_snn_res.0.6")
DimPlot(NK, label=TRUE, split.by = "sample", label.size = 6)


NK <- NormalizeData(NK, assay = "RNA")
NK <- ScaleData(NK, assay = "RNA")

NK <- SetIdent(NK, value = "sample")
DoHeatmap(NK, features = c("CD52", "CD34", "NCAM1", "FCGR3A", "ITGAM", "B3GAT1", "KLRD1", "KLRC2", "IL2RB", "CD27"), assay = "RNA", slot="data", angle = 90) + NoLegend()
nklung_df_stim <- as.data.frame(AverageExpression(NK, verbose = FALSE, slot = "counts", features = c("CD52", "CD34", "NCAM1", "FCGR3A", "ITGAM", "B3GAT1", "KLRD1", "KLRC2", "IL2RB", "CD27") )$RNA)

NK$nk.infection <- paste(Idents(NK), NK$sample, sep = "_")
Idents(NK) <- "nk.infection"
table(NK@meta.data$nk.infection)

my_levels <- c("0_Infected", "0_Mock",  "1_Infected", "1_Mock", "2_Infected", "2_Mock")
# Re-level object@ident
NK@active.ident <- factor(x = NK@active.ident, levels = my_levels)

genes <- read.csv("/storage1/fs1/liang.shan/Active/RITU/HeatMapGenes.csv")
DoHeatmap(NK, features = genes[,1], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[1:50,2], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[51:99,2], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[1:55,3], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[56:111,3], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[,4], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[,5], assay = "RNA", slot="scale.data", angle = 90) 


write.csv(nklung_df, "./NKlung.csv")

NK <- SetIdent(NK, value = "integrated_snn_res.0.6")
NK <- SetIdent(NK, value = "sample")
NKmarks <- FindAllMarkers(NK, assay = "RNA", slot = "data")

nk1 <- WhichCells(NK, idents = c("1"))
NK@meta.data[nk1, "more"] <- "NK1"
nk2 <- WhichCells(NK, idents = c("2"))
NK@meta.data[nk2, "more"] <- "NK2"
nk0 <- WhichCells(NK, idents = c("0"))
NK@meta.data[nk0, "more"] <- "NK0"
NK <- SetIdent(NK, value = "more")


NK$nk.infection <- paste(Idents(NK), NK$sample, sep = "_")
Idents(NK) <- "nk.infection"
table(NK@meta.data$nk.infection)

my_levels <- c("0_Infected", "0_Mock",  "1_Infected", "1_Mock", "2_Infected", "2_Mock")
# Re-level object@ident
NK@active.ident <- factor(x = NK@active.ident, levels = my_levels)

NK0 <- subset(NK, cells = nk0)
NK0 <- ScaleData(NK0, assay = "RNA")
NK0marks <- FindAllMarkers(NK0, assay = "RNA", slot = "data")
write.csv(NK0marks, "/storage1/fs1/liang.shan/Active/RITU/blood0marks.csv")

NK1 <- subset(NK, cells = nk1)
NK1 <- ScaleData(NK1, assay = "RNA")
NK1marks <- FindAllMarkers(NK1, assay = "RNA", slot = "data")
write.csv(NK1marks, "/storage1/fs1/liang.shan/Active/RITU/blood1marks.csv")

NK2 <- subset(NK, cells = nk2)
NK2 <- ScaleData(NK2, assay = "RNA")
NK2marks <- FindAllMarkers(NK2, assay = "RNA", slot = "data")
write.csv(NK2marks, "/storage1/fs1/liang.shan/Active/RITU/blood2marks.csv")


NKcell_blood <- as.data.frame(AverageExpression(NK, verbose = FALSE, slot = "data")$RNA)
# write.csv(NKcell_blood, "/storage1/fs1/liang.shan/Active/NKblood_expression.csv")


bl0 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/blood0marks.csv", header = T)
bl0genes <- bl0$gene
bl0UP <- bl0genes[1:125]
bl0down <- bl0genes[223:278]
up0 <- as.data.frame(AverageExpression(NK0, verbose = FALSE, slot = "data", features = bl0UP)$RNA)
write.csv(up0, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood0UP.csv")
down0 <- as.data.frame(AverageExpression(NK0, verbose = FALSE, slot = "data", features = bl0down)$RNA)
write.csv(down0, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood0Down.csv")

bl1 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/blood1marks.csv", header = T)
bl1genes <- bl1$gene
bl1UP <- bl1genes[1:150]
bl1down <- bl1genes[278:330]
up1 <- as.data.frame(AverageExpression(NK1, verbose = FALSE, slot = "data", features = bl1UP)$RNA)
write.csv(up1, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood1UP.csv")
down1 <- as.data.frame(AverageExpression(NK1, verbose = FALSE, slot = "data", features = bl1down)$RNA)
write.csv(down1, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood1Down.csv")

bl2 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/blood2marks.csv", header = T)
bl2genes <- bl2$gene
bl2UP <- bl2genes[1:200]
bl2down <- bl2genes[215:274]
up2 <- as.data.frame(AverageExpression(NK2, verbose = FALSE, slot = "data", features = bl2UP)$RNA)
write.csv(up2, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood2UP.csv")
down2 <- as.data.frame(AverageExpression(NK2, verbose = FALSE, slot = "data", features = bl2down)$RNA)
write.csv(down2, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/blood2Down.csv")




######################################################## RITU CD8+ T-Cell Analysis
combined.all <- SetIdent(combined.all, value = "annotations")
DimPlot(combined.all, label = T)

CD8 <- subset(combined.all, idents = c("CD8+ Naive/CM", "CD8+ CCL5+", "CD8+ CCL5+ cycle" ))

FeaturePlot(CD8, features = c("CD8A", "CD3G", "CD3E", "CD3D", "CD8B"))

DefaultAssay(CD8) <- "integrated"
CD8 <- RunPCA(CD8, npcs = 30, verbose = FALSE)
CD8 <- RunUMAP(CD8, dims = 1:20) 
CD8 <- FindNeighbors(CD8, reduction = "pca", dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 0.2)
CD8 <- SetIdent(CD8, value = "integrated_snn_res.0.2")
DimPlot(CD8, label=TRUE, split.by = "sample", label.size = 5)

CD8 <- NormalizeData(CD8, assay = "RNA")
CD8 <- ScaleData(CD8, assay = "RNA")

CD8 <- SetIdent(CD8, value = "sample")
DoHeatmap(CD8, features = c("CD3D", "CD3G", "CD3E", "CD8A", "CD8B", "SELL", "CCR7",
                            "TBX21", "IFNG", "CXCR3", "GZMB", "GZMH", "GNLY", "PRF1", "GZMK", "CRTAM", "CD4", "MKI67", "CD40LG"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd8expr <- as.data.frame(AverageExpression(CD8, verbose = FALSE, slot = "counts", features = c("CD3D", "CD3G", "CD3E", "CD8A", "CD8B", "SELL", "CCR7",
                                                                                               "TBX21", "IFNG", "CXCR3", "GZMB, GZMH", "GNLY", "PRF1", "GZMK, CRTAM") )$RNA)

naive <- WhichCells(CD8, idents = c(0,3,5))
ememra <- WhichCells(CD8, idents = c(1,4))
trmem <- WhichCells(CD8, idents = c(2))
CD8@meta.data[naive, "cd8_more"] <- "CD8+ Naive"
CD8@meta.data[ememra, "cd8_more"] <- "CD8+ Tem/emra"
CD8@meta.data[trmem, "cd8_more"] <- "CD8+ Trm/em"
combined.all@meta.data[trmem, "more"] <- "CD8+ Trm/em"
combined.all@meta.data[naive, "more"] <- "CD8+ Naive/CM"
combined.all@meta.data[ememra, "more"] <- "CD8+ Tem/emra"

Idents(CD8) <- "cd8_more"
CD8$cd8.infection <- paste(Idents(CD8), CD8$sample, sep = "_")
Idents(CD8) <- "cd8.infection"
table(CD8@meta.data$cd8.infection)

cd8naive <- subset(CD8, cells = naive)
cd8naive <- ScaleData(cd8naive, assay = "RNA")
cd8naive_marks <- FindAllMarkers(cd8naive, assay = "RNA", slot = "data")
write.csv(cd8naive_marks, "/storage1/fs1/liang.shan/Active/RITU/bloodCD8naive_marks.csv" )

cd8em <- subset(CD8, cells = ememra)
cd8em <- ScaleData(cd8em, assay = "RNA")
cd8em_marks <- FindAllMarkers(cd8em, assay = "RNA", slot = "data")
write.csv(cd8em_marks, "/storage1/fs1/liang.shan/Active/RITU/bloodCD8ememra_marks.csv")

cd8trm <- subset(CD8, cells = trmem)
cd8trm <- ScaleData(cd8trm, assay = "RNA")
cd8trm_marks <- FindAllMarkers(cd8trm, assay = "RNA", slot = "data")
write.csv(cd8trm_marks, "/storage1/fs1/liang.shan/Active/RITU/bloodCD8trmem_marks.csv")


bloodN <- read.csv("/storage1/fs1/liang.shan/Active/RITU/bloodCD8naive_marks.csv", header = T)
bloodNgenes <- bloodN$gene
bloodNUP <- bloodNgenes[1:220]
bloodNdown <- bloodNgenes[221:293]
upN <- as.data.frame(AverageExpression(cd8naive, verbose = FALSE, slot = "data", features = bloodNUP)$RNA)
write.csv(upN, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodNUP.csv")
downN <- as.data.frame(AverageExpression(cd8naive, verbose = FALSE, slot = "data", features = bloodNdown)$RNA)
write.csv(downN, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodNDown.csv")

bloodEM <- read.csv("/storage1/fs1/liang.shan/Active/RITU/bloodCD8ememra_marks.csv", header = T)
bloodemgenes <- bloodEM$gene
bloodEMUP <- bloodemgenes[1:205]
bloodEMdown <- bloodemgenes[206:448]
upEM <- as.data.frame(AverageExpression(cd8em, verbose = FALSE, slot = "data", features = bloodEMUP)$RNA)
write.csv(upEM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodEMUP.csv")
downEM <- as.data.frame(AverageExpression(cd8em, verbose = FALSE, slot = "data", features = bloodEMdown)$RNA)
write.csv(downEM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodEMDown.csv")

bloodRM <- read.csv("/storage1/fs1/liang.shan/Active/RITU/bloodCD8trmem_marks.csv", header = T)
bloodRMgenes <- bloodRM$gene
bloodRMUP <- bloodRMgenes[1:240]
bloodRMdown <- bloodRMgenes[245:354]
upRM <- as.data.frame(AverageExpression(cd8trm, verbose = FALSE, slot = "data", features = bloodRMUP)$RNA)
write.csv(upRM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodRMUP.csv")
downRM <- as.data.frame(AverageExpression(cd8trm, verbose = FALSE, slot = "data", features = bloodRMdown)$RNA)
write.csv(downRM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/bloodRMDown.csv")


























