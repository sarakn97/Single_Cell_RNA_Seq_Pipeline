# scRNA Analysis Bone Marrow (BM)

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
# library(DoubletFinder)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/")
combined.all <- readRDS("./combinedALL.rds")
saveRDS(combined.all, "./combinedALL.rds")

# Read in Cell Calls
BM_Infected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_I/gem_classification.csv")
BM_Infected_mouse <- subset(BM_Infected, BM_Infected$call=="mm10")

BM_NOTinfected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_NI/gem_classification.csv")
BM_NOTinfected_mouse <- subset(BM_NOTinfected, BM_NOTinfected$call=="mm10")

# read in sample
BM_HIV = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_I/filtered_feature_bc_matrix/")
dim(BM_HIV)
#c6888c6  #5356

# remove mouse cells (BM_HIV)  
BM_HIV <- BM_HIV[,!colnames(BM_HIV)%in%BM_Infected_mouse$barcode]
dim(BM_HIV)
#68886   #5318

# Create Sample1 Seurat Object (Infected)
infected_BM <- CreateSeuratObject(counts = BM_HIV, project = "BM_infected", min.cells = 3, min.features = 200)
dim(infected_BM)
#24286 5286
head(infected_BM)

# read sample 2
BM_NI = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_NI/filtered_feature_bc_matrix/")
dim(BM_NI)
#68886  7553

# remove mouse cells (BM_NI) 
BM_NI <- BM_NI[,!colnames(BM_NI)%in%BM_NOTinfected_mouse$barcode]
dim(BM_NI)
#c6888c6 7529

# create sample2 suerat object
NotInf_BM <- CreateSeuratObject(counts = BM_NI, project = "NotInf_BM", min.cells = 3, min.features = 200)
dim(NotInf_BM)
# 29325  7409
head(NotInf_BM)

######################################################################################
# Identify Condition
infected_BM[["sample"]] <- "Infected"
NotInf_BM[["sample"]] <- "Not"

# Merge Samples
merged_seurat <- merge(x=infected_BM, y=NotInf_BM, add.cell.id = c("Infected", "Not"))
# View(merged_seurat@meta.data)

sum(grepl("MT",rownames(merged_seurat)))
# 254
rownames(merged_seurat)[grepl("MT-",rownames(merged_seurat))]

# calculate the proportion of transcripts mapping to mitochondrial genes
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "MT-")

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

#save(merged_seurat, file="./merged_BM_seurat.RData")
#load("./merged_BM_seurat.RData")
max(merged_seurat@meta.data$percent.mt)
#84.3038

# PLot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Filter based on quality plots 
Filter_BM <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt > -Inf & percent.mt < 6)
VlnPlot(Filter_BM, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filter_BM, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Sums post Filtering
sum(Filter_BM$sample=="Infected")
# 50c62
sum(Filter_BM$sample=="Not")
# c6842

# Fix Gene Names
Filter_BM@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",Filter_BM@assays$RNA@counts@Dimnames[[1]])
Filter_BM@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filter_BM@assays$RNA@counts@Dimnames[[1]])

split_seurat <- SplitObject(Filter_BM, split.by = "sample")
split_seurat <- split_seurat[c("Infected", "Not")]

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


################### Find Doublets in Infected BM
# (pk identification) 
sweep.res.list_1 <- paramSweep_v3(split_seurat[[1]], PCs= 1:10, sct=T)
sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
bcmvn_1 <- find.pK(sweep.stats_1)
#pK = 0.15

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
# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[1]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*nrow(split_seurat[[1]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.15, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
# split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.15, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.15_197", sct=T)


# PLot Doublets & Singlets
DimPlot(split_seurat[[1]], reduction = "umap", group.by = c("DF.classifications_0.25_0.15_197"),label = TRUE,repel = TRUE)
sum(split_seurat[[1]]$DF.classifications_0.25_0.15_197=="Doublet") #197

################### Doublet Finder for Non-Infected BM
# (pk identification) 
sweep.res.list_2 <- paramSweep_v3(split_seurat[[2]], PCs= 1:10, sct=T)
sweep.stats_2 <- summarizeSweep(sweep.res.list_2, GT=F)
bcmvn_2 <- find.pK(sweep.stats_2)

# pK = 0.29
######## label Pk chart
pK=as.numeric(as.character(bcmvn_2$pK))
BCmetric=bcmvn_2$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
#######scRNA_sara/scRNA_BM_sara/BM_NI_mixREF

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[2]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0524*nrow(split_seurat[[2]]@meta.data)) # Adjust according to Doublet Rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder ---- edit pK
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
# split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN ="pANN_0.25_0.29_359", sct=T)

# Plot Doublets vs. Singlets
DimPlot(split_seurat[[2]], reduction = "umap", group.by = c("DF.classifications_0.25_0.29_359"),label = TRUE,repel = TRUE)

# REMOVE DOUBLETS #
split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.15_197!="Doublet") 
dim(split_seurat[[1]])
#4865
split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset =DF.classifications_0.25_0.29_359 !="Doublet") 
dim(split_seurat[[2]])
#6483

########################################################## Downsample (4800 cells/sample)
set.seed(111)
split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 4850, replace=F)]
split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 4850, replace=F)]

for (i in 1:length(split_seurat)){
  split_seurat[[i]] <- RunTSNE(split_seurat[[i]], dims = 1:20)
}

########################################################### Integration
# select most variable features
integ_features <- SelectIntegrationFeatures(object.list=split_seurat, nfeatures=6000)

# prepare SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list=split_seurat, anchor.features=integ_features)

# Find Anchors
integ_anchors <- FindIntegrationAnchors(object.list=split_seurat, normalization.method = "SCT",
                                        anchor.features=integ_features)

# Integrate
combined.all <- IntegrateData(anchorset=integ_anchors, normalization.method="SCT")

# Save integrated seurat object & workspace
#saveRDS(combined.all, "/storage1/fs1/leyao.wang/Active/scRNA_BM_sara/SCT_integrated_BM_seurat.rds")

combined.all@meta.data

DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) # reduction = "pca",
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = 0.6)
combined.all <- RunTSNE(object = combined.all)
combined.all <- FindVariableFeatures(combined.all, selection.method = "vst", nfeatures = 2000)
combined.all[["percent.mt"]] <- PercentageFeatureSet(combined.all, pattern = "MT-",assay = "RNA")

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



########## Make tSNE and UMAP plots
DefaultAssay(combined.all) <- "integrated"
DimPlot(combined.all, reduction = "pca", split.by = "sample",label = TRUE,repel = TRUE)

DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "tsne", label = TRUE,repel = TRUE)

### For 2 samples run separately
# BM infected
DimPlot(split_seurat[[1]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[1]], reduction = "umap", label = TRUE,repel = TRUE)
# BM not infected
DimPlot(split_seurat[[2]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[2]], reduction = "umap", label = TRUE,repel = TRUE)

# Use RNA Assay for Visualization, scale and normalize raw counts
DefaultAssay(combined.all) <- "RNA"
NormalizeData(combined.all, assay = "RNA")
ScaleData(combined.all, assay = "RNA")

# Manually Annotate Clusters
combined_markers <- FindAllMarkers(object = combined.all, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25) 

comb_markers <- combined_markers[ , c(c6, 7, 1:5)]
comb_markers <- comb_markers %>%
  dplyr::arrange(cluster, p_val_adj)

# create columns for annotations
combined.all@meta.data$basic <- "NA"
combined.all@meta.data$zoom <- "NA"

#### apply greater resolution ###### sep 9 22
DefaultAssay(combined.all) <- "integrated"
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = 1.0)
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.1")
combined.all <- SetIdent(combined.all, value = "basic")
head(combined.all@meta.data)
DimPlot(combined.all, label=TRUE)
################################################
DefaultAssay(combined.all) <- "integrated"
FeaturePlot(combined.all, features = c("PROM1", "CD34", "MEIS1", "MLLT3", "TESPA1", "HLF", "MSI2", "AVP", "CRHBP"), label = TRUE, cols=c("yellow", "blue")) 
FeaturePlot(combined.all, label=TRUE, features=c("CD24", "RAG1", "VPREB1", "PAX5", "EBF1", "IGHM", "CD19", "CD79A", "DNTT"),cols=c("yellow", "blue"))
FeaturePlot(combined.all, label=TRUE, features=c("SMAD1", "TCF4", "IRF1"),cols=c("yellow", "blue"))
FeaturePlot(combined.all, label=TRUE, features=c("ESRRA", "RARA", "RELA", "STAT3"))
FeaturePlot(combined.all, label=TRUE, features=c("CD34", "CD38", 'MME', "FLT3", "TFRC", "THY1"))
FeaturePlot(combined.all, label=TRUE, features=c("IRF8", "BCL11A", "IRF4", "ESRRA", "RORB", "KLF2"))
FeaturePlot(combined.all, label=TRUE, features=c("CD34", "AVP", "HOPX")) # 11 = HSC
FeaturePlot(combined.all, label=TRUE, features=c("CEBPD","CD34", "MYADM", "KLFc6", "VIM")) #MPP?
FeaturePlot(combined.all, label=TRUE, features=c("DNTT", "CYGB", "MKIc67")) # CLP = 1,c6,1c6 + 5?
FeaturePlot(combined.all, label=TRUE, features=c("NFIL3")) #ILC = 3/5
FeaturePlot(combined.all, label=TRUE, features=c("VPREB3", "DNTT", "RAG1")) # B = 8,4,12,0,2
FeaturePlot(combined.all, label=TRUE, features=c("CEBPD", "ELANE")) # GMP = 9
FeaturePlot(combined.all, label=TRUE, features=c("CEBPD", "MKIc67", "VCAN", "S100A8", "S100A9")) #mono=9/10
FeaturePlot(combined.all, label=TRUE, features=c("BATF3")) #Cdc or mono?
FeaturePlot(combined.all, label=TRUE, features=c("VCAN", "CEBPD", "CLEC10A", "LILRA4", "SIGLECc6")) #14,7,10 =DC
FeaturePlot(combined.all, label=TRUE, features=c("GATA1", "KLF1", 'APOC1', "PF4", "TPSB2", "EPX")) # 13 = myeloid cells


# CLUSTER 10 - MYELOID/ERYTHROID CELLS
c13 <- WhichCells(combined.all, idents = "13")
c13 <- subset(combined.all, cells=c13)
AverageExpression(combined.all, slot="data", assays = "integrated", features = c("GATA1", "APOC1", "KLF1", "PLEK", "BLVRB", "TFRC", "GATA2", "FCER1A", "ITGA2B", "CA1", "HBD", "CD34", "PF4", "TPSB2", "EPX", "SIGLECc6", "MKIc67"))
FeaturePlot(c13, features = c("FCER1A","KLF1", "ITGA2B", "HBG1", "HBG2", "GP9"), label=TRUE, cols=c("yellow", "blue"))
FeaturePlot(c13, features = c( "HDC", "CA1", "GATA2", "PLEK", "BLVRB", "HBD", "GATA1", "TAL1", 'TFRC'), label=TRUE, cols=c("yellow", "blue"))
FeaturePlot(c13, features = c( "APOC1", "GATA1", "PF4", "KLF1", "TPSB2", 'EPX'), label=TRUE, cols=c("yellow", "blue"))
FeaturePlot(c13, features = c( "GATA1", "GATA2","HBD", "BLVRB", "TFRC", 'CA1', "FCER1A", "CLC", "MS4A2"), label=TRUE, cols=c("yellow", "blue"))
FeaturePlot(c13, features = c( "APOC1", "GATA1","KLF1","BLVRB", "TFRC"), label=TRUE, cols=c("yellow", "blue"))
FeaturePlot(c13, features = c( "APOC1", "GATA1","KLF1","CA1", "HBD", "ITGA2B", "FCER1A", "CD14", "CD1c6"), label=TRUE, cols=c("yellow", "blue"))
DoHeatmap(c13, features = c("APOC1", "GATA1","KLF1","CA1", "HBD", "PLEK", "BLVRB", "ITGA2B", "FCER1A", "GATA2",'EPX', "HDC", "TPSB2", "CLC", "MS4A2"), assay = "RNA", angle = 90, slot="data") + NoLegend()

mep <- WhichCells(c13, idents = "1")
memp <- WhichCells(c13, idents = "0")
ebmp <- WhichCells(c13, idents = "2")
tot <- WhichCells(c13, idents = c("1", "0", "2"))
combined.all@meta.data[tot, "basic"] <- "MEMP"
combined.all@meta.data[mep, "zoom"] <- "MEP"
combined.all@meta.data[memp, "zoom"] <- "MEMP"
combined.all@meta.data[ebmp, "zoom"] <- "Eos/Baso/Mast Prog."
DefaultAssay(c13) <- "integrated"
c13 <- RunPCA(c13, npcs = 30, verbose = FALSE)
c13 <- RunUMAP(c13, dims = 1:20) # reduction = "pca",
c13 <- FindNeighbors(c13, reduction = "pca", dims = 1:20)
c13 <- FindClusters(c13, resolution = 0.6)
head(c13@meta.data)
c13<- SetIdent(c13, value = "zoom")
DimPlot(c13, label=TRUE, label.size = 8)

#baso/eoso/mast
DefaultAssay(combined.all) <- "integrated"
FeaturePlot(combined.all, features=c("CLC", "MS4A2", "EPX"), label=TRUE, cols=c("yellow", "blue"))

# NK/ILC
FeaturePlot(combined.all, features=c("NFIL3", "IL2RB", "CD9c6"), label = TRUE) #ILC CELLS
nfil3 = GetAssayData(object = combined.all, 
                     assay = "RNA", slot = "data")["NFIL3",]
cd9c6 = GetAssayData(object = combined.all, 
                     assay = "RNA", slot = "data")["CD9c6",]
tnfa = GetAssayData(object = combined.all, 
                    assay = "RNA", slot = "data")["TNF",]   ####### TNF = macrophage, nk, t
ilc <- names(which(nfil3>3))
combined.all@meta.data[ilc, 19] <- "ILC"
DimPlot(combined.all, cells.highlight=ilc)


saveRDS(combined.all, "./combinedall914.rds")
combined.all <- readRDS("./combinedall914.rds")

# cluster 3
cl4 <- WhichCells(combined.all, idents = 4)
C4 <- subset(combined.all, cells=cl4)
# DefaultAssay(C4) <- "RNA"
DefaultAssay(C4) <- "integrated"
C4 <- RunPCA(C4, npcs = 30, verbose = FALSE)
C4 <- RunUMAP(C4, dims = 1:20) # reduction = "pca",
C4 <- FindNeighbors(C4, reduction = "pca", dims = 1:20)
C4 <- FindClusters(C4, resolution = 0.7)
DimPlot(C4, label=TRUE, label.size = 8)
head(C4@meta.data)
C4<- SetIdent(C4, value = "integrated_snn_res.0.7")
DefaultAssay(C4) <- "RNA"
FeaturePlot(C4, label=TRUE, features=c("MLLT3", "AVP", "HLF", "CRHBP", "CD34", "PROM1"), cols = c("yellow", "blue"), label.size = c6)

FeaturePlot(C4, label=TRUE, features=c("MKIc67", "CENPF", "TUBB", "TUBA1B", "MYADM", "KLFc6", "VIM"))
FeaturePlot(C4, label=TRUE, features=c("TOP2A", "SPINK2", "IRF8", "HLF"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("MLLT3", "MSI2", "BCL11A", "STAT1"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("MPO", "CEBPG", "CEBPA", "CEBPD", "AZU1", "CD38", "CD7", "MME", "SPIB", "LYZ"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("CSF3R", "KIT", "ATP1B3", "CD34", "HOPX", "DNTT", "HLF"), cols = c("yellow", "blue"))
DoHeatmap(C4, features = c("CD38", "CD34", "PROM1", "MLLT3", "MEIS1", "MSI2", "SPINK2", "FLT3", "HOPX",
                           "IRF8",  "STAT1", "BCL11A", "MPO", "LYZ", "VIM", "KIT", "CSF3R", "CD7", "MKIc67",
                           "SPIB", "CEBPG", "CEBPA", "CEBPD", "DBP", "NFIL3",  "HOXA9", "GTF2F1", "TCF3", "TCF12", "CTCF",
                           "ID3", "EGR1", "VPREB1", "CD79A", "DNTT", "CD1c64", "BST2", "KIT", "ICAM3", "TUBB", "TUBA1B", "CENPF", "ADA", "JCHAIN",
                           "ITGAM", "ACY3","EBF1"), assay = "RNA", angle = 90, slot="data") + NoLegend()

FeaturePlot(C4, label=TRUE, features=c("IRF8", "MPO", "LYZ", "CEBPD", "FLT3", "PROM1"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("SPINK2", "PROM1", "CD34", "MPO", "MLLT3"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("SPINK2", "HOPX", "DNTT", "PROM1", "CD34", "SOX4", "HLF", "TOP2A"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c( "NFIL3", "CD7", "SOX4", "MS12", "SPINK2"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("TCF3", "TCF12", "CTCF"), cols = c("yellow", "blue"))
FeaturePlot(C4, label=TRUE, features=c("NFIL3", "ACY3"), cols = c("yellow", "blue"))

c1 <- WhichCells(C4, idents = c("1"))
c2 <- WhichCells(C4, idents = c("2"))
c5 <- WhichCells(C4, idents = c("5"))
c4 <- WhichCells(C4, idents = c("4"))
cc6 <- WhichCells(C4, idents = c("c6"))
c0 <- WhichCells(C4, idents = c("0"))
c3 <- WhichCells(C4, idents = c("3"))
combined.all@meta.data[c1, "basic"] <- "MPP2"
combined.all@meta.data[c2, "basic"] <- "MPP1"
combined.all@meta.data[c3, "basic"] <- "MPP3"
combined.all@meta.data[c5, "basic"] <- "MPP4"
combined.all@meta.data[c4, "basic"] <- "MPP2"
combined.all@meta.data[cc6, "basic"] <- "MPP2"
combined.all@meta.data[c0, "basic"] <- "HSC"
combined.all@meta.data[c1, "zoom"] <- "CMP"
combined.all@meta.data[c2, "zoom"] <- "MPP"
combined.all@meta.data[c3, "zoom"] <- "LMP"
combined.all@meta.data[c5, "zoom"] <- "pre-cDCs"
combined.all@meta.data[c4, "zoom"] <- "CMP"
combined.all@meta.data[cc6, "zoom"] <- "CMP"
combined.all@meta.data[c0, "zoom"] <- "HSC"
combined.all <- SetIdent(combined.all, value = "basic")
combined.all <- SetIdent(combined.all, value = "zoom")
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.1")
DimPlot(combined.all, label=TRUE, label.size = 4)

#CLP 
clp <- WhichCells(combined.all, idents = c("1","7", "17", "5"))
CLPcells <- subset(combined.all, cells=clp)
CLPcells <- RunPCA(CLPcells, npcs = 30, verbose = FALSE)
CLPcells <- RunUMAP(CLPcells, dims = 1:20) # reduction = "pca",
CLPcells <- FindNeighbors(CLPcells, reduction = "pca", dims = 1:20)
CLPcells <- FindClusters(CLPcells, resolution = 0.5)
head(CLPcells@meta.data)
CLPcells<- SetIdent(CLPcells, value = "integrated_snn_res.0.5")
DimPlot(CLPcells, label=TRUE)
DefaultAssay(CLPcells) <- "integrated"
DefaultAssay(CLPcells) <- "RNA"
FeaturePlot(CLPcells, label=TRUE, features = c("CYGB", "TCF3", "TCF12", "CTCF", "SPINK2", "VPREB1", "SPINK2", "CD19", "DNTT", "CD79A"), cols=c("yellow", "blue"))
FeaturePlot(CLPcells, label=TRUE, features = c("CYGB", "FABP5", "ADA"), cols=c("yellow", "blue"))
DoHeatmap(B, features = c("CYGB", "TCF3", "TCF12", "CTCF", "SPINK2", "VPREB1", "SPINK2", 
                          "CD19", "DNTT", "CD79A", "FABP5", "ADA"), assay = "integrated", angle = 90) + NoLegend()



DoHeatmap(combined.all, features = c("IRF8", "MKIc67", "RUNX1", "TCF4", "SPIB", "ZBTB4c6", "IRF4", "ZEB2",
                                     "ESRRA", "KLF2", "RORB", "IRF5", "IRF1", "JUNB", "LTF", 
                                     "RETN", "S100A8", "S100A9", "MPEG1", "CYGB", "TCF3", "TCF12", "CTCF", 
                                     "SPINK2", "VPREB1", "SPINK2", "CD19", "DNTT", "CD79A", "GNLY", "FCGR3A", "STAT1", "CD34", "CD79A", "CD24", "RAG1", "VPREB3", "VPREB1", "DNTT", "PAX5", "IGHM", "MME", "HVCN1",
                                     "SMAD1", "TCF4", "CYGB", "TCF3", "CD19", "TCL1A", "SELL", "BCL2", "CD40",  "EBF1",
                                     "IGHD", "STMN1", "MKIc67", "MS4A1",
                                     "CYGB", "TCF3", "TCF12", "CTCF", "SPINK2", "VPREB1", "SPINK2", 
                                     "CD19", "DNTT", "CD79A", "FABP5", "ADA", "NFIL3", "CD74", "MNDA", "LYZ", "LYST"), assay = "RNA", slot="data", angle = 90) + NoLegend()

####################################################################### Analysis
combined.all <- SetIdent(combined.all, value = "zoom")
combined.all <- SetIdent(combined.all, value = "basic")
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.1")
DimPlot(combined.all, label=TRUE, label.size = 4)
DimPlot(combined.all, cells= c13, label=TRUE, label.size = 7)
table(combined.all@meta.data$basic)
table(combined.all@meta.data$zoom)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
split_seu <- SplitObject(combined.all, split.by = "sample")
split_seu <- split_seu[c("Infected", "Not")]
table(split_seu[[1]]@meta.data$basic)
table(split_seu[[2]]@meta.data$basic)

DimPlot(C4, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
split_c4 <- SplitObject(C3, split.by = "sample")
split_c4 <- split_c3[c("Infected", "Not")]
table(split_c3[[1]]@meta.data$basic)
table(split_c3[[2]]@meta.data$basic)

table(combined.all@meta.data$celltype.stim)

genes_intersted <- readxl::read_xlsx("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/SCcluster_genes.xlsx", col_names = F)
colnames(genes_intersted) <- "Genes"
RidgePlot(combined.all,features= c("CD4", "CD8A"),assay = "RNA", slot= "data", stack = TRUE) + RotatedAxis()
DotPlot(combined.all, assay = "RNA",  features=  genes_intersted$Genes[1:23], cols = c("blue", "red"), dot.scale=8, split.by = "sample") + RotatedAxis()
DotPlot(combined.all, assay = "RNA",  features = genes_intersted$Genes[24:49], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()

DoHeatmap(combined.all, features = genes_intersted$Genes, assay = "RNA", slot="data", angle = 90) + NoLegend()

combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all@meta.data$celltype.stim)

######################################## Fine Differential EXpressed Genes ############################################################
proGran_DEgenes <- FindMarkers(combined.all, ident.1 = "pro-Granulocytes_Infected", ident.2 = "pro-Granulocytes_Not", verbose = FALSE)

myeDC_DEgenes <- FindMarkers(combined.all, ident.1 = "Mye-DC_Infected", ident.2 = "Mye-DC_Not", verbose = FALSE)

premono_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-Monocytes_Infected", ident.2 = "pre-Monocytes_Not", verbose = FALSE)

mono_DEgenes <- FindMarkers(combined.all, ident.1 = "Monocytes_Infected", ident.2 = "Monocytes_Not", verbose = FALSE)

predcs_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-DC Cycle_Infected", ident.2 = "pre-DC Cycle_Not", verbose = FALSE)

predcsCYC_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-DCs_Infected", ident.2 = "pre-DCs_Not", verbose = FALSE)

MEMP_DEgenes <- FindMarkers(combined.all, ident.1 = "MEMP_Infected", ident.2 = "MEMP_Not", verbose = FALSE)

seven_DEgenes <- FindMarkers(combined.all, ident.1 = "7_Infected", ident.2 = "7_Not", verbose = FALSE)

mpp1_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP1_Infected", ident.2 = "MPP1_Not", verbose = FALSE)

mpp2_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP2_Infected", ident.2 = "MPP2_Not", verbose = FALSE)

mpp3_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP3_Infected", ident.2 = "MPP3_Mock", verbose = FALSE)

mpp4_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP4_Infected", ident.2 = "MPP4_Mock", verbose = FALSE)

hsc_DEgenes <- FindMarkers(combined.all, ident.1 = "HSC_Infected", ident.2 = "HSC_Mock", verbose = FALSE)

cdp_DEgenes <- FindMarkers(combined.all, ident.1 = "CDP_Infected", ident.2 = "CDP_Mock", verbose = FALSE)

lmpp_DEgenes <- FindMarkers(combined.all, ident.1 = "LMPP_Infected", ident.2 = "LMPP_Mock", verbose = FALSE)

clp1_DEgenes <- FindMarkers(combined.all, ident.1 = "CLP.1_Infected", ident.2 = "CLP.1_Not", verbose = FALSE)

clp2_DEgenes <- FindMarkers(combined.all, ident.1 = "CLP.2_Infected", ident.2 = "CLP.2_Not", verbose = FALSE)

proB_DEgenes <- FindMarkers(combined.all, ident.1 = "Pro-B_Infected", ident.2 = "Pro-B_Not", verbose = FALSE)

preProB_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-Pro-B_Infected", ident.2 = "pre-Pro-B_Not", verbose = FALSE)

unident_DEgenes <- FindMarkers(combined.all, ident.1 = "unidentified_Infected", ident.2 = "unidentified_Mock", verbose = FALSE)
unident_marks <- FindMarkers(combined.all, ident.1 = "unidentified")
saveRDS(unident_marks, "./unident_cluster_marks.rds")


##################################################################################################
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
preM <- subset(combined.all, idents = "pre-Monocytes")
Idents(preM) <- "sample"
avg.preM <- as.data.frame(AverageExpression(preM, verbose = FALSE, slot = "data")$RNA)
avg.preM$gene <- rownames(avg.preM)
genes.to.label = c("CXCL10", "IFI27", "B2M", "ISG15", "IFITM3", "GBP1", "IF1c6", "LYc6E", "STAT1", "MX1", "SNRPG", "SNHG25")
p1 <- ggplot(avg.preM, aes(Infected, Not)) + geom_point() + ggtitle("Pre-Monocytes")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

Mono <- subset(combined.all, idents = "Monocytes")
Idents(Mono) <- "sample"
avg.Mono <- as.data.frame(AverageExpression(Mono, verbose = FALSE, slot="data")$RNA)
avg.Mono$gene <- rownames(avg.Mono)
genes.to.label = c("CXCL10", "IFI27", "B2M", "ISG15", "ISG20", 'IRF8', "GBP1", "CXCL9", "SCT", "PSME2")
p1 <- ggplot(avg.Mono, aes(Infected, Not)) + geom_point() + ggtitle("Monocytes")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

preDC <- subset(combined.all, idents = "pre-DCs")
Idents(preDC) <- "sample"
avg.preDC <- as.data.frame(AverageExpression(preDC, verbose = FALSE, slot = "data")$RNA)
avg.preDC$gene <- rownames(avg.preDC)
genes.to.label = c("PSMB9", "IFI27", "B2M", "IFI44L", "ISG15", "HNRNPH1", "ISG20", "LYc6E", "TYMP", "BST2")
p1 <- ggplot(avg.preDC, aes(Infected, Not)) + geom_point() + ggtitle("preDC")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


proB <- subset(combined.all, idents = "Pro-B")
Idents(proB) <- "sample"
avg.proB <- as.data.frame(AverageExpression(proB, verbose = FALSE, slot = "counts", assay="RNA")$RNA)
avg.proB$gene <- rownames(avg.proB)
genes.to.label = c("MALAT1", "EBF1", "BACH2", "SSBP2", "IFI27", "MT2A", "ISG15", "TEX14", "MX1")
p1 <- ggplot(avg.proB, aes(Infected, Not)) + geom_point() + ggtitle("proB")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

proB <- subset(combined.all, idents = "Pro-B")
Idents(proB) <- "sample"
avg.proB <- as.data.frame(log10(AverageExpression(proB, verbose = FALSE)$RNA))
avg.proB$gene <- rownames(avg.proB)
genes.to.label = c("MALAT1", "EBF1", "BACH2", "SSBP2", "IFI27", "MT2A", "ISG15", "TEX14", "MX1")
p1 <- ggplot(avg.proB, aes(Infected, Not)) + geom_point() + ggtitle("proB")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

Memp <- subset(combined.all, idents = "MEMP")
Idents(Memp) <- "sample"
avg.Memp <- as.data.frame(AverageExpression(Memp, verbose = FALSE, slot = "data")$RNA)
avg.Memp$gene <- rownames(avg.Memp)
genes.to.label = c("IFI27", "HBA2", "HBA1", "HBG1", "HBD", "IFITM3", "ISG15", "IF1c6", "MT2A", 'TMSB4X', "PSMB9")
p1 <- ggplot(avg.Memp, aes(Infected, Not)) + geom_point() + ggtitle("Memp")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

na <- subset(combined.all, idents = "NA")
Idents(na) <- "sample"
avg.na <- as.data.frame(log10(AverageExpression(na, verbose = FALSE)$RNA))
avg.na$gene <- rownames(avg.na)
genes.to.label = c("TMSB4X", "B2M", "ISG15", "HMGB1", "FCHSG2")
p1 <- ggplot(avg.na, aes(Infected, Not)) + geom_point() + ggtitle("na")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

proGran <- subset(combined.all, idents = "pro-Granulocytes")
Idents(proGran) <- "sample"
avg.proGran <- as.data.frame(AverageExpression(proGran, verbose = FALSE, slot = "data")$RNA)
avg.proGran$gene <- rownames(avg.proGran)
genes.to.label = c("IFI27", "ISG15", "MT2A", "IF1c6", "IFI44L", "IFITM3", "STAT1", "PSMB9", "KLHDC7B", "CXCL10", "LYc6E", "IFIT3")
p1 <- ggplot(avg.proGran, aes(Infected, Not)) + geom_point() + ggtitle("pro-Granulocyte")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

hsc <- subset(combined.all, idents = "HSC")
Idents(hsc) <- "sample"
avg.hsc <- as.data.frame(AverageExpression(hsc, verbose = FALSE, slot = "data")$RNA)
avg.hsc$gene <- rownames(avg.hsc)
genes.to.label = c("MT2A", "STAT1", "IFI27", "IF1c6", "CXCL10", "ISG15", "SAT1", "IFIT3", "IFITM3", "ISG20", "CXCL11", "IRF7")
p1 <- ggplot(avg.hsc, aes(Infected, Not)) + geom_point() + ggtitle("HSC")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

mpp1 <- subset(combined.all, idents = "MPP1")
Idents(mpp1) <- "sample"
avg.mpp1 <- as.data.frame(AverageExpression(mpp1, verbose = FALSE, slot="data")$RNA)
avg.mpp1$gene <- rownames(avg.mpp1)
genes.to.label = c("IFI27", "ISG15", "STAT1", "IF1c6", "CXCL10", "MT2A", "B2M", "IFIT3",
                   "IFITM3",  "GBP1", "IFI44L", "MX1")
p1 <- ggplot(avg.mpp1, aes(Infected, Not)) + geom_point() + ggtitle("MPP1")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

mpp2 <- subset(combined.all, idents = "MPP2")
Idents(mpp2) <- "sample"
avg.mpp2 <- as.data.frame(AverageExpression(mpp2, verbose = FALSE, slot="data")$RNA)
avg.mpp2$gene <- rownames(avg.mpp2)
genes.to.label = c("IFI27", "ISG15", "MT2A", "IF1c6", "CXCL10", "GBP1", "IFITM3", "B2M", "STAT1", "IFIT3", 'PSMB9', "HNRNPH1", "IFI44L")
p1 <- ggplot(avg.mpp2, aes(Infected, Not)) + geom_point() + ggtitle("MPP2")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

mpp3 <- subset(combined.all, idents = "MPP3")
Idents(mpp3) <- "sample"
avg.mpp3 <- as.data.frame(AverageExpression(mpp3, verbose = FALSE)$RNA)
avg.mpp3$gene <- rownames(avg.mpp3)
genes.to.label = c("IFI27", "ISG15", "STAT1", "IF1c6", "CXCL10", "IFIT3", "IFI44L", "IFITM3", "IFIT3", "IFIT1", "MT2A", "PSMB9")
p1 <- ggplot(avg.mpp3, aes(Infected, Mock)) + geom_point() + ggtitle("MPP3")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

mpp4 <- subset(combined.all, idents = "MPP4")
Idents(mpp4) <- "sample"
avg.mpp4 <- as.data.frame(AverageExpression(mpp4, verbose = FALSE)$RNA)
avg.mpp4$gene <- rownames(avg.mpp4)
genes.to.label = c("IFI27", "ISG15", "MT2A", "IF1c6", "CXCL10", "ISG20", "SNHG25", "HLA-DRB1", "STAT1")
p1 <- ggplot(avg.mpp4, aes(Infected, Mock)) + geom_point() + ggtitle("MPP4")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

cdp <- subset(combined.all, idents = "CDP")
Idents(cdp) <- "sample"
avg.cdp <- as.data.frame(AverageExpression(cdp, verbose = FALSE)$RNA)
avg.cdp$gene <- rownames(avg.cdp)
genes.to.label = c("IFI27", "ISG15", "IFI44L", "IF1c6", "PARP14", "IFITM3", "STAT1", "IFI44", "MX1")
p1 <- ggplot(avg.cdp, aes(Infected, Mock)) + geom_point() + ggtitle("CDP")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

clp1 <- subset(combined.all, idents = "CLP.1")
Idents(clp1) <- "sample"
avg.clp1 <- as.data.frame(log10(AverageExpression(clp1, verbose = FALSE)$RNA))
avg.clp1$gene <- rownames(avg.clp1)
genes.to.label = c("IFI27", "ISG15", "MALAT1", "IF1c6", "CXCL10")
p1 <- ggplot(avg.clp1, aes(Infected, Not)) + geom_point() + ggtitle("CLP.1")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

clp2 <- subset(combined.all, idents = "CLP.2")
Idents(clp2) <- "sample"
avg.clp2 <- as.data.frame(log10(AverageExpression(clp2, verbose = FALSE)$RNA))
avg.clp2$gene <- rownames(avg.clp2)
genes.to.label = c("IFI27", "ISG15", "MALAT1", "IF1c6", "CXCL10")
p1 <- ggplot(avg.clp2, aes(Infected, Not)) + geom_point() + ggtitle("CLP.2")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

clp3 <- subset(combined.all, idents = "CLP.3")
Idents(clp3) <- "sample"
avg.clp3 <- as.data.frame(log10(AverageExpression(clp3, verbose = FALSE)$RNA))
avg.clp3$gene <- rownames(avg.clp3)
genes.to.label = c("IFI27", "ISG15", "MALAT1", "IF1c6", "CXCL10")
p1 <- ggplot(avg.clp3, aes(Infected, Not)) + geom_point() + ggtitle("CLP.3")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

#######################################################################################

# 10/3
saveRDS(combined.all, "./combinedALL.rds")
combined.all <- readRDS("./combinedALL.rds")


######################################################################################################### Modifications 10/20/22
# CLP.1 
clp1 <- WhichCells(combined.all, idents = c("CLP.1"))
combined.all@meta.data[clp1, "basic"] <- "LMPP"
combined.all <- SetIdent(combined.all, value = "zoom")
DimPlot(combined.all, label=T)

# MEMP
memp <- WhichCells(combined.all, idents = "MEMP")
Memp <- subset(combined.all, cells=memp)
Memp<- SetIdent(Memp, value = "zoom")
DimPlot(Memp, label=TRUE, label.size = 8)

#MEMP & MPP1
memp_mpp1 <- WhichCells(combined.all, idents = c("MEMP", "MPP1"))
Memp_mpp1<- SetIdent(Memp_mpp1, value = "zoom")
Memp_mpp1 <- subset(combined.all, cells = memp_mpp1)
DoHeatmap(Memp_mpp1, features = c("CD38", "CD34", "PROM1", "MLLT3", "MEIS1", "MSI2", "SPINK2",
                                  "APOC1", "GATA1","KLF1","CA1", "HBD", "PLEK", "BLVRB", "ITGA2B",
                                  "FCER1A", "GATA2","HDC", "TPSB2", "CLC", "MS4A2", "EPX"), assay = "RNA", angle = 90, slot="data") + NoLegend()
DefaultAssay(Memp_mpp1) <- "integrated"
Memp_mpp1 <- RunPCA(Memp_mpp1, npcs = 30, verbose = FALSE)
Memp_mpp1 <- RunUMAP(Memp_mpp1, dims = 1:20) # reduction = "pca",
Memp_mpp1 <- FindNeighbors(Memp_mpp1, reduction = "pca", dims = 1:20)
Memp_mpp1 <- FindClusters(Memp_mpp1, resolution = 0.6)
Memp_mpp1 <- SetIdent(Memp_mpp1, value="integrated_snn_res.0.6")
DimPlot(Memp_mpp1, label = T)
memp_mpp1_marks <- FindAllMarkers(Memp_mpp1)

# pre-Pro-B.1
bcells <- WhichCells(combined.all, idents = c("pre-Pro-B cycle", "pre-Pro-B.1", "pre-Pro-B.2", "Pro-B"))
PROB <- subset(combined.all, cells = bcells)
DimPlot(PROB, label=T)
PROB_marks <- FindAllMarkers(PROB)
write.csv(PROB_marks, "./PROBmarks.csv")



DefaultAssay(PROB) <- "integrated"
PROB <- RunPCA(PROB, npcs = 30, verbose = FALSE)
PROB <- RunUMAP(PROB, dims = 1:20) # reduction = "pca",
PROB <- FindNeighbors(PROB, reduction = "pca", dims = 1:20)
PROB <- FindClusters(PROB, resolution = 0.3)
PROB <- SetIdent(PROB, value="integrated_snn_res.0.3")
DimPlot(PROB, label = T, label.size=8)
PROB_marks <- FindAllMarkers(PROB)

DoHeatmap(PROB, features = c("CD79A", "CD24", "RAG1", "VPREB1", "DNTT", "PAX5", "IGHM", "MME", "HVCN1",
                             "SMAD1", "TCF4", "CYGB", "TCF3", "CD19", "TCL1A", "SELL", "CD40", "CD19", "EBF1",
                             "IGHD", "STMN1", "MKIc67", "MS4A1", "IRF4", 'ACSM3', "LTB", "HMGB1", 'SPINK2',
                             "CD34", "RAG2"), assay = "RNA", slot="data", angle = 90) + NoLegend()


# Prior: "combinedALL.rds"
# Oct 24
saveRDS(combined.all, "./combined_all.rds")
combined.all <- readRDS("./combined_all.rds")

# edits at Shuai request
combined.all <- SetIdent(combined.all, value = "basic")
mpp3_4 <- WhichCells(combined.all, idents = "MPP3")
mpp4_3 <- WhichCells(combined.all, idents = "MPP4")
combined.all@meta.data[mpp3_4, "basic"] <- "MPP4"
combined.all@meta.data[mpp4_3, "basic"] <- "MPP3"
lmpp <- WhichCells(combined.all, idents = "CLP.1")
combined.all@meta.data[lmpp, "basic"] <- "LMPP"
clp1 <- WhichCells(combined.all, idents = "CLP.2")
combined.all@meta.data[clp1, "basic"] <- "CLP.1"
clp2 <- WhichCells(combined.all, idents = "CLP.3")
combined.all@meta.data[clp2, "basic"] <- "CLP.2"


DimPlot(combined.all, label = T, pt.size = 1.5, cols = c("aquamarine2", "blue1", "cadetblue1", "chartreuse", "paleturquoise3", "cyan2", "darkcyan", "darkgoldenrod1", 
                                                         "darkolivegreen2", "darkorchid", "deeppink1", "deepskyblue1", "forestgreen", "firebrick1", "mediumpurple1", 
                                                         "lightgoldenrod1", "lightpink1", "lightsalmon2", "lightslateblue", "plum2", "chocolate1", "yellow3", "violet"))

allmarks <- FindAllMarkers(combined.all, assay = "RNA", slot = "data")

# Use log normalized expression for visualization and DESeq
combined.all <- NormalizeData(combined.all, assay = "RNA")
combined.all <- ScaleData(combined.all, assay = "RNA")

GENEexpr <- as.data.frame(AverageExpression(combined.all, features = c("IFNGR1", "IFNAR1", "CD34", "CD38", "CD7", "MME", "CD19", "THY1", "FLT3", "TFRC"), verbose = FALSE, slot = "data", assays = "RNA"))
write.csv(GENEexpr, "./Marker_Gene_Expr_unseparated.csv")
cd1c6c6 <- as.data.frame(AverageExpression(combined.all, features = c("ALCAM"), verbose = FALSE, slot = "data", assays = "RNA"))
write.csv(cd1c6c6, "./cd1c6c6_expression.csv")

# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c("Pro-B_Infected", "Pro-B_Mock", "pre-Pro-B.1_Infected", "pre-Pro-B.1_Mock", "pre-Pro-B.2_Infected", "pre-Pro-B.2_Mock",
               "pre-Pro-B cycle_Infected", "pre-Pro-B cycle_Mock", "CLP.1_Infected", "CLP.1_Mock", "CLP.2_Infected","CLP.2_Mock",
               "CLP cycle_Infected", "CLP cycle_Mock", "LMPP_Infected", "LMPP_Mock", "MPP4_Infected", "MPP4_Mock", "MPP3_Infected", "MPP3_Mock",
               "MPP2_Infected", "MPP2_Mock", "MPP1_Infected", "MPP1_Mock", "HSC_Infected", "HSC_Mock", "CDP_Infected", "CDP_Mock",
               "MEMP_Infected", "MEMP_Mock", "Pre-Monocytes_Infected", "Pre-Monocytes_Mock", "MP_Infected", "MP_Mock", 
               "Mye-DC_Infected", "Mye-DC_Mock", "pre-DCs_Infected", "pre-DCs_Mock", "pre-DC Cycle_Infected", "pre-DC Cycle_Mock",
               "GMP_Infected", "GMP_Mock", "unidentified_Infected", "unidentified_Mock")

# Re-level object@ident
combined.all@active.ident <- factor(x = combined.all@active.ident, levels = my_levels)
DoHeatmap(combined.all, features = c("ALCAM"), assay = "RNA", angle = 90, slot="scale.data") + NoLegend()

DimPlot(combined.all, label = F, pt.size = 1)

# Change names
NAs <- WhichCells(combined.all, idents = "NA")
combined.all@meta.data[NAs, "basic"] <- "unidentified"

combined.all <- SetIdent(combined.all, value = "main")
not <- WhichCells(combined.all, idents = "Not")
combined.all@meta.data[not, "sample"] <- "Mock"
FeaturePlot(combined.all, features = "ALCAM", cols =c("grey", "red"), split.by = "sample", label = TRUE)
RidgePlot(combined.all, features = "ALCAM")


##################################### PLot Cell Counts
ggplot(combined.all@meta.data, aes(x= basic)) + 
  geom_bar(aes(y = (..count..), fill = sample), position = "dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

DefaultAssay(combined.all) <- "RNA"
FeaturePlot(combined.all, features = c("ELANE", "MPO", "AZU1", "CLEC11A"), label = T, label.size = 2, slot = "data")


### Shuai DEG - Feb 1 2023
VlnPlot(combined.all, features = c("CASP1", "CASP3", "CASP7", "CASP8"), split.by = "sample", 
        pt.size = 0, combine = FALSE, assay = "RNA", slot = "scale.data")
FeaturePlot(combined.all, label = T, cols=c("yellow","blue"), features = c("CASP1"))
avgs <- AverageExpression(combined.all, features = c("CASP1", "CASP3", "CASP7", "CASP8", "GSDMD", "GSDME", "RIPK1", "RIPK3", "MLKL", "ZBP1", "FADD", "STAT1", "IL1A"), slot="data", assays = "RNA")
DoHeatmap(combined.all, features = c("CASP1", "CASP3", "CASP7", "CASP8", "GSDMD", "GSDME", "RIPK1", "RIPK3", "MLKL", "ZBP1", "FADD", "STAT1", "IL1A"), assay = "RNA", slot="scale.data", angle = 90) + NoLegend()

write.csv(as.data.frame(avgs), "/storage1/fs1/liang.shan/Active/Average_Gene_Expression.csv")

DoHeatmap(combined.all, features = c("HLF", "AVP", "CD34", "PROM1", "TOP2A", "MKIc67", "S100A8", "S100A9",
                                     "MPO", "PRTN3", "ELANE", "LTB", "JCHAIN", "GATA2", "ITGA2B", 
                                     "KLF1", "TFRC", "HBG2", "AHSP", "PF4", "ITGB3", "KIT", "KRT1"), assay = "RNA", slot="data", angle = 90) + NoLegend()

combined.all <- SetIdent(combined.all, value = "main")
DimPlot(combined.all, label = T, cols = colors, label.size = 5, pt.size = 3, split.by = "sample")
DimPlot(combined.all, label=T)
mep <- WhichCells(combined.all, idents = "MEP")
memp <- WhichCells(combined.all, idents = "MEMP")
ebmp <- WhichCells(combined.all, idents = "Eos/Baso/Mast Prog.")
combined.all$main <- combined.all$basic
combined.all@meta.data[mep, "main"] <- "MEP"
combined.all@meta.data[memp, "main"] <- "MEMP"
combined.all@meta.data[ebmp, "main"] <- "Eos/Baso/Mast Prog."

gmp <- WhichCells(combined.all, idents = "pro-Granulocytes")
mp <- WhichCells(combined.all, idents = "pre-Monocytes")
mon <- WhichCells(combined.all, idents = "Monocytes")
combined.all$main <- combined.all$basic
combined.all@meta.data[gmp, "main"] <- "GMP"
combined.all@meta.data[mp, "main"] <- "MP"
combined.all@meta.data[mon, "main"] <- "Pre-Monocytes"

# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c("HSC", "MPP1", "MPP2", "MPP3", "MPP4", "LMPP", "CLP.1", "CLP.2", "CLP cycle", "pre-Pro-B.1", "pre-Pro-B.2",
               "pre-Pro-B cycle", "Pro-B", "CDP", "Mye-DC", "pre-DCs", "pre-DC Cycle",   "GMP", "MP",
               "pre-Monocytes", "MEMP", "unidentified")

# Re-level object@ident
combined.all@active.ident <- factor(x = combined.all@active.ident, levels = my_levels)

DoHeatmap(combined.all, features = c('CD38','CD34','AVP', 'PROM1', 'MPO','CEBPD',
                                     'HOPX', 'DNTT',
                                     'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKIc67', 
                                     'VPREB3','RAG1','HVCN1',
                                     'ACY3', 'IRF8','BATF3','LILRA4', 'FCER1A','ELANE', 'S100A8', 'S100A9','VCAN',
                                     'GATA1', 'GATA2','KLF1','APOC1','TPSB2','EPX'), assay = "RNA", slot="data", angle = 90) + NoLegend()

colors = c('#E0E0E0', '#CC9966', '#DC143C', '#FF7F50', '#FFA500', '#FFFF00', '#9ACD32', '#6B8E23',
           '#98FB98', '#20B2AA', '#00FFFF','#1E90FF', '#87CEFA', '#0000CD', '#8A2BE2', '#7B68EE', '#BA55D3',
           '#DB7093', '#FF69B4', '#FFB6C1', '#FAEBD7', '#A0522D', '#BC8F8F', '#B0C4DE')


p = VlnPlot(combined.all,c('CD38','CD34','AVP', 'PROM1', 'MPO','CEBPD', 'DNTT',
                           'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKI67', 
                           'VPREB3','RAG1', 'ACY3', 'IRF8','BATF3','LILRA4', 'FCER1A','ELANE', 'S100A8', 'S100A9','VCAN',
                           'GATA1', 'GATA2','APOC1','EPX'),
            pt.size=0,stack=T,flip=T,fill.by='ident',assay='RNA')+
  NoLegend()+
  scale_x_discrete(position='top')+
  scale_fill_manual(values=colors)+
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


test = subset(data,idents=c('MPP.c1','GMP'))
markers = FindAllMarkers(test,only.pos=T,assay='RNA')
top = markers %>% group_by(cluster) %>% top_n(n=20,wt=avg_log2FC)
DoHeatmap(test,features=top$gene,assay='RNA')+NoLegend()

# SCpubR
library(SCpubr)
SCpubr::do_DimPlot(sample = combined.all, pt.size =1.5, label = T, label.size = 3, label.box = T, split.by = "sample", group.by = "main")
SCpubr::do_FeaturePlot(sample = combined.all,
                       features = c("AVP", "MPO"), label.color = "black", label = T)
combined.all <- SetIdent(combined.all, value = "sample")
de_genes <- FindMarkers(combined.all, ident.1 = "Infected", ident.2 = "Mock")
SCpubr::do_VolcanoPlot(sample = combined.all,
                       de_genes = de_genes,
                       n_genes = 15,
                       pval_cutoff = 1.3,
                       FC_cutoff = 1.5,
                       colors.use = c("red"),
                       font.size=12)
combined.all <- SetIdent(combined.all, value = combined.all$celltype.stim)
SCpubr::do_ExpressionHeatmap(sample = combined.all,
                             features = c("MKI67"), 
                             viridis_direction = -1,
                             # group.by = "sample",
                             flip= T, cell_size = 7)

##### MONACLE3
# TUTORIAL : https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
# devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(tidyverse)

tol_high_contrast_palette <- c("#DDAA33", "#BB55c6c6", "#004488")
tol_vibrant_palette <- c("#0077BB", "#33BBEE", "#009988",
                         "#EE7733", "#CC3311", "#EE3377",
                         "#BBBBBB")
tol_muted_palette <- c("#332288", "#88CCEE", "#44AA99",
                       "#117733", "#999933", "#DDCC77",
                       "#CCc6c677", "#882255", "#AA4499")

cds <- as.cell_data_set(combined.all)
cds <- reduce_dimension(cds)
plot_cells(cds)
cds <- cluster_cells(cds, resolution=1e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds1 <- cds
cds1 <- preprocess_cds(cds1, method = "PCA", norm_method = "none")
cds1 <- reduce_dimension(cds1, reduction_method = "tSNE", preprocess_method = "PCA")
cds1 <- cluster_cells(cds1, reduction_method = "tSNE")
p1 <- plot_cells(cds1, color_cells_by = "cluster", show_trajectory_graph = FALSE, reduction_method = "tSNE")
p2 <- plot_cells(cds1, color_cells_by = "main", show_trajectory_graph = F, reduction_method = "tSNE")
wrap_plots(p1, p2)

cds1 <- reduce_dimension(cds1, reduction_method = "PCA", preprocess_method = "PCA")
cds1 <- cluster_cells(cds1, reduction_method = "PCA")
p1 <- plot_cells(cds1, color_cells_by = "cluster", show_trajectory_graph = FALSE, reduction_method = "PCA")
p2 <- plot_cells(cds1, color_cells_by = "partition", show_trajectory_graph = FALSE, reduction_method = "PCA")
wrap_plots(p1, p2)


integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)


gene_modules <- find_gene_modules(cds[deg_ids,],
                                  resolution=c(10^seq(-c6,-1)))
table(gene_modules$module)

cell_groups <- data.frame(cell = row.names(colData(cds)),
                          cell_group = colData(cds)$orig.ident)
agg_mat <- aggregate_gene_expression(cds,
                                     gene_group_df = gene_modules,
                                     cell_group_df = cell_groups)
dim(agg_mat)

row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   treeheight_row = 0,
                   treeheight_col = 0,
                   clustering_method="ward.D2")

cluster_groups <- data.frame(cell = row.names(colData(cds)),
                             cluster_group = cds@clusters$UMAP[[2]])
agg_mat2 <- aggregate_gene_expression(cds, gene_modules, cluster_groups)
dim(agg_mat2)


row.names(agg_mat2) <- paste0("Module ", row.names(agg_mat2))
pheatmap::pheatmap(agg_mat2,
                   scale="column",
                   treeheight_row = 0,
                   treeheight_col = 0,
                   clustering_method="ward.D2")


gm <- gene_modules[which(gene_modules$module %in% c(8, 18)),]
plot_cells(cds,
           genes=gm,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           trajectory_graph_color = "greyc60")



## Slingshot Trajectory Analysis
# https://bustools.github.io/BUS_notebooks_R/slingshot.html#qc
library(slingshot)
library(BiocParallel)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(Matrix)

combined.all <- SetIdent(combined.all, value = "basic")
DimPlot(combined.all, label = T, split = "sample")

DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 100, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) 
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = c(0.6, 0.8,1.0, 1.3, 1.6, 1.8, 2))
#combined.all <- RunTSNE(object = combined.all)
combined.all <- FindVariableFeatures(combined.all, selection.method = "vst", nfeatures = 2000)
#combined.all <- SetIdent(combined.all, value = combined.all$integrated_snn_res.2)
DimPlot(combined.all, label = T)

c16 <- WhichCells(combined.all, idents = "16")

combined.all <- SetIdent(combined.all, value = combined.all$main)
DimPlot(combined.all, label = T)

hsc <- WhichCells(combined.all, idents = "HSC")

combined.all <- SetIdent(combined.all, value = combined.all$integrated_snn_res.2)
DimPlot(combined.all, label = T)

combined.all$traj.clusters <- as.character(combined.all$integrated_snn_res.0.6)

combined.all@meta.data[c16, "traj.clusters"] <- "14"
combined.all@meta.data[hsc, "traj.clusters"] <- "15"

combined.all$traj.clusters <- as.factor(combined.all$traj.clusters)
combined.all <- SetIdent(combined.all, value = combined.all$traj.clusters)
DimPlot(combined.all, label = T)

mp <- WhichCells(combined.all, idents = "pre-Monocytes")
premono <- WhichCells(combined.all, idents = "Monocytes")
clp2 <- WhichCells(combined.all, idents = "CLP.3")
clp1 <- WhichCells(combined.all, idents = "CLP.2")
lmpp <- WhichCells(combined.all, idents = "CLP.1")
unident <- WhichCells(combined.all, idents = "NA")
gmp <- WhichCells(combined.all, idents = "pro-Granulocytes")

combined.all@meta.data[mp, "basic"] <- "MP"
combined.all@meta.data[premono, "basic"] <- "pre-Monocytes"
combined.all@meta.data[clp2, "basic"] <- "CLP.2"
combined.all@meta.data[clp1, "basic"] <- "CLP.1"
combined.all@meta.data[lmpp, "basic"] <- "LMPP"
combined.all@meta.data[unident, "basic"] <- "unidentified"
combined.all@meta.data[gmp, "basic"] <- "GMP"

combined.all$trajectory_clusters <- combined.all$integrated_snn_res.1.3
mpp1 <- WhichCells(combined.all, idents = "MPP1")
mpp2 <- WhichCells(combined.all, idents = "MPP2")
mpp3 <- WhichCells(combined.all, idents = "MPP3")
mpp4 <- WhichCells(combined.all, idents = "MPP4")
hcs <- WhichCells(combined.all, idents = "HSC")
lmpp <- WhichCells(combined.all, idents = "LMPP")
combined.all$trajectory_clusters <- as.character(combined.all$trajectory_clusters)
combined.all@meta.data[mpp1, "trajectory_clusters"] <- "mpp1"
combined.all@meta.data[mpp2, "trajectory_clusters"] <- "mpp2"
combined.all@meta.data[mpp3, "trajectory_clusters"] <- "mpp3"
combined.all@meta.data[mpp4, "trajectory_clusters"] <- "mpp4"
combined.all@meta.data[hcs, "trajectory_clusters"] <- "hsc"
combined.all@meta.data[lmpp, "trajectory_clusters"] <- "LMPP"
combined.all$trajectory_clusters <- as.factor(combined.all$trajectory_clusters )
combined.all <- SetIdent(combined.all, value = combined.all$trajectory_clusters)
DimPlot(combined.all, label = T)
table(combined.all$trajectory_clusters)

# twelve <- subset(combined.all, idents = 12)
# t12 <- WhichCells(combined.all, idents = 12)
# 
# combined.all@meta.data[t12, "trajectory_clusters"] <- 0


sds <- slingshot(Embeddings(combined.all, "umap"), clusterLabels = combined.all$main, 
                 start.clus = "HSC", stretch = 0)

sds <- SlingshotDataSet(sds)
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

cell_colors_clust <- cell_pal(combined.all$integrated_snn_res.2, hue_pal())

plot(reducedDim(sce), col = cell_colors_clust, pch = 10, cex = 0.4)
lines(sce, lwd = 2, type = 'lineages', col = 'black', show.constraints = T)

plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')

## Differential Expression 
# Get top highly variable genes
top_hvg <- HVFInfo(combined.all) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(300, residual_variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(combined.all, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use) # Do curve 2, so 2nd columnn
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
  plot(reducedDim(sds), col = colors, 
       pch = 16, cex = 0.5, main = top_genes[i])
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}


save.image("./slingshot_traj_analysis.RData")


BiocManager::install("tradeSeq")
library(tradeSeq)










#### Slingshot https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#loading_data
# Save the objects as separate matrices for input in slingshot
dimred <- combined.all@reductions$umap@cell.embeddings
clustering <- combined.all$integrated_snn_res.2
counts <- as.matrix(combined.all@assays$RNA@counts[combined.all@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)

lineages <- SlingshotDataSet(lineages)

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
# Plot the lineages
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")

plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')

nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}




########################## TRAJ ANALYSIS MARCH 30, 2023 - Apr 10
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#:~:text=In%20this%20workshop%2C%20we%20will%20explore%20methods%20for,can%20detect%20large-scale%20changes%2C%20indicative%20of%20differential%20progression.
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

################## GSEA along lineages

pt1 <- which(complete.cases(combined.all$pt1))
combined.all$pseudotime1 <- "no"
combined.all@meta.data[pt1, "pseudotime1"] <- "yes"
combined.all <- SetIdent(combined.all, value = "pseudotime1")
pt1 <- subset(combined.all,  idents = "yes")

DefaultAssay(pt1) <- "integrated"
pt1 <- RunUMAP(pt1, dims = 1:20) 
pt1 <- FindNeighbors(pt1, reduction = "pca", dims = 1:20)
pt1 <- FindClusters(pt1, resolution = c(0.5, 0.8, 1.0))
pt1 <- SetIdent(pt1, value = "integrated_snn_res.0.5")
DefaultAssay(pt1) <- "RNA"
DimPlot(pt1, label = T)
pt1 <- SetIdent(pt1, value = "main")
DimPlot(pt1, label = T)

### lineage 1
library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(hexbin)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

table(Idents(pt1), pt1$sample)

# cluster 3 Cells 
c3 <- subset(pt1, idents = "3")
Idents(c3) = c3$sample
DefaultAssay(c3) = 'RNA'
c3 = NormalizeData(c3) %>% ScaleData()

mtx = as.matrix(GetAssayData(c3,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c3$sample)
c3$sample = factor(c3$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c3$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

# load libraries
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

c3res<- pval
c3res[,"symbol"] <- rownames(c3res)

# Add Annotations
c3res$entrez = mapIds(org.Hs.eg.db,
                      keys=c3res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c3res$ensg = mapIds(org.Hs.eg.db,
                    keys=c3res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c3res$logFC
names(foldchanges) = c3res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# Modify barplot
# stat. mean = mean of the individual statistics from multiple single array based gene set tests.
#Normally, its absolute value measures the magnitude of gene-set level changes,
#and its sign indicates direction of the changes. When saaTest=gs.KSTest, stat.mean
#is always positive.

df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(-log(df$q.val[1:10]), names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")

############ cluster 5
c5 <- subset(pt1, idents = "5")
Idents(c5) = c5$sample
DefaultAssay(c5) = 'RNA'
c5 = NormalizeData(c5) %>% ScaleData()

mtx = as.matrix(GetAssayData(c5,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c5$sample)
c5$sample = factor(c5$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c5$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c5res<- pval
c5res[,"symbol"] <- rownames(c5res)

# Add Annotations
c5res$entrez = mapIds(org.Hs.eg.db,
                      keys=c5res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c5res$ensg = mapIds(org.Hs.eg.db,
                    keys=c5res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c5res$logFC
names(foldchanges) = c5res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:4], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")

############ cluster 0
c0 <- subset(pt1, idents = "0")
Idents(c0) = c0$sample
DefaultAssay(c0) = 'RNA'
c0 = NormalizeData(c0) %>% ScaleData()

mtx = as.matrix(GetAssayData(c0,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c0$sample)
c0$sample = factor(c0$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c0$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c0res<- pval
c0res[,"symbol"] <- rownames(c0res)

# Add Annotations
c0res$entrez = mapIds(org.Hs.eg.db,
                      keys=c0res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c0res$ensg = mapIds(org.Hs.eg.db,
                    keys=c0res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c0res$logFC
names(foldchanges) = c0res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:12], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")



############ cluster 6
c6 <- subset(pt1, idents = "6")
Idents(c6) = c6$sample
DefaultAssay(c6) = 'RNA'
c6 = NormalizeData(c6) %>% ScaleData()

mtx = as.matrix(GetAssayData(c6,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c6$sample)
c6$sample = factor(c6$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c6$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c6res<- pval
c6res[,"symbol"] <- rownames(c6res)

# Add Annotations
c6res$entrez = mapIds(org.Hs.eg.db,
                      keys=c6res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c6res$ensg = mapIds(org.Hs.eg.db,
                    keys=c6res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c6res$logFC
names(foldchanges) = c6res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:4], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")


############ cluster 4
c4 <- subset(pt1, idents = "4")
Idents(c4) = c4$sample
DefaultAssay(c4) = 'RNA'
c4 = NormalizeData(c4) %>% ScaleData()

mtx = as.matrix(GetAssayData(c4,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c4$sample)
c4$sample = factor(c4$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c4$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c4res<- pval
c4res[,"symbol"] <- rownames(c4res)

# Add Annotations
c4res$entrez = mapIds(org.Hs.eg.db,
                      keys=c4res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c4res$ensg = mapIds(org.Hs.eg.db,
                    keys=c4res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c4res$logFC
names(foldchanges) = c4res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:6], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")


############ cluster 1
c1 <- subset(pt1, idents = "1")
Idents(c1) = c1$sample
DefaultAssay(c1) = 'RNA'
c1 = NormalizeData(c1) %>% ScaleData()

mtx = as.matrix(GetAssayData(c1,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c1$sample)
c1$sample = factor(c1$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c1$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c1res<- pval
c1res[,"symbol"] <- rownames(c1res)

# Add Annotations
c1res$entrez = mapIds(org.Hs.eg.db,
                      keys=c1res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c1res$ensg = mapIds(org.Hs.eg.db,
                    keys=c1res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c1res$logFC
names(foldchanges) = c1res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:6], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")


############ cluster 2
c2 <- subset(pt1, idents = "2")
Idents(c2) = c2$sample
DefaultAssay(c2) = 'RNA'
c2 = NormalizeData(c2) %>% ScaleData()

mtx = as.matrix(GetAssayData(c2,slot='counts'))
filter = rowSums(mtx>0) > 20
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=c2$sample)
c2$sample = factor(c2$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + c2$sample)
colnames(design) = c('Control','Infected')
contrast.matrix = makeContrasts(Infected-Control,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)

fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

deg_genes_up = rownames(pval)[which(pval$fdr<0.005 & pval$logFC>0.5)]
deg_genes_down = rownames(pval)[which(pval$fdr<0.005 & pval$logFC<(-0.5))]

c2res<- pval
c2res[,"symbol"] <- rownames(c2res)

# Add Annotations
c2res$entrez = mapIds(org.Hs.eg.db,
                      keys=c2res$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
c2res$ensg = mapIds(org.Hs.eg.db,
                    keys=c2res$symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")


# Prepare Input for GAGE - all log fold change values and entrez IDs
foldchanges = c2res$logFC
names(foldchanges) = c2res$entrez

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 
msigres = gage(foldchanges, gsets=msigdb, same.dir=TRUE)
#lapply(msigres, head)

############# barplot
df1 <- as.data.frame(msigres$greater) %>% filter(q.val <= 0.05)
df2 <- as.data.frame(msigres$less) %>% filter(q.val <= 0.05)
df <- rbind(df1,df2)
par(mar=c(5,22,4,1)+.1)
barplot(df$stat.mean[1:6], names.arg = rownames(df), col="red", horiz=TRUE, cex.names=0.9, las =2, main = "MSigDB Biological Processes", xlab = "-log(adj.p.val")



##### subset to HSC, MPP group

early <- subset(combined.all, idents = c("HSC", "MPP1", "MPP2", "MPP3", "MPP4"))
DefaultAssay(early) <- "integrated"
early <- RunPCA(early, npcs = 30, verbose = FALSE)
early <- RunUMAP(early, dims = 1:20) # reduction = "pca",
early <- FindNeighbors(early, reduction = "pca", dims = 1:20)
early <- FindClusters(early, resolution = 0.7)
DimPlot(early, label=TRUE, label.size = 8)
head(early@meta.data)
early<- SetIdent(early, value = early$main)
DefaultAssay(early) <- "RNA"

DefaultAssay(early) <- "RNA"
early <- FindVariableFeatures(early, nfeatures = 500)
var.features <- early@assays$RNA@var.features

sce <- as.SingleCellExperiment(early, assay = "RNA")
sce <- sce[var.features,]

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


sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$integrated_snn_res.0.7, 
                 start.clus = '0', approx_points = 300)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$integrated_snn_res.2, dist.method = "mnn",
                 start.clus = '12', approx_points = 300)
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

cell_colors_clust <- cell_pal(early$main, hue_pal())
# inf_cell_colors_clust <- cell_pal(infected$main, hue_pal())
# mock_cell_colors_clust <- cell_pal(mock$main, hue_pal())
cols <- unique(cell_colors_clust)
cols <- c("#F8766D", "#EC8239", "#DB8E00", "#C79800", "#AEA200", "#8FAA00", "#64B200", "#00B81B",
          "#00BD5C", "#00C085", "#00C1A7", "#00BFC4", "#00BADE", "#00B2F3", "#00A6FF", "#7C96FF", 
          "#B385FF", "#D874FD", "#EF67EB", "#FD61D3", "#FF63B6", "#FF6B94")
vars <- as.data.frame(table(early$main))
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
par(mar = c(3,4,3,10))
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 20, bty = "n", cex = 1,  inset = c(-0.23, 0), xpd = T)
lines(SlingshotDataSet(crv1), lwd = 3,  col = 'black')

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 1],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 1])



# FINAL VISUALIZATIONS
# Dimensional Reduction UMAP
combined.all <- SetIdent(combined.all, value = combined.all$main)
colors = c('#87CEFA', '#20B2AA', '#00FFFF', '#DB7093',  '#FFB6C1', '#1E90FF',  '#9ACD32',  '#A0522D',
           '#CC9966',  '#0000CD', '#98FB98', '#DC143C', '#E0E0E0',  '#FF69B4', '#FF7F50', '#7B68EE',
           '#6B8E23',  '#BC8F8F', '#8A2BE2',  '#FFFF00', '#FFA500',  '#BA55D3', '#B0C4DE')
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "three", cols = colors) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "sample", cols = colors) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))

p <- Seurat::DimPlot(combined.all, label = F, pt.size = 2.4, label.size = 4.5, repel = T) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", size = 4.5)

# Feature Plots
DefaultAssay(combined.all) <- "RNA"
FeaturePlot(combined.all, features = c("CD34", "CD38", "AVP", "MPO"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)
FeaturePlot(combined.all, features = c("DNTT", "VPREB3","RAG1", "MKI67"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)
FeaturePlot(combined.all, features = c("IRF8", "BATF3", "LILRA4", "S100A8"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)
FeaturePlot(combined.all, features = c("ELANE", "CEBPD", "CYGB", "HVCN1"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)
FeaturePlot(combined.all, features = c("APOC1", "GATA1","GATA2", "HVCN1"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)
FeaturePlot(combined.all, features = c("VCAN", "ACY3", "FCER1A", "HOPX"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = T, repel=T)


# Dot Plot
Seurat::DotPlot(combined.all, assay = "RNA", features = c('CD38','CD34','AVP', 'PROM1', 'MPO','CEBPD', 'HOPX',
                                                          'DNTT', 'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKI67', 
                                                          'RAG1','HVCN1', 'ACY3', 'IRF8','BATF3','LILRA4',
                                                          'FCER1A', 'ELANE', 'S100A8', 'S100A9','VCAN', 
                                                          'GATA1', 'GATA2','KLF1','APOC1','TPSB2','EPX'), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
# Volcano Plot
Idents(combined.all) <- "main"
table(combined.all$main)
MPP <- WhichCells(combined.all, idents = c("MPP1", "MPP2", "MPP3", "MPP4"))
dc <- WhichCells(combined.all, idents = c("CDP", "pre-DCs", "pre-DC Cycle"))
mono <- WhichCells(combined.all, idents = c("MP", "pre-Monocytes", 'Mye-DC'))
clpb <- WhichCells(combined.all, idents = c("LMPP", "CLP.1", "CLP.2", "CLP cycle", "pre-Pro-B cycle", "pre-Pro-B.1", "pre-Pro-B.2", "Pro-B"))

combined.all$lin_main <- combined.all$main
combined.all@meta.data[MPP, "lin_main"] <- "MPP"
combined.all@meta.data[dc, "lin_main"] <- "DC-Lineage"
combined.all@meta.data[mono, "lin_main"] <- "Monocyte-Lineage"
combined.all@meta.data[clpb, "lin_main"] <- "B-Cell_Lineage"
Idents(combined.all) <- "lin_main"
combined.all$celltype.stim.combined <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim.combined"
table(combined.all@meta.data$celltype.stim.combined)

combined.all <- SetIdent(combined.all, value = "three")

HSC_de_genes <- FindMarkers(combined.all, ident.1 = "HSC_Infected", ident.2 = "HSC_Mock")
MEMP_de_genes <- FindMarkers(combined.all, ident.1 = "MEMP_Infected", ident.2 = "MEMP_Mock")
write.csv(HSC_de_genes, "./HSC_DE_GENES.csv")
write.csv(MEMP_de_genes, "./MEMP_DE_GENES.csv")
MPP_de_genes <- FindMarkers(combined.all, ident.1 = "MPP_Infected", ident.2 = "MPP_Mock")
write.csv(MPP_de_genes, "./MPP_DE_GENES.csv")
DClin_de_genes <- FindMarkers(combined.all, ident.1 = "DC-Lineage_Infected", ident.2 = "DC-Lineage_Mock")
write.csv(DClin_de_genes, "./DClin_DE_GENES.csv")
Monolin_de_genes <- FindMarkers(combined.all, ident.1 = "Monocyte-Lineage_Infected", ident.2 = "Monocyte-Lineage_Mock")
write.csv(Monolin_de_genes, "./Monolin_DE_GENES.csv")
CLPBlin_de_genes <- FindMarkers(combined.all, ident.1 = "B-Cell_Lineage_Infected", ident.2 = "B-Cell_Lineage_Mock")
write.csv(CLPBlin_de_genes, "./CLPBlin_DE_GENES.csv")
GMP_de_genes <- FindMarkers(combined.all, ident.1 = "GMP_Infected", ident.2 = "GMP_Mock")
write.csv(GMP_de_genes, "./GMP_DE_GENES.csv")

#### many mm10 genes found upregulated --------
sum(grepl("mm10",rownames(combined.all))) #97
mouse_genes <- rownames(combined.all)[grepl("mm10-",rownames(combined.all))]
DefaultAssay(combined.all) <- "RNA"
mouse_cells <- FetchData(combined.all, vars = mouse_genes, slot = "counts")
mouse_cells$sums <- rowSums(mouse_cells)

mm10 <- GetAssayData(object = combined.all, 
                     assay = "RNA", slot = "counts")[]
mouse <- names(which(mm10>0))
### -------------------------------------------

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

clpb = subset(combined.all,idents='B-Cell_Lineage')
Idents(clpb) = clpb$sample
DefaultAssay(clpb) = 'RNA'
clpb = NormalizeData(clpb) %>% ScaleData()
mtx = as.matrix(GetAssayData(clpb,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=clpb$sample)
clpb$sample = factor(clpb$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + clpb$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./clpb_pval.csv")

MPP = subset(combined.all,idents='MPP')
Idents(MPP) = MPP$sample
DefaultAssay(MPP) = 'RNA'
MPP = NormalizeData(MPP) %>% ScaleData()
mtx = as.matrix(GetAssayData(MPP,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=MPP$sample)
MPP$sample = factor(MPP$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + MPP$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./MPP_pval.csv")

GMP = subset(combined.all,idents='GMP')
Idents(GMP) = GMP$sample
DefaultAssay(GMP) = 'RNA'
GMP = NormalizeData(GMP) %>% ScaleData()
mtx = as.matrix(GetAssayData(GMP,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=GMP$sample)
GMP$sample = factor(GMP$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + GMP$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./GMP_pval.csv")

dclin = subset(combined.all,idents='DC-Lineage')
Idents(dclin) = dclin$sample
DefaultAssay(dclin) = 'RNA'
dclin = NormalizeData(dclin) %>% ScaleData()
mtx = as.matrix(GetAssayData(dclin,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=dclin$sample)
dclin$sample = factor(dclin$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + dclin$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./dclin_pval.csv")


Monolin = subset(combined.all,idents='Monocyte-Lineage')
Idents(Monolin) = Monolin$sample
DefaultAssay(Monolin) = 'RNA'
Monolin = NormalizeData(Monolin) %>% ScaleData()
mtx = as.matrix(GetAssayData(Monolin,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=Monolin$sample)
Monolin$sample = factor(Monolin$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + Monolin$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./Monolin_pval.csv")

MEMP = subset(combined.all,idents='MEMP')
Idents(MEMP) = MEMP$sample
DefaultAssay(MEMP) = 'RNA'
MEMP = NormalizeData(MEMP) %>% ScaleData()
mtx = as.matrix(GetAssayData(MEMP,slot='counts'))
filter = rowSums(mtx>0) > 10
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=MEMP$sample)
MEMP$sample = factor(MEMP$sample,levels=c('Mock','Infected'))
design = model.matrix(~0 + MEMP$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM')
dge = estimateDisp(dge,design=design)
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)
write.csv(pval, "./MEMP_pval.csv")


SCpubr::do_VolcanoPlot(sample = combined.all,
                       de_genes = de_genes,
                       n_genes = 15,
                       pval_cutoff = 0.05,
                       FC_cutoff = 1,
                       colors.use = c("red"),
                       font.size=12)

### PROPORTIONS
prop = prop.table(table(Idents(combined.all),combined.all$sample),1) %>% data.frame()
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Control',],-prop)$cluster)
prop$cluster = factor(prop$cluster,levels=ranked)
counts = table(Idents(combined.all)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster,levels=ranked)

prop_plot = ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='Infected',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='Control',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,0.6))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('red3','royalblue4'),labels=c('Infected','Mock'))+
  labs(y='Proportion (%)')+
  theme(panel.background=element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=16,hjust=0.5,color='black'),
        axis.text.x=element_text(size=18,hjust=1,color='black',angle=45,vjust=1, face = "bold"),
        legend.position='top',
        legend.text=element_text(size=18,color='black', face = "bold"),
        legend.title=element_blank(),
        axis.title.y=element_text(size=16,hjust=0.2,color='black',margin=margin(r=10)),
        axis.title.x=element_blank())

### violin plot
colors = c('#E0E0E0', '#CC9966', '#DC143C', '#FF7F50', '#FFA500', '#FFFF00', '#9ACD32', '#6B8E23',
           '#98FB98', '#20B2AA', '#00FFFF','#1E90FF', '#87CEFA', '#0000CD', '#8A2BE2', '#7B68EE', '#BA55D3',
           '#DB7093', '#FF69B4', '#FFB6C1', '#FAEBD7', '#A0522D', '#BC8F8F', '#B0C4DE')


p = VlnPlot(combined.all,c('CD34','AVP', 'PROM1', 'MPO','CEBPD', 'DNTT',
                           'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKI67', 
                           'VPREB3','RAG1', 'ACY3', 'IRF8','BATF3','LILRA4', 'FCER1A','ELANE', 'S100A8', 'S100A9','VCAN',
                           'GATA1', 'GATA2','APOC1','EPX'),
            pt.size=0,stack=T,flip=T,fill.by='ident',assay='RNA')+
  NoLegend()+
  scale_x_discrete(position='top')+
  scale_fill_manual(values=colors)+
  theme(panel.background=element_rect(color='black',size=0.8),
        panel.spacing=unit(0,'lines'),
        strip.text.y.right=element_text(size=17,hjust=0),
        axis.text.x=element_text(size=16,angle=45,hjust=0, face= "bold"),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=element_line(size=1.2),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title=element_blank())+
  coord_cartesian(clip='off')
