# scRNA Analysis Spleen (Spleen)

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Spleen_sara/")
# saveRDS(combined.all, "./combined_all.rds")
combined.all <- readRDS("./combined_all.rds")


# Read in Cell Calls
Spleen_Infected <- read.csv("/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_Spleen_sara/spleen_infected_mixREF/gem_classification.csv")
Spleen_Infected_mouse <- subset(Spleen_Infected, Spleen_Infected$call=="mm10")

Spleen_NOTinfected <- read.csv("/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_Spleen_sara/spleen_NI_mixREF/gem_classification.csv")
Spleen_NOTinfected_mouse <- subset(Spleen_NOTinfected, Spleen_NOTinfected$call=="mm10")

# read in sample
Spleen_HIV = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_Spleen_sara/spleen_infected_mixREF/filtered_feature_bc_matrix/")
dim(Spleen_HIV)
#68886  #4062

# remove mouse cells (Spleen_HIV)  
Spleen_HIV <- Spleen_HIV[,!colnames(Spleen_HIV)%in%Spleen_Infected_mouse$barcode]
dim(Spleen_HIV)
#68886   #3987

# Create Sample1 Seurat Object (Infected)
infected_Spleen <- CreateSeuratObject(counts = Spleen_HIV, project = "Spleen_infected", min.cells = 3, min.features = 200)
dim(infected_Spleen)
#24881 3967
head(infected_Spleen)

# read sample 2
Spleen_NI = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_Spleen_sara/spleen_NI_mixREF/filtered_feature_bc_matrix/")
dim(Spleen_NI)
#68886  8757

# remove mouse cells (Spleen_NI) 
Spleen_NI <- Spleen_NI[,!colnames(Spleen_NI)%in%Spleen_NOTinfected_mouse$barcode]
dim(Spleen_NI)
#68886 8564

# create sample2 suerat object
NotInf_Spleen <- CreateSeuratObject(counts = Spleen_NI, project = "NotInf_Spleen", min.cells = 3, min.features = 200)
dim(NotInf_Spleen)
# 28090 8239
head(NotInf_Spleen)

######################################################################################
# Identify Condition
infected_Spleen[["sample"]] <- "Infected"
NotInf_Spleen[["sample"]] <- "Not"

# Merge Samples
merged_seurat <- merge(x=infected_Spleen, y=NotInf_Spleen, add.cell.id = c("Infected", "Not"))
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

#save(merged_seurat, file="./merged_Spleen_seurat.RData")
#load("./merged_Spleen_seurat.RData")
max(merged_seurat@meta.data$percent.mt)
#84.3038

# PLot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Filter based on quality plots 
Filter_Spleen <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt > -Inf & percent.mt <10)
VlnPlot(Filter_Spleen, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filter_Spleen, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Sums post Filtering
sum(Filter_Spleen$sample=="Infected")
# 3776
sum(Filter_Spleen$sample=="Not")
# 7569

# Fix Gene Names
Filter_Spleen@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",Filter_Spleen@assays$RNA@counts@Dimnames[[1]])
Filter_Spleen@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filter_Spleen@assays$RNA@counts@Dimnames[[1]])

split_seurat <- SplitObject(Filter_Spleen, split.by = "sample")
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


################### Find Doublets in Infected Spleen
# (pk identification) 
sweep.res.list_1 <- paramSweep_v3(split_seurat[[1]], PCs= 1:10, sct=T)
sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
bcmvn_1 <- find.pK(sweep.stats_1)
#pK = 0.16

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

load("./firstsweep_Spleen.RData")
# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[1]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.029*nrow(split_seurat[[1]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.16, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.16, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.16_110", sct=T)


# PLot Doublets & Singlets
DimPlot(split_seurat[[1]], reduction = "umap", group.by = c("DF.classifications_0.25_0.16_110",
                                                            "DF.classifications_0.25_0.16_96"),label = TRUE,repel = TRUE)

################### Doublet Finder for Non-Infected Spleen
# (pk identification) 
sweep.res.list_2 <- paramSweep_v3(split_seurat[[2]], PCs= 1:10, sct=T)
sweep.stats_2 <- summarizeSweep(sweep.res.list_2, GT=F)
bcmvn_2 <- find.pK(sweep.stats_2)
# pK = 0.1
######## label Pk chart
pK=as.numeric(as.character(bcmvn_2$pK))
BCmetric=bcmvn_2$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
#######scRNA_sara/scRNA_Spleen_sara/Spleen_NI_mixREF

# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[2]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.051*nrow(split_seurat[[2]]@meta.data)) # Adjust according to Doublet Rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder ---- edit pK
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.1, nExp = nExp_poi.adj, reuse.pANN ="pANN_0.25_0.1_386", sct=T)


# REMOVE DOUBLETS #
split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.16_110!="Doublet") 
dim(split_seurat[[1]])
#3666
split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset =DF.classifications_0.25_0.1_386 !="Doublet") 
dim(split_seurat[[2]])
#7183

########################################################## Downsample (3600 cells/sample)
set.seed(111)
split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 3600, replace=F)]
split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 3600, replace=F)]

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


######## Save files
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(combined.all, assay='RNA', slot='counts')
writeMM(counts_matrix, file='/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_Spleen_sara/counts_Spleen_organoid.mtx')
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



########## Make tSNE and UMAP plots
DefaultAssay(combined.all) <- "integrated"
DimPlot(combined.all, reduction = "tsne", split.by = "sample",label = TRUE,repel = TRUE)

DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", label = TRUE,repel = TRUE)

### For 2 samples run separately
# Spleen infected
DimPlot(split_seurat[[1]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[1]], reduction = "umap", label = TRUE,repel = TRUE)
# Spleen not infected
DimPlot(split_seurat[[2]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[2]], reduction = "umap", label = TRUE,repel = TRUE)



DefaultAssay(combined.all) <- "RNA"

# Manually Annotate Clusters
combined_markers <- FindAllMarkers(object = combined.all, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25) 

comb_markers <- combined_markers[ , c(6, 7, 1:5)]
comb_markers <- comb_markers %>%
  dplyr::arrange(cluster, p_val_adj)

# create columns for annotations
# combined.all@meta.data$basic <- "NA"
# combined.all@meta.data$zoom <- "NA"


## October 20 2022
saveRDS(combined.all, "./combined_all_corrected.rds")
save.image("./CA_101622.RData")
load("./CA_101622.RData")
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.7")

FeaturePlot(combined.all, label=T, cols=c("yellow", "blue"), features = c("LYZ", "AIF1", "COTL1", "CD24", "RAG1", 'MS4A1'))
FeaturePlot(combined.all, label=T, cols=c("yellow", "blue"), features = c("IRF8", "IRF4", "IL23R", "IL1R1", "IL7R", "IL4I1"))
FeaturePlot(combined.all, label=T, cols=c("yellow", "blue"), features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD4", "CD40LG"))
combined.all@meta.data$main <- "NA"
bcells <- WhichCells(combined.all, idents = c("0", "1","11", "13", "3", "7", "2", "12", "14"))
tcells <- WhichCells(combined.all, idents =c("4", "5", '9'))
dc <- WhichCells(combined.all, idents =c("17"))
mait <- WhichCells(combined.all, idents =c("15"))
myeloid <- WhichCells(combined.all, idents =c("10", '16'))
nkcells <- WhichCells(combined.all, idents =c("6", "8"))

combined.all <- SetIdent(combined.all, value = "more")
DimPlot(combined.all, label = T)
combined.all@meta.data[bcells, "main"] <- "B Cells"
combined.all@meta.data[tcells, "main"] <- "T Cells"
combined.all@meta.data[nkcells, "main"] <- "NK Cells"
combined.all@meta.data[dc, "more"] <- "pDC Cells"
combined.all@meta.data[myeloid, "main"] <- "Mye Cells"
combined.all@meta.data[mait, "more"] <- "MAIT Cells"

combined.all <- SetIdent(combined.all, value = "more")
combined.all <- SetIdent(combined.all, value = "seurat_clusters")
DimPlot(combined.all, label = T)

mait_marks <- FindMarkers(combined.all, ident.1 = "MAIT Cells", only.pos = T)
c16_marks <- FindMarkers(combined.all, ident.1 = "16", only.pos = T)
# Normalize RNA Counts for Analysis because SCT was performed on a sample wise basis and would create batch effects
DefaultAssay(combined.all) <- "RNA"
NormalizeData(combined.all, assay = "RNA")
ScaleData(combined.all, assay = "RNA")


# ANALYSIS
# load("./combinedALL.RData")
# saveRDS(combined.all, "./combined_all.rds")
combined.all <- readRDS("./combined_all.rds")

combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.6")
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.7")
combined.all <- SetIdent(combined.all, value = "main")
combined.all <- SetIdent(combined.all, value = "mas")

DefaultAssay(combined.all) <- "SCT"
head(combined.all@meta.data)
DimPlot(combined.all, label=TRUE, label.size = 5)
FeaturePlot(combined.all, label = T, cols=c("yellow", "blue"), features = c("FCGR3A", "NCAM1", "KLRB1", "GNLY", "NKG7", "CD14"))
combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"

cd16_DEgenes <- FindMarkers(combined.all, ident.1 = "CD16+ Mono_Infected", ident.2 = "CD16+ Mono_Not", verbose = FALSE)
head(cd16_DEgenes, n = 15)

FeaturePlot(combined.all, cols= c("yellow", "blue"), label=T, features = c("CD8A", "CD8B", "CD4", "CD40LG", "CD3D", "CD3E"))


# T CELLS
TCells <- WhichCells(combined.all, idents = c("9","5", "4"))
tcells <- subset(combined.all, cells = TCells)
DefaultAssay(tcells) <- "integrated"
tcells <- RunPCA(tcells, npcs = 30, verbose = FALSE)
tcells <- RunUMAP(tcells, dims = 1:20) # reduction = "pca",
tcells <- FindNeighbors(tcells, reduction = "pca", dims = 1:20)
tcells <- FindClusters(tcells, resolution = 0.8)
head(tcells@meta.data)
tcells<- SetIdent(tcells, value = "integrated_snn_res.0.8")
DimPlot(tcells, label=T, label.size = 6)
table(tcells$integrated_snn_res.0.8)
FeaturePlot(tcells, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD8A", "CD40LG", "CD3D", "CD3E"))
DefaultAssay(tcells) <- "SCT"
cd4t <- WhichCells(tcells, idents = c("2", "4")) # 328 Cells
cd8t <- WhichCells(tcells, idents = c("1", "0", "3", "5", "6")) # 999 Cells
cd4 <- subset(tcells, cells = cd4t)
cd8 <- subset(tcells, cells = cd8t)

# save.image("./combined_all_HB.RData")
load("./combined_all_HB.RData")

# CD4 
DefaultAssay(cd4) <- "integrated"
cd4 <- RunPCA(cd4, npcs = 30, verbose = FALSE)
cd4 <- RunUMAP(cd4, dims = 1:20) # reduction = "pca",
cd4 <- FindNeighbors(cd4, reduction = "pca", dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 1.2)
cd4<- SetIdent(cd4, value = "integrated_snn_res.1.2")
DimPlot(cd4, label=T, label.size = 6)
table(cd4$integrated_snn_res.1.2)

DefaultAssay(cd4) <- "RNA"
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "CD27", "LYAR"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("TIMP1", "LGALS1", "GZMK", "CST7", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("FOXP3", "IL2RA", "IL1R1", "CD59"))

DoHeatmap(cd4, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                               "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B"), assay = "RNA", slot="data", angle = 90) + NoLegend()
table(cd4$integrated_snn_res.1)

# cd8 
DefaultAssay(cd8) <- "integrated"
cd8 <- RunPCA(cd8, npcs = 30, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:20) # reduction = "pca",
cd8 <- FindNeighbors(cd8, reduction = "pca", dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 0.8)
cd8<- SetIdent(cd8, value = "integrated_snn_res.0.8")
DimPlot(cd8, label=T)

# B CELLS
DefaultAssay(Bcells) <- "integrated"
#bcells <-WhichCells(combined.all, idents = c("4", "1", "6", "7", "10", "0", "3", "12"))
Bcells <- subset(combined.all, cells=bcells)
Bcells <- RunPCA(Bcells, npcs = 30, verbose = FALSE)
Bcells <- RunUMAP(Bcells, dims = 1:20) # reduction = "pca",
Bcells <- FindNeighbors(Bcells, reduction = "pca", dims = 1:20)
Bcells <- FindClusters(Bcells, resolution = 0.6)
head(Bcells@meta.data)
Bcells<- SetIdent(Bcells, value = "integrated_snn_res.0.6")
DimPlot(Bcells, label=T)
DoHeatmap(Bcells, features = c("CR2", "IGHD", "IGHM", "CD24", "CD19", "CD79A", "CYGB", "RAG1", "MME", "MS4A1", 
                               "IL2RA", "BANK1", "BLK", "FCER2", "CD27", "MZB1", "JCHAIN", "PRDM1", "MKI67",
                               "CDK1", "TOP2A", "TCL1A", "SLAMF7", "CD38", "CD40", "DNTT", "CD22", "SDC1",
                               "CD84", "TNFRSF13B", "CD1C", "CD80"), assay = "RNA", slot="data", angle = 90) + NoLegend()

tranb1 <- WhichCells(Bcells, idents = c("12"))
tranb2 <- WhichCells(Bcells, idents = c("2", "5"))
tranb3 <- WhichCells(Bcells, idents = c("6"))
tranbcyc <- WhichCells(Bcells, idents = c("11"))
follb <- WhichCells(Bcells, idents=c("0", "1", "3", "4", "7", "8", "9", "10", "13"))
combined.all@meta.data[tranb1, "more"] <- "Transitional B.1"
combined.all@meta.data[tranb2, "more"] <- "Transitional B.2"
combined.all@meta.data[tranb3, "more"] <- "Transitional B.3"
combined.all@meta.data[tranbcyc, "more"] <- "Transitional Cycle"
combined.all@meta.data[follb, "more"] <- "Follicular B Cells"

twelve <- WhichCells(combined.all, idents = "12")
combined.all@meta.data[twelve, "mas"] <- "transitional B cycle"
t12 <- subset(combined.all, cells = twelve)
t12 <- RunPCA(t12, npcs = 30, verbose = FALSE)
t12 <- RunUMAP(t12, dims = 1:20) # reduction = "pca",
t12 <- FindNeighbors(t12, reduction = "pca", dims = 1:20)
t12 <- FindClusters(t12, resolution = 0.7)
t12<- SetIdent(t12, value = "integrated_snn_res.0.7")
DimPlot(t12, label=T)
DoHeatmap(t12, features = c("CR2", "IGHD", "IGHM", "CD24", "CD19", "CD79A", "CYGB", "RAG1", "MME", "MS4A1", 
                               "IL2RA", "BANK1", "BLK", "FCER2", "CD27", "MZB1", "JCHAIN", "PRDM1", "MKI67",
                               "CDK1", "TOP2A", "TCL1A", "SLAMF7", "CD38", "CD40", "DNTT", "NCAM1", "KLRC1", "CX3CR1", "IL2RB", 
                            "SELL", "CXCR3", "IL2RB", "ITGB2", 'PRF1', "CD21", "EBF1", "PAX5"), assay = "RNA", slot="data", angle = 90) + NoLegend()

trbcyc <- WhichCells(t12, idents = "0")
combined.all@meta.data[trbcyc, "mas"] <- "B-cycle"
nkcyc <- WhichCells(t12, idents = c("1"))
combined.all@meta.data[nkcyc, "mas"] <- "NK cycle"
nkb <- WhichCells(t12, idents = c("2"))
combined.all@meta.data[nkb, "mas"] <- "Plasma"


# one NA cell
na1 <- WhichCells(combined.all, idents = "NA")
combined.all@meta.data[na1, "mas"] <- "B-cycle"

# Cluster 13
DefaultAssay(t13) <- "RNA"
thirteen <- WhichCells(combined.all, idents = "13")
t13 <- subset(combined.all, cells = thirteen)
FeaturePlot(t13, label=TRUE, features = c("SELL", "KLRD1", "CD27", "KLRB1", "CD3D", "CD3E", "CD3G", "CD4", "CD40LG", "CD8A", "CD8B",
                                          "CD44", "CCR6", "CXCR6", "CD27"), cols=c("yellow", "blue"))

FeaturePlot(t13, label=TRUE, features = c("KLRB1", "IL7R", "SLC4A10", "CCR6"), cols=c("yellow", "blue"))


# NA/14
DefaultAssay(f14) <- "RNA"
fourteen <- WhichCells(combined.all, idents = "14")
f14 <- subset(combined.all, cells = fourteen)
FeaturePlot(f14, label=TRUE, features = c("CD4", "IRF8", "HLA-DRA"), cols=c("yellow", "blue"))
c14_marks <- FindMarkers(combined.all, ident.1 = 14) # pDCs
combined.all@meta.data[fourteen, "mas"] <- "pDCs"

# All cells HeatMap
DoHeatmap(combined.all, features = c("CR2", "IGHD", "IGHM", "CD24", "CD19", "CD79A", "CYGB", "RAG1", "MME", "MS4A1", 
                               "IL2RA", "BANK1", "BLK", "FCER2", "CD27", "MZB1", "JCHAIN", "MKI67",
                               "CDK1", "TOP2A", "TCL1A", "SLAMF7", "CD38", "CD40", "DNTT", "CD22", "CD27", "CD8A", "CD8B", "CD40LG",
                               "CCR7", "IL7R", "SELL", "DTHD1", "LYST", "LEF1","GNLY", "NKG7", "NCAM1", "CCR4", "CD4", "GZMH", "GZMK", "FCGR3A", "CD14", "ITGAM",
                               "KLRG1", "KLRB1", "CSF2", "CXCR3", "IRF8", "KLRC1", "IL2RB", "CX3CR1", "ITGB2", "PRF1"), assay = "RNA", slot="data", angle = 90) + NoLegend()



#ANALYSIS
genes_intersted <- readxl::read_xlsx("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara/interested_genes_scRNA.xlsx", col_names = F)
colnames(genes_intersted) <- "Genes"

##### violin plots for genes interested
combined.all <- SetIdent(combined.all, value = "more")
DefaultAssay(combined.all) <- "RNA"
p.violin <- VlnPlot(object = combined.all, 
                    features = genes_intersted$Genes[2:3],
                    split.by = "sample",
                    pt.size = 0.2,
                    combine = TRUE) + theme(legend.position = 'right')

DotPlot(combined.all, assay = "RNA",  features=  genes_intersted$Genes[1:23], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()
DotPlot(combined.all, assay = "RNA",  features=  genes_intersted$Genes[24:49], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()
###########################################

combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
table(combined.all@meta.data$celltype.stim)


cd16_DEgenes <- FindMarkers(combined.all, ident.1 = "CD16+ Mono_Infected", ident.2 = "CD16+ Mono_Not", verbose = FALSE)
cd14_DEgenes <- FindMarkers(combined.all, ident.1 = "CD14+ Mono_Infected", ident.2 = "CD14+ Mono_Not", verbose = FALSE)
cd4naive_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ Naive T_Infected", ident.2 = "CD4+ Naive T_Not", verbose = FALSE)
cd4treg_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ TREG_Infected", ident.2 = "CD4+ TREG_Not", verbose = FALSE)
follb_DEgenes <- FindMarkers(combined.all, ident.1 = "Follicular B Cell_Infected", ident.2 = "Follicular B Cell_Not", verbose = FALSE)
tranB_DEgenes <- FindMarkers(combined.all, ident.1 = "Transitional B Cell_Infected", ident.2 = "Transitional B Cell_Not", verbose = FALSE)
dimNK_DEgenes <- FindMarkers(combined.all, ident.1 = "CD56-dim NK_Infected", ident.2 = "CD56-dim NK_Not", verbose = FALSE)
brightNK_DEgenes <- FindMarkers(combined.all, ident.1 = "CD56-bright NK_Infected", ident.2 = "CD56-bright NK_Not", verbose = FALSE)


cd16 <- subset(combined.all, idents = "CD16+ Mono")
Idents(cd16) <- "sample"
avg.cd16 <- as.data.frame(log10(AverageExpression(cd16, verbose = FALSE)$RNA))
avg.cd16$gene <- rownames(avg.cd16)
genes.to.label = c("LIPC", "TCF7L1", "ABCB11", "mm10---Cox7a2")
p1 <- ggplot(avg.cd16, aes(Infected, Not)) + geom_point() + ggtitle("CD16")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

cd4 <- subset(combined.all, idents = "CD4+ Naive T")
Idents(cd4) <- "sample"
avg.cd4 <- as.data.frame(log10(AverageExpression(cd4, verbose = FALSE)$RNA))
avg.cd4$gene <- rownames(avg.cd4)
genes.to.label = c("LAMB4", "IFI27", "AL163973.2", "AF038458.3", "BRDT")
p1 <- ggplot(avg.cd4, aes(Infected, Not)) + geom_point() + ggtitle("cd4")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


#CD11C
FeaturePlot(combined.all, label = T, cols=c("yellow","blue"), features = "ITGAM")
AverageExpression(combined.all, features = "ITGAX", assays="RNA", slot = "data")
DoHeatmap(combined.all, features = c("ITGAX"), assay = "RNA", slot="data", angle = 90) + NoLegend()











################### Hongbo T-Cell Analysis ######################################
# ANALYSIS
# save.image("./combined_all_HB.RData")
load("./combined_all_HB.RData")
# saveRDS(combined.all, "./hbCombined.all.rds")

combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.6")
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.7")
combined.all <- SetIdent(combined.all, value = "main")
combined.all <- SetIdent(combined.all, value = "mas")
combined.all <- SetIdent(combined.all, value = "more")
# head(combined.all@meta.data)
na <- names(which(is.na(combined.all$more)))
DimPlot(combined.all, label=TRUE, label.size = 5)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE, label.size =4, repel = TRUE)

myeDC <- WhichCells(combined.all, idents = "Myeloid")
combined.all@meta.data[myeDC, "more"] <- "cDC"

combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all@meta.data$celltype.stim)

cd16_DEgenes <- FindMarkers(combined.all, ident.1 = "CD16+ Mono_Infected", ident.2 = "CD16+ Mono_Not", verbose = FALSE)
head(cd16_DEgenes, n = 15)


# Use log normalized expression for visualization and DESeq
combined.all <- NormalizeData(combined.all, assay = "RNA")
combined.all <- ScaleData(combined.all, assay = "RNA")

FeaturePlot(combined.all, cols= c("yellow", "blue"), features = c("CD1C", "FCGR3B", "IL7R"))


# T CELLS
TCells <- WhichCells(combined.all, idents = c("2","8"))
tcells <- subset(combined.all, cells = TCells)
DefaultAssay(tcells) <- "RNA"
tcells <- RunPCA(tcells, npcs = 30, verbose = FALSE)
tcells <- RunUMAP(tcells, dims = 1:20) # reduction = "pca",
tcells <- FindNeighbors(tcells, reduction = "pca", dims = 1:20)
tcells <- FindClusters(tcells, resolution = 0.6)
head(tcells@meta.data)
tcells<- SetIdent(tcells, value = "integrated_snn_res.0.6")
DimPlot(tcells, label=T, label.size = 6)
table(tcells$integrated_snn_res.0.6)
FeaturePlot(tcells, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD8A", "CD40LG", "CD3D", "CD3E"))

cd4t <- WhichCells(tcells, idents = c("2", "4")) # 293 Cells
cd8t <- WhichCells(tcells, idents = c("1", "0", "3", "5")) # 1198 Cells
cd4 <- subset(tcells, cells = cd4t)
cd8 <- subset(tcells, cells = cd8t)

# save.image("./combined_all_HB.RData")
load("./combined_all_HB.RData")

# CD4 
DefaultAssay(cd4) <- "integrated"
cd4 <- RunPCA(cd4, npcs = 30, verbose = FALSE)
cd4 <- RunUMAP(cd4, dims = 1:20) # reduction = "pca",
cd4 <- FindNeighbors(cd4, reduction = "pca", dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 1.7)
cd4<- SetIdent(cd4, value = "integrated_snn_res.1.7")
DimPlot(cd4, label=T, label.size = 6)
table(cd4$integrated_snn_res.1.7)

DefaultAssay(cd4) <- "RNA"
DefaultAssay(cd4) <- "SCT"
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "CD27", "LYAR"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("TIMP1", "LGALS1", "GZMK", "CST7", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("S100A4", "CCR7", "IL32", "ISG15"))


DoHeatmap(cd4, features = c("CD3D", "CD3E", "CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CCR5"), assay = "RNA", slot="data", angle = 90) + NoLegend()
table(cd4$integrated_snn_res.1)

cd4$celltype <- paste(Idents(cd4), cd4$sample, sep = "_")
DimPlot(cd4, label=T, split.by = "sample")

cd4treg <- WhichCells(cd4, idents= c("1"))
combined.all@meta.data[cd4treg, "more"] <- "CD4+ TREG"
cd4emtm <- WhichCells(cd4, idents= c("3"))
combined.all@meta.data[cd4emtm, "more"] <- "CD4+ EM/TM"
cd4n <- WhichCells(cd4, idents= c("0","4", "2", "5", "6"))
combined.all@meta.data[cd4n, "more"] <- "CD4+ "

cd4$celltype.inf <- paste(Idents(cd4), cd4$sample, sep = "_")
Idents(cd4) <- "celltype.inf"
table(cd4$celltype.inf)

# CD4+ CCR5 & TCF7 (TCF-1) expression
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("CCR5", "TCF7"))
DoHeatmap(cd4, features = c("CCR5", "TCF7"), assay = "RNA", slot="data", angle = 90) + NoLegend()

cd4marks <- FindAllMarkers(cd4)

# cd8 
DefaultAssay(cd8) <- "integrated"
cd8 <- RunPCA(cd8, npcs = 30, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:20) # reduction = "pca",
cd8 <- FindNeighbors(cd8, reduction = "pca", dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 1.2)
cd8<- SetIdent(cd8, value = "integrated_snn_res.1.2")
DimPlot(cd8, label=T, label.size = 6)
table(cd8$integrated_snn_res.1.2)

DefaultAssay(cd8) <- "RNA"
DefaultAssay(cd8) <- "SCT"
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("TRADD", "CXCR3", "GZMA", "GZMK", "CCL5", "CD27"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("FGFBP2", "GZMH", "GZMB", "GNLY"))

FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CD8A", "GZMK", "CCL5"))
DoHeatmap(cd8, features = c("CD4", "CD40LG","CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                            "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY", "NCAM1", "FCGR3A", "MKI67"), assay = "RNA", slot="data", angle = 90) + NoLegend()

naive8 <- WhichCells(cd8, idents= c("1", "2", "0", "3"))
ccl5 <- WhichCells(cd8, idents= c("4","5"))
nk_cycle <- WhichCells(cd8, idents= c("6"))
combined.all@meta.data[ccl5, "more"] <- "CD8+ CCL5+"
combined.all@meta.data[naive8, "more"] <- "CD8+ Naive/CM"
combined.all@meta.data[nk_cycle, "more"] <- "NK-Cycle"

cd86 <- subset(cd8, cells = nk_cycle)
DefaultAssay(cd86) <- "integrated"
cd86 <- RunPCA(cd86, npcs = 30, verbose = FALSE)
cd86 <- RunUMAP(cd86, dims = 1:20) # reduction = "pca",
cd86 <- FindNeighbors(cd86, reduction = "pca", dims = 1:20)
cd86 <- FindClusters(cd86, resolution = 1.0)
cd86<- SetIdent(cd86, value = "integrated_snn_res.1")
DimPlot(cd86, label = T)
DoHeatmap(cd86, features = c("CD4", "CD40LG","CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                             "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY", "NCAM1", "FCGR3A", "MKI67"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd8cyc <- WhichCells(cd86, idents= c("1"))
combined.all@meta.data[cd8cyc, "more"] <- "CD8+ CCL5+ cycle"
nkcyc <- WhichCells(cd86, idents= c("0"))
combined.all@meta.data[nkcyc, "more"] <- "NK-cycle"

cd8$celltype.inf <- paste(Idents(cd8), cd8$sample, sep = "_")
Idents(cd8) <- "celltype.inf"
table(cd8$celltype.inf)
DimPlot(cd8, label=T, split.by = "sample")

CD8_total <- WhichCells(combined.all, idents= c("CD8+ CCL5+", "CD8+ CCL5+ cycle", "CD8+ Naive/CM"))
cd8TOT <- subset(combined.all, cells = CD8_total)
DoHeatmap(cd8TOT, features = c("CD4", "CD40LG","CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                             "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY", "NCAM1", "FCGR3A", "MKI67"), assay = "RNA", slot="data", angle = 90) + NoLegend()

# NK
# nk <- WhichCells(combined.all, idents = c("CD56-dim NK", "CD56-bright NK", "NK cycle"))
NK <- subset(combined.all, cells= nkcells)
DefaultAssay(NK) <- "integrated"
NK <- RunPCA(NK, npcs = 30, verbose = FALSE)
NK <- RunUMAP(NK, dims = 1:20) # reduction = "pca",
NK <- FindNeighbors(NK, reduction = "pca", dims = 1:20)
NK <- FindClusters(NK, resolution = 1.0)
NK<- SetIdent(NK, value = "integrated_snn_res.1")
DimPlot(NK, label=T, label.size = 6)
table(NK$integrated_snn_res.0.8)
DefaultAssay(NK) <- "SCT"
FeaturePlot(NK, label=T, cols=c("yellow", "blue"), features = c("CD3D", "CD3E", "GNLY", "NKG7", "CD247", "NCAM1", "FCRG3A"))
DoHeatmap(NK, features = c("CD3D", "CD3E", "GNLY", "NKG7", "CD247", "NCAM1", "FCGR3A", "CXRCR1", "IL2RB"), assay = "RNA", slot="data", angle = 90) + NoLegend()

NK$celltype <- paste(Idents(NK), NK$sample, sep = "_")
Idents(NK) <- "celltype"
table(NK$celltype)
DimPlot(NK, label=T, split.by = "sample")

# "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"
DoHeatmap(NK, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                           "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"), assay = "RNA", slot="data", angle = 90) + NoLegend()

bright_ncam <- WhichCells(NK, idents= c("1", "4", "0", "3"))
dim_ncam <- WhichCells(NK, idents= c("2"))
combined.all@meta.data[bright_ncam, "more"] <- "CD56-bright NK"
combined.all@meta.data[dim_ncam, "more"] <- "CD56-dim NK"


# Myeloid
myeloid <- subset(combined.all, cells = myeloid)
DefaultAssay(myeloid) <- "integrated"
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:20) # reduction = "pca",
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.6)
head(myeloid@meta.data)
myeloid<- SetIdent(myeloid, value = "integrated_snn_res.0.6")
DimPlot(myeloid, label=T, label.size = 6)
table(myeloid$integrated_snn_res.0.6)
DefaultAssay(myeloid) <- "RNA"
FeaturePlot(myeloid, label=T, cols=c("yellow", "blue"), features = c("FCGR3A", "CD14", "CD1C", "IDO1", "IFI27"))
FeaturePlot(myeloid, label=T, cols=c("yellow", "blue"), features = c("FCGR3A", "CD14", "CST3", "LYZ", "HLA-DPA", "HLA-DPB1"))
myeloid$celltype <- paste(Idents(myeloid), myeloid$sample, sep = "_")
Idents(myeloid) <- "celltype"
table(myeloid$celltype)
DimPlot(myeloid, label=T, split.by = "sample")
mye2_marks <- FindMarkers(myeloid, ident.1 = "2")

DoHeatmap(myeloid, features = c("FCGR3A", "CD14", "CST3", "LYZ", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "IDO1",
                           "HLA-DRA", "CLEC10A", "CLEC9A", "S100A8", "S100A9", "CD74", "FTH1", "CD1C", "CD8A", "CD4", "CD40LG", "CD8B",
                           "COTL1", "AIF1", "IRF8", "IRF4"), assay = "RNA", slot="data", angle = 90) + NoLegend()


mye_2 <- WhichCells(myeloid, idents = "2")
mye_2 <- subset(combined.all, cells = mye_2)
DefaultAssay(mye_2) <- "integrated"
mye_2 <- RunPCA(mye_2, npcs = 30, verbose = FALSE)
mye_2 <- RunUMAP(mye_2, dims = 1:20) # reduction = "pca",
mye_2 <- FindNeighbors(mye_2, reduction = "pca", dims = 1:20)
mye_2 <- FindClusters(mye_2, resolution = 1.0)
mye_2<- SetIdent(mye_2, value = "integrated_snn_res.1")
DimPlot(mye_2, label=T, label.size = 6)
mye2_markers <- FindAllMarkers(mye_2)

myeloid$celltype <- paste(Idents(myeloid), myeloid$sample, sep = "_")
Idents(myeloid) <- "celltype"
table(myeloid$celltype)
DimPlot(myeloid, label=T, split.by = "sample")

cd14 <- WhichCells(myeloid, idents= c("1"))
cd16 <- WhichCells(myeloid, idents= c("0"))
c2 <- WhichCells(myeloid, idents= c("2"))
combined.all@meta.data[cd16, "more"] <- "CD16 MONO"
combined.all@meta.data[cd14, "more"] <- "CD14 MONO"
combined.all@meta.data[c2, "more"] <- "Myeloid"


# HONGBO
FeaturePlot(NK, cols=c("grey", "red"), label=T, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                  "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"))
DoHeatmap(NK, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"), assay = "RNA", slot="data", angle = 90) + NoLegend()

combined.all <- SetIdent(combined.all, value = "more")
CD4 <- subset(combined.all, idents = c("CD4+ ","CD4+ EM/TM", "CD4+ TREG"))
CD4 <- ScaleData(CD4, assay = "RNA")
DefaultAssay(CD4) <- "RNA"
FeaturePlot(CD4, cols= c("grey", "red"), label=F, features = c("CX3CR1", "BACH2", "TCF7", "CCR5"))
CD4 <- SetIdent(CD4, value = "sample")
DoHeatmap(CD4, features = c("CCR5", "TCF7", "CX3CR1", "BACH2"), assay = "RNA", slot="scale.data", angle = 90) + NoLegend()


CCR5 = GetAssayData(object = CD4, 
                    assay = "RNA", slot = "data")["CCR5",]
TCF7 = GetAssayData(object = CD4, 
                    assay = "RNA", slot = "data")["TCF7",]
Cd4 = GetAssayData(object = CD4, 
                   assay = "RNA", slot = "data")["CD4",]

names(which(TCF7>0 & CCR5<=0)) #281
names(which(CCR5>0)) #12
names(which(TCF7<=0 & CCR5<=0)) # 35
dp <- names(which(CCR5>0 & TCF7>0)) # 4 (ALL MOCK)
ccr5 <- names(which(CCR5>0 & TCF7<=0)) #8 (1 INF)

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
write.csv(ccr5_df, "/storage1/fs1/liang.shan/Active/scRNA_sara/CCR5_expr_spleen.csv")

DoHeatmap(CCR5, features = c("BACH2", "CX3CR1", "PDCD1", "TOX", "CXCL13", "TIGIT", 
                             "CTLA4", "TNFRSF9", "HAVCR2", "LAG3", "CST7", "GZMK", 
                             "GZMA", "NKG7", "IFNG", "PRF1", "GZMB", "GNLY", "PRDM1", 
                             "BATF", "CD28", "GZMH", "KLRK1", "KLRB1", "CTSW", "BTLA", 
                             "FOXP3", "TNFRSF4", "EEF1A1", "SELL", "CCR7", "IL6R", 
                             "IGFBP4", "IGFL2I", "IL7R"), assay = "RNA", slot="data", angle = 90) 

# NK Cells
nkcells <- WhichCells(combined.all, idents = c("CD56-dim NK", "CD56-bright NK", "NK-cycle"))
NKcells <- subset(combined.all, idents = c("CD56-dim NK", "CD56-bright NK", "NK-cycle"))
NKcells<- ScaleData(NKcells, assay = "RNA")
Idents(NKcells) <- "sample"

RidgePlot(NKcells, features = c("KIR2DL1",  "KIR2DL3"))

FeaturePlot(NKcells, cols=c("grey", "red"), label=T, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                  "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"))
DoHeatmap(NKcells, features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3"), assay = "RNA", slot="scale.data", angle = 90) + NoLegend()

cd158_df <- as.data.frame(AverageExpression(NKcells, verbose = FALSE, slot = "counts", features = c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2",
                                                                                                    "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR3DL2", "KIR3DL3") )$RNA)


save.image("./HONGBO_interest.RData")
load("/HONGBO_interest.RData")


















combined.all <- NormalizeData(combined.all, assay = "RNA")
combined.all <- ScaleData(combined.all, assay = "RNA")




# CD166
Idents(combined.all) <- "celltype.stim"
table(combined.all$celltype.stim)
Idents(combined.all) <- "more"
cd166 <- as.data.frame(AverageExpression(combined.all, features = c("ALCAM"), verbose = FALSE, slot = "data", assays = "RNA"))
write.csv(cd166, "./cd166_expression.csv")

interest <- WhichCells(combined.all, idents = c("CD14 MONO_Infected", "CD14 MONO_Not", "CD16 MONO_Infected", "CD16 MONO_Not", "CD4+ _Infected", "CD4+ _Not"))
INTEREST <- subset(combined.all, cells = interest)
inf14 <- WhichCells(combined.all, idents = c("CD14 MONO_Infected"))
mock14 <- WhichCells(combined.all, idents = c("CD14 MONO_Not"))
inf16 <- WhichCells(combined.all, idents = c("CD16 MONO_Infected"))
mock16 <- WhichCells(combined.all, idents = c("CD16 MONO_Not"))

combined.all <- SetIdent(combined.all, value = "sample")
not <- WhichCells(combined.all, idents = "Not")
combined.all@meta.data[not, "sample"] <- "Mock"
FeaturePlot(combined.all, features = "ALCAM", cols =c("grey", "red"), split.by = "sample", label = TRUE)
RidgePlot(combined.all, features = "ALCAM", idents = c("CD14 MONO_Infected", "CD14 MONO_Not", "CD16 MONO_Infected", "CD16 MONO_Not", "CD4+ _Infected", "CD4+ _Not"), slot = "data")

FetchData(combined.all, cells = inf16,  vars =c("ALCAM"), slot="data")
FetchData(combined.all, cells = inf14,  vars =c("ALCAM"), slot="counts")
mock16 <- as.data.frame(FetchData(combined.all, cells = mock16,  vars =c("ALCAM"), slot="counts"))
mock14 <- as.data.frame(FetchData(combined.all, cells = mock14,  vars =c("ALCAM"), slot="counts"))


my_levels <- c("CD14 MONO_Infected", "CD14 MONO_Mock", "CD16 MONO_Infected", "CD16 MONO_Mock", "cDC_Infected", "cDC_Mock",
             "pDC Cells_Infected", "pDC Cells_Mock", "CD4+ _Infected", "CD4+ _Mock", "CD4+ EM/TM_Infected", "CD4+ EM/TM_Mock", "CD4+ TREG_Infected", "CD4+ TREG_Mock",
             "MAIT Cells_Infected", "MAIT Cells_Mock", "CD56-bright NK_Infected", "CD56-bright NK_Not","CD56-dim NK_Infected",  "CD56-dim NK_Not", "CD8+ CCL5+ cycle_Infected",
             "CD8+ CCL5+ cycle_Not", "CD8+ CCL5+_Infected","CD8+ CCL5+_Not","CD8+ Naive/CM_Infected",  "CD8+ Naive/CM_Not", "Follicular B Cells_Infected", "Follicular B Cells_Mock",
             "NK-cycle_Infected", "NK-cycle_Mock", "Transitional B.1_Infected", "Transitional B.1_Mock",    "Transitional B.2_Infected", "Transitional B.2_Mock", "Transitional B.3_Infected",
             "Transitional B.3_Mock", "Transitional Cycle_Infected",  "Transitional Cycle_Mock")
combined.all@active.ident <- factor(x = combined.all@active.ident, levels = my_levels)
FeaturePlot(combined.all, features = "ALCAM", cols = c("grey", "red"), label = TRUE)

VlnPlot(combined.all, features = c("ALCAM"), split.by = "sample", 
        pt.size = 0, combine = FALSE, assay = "RNA", idents = )

INTEREST <- ScaleData(INTEREST, assay = "RNA")
alcam = GetAssayData(object = INTEREST, 
                      assay = "RNA", slot = "data")["ALCAM",]
AverageExpression(INTEREST, verbose = FALSE, slot = "data", assays = "RNA", features = "ALCAM")


my_levels <- c("CD14 MONO_Infected", "CD14 MONO_Not", "CD16 MONO_Infected", "CD16 MONO_Not", "CD4+ _Infected", "CD4+ _Not")
# Re-level object@ident
INTEREST@active.ident <- factor(x = INTEREST@active.ident, levels = my_levels)

RidgePlot(INTEREST, features = "ALCAM", slot = "counts")

Idents(INTEREST) <- "more"
VlnPlot(INTEREST, features = c("ALCAM"), split.by = "sample", 
        pt.size = 0, combine = FALSE, assay = "RNA", slot = "counts")




# SCMAP
library(SingleCellExperiment)
library(scmap)

#convert to SCE
DefaultAssay(combined.all) <- "RNA"
combined.all_Diet <- DietSeurat(combined.all, graphs = "umap")
sce_combined.all <- as.SingleCellExperiment(combined.all_Diet)
rowData(sce_combined.all)$feature_symbol <- rownames(sce_combined.all)
sce <- selectFeatures(sce_combined.all, suppress_plot = FALSE)



########################################################### RITU NK CELL ANALYSIS
combined.all <- SetIdent(combined.all, value = "more")
DimPlot(combined.all, label = T)

FeaturePlot(combined.all, features = c("GNLY", "PRF1", "TYROBP", "TRDC"), label =T)

NK <- subset(combined.all, idents = c("NK-cycle", "CD56-dim NK", "CD56-bright NK"))
FeaturePlot(NK, features = c("GNLY", "PRF1", "TYROBP", "TRDC"))

DefaultAssay(NK) <- "integrated"
NK <- RunPCA(NK, npcs = 30, verbose = FALSE)
NK <- RunUMAP(NK, dims = 1:20) 
NK <- FindNeighbors(NK, reduction = "pca", dims = 1:20)
NK <- FindClusters(NK, resolution = 0.4)
NK <- SetIdent(NK, value = "integrated_snn_res.0.4")
DimPlot(NK, label=TRUE, split.by = "sample")


NK <- NormalizeData(NK, assay = "RNA")
NK <- ScaleData(NK, assay = "RNA")

NK <- SetIdent(NK, value = "sample")
DoHeatmap(NK, features = c("CD52", "CD34", "NCAM1", "FCGR3A", "ITGAM", "B3GAT1", "KLRD1", "KLRC2", "IL2RB", "CD27", "KLRC1", "KLRK1", "KIR2DL2", "KIR2DL3", "KIR2DL4", 
                           "KIR3DL1", "KIR3DL2", "CCR7", "CXCR1", "CX3CR1", "CD7", "LILRB1", "CD69"), assay = "RNA", slot="scale.data", angle = 90) 
nklung_df_stim <- as.data.frame(AverageExpression(NK, verbose = FALSE, slot = "counts", features = c("CD52", "CD34", "NCAM1", "FCGR3A", "ITGAM", "B3GAT1", "KLRD1", "KLRC2", "IL2RB", "CD27") )$RNA)


genes <- read.csv("/storage1/fs1/liang.shan/Active/RITU/HeatMapGenes.csv")
DoHeatmap(NK, features = genes[,1], assay = "RNA", slot="scale.data", angle = 90)
DoHeatmap(NK, features = genes[1:50,2], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[51:99,2], assay = "RNA", slot="scale.data", angle = 90)
DoHeatmap(NK, features = genes[1:55,3], assay = "RNA", slot="scale.data", angle = 90)
DoHeatmap(NK, features = genes[56:111,3], assay = "RNA", slot="scale.data", angle = 90)
DoHeatmap(NK, features = genes[,4], assay = "RNA", slot="scale.data", angle = 90) 
DoHeatmap(NK, features = genes[,5], assay = "RNA", slot="scale.data", angle = 90) 



write.csv(nklung_df, "./NKlung.csv")

NK <- SetIdent(NK, value = "integrated_snn_res.0.4")
NK <- SetIdent(NK, value = "sample")
NKmarks <- FindAllMarkers(NK, assay = "RNA", slot = "data")

nk1 <- WhichCells(NK, idents = c("1"))
NK@meta.data[nk1, "more"] <- "NK1"
nk2 <- WhichCells(NK, idents = c("2"))
NK@meta.data[nk2, "more"] <- "NK2"
nk3 <- WhichCells(NK, idents = c("3"))
NK@meta.data[nk3, "more"] <- "NK3"
nk4 <- WhichCells(NK, idents = c("4"))
NK@meta.data[nk4, "more"] <- "NK4"
nk0 <- WhichCells(NK, idents = c("0"))
NK@meta.data[nk0, "more"] <- "NK0"
NK <- SetIdent(NK, value = "more")


NK$nk.infection <- paste(Idents(NK), NK$sample, sep = "_")
Idents(NK) <- "nk.infection"
table(NK@meta.data$nk.infection)

my_levels <- c("0_Infected", "0_Not",  "1_Infected", "1_Not", "2_Infected", "2_Not", "3_Infected", "3_Not", "4_Infected", "4_Not")
# Re-level object@ident
NK@active.ident <- factor(x = NK@active.ident, levels = my_levels)

NK0 <- subset(NK, cells = nk0)
NK0 <- ScaleData(NK0, assay = "RNA")
NK0marks <- FindAllMarkers(NK0, assay = "RNA", slot = "data")
write.csv(NK0marks, "/storage1/fs1/liang.shan/Active/RITU/spleen0marks.csv")

NK1 <- subset(NK, cells = nk1)
NK1 <- ScaleData(NK1, assay = "RNA")
NK1marks <- FindAllMarkers(NK1, assay = "RNA", slot = "data")
write.csv(NK1marks, "/storage1/fs1/liang.shan/Active/RITU/spleen1marks.csv")


NK2 <- subset(NK, cells = nk2)
NK2 <- ScaleData(NK2, assay = "RNA")
NK2marks <- FindAllMarkers(NK2, assay = "RNA", slot = "data")
write.csv(NK2marks, "/storage1/fs1/liang.shan/Active/RITU/spleen2marks.csv")


NK3 <- subset(NK, cells = nk3)
NK3 <- ScaleData(NK3, assay = "RNA")
NK3marks <- FindAllMarkers(NK3, assay = "RNA", slot = "data")
write.csv(NK3marks, "/storage1/fs1/liang.shan/Active/RITU/spleen3marks.csv")


NK4 <- subset(NK, cells = nk4)
NK4 <- ScaleData(NK4, assay = "RNA")
NK4marks <- FindAllMarkers(NK4, assay = "RNA", slot = "data")
write.csv(NK4marks, "/storage1/fs1/liang.shan/Active/RITU/spleen4marks.csv")


NKcell_spleen <- as.data.frame(AverageExpression(NK, verbose = FALSE, slot = "data")$RNA)
write.csv(NKcell_spleen, "/storage1/fs1/liang.shan/Active/NKspleen_expression.csv")


spl0 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleen0marks.csv", header = T)
spl0genes <- spl0$gene
spl0UP <- spl0genes[1:200]
spl0down <- spl0genes[207:406]
up0 <- as.data.frame(AverageExpression(NK0, verbose = FALSE, slot = "data", features = spl0UP)$RNA)
write.csv(up0, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen0UP.csv")
down0 <- as.data.frame(AverageExpression(NK0, verbose = FALSE, slot = "data", features = spl0down)$RNA)
write.csv(down0, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen0Down.csv")

spl1 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleen1marks.csv", header = T)
spl1genes <- spl1$gene
spl1UP <- spl1genes[1:200]
spl1down <- spl1genes[299:498]
up1 <- as.data.frame(AverageExpression(NK1, verbose = FALSE, slot = "data", features = spl1UP)$RNA)
write.csv(up1, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen1UP.csv")
down1 <- as.data.frame(AverageExpression(NK1, verbose = FALSE, slot = "data", features = spl1down)$RNA)
write.csv(down1, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen1Down.csv")

spl2 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleen2marks.csv", header = T)
spl2genes <- spl2$gene
spl2UP <- spl2genes[1:200]
spl2down <- spl2genes[203:402]
up2 <- as.data.frame(AverageExpression(NK2, verbose = FALSE, slot = "data", features = spl2UP)$RNA)
write.csv(up2, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen2UP.csv")
down2 <- as.data.frame(AverageExpression(NK2, verbose = FALSE, slot = "data", features = spl2down)$RNA)
write.csv(down2, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen2Down.csv")

spl3 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleen3marks.csv", header = T)
spl3genes <- spl3$gene
spl3UP <- spl3genes[1:130]
spl3down <- spl3genes[158:287]
up3 <- as.data.frame(AverageExpression(NK3, verbose = FALSE, slot = "data", features = spl3UP)$RNA)
write.csv(up3, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen3UP.csv")
down3 <- as.data.frame(AverageExpression(NK3, verbose = FALSE, slot = "data", features = spl3down)$RNA)
write.csv(down3, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen3Down.csv")

spl4 <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleen4marks.csv", header = T)
spl4genes <- spl4$gene
spl4UP <- spl4genes[1:65]
spl4down <- spl4genes[76:129]
up4<- as.data.frame(AverageExpression(NK4, verbose = FALSE, slot = "data", features = spl4UP)$RNA)
write.csv(up4, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen4UP.csv")
down4 <- as.data.frame(AverageExpression(NK4, verbose = FALSE, slot = "data", features = spl4down)$RNA)
write.csv(down4, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleen4Down.csv")


######################################################## RITU CD8+ T-Cell Analysis
combined.all <- SetIdent(combined.all, value = "more")
DimPlot(combined.all, label = T)

CD8 <- subset(combined.all, idents = c("CD8+ Naive/CM", "CD8+ CCL5+", "CD8+ CCL5+ cycle" ))

FeaturePlot(CD8, features = c("CD8A", "CD3G", "CD3E", "CD3D", "CD8B"))

DefaultAssay(CD8) <- "integrated"
CD8 <- RunPCA(CD8, npcs = 30, verbose = FALSE)
CD8 <- RunUMAP(CD8, dims = 1:20) 
CD8 <- FindNeighbors(CD8, reduction = "pca", dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 0.4)
CD8 <- SetIdent(CD8, value = "integrated_snn_res.0.4")
DimPlot(CD8, label=TRUE, split.by = "sample", label.size = 5)

CD8 <- NormalizeData(CD8, assay = "RNA")
CD8 <- ScaleData(CD8, assay = "RNA")

CD8 <- SetIdent(CD8, value = "sample")
DoHeatmap(CD8, features = c("CD3D", "CD3G", "CD3E", "CD8A", "CD8B", "SELL", "CCR7",
                            "TBX21", "IFNG", "CXCR3", "GZMB", "GZMH", "GNLY", "PRF1", "GZMK", "CRTAM", "CD4", "MKI67", "CD40LG"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd8expr <- as.data.frame(AverageExpression(CD8, verbose = FALSE, slot = "counts", features = c("CD3D", "CD3G", "CD3E", "CD8A", "CD8B", "SELL", "CCR7",
                                                                                               "TBX21", "IFNG", "CXCR3", "GZMB, GZMH", "GNLY", "PRF1", "GZMK, CRTAM") )$RNA)

naive <- WhichCells(CD8, idents = c(0,1))
ememra <- WhichCells(CD8, idents = c(2))
trmem <- WhichCells(CD8, idents = c(3))
CD8@meta.data[naive, "cd8_more"] <- "CD8+ Naive"
CD8@meta.data[ememra, "cd8_more"] <- "CD8+ Tem/emra"
CD8@meta.data[trmem, "cd8_more"] <- "CD8+ Trm/em"
combined.all@meta.data[naive, "more"] <- "CD8+ Trm/em"
combined.all@meta.data[naive, "more"] <- "CD8+ Naive/CM"
combined.all@meta.data[ememra, "more"] <- "CD8+ Tem/emra"

Idents(CD8) <- "cd8_more"
CD8$cd8.infection <- paste(Idents(CD8), CD8$sample, sep = "_")
Idents(CD8) <- "cd8.infection"
table(CD8@meta.data$cd8.infection)

cd8naive <- subset(CD8, cells = naive)
cd8naive <- ScaleData(cd8naive, assay = "RNA")
cd8naive_marks <- FindAllMarkers(cd8naive, assay = "RNA", slot = "data")
write.csv(cd8naive_marks, "/storage1/fs1/liang.shan/Active/RITU/spleenCD8naive_marks.csv" )

cd8em <- subset(CD8, cells = ememra)
cd8em <- ScaleData(cd8em, assay = "RNA")
cd8em_marks <- FindAllMarkers(cd8em, assay = "RNA", slot = "data")
write.csv(cd8em_marks, "/storage1/fs1/liang.shan/Active/RITU/spleenCD8ememra_marks.csv")

cd8trm <- subset(CD8, cells = trmem)
cd8trm <- ScaleData(cd8trm, assay = "RNA")
cd8trm_marks <- FindAllMarkers(cd8trm, assay = "RNA", slot = "data")
write.csv(cd8trm_marks, "/storage1/fs1/liang.shan/Active/RITU/spleenCD8trmem_marks.csv")



spleenN <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleenCD8naive_marks.csv", header = T)
spleenNgenes <- spleenN$gene
spleenNUP <- spleenNgenes[1:200]
spleenNdown <- spleenNgenes[226:349]
upN <- as.data.frame(AverageExpression(cd8naive, verbose = FALSE, slot = "data", features = spleenNUP)$RNA)
write.csv(upN, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenNUP.csv")
downN <- as.data.frame(AverageExpression(cd8naive, verbose = FALSE, slot = "data", features = spleenNdown)$RNA)
write.csv(downN, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenNDown.csv")

spleenEM <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleenCD8ememra_marks.csv", header = T)
spleenemgenes <- spleenEM$gene
spleenEMUP <- spleenemgenes[1:240]
spleenEMdown <- spleenemgenes[243:384]
upEM <- as.data.frame(AverageExpression(cd8em, verbose = FALSE, slot = "data", features = spleenEMUP)$RNA)
write.csv(upEM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenEMUP.csv")
downEM <- as.data.frame(AverageExpression(cd8em, verbose = FALSE, slot = "data", features = spleenEMdown)$RNA)
write.csv(downEM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenEMDown.csv")

spleenRM <- read.csv("/storage1/fs1/liang.shan/Active/RITU/spleenCD8trmem_marks.csv", header = T)
spleenRMgenes <- spleenRM$gene
spleenRMUP <- spleenRMgenes[1:80]
spleenRMdown <- spleenRMgenes[81:155]
upRM <- as.data.frame(AverageExpression(cd8trm, verbose = FALSE, slot = "data", features = spleenRMUP)$RNA)
write.csv(upRM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenRMUP.csv")
downRM <- as.data.frame(AverageExpression(cd8trm, verbose = FALSE, slot = "data", features = spleenRMdown)$RNA)
write.csv(downRM, "/storage1/fs1/liang.shan/Active/RITU/genes_lists/spleenRMDown.csv")



















