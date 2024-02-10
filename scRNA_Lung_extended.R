# SCRNA LUNG 

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
# library(DoubletFinder)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# set working directory
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara")

# saveRDS(combined.all, "./combined_all_2.rds")
combined.all <- readRDS("./combined_all_2.rds")


# Read in Cell Calls
lung_Infected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara/lung_Inf_mixREF/gem_classification.csv")
lung_Infected_mouse <- subset(lung_Infected, lung_Infected$call=="mm10")

lung_NOTinfected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara/lung_NI_mixREF/gem_classification.csv")
lung_NOTinfected_mouse <- subset(lung_NOTinfected, lung_NOTinfected$call=="mm10")

# read in sample
lung_HIV = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara/lung_Inf_mixREF/filtered_feature_bc_matrix/")
dim(lung_HIV)
#68886  #5163

# remove mouse cells (lung_HIV)  
lung_HIV <- lung_HIV[,!colnames(lung_HIV)%in%lung_Infected_mouse$barcode]
dim(lung_HIV)
#68886   #5140

# Create Sample1 Seurat Object (Infected)
infected_lung <- CreateSeuratObject(counts = lung_HIV, project = "lung_infected", min.cells = 3, min.features = 200)
dim(infected_lung)
#24364  5096
head(infected_lung)

# read sample 2
lung_NI = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara/lung_NI_mixREF/filtered_feature_bc_matrix/")
dim(lung_NI)
#68886  6923

# remove mouse cells (lung_NI) 
lung_NI <- lung_NI[,!colnames(lung_NI)%in%lung_NOTinfected_mouse$barcode]
dim(lung_NI)
#68886     6870

# create sample2 suerat object
NotInf_lung <- CreateSeuratObject(counts = lung_NI, project = "NotInf_lung", min.cells = 3, min.features = 200)
dim(NotInf_lung)
# 28719  6727
head(NotInf_lung)

######################################################################################
# Identify Condition
infected_lung[["sample"]] <- "Infected"
NotInf_lung[["sample"]] <- "Not"

# Merge Samples
merged_seurat <- merge(x=infected_lung, y=NotInf_lung, add.cell.id = c("Infected", "Not"))

sum(grepl("MT",rownames(merged_seurat)))
#252
rownames(merged_seurat)[grepl("MT-",rownames(merged_seurat))]

# calculate the proportion of transcripts mapping to mitochondrial genes
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "MT-")

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


#save(merged_seurat, file="./merged_Lung_seurat.RData")
#load("./merged_Lung_seurat.RData")
max(merged_seurat@meta.data$percent.mt)
#87.24

# PLot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Filter based on quality plots 
Filter_lung <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt > -Inf & percent.mt <15)
VlnPlot(Filter_lung, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filter_lung, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Sums post Filtering
sum(Filter_lung$sample=="Infected")
#4903
sum(Filter_lung$sample=="Not")
#6369

Filter_lung@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",Filter_lung@assays$RNA@counts@Dimnames[[1]])
Filter_lung@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filter_lung@assays$RNA@counts@Dimnames[[1]])

######### Visualize Mitochondria
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
############

split_seurat <- SplitObject(Filter_lung, split.by = "sample")
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
################### Find Doublets in Infected Lung
# (pk identification) 
sweep.res.list_1 <- paramSweep_v3(split_seurat[[1]], PCs= 1:10, sct=T)
sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
bcmvn_1 <- find.pK(sweep.stats_1)
#pK = 0.07

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
nExp_poi <- round(0.038*nrow(split_seurat[[1]]@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder   -----edit pK
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[1]] <- doubletFinder_v3(split_seurat[[1]], PCs=1:10, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.07_186", sct=T)

# PLot Doublets & Singlets
DimPlot(split_seurat[[1]], reduction = "umap", group.by = c("DF.classifications_0.25_0.07_186",
                                                            "DF.classifications_0.25_0.07_165"),label = TRUE,repel = TRUE)
sum(split_seurat[[1]]$DF.classifications_0.25_0.07_186=="Doublet") #186

################### Doublet Finder for Non-Infected Lung
# (pk identification) 
sweep.res.list_2 <- paramSweep_v3(split_seurat[[2]], PCs= 1:10, sct=T)
sweep.stats_2 <- summarizeSweep(sweep.res.list_2, GT=F)
bcmvn_2 <- find.pK(sweep.stats_2)
# pK = 0.11
######## label Pk chart
pK=as.numeric(as.character(bcmvn_2$pK))
BCmetric=bcmvn_2$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
#######
# Homotypic Doublet Proportion Estimate
annotations <- split_seurat[[2]]@meta.data$SCT_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.049*nrow(split_seurat[[2]]@meta.data)) # Adjust according to Doublet Rate (8700 cells) 3.9%
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run Doublet Finder ---- edit pK
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
split_seurat[[2]] <- doubletFinder_v3(split_seurat[[2]], PCs=1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi.adj, reuse.pANN ="pANN_0.25_0.11_312", sct=T)

# Plot Doublets vs. Singlets
DimPlot(split_seurat[[2]], reduction = "umap", group.by = c("DF.classifications_0.25_0.11_312",
                                                            "DF.classifications_0.25_0.11_280"),label = TRUE,repel = TRUE)

# REMOVE DOUBLETS #
split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.07_186!="Doublet") 
dim(split_seurat[[1]])
#4717
split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset =DF.classifications_0.25_0.11_312 !="Doublet") 
dim(split_seurat[[2]])
#6057

########################################################## Downsample (4700 cells/sample)
set.seed(111)
split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 4700, replace=F)]
split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 4700, replace=F)]

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
save.image("./Combined_all_rdata.RData")
load("./Combined_all_rdata.RData")
#saveRDS(combined.all, "/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/SCT_integrated_lung_seurat.rds")

head(combined.all@meta.data)

# Dimensional Reduction
DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) # reduction = "pca",
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = 0.6)
combined.all <- RunTSNE(object = combined.all)
combined.all <- FindVariableFeatures(combined.all, selection.method = "vst", nfeatures = 2000)
combined.all[["percent.mt"]] <- PercentageFeatureSet(combined.all, pattern = "MT-",assay = "RNA")

VlnPlot(object = combined.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(combined.all@meta.data, 5)
rownames(combined.all)

########## Make tSNE and UMAP plots
DimPlot(combined.all, reduction = "tsne", split.by = "sample",label = TRUE,repel = TRUE)

DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", label = TRUE,repel = TRUE)

### For 2 samples run separately
# Lung infected
DimPlot(split_seurat[[1]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[1]], reduction = "umap", label = TRUE,repel = TRUE)
# Lung not infected
DimPlot(split_seurat[[2]], reduction = "tsne", label = TRUE,repel = TRUE)
DimPlot(split_seurat[[2]], reduction = "umap", label = TRUE,repel = TRUE)


# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(combined.all, assay='RNA', slot='counts')
writeMM(counts_matrix, file='/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/counts_Lung_organoid.mtx')
rownames(counts_matrix)

#Lung_infected_seurat <- split_seurat[[1]]
#Lung_NOTinfected_seurat <- split_seurat[[2]]
#saveRDS(Lung_infected_seurat, "/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/Lung_infected_QCdone_SCTmore.rds")
#saveRDS(Lung_NOTinfected_seurat, "/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/Lung_NOTinfected_QCdone_SCTmore.rds")

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

saveRDS(combined.all, "./combined_all_2.rds")
combined.all <- readRDS("./combined_all_2.rds")

# MANUALLY ANNOTATE
FeaturePlot(combined.all, cols= c("yellow", "blue"), label=T, features = c("CD8A", "CD8B", "CD4", "CD40LG", "CD3D", "CD3E", "CD19", "FCGR3A", "CD14"))


# Normalize and SCale RNA slot
DefaultAssay(combined.all) <- "RNA"
NormalizeData(combined.all, assay = "RNA")
ScaleData(combined.all, assay = "RNA")

DimPlot(combined.all, label=TRUE)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
################################ 10/13
# subcluster main annotations
combined.all <- SetIdent(combined.all, value = "seurat_clusters" )
DimPlot(combined.all, label=TRUE)
head(combined.all@meta.data)

nKcells <- WhichCells(combined.all, idents= c("1", "7", "4"))
Myeloid <- WhichCells(combined.all, idents= c("0", "8", "11", "12"))
Bcells <- WhichCells(combined.all, idents= "2")
C13 <- WhichCells(combined.all, idents= "13")
c13 <- subset(combined.all, cells = C13)
Tcells <- WhichCells(combined.all, idents= c("3", "6", "5", "9", "10"))

Tcell <- subset(combined.all, cells=Tcells)
Bcell <- subset(combined.all, cells=Bcells)
myeloid <- subset(combined.all, cells= Myeloid)
NKcells <- subset(combined.all, cells=nKcells)

combined.all@meta.data[nKcells, "more"] <- "NK"
combined.all@meta.data[Myeloid, "main"] <- "Myeloid"
combined.all@meta.data[Bcells, "more"] <- "B Cells"
combined.all@meta.data[C13, "main"] <- "C13"
combined.all@meta.data[Tcells, "main"] <- "T Cells"

# c13 - ambiguous? -> MAST
c13_marks <- FindMarkers(combined.all, ident.1 = 13)
FeaturePlot(c13, features = c("TPSB2", "CTSG", "TPSAB1", "MS4A2", 'HBD', "CPA3", "SLC18A2", "HPGDS", "KIT", "CD33", "FCGR2A"), cols=c("yellow", "blue"))+ NoLegend()
DoHeatmap(c13, features = c("TPSB2", "CTSG", "TPSAB1", "MS4A2", 'HBD', "CPA3", "SLC18A2", "HPGDS", "KIT", "CD33", "FCGR2A"), assay = "RNA", slot="data", angle = 90) + NoLegend()
combined.all@meta.data[C13, "more"] <- "MAST"

# t-cells
DefaultAssay(Tcell) <- "RNA"
Tcell <- RunPCA(Tcell, npcs = 30, verbose = FALSE)
Tcell <- RunUMAP(Tcell, dims = 1:20) 
Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:20)
Tcell <- FindClusters(Tcell, resolution = 1.0)
head(Tcell@meta.data)
Tcell <- SetIdent(Tcell, value = "integrated_snn_res.1")
DimPlot(Tcell, label=TRUE)

FeaturePlot(Tcell, cols= c("yellow", "blue"), label=T, features = c("CD8A", "CD8B", "CD4", "CD40LG"))



cd8naive <- WhichCells(Tcell, idents = c("0", "1", "6"))
cd4treg <- WhichCells(Tcell, idents = c("5"))
cd8em <- WhichCells(Tcell, idents = c("4", "2", "9"))
cd4eff <- WhichCells(Tcell, idents = c("8"))
cd4naive <- WhichCells(Tcell, idents = c("3"))
cd4cm <- WhichCells(Tcell, idents = c("7"))
NKgdT <- WhichCells(Tcell, idents = c("10"))
# combined.all$annotations <- "NA"
combined.all@meta.data[cd8naive, "annotations"] <- "CD8+ Naive"
combined.all@meta.data[cd8em, "annotations"] <- "CD8+ Eff/Mem"
combined.all@meta.data[cd4naive, "annotations"] <- "CD4+ Naive"
combined.all@meta.data[cd4cm, "annotations"] <- "CD4+ CM"
combined.all@meta.data[cd4treg, "annotations"] <- "CD4+ Treg"
combined.all@meta.data[cd4eff, "annotations"] <- "CD4+ EM"
combined.all@meta.data[NKgdT, "annotations"] <- "NK/gdT"
combined.all<- SetIdent(combined.all, value = "annotations")
DimPlot(combined.all, label=TRUE)

nkt_gdt <- WhichCells(Tcell, idents = "10")
cd4t <- WhichCells(Tcell, idents = c("4", "6", "7", "8", "9")) # 293 Cells
cd8t <- WhichCells(Tcell, idents = c("0", "1", "2", "3", "5", "11")) # 1198 Cells
cd4 <- subset(Tcell, cells = cd4t)
cd8 <- subset(Tcell, cells = cd8t)
NKT.GDT <- subset(Tcell, cells = nkt_gdt)

combined.all@meta.data[nkt_gdt, "more"] <- "NKT"

################################################### CD4
DefaultAssay(cd4) <- "integrated"
cd4 <- RunPCA(cd4, npcs = 30, verbose = FALSE)
cd4 <- RunUMAP(cd4, dims = 1:20) # reduction = "pca",
cd4 <- FindNeighbors(cd4, reduction = "pca", dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 1.0)
cd4<- SetIdent(cd4, value = "more")
DimPlot(cd4, label=T, label.size = 6)


DefaultAssay(cd4) <- "RNA"
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "S100A4", "IL32", "USP10", "TIMP1", "CD27", "LYAR"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("TIMP1", "LGALS1", "GZMK", "CST7", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("S100A4", "CCR7", "IL32", "ISG15"))
DoHeatmap(cd4, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B", "CCR5", "S100A4", "IL32", "ISG15"), assay = "RNA", slot="data", angle = 90) + NoLegend()
cd4treg <- WhichCells(cd4, idents = c("2"))
cd4cmn <- WhichCells(cd4, idents = c("0", "5", "7"))
cd4emtm <- WhichCells(cd4, idents = c("1", "3", '6'))
cd4temra <- WhichCells(cd4, idents = c("4"))
cd47 <- WhichCells(cd4, idents = c("7"))
combined.all@meta.data[cd4treg, "more"] <- "CD4+ TREG"
combined.all@meta.data[cd4temra, "more"] <- "CD4+ TEMRA"
combined.all@meta.data[cd4cmn, "more"] <- "CD4+ Naive/CM"
combined.all@meta.data[cd4emtm, "more"] <- "CD4+ EM/TM"
combined.all@meta.data[cd47, "more"] <- "CD4+ c7"

combined.allT <- combined.all
cd4naive <- WhichCells(cd4, idents = c("1", "4"))
cm <- WhichCells(cd4, idents = c("6"))
tm.em <- WhichCells(cd4, idents = c("2", "3"))
temra <- WhichCells(cd4, idents = c("5"))
treg <- WhichCells(cd4, idents = c("0"))
# combined.all$annotations <- "NA"
# Tcell$annotate <- "NA"
cd4@meta.data[temra, "annotate"] <- "CD4+ TEMRA"
cd4@meta.data[treg, "annotate"] <- "CD4+ TREG"
cd4@meta.data[cd4naive, "annotate"] <- "CD4+ Naive"
cd4@meta.data[cm, "annotate"] <- "CD4+ CM"
cd4@meta.data[tm.em, "annotate"] <- "CD4+ EM/TM"
cd4 <- SetIdent(cd4, value = "annotate")
DimPlot(cd4, label=T)
# cd4@meta.data <- cd4@meta.data[, -which(colnames(cd4@meta.data) %in% c('integrated_snn_res.0.6', "integrated_snn_res.1", "main", "finer"))]
cd4$celltype.cd4 <- paste(Idents(cd4), cd4$sample, sep = "_")
Idents(cd4) <- "celltype.cd4"
table(cd4@meta.data$celltype.cd4)
DimPlot(cd4, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

combined.allT@meta.data[temra, "annotations"] <- "CD4+ TEMRA"
combined.allT@meta.data[treg, "annotations"] <- "CD4+ TREG"
combined.allT@meta.data[naive, "annotations"] <- "CD4+ Naive"
combined.allT@meta.data[cm, "annotations"] <- "CD4+ CM"
combined.allT@meta.data[tm.em, "annotations"] <- "CD4+ EM/TM"
combined.allT <- SetIdent(combined.allT, value = "annotations")
DimPlot(combined.allT, label=T)


########################################### add annotatiions to combined.all
combined.all@meta.data[temra, "annotations"] <- "CD4+ TEMRA"
combined.all@meta.data[treg, "annotations"] <- "CD4+ TREG"
combined.all@meta.data[cd4naive, "annotations"] <- "CD4+ Naive"
combined.all@meta.data[cm, "annotations"] <- "CD4+ CM"
cdombined.all@meta.data[tm.em, "annotations"] <- "CD4+ EM/TM"
combined.all@meta.data[naive, "annotations"] <- "CD8+ Naive/CM"
combined.all@meta.data[tm_em_temra, "annotations"] <- "CD8+ EM/TM/TEMRA"

################################################################################################ cd8 
DefaultAssay(cd8) <- "integrated"
cd8 <- RunPCA(cd8, npcs = 30, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:20) # reduction = "pca",
cd8 <- FindNeighbors(cd8, reduction = "pca", dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 0.8)
cd8<- SetIdent(cd8, value = "more")
DimPlot(cd8, label=T, label.size = 6)

DefaultAssay(cd8) <- "RNA"
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3", "GZMA"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("FGFBP2", "GZMH", "GZMB", "GNLY"))

FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CD8A", "GZMK", "CCL5", "CCR5", "CD27", "FCGR3A"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("KLRC2", "TRDC", "ZNF683", "TRDV1", "CD3D", "CD3E"))
DoHeatmap(cd8, features = c("CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                            "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY",
                            "FGFBP2", "NKG7", "CCL4", "KLRC2", "TRDC", "ZNF683", "TRDV1", "CD3D", "CD3E"), assay = "RNA", slot="data", angle = 90) + NoLegend()



cd8naive <- WhichCells(cd8, idents = c("0", "1", "4"))
tm_em_temra <- WhichCells(cd8, idents = c("2", "3", '5'))
gdt <- WhichCells(cd8, idents = "6")
combined.all@meta.data[gdt, "more"] <- "gdT"
combined.all@meta.data[cd8naive, "more"] <- "CD8+ Naive/CM"
combined.all@meta.data[tm_em_temra, "more"] <- "CD8+ CCL5+"
cd8 <- SetIdent(cd8, value = "more")
DimPlot(cd8, label=T, label.size = 6)

# subset ccl5 cluster
cd8tem <- subset(cd8, cells = tm_em_temra)
DefaultAssay(cd8tem) <- "integrated"
cd8tem <- RunPCA(cd8tem, npcs = 30, verbose = FALSE)
cd8tem <- RunUMAP(cd8tem, dims = 1:20) # reduction = "pca",
cd8tem <- FindNeighbors(cd8tem, reduction = "pca", dims = 1:20)
cd8tem <- FindClusters(cd8tem, resolution = 0.8)
cd8tem<- SetIdent(cd8tem, value = "integrated_snn_res.0.8")
DimPlot(cd8tem, label=T, label.size = 6)
DoHeatmap(cd8tem, features = c("CD4", "CD40LG","CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                            "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY"), assay = "RNA", slot="data", angle = 90) + NoLegend()

######

cd8$celltype.cd8 <- paste(Idents(cd8), cd8$sample, sep = "_")
Idents(cd8) <- "celltype.cd8"
table(cd8@meta.data$celltype.cd8)
DimPlot(cd8, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)

#NKT/gdT
DefaultAssay(NKT.GDT) <- "RNA"
NKT.GDT <- RunPCA(NKT.GDT, npcs = 30, verbose = FALSE)
NKT.GDT <- RunUMAP(NKT.GDT, dims = 1:20) # reduction = "pca",
NKT.GDT <- FindNeighbors(NKT.GDT, reduction = "pca", dims = 1:20)
NKT.GDT <- FindClusters(NKT.GDT, resolution = 0.3)
NKT.GDT<- SetIdent(NKT.GDT, value = "nkt")
DimPlot(NKT.GDT, label=T, label.size = 6)

FeaturePlot(NKT.GDT, cols=c("yellow", "blue"), label=T, features = c("FGFBP2", "NKG7", "CD3E", "CD3D", "CD3G", "CCL5", 
                                                                     "CCL4", "GZMH", "CD8A", "GNLY", "S1PR1"))
DoHeatmap(NKT.GDT, features = c("FGFBP2", "NKG7", "CD3E", "CD3D", "CD3G", "CCL5", 
                                "CCL4", "GZMH", "CD8A", "GNLY", "IL7R", "XCL1", "KLRC2", "GZMB", "KLF2"), assay = "RNA", slot="data", angle = 90) + NoLegend()

NKT.GDT$nkt <- "NKT"
NKT.GDT$celltype <- paste(Idents(NKT.GDT), NKT.GDT$sample, sep = "_")
Idents(NKT.GDT) <- "celltype"
table(NKT.GDT@meta.data$celltype)
DimPlot(NKT.GDT, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)


################################################################################################# NK
DefaultAssay(NKcells) <- "SCT"
NKcells <- RunPCA(NKcells, npcs = 30, verbose = FALSE)
NKcells <- RunUMAP(NKcells, dims = 1:20) # reduction = "pca",
NKcells <- FindNeighbors(NKcells, reduction = "pca", dims = 1:20)
NKcells <- FindClusters(NKcells, resolution = 0.5)
NKcells<- SetIdent(NKcells, value = "integrated_snn_res.0.5")
DimPlot(NKcells, label=T, label.size = 6)

FeaturePlot(NKcells, cols=c("yellow", "blue"), label=T, features = c("NCAM1", "GNLY", "NKG7", "EOMES", "GZMB",
                                                                     "CD3D", "CD3E", "CCL3", "FGFBP2", "KIR2DL4", "IL7R", "XCL1"))
DoHeatmap(NKcells, features = c("NCAM1", "GNLY", "NKG7", "EOMES", "GZMB",
                            "CD3D", "CD3E", "CCL3", "FGFBP2", "KIR2DL4", "IL7R", "XCL1"), assay = "RNA", slot="data", angle = 90) + NoLegend()

NKcells$nk <- "NK"
NKcells$celltype <- paste(Idents(NKcells), NKcells$sample, sep = "_")
Idents(NKcells) <- "celltype"
table(NKcells@meta.data$celltype)
DimPlot(NKcells, reduction = "umap", split.by = "sample",repel = TRUE)


################################################################################################## B cells
DefaultAssay(Bcell) <- "integrated"
Bcell <- RunPCA(Bcell, npcs = 30, verbose = FALSE)
Bcell <- RunUMAP(Bcell, dims = 1:20) # reduction = "pca",
Bcell <- FindNeighbors(Bcell, reduction = "pca", dims = 1:20)
Bcell <- FindClusters(Bcell, resolution = 0.5)
Bcell<- SetIdent(Bcell, value = "integrated_snn_res.0.5")
DimPlot(Bcell, label=T, label.size = 6)

FeaturePlot(Bcell, cols=c("yellow", "blue"), label=T, features = c("CD19", "CD79A", "MS4A1", "VPREB3", "ITGAX", "SELL"))
DoHeatmap(Bcell, features = c("CD19", "CD79A", "MS4A1", "VPREB3", "ITGAX", "SELL"), assay = "RNA", slot="data", angle = 90) + NoLegend()

Bcell$b <- "B"
Bcell<- SetIdent(Bcell, value = "b")
DimPlot(Bcell, reduction = "umap", split.by = "sample",repel = TRUE)
Bcell$celltype <- paste(Idents(Bcell), Bcell$sample, sep = "_")
Idents(Bcell) <- "celltype"
table(Bcell@meta.data$celltype)

# MYELOID
DefaultAssay(myeloid) <- "integrated"
DefaultAssay(myeloid) <- "RNA"
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:20) 
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.6)
head(myeloid@meta.data)
myeloid <- SetIdent(myeloid, value = "more")
DimPlot(myeloid, label=TRUE)
FeaturePlot(myeloid, features = c( "RGS1", "PLD4", "BIRC3", "IDO1", "CDKN1C", "ADAM28", "FLT3", "SNCA"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features = c( "C1QA", "C1QB", "C1QC", "APOE", "APOC1"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features = c( "FCGR3A", "CD14", "CCR7", "JCHAIN", "CD1C", "CLEC9A", "CLEC10A", "FCN1", "LILRA4", "CX3CR1", "S100A12"), cols= c("yellow", "blue"), label = TRUE)
DoHeatmap(myeloid, features = c("FCGR3A", "CD14", "CCR7", "LILRA4",
                              "CX3CR1", "S100A12","RGS1", "PLD4", "BIRC3", "IDO1", "CDKN1C", "ADAM28",
                              "CX3CR1", "FLT3", "HLA-DRA", "CD86",  "LYZ", "S100A8", "S100A10", "CD74",
                              "IFI30", 'HLA-DPB1', "ITGAX", "CD33", "IL3RA", "NRP1", 'SELL', "CXCR2", "IRF8", "IL3RA", "ITGAX", "ITGAE", 
                              "HLA-DRA", "CCR2", "THBD", "ITGAM", "FCGR3A", "CD14", "CCR7", "BATF3", "IRF4","IRF7",  "KLF4", "ID2", "CLEC9A", 
                              "CLEC10A", "CD1A", "CD1C", "CD2", "FCER1A", "ITGAM", "MRC1", "BTLA", "CADM1", "CD68", "CSF3R","SIRPA", "CD4"), assay = "RNA", slot="data", angle = 90) + NoLegend()
FeaturePlot(myeloid, features = c( "HLA-DRA", "CD86", "LYZ", "CD74",
                                   "IFI30", 'HLA-DPB1', 'CPVL'), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features = c( "ITGAM", "ITGAX", "CD13", "CD33", "CD1C", "THBD", "BATF3"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features = c("NRP1", 'CD1C', "IL3RA", "ITGAX","FCN1", "CLEC10A", "FCGR3A", "CD14"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features =c("IRF8", "IL3RA", "ITGAX", "ITGAE", "HLA-DRA", "CCR2", "THBD", "ITGAM", "FCGR3A", "CD14", "CCR7", "BATF3", "IRF4", "KLF4", "ID2", "IRF7"),cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(myeloid, features =c("CD14", "CLEC9A", "ID2", "CLEC10A", "CD1A", 'CD1C', "CD2", "FCER1A", "ITGAM", "MRC1", "BTLA", "CADM1"),cols= c("yellow", "blue"), label = TRUE)

myeloid6_markers <- FindMarkers(object = myeloid, ident.1 = "6")
myeloid8_markers <- FindMarkers(object = myeloid, ident.1 = "8")


cd14 <- WhichCells(myeloid, idents = c("3", "2"))
cd16 <- WhichCells(myeloid, idents=c("0", "1"))
c8 <- WhichCells(myeloid, idents=c("8"))
mdc <- WhichCells(myeloid, idents=c("7"))
imono <- WhichCells(myeloid, idents=c("4", "5"))
neut <- WhichCells(myeloid, idents=c("6"))
combined.all@meta.data[cd14, "more"] <- "CD14 MONO"
combined.all@meta.data[cd16, "more"] <- "CD16 MONO"
combined.all@meta.data[c8, "more"] <- "CD16 MONO"
combined.all@meta.data[imono, "more"] <- "CD14+/16+ MONO"
combined.all@meta.data[mdc, "more"] <- "Mye DC"
combined.all@meta.data[neut, "more"] <- "CD16 MONO"

DimPlot(myeloid, label=T, label.size = 5, split.by = "sample")
myeloid$celltype <- paste(Idents(myeloid), myeloid$sample, sep = "_")
Idents(myeloid) <- "celltype"
table(myeloid@meta.data$celltype)

#myeloid@meta.data$annotate <- "NA"
myeloid@meta.data[cd14, "annotate"] <- "CD14 MONO"
myeloid@meta.data[cd16, "annotate"] <- "CD16 MONO"
myeloid@meta.data[dc, "annotate"] <- "cDC"
myeloid@meta.data[imono, "annotate"] <- "CD14+/16+ MONO"
myeloid@meta.data[mdc, "annotate"] <- "CD16 MONO"
myeloid@meta.data[neut, "annotate"] <- "CD16 MONO"
myeloid <- SetIdent(myeloid, value = "annotate")
DimPlot(myeloid, label=T)
myeloid$celltype.myeloid <- paste(Idents(myeloid), myeloid$sample, sep = "_")
Idents(myeloid) <- "celltype.myeloid"
table(myeloid@meta.data$celltype.myeloid)
DimPlot(myeloid, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)


# CLuster 11
C11 <- WhichCells(combined.all, idents= "11")
c11_marks <- FindMarkers(combined.all, ident.1 = 11)
write.csv(c11_marks, "./c11_markers.csv")
c11 <- subset(combined.all, cells = C11)
# GATA2 stat1
DoHeatmap(combined.all, features = c("FCGR3A", "CD14", "CCR7", "JCHAIN", "CD1C", "CLEC10A", "LILRA4",
                                "CX3CR1", "S100A12","CD34", "CD38", "PAX5", "GATA1", "GATA2", "STAT5", "CD4", "CD3D", "CD3E", 
                                "GNLY", "NKG7", "CD8A", "CD8B", "CCR7", "CD24", "CD19", "KLRG1", "IL7R", "CD34", "STAT1", "BLVRB", "PLEK", "TAL1",
                            "HBG1", "HBG2", "ITGA2B", "CA1", "HBD", "TFRC", "TPSB2", "CTSG", "HBGD", "TPSAB1", "CPA3", "KIT"), assay = "RNA", slot="data", angle = 90) + NoLegend()
DoHeatmap(c11, features = c("TPSB2", "CTSG", "HBGD", "TPSAB1", "CPA3", "KIT", "MS4A2", "HBD", "SLC18A2", "HPGDS", "CPA3"), assay = "RNA", slot="data", angle = 90) + NoLegend()
FeaturePlot(c11, label=T, features = c("TPSB2", "CTSG","TPSAB1", "CPA3", "HBD", "MS4A2","SLC18A2", "HPGDS", "CPA3"), cols=c("yellow", "blue"))


combined.all@meta.data[C11, "annotations"] <- "MAST"









######################################################################## Annotate Clusters using SINGLE R
# get correct matrixStats version
#require(devtools)
#install_version("matrixStats", version = "0.60.1", repos = "http://cran.us.r-project.org")

# Install SingleR
#BiocManager::install("SingleR")
library(SingleR)
library(celldex)

# establish reference
dice.ref <- DatabaseImmuneCellExpressionData()
sce <- as.SingleCellExperiment(combined.all)
combined.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
combined.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
table(combined.main$pruned.labels)
table(combined.fine$pruned.labels)
combined.all@meta.data$combined.main <- combined.main$pruned.labels
combined.all <- SetIdent(combined.all, value = "combined.main")
DimPlot(combined.all, label = T , repel = T, label.size = 3) + NoLegend()
combined.all <-  subset(combined.all, combined.fine != "")

########## Make tSNE and UMAP plots
combined.all <- SetIdent(combined.all, value = "combined.main")

DimPlot(combined.all, reduction = "tsne", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap",label = TRUE,repel = TRUE)

##############################################################################################################
######## read in interested genes
genes_intersted <- readxl::read_xlsx("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/interested_genes_scRNA.xlsx", col_names = F)
colnames(genes_intersted) <- "Genes"

##### violin plots for genes interested
combined.all <- SetIdent(combined.all, value = "annotations")
DefaultAssay(combined.all) <- "RNA"
p.violin <- VlnPlot(object = combined.all, 
                    features = genes_intersted$Genes[48:50],
                    split.by = "sample",
                    pt.size = 0.2,
                    combine = TRUE) + theme(legend.position = 'right')

# More Visualizations
RidgePlot(combined.all,features=  genes_intersted$Genes[1:15], stack = TRUE)
DotPlot(combined.all, assay = "RNA",  features=  genes_intersted$Genes[1:23], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()
DotPlot(combined.all, assay = "RNA",  features=  genes_intersted$Genes[24:49], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()

DoHeatmap(combined.all, features = genes_intersted$Genes, assay = "RNA", slot="data", angle = 90) + NoLegend()



######################################################################################################## EXTRA
##### Dot PLot to visualize genes of interest across conditions
DotPlot(combined.all, features = genes_intersted$Genes[1:23], cols = c("blue", "red"), dot.scale = 8, split.by = "sample") +
  RotatedAxis()
DotPlot(combined.all, features = genes_intersted$Genes[24:49], cols = c("blue", "red"), dot.scale = 8, split.by = "sample") +
  RotatedAxis()

#############################################################################################
Idents(combined.all) <- "more"
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", label = TRUE,repel = TRUE, cells = Myeloid, label.size = 6)

combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all@meta.data$celltype.stim)

NormalizeData(combined.all, assay = "RNA")
ScaleData(combined.all, assay = "RNA")

cd4EM_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ EM_Infected", ident.2 = "CD4+ EM_Not", verbose = FALSE)
head(cd4EM_DEgenes, n = 15)


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
cd4EM <- subset(combined.all, idents = "CD4+ EM")
Idents(cd4EM) <- "sample"
avg.cd4EM <- as.data.frame(log10(AverageExpression(cd4EM, verbose = FALSE)$RNA))
avg.cd4EM$gene <- rownames(avg.cd4EM)
genes.to.label = c("MALAT1", "S100A6", "IFI27", "CXCL13", "S100A4")
p1 <- ggplot(avg.cd4EM, aes(Infected, Not)) + geom_point() + ggtitle("CD4+ EM")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd14_DEgenes <- FindMarkers(combined.all, ident.1 = "CD14 MONO_Infected", ident.2 = "CD14 MONO_Not", verbose = FALSE)

# Plot Average Expression
cd14M <- subset(combined.all, idents = "CD14 MONO")
Idents(cd14M) <- "sample"
avg.cd14M <- as.data.frame(log10(AverageExpression(cd14M, verbose = FALSE)$RNA))
avg.cd14M$gene <- rownames(avg.cd14M)
genes.to.label = c("IFI27", "CXCL10", "B2M")
p1 <- ggplot(avg.cd14M, aes(Infected, Not)) + geom_point() + ggtitle("CD14 MONO")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


############################## April Visualizations
combined.all <- SetIdent(combined.all, value = "main")
DimPlot(combined.all, label=T, label.size =5, pt.size = 1)
table(combined.all$sample)
inf <- subset(combined.all, ident = "Infected")
not <- subset(combined.all, ident = "Not")
table(inf$main)
table(inf$more)
table(not$main)
table(not$more)
mast <- WhichCells(combined.all, idents = "C13")
combined.all@meta.data[mast, "main"] <- "Mast"

DotPlot(combined.all, assay = "RNA",  features=  inf_deg$gene[1:23], cols = c("blue", "red"), dot.scale=8, split.by = "sample") + RotatedAxis()
DotPlot(combined.all, assay = "RNA",  features=  inf_deg$gene[24:49], cols = c("blue", "red"), dot.scale=8) + RotatedAxis()
library(SCpubr)
combined.all <- SetIdent(combined.all, value = "main")
de_genes <- FindMarkers(combined.all, ident.1 = "Infected", ident.2 = "Not", assay = "RNA", slot = "data")
SCpubr::do_VolcanoPlot(sample = combined.all,
                       de_genes = de_genes,
                       n_genes = 5,
                       pval_cutoff = 1.3,
                       FC_cutoff = 1.5,
                       colors.use = c("red"),
                       font.size=12)
SCpubr::do_ExpressionHeatmap(sample = combined.all,
                             features = row.names(deg_b)[1:20],
                             viridis_direction = -1,
                             group.by = "sample",
                             flip= T)

# Find DEGs
combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all@meta.data$celltype.stim)
deg_b <- FindMarkers(combined.all, ident.1 = "B Cells_Infected", ident.2 = "B Cells_Not", assay = "RNA", slot = "data")
deg_mye <- FindMarkers(combined.all, ident.1 = "Myeloid_Infected", ident.2 = "Myeloid_Not", assay = "RNA", slot = "data")
deg_nk <- FindMarkers(combined.all, ident.1 = "NK_Infected", ident.2 = "NK_Not", assay = "RNA", slot = "data")
deg_t <- FindMarkers(combined.all, ident.1 = "T Cells_Infected", ident.2 = "T Cells_Not", assay = "RNA", slot = "data")
deg_nkt <- FindMarkers(combined.all, ident.1 = "NKT_Infected", ident.2 = "NKT_Not", assay = "RNA", slot = "data")




