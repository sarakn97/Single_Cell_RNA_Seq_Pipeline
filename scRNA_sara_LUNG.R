# SCRNA LUNG 

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# set working directory
setwd("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara")

# Read in Cell Calls
lung_Infected <- read.csv("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/lung_Inf_mixREF/gem_classification.csv")
lung_Infected_mouse <- subset(lung_Infected, lung_Infected$call=="mm10")

lung_NOTinfected <- read.csv("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/lung_NI_mixREF/gem_classification.csv")
lung_NOTinfected_mouse <- subset(lung_NOTinfected, lung_NOTinfected$call=="mm10")

# read in sample
lung_HIV = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/lung_Inf_mixREF/filtered_feature_bc_matrix/")
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
lung_NI = Read10X("/storage1/fs1/leyao.wang/Active/scRNA_Lung_sara/lung_NI_mixREF/filtered_feature_bc_matrix/")
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

DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE)
combined.all <- RunUMAP(combined.all, dims = 1:20) # reduction = "pca",
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20)
combined.all <- FindClusters(combined.all, resolution = 0.5)
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




#####################################################################################################################
FeaturePlot(combined.all, features = "CD8A")
table(combined.all$integrated_snn_res.0.5)
combined_markers<- FindAllMarkers(object = combined.all, 
                                        only.pos = TRUE,
                                        logfc.threshold = 0.25)
combined_markers<- combined_markers[ , c(6, 7, 1:5)]
combined_markers <- combined_markers %>%
  dplyr::arrange(cluster, p_val_adj)


cluster5 <- names(which(combined.all@meta.data$integrated_snn_res.0.5 == "5"))
C5 <- subset(combined.all, integrated_snn_res.0.5 == "5")
DimPlot(combined.all, cells.highlight = cluster5)

DefaultAssay(combined.all) <- "integrated"
FetchData(combined.all, vars=c("CD8A", "CD3D", ""), slot = "data")
# myeloid cells
AverageExpression(object= combined.all, assays = "RNA", features= c("LYZ", "AIF1", "COTL1"))
# B Cells
AverageExpression(object= combined.all, assays = "RNA", features= c("CD19", "CD20", "CD34", "CD27", "MS4A1", "MZB1"))

NKcells <- WhichCells(combined.all, idents= c("1", "7", "9"))
Monocytes <- WhichCells(combined.all, idents= c("0", "6"))
Bcells <- WhichCells(combined.all, idents= "2")
combined.all@meta.data[Bcells, 19] <- "B Cells"
combined.all@meta.data$fine <- "" #21

combined.all@meta.data[NKcells, 20] <- "NK"

# T-CELLS
combined.all<- SetIdent(combined.all, value = "seurat_clusters")
AverageExpression(object= combined.all, assays = "RNA", features= c("CD3D", "CD4", "CD8A", "CCR6", "CD8B", "CD40LG"))
Tcells <- WhichCells(combined.all, idents= c("3", "4", "5", "8", "10"))
combined.all@meta.data[Tcells, 19] <- "T Cells"
tcells <- subset(combined.all, cells= Tcells)
CD8A_expression = GetAssayData(object = tcells, 
                               assay = "RNA", slot = "data")["CD8A",]
CD4_expression = GetAssayData(object = tcells, 
                              assay = "RNA", slot = "data")["CD4",]
both <- names(which(CD4_expression >0 & CD8A_expression >0)) # 35 express both
neither <- names(which(CD4_expression==0 & CD8A_expression==0)) #1048
cd8 <- names(which(CD8A_expression >0))
cd4 <- names(which(CD4_expression >0))
cd8tcells <- cd8[!(cd8 %in% both)] 
cd4tcells <- cd4[!(cd4 %in% both)] 
CD4 <- subset(combined.all, cells= cd4tcells)
CD8 <- subset(combined.all, cells= cd8tcells)
combined.all@meta.data[cd4tcells, 20] <- "CD4+ T-Cell"
DimPlot(combined.all, cells.highlight = cd4tcells)

# BOTH
BOTH <- subset(combined.all, cells=both)
combined.all@meta.data[both, 20] <- "CD4+/CD8+ T-Cell"
FetchData(combined.all, slot = "data", cells = both, vars=c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB", "GNLY", "GZMB", "CD7", "DUSP2"))
DimPlot(combined.all, cells.highlight = both)
sell = GetAssayData(object = BOTH, 
                    assay = "RNA", slot = "data")["SELL",]
lef1 = GetAssayData(object = BOTH, 
                    assay = "RNA", slot = "data")["LEF1",]
ltb = GetAssayData(object = BOTH, 
                   assay = "RNA", slot = "data")["LTB",]
CD4N <- names(which(sell>0 & lef1>0 & ltb>0))
combined.all@meta.data[CD4N, 20] <- "CD4+ T-Cell"
combined.all@meta.data[CD4N, 21] <- "CD4+ naive"
both <- both[!(both %in% CD4N)]
gnly = GetAssayData(object = BOTH, 
                    assay = "RNA", slot = "data")["GNLY",]
nkg7 = GetAssayData(object = BOTH, 
                    assay = "RNA", slot = "data")["NKG7",]
cd8a = GetAssayData(object = BOTH, 
                   assay = "RNA", slot = "data")["CD8A",]
cd8 <- names(which(gnly>6 & nkg7 >6 & cd8a >4))
both <- both[!(both %in% cd8hh)]

#classify remaining

#NEITHER - SEptember 8
NEITH <- subset(combined.all, cells= neither)
NEITH <- SetIdent(NEITH, value = "seurat_clusters")
DimPlot(combined.all, cells.highlight = neither)
AverageExpression(object= NEITH, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB"))
GZMH <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["GZMH",]
ITGA1 <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["ITGA1",]
GNLY <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["GNLY",]
CD8B <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["CD8B",]
cd3d <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["CD3D",]
cd3e <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["CD3E",]
cd3g <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["CD3G",]
names(which(CD8B>0)) #192
names(which(ITGA1>0)) #190
names(which(GNLY >0)) #360
names(which(GZMH >0)) #206
CD3d <- names(which(cd3d>0)) #770
CD3e <- names(which(cd3e>0)) #795
CD3g <- names(which(cd3g>0)) #757
names(which(cd3d==0 & cd3e ==0 & cd3g==0)) #109
Lef1 = GetAssayData(object = NEITH, 
                    assay = "RNA", slot = "data")["LEF1",]
cd28 = GetAssayData(object = NEITH, 
                    assay = "RNA", slot = "data")["CD28",]
names(which(Lef1>0)) #608
names(which(cd28>0)) #487
CD4NorTREG <- names(which(Lef1>0 & cd28>0)) #327 CD4 T-Cells (Naive or TREG)
combined.all@meta.data[CD4NorTREG, 20] <- "CD4+ T-Cell"
neither <- neither[!(neither %in% CD4NorTREG)]
NEITH <- subset(combined.all, cells= neither) # 721 remaining
cd40lg = GetAssayData(object = NEITH, 
                      assay = "RNA", slot = "data")["CD40LG",]
cd40lgenes <- names(which(cd40lg>0)) #111
cd4ORmait <- subset(NEITH, cells=cd40lgenes)
AverageExpression(cd4ORmait, assays="RNA", features= c("CCR6", "LEF1", "IL7R", "CD8B", "SLC4A10", "KLRB1", "CCR7", "TRAV12")) #
ccr7 = GetAssayData(object = cd4ORmait, 
                      assay = "RNA", slot = "data")["CCR7",] 
klrb1 = GetAssayData(object = cd4ORmait, 
                      assay = "RNA", slot = "data")["KLRB1",]
mait <- names(which(klrb1>0 & ccr7==0)) #62
DimPlot(combined.all, cells.highlight = mait)
combined.all@meta.data[mait, 20] <- "MAIT T-Cell"
neither <- neither[!(neither %in% mait)]
NEITH <- subset(combined.all, cells= neither) # 677 remaining
foxp3 = GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["FOXP3",]
ccr4 = GetAssayData(object = NEITH, 
                    assay = "RNA", slot = "data")["CCR4",]
il2ra = GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["IL2RA",]
names(which(foxp3>0)) #38
names(which(il2ra>0)) #60
cd4cells <- names(which(il2ra>0 | foxp3>0 | ccr4>0)) #98 T-CELLS CD4 (TREG)
DimPlot(combined.all, cells.highlight = cd4cells)
combined.all@meta.data[cd4cells, 20] <- "CD4+ T-Cell"
neither <- neither[!(neither %in% cd4cells)]
NEITH <- subset(combined.all, cells= neither) # 579 remaining

klrc2 = GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["KLRC2",]
trdv1 = GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["TRDV1",]
NK_or_gdT <- names(which(klrc2>0)) #98
NK_gdT <- subset(NEITH, cells= NK_or_gdT)
znf683 <- GetAssayData(object = NK_gdT, 
                       assay = "RNA", slot = "data")["ZNF683",]
GDTznf <- names(which(znf683>0))
combined.all@meta.data[GDTznf, 20] <- "gdT"
NK_or_gdT <- NK_or_gdT[!(NK_or_gdT %in% GDTznf)] #90
NArem <- names(which(NK_gdT$finer != "gdT")) #remove those already identified as gDT from september 2 below
NK_gdT<- subset(combined.all, cells= NArem) 
NK_gdT <- SetIdent(NK_gdT, value = "seurat_clusters")
AverageExpression(NK_gdT, slot="data", assays= "integrated", features = c("CCL5", "KLRC2", "TRDV1", "GNLY", "CXCR31", "KLRB2", "DUSP2", "NCR1"))
DimPlot(combined.all, cells.highlight = NArem)
gdT2 <- WhichCells(NK_gdT, idents= "10")
combined.all@meta.data[gdT2, 20] <- "gdT"
NArem <- NArem[!(NArem %in% gdT2)] #28
NK_gdT<- subset(combined.all, cells= NArem) # 28 cells
AverageExpression(NK_gdT, slot="data", assays= "integrated", features = c("CCL5", "CD8B", "KLRC2", "TRDV1", "GNLY", "CXCR31", "KLRB2", "DUSP2", "NCR1"))
ccl5 = GetAssayData(object = NK_gdT, 
                     assay = "RNA", slot = "data")["CCL5",]
klrc2 = GetAssayData(object = NK_gdT, 
                     assay = "RNA", slot = "data")["KLRC2",]
gdt3 <- names(which(ccl5>0 & klrc2>0))
combined.all@meta.data[gdt3, 20] <- "gdT"
NArem <- NArem[!(NArem %in% gdt3)] # 3 left
NK_gdT<- subset(combined.all, cells= NArem) # 3 cells - NK or gdT ???

neither <- neither[!(neither %in% NK_or_gdT)]
NEITH <- subset(combined.all, cells= neither) #481
DimPlot(combined.all, cells.highlight= neither)



# Gamma Delta - September 2
klrc2 <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["KLRC2",]
trdv1 <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["TRDV1",]
TRDC <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["TRDC",]
GDT <- names(which((klrc2>0 & trdv1>0) |(klrc2>0 & TRDC>0) |(trdv1>0 & TRDC>0))) #110
DimPlot(combined.all, cells.highlight = GDT)
gdt <- subset(combined.all, cells=GDT)
#NEITHER
gdt <- SetIdent(gdt, value = "seurat_clusters")
AverageExpression(object= gdt, assays = "RNA", features= c("CCL5", "GNLY", "CD7", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "KLRC2", "TRDV1", "TRDC"))
names(which(trdv1>0)) #30
names(which(TRDC>0))  #192
WhichCells(gdt, idents = "3")
cd40lg <- GetAssayData(object =gdt, 
                     assay = "RNA", slot = "data")["CD40LG",]
names(which(cd40lg>0))
gdt5 <- WhichCells(gdt, idents = "5")
GDT5 <- subset(gdt, cells= gdt5)
AverageExpression(object= GDT5, assays = "RNA", features= c("CCL5", "GNLY", "CD7", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "KLRC2", "TRDV1", "TRDC"))
gdt10 <- WhichCells(gdt, idents = "10")
combined.all@meta.data[gdt10, 20] <- "gdT"
gdtmaybe <- GDT[!(GDT %in% gdt10)]

GDT8 <- subset(gdt, cells= gdt8)
AverageExpression(object= GDT8, assays = "RNA", features= c("CCL5", "GNLY", "CD7", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "KLRC2", "TRDV1", "TRDC"))
FetchData(gdt, cells = gdt10,  vars =c("CCL5", "GNLY", "CD7", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "KLRC2", "TRDV1", "TRDC"), slot="data" )
cd40lg <- GetAssayData(object = NEITH, 
                     assay = "RNA", slot = "data")["CD40LG",]
cd40posNEITH <- names(which(cd40lg>0))
C40PN <- subset(combined.all, cells=cd40posNEITH)
AverageExpression(object= C40PN, assays = "RNA", features=c("CD3E", "CD3D", "CXCR6",
                                                            "CCL5", "PFN1", "LTB","IL7R", "IL32",
                                                            "CCL4", "GZMB", "GZMH", 'NKG7'))

combined.all@meta.data[cd40posNEITH, 20] <- "CD4+ T-Cell"

unident <- WhichCells(combined.all, idents = "NA")
table(combined.all@meta.data$finer)

NAmarks <-FindConservedMarkers(combined.all, ident.1 = "NA", grouping.var = "sample", verbose= FALSE)







#saveRDS(combined.all, "./combinedall91.rds")
#saveRDS(NAmarks, "./NAmarks.rds")
# separate infected vs. non-infected
# run blood scRNA
# FindAllMarkers()
# separate NK
NAmarks <- readRDS("./NAmarks.rds")


#### CD4
AverageExpression(object= CD4, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB"))
AverageExpression(object= CD8, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB"))

sell = GetAssayData(object = CD4, 
                               assay = "RNA", slot = "data")["SELL",]
lef1 = GetAssayData(object = CD4, 
                              assay = "RNA", slot = "data")["LEF1",]
ltb = GetAssayData(object = CD4, 
                    assay = "RNA", slot = "data")["LTB",]

CD4naive <- names(which(sell >0 & ltb >0 & lef1 >0))
combined.all@meta.data[CD4naive, 21] <- "CD4+ naive"
cd4rem <- cd4tcells[!(cd4tcells %in% CD4naive)] 
CD4rem <- subset(CD4, cells = cd4rem)
AverageExpression(object= CD4rem, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB"))
sell = GetAssayData(object = CD4, 
                    assay = "RNA", slot = "data")["SELL",]
lef1 = GetAssayData(object = CD4, 
                    assay = "RNA", slot = "data")["LEF1",]
ltb = GetAssayData(object = CD4, 
                   assay = "RNA", slot = "data")["LTB",]

CD4remC5 <- WhichCells(CD4rem, idents = "5")
FetchData(combined.all, vars = c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB", "IL32", "PFN1", "CCL5"), slot = "data", cells = CD4remC5)
CD4trm <- CD4remC5
combined.all@meta.data[CD4trm, 21] <- "CD4+ TRM"
cd4rem <- cd4rem[!(cd4rem %in% CD4remC5)] 
CD4rem <- subset(CD4rem, cells = cd4rem)
AverageExpression(object= CD4rem, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB"))
sell = GetAssayData(object = CD4rem, 
                    assay = "RNA", slot = "data")["SELL",]
CD4remSELL <- names(which(sell >0))
FetchData(combined.all, vars = c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB", "IL32", "PFN1", "CCL5", "CXCR6", "ITGAE"), slot = "data", cells = CD4remSELL)
CD4naive2 <- CD4remSELL
combined.all@meta.data[CD4naive2, 21] <- "CD4+ naive"
cd4rem <- cd4rem[!(cd4rem %in% CD4naive2)] 
CD4rem <- subset(CD4rem, cells = cd4rem)
AverageExpression(object= CD4rem, assays = "RNA", features= c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB", "IL32", "PFN1", "CCL5", "CXCR6", "ITGAE"))
FetchData(combined.all, vars = c("CD4", "CD8A", "CD8B", "CD40LG", "SRGN", "CXCR6", "GZMH", "NKG7", "ZNF683", "CD28", "KLF2", "S1PR1", "IL7R", "SELL", "LEF1", "LTB", "IL32", "PFN1", "CCL5", "CXCR6", "CCR4", "CCR6", "CCR5"), slot = "data", cells = cd4rem)
CD4em <- WhichCells(CD4rem, idents = "8")
combined.all@meta.data[CD4em, 21] <- "CD4+ effector"
cd4rem <- cd4rem[!(cd4rem %in% trm2)] 
combined.all@meta.data[jajajja, 21] <- "CD4+ naive"

#treg
naive <- names(which(CD4$fine == "CD4+ naive"))
cd4Naive <- subset(combined.all, cells = naive)
AverageExpression(object = cd4Naive, assays = "RNA", features=c( "CCR4", "CTLA4", "IL2RA", "TIGIT", "CCL5", "CD40LG"))
foxp3 = GetAssayData(object = cd4Naive, 
                    assay = "RNA", slot = "data")["FOXP3",]
treg1 <- names(which(foxp3>0)) ############################
naive <- naive[!(naive %in% treg1)] 
CD4naiveC3 <- WhichCells(cd4Naive, idents = "3")
naive <- naive[!(naive %in% CD4naiveC3)] 
il2ra = GetAssayData(object = cd4Naive, 
                     assay = "RNA", slot = "data")["IL2RA",]
ctla4 = GetAssayData(object = cd4Naive, 
                     assay = "RNA", slot = "data")["CTLA4",]
treg2 <- names(which(il2ra >0 & ctla4 >0 )) #################################
naive <- naive[!(naive %in% treg2)]  
ccr4 = GetAssayData(object = cd4Naive, 
                     assay = "RNA", slot = "data")["CCR4",]
tigit = GetAssayData(object = cd4Naive, 
                    assay = "RNA", slot = "data")["TIGIT",]
treg3 <- names(which((il2ra >0 & ccr4 >0) | (tigit >0 & ctla4 >0))) #########
naive <- naive[!(naive %in% treg3)]  
CD4naive4 <- WhichCells(cd4Naive, idents = "4")
FetchData(cd4Naive, vars= c("CCR4", "CTLA4", "IL2RA", "TIGIT", "CCL5", "CD40LG"), slot = "data", cells = CD4naive4)
CD4naive8 <- WhichCells(cd4Naive, idents = "8")
FetchData(cd4Naive, vars= c("CCR4", "CTLA4", "IL2RA", "TIGIT", "CCL5", "CD40LG", "CD62L", "CCR7"), slot = "data", cells = CD4naive8)


combined.all <- SetIdent(combined.all, value= "seurat_clusters")
DimPlot(combined.all, cells.highlight = cd8tcells)
# save.image("./clusterannotations831.RData")

# CD8 T-CELL
head(combined.all@meta.data)
cd8tcells <- cd8[!(cd8 %in% both)] 
CD8 <- subset(combined.all, cells= cd8tcells)
combined.all@meta.data[cd8tcells, 20] <- "CD8+ T-Cell"
AverageExpression(CD8, assays="RNA", features=c("CD8A", "CD8B", "CXCR6", "ITGA1", "ITGAE", "IL7R",
                                                "GZMH", "KLRG1", "GNLY", "NKG7", "CCL5", "EOMES", "TIGIT", "TNIP3", "LYST", 
                                                "CRTAM", "GZMA", "GZMB", "GZMK", "ZNF683", "PFN1", "CD7", "SRGN"))
# CLUSTER 3 LOOKS LIKE EFECTOR CD8










####### B CELLS
Bcells <- WhichCells(combined.all, idents= "2")
combined.all@meta.data[Bcells, 19] <- "B Cells"
B <- subset(combined.all, cells= Bcells)
CD19 = GetAssayData(object = B, 
                    assay = "RNA", slot = "data")["CD19",]
CD479A = GetAssayData(object = B, 
                      assay = "RNA", slot = "data")["CD79A",]
IGHD = GetAssayData(object = B, 
                    assay = "RNA", slot = "data")["IGHD",]
FCER2 = GetAssayData(object = B, 
                     assay = "RNA", slot = "data")["FCER2",]
TCL1A = GetAssayData(object = B, 
                     assay = "RNA", slot = "data")["TCL1A",]
MS4A1 = GetAssayData(object = B, 
                     assay = "RNA", slot = "data")["MS4A1",]

bnaive <- names(which(IGHD >0 & TCL1A >0)) #812
Bnotnaive <- Bcells[!(Bcells %in% bnaive)] #323
BnotNaive <- subset(B, cells = Bnotnaive)
MZB1 = GetAssayData(object = BnotNaive, 
                    assay = "RNA", slot = "data")["MZB1",]
plasmaB <- names(which(MZB1 >0)) #178
PlasmaB <- subset(combined.all, cells=plasmaB)
bmem <-  Bnotnaive[!(Bnotnaive %in% plasmaB)] #145
Bmem <- subset(combined.all, cells= bmem)
AverageExpression(PlasmaB, assays = "RNA", features = c("MZB1", "JCHAIN", "IGHA1", "CD27", "TNFRSF17", "CD19", "IGHD", "TCL1A", "MS4A1", "BANK1"))
MZB1 = GetAssayData(object = PlasmaB, 
                      assay = "RNA", slot = "data")["MZB1",]
TCL1A = GetAssayData(object = PlasmaB, 
                     assay = "RNA", slot = "data")["TCL1A",]
IGHD = GetAssayData(object = PlasmaB, 
                     assay = "RNA", slot = "data")["IGHD",]
names(which(MZB1>0 ))
names(which(IGHD>0))
names(which(TCL1A>0))
FetchData(object= PlasmaB, slot="data", vars= c("TCL1A", "IGHD"))
DimPlot(combined.all, cells.highlight = bnaive, reduction = "umap",label = TRUE,repel = TRUE)
AverageExpression(Bmem, assays = "RNA", features = c("MZB1", "JCHAIN", "IGHA1", "CD27", "TNFSF13B", "CD19", "IGHD", "TCL1A", "CD79A", "MS4A1", "BANK1"))
combined.all@meta.data[bnaive, 20] <- "B naive"
combined.all@meta.data[bmem, 20] <- "B mem"
combined.all@meta.data[plasmaB, 20] <- "B plasma"


############### MYELOID
MYELOID <- WhichCells(combined.all, idents= c("0", "6"))
MY <- subset(combined.all, cells= MYELOID)

CLEC10A = GetAssayData(object = MY, 
                       assay = "RNA", slot = "data")["CLEC10A",]
CLEC9A = GetAssayData(object = MY, 
                      assay = "RNA", slot = "data")["CLEC9A",]
CD1C = GetAssayData(object = MY, 
                     assay = "RNA", slot = "data")["CD1C",]
CCR7 = GetAssayData(object = MY, 
                    assay = "RNA", slot = "data")["CCR7",]
JCHAIN = GetAssayData(object = MY, 
                      assay = "RNA", slot = "data")["JCHAIN",]
LILRA4 = GetAssayData(object = MY, 
                      assay = "RNA", slot = "data")["LILRA4",]
names(which(CCR7 >0))
names(which(JCHAIN >0))
names(which(CD1C >0))
names(which(CLEC9A >0))
names(which(CLEC10A >0))
names(which(FCN1>0))
names(which(LILRA4>0))
DC <- names(which(LILRA4>0 | JCHAIN >0 | CD1C>0 | CLEC9A>0))
combined.all@meta.data[DC, 20] <- "DC"
mye <- MYELOID[!(MYELOID %in% DC)]
my <- subset(MY, cells= mye)
FCGR3A = GetAssayData(object = my, 
                      assay = "RNA", slot = "data")["FCGR3A",]
LILRB2 = GetAssayData(object = my, 
                      assay = "RNA", slot = "data")["LILRB2",]
names(which(LILRB2>0))
names(which(FCGR3A>0))
CDKN1C = GetAssayData(object = my, 
                      assay = "RNA", slot = "data")["CDKN1C",]
names(which(CDKN1C>0))
monoCD16 <- names(which(LILRB2>0 & FCGR3A>0 & CDKN1C>0))
combined.all@meta.data[monoCD16, 20] <- "CD16+ mono"
myNotCD16 <- MYELOID[!(MYELOID %in% monoCD16)]
notCD16 <- subset(MY, cells= myNotCD16)
FCN1 = GetAssayData(object = notCD16, 
                    assay = "RNA", slot = "data")["FCN1",]
S100A12 = GetAssayData(object = notCD16, 
                       assay = "RNA", slot = "data")["S100A12",]
CD14 = GetAssayData(object = notCD16, 
                    assay = "RNA", slot = "data")["CD14",]
monoCD14 <- names(which(CD14 >0 & S100A12 >1 & FCN1 >0))
combined.all@meta.data[monoCD14, 20] <- "CD14+ mono"
myREM <- myNotCD16[!(myNotCD16 %in% monoCD14)]
remainMY <- subset(notCD16, cells= myREM)
MARCO = GetAssayData(object = remainMY, 
                     assay = "RNA", slot = "data")["MARCO",]
macro <- names(which(MARCO>0))

# investigate this later
combined.all@meta.data[myREM, 20] <- "macrophages"

DimPlot(combined.all, cells.highlight = monoCD16)

FeaturePlot(combined.all, features = c("KIT", "SOX4", "TNFRSF18"))
# TNFRSF18 = nk & t-reg

FeaturePlot(combined.all, label=T, cols=c("yellow", "blue"), features = c("LYZ", "AIF1", "COTL1", "GNLY", "NKG7", 'MS4A1'))
#########################################################################################################################




# Normalize and SCale RNA slot
DefaultAssay(combined.all) <- "RNA"
NormalizeData(combined.all, assay = "RNA")
ScaleData(combined.all, assay = "RNA")


combined.all <- SetIdent(combined.all, value = "annotations" )
combined.all <- SetIdent(combined.all, value = "seurat_clusters" )
combined.all <- SetIdent(combined.all, value = "main" )
DimPlot(combined.all, label=TRUE)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
################################ 10/13
saveRDS(combined.all, "./combined_all_1013.rds")
combined.all <- readRDS( "./combined_all_1013.rds")

# subcluster main annotations
combined.all <- SetIdent(combined.all, value = "seurat_clusters" )
DimPlot(combined.all, label=TRUE)
head(combined.all@meta.data)

NKcells <- WhichCells(combined.all, idents= c("1", "7", "9"))
Myeloid <- WhichCells(combined.all, idents= c("0", "6"))
Bcells <- WhichCells(combined.all, idents= "2")
C11 <- WhichCells(combined.all, idents= "11")
Tcells <- WhichCells(combined.all, idents= c("3", "4", "5", "8", "10"))
Tcell <- subset(combined.all, cells=Tcells)
Bcell <- subset(combined.all, cells=Bcells)
myeloid <- subset(combined.all, cells= Myeloid)
NKcells <- subset(combined.all, cells=NKcells)

combined.all@meta.data[NKcells, "main"] <- "NK"
combined.all@meta.data[Myeloid, "main"] <- "Myeloid"
combined.all@meta.data[Bcells, "main"] <- "B Cells"
combined.all@meta.data[C11, "main"] <- "C11"
combined.all@meta.data[Tcells, "main"] <- "T Cells"


# t-cells
DefaultAssay(Tcell) <- "integrated"
Tcell <- RunPCA(Tcell, npcs = 30, verbose = FALSE)
Tcell <- RunUMAP(Tcell, dims = 1:20) 
Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:20)
Tcell <- FindClusters(Tcell, resolution = 1.0)
head(Tcell@meta.data)
Tcell <- SetIdent(Tcell, value = "integrated_snn_res.1")
DimPlot(Tcell, label=TRUE)
# 5 + TREG CD4
# 4 = CD8 EM/EMRA
# 2 = TRM/EM CD8
# 8 = CD4 EFFECTOR
# 3 = CD4 NAIVE
# 0,1,6 = CD8 NAIVE
# 7 ?
DoHeatmap(Tcell, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B", "CCR5", "CCR7", "LEF1", "SELL",
                            "CD8B", "CD69", "CD8A", "TRADD", "CXCR3", "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY", "GZMA"), assay = "RNA", slot="data", angle = 90) + NoLegend()

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


cd4t <- WhichCells(Tcell, idents = c("3", "5", "7", "8")) # 293 Cells
cd8t <- WhichCells(Tcell, idents = c("9", "2", "4", "0", "1","6")) # 1198 Cells
cd4 <- subset(Tcell, cells = cd4t)
cd8 <- subset(Tcell, cells = cd8t)
################################################### CD4
DefaultAssay(cd4) <- "integrated"
cd4 <- RunPCA(cd4, npcs = 30, verbose = FALSE)
cd4 <- RunUMAP(cd4, dims = 1:20) # reduction = "pca",
cd4 <- FindNeighbors(cd4, reduction = "pca", dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 0.8)
cd4<- SetIdent(cd4, value = "integrated_snn_res.0.8")
DimPlot(cd4, label=T, label.size = 8)


DefaultAssay(cd4) <- "RNA"
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "S100A4", "IL32", "USP10", "TIMP1", "CD27", "LYAR"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("TIMP1", "LGALS1", "GZMK", "CST7", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("FOXP3", "IL2RA", "IL1R1", "CD59"))
FeaturePlot(cd4, label=T, cols=c("yellow", "blue"), features = c("S100A4", "CCR7", "IL32", "ISG15"))
DoHeatmap(cd4, features = c("CCR7", "CD4", "CD40LG", "LEF1", "SELL", "USP10", "TIMP1", "LGALS1","CD27", "GZMK", "CST7",
                            "LYAR", "CCL5", "GZMA", "FOXP3", "IL2RA", "IL1R1", "CD59", "CD8A", "CD8B", "CCR5", "S100A4", "IL32", "ISG15"), assay = "RNA", slot="data", angle = 90) + NoLegend()

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


###########################################add annotatiions to combined.all
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
cd8<- SetIdent(cd8, value = "integrated_snn_res.0.8")
DimPlot(cd8, label=T, label.size = 6)

DefaultAssay(cd8) <- "RNA"
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3", "GZMA"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY"))
FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("FGFBP2", "GZMH", "GZMB", "GNLY"))

FeaturePlot(cd8, label=T, cols=c("yellow", "blue"), features = c("CD8A", "GZMK", "CCL5"))
DoHeatmap(cd8, features = c("CD4", "CD40LG","CCR7", "LEF1", "SELL", "CD8B", "CD69", "CD8A", "TRADD", "CXCR3",
                            "GZMA", "GZMK", "CCL5", "CD27", "FGFBP2", "GZMH", "GZMB", "GNLY"), assay = "RNA", slot="data", angle = 90) + NoLegend()


cd8naive <- WhichCells(cd8, idents = c("0", "4", "2"))
tm_em_temra <- WhichCells(cd8, idents = c("1", "3", '5'))
# cd8$annotate <- "NA"
cd8@meta.data[naive, "annotate"] <- "CD8+ Naive/CM"
cd8@meta.data[tm_em_temra, "annotate"] <- "CD8+ EM/TM/TEMRA"
cd8 <- SetIdent(cd8, value = "annotate")
DimPlot(Tcell, label=T, label.size = 6)

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

#############################################################################################################################
# October 14: Re-Run CD4 annotations and finish lung powerpoint (CD8)


CD8A = GetAssayData(object = Tcell, 
                     assay = "RNA", slot = "data")["CD8A",]
CD4 = GetAssayData(object = Tcell, 
                     assay = "RNA", slot = "data")["CD4",]

cd4 <- names(which(CD4>0))
cd8 <- names(which(CD8A>0))
DimPlot(Tcell, cells.highlight = cd8)
DefaultAssay(Tcell) <- "RNA"
tcell_markers <- FindAllMarkers(object = Tcell, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25) 

tcell_markers <- tcell_markers[ , c(6, 7, 1:5)]
tcell_markers <- tcell_markers %>%
  dplyr::arrange(cluster, p_val_adj)
top5_tcell <- tcell_markers %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)

# tcell_markers <- readRDS("./tcell_markers.rds")
DoHeatmap(Tcell, features = c('CD40LG', "CD4","CD8A", "CD8B", "CCR7", "IL7R", "SELL", "DTHD1", "LYST", "LEF1", "CCR4", "GZMH", "GZMK",
                              "GNLY", "NKG7", 'TNIP3', "CCL5", "FOXP3", "TNFRSF18", "CTLA4", "CXCR6", "CXCR3", 
                               "CCL4", "CCL3", "CMC1", "KLRC2", "FGFBP2", "PDCD1", "TIGIT", "HAVCR2", "LAG3", 
                               "PRDM1", "CD28", "PFN1", "CCR6", "IL2RA",  "IL2RB", "IFNG", "EOMES", "KLRF1",  "CRTAM"), assay = "RNA", slot = "data", angle = 90) + NoLegend()
DefaultAssay(Tcell) <- "integrated"
DefaultAssay(Tcell) <- "RNA"
FeaturePlot(Tcell, features = c("CD69", "HAVCR2", "LAG3", "KLRG1", "PRDM1", "RUNX3", "CTSG", "KIT", "TTN"), cols= c("yellow", "blue"), label=TRUE)
FeaturePlot(Tcell, features = c("CD8A", "CD4", "IL2RB"), cols= c("yellow", "blue"), label=TRUE)
FeaturePlot(Tcell, features = c("LEF1", "CD40LG", "LTB", "SELL", "CD28", "KLF2"), cols= c("yellow", "blue"), label=TRUE)
FeaturePlot(Tcell, features = c("CD8A", "CCR7"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("CD28", "LEF1", "SELL", "ZNF683"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("CD40LG", "CTLA4", "IL2RA", "LTB"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features=c("STAT1", "STAT4", "STAT3", "STAT6"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features=c("GATA3", "TBX21", "IRF4", "BATF"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features=c("CCR1", "CCR5", "CCR3", "CCR8"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features=c("CCR6", "GZMH"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("TNFRSF18", "FOXP3", "CCR4", "CTLA4", "CD4", "CD40LG", "IL7R", "IL2RA", "CD28", "LTB"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("FGFBP2", "NKG7", "CD44"), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(Tcell, features = c("CRTAM", "LYST", 'EOMES',"DTHD1", "TNIP3", 'GZMK', "GZMH", "KLRG1"), cols= c("yellow", "blue"), label=TRUE)


DimPlot(Tcell, cells.highlight = cd8, label=TRUE)
FeaturePlot(Tcell, features = c("ITGA1", "ITGAE", "GZMH", "KLRG1"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("GZMK", "GZMH", "DTHD1", "GNLY"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("TIGIT", "NKG7", "LYST", "TNIP3"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("GZMA", "GZMB", "IL7R", "CD8B"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c( "KLRB1", "CXCR3", "CCR5", "CXCR6"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c( "KLRF1", "FGFBP2", "EOMES", "KLFRC2", "ZNF683"), cols= c("yellow", "blue"), label = TRUE)
FeaturePlot(Tcell, features = c("CD44", "KLRG1", "CD57", "EOMES", "CCR7", "CD127", "CD62L", "CXCR3"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("SELL", "IL7R","IL2RA", "IL2RB", 'CD69', "OX40", "LAG3", "ICOS"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("BLIMP1", "GITRL", "FOXP3"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("IFNG", "IL10", "IL2", 'IL4', "GZMB", "GZMA", 'PRF1', 'CCL3', "CCL4"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("IFNG", "TNF", "IL12", "IL4", "IL5", 'IL13'), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("IL9", "TGFB1", "IL17", "IL21", "IL22", "IL25", "IL26", "IL23"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("IL10", "IL6", "IL1", "IL2", "IL21"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5
FeaturePlot(Tcell, features = c("CD69", "CXCR6", "CD49A","PD1", "IFNG", "IL2", "IL17A", "TNF"), cols= c("yellow", "blue"), label=TRUE) #CD4 TREG=5

FeaturePlot(Tcell, features = c("IL7R", "CCR7", "SELL", "CXCR3", "CD8A", "CD8B", "CD40LG"), cols= c("yellow", "blue"), label=TRUE) #CD8 NAIVE : 0/1/6
FeaturePlot(Tcell, features = c("CD27", "FAS", "ITGAL", "ITGB2", "ITGB1", 'IL2', "IFNG"), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(Tcell, features = c("PDCD1", "IL7R", "CD8A", "CD8B", 'IL32', "CD4", 'DTHD1'), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(Tcell, features = c("TIGIT", "TNFRSF4", "CCR7", 'IL7R', "CMC1", "GZMH", "KLRG1", 'NKG7'), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(Tcell, features = c("S1PR1", 'KLF2', "KLRF1", 'GNLY', "CCL4", "CCL5", "CD40LG", 'LEF1', "CD27"), cols= c("yellow", "blue"), label=TRUE)


marks2vs4 <- FindMarkers(Tcell, ident.1 = "2", ident.2 = "4")
# closer look at T-cell clusters 4,2, & 9
cd8_429 <- WhichCells(Tcell, idents = c("4", "2", "9"))
CD8429 <- subset(Tcell, cells = cd8_429)
DefaultAssay(CD8429) <- "RNA"
CD8429 <- RunPCA(CD8429, npcs = 30, verbose = FALSE)
CD8429 <- RunUMAP(CD8429, dims = 1:20) 
CD8429 <- FindNeighbors(CD8429, reduction = "pca", dims = 1:20)
CD8429 <- FindClusters(CD8429, resolution = 0.5)
head(CD8429@meta.data)
CD8429 <- SetIdent(CD8429, value = "integrated_snn_res.0.5")
DimPlot(CD8429, label=TRUE)
FeaturePlot(CD8429, features = c("GZMK", "TIGIT", "LYST", "GZMH", "NKG7", "GNLY", "DTHD1",
                                 "TNIP3", "EOMES", "CRTAM", "CX3CR1", "CD27", "IL7R", "KLRG1", "CMC1", "CTLA4"), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(CD8429, features = c("KLRG1", "CMC1", "IL7R", "TNIP3", "EOMES", "CRTAM", "CX3CR1", "CD27"), cols= c("yellow", "blue"), label=TRUE) 
FeaturePlot(CD8429, features = c("PDCD1", "ZNF683", "IL7R"), cols= c("yellow", "blue"), label=TRUE) 
DoHeatmap(CD8429, features = c("GZMK", "KLRG1", "TIGIT", "LYST", "GZMH", "NKG7", "GNLY", "DTHD1", 
                               "TNF", "TBX21", "GATA3", "IRF4", "KLRB1", "CCL5", "CXCR3", 
                               "FASLG", "CCL4", "CD27", "IL7R", "CCR7", "CTLA4"), assay = "RNA", slot="data", angle = 90) + NoLegend()




# MYELOID
DefaultAssay(myeloid) <- "integrated"
DefaultAssay(myeloid) <- "RNA"
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:20) 
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.6)
head(myeloid@meta.data)
myeloid <- SetIdent(myeloid, value = "integrated_snn_res.0.6")
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


myeloid9_markers <- FindMarkers(object = myeloid, ident.1 = "9")
write.csv(myeloid9_markers, "./nine_marks.csv")

myeloid_markers <- myeloid_markers[ , c(6, 7, 1:5)]
myeloid_markers <- myeloid_markers %>%
  dplyr::arrange(cluster, p_val_adj)

cd14 <- WhichCells(myeloid, idents = c("1", "2"))
cd16 <- WhichCells(myeloid, idents=c("3", "6"))
dc <- WhichCells(myeloid, idents=c("8"))
mdc <- WhichCells(myeloid, idents=c("5"))
imono <- WhichCells(myeloid, idents=c("0", "4", "7"))
neut <- WhichCells(myeloid, idents=c("9"))
combined.all@meta.data[cd14, "annotations"] <- "CD14 MONO"
combined.all@meta.data[cd16, "annotations"] <- "CD16 MONO"
combined.all@meta.data[dc, "annotations"] <- "cDC"
combined.all@meta.data[imono, "annotations"] <- "CD14+/16+ MONO"
combined.all@meta.data[mdc, "annotations"] <- "CD16 MONO"
combined.all@meta.data[neut, "annotations"] <- "CD16 MONO"

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





# NK CELLS
NKcells <- RunPCA(NKcells, npcs = 30, verbose = FALSE)
NKcells <- RunUMAP(NKcells, dims = 1:20) 
NKcells <- FindNeighbors(NKcells, reduction = "pca", dims = 1:20)
NKcells <- FindClusters(NKcells, resolution = 0.5)
head(NKcells@meta.data)
NKcells <- SetIdent(NKcells, value = "integrated_snn_res.0.5")
DimPlot(NKcells, label=TRUE)

NKcells_markers <- FindAllMarkers(object = NKcells, 
                                  only.pos = TRUE,
                                  logfc.threshold = 0.25) 

NKcells_markers <- NKcells_markers[ , c(6, 7, 1:5)]
NKcells_markers <- NKcells_markers %>%
  dplyr::arrange(cluster, p_val_adj)

FeaturePlot(NKcells, features = c("CD40LG", "CCR6", "SLC4A10", "GNLY", "NKG7", "FGFBP2", "ITGAD", "KLRC3"), cols= c("yellow", "blue"), label = TRUE)



# nk & b
nk<- WhichCells(combined.all, idents = c("1", "7", "9"))
bcell <- WhichCells(combined.all, idents=c("2"))
combined.all@meta.data[nk, "annotations"] <- "NK"
combined.all@meta.data[bcell, "annotations"] <- "B Cells"


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
combined.all <- SetIdent(combined.all, value = "seurat_clusters")
combined.all <- SetIdent(combined.all, value = "combined.main")

DimPlot(combined.all, reduction = "tsne", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap",label = TRUE,repel = TRUE)
FeaturePlot(combined.all, "CD8A") 





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
Idents(combined.all) <- "annotations"
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", label = TRUE,repel = TRUE)

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


cd16_DEgenes <- FindMarkers(combined.all, ident.1 = "CD16 MONO_Infected", ident.2 = "CD16 MONO_Not", verbose = FALSE)

# Plot Average Expression
cd16M <- subset(combined.all, idents = "CD16 MONO")
Idents(cd16M) <- "sample"
avg.cd16M <- as.data.frame(log10(AverageExpression(cd16M, verbose = FALSE)$RNA))
avg.cd16M$gene <- rownames(avg.cd16M)
genes.to.label = c("IFI27", "HLA-DRA", "B2M", "HLA-DPA1", "RBM47")
p1 <- ggplot(avg.cd16M, aes(Infected, Not)) + geom_point() + ggtitle("CD16 MONO")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


cd4toreg_DEgenes <- FindMarkers(combined.all, ident.1 = "CD4+ Treg_Infected", ident.2 = "CD4+ Treg_Not", verbose = FALSE)

# Plot Average Expression
cd4Treg <- subset(combined.all, idents = "CD4+ Treg")
Idents(cd4Treg) <- "sample"
avg.cd4Treg <- as.data.frame(log10(AverageExpression(cd4Treg, verbose = FALSE)$RNA))
avg.cd4Treg$gene <- rownames(avg.cd4Treg)
genes.to.label = c("IFI27", "CRIP1", "CCL4", "MT-CO1", "MT-ND4", "GNLY")
p1 <- ggplot(avg.cd4Treg, aes(Infected, Not)) + geom_point() + ggtitle("CD4 TREG")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

cd8em_DEgenes <- FindMarkers(combined.all, ident.1 = "CD8+ Eff/Mem_Infected", ident.2 = "CD8+ Eff/Mem_Not", verbose = FALSE)

# Plot Average Expression
cd8em <- subset(combined.all, idents = "CD8+ Eff/Mem")
Idents(cd8em) <- "sample"
avg.cd8em <- as.data.frame(log10(AverageExpression(cd8em, verbose = FALSE)$RNA))
avg.cd8em$gene <- rownames(avg.cd8em)
genes.to.label = c("IFI27", "GZMB", "CCL4", "GNLY")
p1 <- ggplot(avg.cd8em, aes(Infected, Not)) + geom_point() + ggtitle("CD8+ Eff/Mem")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

nk_DEgenes <- FindMarkers(combined.all, ident.1 = "NK_Infected", ident.2 = "NK_Not", verbose = FALSE)

# Plot Average Expression
NK <- subset(combined.all, idents = "NK")
Idents(NK) <- "sample"
avg.NK <- as.data.frame(log10(AverageExpression(NK, verbose = FALSE)$RNA))
avg.NK$gene <- rownames(avg.NK)
genes.to.label = c("IFI27", "GNLY", "CCL4", "CSF2", "CCL4L2")
p1 <- ggplot(avg.NK, aes(Infected, Not)) + geom_point() + ggtitle("NK Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)

#CD11C
FeaturePlot(combined.all, label = T, cols=c("yellow","blue"), features = c("ITGAX"))
AverageExpression(combined.all, features = "ITGAX", slot="data", assays = "RNA")
DoHeatmap(combined.all, features = c("ITGAX"), assay = "RNA", slot="data", angle = 90) + NoLegend()




saveRDS(combined.all, "./combined_all_origANNOTATE.rds")
saveRDS(combined.all, "./combined_all_new.rds")




