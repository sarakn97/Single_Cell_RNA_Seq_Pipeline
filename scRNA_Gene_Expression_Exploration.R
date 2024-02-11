######### SHAN LAB scRNA #######################################################
#### Date : August 2023
#### Sara Nicholson
#### Description: File for Shan Lab Members to investigate & visualize Gene Expression 
# in our 4 Tissue Samples (Bone Marrow, Spleen, Blood, Lung) with HIV-Infected & Control Samples
################################################################################

# If you want to make many changes to the script, you can copy it and change the name 
# so you are free to make any changes you would like. Small Changes like changing
# gene names dont require creating a copy file


##### INITIALIZATIONS -Always run these lines 10-15
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(edgeR)

# Here you will load into your working environment the data that you would like to analyze. 
# You highlight the chunk of code under the tissue that you would like to analyze and select "Run"
# Under the environment tab, you should now see "combined.all" as it is now loaded into the workspace
# The "combined.all" object has been pre-processed including alignment to human and mouse genome,
# mouse cells removed, and genes per cell quantified and normalized.

#-------------------------------Lung-------------------------------------------#
# set working directory
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Lung_sara") # Set path to the data-file
# saveRDS(combined.all, "./combined_all_ShanLab.rds") # only run this to save new resolutions or annotations
combined.all <- readRDS("./combined_all_ShanLab.rds")

#------------------------------Spleen------------------------------------------#
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Spleen_sara/")
# saveRDS(combined.all, "./combined_all_ShanLab_Spleen.rds")
combined.all <- readRDS("./combined_all_ShanLab_Spleen.rds")

#------------------------------Blood-------------------------------------------#
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_Blood")
# saveRDS(combined.all, "./combined_all_ShanLab_Blood.rds")
combined.all <- readRDS("./combined_all_ShanLab_Blood.rds")

#----------------------------Bone-Marrow---------------------------------------#
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/")
# saveRDS(combined.all, "./combined_all_ShanLab_BM.rds")
combined.all <- readRDS("./combined_all_ShanLab_BM.rds")

# This shows the genes and counts per cell stored in "combined.all", your seurat object.
# Genes are rows and cell names are columns
head(combined.all@assays$RNA@counts)
# simply running the name of the seurat object will show you the total gene(feature) counts and number of cells (samples)
combined.all 
# you can view the meta.data columns below, this data contains information for each individual cell such as what type of cell it is, 
# what cluster it belongs to and from with sample - infected or control - it is derived. You can always add new columns such as a new 
# set of annotations for the cells to the meta data
head(combined.all@meta.data)

############################# Run UMAP ##########################################
# DimPlot() command shows the UMAP dimensional reduction plot for the object
DimPlot(combined.all, label=T, label.size =3, pt.size = 1)
# Split UMAP by Infection Status
DimPlot(combined.all, label=F, label.size =3, pt.size = 1, split.by="sample")

# Pre-defined Annotations by Sara
colnames(combined.all@meta.data) # see identity options
combined.all <- SetIdent(combined.all, value = "main") # clusters annotated by cell-type 
combined.all <- SetIdent(combined.all, value = "more") # more specific annotations (Bone-Marrow does not have "more")

DimPlot(combined.all, label=T, label.size =5, pt.size = 1) # UMAP
table(Idents(combined.all), combined.all$sample) # Cell counts by Infection and Seurat Object Identity

################################ If you need a different resolution for clustering run the following:
# run all lines
DefaultAssay(combined.all) <- "integrated"
combined.all <- RunPCA(combined.all, npcs = 30, verbose = FALSE) # Run PCA dimensional Reduction
combined.all <- RunUMAP(combined.all, dims = 1:20)               # Run UMAP dimensional Reduction
combined.all <- FindNeighbors(combined.all, reduction = "pca", dims = 1:20) # Run KNN clustering 
combined.all <- FindClusters(combined.all, resolution = c(0.3, 0.5, 0.6, 0.8, 1.0))   # change # to new resolutions if desired, higher = more clusters
DefaultAssay(combined.all) <- "RNA"
combined.all <- ScaleData(combined.all)

colnames(combined.all@meta.data) # lists new identities to choose from
# "integrated.snn.res.##" are differnt clustering resolutions
combined.all <- SetIdent(combined.all, value = "integrated_snn_res.0.3")
DimPlot(combined.all, label=T, label.size =3, pt.size = 1)

################################################### Subset Combined.all to Infected or Non-Infected
# All cells are labeled as "Infected" or "Control" in the "sample" column of the meta.data table
table(combined.all$sample) # print cell counts for each category
# If you only want to explore the Non-Infected cells, you can subset your seurat object to only those cells and save it under a new name
# for example, I will subset only the NOn-INfected cells and save it under an object called "control"
# You must first set the seurat identity to "sample" so it knows which category under meta.data that you will use to subset the cells

combined.all <- SetIdent(combined.all, value = "sample") # set identity
control <- subset(combined.all, idents = "Control") # subset

# You now have a seurat object that only has the Control cells and all of the meta-data that belongs to those cells
# All of the same visualizations and functions can be performed on this subset, it just has a different name : "control"
# You will, however, need to re-run Dimensional Reduction and Clustering as shown below:
DefaultAssay(control) <- "integrated"
control <- RunPCA(control, npcs = 30, verbose = FALSE) # Run PCA dimensional Reduction
control <- RunUMAP(control, dims = 1:20)               # Run UMAP dimensional Reduction
control <- FindNeighbors(control, reduction = "pca", dims = 1:20) # Run KNN clustering 
control <- FindClusters(control, resolution = c(0.3, 0.5, 0.6, 0.8, 1.0))   # change # to new resolutions if desired, higher = more clusters
DefaultAssay(control) <- "RNA"
control <- ScaleData(control) # Re-Scale data

# You can re-establish the pre-established identities of the cells as the above command "FindClusters" 
# will automatically set the Identity to the new numbered clusters you created
combined.all <- SetIdent(combined.all, value = "more")
DimPlot(combined.all, label=T, label.size =3, pt.size = 1) # View UMAP

################################################################################# Find Cluster Markers
# This command will find markers of the clusters defined by the current identity of the Seurat object
marks <- FindAllMarkers(combined.all)
# write.csv(marks, "./seurat_lung_cluster_marks.csv") # This will write a csv of the markers that will be stored in the storage. 
# You can also go to the environment tab and click on marks to view it here in R


############################ SUBSET SEURAT OBJECT ############################################################################
# You can also subset to different cell types based on the annotations you provide, as we did before to subset to Control cells
combined.all <- SetIdent(combined.all, value = "main")
table(Idents(combined.all)) # View cells to choose from, will choose T Cells but can also create a list of cells to choose and use different annotations
Tcells <- subset(combined.all, idents = c("T Cells"))

# Re-Run dimensional reduction and clustering
DefaultAssay(Tcells) <- "integrated"
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE) # Run PCA dimensional Reduction
Tcells <- RunUMAP(Tcells, dims = 1:20)               # Run UMAP dimensional Reduction
Tcells <- FindNeighbors(Tcells, reduction = "pca", dims = 1:20) # Run KNN clustering 
Tcells <- FindClusters(Tcells, resolution = c(0.5, 0.6, 0.8, 1.0))   # change # to new resolutions if desired, higher = more clusters
DefaultAssay(Tcells) <- "RNA"
Tcells <- ScaleData(Tcells)

colnames(Tcells@meta.data) # lists new identities to choose from
Tcells <- SetIdent(Tcells, value = "more")
DimPlot(Tcells, label=T, label.size =5, pt.size = 1)

### Now you would put "Tcells" in all other functions instead of combined.all ###

##############################################################################################################################
############################ VISUALIZATIONS/EXPLORATION ######################################################################
##############################################################################################################################

############################ FEATURE PLOTS
# https://www.rdocumentation.org/packages/Seurat/versions/4.3.0.1/topics/FeaturePlot : link for options to change feature plots
FeaturePlot(Tcells, label=T, features = c("CD4", "CD8A"), cols = c("yellow", "red"))
FeaturePlot(combined.all, label=T, features = c("IL18"))
############################ Highlighting cells by expression
# change the genes as desired
IL18 = GetAssayData(object = combined.all, 
                    assay = "RNA", slot = "data")["IL18",]
CCR5 = GetAssayData(object = Tcells, 
                    assay = "RNA", slot = "data")["CCR5",]
bach2 <- GetAssayData(object = Tcells, 
                      assay = "RNA", slot = "data")["BACH2",]

# select cells that express that gene at any level (or at a certain level if you change 0 to a higher number)
il18 <- names(which(IL18>0)) 
ccr5 <- names(which(CCR5>0)) 
bach2 <- names(which(bach2>0)) 
length(cd4)
# can create selections for genes expressed together as shown below or use this symbol : "|" to signify "or"
# ">" greater than, "<" less than , "==" equals, "!=" not equals
bothCD4_CD8 <- names(which(bach2>0 & CD8A >0)) 

############# Plot highlights
# link to DimPLot() options: https://www.rdocumentation.org/packages/Seurat/versions/4.3.0.1/topics/DimPlot
DimPlot(combined.all, label = TRUE, cells.highlight = bothCD4_CD8, pt.size = 0.6, label.size = 6, sizes.highlight = 0.6)
DimPlot(combined.all, label = TRUE, cells.highlight = il18,  pt.size = 2, label.size = 6)
DimPlot(Tcells, label = T, cells.highlight = list(cd4, cd8a, ccr5), cols.highlight = c("red", "blue", "green"), cols = "grey",  pt.size = 0.6, label.size = 5, sizes.highlight = 0.6)

########### HEAT MAP
DoHeatmap(Tcells, features = c("CD3D", "CD3G", "CD4", "CD40LG", "CD8A", "CCR7", "SELL", "LEF1", "TCF7", "LGALS1",
                               "USP10", "ANXA1", "TIMP1", "GZMA", "GZMK", "CCL5", "LYAR", "FOXP3","MKI67"), assay = "RNA", slot="data", angle = 90) + NoLegend()


############ VIOLIN PLOT
colors = c('#E0E0E0', '#CC9966', '#DC143C', '#FF7F50', '#FFA500', '#FFFF00', '#9ACD32', '#6B8E23',
           '#98FB98', '#20B2AA', '#00FFFF','#1E90FF', '#87CEFA', '#0000CD', '#8A2BE2', '#7B68EE', '#BA55D3',
           '#DB7093', '#FF69B4', '#FFB6C1', '#FAEBD7', '#A0522D', '#BC8F8F', '#B0C4DE')

# genes <- 'MS4A1', "CD14", "FCGR3A", "IRF8", "KIT", "IL7R","NKG7", "NCAM1", "ANXA1", "LGALS1", "FOXP3", "GZMA", "GZMK", 'GZMB', "KLRG1", "CCR7", "SELL", "CD3G", "CD3D","CD4", "CD8A
VlnPlot(combined.all,c("IL18", "CD4"),
        pt.size=0,stack=T,flip=T,fill.by='ident',assay='RNA')+
  NoLegend()+
  scale_x_discrete(position='top')+
  scale_fill_manual(values=colors)+
  theme(panel.background=element_rect(color='black',size=0.8),
        panel.spacing=unit(0,'lines'),
        strip.text.y.right=element_text(size=12,face='italic',hjust=0),
        axis.text.x=element_text(size=12,angle=50,hjust=0),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=element_line(size=0.8),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=margin(2,2,2,2,'pt'),
        plot.title=element_blank())+
  coord_cartesian(clip='off')


################################ Average Expression
AverageExpression(combined.all, features = c("IL18"), assay = "RNA", slot = "data")


################################################################## PROPORTIONS PLOT
combined.all <- SetIdent(combined.all, value = "more") # Set Identity
prop = prop.table(table(Idents(combined.all),combined.all$sample),1) %>% data.frame()
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Control',],-prop)$cluster)
prop$cluster = factor(prop$cluster,levels=ranked)
counts = table(Idents(combined.all)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster,levels=ranked)

ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='Infected',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='Control',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,0.6))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('royalblue4','red3'),labels=c('Control','Infected'))+
  labs(y='Proportion (%)')+
  theme(panel.background=element_rect(),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=12,hjust=0.5,color='black'),
        axis.text.x=element_text(size=12,hjust=1,color='black',angle=45,vjust=1),
        legend.position='top',
        legend.text=element_text(size=16,color='black'),
        legend.title=element_blank(),
        axis.title.y=element_text(size=16,hjust=0.2,color='black',margin=margin(r=10)),
        axis.title.x=element_blank())


##### HOW TO CREATE NEW ANNOTATIONS

# For example, if you had a subset of T-Cells and then clustered them into 6 clusters, 
# you could annotate them as follows:

cd8_naive <- WhichCells(Tcells, idents = c(0,1))
cd4_naive <- WhichCells(Tcells, idents = c(2))
treg <- WhichCells(Tcells, idents = c(4))
cd8temrm <- WhichCells(Tcells, idents = c(5))
cd8temra <- WhichCells(Tcells, idents <- c(3,6))

combined.all@meta.data[cd8_naive, "more"] <- "CD8+ Naive/CM"
combined.all@meta.data[cd4_naive, "more"] <- "CD4+ Naive/CM"
combined.all@meta.data[treg, "more"] <- "CD4+ Treg"
combined.all@meta.data[cd8temrm, "more"] <- "CD8+ Tem/rm"
combined.all@meta.data[cd8temra, "more"] <- "CD8+ Tem/emra"

############################### Differential Expression ##########################
# non-parametric Wilcoxon rank sum test, adj.p.val is Bonferroni corrected
combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim" 
table(combined.all@meta.data$celltype.stim)

DEgenes <- FindMarkers(combined.all, ident.1 = "Treg_Infected", ident.2 = "Treg_Control", verbose = FALSE)
head(DEgenes, n = 15) # show first 15 DEGs

########################## Differential Expression Analysis with edgeR
# Subset Seurat Object to the cluster you would like to perform DE-Seq between Infected and Control
combined.all <- SetIdent(combined.all, value = "main") 
DimPlot(combined.all, label = T)

tcells = subset(combined.all,idents='T Cells') # fill in Identity of your choice
Idents(tcells) = tcells$sample
DefaultAssay(tcells) = 'RNA'
tcells = NormalizeData(tcells) %>% ScaleData()
mtx = as.matrix(GetAssayData(tcells,slot='counts'))
filter = rowSums(mtx>0) > 10 # filter to genes with atleast 10 counts across cells
dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=tcells$sample) # create edgeR DGElist object
tcells$sample = factor(tcells$sample,levels=c('Control','Infected'))
design = model.matrix(~0 + tcells$sample)
colnames(design) = c('Mock','Infected')
contrast.matrix = makeContrasts(Infected-Mock,levels=design)
dge = calcNormFactors(dge,method='TMM') # TMM normalize
dge = estimateDisp(dge,design=design)

# Negative Binomial Generalized Linear Model
fit = glmQLFit(dge,design=design)
res = glmQLFTest(fit,contrast=contrast.matrix)
pval = res$table
pval$fdr = p.adjust(pval$PValue,method='fdr')
pval = arrange(pval,fdr)

# write results to csv, change name of file (can also change file path) or simply view in RStudio
write.csv(pval, "./tcells_pval.csv")

#### import results
tcells_res <- read.csv("./tcells_pval.csv")  
head(tcells_res)
tcells_res_significant <- tcells_res %>% filter(fdr < 0.05 & logFC > 1) # filter to significant only at indicated levels
dim(tcells_res_significant) # 74 significant

##################### VOLCANO PLOTS
data <- tcells_res # change to desired subset 

diffexpr <- which(data$fdr <= 0.05)
data$FDR_Significance <- "Not Significant"

# Remove cases with NA Values
de <- data[complete.cases(data), ]
# if log2Foldchange > 1 and pvalue < 0.05, set as up-regulated ## CHANGE AS NEEDED
de$FDR_Significance[de$logFC > 1 & de$fdr < 0.05] <- "Significantly Up-Regulated"
# if log2Foldchange < -1 and pvalue < 0.05, set as down-regulated
de$FDR_Significance[de$logFC < -1 & de$fdr < 0.05] <- "Significantly Down-Regulated"
# Labels for Significant
de$delabel <- NA
de$delabel[de$FDR_Significance != "Not Significant"] <- de$X[de$FDR_Significance != "Not Significant"]


library(ggrepel)
options(ggrepel.max.overlaps = Inf)

ggplot(data=de, aes(x=logFC, y=-log10(fdr), col=FDR_Significance, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size=3) + # change font size of labels
  scale_color_manual(values=c("black", "blue",  "red", "green",  "purple")) + # change colors
  geom_vline(xintercept=c(-1, 1), col="blue", linetype=2) + # change boundaries
  geom_hline(yintercept=-log10(0.05), col="blue", linetype=2) + 
  labs(y = "-log10(Adjusted.Pvalue)", x = "Log2 Fold Change") + 
  theme(axis.title.x = element_text(size=16, face="bold"), axis.title.y = element_text(size=16, face="bold"), legend.text = element_text(size=), legend.title = element_text(size=12)) + ylim(0,200) + xlim(-3.5, 4) 
# may need to change +xlim() + ylim()








