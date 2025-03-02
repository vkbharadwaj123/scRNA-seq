### adapted from Bioinformagician Youtube channel

# load libraries
library(Seurat)
library(tidyverse)

# Load the non small cell lung cancer dataset
nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m) # prints the different modalities
cts <-  nsclc.sparse.m$`Gene Expression` # we want just the gene expression


# Initialize the Seurat object with raw data
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200) # keep all cells that have 200 genes
str(nsclc.seurat.obj) 
nsclc.seurat.obj


# pre processing
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

# violin plot 
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# feature scatter plot -- 
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# keeping the cells with < 5% mitochondrial RNA, < 200 < 2500 genes
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                          percent.mt < 5)

# normalizing data
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj) # uses log normalization
str(nsclc.seurat.obj)


#Identify highly variable features 
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# scaling -- needed for PCA
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# principal component analysis = linear dimensional reduction
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# determine # of PCs to use via elbow plot
ElbowPlot(nsclc.seurat.obj)

# clustering using Shared Nearest-Neighbors clustering
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15) # using 15 for the sake of this analysis 

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
Idents(nsclc.seurat.obj)

# UMAP = non-linear dimensional reduction
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")

