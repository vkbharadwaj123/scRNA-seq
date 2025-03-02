# practicing working with multiple conditions (control vs treatment)
### adapted from DIY Transcriptomics

library(tidyverse)
library(Seurat)

load("spleen.naive.seurat")
DimPlot(spleen.naive.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

load("spleen.toxoInfected.seurat")
DimPlot(spleen.toxoInfected.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

targets <- read_tsv("studyDesign_singleCell.txt")

sampleID <- targets$sampleID
treatment <- targets$treatment

# annote seurat objects
spleen.naive.seurat$treatment <- treatment[1]
spleen.toxoInfected.seurat$treatment <- treatment[2]

spleen.toxoInfected.seurat@meta.data$treatment

spleen_features <- SelectIntegrationFeatures(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat))
spleen_anchors <- FindIntegrationAnchors(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat), 
                                         anchor.features = spleen_features)
spleen_integrated <- IntegrateData(anchorset = spleen_anchors)

spleen_integrated <- ScaleData(spleen_integrated, verbose = FALSE)
spleen_integrated <- RunPCA(spleen_integrated, npcs = 30, verbose = FALSE)
spleen_integrated <- RunUMAP(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindNeighbors(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindClusters(spleen_integrated, resolution = 0.5)
DimPlot(spleen_integrated, reduction = "umap", label = TRUE)

# splitting the clusters in the plot
DimPlot(spleen_integrated, reduction = "umap", split.by = "treatment",
        group.by = "seurat_clusters", label = TRUE)


# plotting certain genes on the umap
FeaturePlot(spleen_integrated, reduction = "umap", features = c("Cd4", "Cd8a"), # some immune markers
            pt.size = 0.4, order = TRUE, split.by = "treatment", min.cutoff = "q10", label = FALSE)


# using SingleR to identify the clusters
spleen_integrated.sce = as.SingleCellExperiment(spleen_integrated)

ref <- MouseRNAseqData(ensembl = FALSE)

predictions = SingleR(test = spleen_integrated.sce, assay.type.test = 1, ref = ref, labels = ref$label.main)

spleen_integrated.sce[["labels"]] <- predictions$labels

spleen_integrated.sce <- as.Seurat(spleen_integrated.sce, counts = NULL)
DimPlot(spleen_integrated.sce, reduction = "UMAP", split.by = "treatment", group.by = "labels", label = TRUE)



