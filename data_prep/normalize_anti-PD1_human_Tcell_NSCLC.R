# this is for normalize data "anti-PD1_human_Tcell_NSCLC" in Lewis:
# conda activate r-mofa

library(Seurat)

setwd("/storage/htc/joshilab/Su_Li/Alg_development/scRNA_TCRBCR_surfaceProtein/data_collection/anti-PD1_human_Tcell_NSCLC/")

pbmc.data = readRDS("GSE179994_all.Tcell.rawCounts.rds")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

