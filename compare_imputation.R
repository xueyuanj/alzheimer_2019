library(dplyr)
library(Linnorm)
library(scater)
library(Seurat)
library(SingleCellExperiment)
library(stringr)

# strct filtering, ~1900 genes left
 
trem2_imputed = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.imputed.rds")
trem2_unimputed = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.unimputed.rds")

head(as.matrix(GetAssayData(trem2_unimputed, slot = "counts"))[,1:4])
head(as.matrix(GetAssayData(trem2_imputed, slot = "counts"))[,1:4])

DimPlot(trem2_unimputed, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(trem2_imputed, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# gene markers for imputed data
FeaturePlot(trem2_unimputed, features = c("Cx3cr1", "Hexb",  "Trem2", "Cst3"))
FeaturePlot(trem2_unimputed, features = c("B2m", "Apoe", "Fth1", "Tyrobp"))
FeaturePlot(trem2_unimputed, features = c("Lpl", "Cst7", "Apoe"))

# gene markers for un-imputed data
FeaturePlot(trem2_imputed, features = c("Cx3cr1", "Hexb",  "Trem2", "Cst3"))
FeaturePlot(trem2_imputed, features = c("B2m", "Apoe", "Fth1", "Tyrobp"))
FeaturePlot(trem2_imputed, features = c("Lpl", "Cst7", "Apoe"))




