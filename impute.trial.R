library(dplyr)
library(DrImpute)
library(Linnorm)
library(SAVER)
library(scater)
library(SC3)
library(Seurat)
library(SingleCellExperiment)
library(stringr)

# set working directory
ad_wt_clean = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.ad_wt.umi.merged.clean2.rds")

# normalize with Linnorm
ad_wt_log = Linnorm(counts(ad_wt_clean))
ad_wt_linnorm = Linnorm.Norm(counts(ad_wt_clean), output = "XPM")
ad_wt_linnorm_log = log(ad_wt_linnorm+1)

# apply Seurat
ad_wt_seurat = CreateSeuratObject(counts = ad_wt_linnorm_log, project = "dam", min.cells = 3, min.features = 200)
#ad_wt_seurat[["percent.ercc"]] = PercentageFeatureSet(ad_wt_seurat, pattern ="^ERCC-" )

VlnPlot(ad_wt_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#cor_plot1 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "percent.ercc")
cor_plot2 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(cor_plot2))

# this step is no longer needed
ad_wt_seurat_filter = subset(ad_wt_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# default Seurat normalization
#ad_wt_seurat_filter = NormalizeData(ad_wt_seurat_filter)


# impute the drop-outs
# DrImpute
#ad_wt_imputed = DrImpute(as.matrix(GetAssayData(ad_wt_seurat_filter, slot = "counts")))

# SAVER
#ad_wt_saver = saver(as.matrix(GetAssayData(ad_wt_seurat_filter, slot = "counts")), ncores = 12, estimates.only = T)
 

# select the 500 most variable features 
ad_wt_seurat_filter = FindVariableFeatures(ad_wt_seurat_filter, selection.method = "disp", nfeatures = 500)

seurat_top20 = head(VariableFeatures(ad_wt_seurat_filter), 20)

var_plot1 = VariableFeaturePlot(ad_wt_seurat_filter)
var_plot2 = LabelPoints(plot = var_plot1, points = seurat_top20, repel = T)
CombinePlots(plots = list(var_plot2))

# scale the data, perform dimension reduction
all.genes = rownames(ad_wt_seurat_filter)
ad_wt_seurat_filter = ScaleData(ad_wt_seurat_filter, features = all.genes)

ad_wt_seurat_filter = RunPCA(ad_wt_seurat_filter, features = VariableFeatures(object = ad_wt_seurat_filter))

DimPlot(ad_wt_seurat_filter, reduction = "pca")

DimHeatmap(ad_wt_seurat_filter, dims = 1, cells = 500, balanced = T)

ad_wt_seurat_filter = FindNeighbors(ad_wt_seurat_filter, dims = 1:40, k.param = 30)
ad_wt_seurat_filter = FindClusters(ad_wt_seurat_filter, resolution = 0.3)
ad_wt_seurat_filter = RunUMAP(ad_wt_seurat_filter, features = VariableFeatures(object = ad_wt_seurat_filter))

DimPlot(ad_wt_seurat_filter, reduction = "umap")

# to identify clusters, first get the markders of each cluster

# microglial 1
#VlnPlot(ad_wt_seurat_filter, features = c("Cst3", "Hexb") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cst3", "Hexb", "Cst3", "Hexb"))

# microglial 2
# microglial 3
#VlnPlot(ad_wt_seurat_filter, features = c("Apoe", "Cst7", "Lpl", "Cd9", "Trem2") )
FeaturePlot(ad_wt_seurat_filter, features = c("Apoe", "Cst7", "Lpl", "Cd9"))

# perviascular Mf
#VlnPlot(ad_wt_seurat_filter, features = c("Cd74", "Cd163", "Mrc1") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cd74", "Cd163", "Mrc1"))

# Monocytes
#VlnPlot(ad_wt_seurat_filter, features = c("S100a4", "Cd74") )
FeaturePlot(ad_wt_seurat_filter, features = c("S100a4", "Cd74", "S100a4", "Cd74"))

# Mature B-cells and immature B-cells
#VlnPlot(ad_wt_seurat_filter, features = c("Cd74", "Cd79b", "Rag1") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cd74", "Cd79b", "Rag1"))

# T/NK cells
#VlnPlot(ad_wt_seurat_filter, features = c("Nkg7", "Trbc2") )
FeaturePlot(ad_wt_seurat_filter, features = c("Nkg7", "Trbc2", "Nkg7", "Trbc2"))

# Granulocytes 1 &2
#VlnPlot(ad_wt_seurat_filter, features =c("Camp", "S100a9")  )
FeaturePlot(ad_wt_seurat_filter, features = c("Camp", "S100a9", "Camp", "S100a9"))

new_cluster_ids = c("microglia1", "T/NK; monocyte", "perviascular Mf; B cell",  "microglia2/3",  "granulocytes")


names(new_cluster_ids) = levels(ad_wt_seurat_filter)
ad_wt_seurat_filter = RenameIdents(ad_wt_seurat_filter, new_cluster_ids)
DimPlot(ad_wt_seurat_filter, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# print out the top cluster markers 
for (i in 1:length(new_cluster_ids)) {
  cluster.markers = FindMarkers(ad_wt_seurat_filter, ident.1 = new_cluster_ids[i], min.pct = 0.25, only.pos = T,logfc.threshold = 0.25)
  print(new_cluster_ids[i])
  print(head(cluster.markers, n = 30))
}
