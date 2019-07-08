library(dplyr)
library(GEOquery)
library(reshape)
library(scater) 
library(Seurat)
library(stringr)
library(tidyverse)

setwd("~/Documents/2019_summer_eisai/raw_data/GSE121654_RAW")

# get the data information
gse121654 = getGEO("GSE121654")

# extract the injured mouse information
design_infor = pData(phenoData(gse121654[[1]]))
injurmouse = dplyr::filter(design_infor, grepl("100", title))%>%dplyr::filter(grepl("W", title))

for (i in 1:ncol(injurmouse)) {
  print(unique(injurmouse[,i]))
}

injured_urls = select(injurmouse, supplementary_file_1)%>%unlist()%>%as.character()

# download the processed matrix
sapply(injured_urls, function(x) download.file(x, unlist(strsplit(x, "/"))[9] ))

# merge the data into one big table
inj_files = sapply(injured_urls, function(x) unlist(strsplit(x, "/"))[9], USE.NAMES = F ) 
injured_df = lapply(inj_files, function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

injured_merge = injured_df%>%reduce(left_join, by = "GENE")
injured_merge = na.omit(injured_merge)

# save a copy of merged dataset
saveRDS(injured_merge, file = "~/Documents/2019_summer_eisai/processed_data/irm.umi.merged.rds")

# read the saved file, and get to the analysis
injured_mg = readRDS("~/Documents/2019_summer_eisai/processed_data/irm.umi.merged.rds")

# check ERCC and MT genes
totalgenes = injured_mg$GENE

# no ERCC, 11 MT genes
mtgenes = totalgenes[grepl("^mt", totalgenes)]
#psgenes = totalgenes[grepl("\\.ps", totalgenes)]
#endogenes.t = totalgenes[!totalgenes %in% mtgenes ]
#endogenes = endogenes.t[!endogenes.t %in% psgenes]
endogenes = totalgenes[!totalgenes %in% mtgenes]

rownames(injured_mg) = totalgenes
injured_mg_input = injured_mg[, 2:ncol(injured_mg)]
# 13316 cells

injure_seurat = CreateSeuratObject(counts = injured_mg_input, assay = "RNA", project = "irm", min.cells = 3, min.features = 200)
injure_seurat[["percent.mt"]] = PercentageFeatureSet(injure_seurat, pattern ="^mt" )

VlnPlot(injure_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cor_plot1 = FeatureScatter(injure_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
cor_plot2 = FeatureScatter(injure_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(cor_plot1, cor_plot2))

# QC of the dataset
injure_seurat_filter = subset(injure_seurat, subset= nFeature_RNA > 650 & nFeature_RNA < 2500 & percent.mt < 10)
injure_seurat_filter = injure_seurat_filter[endogenes,]
# 10672 genes, 11586 cells

injure_seurat_filter = NormalizeData(injure_seurat_filter)

# feature selection
injure_seurat_filter = FindVariableFeatures(injure_seurat_filter, selection.method = "disp", nfeatures = 500)

injure_top20 = head(VariableFeatures(injure_seurat_filter), 20)
var_plot1 = VariableFeaturePlot(injure_seurat_filter)
var_plot2 = LabelPoints(plot = var_plot1, points = injure_top20, repel = T)
CombinePlots(plots = list(var_plot2))

# dimention reduction
injure_seurat_filter = ScaleData(injure_seurat_filter)

# PCA
injure_seurat_filter = RunPCA(injure_seurat_filter, features = VariableFeatures(object = injure_seurat_filter))
DimPlot(injure_seurat_filter, reduction = "pca")

# tSNE
injure_seurat_filter = FindNeighbors(injure_seurat_filter, dims = 1:15, k.param = 30)
injure_seurat_filter = FindClusters(injure_seurat_filter, resolution = 0.125)
injure_seurat_filter = RunUMAP(injure_seurat_filter, features = VariableFeatures(object = injure_seurat_filter))

DimPlot(injure_seurat_filter, reduction = "umap")

# save a copy of the Seurat object
saveRDS(injure_seurat_filter, file = "~/Documents/2019_summer_eisai/processed_data/irm.seurat.v1.rds")

# read the saved Seurat data
injure_seurat_filter = readRDS("~/Documents/2019_summer_eisai/processed_data/irm.seurat.v1.rds")

# try tSNE
injure_seurat_filter = RunTSNE(injure_seurat_filter, features = VariableFeatures(object = injure_seurat_filter))
DimPlot(injure_seurat_filter, reduction = "tsne")

# deal with batch effect
# split the data into two parts:
# 1) batch1 -- cluster 0, 2
# 2) batch2 -- cluster 1, 3, 4

irm_cluster_id = as.data.frame(injure_seurat_filter$seurat_clusters)
irm_cluster_id[,"cell"] = rownames(irm_cluster_id)
colnames(irm_cluster_id) = c("cluster", "cell")

batch_1 = dplyr::filter(irm_cluster_id, cluster =="0"|cluster == "2")%>%select(cell)%>%unlist()%>%as.vector()
batch_2 = dplyr::filter(irm_cluster_id, cluster =="1"|cluster == "3"|cluster =="4")%>%select(cell)%>%unlist()%>%as.vector()

# split the seurat object into batch_1 and batch_2
irm_batch1 = subset(injure_seurat_filter, cells = batch_1)
irm_batch2 = subset(injure_seurat_filter, cells = batch_2)

irm_batch1 = NormalizeData(irm_batch1) 
irm_batch2 = NormalizeData(irm_batch2)

irm_batch1 = FindVariableFeatures(irm_batch1, selection.method = "disp", nfeatures = 500)
irm_batch2 = FindVariableFeatures(irm_batch2, selection.method = "disp", nfeatures = 500)

irm_batch1 = ScaleData(irm_batch1)
irm_batch2 = ScaleData(irm_batch2)

irm_anchors = FindIntegrationAnchors(object.list = c(irm_batch1, irm_batch2), dims = 1:15)

irm_integrated = IntegrateData(anchorset = irm_anchors, dims = 1:15)

DefaultAssay(irm_integrated) = "integrated"

irm_integrated = ScaleData(irm_integrated)
irm_integrated = RunPCA(irm_integrated)
DimPlot(irm_integrated, reduction = "pca")


irm_integrated = FindNeighbors(irm_integrated, dims = 1:15, k.param = 30)
irm_integrated = FindClusters(irm_integrated , resolution = 0.125)
irm_integrated = RunTSNE(irm_integrated)
DimPlot(irm_integrated, reduction = "tsne")

irm_integrated = RunUMAP(irm_integrated, dims = 1:15)
DimPlot(irm_integrated, reduction = "umap")


# save a copy of batch-corrected data
saveRDS(irm_integrated, file = "~/Documents/2019_summer_eisai/processed_data/irm.seurat.batch_corrected.rds")



# plot cluster markers
# dam markers
FeaturePlot(irm_integrated, features = c("Apoe", "Cst7", "Lpl", "Spp1"), reduction = "tsne")

# ir1 markers
FeaturePlot(irm_integrated, features = c("P2ry12", "Selplg", "Cx3cr1"), reduction = "tsne")

# ir2 overall markers
FeaturePlot(irm_integrated, features = c("Ifi27l2a", "Ifi204", "Cxcl10", "Ccl4"), reduction = "tsne")

# get de genes for all 5 clusters
for (i in 0:4) {
  cluster.markers = FindMarkers(injure_seurat_filter, ident.1 = i, min.pct = 0.25, only.pos = T,logfc.threshold = 0.25)
  print(i)
  print(head(cluster.markers, n = 30))
}


# subset the data, highlight each condition on the plot
# creat vectors of well IDs that belong to each of the three conditions

controwb = dplyr::filter(injurmouse, grepl("contro", title))
salinewm = dplyr::filter(injurmouse, grepl("saline", title))
lpcwm = dplyr::filter(injurmouse, grepl("LPC", title))

# a wrapper that extract well IDs (not provided by the data table)
selectwell = function(x){
  injured_urls = select(x, supplementary_file_1)%>%unlist()%>%as.character()
  
  # merge the data into one big table
  inj_files = sapply(injured_urls, function(x) unlist(strsplit(x, "/"))[9], USE.NAMES = F ) 
  injured_df = lapply(inj_files, function(x){
    read.csv(list.files(pattern = x), sep="\t", header=TRUE)})
  
  injured_merge = injured_df%>%reduce(left_join, by = "GENE")
  injured_merge = na.omit(injured_merge)
  return(colnames(injured_merge)[2:ncol(injured_merge)])
  
}

contrlwb_wells = selectwell(controwb)
salinewm_wells = selectwell(salinewm)
lpcwm_wells = selectwell(lpcwm)


# first plot cell source 

DimPlot(object = irm_integrated, cells.highlight = contrlwb_wells, reduction = "tsne")
DimPlot(object = irm_integrated, cells.highlight = salinewm_wells, reduction = "tsne")
DimPlot(object = irm_integrated, cells.highlight = lpcwm_wells, reduction = "tsne")

DimPlot(object = irm_integrated, cells.highlight = contrlwb_wells, reduction = "umap")
DimPlot(object = irm_integrated, cells.highlight = salinewm_wells, reduction = "umap")
DimPlot(object = irm_integrated, cells.highlight = lpcwm_wells, reduction = "umap")

# select cluster 0 and 2, redo the plot
irm_try = subset(injure_seurat_filter, idents = c(0,2))
DimPlot(irm_try, reduction = "umap",label = TRUE, pt.size = 0.5)+NoLegend()




















