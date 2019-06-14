library(dplyr)
library(Linnorm)

setwd("~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW")

# get the design matrix
design_matrix = read.csv("GSE98969_experimental_design_f.txt",sep = "\t", header = T)
# fix the inconsistent format in Batch_desc column
design_matrix$Batch_desc = gsub(" ", "_", design_matrix$Batch_desc)
# check the number of cells per well
table(design_matrix$Number_of_cells)

# to re-create figure one, select the 6 six-month old AD and WT mouses (there is a mouse 4 in AD)
admouseinfor = filter (design_matrix, Mouse_ID=="5XFAD")%>%filter(grepl("^AD6m_mouse", Batch_desc))%>%filter(!grepl("mouse4", Batch_desc))
wtmouseinfor = filter (design_matrix, Mouse_ID=="C57BL/6")%>%filter(grepl("^WT6m_mouse", Batch_desc))

# get the files using sample information
adfiles = unique(droplevels(admouseinfor)$Amp_batch_ID)
wtfiles = unique(droplevels(wtmouseinfor)$Amp_batch_ID)

# read the data and merge them into one data frame
adtables = lapply(as.character(adfiles), function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

admegatable = data.frame(adtables)

wttables = lapply(as.character(wtfiles), function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

wtmegatable = data.frame(wttables)

ad_wt_table = data.frame(admegatable, wtmegatable)

# save a copy of this file
saveRDS(ad_wt_table, file = "./pbmc3k_final.rds")

ad_wt_table = readRDS("dam.ad_wt.umi.merged.rds")


# Seurat version

ad_wt_seurat = CreateSeuratObject(counts = ad_wt_table, project = "dam", min.cells = 1, min.features = 200)
ad_wt_seurat[["percent.ercc"]] = PercentageFeatureSet(ad_wt_seurat, pattern ="^ERCC-" )

VlnPlot(ad_wt_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.ercc"), ncol = 3)

cor_plot1 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "percent.ercc")
cor_plot2 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(cor_plot1, cor_plot2))

ad_wt_seurat_filter = subset(ad_wt_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

ad_wt_seurat_filter = NormalizeData(ad_wt_seurat_filter)

# select the 220 most variable features 
ad_wt_seurat_filter = FindVariableFeatures(ad_wt_seurat_filter, selection.method = "vst", nfeatures = 220)

seurat_top10 = head(VariableFeatures(ad_wt_seurat_filter), 10)

var_plot1 = VariableFeaturePlot(ad_wt_seurat_filter)
var_plot2 = LabelPoints(plot = var_plot1, points = seurat_top10, repel = T)
CombinePlots(plots = list(var_plot1, var_plot2))

# scale the data, perform dimension reduction
all.genes = rownames(ad_wt_seurat_filter)
ad_wt_seurat_filter = ScaleData(ad_wt_seurat_filter, features = all.genes)

ad_wt_seurat_filter = RunPCA(ad_wt_seurat_filter, features = VariableFeatures(object = ad_wt_seurat_filter))

DimPlot(ad_wt_seurat_filter, reduction = "pca")

DimHeatmap(ad_wt_seurat_filter, dims = 1, cells = 500, balanced = T)


ad_wt_seurat_filter = FindNeighbors(ad_wt_seurat_filter, dims = 1:10)
ad_wt_seurat_filter = FindClusters(ad_wt_seurat_filter, resolution = 0.5)
ad_wt_seurat_filter = RunUMAP(ad_wt_seurat_filter, dims = 1:10)
DimPlot(ad_wt_seurat_filter, reduction = "umap")




