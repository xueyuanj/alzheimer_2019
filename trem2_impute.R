library(dplyr)
library(DrImpute)
library(Linnorm)
library(scater)
library(Seurat)
library(SingleCellExperiment)
library(stringr)

setwd("~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW")

# get the design matrix
design_matrix = read.csv("GSE98969_experimental_design_f.txt",sep = "\t", header = T)
# fix the inconsistent format in Batch_desc column
design_matrix$Batch_desc = gsub(" ", "_", design_matrix$Batch_desc)

# get wells with cells
cell_well = dplyr::filter(design_matrix, Number_of_cells == "1")%>%select(Well_ID)%>%unlist()

# Trem2 mouse designs
trem2_mouse = filter(design_matrix, grepl("^Trem2", Batch_desc))

clean_trem2_matrix = trem2_mouse[(which(trem2_mouse$Well_ID%in%cell_well)),]


# start from saved file
trem2df = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.merged.rds")

# remove empty wells from this dataset
trem2df = trem2df[,(which(colnames(trem2df)%in%cell_well))]
# 4934 cells

# create the data
trem2_qc = SingleCellExperiment(assays = list(counts = as.matrix(trem2df)), colData = clean_trem2_matrix)

trem_total_genes = rownames(trem2_qc)

# check data quality
isSpike(trem2_qc, "ERCC") = grepl("^ERCC-", trem_total_genes)
isSpike(trem2_qc, "MT") = grepl("^mt-", trem_total_genes)

trem2_qc = calculateQCMetrics(trem2_qc,feature_controls = list(ERCC = isSpike(trem2_qc, "ERCC"), 
                                                               MT = isSpike(trem2_qc, "MT")))
plotColData(trem2_qc,x = "total_features_by_counts",y = "pct_counts_MT",colour = "Amp_batch_ID")
plotColData(trem2_qc, x = "total_features_by_counts",y = "pct_counts_ERCC",colour = "Amp_batch_ID")

# remove ERCC genes and mt genes
trem2_qc = trem2_qc[!(isSpike(trem2_qc, "ERCC")),]
trem2_qc = trem2_qc[!(isSpike(trem2_qc, "MT")),]

# remove the low-quality cells
# keep cells with 500 < UMI < 2500
umi_lower = trem2_qc$total_counts_endogenous > 500
umi_upper = trem2_qc$total_counts_endogenous < 2500

# keep cells with expressed genes > 300
feature_lower = trem2_qc$total_features_by_counts > 300

trem2_filter_by_cell = trem2_qc[,umi_lower&umi_upper&feature_lower]
# 1965 cells left

filter_by_cell_counts = counts(trem2_filter_by_cell)
gene_lower = rowSums(filter_by_cell_counts) > (0.005* ncol(filter_by_cell_counts))
#table(gene_lower)
# 10919 genes left

trem2_filter_by_cell_by_gene = trem2_filter_by_cell[gene_lower,]

# save a copy of the clean Trem2 data
#saveRDS(trem2_filter_by_cell_by_gene, file = "~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.clean.1.rds")

# get the saved data
#trem2_clean = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.clean.1.rds")

trem2_seurat = CreateSeuratObject(counts = counts(trem2_filter_by_cell_by_gene), project = "trem2", min.cells = 3, min.features = 200)


# normalization
trem2_seurat_filter = NormalizeData(trem2_seurat)

# impute with DrImpute
trem2_seurat_imputed = DrImpute(as.matrix(GetAssayData(trem2_seurat_filter, slot = "counts")))
colnames(trem2_seurat_imputed) = colnames(counts(trem2_filter_by_cell_by_gene))
trem2_imputed_input = CreateSeuratObject(counts = trem2_seurat_imputed, project = "trem2", min.cells = 3, min.features = 200)

# select the 500 most variable features 
trem2_seurat_filter = FindVariableFeatures(trem2_imputed_input, selection.method = "disp", nfeatures = 500)

trem2_seurat_top20 = head(VariableFeatures(trem2_seurat_filter), 20)

var_plot1 = VariableFeaturePlot(trem2_seurat_filter)
var_plot2 = LabelPoints(plot = var_plot1, points = trem2_seurat_top20, repel = T)
CombinePlots(plots = list(var_plot2))


# scale the data, perform dimension reduction
all.genes = rownames(trem2_seurat_filter)
trem2_seurat_filter = ScaleData(trem2_seurat_filter, features = all.genes)

trem2_seurat_filter = RunPCA(trem2_seurat_filter, features = VariableFeatures(object = trem2_seurat_filter))

DimPlot(trem2_seurat_filter, reduction = "pca")

#DimHeatmap(trem2_seurat_filter, dims = 1, cells = 500, balanced = T)

trem2_seurat_filter = FindNeighbors(trem2_seurat_filter, dims = 1:40, k.param = 30)
trem2_seurat_filter = FindClusters(trem2_seurat_filter, resolution = 0.3)
trem2_seurat_filter = RunUMAP(trem2_seurat_filter, features = VariableFeatures(object = trem2_seurat_filter))

# tSNE map without annotation
DimPlot(trem2_seurat_filter, reduction = "umap")

# highlight the marker genes
# microglia1
FeaturePlot(trem2_seurat_filter, features = c("Cx3cr1", "Hexb",  "Trem2", "Cst3"))

# microglia2
FeaturePlot(trem2_seurat_filter, features = c("B2m", "Apoe", "Fth1", "Tyrobp"))

# microglia3
FeaturePlot(trem2_seurat_filter, features = c("Lpl", "Cst7", "Apoe"))

# save a copy of un-imputed file
saveRDS(trem2_seurat_filter, file="~/Documents/2019_summer_eisai/processed_data/dam.trem2.unimputed.rds")

# save a copy of the imputed file
saveRDS(trem2_seurat_filter, file="~/Documents/2019_summer_eisai/processed_data/dam.trem2.imputed.rds")

#new_cluster_ids = c("homeostatic microglia", "stage 2 DAM", "stage 3 DAM", "unknown 1", "unknown 2")




