library(dplyr)
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

trem2_files = unique(droplevels(trem2_mouse)$Amp_batch_ID)

# read the data and merge all the Trem2 ones
tremdf = lapply(as.character(trem2_files), function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

trem2megatable = data.frame(tremdf)

# keep a copy of the Trem2 data
saveRDS(trem2megatable, file = "~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.merged.rds")

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
umi_lower = trem2_qc$total_counts_endogenous > 200
#umi_upper = trem2_qc$total_counts_endogenous < 2500

# keep cells with expressed genes > 300
#feature_lower = trem2_qc$total_features_by_counts > 300

trem2_filter_by_cell = trem2_qc[,umi_lower]
# 1978 cells left

filter_by_cell_counts = counts(trem2_filter_by_cell)
gene_lower = rowSums(filter_by_cell_counts) > (0.005* ncol(filter_by_cell_counts))
#table(gene_lower)
# 10925 genes left

trem2_filter_by_cell_by_gene = trem2_filter_by_cell[gene_lower,]

# save a copy of the clean Trem2 data
saveRDS(trem2_filter_by_cell_by_gene, file = "~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.clean.1.rds")

# get the saved data
trem2_clean = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.umi.clean.1.rds")

trem2_seurat = CreateSeuratObject(counts = counts(trem2_clean), project = "trem2", min.cells = 3, min.features = 200)

VlnPlot(trem2_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
cor_plot2 = FeatureScatter(trem2_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(cor_plot2))

# normalization
trem2_seurat_filter = NormalizeData(trem2_seurat)

# select the 500 most variable features 
trem2_seurat_filter = FindVariableFeatures(trem2_seurat_filter, selection.method = "disp", nfeatures = 500)

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
trem2_seurat_filter = RunTSNE(trem2_seurat_filter, features = VariableFeatures(object = trem2_seurat_filter))

# tSNE map without annotation
DimPlot(trem2_seurat_filter, reduction = "tsne")
# homeostatic microglia
FeaturePlot(trem2_seurat_filter, features = c("Cx3cr1", "Hexb",  "Trem2", "Cst3"), reduction = "tsne")
# stage 2 microglia
FeaturePlot(trem2_seurat_filter, features = c("B2m", "Apoe", "Fth1", "Tyrobp"), reduction = "tsne")
# stage 3 microglia
FeaturePlot(trem2_seurat_filter, features = c("Lpl", "Cst7", "Apoe"), reduction = "tsne")

# UMAP map without annotation
DimPlot(trem2_seurat_filter, reduction = "umap")

# homeostatic microglia
FeaturePlot(trem2_seurat_filter, features = c("Cx3cr1", "Hexb",  "Trem2", "Cst3"))
# stage 2 microglia
FeaturePlot(trem2_seurat_filter, features = c("B2m", "Apoe", "Fth1", "Tyrobp"))
# stage 3 microglia
FeaturePlot(trem2_seurat_filter, features = c("Lpl", "Cst7", "Apoe"))

new_cluster_ids = c("homeostatic microglia", "stage 1 DAM", "stage 2 DAM", "unknown 1", "unknown 2")

names(new_cluster_ids) = levels(trem2_seurat_filter)
trem2_seurat_filter = RenameIdents(trem2_seurat_filter, new_cluster_ids)
DimPlot(trem2_seurat_filter, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(trem2_seurat_filter, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

# plot the source of cells
ad = dplyr::filter(clean_trem2_matrix, Mouse_ID == "5XFAD")%>%select(Well_ID)%>%unlist()%>%as.vector()
trem2ad=dplyr::filter(clean_trem2_matrix, Mouse_ID == "Trem2KO_5XFAD")%>%select(Well_ID)%>%unlist()%>%as.vector()
wt = dplyr::filter(clean_trem2_matrix, Mouse_ID == "C57BL/6")%>%select(Well_ID)%>%unlist()%>%as.vector()
trem2ko= dplyr::filter(clean_trem2_matrix, Mouse_ID == "Trem2KO")%>%select(Well_ID)%>%unlist()%>%as.vector()


UMAPPlot(object = trem2_seurat_filter, cells.highlight = ad)
UMAPPlot(object = trem2_seurat_filter, cells.highlight = trem2ad)
UMAPPlot(object = trem2_seurat_filter, cells.highlight = wt)
UMAPPlot(object = trem2_seurat_filter, cells.highlight = trem2ko)


TSNEPlot(object = trem2_seurat_filter, cells.highlight = ad)
TSNEPlot(object = trem2_seurat_filter, cells.highlight = trem2ad)
TSNEPlot(object = trem2_seurat_filter, cells.highlight = wt)
TSNEPlot(object = trem2_seurat_filter, cells.highlight = trem2ko)

# extract the data, then plot genes by cell design batch
clean_counts = as.data.frame(as.matrix(GetAssayData(trem2_seurat_filter, slot = "data")))

ad_set = t(clean_counts[, names(clean_counts)%in%ad])
trem2ad_set = t(clean_counts[, names(clean_counts)%in%trem2ad])
wt_set = t(clean_counts[, names(clean_counts)%in%wt])
trem2ko_set = t(clean_counts[, names(clean_counts)%in%trem2ko])

par(mar=c(5.1,4.1,4.1,2.1))
par(mai=c(1.02,0.82,0.82,0.42))

plotbox = function(x){
  boxplot(wt_set[,x],trem2ko_set[,x], trem2ad_set[, x], ad_set[,x], notch=T, 
          outline=F, names= c("WT", "Trem2-/-", "Trem2-/-AD", "AD"), col= c("yellow", "yellow", "orange", "red"),
          boxlwd = 2,main = x)}

genemarkers = c("Trem2", "Hexb", "Cx3cr1", "B2m", "Apoe", "Fth1", "Cst7", "Lpl")

sapply(genemarkers, plotbox)

# save a copy of this Seurat object
saveRDS(trem2_seurat_filter, file = "~/Documents/2019_summer_eisai/processed_data/dam.trem2.seurat.optimal.rds")

# write down the marker genes of homeostatic microglia, stage 1 DAM and stage 2 DAM

# up in homeo to stage 1
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = "stage 1 DAM", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in homeo to stage 2
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = "stage 2 DAM", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in hemeo to stage 1 &2
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = c("stage 1 DAM", "stage 2 DAM"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_12.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to homeo
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = "homeostatic microglia", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_homeo.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to stage 2
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = "stage 2 DAM", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to homeo and stage 2
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = c("homeostatic microglia", "stage 2 DAM"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_homeo2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to homeo
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = "homeostatic microglia", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to stage 1
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = "stage 1 DAM", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to homeo and stage 1
write.table(FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = c("homeostatic microglia", "stage 1 DAM"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

trem2_try = subset(trem2_seurat_filter, idents = c("homeostatic microglia", "stage 1 DAM", "stage 2 DAM"))
DimPlot(trem2_try, reduction = "umap",label = TRUE, pt.size = 0.5)+NoLegend()

