library(Seurat)
library(dplyr)

###################################################################
### 2700 PBMCs
###################################################################
pbmc.data = Read10X(data.dir = "~/Documents/2019_summer_eisai/filtered_gene_bc_matrices/hg19")
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

CombinePlots(plots=list(plot1, plot2))

pbmc = subset(pbmc, subset = nFeature_RNA>200 & nFeature_RNA <2500 & percent.mt < 5)

pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc = NormalizeData(pbmc)

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(pbmc), 10)

plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)

CombinePlots(plots= list(plot1, plot2) )

all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
  
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims=1:5, nfeatures = 5)  

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)  

pbmc = JackStraw(pbmc, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

pbmc = RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "./pbmc_tutorial.rds")

cluster1.markers = FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n= 5)

cluster5.markers = FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n=5)

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_logFC)

cluster1.markers = FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

top10 = pbmc.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene)

new.cluster.ids = c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) = levels(pbmc)
pbmc = RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.25) +NoLegend()

saveRDS(pbmc, file = "./pbmc3k_final.rds")

###################################################################
### simulated vs control PBMCs
###################################################################

library(Seurat)
library(cowplot)

ctrl.data <- read.table(file = "~/Documents/2019_summer_eisai/immune_alignment_expression_matrices/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "~/Documents/2019_summer_eisai/immune_alignment_expression_matrices/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

ctrl = CreateSeuratObject (counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim = "CTRL"

ctrl = subset(ctrl, subset = nFeature_RNA > 500)
ctrl = NormalizeData(ctrl, verbose = FALSE)
ctrl = FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

stim = CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim = "STIM"

stim = subset(stim, subset = nFeature_RNA >500)
stim = NormalizeData(stim, verbose = FALSE)
stim = FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

immune.anchors = FindIntegrationAnchors(object.list = list(ctrl, stim), dim= 1:20)
immune.combined = IntegrateData(anchorset = immune.anchors, dim=1:20)

DefaultAssay(immune.combined) = "integrated"

immune.combined = ScaleData(immune.combined, verbose = FALSE)
immune.combined = RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined = RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined = FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined = FindClusters(immune.combined, resolution = 0.5)

p1 = DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 = DimPlot(immune.combined, reduction = "umap", label = TRUE)

plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

DefaultAssay(immune.combined) = "RNA"
nk.markers = FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9")
immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined, label = TRUE)

Idents(immune.combined) = factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", 
                                                                     "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
                                                                     "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()

t.cells = subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) = "stim"
avg.t.cells = log1p(AverageExpression(t.cells, verbose = F)$RNA)
avg.t.cells$gene = rownames(avg.t.cells)

cd14.mono = subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono)= "stim"
avg.cd14.mono = log1p(AverageExpression(cd14.mono, verbose = F)$RNA)
avg.cd14.mono$gene = rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

p1 = ggplot(avg.t.cells, aes(CTRL, STIM))+
  geom_point()+ ggtitle("CD4 Naive T Cells")
p1 = LabelPoints(plot = p1, points = genes.to.label, repel = T)
p2 = ggplot(avg.cd14.mono, aes(CTRL, STIM))+
  geom_point()+ggtitle("CD14 Monocytes")
p2 = LabelPoints(plot = p2, points = genes.to.label, repel = T)
plot_grid(p1, p2)

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))
plots = VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)


###################################################################
### Clustering of the Microwell-seq Mouse Cell 
###################################################################

library(Seurat)

mca.matrix = readRDS(file = "~/Documents/2019_summer_eisai/MCA/MCA_merged_mat.rds")
mca.metadata = read.csv(file = "~/Documents/2019_summer_eisai/MCA/MCA_All-batch-removed-assignments.csv", row.names=1)

mca = CreateSeuratObject(counts = mca.matrix, meta.data = mca.metadata, project = "MouseCellAtlas")

mca = NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)

mca = FindVariableFeatures(mca)

mca[["percent.mt"]] = PercentageFeatureSet(mca, pattern ="^mt-")
mca = ScaleData(mca, vars.to.regress = "percent.mt")

##Error message: vector memory exhausted (limit reached?)


###################################################################
### Multi-modal analysis
###################################################################

library(Seurat)
library(ggplot2)

cbmc.rna = as.sparse(read.csv(file = "~/Documents/2019_summer_eisai/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",",                             header = TRUE, row.names = 1))
cbmc.rna = CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc.adt = as.sparse(read.csv(file = "~/Documents/2019_summer_eisai/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",                            header = TRUE, row.names = 1))
cbmc.adt = cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]

cbmc = CreateSeuratObject(counts = cbmc.rna)
cbmc = NormalizeData(cbmc)
cbmc = FindVariableFeatures(cbmc)
cbmc = ScaleData(cbmc)
cbmc = RunPCA(cbmc, verbose = F)
ElbowPlot(cbmc, ndims = 40)

cbmc = FindNeighbors(cbmc, dims = 1:25)
cbmc = FindClusters(cbmc, resolution = 0.8)
cbmc = RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")

cbmc.rna.markers = FindAllMarkers (cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = T)

new.cluster.ids = c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs")
names(new.cluster.ids) = levels(cbmc)
cbmc = RenameIdents(cbmc, new.cluster.ids)

DimPlot(cbmc, label = T) + NoLegend()

cbmc[["ADT"]] = CreateAssayObject(counts = cbmc.adt)
cbmc = NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc = ScaleData(cbmc, assay = "ADT")

FeaturePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16", "CD3E", "ITGAX", "CD8A", 
                               "FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

RidgePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16"), ncol = 2)

FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "CD3E")

tcells = subset(cbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")

saveRDS(cbmc, file = "~/Documents/2019_summer_eisai/cbmc.rds")

cbmc.small = subset(cbmc, downsample=300)
adt.markers = FindAllMarkers(cbmc.small, assay = "ADT", only.pos = T)
DoHeatmap(cbmc.small, features = unique(adt.markers$gene), assay = "ADT", angle = 90)+NoLegend()

cbmc = subset(cbmc, idents = c("Multiplets", "Mouse"), invert = T)

DefaultAssay(cbmc) = "ADT"
cbmc = RunPCA(cbmc, features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", verbose = F)
DimPlot(cbmc, reduction = "pca_adt")

adt.data = GetAssayData(cbmc, slot = "data")
adt.dist = dist(t(adt.data))

cbmc[["rnaClusterID"]] = Idents(cbmc)
cbmc[["tsne_adt"]] = RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cbmc[["adt_snn"]] = FindNeighbors(adt.dist)$snn
cbmc = FindClusters(cbmc, resolution = 0.2, graph.name = "adt_snn")

clustering.table = table(Idents(cbmc), cbmc$rnaClusterID)
clustering.table

new.cluster.ids = c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets", 
                    "CD16+ Mono", "pDCs", "B")
names(new.cluster.ids)= levels(cbmc)
cbmc = RenameIdents(cbmc, new.cluster.ids)

tsne_rnaClusters = DimPlot(cbmc, reduction = "tsne_adt", group.by = "rnaClusterID") + NoLegend()
tsne_rnaClusters = tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters = LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 4)

tsne_adtClusters = DimPlot(cbmc, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()
tsne_adtClusters = tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters = LabelClusters(plot = tsne_adtClusters, id = "ident", size = 4)

CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

tcells = subset(cbmc, idents = c("CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "CD4", feature2 = "CD8")
RidgePlot(cbmc, features = c("adt_CD11c", "adt_CD8", "adt_CD16", "adt_CD4", "adt_CD19", "adt_CD14"), 
          ncol = 2)

saveRDS(cbmc, file = "~/Documents/2019_summer_eisai/cbmc_multimodal.rds")



###################################################################
### Multiple dataset integration and label transfer
###################################################################
library(Seurat)

pancreas.data = readRDS(file = "~/Documents/2019_summer_eisai/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata = readRDS(file = "~/Documents/2019_summer_eisai/pancreas_v3_files/pancreas_metadata.rds")

pancreas = CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list = SplitObject(pancreas, split.by = "tech")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] = NormalizeData(pancreas.list[[i]], verbose = F)
  pancreas.list[[i]] = FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}


reference.list = pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors = FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated = IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)

DefaultAssay(pancreas.integrated) = "integrated"
pancreas.integrated = ScaleData(pancreas.integrated, verbose = F)
pancreas.integrated = RunPCA(pancreas.integrated, npcs = 30, verbose = F)
pancreas.integrated = RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)

p1 = DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech") 
p2 = DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", repel = T, label = T)+NoLegend()
plot_grid(p1, p2)

pancreas.query = pancreas.list[["fluidigmc1"]]
pancreas.anchors = FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, dims = 1:30)
predictions = TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, dims = 1:30)
pancreas.query = AddMetaData(pancreas.query, metadata = predictions)

pancreas.query$prediction.match = pancreas.query$predicted.id == pancreas.query$celltype

table(pancreas.query$prediction.match)

VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")



###################################################################
### sctransfer
###################################################################

library(Seurat)
library(ggplot2)
library(sctransform)

pbmc_data = Read10X(data.dir = "~/Documents/2019_summer_eisai/filtered_gene_bc_matrices/hg19")
pbmc = CreateSeuratObject(counts = pbmc_data)

pbmc = PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc = SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = F)

pbmc = RunPCA(pbmc, verbose = F)
pbmc = RunUMAP(pbmc, dims = 1:30, verbose = F)

pbmc = FindNeighbors(pbmc, dims = 1:30, verbose = F)
pbmc = FindClusters(pbmc, verbose = F)
DimPlot(pbmc, label = T)+NoLegend()


VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), 
        pt.size = 0.2, ncol = 4)
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, 
            ncol = 3)
FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2, 
            ncol = 3)






