library(dplyr)
library(Linnorm)
library(RColorBrewer)
library(scater)
library(SC3)
library(Seurat)
library(SingleCellExperiment)
library(stringr)

setwd("~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW")

# get the design matrix
design_matrix = read.csv("GSE98969_experimental_design_f.txt",sep = "\t", header = T)
# fix the inconsistent format in Batch_desc column
design_matrix$Batch_desc = gsub(" ", "_", design_matrix$Batch_desc)
# check the number of cells per well
table(design_matrix$Number_of_cells)

# to re-create figure one, select the 6 six-month old AD and WT mouses (there is a mouse 4 in AD)
admouseinfor = dplyr::filter (design_matrix, Mouse_ID=="5XFAD")%>%dplyr::filter(grepl("^AD6m_mouse", Batch_desc))%>%dplyr::filter(!grepl("mouse4", Batch_desc))
wtmouseinfor = dplyr::filter (design_matrix, Mouse_ID=="C57BL/6")%>%dplyr::filter(grepl("^WT6m_mouse", Batch_desc))

ad_wt_table = readRDS("dam.ad_wt.umi.merged.rds")

# check the empty wells
#empty_well = dplyr::filter(design_matrix, Number_of_cells == "0")%>%select(Well_ID)%>%unlist()
#empty_expr = ad_wt_table[,which(colnames(ad_wt_table)%in%empty_well)]

# 192 empty wells
# non-zero expression values in empty wells. remove them from downstream analysis
cell_well = dplyr::filter(design_matrix, Number_of_cells == "1")%>%select(Well_ID)%>%unlist()
ad_wt_table = ad_wt_table[,(which(colnames(ad_wt_table)%in%cell_well))]

# update design matrix, keep only the non-empty wells
ori_ad_wt_matrix = rbind(admouseinfor, wtmouseinfor)
clean_ad_wt_matrix = ori_ad_wt_matrix[(which(ori_ad_wt_matrix$Well_ID%in%cell_well)),]

# use scater to do filtering and plot QC plots
scater_qc = SingleCellExperiment(assays = list(counts = as.matrix(ad_wt_table)), colData = clean_ad_wt_matrix)

# check ERCC and MT genes
total_genes = rownames(ad_wt_table) 

# check MT genes
#mtlist = total_genes[grepl("^mt", total_genes)]
#head(mtlist,5); length(mtlist)
# 37

# check ERCC genes
#ercclist = total_genes[grepl("^ERCC", total_genes)]
#head(ercclist, 5); length(ercclist)
# 92

# annotate the ERCC and MT genes
isSpike(scater_qc, "ERCC") = grepl("^ERCC-", total_genes)
isSpike(scater_qc, "MT") = grepl("^mt-", total_genes)

# plot the percentage of ERCC and mt genes. 
# filter out the ones with > 10% mitochondrial genes

scater_qc = calculateQCMetrics(scater_qc,feature_controls = list(ERCC = isSpike(scater_qc, "ERCC"), 
    MT = isSpike(scater_qc, "MT")))

plotColData(scater_qc,x = "total_features_by_counts",y = "pct_counts_MT",colour = "Amp_batch_ID")
plotColData(scater_qc, x = "total_features_by_counts",y = "pct_counts_ERCC",colour = "Amp_batch_ID")

#mtdf = ad_wt_table[(which(rownames(ad_wt_table)%in%mtlist)),]
#rowSums(mtdf)

# no MT genes found

# remove ERCC genes and mt genes
scater_qc = scater_qc[!(isSpike(scater_qc, "ERCC")),]
scater_qc = scater_qc[!(isSpike(scater_qc, "MT")),]

# keep cells with 500 < UMI < 2500
umi_lower = scater_qc$total_counts_endogenous > 500
umi_upper = scater_qc$total_counts_endogenous < 2500

#table(umi_lower); table(umi_upper)
#3875 cellls < 500, 274 cells > 2500

# keep cells with expressed genes > 300
feature_lower = scater_qc$total_features_by_counts > 300

table(feature_lower)
# 4053 cells < 300

scater_filter_by_cell = scater_qc[,umi_lower&umi_upper&feature_lower]
# 8284 cells left

# remove genes with mean expression smaller than 0.005 UMIs/cell
filter_by_cell_counts = counts(scater_filter_by_cell)
gene_lower = rowSums(filter_by_cell_counts) > (0.005* ncol(filter_by_cell_counts))
#table(gene_lower)
# 10976 genes left

scater_filter_by_cell_by_gene = scater_filter_by_cell[gene_lower,]

# save a copy of file 
saveRDS(scater_filter_by_cell_by_gene, file = "./dam.ad_wt.umi.merged.clean2.rds")

# get the copy of file
ad_wt_clean = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.ad_wt.umi.merged.clean2.rds")

# apply Seurat
ad_wt_seurat = CreateSeuratObject(counts = counts(ad_wt_clean), project = "dam", min.cells = 3, min.features = 200)
#ad_wt_seurat[["percent.ercc"]] = PercentageFeatureSet(ad_wt_seurat, pattern ="^ERCC-" )

VlnPlot(ad_wt_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#cor_plot1 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "percent.ercc")
#cor_plot2 = FeatureScatter(ad_wt_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(cor_plot2))

# QC
ad_wt_seurat_filter = subset(ad_wt_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

ad_wt_seurat_filter = NormalizeData(ad_wt_seurat_filter)

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

# DimHeatmap(ad_wt_seurat_filter, dims = 1, cells = 500, balanced = T)

ad_wt_seurat_filter = FindNeighbors(ad_wt_seurat_filter, dims = 1:40, k.param = 30)
ad_wt_seurat_filter = FindClusters(ad_wt_seurat_filter, resolution = 0.3)
ad_wt_seurat_filter = RunUMAP(ad_wt_seurat_filter, features = VariableFeatures(object = ad_wt_seurat_filter))
ad_wt_seurat_filter = RunTSNE(ad_wt_seurat_filter, features = VariableFeatures(object = ad_wt_seurat_filter))

# UMAP
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

# tSNE
DimPlot(ad_wt_seurat_filter, reduction = "tsne")
# microglial 1
#VlnPlot(ad_wt_seurat_filter, features = c("Cst3", "Hexb") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cst3", "Hexb", "Cst3", "Hexb"), reduction = "tsne")
# microglial 2
# microglial 3
#VlnPlot(ad_wt_seurat_filter, features = c("Apoe", "Cst7", "Lpl", "Cd9", "Trem2") )
FeaturePlot(ad_wt_seurat_filter, features = c("Apoe", "Cst7", "Lpl"), reduction = "tsne")
# perviascular Mf
#VlnPlot(ad_wt_seurat_filter, features = c("Cd74", "Cd163", "Mrc1") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cd74", "Cd163", "Mrc1"), reduction = "tsne")
# Monocytes
#VlnPlot(ad_wt_seurat_filter, features = c("S100a4", "Cd74") )
FeaturePlot(ad_wt_seurat_filter, features = c("S100a4", "Cd74", "S100a4", "Cd74"), reduction = "tsne")
# Mature B-cells and immature B-cells
#VlnPlot(ad_wt_seurat_filter, features = c("Cd74", "Cd79b", "Rag1") )
FeaturePlot(ad_wt_seurat_filter, features = c("Cd74", "Cd79b", "Rag1"), reduction = "tsne")
# T/NK cells
#VlnPlot(ad_wt_seurat_filter, features = c("Nkg7", "Trbc2") )
FeaturePlot(ad_wt_seurat_filter, features = c("Nkg7", "Trbc2", "Nkg7", "Trbc2"), reduction = "tsne")
# Granulocytes 1 &2
#VlnPlot(ad_wt_seurat_filter, features =c("Camp", "S100a9")  )
FeaturePlot(ad_wt_seurat_filter, features = c("Camp", "S100a9", "Camp", "S100a9"), reduction = "tsne")


new_cluster_ids = c("microglia 1", "T/NK", "microglia2/3", "B-cells", "granulocytes", "perviascular Mf")

names(new_cluster_ids) = levels(ad_wt_seurat_filter)
ad_wt_seurat_filter = RenameIdents(ad_wt_seurat_filter, new_cluster_ids)
saveRDS(ad_wt_seurat_filter, file = "~/Documents/2019_summer_eisai/processed_data/dam.seurat.default_cluster.rds")

# print labled map
DimPlot(ad_wt_seurat_filter, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(ad_wt_seurat_filter, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()

DimPlot(ad_wt_seurat_filter, reduction = "tsne", pt.size = 0.5) 

# pinpoint each cluster, redo the coloring
dam_embeds = as.data.frame(Embeddings(ad_wt_seurat_filter[["tsne"]]))
dam_embeds[,"cell"] = rownames(dam_embeds)

# cluster identity of cells
dam_cluster_id = as.data.frame(ad_wt_seurat_filter$seurat_clusters)
dam_cluster_id[,"cell"] = rownames(dam_cluster_id)
colnames(dam_cluster_id) = c("cluster", "cell")

# stage 2 microglia
# subset cluster 2
dam_cluster2 = dplyr::filter(dam_cluster_id, cluster =="2")%>%select(cell)%>%unlist()%>%as.vector()
# select based on coordinates
dam_stage2_coor = dplyr::filter(dam_embeds, -10<tSNE_1&tSNE_1<37&-10<tSNE_2&tSNE_2<25)%>%select(cell)%>%unlist()%>%as.vector()

# microglia 3
dam_stage2 = intersect(dam_cluster2, dam_stage2_coor)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = dam_stage2, value = "microglia 3")

# microglia 2
dam_stage1 = setdiff(dam_cluster2, dam_stage2)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = dam_stage1, value = "microglia 2")

# monocyte
# subset cluster 1
cluster1 = dplyr::filter(dam_cluster_id, cluster =="1")%>%select(cell)%>%unlist()%>%as.vector()
monocyte_coor = dplyr::filter(dam_embeds, 0<tSNE_1&tSNE_1<15&(-40)<tSNE_2&tSNE_2<(-25) )%>%select(cell)%>%unlist()%>%as.vector()
monocytes = intersect(cluster1, monocyte_coor)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = monocytes, value = "monocytes")

# separate immature and mature B cells
# subset cluster 3
cluster3 = dplyr::filter(dam_cluster_id, cluster =="3")%>%select(cell)%>%unlist()%>%as.vector()
immature_coor = dplyr::filter(dam_embeds, tSNE_2<(-5) )%>%select(cell)%>%unlist()%>%as.vector()
inmatureB = intersect(cluster3, immature_coor)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = inmatureB, value = "immature B-cells")

matureB = setdiff(cluster3, inmatureB)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = matureB, value = "mature B-cells")


# granulocytes 1 & 2
# subset cluster 4
cluster4 = dplyr::filter(dam_cluster_id, cluster =="4")%>%select(cell)%>%unlist()%>%as.vector()
granu2_coor = dplyr::filter(dam_embeds, 30<tSNE_1&tSNE_2<(-26) )%>%select(cell)%>%unlist()%>%as.vector()
granu2 = intersect(cluster4, granu2_coor)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = granu2, value = "granulocytes 2")

granu1 = setdiff(cluster4, granu2)
ad_wt_seurat_filter = SetIdent(object = ad_wt_seurat_filter, cells = granu1, value = "granulocytes 1")

# tester of positions
#DimPlot(object = ad_wt_seurat_filter, cells.highlight = granu2, reduction = "tsne")


DimPlot(ad_wt_seurat_filter, reduction = "tsne", label = F, pt.size = 0.5) 

saveRDS(ad_wt_seurat_filter, file = "~/Documents/2019_summer_eisai/processed_data/dam.seurat.reassign_cluster.rds")


display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=FALSE)
display.brewer.pal(11, "Spectral")
brewer.pal(11, "Spectral")

replot = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.seurat.reassign_cluster.rds")
DimPlot(replot, reduction = "tsne", label = F, pt.size = 0.5, order= rev(c("microglia 1",
        "microglia 2", "microglia 3", "perviascular Mf", "monocytes", "mature B-cells",
        "immature B-cells", "T/NK", "granulocytes 1", "granulocytes 2")), cols = c("#FEE08B",
                                                                                   "#F46D43",  "#9E0142", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2",  "#D53E4F", "#FDAE61", "#E6F598")  ) 

# print out the top cluster markers 

ad_wt.markers = FindAllMarkers(ad_wt_seurat_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ad_wt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


#saveRDS(ad_wt_seurat_filter, file = "~/Documents/2019_summer_eisai/progress_report/dam_seurat_v1.rds")

#ad_wt_seurat_filter = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.seurat.default_cluster.rds")

DimPlot(ad_wt_seurat_filter, reduction = "tsne")

# get the genes differentially expressed in the three microglia clusters
# up in 1 to 2
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = "microglia 2", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)
# up in 1 to 3
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = "microglia 3", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_3.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 1 to stage 2&3
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = c("microglia 2", "microglia 3"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_23.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 1
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = "microglia 1", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 3
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = "microglia 3", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_3.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 1&3
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = c("microglia 1", "microglia 3"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_13.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 1
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = "microglia 1", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 2
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = "microglia 2", min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 1&2
write.table(FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = c("microglia 1", "microglia 2"), min.pct = 0.25), 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_12.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)





