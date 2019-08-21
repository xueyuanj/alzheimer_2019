library(dplyr)
library(GEOquery)
library(Seurat)

setwd("/public_genomic_data_downloads/czhang/sharon/fat_trem2/GSE128518_RAW")

# get the accession number
gse128518 = getGEO("GSE128518")
h_design_infor = Biobase::pData(phenoData(gse128518[[1]]))
m_design_infor = Biobase::pData(phenoData(gse128518[[2]]))

intresting= c("title", "cell type:ch1", "mouse age (weeks):ch1", "selection marker:ch1",
              "strain:ch1",  "time point (weeks):ch1", "treatment (diet):ch1")

sapply(intresting, function(x)
  print(unique(m_design_infor[,x])))

# figure1, wild type NC and HFD mice at 6, 12, 18 weeks
m_nc_hfd = dplyr::filter(m_design_infor, `strain:ch1` == "C57BL/6")
m_nc_hfd_supp = select(m_nc_hfd, supplementary_file_1)%>%unlist()%>%as.character() 
m_nc_hfd_f = sapply(m_nc_hfd_supp, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 


wt_m = lapply(m_nc_hfd_f, function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

wt_m_df = data.frame(wt_m)

# save a copy of UMI file
saveRDS(wt_m_df, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.wt_nc_hfd.umi.rds")


setwd("/public_genomic_data_downloads/czhang/sharon/fat_trem2")

# QC part for the wildtype mouse
wt_nc_hfd_df = readRDS("/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.wt_nc_hfd.umi.rds")

totalgenes = rownames(wt_nc_hfd_df)

# identify the genes to be removed from analysis
mtgenes = totalgenes[grepl("^mt", totalgenes)]
erccgene = totalgenes[grepl("^ERCC", totalgenes)]
rpgene = totalgenes[grepl("^Rp", totalgenes)]
immunoglobin = totalgenes[grepl("^Ig", totalgenes)]

endogenes = totalgenes[!totalgenes %in% mtgenes & !totalgenes %in% erccgene & !totalgenes %in% rpgene & !totalgenes %in% immunoglobin]

# seurat part
wtm_seurat = CreateSeuratObject(counts = wt_nc_hfd_df, assay = "RNA", project = "nc_hfd", min.cells = 3, min.features = 200)

# check MT gene and ERCC
wtm_seurat[["percent.mt"]] = PercentageFeatureSet(wtm_seurat, pattern = "^mt")
wtm_seurat[["percent.ERCC"]] = PercentageFeatureSet(wtm_seurat, pattern = "^ERCC")

VlnPlot(wtm_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(wtm_seurat, features = c("percent.mt", "percent.ERCC"), ncol = 2)


cor_plot1 = FeatureScatter(wtm_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
cor_plot2 = FeatureScatter(wtm_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cor_plot3 = FeatureScatter(wtm_seurat, feature1 = "nCount_RNA", feature2 = "percent.ERCC")

CombinePlots(plots = list(cor_plot1, cor_plot2, cor_plot3))

# filter cells/genes 
wtm_seurat_filter = subset(wtm_seurat, subset = nCount_RNA > 400  & percent.mt < 40 & percent.ERCC < 20)

wtm_seurat_filter = wtm_seurat_filter[endogenes,]
# 25951 features across 20157 samples

# normalize
wtm_seurat_filter = NormalizeData(wtm_seurat_filter)

# feature selection
wtm_seurat_filter = FindVariableFeatures(wtm_seurat_filter, selection.method = "disp", nfeatures = 425)

wtm_top20 = head(VariableFeatures(wtm_seurat_filter), 20)
wt_varplot1 = VariableFeaturePlot(wtm_seurat_filter)
wt_varplot2 = LabelPoints(plot = wt_varplot1, points = wtm_top20, repel = T)
CombinePlots(plots = list(wt_varplot2))

# dimention reduction
wtm_seurat_filter = ScaleData(wtm_seurat_filter)

# PCA
wtm_seurat_filter = RunPCA(wtm_seurat_filter, features = VariableFeatures(object = wtm_seurat_filter))

# tSNE
wtm_seurat_filter = FindNeighbors(wtm_seurat_filter, dims = 1:15, k.param = 30)
wtm_seurat_filter = FindClusters(wtm_seurat_filter, resolution = 0.2)
wtm_seurat_filter = RunTSNE(wtm_seurat_filter, features = VariableFeatures(object = wtm_seurat_filter))
#wtm_seurat_filter = RunUMAP(wtm_seurat_filter, features = VariableFeatures(object = wtm_seurat_filter))

DimPlot(wtm_seurat_filter, reduction = "tsne")
#DimPlot(wtm_seurat_filter, reduction = "umap")
#ori_cluster_id = as.data.frame( wtm_seurat_filter$seurat_clusters)
#ori_cluster_id[,"cell"] = rownames(ori_cluster_id)
#colnames(ori_cluster_id) = c("cluster", "cell")

# where is cluster 4
#cluster4 = dplyr::filter(ori_cluster_id, cluster == as.character(4))%>%select(cell)%>%unlist()%>%as.vector()
#DimPlot(object = wtm_seurat_filter, cells.highlight = cluster4, reduction = "tsne")


cellmarkers= c("Lyz1", "Fn1", "Mrc1", "C1qc", "Apoe", "Cd79b", "Ly6d", "Fcer1a", "S100a9",
               "Irf8", "Cd209a", "Siglech", "Foxp3", "Cd3d", "Gzma")

# monocyte
FeaturePlot(wtm_seurat_filter, features = c("Lyz1", "Fn1", "Mrc1"), reduction = "tsne")

# macrophage
FeaturePlot(wtm_seurat_filter, features = c("Mrc1", "C1qc", "Apoe"), reduction = "tsne")

# B cells
FeaturePlot(wtm_seurat_filter, features = c("Cd79b", "Ly6d", "Cd79b", "Ly6d"), reduction = "tsne")

# mast cell
# neu
# cDC1
# cDC2
FeaturePlot(wtm_seurat_filter, features = c("Fcer1a", "S100a9",
                                            "Irf8", "Cd209a"), reduction = "tsne")
# pDC
# T cells
FeaturePlot(wtm_seurat_filter, features = c("Siglech", "Foxp3", "Cd3d", "Gzma"), reduction = "tsne")

# hard to annotate the clusters
# save of copy of the original tSNE plot

saveRDS(wtm_seurat_filter, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.wt_nc_hfd.seurat.rds")


# highlight the cells of different origin
# obrain the cell metadata
allcellmeta = read.csv("GSE128518_metadata.txt", header = T, sep = "\t")

# HFD 
hfdmouse = dplyr::filter(m_design_infor, `treatment (diet):ch1` == "High-fat diet" & `strain:ch1` == "C57BL/6")

# HFD different age
week6hfd = dplyr::filter(hfdmouse, `time point (weeks):ch1` == "6")%>%select(title)%>%unlist()%>%as.vector()
week12hfd = dplyr::filter(hfdmouse, `time point (weeks):ch1` == "12")%>%select(title)%>%unlist()%>%as.vector()
week18hfd = dplyr::filter(hfdmouse, `time point (weeks):ch1` == "18")%>%select(title)%>%unlist()%>%as.vector()

# select the cells correspond to each treatment
hfd6cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week6hfd)%>%select(Well_ID)%>%unlist()%>%as.vector()
hfd12cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week12hfd)%>%select(Well_ID)%>%unlist()%>%as.vector()
hfd18cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week18hfd)%>%select(Well_ID)%>%unlist()%>%as.vector()

# plot those on the tSNE
DimPlot(object = wtm_seurat_filter, cells.highlight = hfd6cells, reduction = "tsne")
DimPlot(object = wtm_seurat_filter, cells.highlight = hfd12cells, reduction = "tsne")
DimPlot(object = wtm_seurat_filter, cells.highlight = hfd18cells, reduction = "tsne")

# NC
ncmouse = dplyr::filter(m_design_infor, `treatment (diet):ch1` == "Normal chow" & `strain:ch1` == "C57BL/6")

# NC different age

# NC different age
week6nc = dplyr::filter(ncmouse, `time point (weeks):ch1` == "6")%>%select(title)%>%unlist()%>%as.vector()
week12nc = dplyr::filter(ncmouse, `time point (weeks):ch1` == "12")%>%select(title)%>%unlist()%>%as.vector()
week18nc = dplyr::filter(ncmouse, `time point (weeks):ch1` == "18")%>%select(title)%>%unlist()%>%as.vector()

# select the cells correspond to each treatment
nc6cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week6nc)%>%select(Well_ID)%>%unlist()%>%as.vector()
nc12cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week12nc)%>%select(Well_ID)%>%unlist()%>%as.vector()
nc18cells = dplyr::filter(allcellmeta, Amp_batch_ID%in%week18nc)%>%select(Well_ID)%>%unlist()%>%as.vector()

# plot those on the tSNE
DimPlot(object = wtm_seurat_filter, cells.highlight = nc6cells, reduction = "tsne")
DimPlot(object = wtm_seurat_filter, cells.highlight = nc12cells, reduction = "tsne")
DimPlot(object = wtm_seurat_filter, cells.highlight = nc18cells, reduction = "tsne")


# aliquot the 2 major cell types for downstream analysis
wtm_2celltype = subset(wtm_seurat_filter, idents = c("0",  "3"))

# redo the KNN and tSNE plot
wtm_2celltype =  FindVariableFeatures(wtm_2celltype, selection.method = "disp", nfeatures = 425)

wtm_2celltype_top20 = head(VariableFeatures(wtm_2celltype), 20)
wt_2celltype_varplot1 = VariableFeaturePlot(wtm_2celltype)
wt_2celltype_varplot2 = LabelPoints(plot = wt_2celltype_varplot1, points =wtm_2celltype_top20, repel = T)
CombinePlots(plots = list(wt_2celltype_varplot2))

# dimention reduction
wtm_2celltype = ScaleData(wtm_2celltype)

# PCA
wtm_2celltype = RunPCA(wtm_2celltype, features = VariableFeatures(object = wtm_2celltype))

# tSNE
wtm_2celltype = FindNeighbors(wtm_2celltype, dims = 1:15, k.param = 30)
wtm_2celltype = FindClusters(wtm_2celltype, resolution = 0.125)
wtm_2celltype = RunTSNE(wtm_2celltype, features = VariableFeatures(object = wtm_2celltype))
#wtm_seurat_filter = RunUMAP(wtm_seurat_filter, features = VariableFeatures(object = wtm_seurat_filter))

DimPlot(wtm_2celltype, reduction = "tsne")

# plot the markers
# monocytes
FeaturePlot(wtm_2celltype, features = c("Retnla", "Fn1", "Plac8", "Clec4e"), reduction = "tsne")

# macrophage
FeaturePlot(wtm_2celltype, features = c("Lyve1", "Cd209f", "Cd163", "Nceh1"), reduction = "tsne")
FeaturePlot(wtm_2celltype, features = c("Cd9", "Spp1", "Trem2", "Lpl"), reduction = "tsne")

# highlight the cells of different time points
DimPlot(object = wtm_2celltype, cells.highlight = hfd6cells, reduction = "tsne")
DimPlot(object = wtm_2celltype, cells.highlight = hfd12cells, reduction = "tsne")
DimPlot(object = wtm_2celltype, cells.highlight = hfd18cells, reduction = "tsne")

# annotate cells
mon_mac_coor = as.data.frame(Embeddings(wtm_2celltype[["tsne"]]))
mon_mac_coor[,"cell"] = rownames(mon_mac_coor)

mon_mac_cluster_id = as.data.frame(wtm_2celltype$seurat_clusters)
mon_mac_cluster_id[,"cell"] = rownames(mon_mac_cluster_id)
colnames(mon_mac_cluster_id) = c("cluster", "cell")

# monocytes
mono_cells = dplyr::filter(mon_mac_cluster_id, cluster =="1")%>%select(cell)%>%unlist()%>%as.vector()

# macrophage 1
mac_1_cells = dplyr::filter(mon_mac_cluster_id, cluster =="2")%>%select(cell)%>%unlist()%>%as.vector()

# macrophage 2
mac23cluster = dplyr::filter(mon_mac_cluster_id, cluster =="0")%>%select(cell)%>%unlist()%>%as.vector()
mac2_cell_range = dplyr::filter(mon_mac_coor, tSNE_1<12)%>%select(cell)%>%unlist()%>%as.vector()
mac2_cells = intersect(mac23cluster, mac2_cell_range)

# macrophage 3
mac3_cells = setdiff(mac23cluster, mac2_cells)

wtm_2celltype = SetIdent(object = wtm_2celltype, cells = mono_cells, value = "monocytes")
wtm_2celltype = SetIdent(object = wtm_2celltype, cells = mac_1_cells, value = "macrophage 1")
wtm_2celltype = SetIdent(object = wtm_2celltype, cells = mac2_cells, value = "macrophage 2")
wtm_2celltype = SetIdent(object = wtm_2celltype, cells = mac3_cells, value = "macrophage 3")

DimPlot(wtm_2celltype, reduction = "tsne")

# save a copy of the data
saveRDS(wtm_2celltype, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.mon.mac.rds")

# get DEG between Mac1 and Mac 3
mac1_3 = FindMarkers(wtm_2celltype, ident.1 = "macrophage 3", ident.2 ="macrophage 1", min.pct = 0.25)
mac1_3[,"gene"] = rownames(mac1_3)
mac1_3_sig = dplyr::filter(mac1_3, p_val_adj < 0.001)

# keep a copy of the file
write.table(mac1_3_sig, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/macrophage.1.3.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

upinmac3 =  dplyr::filter(mac1_3_sig,avg_logFC>0.5)
upinmac1 = dplyr::filter(mac1_3_sig,avg_logFC<0)


write.table(head(upinmac3[order(upinmac3$avg_logFC, decreasing = T),], 20),
            file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/macrophage.up3to1.top20.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
head(upinmac1, 20)

# get DEG between Mac2 and Mac 3
mac2_3 = FindMarkers(wtm_2celltype, ident.1 = "macrophage 3", ident.2 ="macrophage 2", min.pct = 0.25)
mac2_3[,"gene"] = rownames(mac2_3)
mac2_3_sig = dplyr::filter(mac2_3, p_val_adj < 0.001)

# keep a copy of the file
write.table(mac2_3_sig, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/macrophage.2.3.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

upinmac3 =  dplyr::filter(mac2_3_sig,avg_logFC>0)
upinmac1 = dplyr::filter(mac2_3_sig,avg_logFC<0)

head(upinmac1, 20)
head(upinmac3, 20)



# get DEG between Mac1 and Mac 2
mac1_2 = FindMarkers(wtm_2celltype, ident.1 = "macrophage 2", ident.2 ="macrophage 1", min.pct = 0.25)
mac1_2[,"gene"] = rownames(mac1_2)
mac1_2_sig = dplyr::filter(mac1_2, p_val_adj < 0.001)

# keep a copy of the file
write.table(mac1_2_sig, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/macrophage.1.2.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

upinmac2 =  dplyr::filter(mac1_2_sig,avg_logFC>0)
upinmac1 = dplyr::filter(mac1_2_sig,avg_logFC<0)

head(upinmac1, 20)
head(upinmac2, 20)



# download the data
# sapply(oat_urls, function(x) download.file(x, unlist(strsplit(x, "/"))[9]))

# trem2 KO part

setwd("/public_genomic_data_downloads/czhang/sharon/fat_trem2/GSE128518_RAW")

#hfd_trem_part = rbind(trem2_ko_hfd, litterm_hfd)


# ignore the paper, QC with both Trem2 KO and their WT littermates (both HFD and NC)

trem_part = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 knock-out C57BL/6"|`strain:ch1` == "Trem2 wild-type C57BL/6 Littermates")

trem_supp = select(trem_part, supplementary_file_1)%>%unlist()%>%as.character() 
trem_f = sapply(trem_supp, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 


trem_df = lapply(trem_f, function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

trem_data = data.frame(trem_df)

#> dim(trem_data)
#[1] 52634 14976
# before integrate NC, 10368 cells. Very little NC cells anyway


# QC of the cells from high-fat diet
trem_seurat = CreateSeuratObject(counts = trem_data, assay = "RNA", project = "Trem2 part", min.cells = 3,
                                     min.features = 200)

# save a copy
saveRDS(trem_seurat, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.trem2.seurat.raw.rds")

trem_seurat[["percent.mt"]] = PercentageFeatureSet(trem_seurat, pattern = "^mt")
trem_seurat[["percent.ERCC"]] = PercentageFeatureSet(trem_seurat, pattern = "^ERCC")

VlnPlot(trem_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(trem_seurat, features = c("percent.mt", "percent.ERCC"), ncol = 2)

# no ERCC in this dataset

cor_plot1 = FeatureScatter(trem_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
cor_plot2 = FeatureScatter(trem_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#cor_plot3 = FeatureScatter(trem_seurat, feature1 = "nCount_RNA", feature2 = "percent.ERCC")

CombinePlots(plots = list(cor_plot1, cor_plot2))

# filter out low-quality ones
trem_seurat_filter = subset(trem_seurat, subset =  nCount_RNA > 400  & percent.mt < 40)

trem_seurat_filter = trem_seurat_filter[endogenes,]
# 24218 genes, 10792 cells

# normalize
trem_seurat_filter = NormalizeData(trem_seurat_filter)
trem_seurat_filter = FindVariableFeatures(trem_seurat_filter, selection.method = "disp", nfeatures = 425)

# run tSNE, isolate monocytes and macrophages
trem_seurat_filter = ScaleData(trem_seurat_filter)
trem_seurat_filter = RunPCA(trem_seurat_filter, features = VariableFeatures(object = trem_seurat_filter))

trem_seurat_filter = FindNeighbors(trem_seurat_filter, dims = 1:15, k.param = 30)
trem_seurat_filter = FindClusters(trem_seurat_filter, resolution = 0.2)
trem_seurat_filter = RunTSNE(trem_seurat_filter, features = VariableFeatures(object = trem_seurat_filter))
#wtm_seurat_filter = RunUMAP(wtm_seurat_filter, features = VariableFeatures(object = wtm_seurat_filter))

DimPlot(trem_seurat_filter, reduction = "tsne")

# plot the marker genes
# monocyte
FeaturePlot(trem_seurat_filter, features = c("Lyz1", "Fn1", "Mrc1"), reduction = "tsne")
# cluster 4

# macrophage
FeaturePlot(trem_seurat_filter, features = c("Mrc1", "C1qc", "Apoe"), reduction = "tsne")
# cluster 0

# save a copy
saveRDS(trem_seurat_filter, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.trem2.seurat.clustered.rds")

# isolate cluster 0 and 4
trem2_mon_mac = subset(trem_seurat_filter, idents = c("0", "4"))

# quick check
#DimPlot(trem2_mon_mac, reduction = "tsne")

wtm_2celltype = FindVariableFeatures(wtm_2celltype, selection.method = "disp", nfeatures = 200)
trem2_mon_mac = FindVariableFeatures(trem2_mon_mac, selection.method = "disp", nfeatures = 200)

trem_anchor = FindTransferAnchors(reference = wtm_2celltype, query = trem2_mon_mac, dims = 1:30)
predictions = TransferData(anchorset = trem_anchor, refdata = Idents(wtm_2celltype))
trem2_mon_mac= AddMetaData(trem2_mon_mac, metadata = predictions)

# plot the expression of gene markers
# Mac 1
VlnPlot(trem2_mon_mac, c("C4b", "F13a1"), group.by = "predicted.id")
VlnPlot(trem2_mon_mac, c("Cbr2", "Lyve1"), group.by = "predicted.id")

# Mac 2 and LAM
VlnPlot(trem2_mon_mac, c("Trem2", "Nceh1"), group.by = "predicted.id")
VlnPlot(trem2_mon_mac, c("Cd36", "Cd9"), group.by = "predicted.id")
VlnPlot(trem2_mon_mac, c("Ctsl", "Lpl", "Fabp4"), group.by = "predicted.id")

# save a copy
saveRDS(trem2_mon_mac, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.trem2.projection.rds")

# aliquote cells of different conditions
# select each dataframe
trem2_ko_hfd = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 knock-out C57BL/6"&`treatment (diet):ch1` == "High-fat diet")%>%select(title)%>%unlist()%>%as.vector()
trem2_ko_nc = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 knock-out C57BL/6"&`treatment (diet):ch1` == "Normal chow")%>%select(title)%>%unlist()%>%as.vector()
litterm_hfd = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 wild-type C57BL/6 Littermates"&`treatment (diet):ch1` == "High-fat diet")%>%select(title)%>%unlist()%>%as.vector()
litterm_nc = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 wild-type C57BL/6 Littermates"&`treatment (diet):ch1` == "Normal chow")%>%select(title)%>%unlist()%>%as.vector()

# select cells of each category
cells_trem2_ko_hfd = dplyr::filter(allcellmeta, Amp_batch_ID%in%trem2_ko_hfd)%>%select(Well_ID)%>%unlist()%>%as.vector()
cells_litterm_hfd = dplyr::filter(allcellmeta, Amp_batch_ID%in%litterm_hfd)%>%select(Well_ID)%>%unlist()%>%as.vector()
cells_trem2_ko_nc = dplyr::filter(allcellmeta, Amp_batch_ID%in%trem2_ko_nc)%>%select(Well_ID)%>%unlist()%>%as.vector()
cells_litterm_nc = dplyr::filter(allcellmeta, Amp_batch_ID%in%litterm_nc)%>%select(Well_ID)%>%unlist()%>%as.vector()
  
seurat_trem2_ko_hfd = subset(trem2_mon_mac, cells= cells_trem2_ko_hfd)
seurat_litterm_hfd = subset(trem2_mon_mac, cells= cells_litterm_hfd)
seurat_trem2_ko_nc = subset(trem2_mon_mac, cells= cells_trem2_ko_nc)
seurat_litterm_nc = subset(trem2_mon_mac, cells= cells_litterm_nc)

summarytable = rbind(table(seurat_litterm_nc$predicted.id),
                     table(seurat_trem2_ko_nc$predicted.id),
                     table(seurat_litterm_hfd$predicted.id),
                     table(seurat_trem2_ko_hfd$predicted.id)
                     )
rownames(summarytable) = c("wt_nc", "ko_nc", "wt_hfd", "ko_hfd")
summarytable = t(summarytable)



barplot(summarytable, col = c("darkseagreen1", "slategray1", "mediumpurple", "indianred1")) 
legend(0.2, 2500 ,legend = rev(rownames(summarytable)),
       col = rev(c("darkseagreen1", "slategray1", "mediumpurple", "indianred1")),
       lty=1:2, cex=1)

# get the percentage one
percent_table =t(t(summarytable)/colSums(summarytable)) 

print(percent_table)
barplot(percent_table, col = c("darkseagreen1", "slategray1", "mediumpurple", "indianred1")) 


############################################################################################
# human part


# the single obese individual
# the whole adipose tissue
human_eat = dplyr::filter(h_design_infor, `organ:ch1` == "EAT")

human_matrix = select(human_eat, supplementary_file_3)%>%unlist()%>%as.character() 

human_m_f = sapply(human_matrix, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 

obese_indi = lapply(human_m_f, function(x){
  readMM(list.files(pattern = x))})

# cell barcode
human_barcodes = select(human_eat, supplementary_file_1)%>%unlist()%>%as.character() 
human_bar_f = sapply(human_barcodes, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 

obese_barcodes = lapply(human_bar_f, function(x){
  read.table(list.files(pattern = x), header = F)})


# features
human_features = select(human_eat, supplementary_file_2)%>%unlist()%>%as.character() 
human_feature_f = sapply(human_features, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 
obese_feature = lapply(human_feature_f, function(x){
  read.table(list.files(pattern = x), header = F)})

# add rownames and colnames to the dataframe
colnames(obese_indi[[1]]) = obese_barcodes[[1]][,1]
colnames(obese_indi[[2]]) = obese_barcodes[[2]][,1]

rownames(obese_indi[[1]]) = obese_feature[[1]][,2]
rownames(obese_indi[[2]]) = obese_feature[[2]][,2]

# to matrix are too large to handle as a whole
# remove empy columns
human_part1 = obese_indi[[1]]
human_part1 = human_part1[,colSums(human_part1 != 0) >200]

human_part2 = obese_indi[[2]]
human_part2 = human_part2[,colSums(human_part2 != 0) >200]

human_obese_df = data.frame(human_part1, human_part2)

h_totalgene = rownames(human_obese_df)
h_mtgene = h_totalgene[grepl("^MT", h_totalgene)]
h_erccgene =  h_totalgene[grepl("^ERCC", h_totalgene)]
h_rpgene =  h_totalgene[grepl("^RP", h_totalgene)]
h_immunogene =  h_totalgene[grepl("^IG", h_totalgene)]

h_endogenes = h_totalgene[!h_totalgene %in% h_mtgene & !h_totalgene %in% h_erccgene & !h_totalgene %in% h_rpgene & !h_totalgene %in% h_immunogene]

# merge the two and create a Seurat object 
human_obese_seurat = CreateSeuratObject(counts = human_obese_df, assay = "RNA", project = "obese_individual_1",  min.cells = 3)

# save a copy
saveRDS(human_obese_seurat, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.obese_individual.seurat.raw.rds")

# QC of the data
# MT gene and ERCC
human_obese_seurat[["percent.mt"]] = PercentageFeatureSet(human_obese_seurat, pattern = "^MT")
human_obese_seurat[["percent.ERCC"]] = PercentageFeatureSet(human_obese_seurat, pattern = "^ERCC")

VlnPlot(human_obese_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(human_obese_seurat, features = c("percent.mt", "percent.ERCC"), ncol = 2)

cor_plot1 = FeatureScatter(human_obese_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
cor_plot2 = FeatureScatter(human_obese_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cor_plot3 = FeatureScatter(human_obese_seurat, feature1 = "nCount_RNA", feature2 = "percent.ERCC")

CombinePlots(plots = list(cor_plot1, cor_plot2, cor_plot3))

# filter the genes and cells
human_obese_seurat_filter = subset(human_obese_seurat, subset = nCount_RNA > 400  & percent.mt < 10 & percent.ERCC < 20)

human_obese_seurat_filter = human_obese_seurat_filter[h_endogenes,]
# 17032 genes, 15387 cells

# normalize
human_obese_seurat_filter = NormalizeData(human_obese_seurat_filter)

# feature selection
human_obese_seurat_filter = FindVariableFeatures(human_obese_seurat_filter, selection.method = "disp", nfeatures = 500)

h_obese_top20 = head(VariableFeatures(human_obese_seurat_filter), 20)
h_obese_plot1 = VariableFeaturePlot(human_obese_seurat_filter)
h_obese_plot2 = LabelPoints(plot = h_obese_plot1, points = h_obese_top20, repel = T)
CombinePlots(plots = list(h_obese_plot2))

# dimention reduction
human_obese_seurat_filter = ScaleData(human_obese_seurat_filter)

# PCA
human_obese_seurat_filter = RunPCA(human_obese_seurat_filter, features = VariableFeatures(object = human_obese_seurat_filter))

# tSNE
human_obese_seurat_filter = FindNeighbors(human_obese_seurat_filter, dims = 1:15, k.param = 30)
human_obese_seurat_filter = FindClusters(human_obese_seurat_filter, resolution = 0.2)
human_obese_seurat_filter = RunTSNE(human_obese_seurat_filter, features = VariableFeatures(object = human_obese_seurat_filter))

DimPlot(human_obese_seurat_filter, reduction = "tsne")

# identify the clusters with cell type markers
# macrophages
FeaturePlot(human_obese_seurat_filter, features = c("MRC1", "C1QC", "APOE", "C1QB"), reduction = "tsne")

# LAM
FeaturePlot(human_obese_seurat_filter, features = c("CHIT1", "TREM2", "LIPA", "SPP1"), reduction = "tsne")

# save a copy of this
saveRDS(human_obese_seurat_filter, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.obese_individual.seurat.clustered.rds")

# isolate the macrophage
# annotate the important ones
#h_obese_cell_id = as.data.frame(human_obese_seurat_filter$seurat_clusters)
#h_obese_cell_id[,"cell"] = rownames(h_obese_cell_id)
#colnames(h_obese_cell_id) = c("cluster", "cell")

#h_obese_coordi = as.data.frame(Embeddings(human_obese_seurat_filter[["tsne"]]))
#h_obese_coordi[,"cell"] = rownames(h_obese_coordi)

# then isolate LAM
#DimPlot(human_obese_seurat_filter, reduction = "tsne") + geom_abline(slope = -0.5, intercept = -14)

#human_obese_mac = dplyr::filter(h_obese_cell_id, cluster =="3")%>%select(cell)%>%unlist()%>%as.vector()
#lamrange = dplyr::filter(h_obese_coordi, tSNE_2< (-0.5*tSNE_1 -14) )%>%select(cell)%>%unlist()%>%as.vector()

# isolate the macrophage, subcluster them
h_obese_mac= subset(human_obese_seurat_filter, idents = "3")

h_obese_mac = FindVariableFeatures(h_obese_mac, selection.method = "disp", nfeatures = 500)
h_obese_mac = ScaleData(h_obese_mac)
h_obese_mac = RunPCA(h_obese_mac, features = VariableFeatures(object = h_obese_mac))
h_obese_mac = FindNeighbors(h_obese_mac, dims = 1:15, k.param = 30)
h_obese_mac = FindClusters(h_obese_mac, resolution = 0.2)
h_obese_mac = RunTSNE(h_obese_mac, features = VariableFeatures(object = h_obese_mac))
DimPlot(h_obese_mac, reduction = "tsne")

# macrophage markers
FeaturePlot(h_obese_mac, features = c("MRC1", "C1QC", "APOE", "C1QB"), reduction = "tsne")

# ID the LAM by plotting the markers
FeaturePlot(h_obese_mac, features = c("CHIT1", "TREM2", "LIPA", "SPP1"), reduction = "tsne")

# save a copy
saveRDS(h_obese_mac, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.obese_individual.seurat.macrophage.clustered.rds")

# DEG between LAM and all the other macrophages
lam_othermac = FindMarkers(h_obese_mac, ident.1 = "3" , ident.2 = c("0", "1", "2", "4"), min.pct = 0.25)

# save a copy
write.table(lam_othermac, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/human.lam.all_other_macrophage.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

lam_othermac = read.table("/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/human.lam.all_other_macrophage.deg.txt")
lam_othermac[, "gene"] = rownames(lam_othermac)
lam_sig = dplyr::filter(lam_othermac, p_val_adj < 0.001)

upinlam = dplyr::filter(lam_sig,avg_logFC>0.5)
downinlam =  dplyr::filter(lam_sig,avg_logFC<0)

upinlam = upinlam[order(upinlam$avg_logFC, decreasing = T),]
downinlam = downinlam[order(downinlam$avg_logFC, decreasing = T),]

head(upinlam, 20)

write.table(head(upinlam, 20), file="/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/human.up_in_lam.all_other_macrophage.top20.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
head(downinlam, 20)

# the proportion of LAM out of all macrophages
# isolate the LAM cells 
# (deal with this later)



# QC all the cells from the 6 individual
# shown here are all the immune cells
human_oat = dplyr::filter(h_design_infor, `organ:ch1` == "OAT")
human_supp = select(human_oat, supplementary_file_1)%>%unlist()%>%as.character() 

human_f = sapply(human_supp, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 

human_indiv = lapply(human_f, function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

human_df = data.frame(human_indiv)

#dim(human_df)
#[1] 57874  8832


# cell-individual map
bmi23 = dplyr::filter(human_oat, grepl("23", characteristics_ch1.6))
bmi37 = dplyr::filter(human_oat, grepl("37", characteristics_ch1.6))
bmi38 = dplyr::filter(human_oat, grepl("38", characteristics_ch1.6))
bmi41 = dplyr::filter(human_oat, grepl("41", characteristics_ch1.6))
bmi43 = dplyr::filter(human_oat, grepl("43", characteristics_ch1.6))
bmi44 = dplyr::filter(human_oat, grepl("44", characteristics_ch1.6))

# wrapper to extract the well id
selectwell = function(x){
  human_supp = select(x, supplementary_file_1)%>%unlist()%>%as.character() 
  
  human_f = sapply(human_supp, function(x) {
    gzfile = unlist(strsplit(x, "/"))[9] 
    unzipfile = sub(".gz", "", gzfile)
  },  USE.NAMES = F ) 
  
  human_indiv = lapply(human_f, function(x){
    read.csv(list.files(pattern = x), sep="\t", header=TRUE)})
  
  human_df = data.frame(human_indiv)
  
  return(colnames(human_df))
  
}

# cells belong to each individual
bmi23_well = selectwell(bmi23)
# 1536
bmi37_well = selectwell(bmi37)
# 1536
bmi38_well = selectwell(bmi38)
# 1152
bmi41_well = selectwell(bmi41)
# 1920
bmi43_well = selectwell(bmi43)
# 1536
bmi44_well = selectwell(bmi44)
# 1152

human_other6_seurat = CreateSeuratObject(counts = human_df, assay = "RNA", project = "All other human", min.cells = 3,
                                         min.features = 200)

# save a copy for demo
saveRDS(human_other6_seurat, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.other6.seurat.raw.rds")

human_other6_seurat[["percent.mt"]] = PercentageFeatureSet(human_other6_seurat, pattern = "^MT")
human_other6_seurat[["percent.ERCC"]] = PercentageFeatureSet(human_other6_seurat, pattern = "^ERCC")

VlnPlot(human_other6_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(human_other6_seurat, features = c("percent.mt", "percent.ERCC"), ncol = 2)

cor_plot1 = FeatureScatter(human_other6_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
cor_plot2 = FeatureScatter(human_other6_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cor_plot3 = FeatureScatter(human_other6_seurat, feature1 = "nCount_RNA", feature2 = "percent.ERCC")

CombinePlots(plots = list(cor_plot1, cor_plot2, cor_plot3))

# filter the genes and cells
human_other6_seurat_filter = subset(human_other6_seurat, subset = nCount_RNA > 400  & percent.mt < 30 & percent.ERCC < 25)

human_other6_seurat_filter = human_other6_seurat_filter[h_endogenes,]
# 14338 genes, 4760 cells

# normalize
human_other6_seurat_filter = NormalizeData(human_other6_seurat_filter)
human_other6_seurat_filter = FindVariableFeatures(human_other6_seurat_filter, selection.method = "disp", nfeatures = 500)

# run tSNE
human_other6_seurat_filter = ScaleData(human_other6_seurat_filter)
human_other6_seurat_filter = RunPCA(human_other6_seurat_filter, features = VariableFeatures(object = human_other6_seurat_filter))

human_other6_seurat_filter = FindNeighbors(human_other6_seurat_filter, dims = 1:15, k.param = 30)
human_other6_seurat_filter = FindClusters(human_other6_seurat_filter, resolution = 0.2)
human_other6_seurat_filter = RunTSNE(human_other6_seurat_filter, features = VariableFeatures(object = human_other6_seurat_filter))

DimPlot(human_other6_seurat_filter, reduction = "tsne")

# save a copy
saveRDS(human_other6_seurat_filter, file =  "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.other6.seurat.filtered.rds")

# macrophage and monocytes
# macrophages
FeaturePlot(human_other6_seurat_filter, features = c("MRC1", "C1QC", "APOE", "C1QB"), reduction = "tsne")

# LAM
FeaturePlot(human_other6_seurat_filter, features = c("CHIT1", "TREM2", "LIPA", "SPP1"), reduction = "tsne")

# isolate cluster 1 - macrophage
h_other6_mac= subset(human_other6_seurat_filter, idents = "1")

# subclustering
h_other6_mac = FindVariableFeatures(h_other6_mac, selection.method = "disp", nfeatures = 500)
h_other6_mac = ScaleData(h_other6_mac)
h_other6_mac = RunPCA(h_other6_mac, features = VariableFeatures(object = h_other6_mac))
h_other6_mac = FindNeighbors(h_other6_mac, dims = 1:15, k.param = 30)
h_other6_mac = FindClusters(h_other6_mac, resolution = 0.8)
h_other6_mac = RunTSNE(h_other6_mac, features = VariableFeatures(object = h_other6_mac))

DimPlot(h_other6_mac, reduction = "tsne")

# save a copy
saveRDS(h_other6_mac, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/processed_data/lam.human.other6.mac.rds")

#####??????
# DEG between LAM and all the other macrophages
other6_lam_othermac = FindMarkers(h_other6_mac, ident.1 = "2" , ident.2 = c("0", "1"), min.pct = 0.25)

# save a copy
write.table(other6_lam_othermac, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/results/human.other6.lam.all_other_macrophage.deg.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

other6_lam_othermac[, "gene"] = rownames(other6_lam_othermac)
other6_lam_sig = dplyr::filter(other6_lam_othermac, p_val_adj < 0.001)

other6_upinlam = dplyr::filter(other6_lam_sig,avg_logFC>0)
other6_downinlam =  dplyr::filter(other6_lam_sig,avg_logFC<0)

other6_upinlam = other6_upinlam[order(other6_upinlam$avg_logFC, decreasing = T),]
other6_downinlam = other6_downinlam[order(other6_downinlam$avg_logFC, decreasing = T),]

head(other6_upinlam, 20)
head(other6_downinlam, 20)


# ID the LAM by plotting the markers
FeaturePlot(h_other6_mac, features = c("CHIT1", "TREM2", "LIPA", "SPP1"), reduction = "tsne")

# LAM should be cluster 2

# calculate the percentage of LAM out of macrophage for each individual
# cluster ID for each cell
other6_mac_cluster_id = as.data.frame(h_other6_mac$seurat_clusters)
other6_mac_cluster_id[,"cell"] = rownames(other6_mac_cluster_id)
colnames(other6_mac_cluster_id) = c("cluster", "cell")

other6_mac_cells = select(other6_mac_cluster_id, cell)%>%unlist()%>%as.vector()

dim(dplyr::filter(other6_mac_cluster_id, cluster =="2"))
# 102 LAM cells

mac_cluster2_cells = dplyr::filter(other6_mac_cluster_id, cluster =="2")%>%select(cell)%>%unlist()%>%as.vector()

# percentage of LAM out of macrophage
comput_perctg = function(x){
  n_of_mac = length(intersect(other6_mac_cells, x))
  n_of_lam = length(intersect(mac_cluster2_cells, x))
  return( c(n_of_lam, n_of_mac, n_of_lam/n_of_mac))
}

bmi23_pct = comput_perctg(bmi23_well)
bmi37_pct = comput_perctg(bmi37_well)
bmi38_pct = comput_perctg(bmi38_well)
bmi41_pct = comput_perctg(bmi41_well)
bmi43_pct = comput_perctg(bmi43_well)
bmi44_pct = comput_perctg(bmi44_well)

# the obese individual macrophage
fat_indi_mac_cells = as.data.frame(h_obese_mac$seurat_clusters)
fat_indi_lam = dplyr::filter(fat_indi_mac_cells, h_obese_mac$seurat_clusters == "3")
fat_lam_pc = nrow(fat_indi_lam)/nrow(fat_indi_mac_cells)

macrophage_m = matrix(c(23, 37, 38, 41, 43, 44, 46, 0.4230769, 0.1875, 0.2210526,
                  0.1388889, 0.4411765, 0.2758621, 0.1144997), ncol =2)

plot(macrophage_m[,1], macrophage_m[,2]*100, xlab = "bmi", ylab = "%LAM(macrophage)", pch = 19, col = "red")



# percentage of LAM out of CD45+ (immune) cells
other6_cd45cells_qc = names(human_other6_seurat_filter$seurat_clusters)

comput_perctg_cd45 = function(x){
  n_of_immune = length(intersect(other6_cd45cells_qc, x))
  n_of_lam = length(intersect(mac_cluster2_cells, x))
  return( c(n_of_lam, n_of_immune, n_of_lam/n_of_immune))
}

bmi23_pct_cd45 = comput_perctg_cd45(bmi23_well)
bmi37_pct_cd45 = comput_perctg_cd45(bmi37_well)
bmi38_pct_cd45 = comput_perctg_cd45(bmi38_well)
bmi41_pct_cd45 = comput_perctg_cd45(bmi41_well)
bmi43_pct_cd45 = comput_perctg_cd45(bmi43_well)
bmi44_pct_cd45 = comput_perctg_cd45(bmi44_well)

cd45_m = matrix(c(23, 37, 38, 41, 43, 44, 9.369676e-03, 0.01256983, 0.03825137,
                     0.01698754, 0.03989362, 0.02332362), ncol =2)
colnames(cd45_m) = c("bmi", "%LAM(CD45+)")

plot(cd45_m[,1], cd45_m[,2]*100, xlab = "bmi", ylab = "%LAM(CD45+)", pch = 19, col = "red")
