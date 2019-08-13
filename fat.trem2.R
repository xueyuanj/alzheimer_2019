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

upinmac3 =  dplyr::filter(mac1_3_sig,avg_logFC>0)
upinmac1 = dplyr::filter(mac1_3_sig,avg_logFC<0)

head(upinmac1, 20)
head(upinmac3, 20)

# download the processed human data (OAT)
human_oat = dplyr::filter(design_infor, `organ:ch1`=="OAT")
oat_urls = select(human_oat, supplementary_file_1)%>%unlist()%>%as.character()

# download the data
# sapply(oat_urls, function(x) download.file(x, unlist(strsplit(x, "/"))[9]))

# trem2 KO part

setwd("/public_genomic_data_downloads/czhang/sharon/fat_trem2/GSE128518_RAW")

trem2_ko_hfd = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 knock-out C57BL/6"&`treatment (diet):ch1` == "High-fat diet")
trem2_ko_nc = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 knock-out C57BL/6"&`treatment (diet):ch1` == "Normal chow")

litterm_hfd = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 wild-type C57BL/6 Littermates"&`treatment (diet):ch1` == "High-fat diet")
litterm_nc = dplyr::filter(m_design_infor, `strain:ch1` == "Trem2 wild-type C57BL/6 Littermates"&`treatment (diet):ch1` == "Normal chow")

hfd_trem_part = rbind(trem2_ko_hfd, litterm_hfd)

hfd_trem_supp = select(hfd_trem_part, supplementary_file_1)%>%unlist()%>%as.character() 
hfd_trem_f = sapply(hfd_trem_supp, function(x) {
  gzfile = unlist(strsplit(x, "/"))[9] 
  unzipfile = sub(".gz", "", gzfile)
},  USE.NAMES = F ) 


hfd_trem_df = lapply(hfd_trem_f, function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

hfd_trem_data = data.frame(hfd_trem_df)


# human part



