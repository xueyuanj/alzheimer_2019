library(dplyr)
library(ggplot2)
library(Matrix)
library(Seurat)

# /anaconda2/bin/R

setwd("/public_genomic_data_downloads/czhang/sharon/")

seurat.filtered = readRDS("/public_genomic_data_downloads/czhang/sharon/rosmap.filtered.seurat.raw.rds")

allgenes = rownames(seurat.filtered)

# VlnPlot(seurat.filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seurat.filtered = NormalizeData(seurat.filtered)
seurat.filtered = FindVariableFeatures(seurat.filtered, selection.method = "disp", nfeatures = 3188)


# plot the top variable genes
seurat_top20 = head(VariableFeatures(seurat.filtered), 20)

var_plot1 = VariableFeaturePlot(seurat.filtered)
var_plot2 = LabelPoints(plot = var_plot1, points = seurat_top20, repel = T)

ggsave(var_plot2, filename = "rosmap.top_var_plot.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

#CombinePlots(plots = list(var_plot2))

seurat.filtered = ScaleData(seurat.filtered, features = allgenes)

# pca
seurat.filtered = RunPCA(seurat.filtered, features = VariableFeatures(object = seurat.filtered))
DimPlot(seurat.filtered, reduction = "pca")

# umap and tsne
seurat.filtered = FindNeighbors(seurat.filtered, dims = 1:50, k.param = 30)
seurat.filtered = FindClusters(seurat.filtered, resolution = 0.3)
#seurat.filtered = RunUMAP(seurat.filtered, features = VariableFeatures(object = seurat.filtered))
seurat.filtered = RunTSNE(seurat.filtered, features = VariableFeatures(object = seurat.filtered))

# save a copy of the plot
ori_tsne = DimPlot(seurat.filtered, reduction = "tsne")

ggsave(ori_tsne, filename = "rosmap.tsne.pre_cluster.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# plot the cluster markers
# excitatory neurons
exmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("SLC17A7", "CAMK2A", "NRGN"))
ggsave(exmarker, filename = "rosmap.tsne.ex_neuron.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# inhibitory neurons
inmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("SYT1", "SNAP25", "GAD1", "GAD2"))
ggsave(inmarker, filename = "rosmap.tsne.in_neuron.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# Oligodendrocytes
oligmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("MBP", "MOBP", "PLP1"))
ggsave(oligmarker, filename = "rosmap.tsne.olig.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# Oligodendrocyte progenitors
opcmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("PDGFRA", "VCAN", "CSPG4"))
ggsave(opcmarker, filename = "rosmap.tsne.opc.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# microglia
micmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("CD74", "CSF1R", "C3"))
ggsave(micmarker, filename = "rosmap.tsne.mic.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# astrocytes
astmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("AQP4", "GFAP", "AQP4", "GFAP"))
ggsave(astmarker, filename = "rosmap.tsne.ast.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# endothelial
endmarker = FeaturePlot(seurat.filtered, reduction = "tsne", features = c("FLT1", "CLDN5", "FLT1", "CLDN5"))
ggsave(endmarker, filename = "rosmap.tsne.end.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# due to the complex colors, plot each cluster, then collapse them

seurat.filtered = readRDS("/public_genomic_data_downloads/czhang/sharon/rosmap.filtered.seurat.tsne.rds")

ori_cluster_id = as.data.frame( seurat.filtered$seurat_clusters)
ori_cluster_id[,"cell"] = rownames(ori_cluster_id)
colnames(ori_cluster_id) = c("cluster", "cell")


for (i in 0:18) {
  
  # select each cluster
  thiscluster = dplyr::filter(ori_cluster_id, cluster == as.character(i))%>%select(cell)%>%unlist()%>%as.vector()
  
  # highlight those cells
  cplot = DimPlot(object = seurat.filtered, cells.highlight = thiscluster, reduction = "tsne")
 
   # save the plot
  plotname = paste(c("cluster", as.character(i), "tsne", "jpeg"), collapse = ".")
  
  ggsave(cplot, filename = plotname, width = 10, height = 8, dpi = 150, units = "in", device='jpeg')
  
  }


# verify the identity of microglia, how much they overlap with the results in the paper

# my microglia cells
micglia_c = dplyr::filter(ori_cluster_id, cluster == "11"|cluster == "18")%>%select(cell)%>%unlist()%>%as.vector()
# 1926 cells, about 3% of the total cells

# microglia cells from the paper

paper_tsne = read.delim("filtered_column_metadata.txt")
right_mic = dplyr::filter(paper_tsne, broad.cell.type == "Mic")%>%select(TAG)%>%unlist()%>%as.vector()
# 1926 cells

# check overlap
length(intersect(micglia_c, right_mic))
# 1920 cells
# all the microglia are clustered as microglia in the paper

# now reassign the cluster information, start with straightforward ones

# mic
seurat.filtered = SetIdent(object = seurat.filtered, cells = micglia_c, value = "microglia")

# olig
oligodendrocyte = dplyr::filter(ori_cluster_id, cluster == "0")%>%select(cell)%>%unlist()%>%as.vector() 
seurat.filtered = SetIdent(object = seurat.filtered, cells = oligodendrocyte, value = "oligodendrocyte")

# opc
oligoprogenitor = dplyr::filter(ori_cluster_id, cluster == "9")%>%select(cell)%>%unlist()%>%as.vector() 
seurat.filtered = SetIdent(object = seurat.filtered, cells = oligoprogenitor, value = "oligodendrocyte progenitor")

# ast
astrocyte = dplyr::filter(ori_cluster_id, cluster == "6")%>%select(cell)%>%unlist()%>%as.vector() 
seurat.filtered = SetIdent(object = seurat.filtered, cells = astrocyte, value = "astrocyte")

# nuerons
allnuerons = dplyr::filter(ori_cluster_id, cluster == "1"| cluster =="2"| cluster =="3"| cluster =="4"| cluster =="5"|
                           cluster =="7"| cluster =="8"| cluster == "10"| cluster =="12"| cluster =="13"| cluster == "14"| cluster =="15"| cluster =="16")%>%select(cell)%>%unlist()%>%as.vector() 

# need coordinate information
rosmap_embeds = as.data.frame(Embeddings(seurat.filtered[["tsne"]]))
rosmap_embeds[,"cell"] = rownames(rosmap_embeds)

inhibneurons_coor =  dplyr::filter(rosmap_embeds, (-40)<tSNE_1&tSNE_1< (-10)&(-25)<tSNE_2&tSNE_2<5)%>%select(cell)%>%unlist()%>%as.vector()
inhineuron_c = dplyr::filter(ori_cluster_id, cluster == "4"|cluster == "7"|cluster == "13"|cluster == "15"|cluster == "16")%>%select(cell)%>%unlist()%>%as.vector() 

inneurons = intersect(inhibneurons_coor, inhineuron_c)
seurat.filtered = SetIdent(object = seurat.filtered, cells = inneurons, value = "inhibitory neurons")


exneurons =setdiff(allnuerons, inneurons)
seurat.filtered = SetIdent(object = seurat.filtered, cells = exneurons, value = "excitatory neurons")

setwd("C:/Users/esi17551/human_ad")

seurat.filtered = readRDS("rosmap.miss_end_per.rds")

# endothelial
endopri_cluster = dplyr::filter(ori_cluster_id, cluster == "17")%>%select(cell)%>%unlist()%>%as.vector() 
endocoor =  dplyr::filter(rosmap_embeds, tSNE_2> (tSNE_1 + 60) )%>%select(cell)%>%unlist()%>%as.vector()
endothelial = intersect(endopri_cluster, endocoor)
seurat.filtered = SetIdent(object = seurat.filtered, cells = endothelial, value = "endothelial")

# pericyte
pericyte = setdiff(endopri_cluster, endothelial)
seurat.filtered = SetIdent(object = seurat.filtered, cells = pericyte, value = "pericyte")

# redo the plot with new cluster assignment
newplot = DimPlot(seurat.filtered, reduction = "tsne", label = F, pt.size = 0.5 ) 
ggsave(newplot, filename = "rosemap.tsne.major.cluster.jpeg", width = 10, height = 8, dpi = 150, units = "in", device='jpeg' )

# get the AD and no AD patients infor
clinical_infor = read.csv("supplement_3.csv")
idmap = read.csv("id_mapping.csv")
uniqidmap = unique(idmap[,c("projid", "Subject")])

# table with both subject and project id
mergeid = merge(clinical_infor, uniqidmap, by = "Subject")

earlyp = dplyr::filter(mergeid, pathology.group == "early-pathology")%>%select(projid)%>%unlist()%>%as.vector()
latep = dplyr::filter(mergeid, pathology.group == "late-pathology")%>%select(projid)%>%unlist()%>%as.vector()

adpatient = c(earlyp, latep)
noad =dplyr::filter(mergeid, pathology.group == "no-pathology")%>%select(projid)%>%unlist()%>%as.vector()

# get the cell ID for each group
adcells = dplyr::filter(paper_tsne, projid %in% adpatient)%>%select(TAG)%>%unlist()%>%as.vector()
noadcells = dplyr::filter(paper_tsne, projid %in% noad)%>%select(TAG)%>%unlist()%>%as.vector()


rosmap_cluster_id = as.data.frame(seurat.filtered$seurat_clusters)
rosmap_cluster_id[,"cell"] = rownames(rosmap_cluster_id)
colnames(rosmap_cluster_id) = c("cluster", "cell")

major6types = c("excitatory neurons", "inhibitory neurons", "astrocyte", "oligodendrocyte progenitor",
                "oligodendrocyte", "microglia")


for(i in 1:length(major6types)){
  thiscluster = major6types[i]
  
  # get the cell cluster
  clustercells = SubsetData(seurat.filtered, ident.use = thiscluster)
  
  # set AD and no-AD identity
  clustercells = SetIdent(object =clustercells, cells = adcells, value = "AD")
  clustercells = SetIdent(object =clustercells, cells = noadcells, value = "no-AD")
  
  markers = FindMarkers(clustercells, ident.1 = "AD" , ident.2 = "no-AD", min.pct = 0.25, logfc.threshold = 0.25)
  
  tablename = paste( "AD", "no-AD", thiscluster, "txt", sep = ".")
  write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
}


# early vs no ad
# late vs early
# late vs no ad
earlyadcells = dplyr::filter(paper_tsne, projid %in% earlyp)%>%select(TAG)%>%unlist()%>%as.vector()
lateadcells = dplyr::filter(paper_tsne, projid %in% latep)%>%select(TAG)%>%unlist()%>%as.vector()

for(i in 1:length(major6types)){
  thiscluster = major6types[i]
  
  # get the cell cluster
  clustercells = SubsetData(seurat.filtered, ident.use = thiscluster)
  
  # set AD and no-AD identity
  clustercells = SetIdent(object =clustercells, cells = earlyadcells, value = "early-AD")
  clustercells = SetIdent(object =clustercells, cells = lateadcells, value = "late-AD")
  clustercells = SetIdent(object =clustercells, cells = noadcells, value = "no-AD")
  
  # early vs no AD
  markers = FindMarkers(clustercells, ident.1 = "early-AD" , ident.2 = "no-AD", min.pct = 0.25, logfc.threshold = 0.25)
  tablename = paste( "early-AD", "no-AD", thiscluster, "txt", sep = ".")
  write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
  
  # late vs early AD
  markers = FindMarkers(clustercells, ident.1 = "early-AD" , ident.2 = "late-AD", min.pct = 0.25, logfc.threshold = 0.25)
  tablename = paste( "early-AD", "late-AD", thiscluster, "txt", sep = ".")
  write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
  
  # late vs no
  markers = FindMarkers(clustercells, ident.1 = "late-AD" , ident.2 = "no-AD", min.pct = 0.25, logfc.threshold = 0.25)
  tablename = paste( "late-AD", "no-AD", thiscluster, "txt", sep = ".")
  write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
  
}


