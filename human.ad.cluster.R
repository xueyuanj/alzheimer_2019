library(boot)
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
# 1920 cells

# check overlap
length(intersect(micglia_c, right_mic))
# 1920 cells
# all the microglia in the paper

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
newplot

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


# identify DEG in each cell type
# ad vs noad

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
  clustercells = subset(seurat.filtered, idents = thiscluster)
  
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

# subcluster of the microglia
microgliacluster = subset (seurat.filtered, idents = "microglia")
microgliacluster = NormalizeData(microgliacluster)
microgliacluster = FindVariableFeatures(microgliacluster, selection.method = "disp", nfeatures = 3188)
microgliacluster = ScaleData(microgliacluster, features = allgenes)

# pca
microgliacluster = RunPCA(microgliacluster, features = VariableFeatures(object = microgliacluster))
microgliacluster = FindNeighbors(microgliacluster, dims = 1:50, k.param = 30)
microgliacluster = FindClusters(microgliacluster, resolution = 0.3)
microgliacluster = RunTSNE(microgliacluster, features = VariableFeatures(object = microgliacluster))
microglia_tsne = DimPlot(microgliacluster, reduction = "tsne")
ggsave(microglia_tsne, filename = "rosmap.tsne.microglia.sub.v4.jpeg",width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# save a copy of the microglia subcluster data
saveRDS(microgliacluster, file = "rosmap.microglia.subcluster.rds")
microgliacluster = readRDS("rosmap.microglia.subcluster.rds")

microglia_cluster_id = as.data.frame(microgliacluster$seurat_clusters)
microglia_cluster_id[,"cell"] = rownames(microglia_cluster_id)
colnames(microglia_cluster_id) = c("cluster", "TAG")

# check the original subcluster information
paper_mic_sub = dplyr::filter(paper_tsne, broad.cell.type == "Mic")
compare_sub = merge(microglia_cluster_id, paper_mic_sub, by = "TAG")

compare_sub = merge(microglia_cluster_id, paper_mic_sub, by = "TAG")
table(compare_sub[, c("cluster", "Subcluster")])

#    Mic0  Mic1  Mic2  Mic3
# 0  1122  225   2     10
# 1  13    256   1     0
# 2  2     5     165   1
# 3  4     23    1     90

# plot the gene markers in the paper
papermic1_marker = c("FTL", "RPLP2", "C1QC", "FTH1", "SPP1", "RPLP1", "SLC11A1", "RPL28", "RPS20", "TMEM163")
papermic2_marker = c("SLC38A1", "STAT4", "PRF1", "FYN", "CD247", "SKAP1", "THEMIS", "NKG7", "BCL11B", "ITK")

mic1marker_plot = FeaturePlot(microgliacluster, features = papermic1_marker, reduction = "tsne", ncol = 3)
ggsave(mic1marker_plot, filename = "mic1_feature.tsne.jpeg", width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

mic2marker_plot = FeaturePlot(microgliacluster, features = papermic2_marker, reduction = "tsne", ncol = 3)
ggsave(mic2marker_plot, filename = "mic2_feature.tsne.jpeg", width = 10, height = 8, dpi = 150, units = "in", device='jpeg')

# run tests to identify DEG in one cluster against all others 
fourcluster = as.character(seq(0,3))

#c_pairwise = combn(fourcluster, 2)
# 
# for (i in 1:ncol(c_pairwise)) {
#   c1 = c_pairwise[1, i]; c2 = c_pairwise[2, i]
#   markers = FindMarkers(microgliacluster, ident.1 = c1 , ident.2 = c2, min.pct = 0.25, logfc.threshold = 0.5)
#   tablename = paste( "mic_subcluster", c1, c2,"txt", sep = ".")
#   write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
# }

for (i in 1:length(fourcluster)) {
  c1 = fourcluster[i]; c2 = fourcluster[-i]
  
  markers = FindMarkers(microgliacluster, ident.1 = c1 , ident.2 = c2, min.pct = 0.25, logfc.threshold = 0.5)
  tablename = paste( "mic_subcluster", c1, "txt", sep = ".")
  write.table(markers, file = tablename, quote = F, row.names = T, col.names = T, sep = "\t")
}


# first merge all the available information
cell_patient = paper_tsne[,c("TAG", "projid")]
cluster_patient = merge(cell_patient, microglia_cluster_id, by="TAG")
cell_trait = merge(cluster_patient, mergeid, by= "projid")

# add APOE genotype data
rosmap_clinical = read.csv("ROSMAP_Clinical_2019-05-03.csv")
apoedf = rosmap_clinical[,c("projid", "apoe_genotype")]

cell_trait = merge(cell_trait, apoedf, by = "projid")

# save a copy of this file
saveRDS(cell_trait, file = "microglia.subcluster.trait.rds")

cell_trait = readRDS("microglia.subcluster.trait.rds")

cat_trait = c("amyloid.group", "braaksc", "cogdx", "msex", "apoe_genotype")
quan_trait = c("gpath", "nft", "plaq_n", "amyloid", "cogn_global_lv")

# fisher's exact test for categorical traits

cat_p_correct=sapply(cat_trait, function(x){
  testtable = as.matrix(table(cell_trait[, c("cluster", x)]))
  print(testtable)
  ncluster = nrow(testtable); ntrait = ncol(testtable)
  pmatrix = as.data.frame(as.matrix(NA, ncol = ntrait, nrow=ncluster))
  for (i in 1:ncluster) {
    for (j in 1:ntrait) {
      
      fishertable = as.data.frame(as.matrix(NA, ncol=2, nrow=2))
      fishertable[1,1] = testtable[i, j]
      fishertable[1,2] = rowSums(as.matrix(testtable[,-j]))[i]
      fishertable[2,1] = colSums(as.matrix(testtable[-i,]))[j]
      fishertable[2,2] = sum(testtable[-i,-j])
      pvalue = fisher.test(fishertable)
      pmatrix[i,j] = pvalue
    }
  }
  fdr = apply(pmatrix, 2, function(x){
    p.adjust(x, method = "bonferroni", n=4)
  })
  colnames(fdr) = colnames(testtable)
  rownames(fdr) = rownames(testtable)
  print(fdr)
  return(fdr)
})

# save a copy
sink("trait_subpop.cat.association.txt")
print(cat_p_correct)
sink()


#write.table(cat_p_correct, file = "trait_subpop.cat.association.txt", quote = F, col.names = T, sep = "\t")

# boostrap for microglia subclusters and traits
cluster_ids = c("0", "1", "2", "3")

p_uncorrect =  sapply(quan_trait, function(x){
  testtable = cell_trait[, c("cluster", x)]
  sapply(cluster_ids, function(y){
    cluster_t = dplyr::filter(testtable, cluster == y)
    size = nrow(cluster_t)
   
    # get the distribution of bootrapping
    nulldistr = replicate(10000, {
      onesample = sample(testtable[,x], size, replace = F)
      mean(onesample)
    })
    # trait mean for this subcluster
    clustermean = mean(cluster_t[,x])
    
    # calculate z-score
    pop_sd = sd(nulldistr)*sqrt((length(nulldistr)-1)/(length(nulldistr)))
    pop_mean = mean(nulldistr)
    z = (clustermean - pop_mean) / pop_sd
    
    
    # p-value
    p_raw = 1 -pnorm(z)
    print(p_raw)
    return(c(z,p_raw))
  } )
} )
print(p_uncorrect)

# save a copy of the results
write.table(p_uncorrect, file = "trait_subpop.quan.association.txt", quote = F, col.names = T, sep = "\t")

# 
