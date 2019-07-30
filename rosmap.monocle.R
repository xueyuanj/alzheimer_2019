library(dplyr)
library(monocle3)

seurat.filtered = readRDS("/public_genomic_data_downloads/czhang/sharon/rosmap.miss_end_per.rds")
microgliacluster = subset (seurat.filtered, idents = "microglia")
rosmap_m_count = as.matrix(GetAssayData(microgliacluster, slot = "counts"))

genemeta = cbind(rownames(rosmap_m_count), rownames(rosmap_m_count))
genemeta = as.data.frame(genemeta)
colnames(genemeta) = c("id", "gene_short_name")
rownames(genemeta) = genemeta$id

dfmeta = readRDS("/public_genomic_data_downloads/czhang/sharon/microglia.subcluster.trait.rds")
rownames(dfmeta) = dfmeta$TAG
colnames(dfmeta) = c("projid", "TAG", "mic_cluster", "Subject", "amyloid", "plaq_n", "nft", "tangles",
                     "cogn_global_lv", "gpath", "gpath_3neocort", "amyloid.group", "caa_4gp", "ceradsc",
                     "braaksc", "niareagansc", "cogdx", "msex", "pathology.group", "apoe_genotype")

# reorder the count matrix according to the rows in cell metadata

mic_count_input = rosmap_m_count[,rownames(dfmeta)]

microglia_monocle = new_cell_data_set(mic_count_input, cell_metadata = dfmeta, gene_metadata = genemeta)

# start the analysis
microglia_monocle = preprocess_cds(microglia_monocle, num_dim = 100)
microglia_monocle = reduce_dimension(microglia_monocle)

microglia_monocle = cluster_cells(microglia_monocle)
microglia_monocle = learn_graph(microglia_monocle)

# save a copy of the microglia results
saveRDS(microglia_monocle, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/rosmap.my.mic.monocle.rds")


# the microglia data is of low dimension, try the original dataset from the paper

seurat.filtered = readRDS("/public_genomic_data_downloads/czhang/sharon/rosmap.filtered.seurat.raw.rds")
allcell = as.matrix(GetAssayData(seurat.filtered, slot = "counts"))


paper_tsne = read.delim("/public_genomic_data_downloads/czhang/sharon/filtered_column_metadata.txt")
clinical_infor = read.csv("/public_genomic_data_downloads/czhang/sharon/supplement_3.csv")
idmap = read.csv("/public_genomic_data_downloads/czhang/sharon/id_mapping.csv")
uniqidmap = unique(idmap[,c("projid", "Subject")])
mergeid = merge(clinical_infor, uniqidmap, by = "Subject")
# gather the information for each cell
cell_paperinfor = merge(paper_tsne, mergeid, by = "projid")
rownames(cell_paperinfor) = cell_paperinfor$TAG
all_cell_meta = cell_paperinfor[colnames(allcell) ,]

all_cell_gene = cbind(rownames(allcell), rownames(allcell))
colnames(all_cell_gene) = c("id", "gene_short_name")
rownames(all_cell_gene) = all_cell_gene[,"id"]

# gather all information
rosmap_all = new_cell_data_set(allcell, cell_metadata = all_cell_meta, gene_metadata = all_cell_gene)

# perform trajectory analysis
rosmap_all = preprocess_cds(rosmap_all, num_dim = 100)
rosmap_all = reduce_dimension(rosmap_all)

rosmap_all = cluster_cells(rosmap_all)
rosmap_all = learn_graph(rosmap_all)

saveRDS(rosmap_all, file = "/public_genomic_data_downloads/czhang/sharon/fat_trem2/rosmap.paper_infor.monocle.rds")

microglia_monocle = readRDS("/public_genomic_data_downloads/czhang/sharon/fat_trem2/rosmap.my.mic.monocle.rds")
rosmap_all = readRDS("/public_genomic_data_downloads/czhang/sharon/fat_trem2/rosmap.paper_infor.monocle.rds")

alltraitsofinterest = c("broad.cell.type", "Subcluster", "amyloid.group", "ceradsc", "braaksc", "cogdx", "msex","pathology.group")

mictraitofinterest = c("mic_cluster", "amyloid.group", "ceradsc", "braaksc", "cogdx", "msex","pathology.group", "apoe_genotype")
# all cells
traitplot = function(x){
  plot_cells(rosmap_all, color_cells_by = x,label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE) 
}
lapply(alltraitsofinterest, traitplot)

plot_cells(rosmap_all, 
           label_cell_groups = F,
           show_trajectory_graph=FALSE,
           genes = c("FTL", "SPP1", "APOE", "TREM2")) 

# microglias
traitplot2 = function(x){
  plot_cells(microglia_monocle, color_cells_by = x,label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE) 
}
lapply(mictraitofinterest, traitplot2)

plot_cells(microglia_monocle, 
           label_cell_groups = F,
           show_trajectory_graph=FALSE,
           genes = c("FTL", "SPP1", "APOE", "TREM2"))



