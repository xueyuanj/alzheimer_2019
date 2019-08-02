library(dplyr)
library(Linnorm)
library(RColorBrewer)
library(scater)
library(SC3)
library(Seurat)
library(SingleCellExperiment)
library(stringr)
setwd("~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW")

ad_wt_table = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.ad_wt.umi.merged.rds")

# get the design matrix
design_matrix = read.csv("GSE98969_experimental_design_f.txt",sep = "\t", header = T)
# fix the inconsistent format in Batch_desc column
design_matrix$Batch_desc = gsub(" ", "_", design_matrix$Batch_desc)
# check the number of cells per well
table(design_matrix$Number_of_cells)

# to re-create figure one, select the 6 six-month old AD and WT mouses (there is a mouse 4 in AD)
admouseinfor = dplyr::filter (design_matrix, Mouse_ID=="5XFAD")%>%dplyr::filter(grepl("^AD6m_mouse", Batch_desc))%>%dplyr::filter(!grepl("mouse4", Batch_desc))
wtmouseinfor = dplyr::filter (design_matrix, Mouse_ID=="C57BL/6")%>%dplyr::filter(grepl("^WT6m_mouse", Batch_desc))

cell_well = dplyr::filter(design_matrix, Number_of_cells == "1")%>%select(Well_ID)%>%unlist()
ad_wt_table = ad_wt_table[,(which(colnames(ad_wt_table)%in%cell_well))]

# update design matrix, keep only the non-empty wells
ori_ad_wt_matrix = rbind(admouseinfor, wtmouseinfor)
clean_ad_wt_matrix = ori_ad_wt_matrix[(which(ori_ad_wt_matrix$Well_ID%in%cell_well)),]

# use scater to do filtering and plot QC plots
scater_qc = SingleCellExperiment(assays = list(counts = as.matrix(ad_wt_table)), colData = clean_ad_wt_matrix)

# check ERCC and MT genes
total_genes = rownames(ad_wt_table) 

isSpike(scater_qc, "ERCC") = grepl("^ERCC-", total_genes)
isSpike(scater_qc, "MT") = grepl("^mt-", total_genes)

# plot the percentage of ERCC and mt genes. 
# filter out the ones with > 10% mitochondrial genes

scater_qc = calculateQCMetrics(scater_qc,feature_controls = list(ERCC = isSpike(scater_qc, "ERCC"), 
                                                                 MT = isSpike(scater_qc, "MT")))
plotColData(scater_qc,x = "total_features_by_counts",y = "pct_counts_MT",colour = "Amp_batch_ID")
plotColData(scater_qc, x = "total_features_by_counts",y = "pct_counts_ERCC",colour = "Amp_batch_ID")

# remove ERCC percentage > 75%

erccfilter = scater_qc$pct_counts_ERCC< 50
#table(erccfilter)
#erccfilter
#FALSE  TRUE 
#1104 11759 

# keep cells with 500 < UMI < 2500
umi_lower = scater_qc$total_counts_endogenous > 500
umi_upper = scater_qc$total_counts_endogenous < 2500

#table(umi_lower); table(umi_upper)
#3875 cellls < 500, 274 cells > 2500

# keep cells with expressed genes > 300
feature_lower = scater_qc$total_features_by_counts > 300

scater_filter_by_cell = scater_qc[, erccfilter&umi_lower&umi_upper&feature_lower]
dim(scater_filter_by_cell)
