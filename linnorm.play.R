library(dplyr)
library(Linnorm)

setwd("~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW")

# get the design matrix
design_matrix = read.csv("GSE98969_experimental_design_f.txt",sep = "\t", header = T)
# fix the inconsistent format in Batch_desc column
design_matrix$Batch_desc = gsub(" ", "_", design_matrix$Batch_desc)
# check the number of cells per well
table(design_matrix$Number_of_cells)

# to re-create figure one, select the 6 six-month old AD and WT mouses (there is a mouse 4 in AD)
admouseinfor = filter (design_matrix, Mouse_ID=="5XFAD")%>%filter(grepl("^AD6m_mouse", Batch_desc))%>%filter(!grepl("mouse4", Batch_desc))
wtmouseinfor = filter (design_matrix, Mouse_ID=="C57BL/6")%>%filter(grepl("^WT6m_mouse", Batch_desc))

# get the files using sample information
adfiles = unique(droplevels(admouseinfor)$Amp_batch_ID)
wtfiles = unique(droplevels(wtmouseinfor)$Amp_batch_ID)

# read the data and merge them into one data frame
adtables = lapply(as.character(adfiles), function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

admegatable = data.frame(adtables)

wttables = lapply(as.character(wtfiles), function(x){
  read.csv(list.files(pattern = x), sep="\t", header=TRUE)})

wtmegatable = data.frame(wttables)

ad_wt_table = data.frame(admegatable, wtmegatable)

# save a copy of this file
saveRDS(ad_wt_table, file = "./pbmc3k_final.rds")

# get the processed data (read count)
countfile = list.files(path = "~/Documents/2019_summer_eisai/raw_data/GSE98969_RAW", pattern = "GSM" )

testdf = read.csv(list.files(pattern = as.character(adfiles[1])), sep="\t", header=TRUE)

table(grepl("^ERCC-", rownames(testdf)))
table(rownames(testdf) %in% 
        c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
          "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
          "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
          "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
          "ENSG00000198840"))
# save a copy of the mega matrix
