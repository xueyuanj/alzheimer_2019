library(Seurat)

ad_wt_seurat_filter = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.seurat.reassign_cluster.rds")
DimPlot(ad_wt_seurat_filter, reduction = "tsne")

# get the genes differentially expressed in the three microglia clusters
# up in 1 to 2
up1to2 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = "microglia 2", 
                     min.pct = 0.25, only.pos = T, logfc.threshold = 0.1)
write.table(up1to2, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 1 to 3
up1to3 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = "microglia 3",
                     min.pct = 0.25, only.pos = T, logfc.threshold = 0.1)
write.table(up1to3, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_3.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 1 to stage 2&3
up1to23 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 1", ident.2 = c("microglia 2", "microglia 3"), 
                      min.pct = 0.25, only.pos = T, logfc.threshold = 0.1)
write.table(up1to23, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_23.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 1
up2to1 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = "microglia 1", min.pct = 0.25,
                     only.pos = T, logfc.threshold = 0.1)
write.table(up2to1, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 3
up2to3 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = "microglia 3", min.pct = 0.25, 
                     only.pos = T, logfc.threshold = 0.1)
write.table(up2to3, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_3.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 2 to 1&3
up2to13 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 2" , ident.2 = c("microglia 1", "microglia 3"), min.pct = 0.25,
                      only.pos = T, logfc.threshold = 0.1)
write.table(up2to13, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_13.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 1
up3to1 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = "microglia 1", min.pct = 0.25,
                     only.pos = T, logfc.threshold = 0.1)
write.table(up3to1, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 2
up3to2 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = "microglia 2", min.pct = 0.25,
                     only.pos = T, logfc.threshold = 0.1)
write.table(up3to2, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in 3 to 1&2
up3to12 = FindMarkers(ad_wt_seurat_filter, ident.1 = "microglia 3" , ident.2 = c("microglia 1", "microglia 2"), min.pct = 0.25,
                      only.pos = T, logfc.threshold = 0.1)
write.table(up3to12, 
            file="~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_12.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)


####################################################################################################

trem2_seurat_filter = readRDS("~/Documents/2019_summer_eisai/processed_data/dam.trem2.seurat.optimal.rds")

# up in homeo to stage 1
uphomeoto1 = FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = "stage 1 DAM", min.pct = 0.25,
                         only.pos = T, logfc.threshold = 0.1)
write.table(uphomeoto1, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in homeo to stage 2
uphomeoto2 = FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = "stage 2 DAM", min.pct = 0.25,
                         only.pos = T, logfc.threshold = 0.1)
write.table(uphomeoto2, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in hemeo to stage 1 &2
uphomeoto12 = FindMarkers(trem2_seurat_filter, ident.1 = "homeostatic microglia", ident.2 = c("stage 1 DAM", "stage 2 DAM"), min.pct = 0.25,
                          only.pos = T, logfc.threshold = 0.1)
write.table(uphomeoto12, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_12.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to homeo
up1tohomeo = FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = "homeostatic microglia", min.pct = 0.25,
                         only.pos = T, logfc.threshold = 0.1)
write.table(up1tohomeo, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_homeo.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to stage 2
up1to2 = FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = "stage 2 DAM", min.pct = 0.25,
                     only.pos = T, logfc.threshold = 0.1)
write.table(up1to2, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 1 to homeo and stage 2
up1tohomeo2 = FindMarkers(trem2_seurat_filter, ident.1 = "stage 1 DAM" , ident.2 = c("homeostatic microglia", "stage 2 DAM"), min.pct = 0.25,
                          only.pos = T, logfc.threshold = 0.1)
write.table(up1tohomeo2, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_homeo2.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to homeo
up2tohomeo = FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = "homeostatic microglia", min.pct = 0.25,
                         only.pos = T, logfc.threshold = 0.1)
write.table(up2tohomeo, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to stage 1
up2to1 = FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = "stage 1 DAM", min.pct = 0.25,
                     only.pos = T, logfc.threshold = 0.1)
write.table(up2to1, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)

# up in stage 2 to homeo and stage 1
up2tohomeo1 = FindMarkers(trem2_seurat_filter, ident.1 = "stage 2 DAM" , ident.2 = c("homeostatic microglia", "stage 1 DAM"), min.pct = 0.25,
                          only.pos = T, logfc.threshold = 0.1)
write.table(up2tohomeo1, 
            file="~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo1.txt",
            quote = F, sep = "\t", col.names = F, row.names = T)




