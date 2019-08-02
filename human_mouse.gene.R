library(dplyr)

setwd("~/Documents/2019_summer_eisai/progress_report/20190729")

# DAM versus homeostatic microglia
dam3to1 = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_1.txt")
colnames(dam3to1) = c("gene", "p_val", "avg_log_FC", "pct.1", "pct.2", "p_val_adj")

# Trem2-/- dataset, stage 2 (Trem2 dependent) versus stage 1 (Trem2 independent)
trem3to2 = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_1.txt")
colnames(trem3to2) = c("gene", "p_val", "avg_log_FC", "pct.1", "pct.2", "p_val_adj")

# ROSMAP Mic1
mic1 = read.table("~/Documents/2019_summer_eisai/rosmap.candidates/mic_subcluster.1.txt", header = T)
mic1[,"gene"] = rownames(mic1)

updam3to1 = dplyr::filter(dam3to1, avg_log_FC > 0.25 & p_val_adj < 0.05)
uptrem3to2 = dplyr::filter(trem3to2, avg_log_FC > 0.25 & p_val_adj < 0.05)
upmic1 = dplyr::filter(mic1, avg_logFC > 0.25 & p_val_adj < 0.05)

updam3to1$gene = toupper(updam3to1$gene); uptrem3to2$gene = toupper(uptrem3to2$gene)

# test the significance of overlap between human and mouse list
# human 17775 genes, mouse 10976 genes, the total genome = 10976
dampart_m = matrix(c(10976-length(union(updam3to1$gene, upmic1$gene)), length(setdiff(updam3to1$gene, upmic1$gene)), 
                     length (setdiff(upmic1$gene, updam3to1$gene)), length(intersect(updam3to1$gene, upmic1$gene)) ),
                   nrow = 2)

fisher.test(dampart_m, alternative = "greater")

# show the gene overlap
print(intersect(upmic1$gene,updam3to1$gene))

trempart_m = matrix(c(10976-length(union(uptrem3to2$gene, upmic1$gene)), length(setdiff(uptrem3to2$gene, upmic1$gene)), 
                     length (setdiff(upmic1$gene, uptrem3to2$gene)), length(intersect(uptrem3to2$gene, upmic1$gene)) ),
                   nrow = 2)

fisher.test(trempart_m, alternative = "greater")


print(intersect(upmic1$gene, uptrem3to2$gene))
# print(intersect(updam3to1$gene, uptrem3to2$gene))


# check the overlap between "homeostatic" microglia in human and mouse
# due to the lack of intepretation of human microglia cluster 0, the DEG were obtained
# using all other three clusters as outgroup
# therefore, the mouse ones also used the two other groups as outgroup, to make it consistent
mic0 = read.table("~/Documents/2019_summer_eisai/rosmap.candidates/mic_subcluster.0.txt", header = T)
mic0[,"gene"] = rownames(mic0)
mic0sig = dplyr::filter(mic0, avg_logFC > 0.25 & p_val_adj < 0.05)


damhomeo = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_1_to_23.txt")
colnames(damhomeo) = c("gene", "p_val", "avg_log_FC", "pct.1", "pct.2", "p_val_adj")
damhomeosig = dplyr::filter(damhomeo, avg_log_FC > 0.25 & p_val_adj < 0.05)
damhomeosig$gene = toupper(damhomeosig$gene)


dampart_homeo = matrix(c(10976-length(union(damhomeosig$gene, mic0sig$gene)), length(setdiff(damhomeosig$gene, mic0sig$gene)), 
                     length (setdiff(mic0sig$gene, damhomeosig$gene)), length(intersect(damhomeosig$gene, mic0sig$gene)) ),
                   nrow = 2)

fisher.test(dampart_homeo, alternative = "greater")
print(intersect(mic0sig$gene, damhomeosig$gene))

tremhomeo = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_homeo_to_12.txt")
colnames(tremhomeo) = c("gene", "p_val", "avg_log_FC", "pct.1", "pct.2", "p_val_adj")
tremhomeosig =  dplyr::filter(tremhomeo, avg_log_FC > 0.25 & p_val_adj < 0.05)
tremhomeosig$gene = toupper(tremhomeosig$gene)

trempart_homeo = matrix(c(10976-length(union(tremhomeosig$gene, mic0sig$gene)), length(setdiff(tremhomeosig$gene, mic0sig$gene)), 
                         length (setdiff(mic0sig$gene, tremhomeosig$gene)), length(intersect(tremhomeosig$gene, mic0sig$gene)) ),
                       nrow = 2)

fisher.test(trempart_homeo, alternative = "greater")

print(intersect(mic0sig$gene, tremhomeosig$gene))

