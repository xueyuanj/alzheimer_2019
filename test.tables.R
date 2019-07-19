library(dplyr)

setwd("C:/Users/esi17551/human_ad")

# AD versus no AD

interestlist = c("AD.no-AD", "early-AD.late-AD", "early-AD.no-AD", "late-AD.no-AD", "mic_subcluster")

# AD progression: min.pct = 0.25, logfc.threshold = 0.25
# microglia subpopulation: min.pct = 0.25, logfc.threshold = 0.5

showresult = function(p){
  adnoad = list.files(pattern = p)
  print(adnoad)
  print(p)
  for (i in 1:length(adnoad)) {
    thismarkert = read.csv(adnoad[i], sep = "\t")
    print(adnoad[i])
    
    # adj.p < 0.01
    adjusted = thismarkert[which(thismarkert$p_val_adj < 0.01),]
    
    # up  
    print("up-regulated")
    print(na.omit(adjusted[which(thismarkert$avg_logFC > 0),]))
    
    # down
    print("down-regulated")
    print(na.omit(adjusted[which(thismarkert$avg_logFC < 0),]))
  }
}

sapply(interestlist, showresult)




