# 2 versus 1
dam2to1 = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_2_to_1.txt")
trem2to1 = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_1_to_homeo.txt")

damup2genes = dam2to1[,1]; tremup2genes = trem2to1[,1]

length(damup2genes); length(tremup2genes)

common2to1 = intersect(damup2genes, tremup2genes)

common2to1

dam2uniq = damup2genes[!(damup2genes%in%common2to1)]
trem2uniq = tremup2genes[!(tremup2genes%in%common2to1)]


# 3 versus 1
dam3to1 = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_1.txt")
trem3to1 = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo.txt")

damup3genes = dam3to1[which(dam3to1$V6 < 0.05 & dam3to1$V3 > 0.4),1]
tremup3genes = trem3to1[which(trem3to1$V6 < 0.05 & trem3to1$V3 > 0.4),1]

length(damup3genes); length(tremup3genes)

common3to1 = intersect(damup3genes, tremup3genes)

common3to1

dam3uniq = damup3genes[!(damup3genes%in%common3to1)]
trem3uniq = tremup3genes[!(tremup3genes%in%common3to1)]

# 3 versus 2
dam3to2 = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_2.txt")
trem3to2 = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_1.txt")

damup3genes = dam3to2[which(dam3to2$V6 < 0.05),1] 
# change to log fold change 0.2 result in only 22 genes

tremup3genes = trem3to2[which(trem3to2$V6 < 0.05),1]

length(damup3genes); length(tremup3genes)

common3to2 = intersect(damup3genes, tremup3genes)

common3to2

dam3uniq = damup3genes[!(damup3genes%in%common3to2)]
trem3uniq = tremup3genes[!(tremup3genes%in%common3to2)]






# 3 versus 12
# dam3to12 = read.table("~/Documents/2019_summer_eisai/results/dam.genelist.up_3_to_12.txt")
# trem3to12 = read.table("~/Documents/2019_summer_eisai/results/trem2.genelist.up_2_to_homeo1.txt")
# 
# damup3genes = dam3to12[,1]; tremup3genes = trem3to12[,1]
# 
# common3to12 = intersect(damup3genes, tremup3genes)
# dam3uniq = damup3genes[!(damup3genes%in%common3to12)]
# trem3uniq = tremup3genes[!(tremup3genes%in%common3to12)]





