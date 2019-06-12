# alzheimer_2019

##


## 1.single-cell RNA-seq identifies disease-associated microglial (DAM)

Keren-Shaul, et al. 2017. Cell. (GSE98969 single-cell RNA-seq)

##

#Obtain the processed sc RNA-seq data (read counts) from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98969/suppl/

#Obtain mouse reference genome from GENCODE https://www.gencodegenes.org/mouse/

#Parse the reference genome, select annotation only for genes. Get the length of each gene

$python parse.ref.py

--> gencode.vM21.annotation.gene.txt

#merge the 97 samples, and calculate TPM
