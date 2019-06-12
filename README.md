# alzheimer_2019

##


## 1.single-cell RNA-seq identifies disease-associated microglial (DAM)

Keren-Shaul, et al. 2017. Cell. (GSE98969 single-cell RNA-seq)

##

1)Obtain the processed sc RNA-seq data (read counts) from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98969/suppl/

2)Obtain mouse reference genome from GENCODE https://www.gencodegenes.org/mouse/

3)Parse the reference genome, select annotation only for genes. Get the length of each gene

$python parse.ref.py

--> gencode.vM21.annotation.gene.txt

4)merge the 97 samples, and calculate TPM
