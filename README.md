# alzheimer_2019

##


## 1.single-cell RNA-seq identifies disease-associated microglial (DAM)

Keren-Shaul, et al. 2017. Cell. (GSE98969 single-cell RNA-seq)

##

1) obtain the processed sc RNA-seq data (read counts) from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98969/suppl/
##
2) obtain mouse reference genome from GENCODE https://www.gencodegenes.org/mouse/
##
3) parse the reference genome, select annotation only for genes. Get the length of each gene

  $python parse.ref.py

  --> gencode.vM21.annotation.gene.txt
##
4) merge the AD and WT mice (3 each)

  --> dam.ad_wt.umi.merged.rds
##
5) clean the data

   left with 8284 cells and 10976 genes
  --> dam.ad_wt.umi.merged.clean2.rds
##
6) clustering analysis

   $ seurat.v1.R
   
  --> dam_seurat_v1.rds
  
##
7) similar analysis of single-cell RNA-seq from Trem2-/- mouse

   $ dam.trem.R
##
8) explored Linnorm and Imputation (DrImpute)

   $ impute.trial.R
   $ compare.imputation.R


##


## 2.Injury-responsive microglial (IRM)

Timothy Hammond, et al. 2019. Immunity. (GSE121654 single-cell RNA-seq)

##

The data cleaning steps are as usual

However, the data has batch effect while batch information can't be found

Manually split clusters into batches, then use Seurat data integration function to integrate the batch data

$ irm.R

##

Focus on the gene candidates that are differentially expressed in the diseased cluster

$ gene.candidates.R

$ compare.genelist.R

## 3.ROSMAP scRNA-seq analysis

Hansruedi Mathys, et al. 2019. Nature. (limited access to the dataset)

##
All the analysis are in one script

$ human.ad.cluster.R

General steps:

1)

2)

A separate file for subcluster and categorical trait association test

$ categorical.test.Rmd

