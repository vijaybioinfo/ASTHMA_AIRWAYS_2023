#!/bin/bash
############################
# Single-Cell: Clustering ##
############################

# Wrapper for script seurat/dimReduction_v#.R

# Call this script as:
# . /path/to/script/thisScript.sh

##### Default parameters ####
# directories need / in the end except OUTDIR
# 10X folder PATH (before outs), CSV or Seurat object (same as MYCELLSF)
SOURCE_FILE=/home/ciro/large/hayley/raw/sever_asthma3/COUNTS_hg19/AGGR/CD4US_filt_mapped
#=Search # annotation/meta.data file; it can be a table with the libraries and their groups Library, Condition_1, Condition_2, etc
ANNOTF="/mnt/BioHome/ciro/asthma_pjs/info/airways_annot3batches.csv /home/ciro/large/asthma_pjs/results/demuxlet/sever_asthma/resting_cd4.rdata"
OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat #=/mnt/BioScratch/cramirez # dont put "/" at the end
PJNAME=resting_cd4 #=10XData # project name
# Filtering: file with cell barcodes to take; or mycolumn~myclass; an expression column == class | column2 == class
# or a "list(c('column1','class2'), c('column14','-class1'))"
# CELLSF #=ningun_file


#### Workflow ####
WFLOW=$1 # you can give this from outside
WFLOW=${WFLOW:-1111} # string with loical gates for: filtering, pca, dim. reduction and markers


# MTPATTERN #="c('^MT-', '^mt-', '^M-', '^m-', '^hg19_MT-', '^mm10_mt-')" # mitochondrial pattern find in genes
# Normalisation type, number of features and percentage of variance separated by underscores
# example LogNormalize_vst2000_30; LogNormalize, vst highly variable genes selection and 30% of variance explained
NORM_TYPE=LogNormalize_2000_15p
# SUBS_NAMES #="c('nGene', 'nUMI', 'percent.mito')" # variable to filter
LOW_FILT_CELLS="c(200, -Inf, -Inf)" # low bound of previous varaibles, -2 for default and -1 for boundaries
HIGH_FILT_CELLS="c(6000, 30000, 0.15)" # upper bound
# REGRESS_VAR #=NULL # variables to regress out (alredy including nUMI and percent.mito)
X_VARGS="c(low = 0.01, high = 8)" # mean cutoffs; only for mvp
# Y_VARGS #="c(low = 1, high = Inf)" # dispersion cutoffs; only for mvp

PCS_COMP=40 # number of PCs to calculate
# REDTYPE #="c('pca', 'tsne', 'umap')" # reduction types to use
CHOOSPC=0 # (elbow) 0 choose based on elbow, 1 max elbow or significance, > 1 that number of PCs
PXITIES="c(75)" #="c(5, 30, 50, 75)" # Perplexities to revise
RES_EXP="c(0.4, 0.6)" #="c(0.4, 0.6, 1.2)" # resolutions to calculate
N_NEIS="c(15)" #="c(5, 10, 15, 30)" # No. neighbors'

# DE_TEST #=MAST # type of Diff. Expr. test
# NTOPG #=16 # plot no. of top genes
# LOGFC #=0.25 # log fold Change cutoff
# QVAL #=0.05 # adjusted p-value cutoff
GCOLOUR=/mnt/BioHome/ciro/asthma_pjs/info/airways_global_colours.csv #=ningun_file # file for colours
# VISF #=data # type of data to Visualise and add stats to markers
#####

## Taking singlets, and removing 7 (low quality) and 9 (macros)
ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p/clustering/zetInfo/metadata_19PCs_30Ks_0.06667JD.RData
PJNAME=resting_cd4_15p_sng_nomacro_x7n9
CELLSF="list(c('orig.class', 'SNG'), c('RNA_snn_res.0.4', '-7', '-9'))"

## Taking singlets, and removing 4, 7, and 10 (low quality), and 9 (macros)
ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p/clustering/zetInfo/metadata_19PCs_30Ks_0.06667JD.RData
PJNAME=resting_cd4_15p_sng_nomacro_x4n7n9n10
CELLSF="list(c('orig.class', 'SNG'), c('RNA_snn_res.0.4', '-7', '-9', '-4', '-10'))"
CHOOSPC=20

# ## Subclustering: 0
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# CELLSF="list(c('RNA_snn_res.0.4', '0'))"
# PJNAME=cluster0
# NORM_TYPE=LogNormalize_2000_15p

# ## Subclustering: 1
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# CELLSF="list(c('RNA_snn_res.0.4', '1'))"
# PJNAME=cluster1
# NORM_TYPE=LogNormalize_2000_15p

# ## Subclustering: 3
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# CELLSF="list(c('RNA_snn_res.0.4', '3'))"
# PJNAME=cluster3
# NORM_TYPE=LogNormalize_2000_20p

# ## Subclustering: 4
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# CELLSF="list(c('RNA_snn_res.0.4', '4'))"
# PJNAME=cluster4
# NORM_TYPE=LogNormalize_2000_15p

# ## Subclustering: 0 and 1
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# CELLSF="list(c('RNA_snn_res.0.4', '1', '0'))"
# PJNAME=cluster0n1
# NORM_TYPE=LogNormalize_2000_15p
# RES_EXP="c(0.2, 0.4, 0.6)"

# ## Subclustering: 0 and 1 - imputation
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData
# OUTDIR=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9
# PJNAME=cluster0n1_imput
# NORM_TYPE=LogNormalize_2000_20p
# RES_EXP="c(0.2, 0.4, 0.6)"
# SOURCE_FILE=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/umi_matrix20PCs_30Ks_0_saver_estimate.rds

WT=30:00:00
MEM=60gb
# MFID=3552100

. /mnt/BioHome/ciro/scripts/seurat/jobConstSeu.sh
