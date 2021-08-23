#!/bin/bash
############################
# Single-Cell: Clustering ##
############################

# Wrapper for script seurat/dimReduction_v#.R

# Call this script as:
# . /path/to/script/thisScript.sh

COMDIR=/mnt/BioAdHoc/Groups/vd-vijay/cramirez
##### Default parameters ####
# directories need / in the end except OUTDIR
SOURCE_FILE=${COMDIR}/hayley/raw/sever_asthma3/COUNTS_hg19/AGGR/ovCD8_aggr_mapped # 10X folder PATH (before outs), CSV or Seurat object (same as MYCELLSF)
ANNOTF=/mnt/BioHome/ciro/asthma_pjs/info/sever_asthma3batches_annot.csv #=Search # annotation/meta.data file; it can be a table with the libraries and their groups Library, Condition_1, Condition_2, etc
OUTDIR=${COMDIR}/asthma_pjs/results/clust_seurat #=/mnt/BioScratch/cramirez # dont put "/" at the end
PJNAME=restim_cd8cd3 #=10XData # project name
# Filtering: file with cell barcodes to take; or mycolumn~myclass; an expression column == class | column2 == class
# or a "list(c('column1','class2'), c('column14','-class1'))"
# CELLSF #=ningun_file


#### Workflow ####
WFLOW=$1 # you can give this from outside
WFLOW=${WFLOW:-1111} # string with loical gates for: filtering, pca, dim. reduction and markers


# MTPATTERN #="c('^MT-', '^mt-', '^M-', '^m-', '^hg19_MT-', '^mm10_mt-')" # mitochondrial pattern find in genes
# Normalisation type, number of features and percentage of variance separated by underscores
# example LogNormalize_vst2000_30; LogNormalize, vst highly variable genes selection and 30% of variance explained
NORM_TYPE=LogNormalize_2000_30p
# SUBS_NAMES #="c('nGene', 'nUMI', 'percent.mito')" # variable to filter
LOW_FILT_CELLS="c(200, -Inf, -Inf)" # low bound of previous varaibles, -2 for default and -1 for boundaries
HIGH_FILT_CELLS="c(6000, 30000, 0.15)" # upper bound
# REGRESS_VAR #=NULL # variables to regress out (alredy including nUMI and percent.mito)
X_VARGS="c(low = 0.01, high = 8)" # mean cutoffs; only for mvp
# Y_VARGS #="c(low = 1, high = Inf)" # dispersion cutoffs; only for mvp

PCS_COMP=40 # number of PCs to calculate
# REDTYPE #="c('pca', 'tsne', 'umap')" # reduction types to use
CHOOSPC=0 # 0 choose based on elbow, 1 max elbow or significance, > 1 that number of PCs
PXITIES="c(100)" #="c(5, 30, 50, 75)" # Perplexities to revise
RES_EXP="c(0.4, 0.6)" #="c(0.4, 0.6, 1.2)" # resolutions to calculate
N_NEIS="c(15)" #="c(5, 10, 15, 30)" # No. neighbors'

# DE_TEST #=MAST # type of Diff. Expr. test
# NTOPG #=16 # plot no. of top genes
# LOGFC #=0.25 # log fold Change cutoff
# QVAL #=0.05 # adjusted p-value cutoff
GCOLOUR=/mnt/BioHome/ciro/asthma_pjs/info/sever_asthma_colors.csv #=ningun_file # file for colours
# VISF #=data # type of data to Visualise and add stats to markers
#####

## Only resting - singlets
ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/restim_cd8cd3/clustering/zetInfo/demuxlet_annotation_resting.Rdata
CELLSF=ningun_file
PJNAME=resting_cd8_15p_sng
NORM_TYPE=LogNormalize_2000_15p
CHOOSPC=13 # 18 (elbow)

# ## Only resting and regress donor - singlets
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/restim_cd8cd3/clustering/zetInfo/demuxlet_annotation_resting.Rdata
# CELLSF=ningun_file
# PJNAME=resting_cd8_15p_sng_reg
# NORM_TYPE=LogNormalize_2000_15p
# CHOOSPC=15 # 20 (elbow)
# REGRESS_VAR=orig.donor

# ## Only resting - singlets - donor and macro removed
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng/summary/donor_and_macro_removed.Rdata
# PJNAME=resting_cd8_15p_sng_noMacro_ds11
# NORM_TYPE=LogNormalize_2000_15p
# CHOOSPC=19 # 22 (elbow)
# REGRESS_VAR=NULL

# ## Only resting - singlets - donor and macro removed - regress donor
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng/summary/donor_and_macro_removed.Rdata
# PJNAME=resting_cd8_15p_sng_noMacro_ds11_reg
# NORM_TYPE=LogNormalize_2000_15p
# CHOOSPC=23 # 16 (elbow)
# REGRESS_VAR=orig.donor

## Only resting - singlets - donor and macro removed
ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng/summary/macro_removed.Rdata
PJNAME=resting_cd8_15p_sng_noMacro
NORM_TYPE=LogNormalize_2000_15p
CHOOSPC=0 #  (elbow)
REGRESS_VAR=NULL

# ## Only resting - singlets - donor and macro removed
# ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng/summary/macro_removed.Rdata
# PJNAME=resting_cd8_15p_sng_noMacro_reg
# NORM_TYPE=LogNormalize_2000_15p
# CHOOSPC=0 #  (elbow)
# REGRESS_VAR=orig.donor

## Only resting - singlets - donor and macro removed
ANNOTF=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng_noMacro/clustering/zetInfo/metadata_18PCs_30Ks_0.06667JD.RData
PJNAME=resting_cd8_15p_sng_nomacro_x8
NORM_TYPE=LogNormalize_2000_15p
CHOOSPC=25 # 19 (elbow)
REGRESS_VAR=NULL
CELLSF="list(c('RNA_snn_res.0.6', '-8'))"

WT=30:00:00
MEM=40gb

# # Sequennce of filtering
# restim_cd8cd3
# resting_cd8_15p_sng
# resting_cd8_15p_sng_noMacro
# resting_cd8_15p_sng_nomacro_x8

. /mnt/BioHome/ciro/scripts/seurat/jobConstSeu.sh
