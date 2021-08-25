#!/share/apps/R/3.5/bin/R

############
# Figure 2 #
############

# This script will create the dimentional reduction plots for figure 2
# of the airways cd4 paper

libs = '~/R/newer_packs_library/3.5'
if(dir.exists(libs)).libPaths(libs)
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(cowplot)
  library(RcppParallel)
  library(uwot)
  library(factoextra)
  library(patchwork)
}); theme_set(theme_cowplot())

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/file_reading.R')
# readfile
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/filters.R')
# filters_subset_df
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/utilities.R')
# is.file.finished, v2cols
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/stats_var_features.R')
# getMostVariableGenes
source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plots_correlation.R")
# nclust_optimal
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plots.R')
# custom_heatmap
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/stats_summary_table.R')
# stats_summary_table (called within custom_heatmap)
source("https://raw.githubusercontent.com/vijaybioinfo/ASTHMA_AIRWAYS_2021/main/figures/dim_reduction_fun.R")
# dimred_hvg

fig_dir = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/results/figures_cd4'
if(dir.exists(fig_dir)) setwd(fig_dir)

# Global variables
colsname <- "/home/ciro/asthma_airways/info/airways_global_colours.csv"
grcols <- read.csv(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1);

### 2.X Dim reduction plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fannot = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/raw/mapping_2020_08_27/metadata_filtered_airways_varsubset.csv"
edata_counts_f = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/raw/mapping_2020_08_27/raw.Rdata"
mdata_init <- readfile(fannot, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
edata_counts_init <- as.matrix(readfile(edata_counts_f, check.names = FALSE, row.names = 1))

# Batch corrected - resting
pcs = c(7:9); pxs = c(10, 12)
edata_norm_f = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/raw/batchcor_cd4/combat_corrected_20517genes_231samples_Sex-plateCorrected.Rdata"
output_dir_i = 'umap_resting_batch_corr'
sample_filter_i = filters_subset_df(list(c("Stim", "Unstim"), c("Cell_type", "-CD8-TRM", "-CD4-TFH")), mdata_init, v = T)
# Batch corrected - stim
pcs = c(7:9); pxs = c(10, 12)
edata_norm_f = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/raw/batchcor_cd4/combat_corrected_20517genes_231samples_Sex-plateCorrected.Rdata"
output_dir_i = 'umap_stim_batch_corr'
sample_filter_i = filters_subset_df(list(c("Stim", "Stim"), c("Cell_type", "-CD8-TRM", "-CD4-TFH")), mdata_init, v = T)

edata_norm_init <- readfile(edata_norm_f)
# mdata$Cell_type <- factor(mdata$Cell_type, c("CD4-TRM", "CD4-TRM-like", "CD4-Teff", "CD4-Treg"))

dimred_hvg(
  mdata = mdata_init,
  edata_counts = edata_counts_init,
  edata_norm = edata_norm_init,
  pc_n = pxs,
  perplex_or_neigh = pcs,
  showgenes = c("ITGAE", "ITGA1", "FOXP3"),
  cnames = c("Stim", "Disease", "Cell_type", "Sex"),
  output_dir = output_dir_i,
  samples_filter = sample_filter_i,
  features_filter = rownames(edata_norm_init),
  take_sig = FALSE,
  padjthr = 0.05,
  fitv = 0.5,
  top_n = 800,
  couls = grcols,
  show_rownames = FALSE
)
