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
})
theme_set(theme_cowplot())

source('/home/ciro/scripts/handy_functions/devel/file_reading.R')
source('/home/ciro/scripts/handy_functions/devel/filters.R')
source('/home/ciro/scripts/handy_functions/devel/utilities.R')

dirfig = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/final_figures_cd4'
dir.create(dirfig); setwd(dirfig)

# Global variables
colsname <- "/home/ciro/asthma_pjs/info/airways_global_colours.csv"
grcols <- read.csv(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1);

### 2.X Dim reduction plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Try UMAP with uwot
## Inputs
fannot = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/mapping_2020_08_27/metadata_filtered_airways_varsubset.csv"
fedataraw = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/mapping_2020_08_27/raw.Rdata"
fedata = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/mapping_2020_08_27/TPM.Rdata"
# 7PCs_12PX; 8PCs_10PX; 9PCs_10PX; 9PCs_12PX
pcs = c(2:3, 5, 7:9, 10)
pxs = c(7, 10, 12)
showgenes = c("ITGAE", "ITGA1", "FOXP3")
couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
cnames = c("Stim", "Disease", "Cell_type", "Plate", "Sex", "Sex_Disease", "Study_ID")

# Batch corrected - resting
pcs = c(7:9); pxs = c(10, 12)
fedata = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/batchcor_cd4/combat_corrected_20517genes_231samples_Sex-plateCorrected.Rdata"
outdir = "tsnes_resting_batch_corr"
outdir = 'umap_resting_batch_corr'
hvg_n = 800
hvg_n = 1313
selectss = list(c("Stim", "Unstim"), c("Cell_type", "-CD8-TRM", "-CD4-TFH"))

# Batch corrected - stim
pcs = c(7:9); pxs = c(10, 12)
fedata = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/batchcor_cd4/combat_corrected_20517genes_231samples_Sex-plateCorrected.Rdata"
outdir = 'tsnes_stim_batch_corr'
outdir = 'umap_stim_batch_corr'
hvg_n = 800
hvg_n = 2265
selectss = list(c("Stim", "Stim"), c("Cell_type", "-CD8-TRM", "-CD4-TFH"))
# outdir = 'umap_stim_batch_corr_xoutlier'

## Reading
annot <- readfile(fannot, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
if(grepl("outlier", outdir)) annot <- annot[!rownames(annot) %in% "071_NIHMA_233_KA_CD4_Teff_Stim", ]
expdata_counts <- as.matrix(readfile(fedataraw, check.names = FALSE, row.names = 1))
expdata <- readfile(fedata)
samples <- filters_subset_df(selectss, annot, v = T)
mygenes <- rownames(expdata) # in case you wanna filter

dir.create(outdir)
source('/mnt/BioHome/ciro/scripts/functions/hvg.R')
sufix <- if(hvg_n != 800) hvg_n
fname <- paste0(outdir, '/_highly_variable_genes', sufix, '.pdf')
if(!file.exists(fname)) pdf(fname)
var_df <- getMostVariableGenes(
  counts = expdata_counts[mygenes, samples],
  normalise = TRUE,
  padjthr = 0.05,
  fitv = 0.5,
  top_n = hvg_n,
  take_sig = FALSE,
  plot = !is.file.finished(fname),
  verbose = TRUE
)
graphics.off()
keep_genes <- names(var_df$means)[var_df$top_n]
expdatass <- log2(expdata[keep_genes, samples] + 1)
summary(rowMeans(expdatass[keep_genes, samples]))
cnames <- unique(c(cnames, grep(paste0(showgenes, collapse = "|"), rownames(expdata), value = TRUE)))

# Clustering
df <- stats::dist(t(scale(expdatass))) # t(scale(expdatass)) same result
cluster_prefix = "hcluster"
summary(c(abs(df)))
str(df)
source("/home/ciro/scripts/handy_functions/devel/plots_correlation.R")
sub_grps <- nclust_optimal(mat = df, prefix = paste0(outdir, '/_', cluster_prefix, sufix))

annot <- annot[, grep(cluster_prefix, colnames(annot), invert = TRUE)]
for(i in names(sub_grps)){
  clusters <- sub_grps[[i]]
  k_clust = gsub("N", "", i)
  annot$hcluster = unname(clusters[rownames(annot)])
  annot$hcluster[!is.na(annot$hcluster)] <- paste0(cluster_prefix, "_", annot$hcluster[!is.na(annot$hcluster)])
  cnames <- unique(c(cnames, paste0(cluster_prefix, k_clust)))
  colnames(annot) <- gsub("hcluster$", paste0(cluster_prefix, k_clust), colnames(annot))
}
table(annot$Cell_type, annot$hcluster3)

pca_res <- prcomp(expdatass, center = TRUE)
pca_df <- data.frame(pca_res$rotation, stringsAsFactors = FALSE, check.names = FALSE)
annot$Sex_Disease <- paste0(annot$Sex, "_", annot$Disease)

fname <- paste0(outdir, '/_elbow.pdf')
if(!is.file.finished(fname)){
  pdf(fname)
  plot(1:length(pca_res$sdev), pca_res$sdev)
  dev.off()
}

fname <- paste0(outdir, '/_heatmap.pdf')
if(!file.exists(fname)){
  source('/home/ciro/scripts/handy_functions/devel/plots.R')
  source('/home/ciro/scripts/handy_functions/R/stats_summary_table.R')
  annot$Cell_type <- factor(annot$Cell_type, c("CD4-TRM", "CD4-TRM-like", "CD4-Teff", "CD4-Treg"))
  # tvar <- "/home/ciro/large/asthma_pjs/results/deseq2/airways_platex/summary/resting_celltypes_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"
  # degs <- readfile(tvar, stringsAsFactors = FALSE)
  # keep_genes <- sub("'", "", degs$gene_name[!is.na(degs$group)])
  cnames_i <- if(grepl("map2", fname)) colnames(annot) %in% cnames else c("Cell_type", "Disease")[1]
  pdf(fname, width = 8, height = 8, onefile = FALSE)
  custom_heatmap(
    object = if(grepl("map2", fname)) expdatass else expdata[keep_genes, samples],
    annoc = annot[, colnames(annot) %in% cnames_i, drop = FALSE],
    use_mean = if(grepl("map2", fname)) FALSE else cnames_i,
    orderby = c("Cell_type"), feature_order = "pca", couls = grcols, verbose = TRUE
  )
  graphics.off()
}

pt_size = 4
for(pc in pcs){
  for(px in pxs){
    set.seed(27)
    drdata <- if(grepl("tsne", outdir)){
      dim_axes = c("t-SNE 1", "t-SNE 2")
      Rtsne::Rtsne(dist(pca_df[, 1:pc]), perplexity = px, pca = FALSE, theta = 0, is_distance = TRUE, max_iter = 1000)$Y
    }else{
      dim_axes = c("UMAP 1", "UMAP 2")
      uwot::umap(t(expdatass), pca = pc, n_neighbors = px)
    }
    ddf <- data.frame(x.var = drdata[, 1], y.var = drdata[, 2])
    rownames(ddf) <- rownames(pca_df)
    # for(cname in grep("^Cell_type|hcluster", colnames(annot), value = TRUE)){
    for(cname in cnames){
      fname <- paste0(outdir, '/hv', length(keep_genes), 'genes_', paste0(pc, "PCs_", px), 'PX_', cname, '_size', pt_size, '.pdf')
      cat(fname, '\n'); #if(file.exists(fname)) next
      ddf$col <- if(cname %in% colnames(annot)) factormix(annot[samples, cname]) else "NONE"
      ddf$col <- if(cname %in% rownames(expdata)) log2(c(expdata[cname, rownames(ddf)]) + 1) else ddf$col
      # ddf$col <- if(cname %in% rownames(expdata)) factor(c(expdata[cname, rownames(ddf)]) > 10) else ddf$col
      scolor <- if(is.numeric(ddf$col)){
        scale_colour_gradientn(colours = couls)
      }else{ scale_colour_manual(values = v2cols(levels(ddf$col), grcols)) }

      gg.d <- ggplot(ddf, aes(x = x.var, y = y.var)) +
          geom_point(size = pt_size, shape = 19, alpha = 0.7, aes(colour = col)) +
          scolor + theme_classic() +
          labs(
            title = gsub("_|.pdf", " ", basename(fname)),
            x = dim_axes[1], y = dim_axes[2], colour = "Group"
          ) + theme_bw()+
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(size = 2),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.7),
            legend.position = "right", legend.text = element_text(size = 10)
          )

      pdf(fname, width = 8, height = 7)
      print(gg.d)
      dev.off()
    }
  }
}

# dir.create(paste0(outdir, "/archived"))
# list.files(outdir, pattern = "hv20414genes", full.names=T)
# system(paste0("mv ", outdir, "/hv", length(keep_genes), "gene* ", outdir, "/archived/"))
