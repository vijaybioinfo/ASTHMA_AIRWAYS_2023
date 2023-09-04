#!/usr/bin/R5

###############
# Extra tasks #
###############

# This script will create plots and tables for exploration of our airways cd4 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

## Destiny Folder ##
dirfig <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/a1_final_figures_cd4/qcs'
setwdc(dirfig)
clustnamebk = "none"

# Global variables
colsname <- "/home/ciro/asthma_pjs/info/airways_global_colours.csv"
grcols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1);

### QC per library - class ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p/clustering/zetInfo/clustCells19PCs_30Ks_0.06667JD.RData"
out_dir = paste0(dirfig, '/var_report_')

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
if(!dir.exists(out_dir)) dir.create(out_dir)
source('/mnt/BioHome/ciro/scripts/clustering/cluster_report.R')
void <- cluster_reports(
  object = mycells,
  confounders = c('orig.class'),
  resolutions = c('origlib'),
  reductions = list(tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2')),
  prefix = out_dir,
  gcols = grcols,
  normalise_bar = FALSE,
  v = TRUE
); graphics.off()

## Checking bar plot - it didn't need normalise = TRUE!
source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
dfplot <- mycells@meta.data
pp <- plot_pct(x = dfplot, groups = c("orig.class", "origlib"), normalise = TRUE, return_table = TRUE)
pp$table
pp <- plot_pct(x = dfplot, groups = c("orig.class", "origlib"), normalise = FALSE, return_table = TRUE)
pp$table

p2 <- pp$plot + coord_flip() + mytheme
pdf("var_report_barplot.pdf", width = 5, height = 15)
print(p2)
dev.off()

# Donots are good!! What a relief... for now.
p <- plot_pct(x = dfplot, groups = c("orig.class", "origlib"), type = "donut", return_table = TRUE)
p$table

### Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = c(
  "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/clustering/zetInfo/clustCells19PCs_30Ks_0.06667JD.RData",
  "/home/ciro/large/asthma_pjs/results/clust_seurat/stim_cd4_15p_sng_nomacro/clustering/zetInfo/clustCells18PCs_30Ks_0.06667JD.RData"
)
low.filt.cells = c(200, -Inf, -Inf)
high.filt.cells = c(6000, 30000, 0.15)
subs.names = c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
names(low.filt.cells) <- subs.names
names(high.filt.cells) <- subs.names

## Reading
lmycells <- lapply(clustname, theObjectSavedIn)
mycellscom <- if(length(lmycells) > 1){
  tvar <- lmycells[[1]]
  for(i in 2:length(lmycells)){
    tvar <- merge(tvar, lmycells[[i]])
  }
  tvar
}else{ lmycells }
rm(lmycells)
tvar <- reshape2::melt(table(mycellscom$orig.stim, mycellscom$origlib))
tvar <- tvar[tvar[, 3] > 0, ]
tvar <- tvar[order(tvar[, 1]), ]
mycellscom$Library <- factor(mycellscom$origlib, levels = unique(as.character(tvar[, 2])))
table(mycellscom$Library)

## Operations
mycols <- rep("#2e82ff", length(table(mycellscom$Library)))
names(mycols) <- names(table(mycellscom$Library))
for(qcvar in subs.names[1]){
  p <- VlnPlot(
    object = mycellscom,
    features = qcvar,
    group.by = "Library",
    cols = mycols
  ) + NoLegend()
  fname <- paste0("sf2a_qc_", qcvar)
  pdf(paste0(fname, ".pdf"), width = 12)
  print(rm_layer(p, "point"))
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 12)
  print(shut_it(p, "point|txt"))
  dev.off()
}

### Donor-wise summary ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
pref = "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p/zetInfo/sources/metadata.RData"
pref = "/home/ciro/large/asthma_pjs/results/demuxlet/sever_asthma/resting_cd4_15p_16PCsmetadata.Rdata"
postf = "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/clustering/zetInfo/metadata_22PCs_30Ks_0.06667JD.RData"

## Reading
annot_pre <- theObjectSavedIn(pref)
annot_post <- theObjectSavedIn(postf)

tvar <- table(annot_pre[, c("origlib", "orig.class")], useNA = 'always')
tvar <- tvar[rowSums(tvar) > 0, ]
dbrate <- if("Doublet" %in% colnames(tvar)) round(tvar[, 1] / table(annot_pre$origlib)[rownames(tvar)] * 100, 2) else 0
dprops <- cbind(
  tvar,
  dbrate,
  table(annot_pre$origlib)[rownames(tvar)]
)
colnames(dprops) <- c(colnames(tvar)[-ncol(tvar)], "GexOnly", "Doublet_Percent", "Total")
dprops
write.csv(dprops, file = gsub("PCs.*", "PCs_doublet_rate.csv", postf))
table(annot_pre[, c("origlib", "in_gex")], useNA = 'always')

annot_pre$orig.master <- annot_pre$orig.tmp <- annot_pre[, as_donor]
annot_pre$orig.tmp[is.na(annot_pre$orig.tmp)] <- "void"
annot_pre$orig.master <- paste0(annot_pre$orig.master, "-", annot_pre$orig.virus2)
annot_pre$orig.master[annot_pre$orig.tmp == "void"] <- annot_pre[annot_pre$orig.tmp == "void", ]$orig.virus2

annot_post$orig.master <- annot_post$orig.tmp <- annot_post[, as_donor]
annot_post$orig.tmp[is.na(annot_post$orig.tmp)] <- "void"
annot_post$orig.master <- paste0(annot_post$orig.master, "-", annot_post$orig.virus2)
annot_post$orig.master[annot_post$orig.tmp == "void"] <- annot_post[annot_post$orig.tmp == "void", ]$orig.virus2
table(annot_post[, as_donor], annot_post$orig.hospital, useNA = 'always')

table(annot_pre$orig.virus2, useNA = 'always')
table(annot_post$orig.virus2, useNA = 'always')
tvar <- reshape2::melt(table(annot_post$orig.master, annot_post$origlib)); tvar <- tvar[tvar$value > 0, ]
tvar <- table(annot_pre$orig.master)
tvar2 <- table(annot_post$orig.master)
ddf <- remove.factors(data.frame(
  tvar[names(tvar2)],
  tvar2
)[, c(1, 2, 4)])
colnames(ddf) <- c("Donor", "Pre", "Post")
summddf <- ddf
table(annot_post[, c('orig.HT_ID', 'orig.virus2')])
extra_cols <- extra_cols[extra_cols %in% colnames(annot_post)]
if(length(extra_cols) > 0){
  tvar <- data.frame(sapply(extra_cols, function(x){
    sapply(make_list(x = annot_post, colname = "orig.master", col_objects = x),
      function(x) paste0(unique(x), collapse = "|") )
  }), stringsAsFactors = FALSE)
  summddf <- cbind(summddf, tvar[ddf[, 1], , drop = FALSE])
  summddf <- summddf[, c("Donor", extra_cols, "Pre", "Post")]
}
str(summddf)
pcs <- gsub(".*_([0-9]{1,}PCs_).*", "\\1", postf)
for(nres in grep("RNA_snn_res", colnames(annot_post), value = TRUE)){
  cat(nres, "\n")
  annot_post$tmp <- factormix(annot_post[, nres])
  tmp <- as.data.frame.matrix(table(annot_post[, c("orig.master", 'tmp')]))
  tvar <- cbind(summddf, tmp[summddf[, 1], , drop = FALSE])
  data.table::setorderv(x = tvar, cols = rev(extra_cols))
  sums <- lapply(extra_cols, function(caty){
    y <- plyr::ddply(.data = tvar[tvar[, caty] != "NA", ], .variable = caty, function( slice ){
        sapply(c("Pre", "Post", colnames(tmp)), function(x) sum(slice[[x]]) )
    })
    y[!grepl("^HT", y[, 1]), ]
  })
  tvar2 <- data.table::rbindlist(c(list(tvar), sums), fill = TRUE)
  tvar2[is.na(tvar2)] <- ""
  fname <- paste0(dirnamen(postf, 2), "/", pcs, nres, "/proportions_a1_filtering_process.csv")
  write.csv(tvar2, file = fname, row.names = FALSE)
}
