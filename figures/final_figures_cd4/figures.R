#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2020-12-09
# ---

# This script will create plots for any figure from our cd4 data
# Each figure's code will then probably be moved to its correct file

# include = "sc_cd4resting"
source("/home/ciro/asthma_airways/scripts/final_figures_cd4/global.R")

### Dimentional reduction plot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = paste0("f3a_", names(redu))
myobject = "sc_cd4resting"
myobject = "sc_cd4stim"
scells = colnames(eval(parse(text = myobject)))
mycols = eval(parse(text = paste0(myobject, "_ident")))$colours
group_by = eval(parse(text = paste0(myobject, "_clust")))

dir.create("umaps")
group_by = "orig.asthma"
mycols = c(MA = "#0000FF", SA = "#FB0207")
result_id = paste0("umaps/", myobject, "_disease_all")
# result_id = paste0("umaps/", myobject, "_disease_female")
# scells <- filters_subset_df(c("orig.Sex", "Female"), eval(parse(text = myobject))@meta.data, v = TRUE)
# result_id = paste0("umaps/", myobject, "_disease_male")
# scells <- filters_subset_df(c("orig.Sex", "Male"), eval(parse(text = myobject))@meta.data, v = TRUE)
# scells <- sample_even(annot = eval(parse(text = myobject))@meta.data[scells, ], cname = group_by, v = TRUE)

p <- DimPlot(
  object = eval(parse(text = myobject))[, scells],
  cols = mycols,
  reduction = names(redu),
  group.by = group_by
) + labs(x = redu[[1]][1], y = redu[[1]][2])

pdf(paste0(result_id, ".pdf"))
print(p)
graphics.off()
pdf(paste0(result_id, "_blank.pdf"))
print(plot_blank(p))
graphics.off()

### Disease dimentional reduction plot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "disease_split_"
cname = "orig.asthma"

scells <- sample_even(annot = sc_cd4resting@meta.data, cname = cname, v = TRUE)
p <- DimPlot(
  object = sc_cd4resting[, scells],
  # cells = scells,
  cols = v2cols(sc_cd4resting@meta.data[, sc_cd4resting_clust], sc_cd4resting_ident$colours),
  reduction = names(redu),
  group.by = sc_cd4resting_clust,
  label = TRUE, #pt.size = 0.8
  split.by = cname,
) + labs(x = redu[[1]][1], y = redu[[1]][2]) + theme(legend.position = "none")

fname <- paste0(result_id, names(redu))
pdf(paste0(fname, ".pdf"), width = 12)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 12)
print(plot_blank(p))
graphics.off()

### Molecules: dot-plots/violins ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create("dotplots")
genes = c(
  "IL21", "IL6", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "GZMB", "TNF", "TNFSF14",
  "TNFSF12", "TNFSF10", "AREG", "IL2", "IL3", "CSF2", "CCL20", "IFNG", "IL4",
  "IL5", "IL13", "IL17A", "IL17F", "IL22", "IL10", "TGFB1"
)

expr_thresh = 0
result_id = "dotplots/sc_cd4merged"
if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
sc_cd4merged$Stim = ifelse(sc_cd4merged$orig.stim == "ST", "Stim", "Resting")
myobject = "sc_cd4merged"
columns_i = c("Stim", "orig.asthma")
eval(parse(text = paste0(
  myobject, "@meta.data$Group = ident_combine(", myobject, "@meta.data, columns_i)"
)))
table(eval(parse(text = myobject))@meta.data$Group)
mygenes = show_found(genes, rownames(eval(parse(text = myobject))), v = TRUE)
scells = sample_even(eval(parse(text = myobject))@meta.data, "Group", -3000, v = TRUE)
edata_i = expm1(eval(parse(text = myobject))@assays$RNA@data[mygenes, scells])
mdata_i = eval(parse(text = myobject))@meta.data[scells, ]

expr_thresh = 10
bulk_metadata$Stim1 = ifelse(bulk_metadata$Stim == "Unstim", "Resting", "Stim")
columns_i = c("Stim1", "Disease")
bulk_metadata$Group = ident_combine(bulk_metadata, columns_i)
table(bulk_metadata$Cell_type)
sselect = list(c("Cell_type", "CD4-Teff"))
sselect = list(c("Cell_type", "CD4-TRM-like"))
sselect = list(c("Cell_type", "CD4-TRM"))
result_id = paste0("dotplots/bulk_", sselect[[1]][-1])
ssamples <- filters_subset_df(sselect, bulk_metadata, v = T)
edata_i = bulk_tpm[features_parse_ensembl(rownames(bulk_tpm)) %in% genes, ssamples]
rownames(edata_i) = features_parse_ensembl(rownames(edata_i))
mygenes = show_found(genes, rownames(edata_i), v = TRUE)
mdata_i = bulk_metadata[ssamples, ]

source("/home/ciro/scripts/handy_functions/devel/plots_dotplot.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
for(scale_mean_q in c(FALSE, TRUE)){
  p <- dot_plot(
    edata = edata_i,
    mdata = mdata_i,
    columns = "Group",
    scale_mean = scale_mean_q,
    features = mygenes,
    cols_limits = if(scale_mean_q) c(0, 1.5),
    values_norm = if(scale_mean_q) c else function(x) log2(x+1),
    clust_rows = !TRUE, cols = c('#fffef0', '#ff0000'), values_exprthr = expr_thresh
  )
  fname <- paste0(result_id, ifelse(scale_mean_q, "_mean_scaled", ""))
  cat(fname, "\n")
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(p)
  graphics.off()
  pdf(paste0(fname, "_blank.pdf"), width = 7.5, height = 10)
  print(plot_blank(p))
  graphics.off()
}

pdata = cbind(mdata_i, t(as.matrix(eval(parse(text = myobject))@assays$RNA@data[mygenes, rownames(mdata_i)])))
# "violins_sc_cd4merged/dots_"
for(result_id in c("violins_sc_cd4merged/all_", "violins_sc_cd4merged/chilli_")){
  if(!dir.exists(dirname(result_id))) dir.create(dirname(result_id))
  cat(result_id, "\n")
  for(g in mygenes){
    cat(" -", g, "\n")
    p <- violin(
      dat = pdata,
      xax = "Group",
      yax = g,
      dots = grepl("dots", result_id),
      colour_by = "pct", couls = 3, vsetting = list(trim = FALSE),
      chilli = grepl("chilli", result_id)
    )
    p <- plot_rm_layer(p, "Box")

    fname <- paste0(result_id, g)
    pdf(paste0(fname, ".pdf"))
    print(p)
    graphics.off()
    pdf(paste0(fname, "_blank.pdf"))
    print(plot_blank(p))
    graphics.off()
  }
}

### Dot-plot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "dotplot" #paste0("dotplot", format(Sys.time(), '_%Y_%m_%d'))
genes = c(
  "ITGAE", "ITGA1", "AMICA1", "HOPX", "CLU", "ALOX5AP", "ZNF683", "GZMA", "GZMB",
  "GZMH", "MX2", "LY6E", "IFITM1", "CXCL13", "PDCD1", "CTLA4", "MAF", "TIGIT",
  "FOXP3", "IL2RA", "IKZF2", "ICAM2", "SELL", "CCR7", "TCF7", "S1PR1", "HSPA1A",
  "IFNG", "FOS", "TUBB", "TOP2A", "MKI67"
)

mygenes <- show_found(genes, rownames(sc_cd4resting), v = TRUE)

p <- DotPlot(
  object = sc_cd4resting,
  features = mygenes,
  group.by = "cluster",
  # cols = c('#fff4ba', '#ff0000'), col.min = -1.5, col.max = 1.5,
  cols = c('#fffef0', '#ff0000'), col.min = 0, col.max = 1.5,
  dot.min = 0.1
) + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)

fname <- paste0(result_id, "")
pdf(paste0(fname, ".pdf"), width = 10, height = 10)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 7.5, height = 10)
print(plot_blank(p))
graphics.off()

source("/home/ciro/scripts/handy_functions/devel/plots.R")
p <- dot_plot(
  edata = sc_cd4resting@assays$RNA@data,
  mdata = sc_cd4resting@meta.data,
  columns = "cluster",
  features = mygenes,
  cols_limits = c(0, 1.5),
  values_trans = expm1, clust_rows = TRUE, clust_cols = TRUE, cols = c('#fffef0', '#ff0000')
)
fname <- paste0(result_id, "_hclust")
pdf(paste0(fname, ".pdf"), width = 10, height = 10)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 7.5, height = 10)
print(plot_blank(p))
graphics.off()

### Markers heatmap ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntop = 200
result_id = "heatmap_markers_cd4resting_shared_mean"
marknames = paste0(
  dirnamen(global_objects_f['sc_cd4resting'], 3), "/markers/20PCs_",
  sc_cd4resting_clust,
  "_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
)
ident_info_i = ident_order(sc_cd4resting_ident)
myobject = "sc_cd4resting"
master_column = "cluster"

result_id = "heatmap_markers_cd4resting0n1_shared_mean"
marknames = paste0(
  dirnamen(sc_cd4resting_0n1_f, 3), "/markers/20PCs_",
  sc_cd4resting_0n1_clust,
  "_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
)
ident_info_i = list()
myobject = "sc_cd4resting_0n1"
master_column = sc_cd4resting_0n1_clust

result_id = "heatmap_markers_cd4stim0.4_shared_mean"; sc_cd4stim_clust = "RNA_snn_res.0.4"
result_id = "heatmap_markers_cd4stim0.6_shared_mean"; sc_cd4stim_clust = "RNA_snn_res.0.6"
marknames = paste0(
  dirnamen(sc_cd4stim_f, 3), "/markers/17PCs_",
  sc_cd4stim_clust,
  "_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
)
ident_info_i = list()
myobject = "sc_cd4stim"
master_column = sc_cd4stim_clust

mygenes <- readfile(marknames, stringsAsFactors = FALSE, row.names = 1)

if(!grepl("uniq", result_id)){
  mygenesl <- make_list(mygenes, colname = "cluster", col_objects = "gene")
  str(mygenesl)
  uddf <- as.data.frame.matrix(table(mygenes[, c("gene", "cluster")]))
  if(length(ident_info_i) > 0){
    uddf <- uddf[, names(ident_info_i[[1]])]
    colnames(uddf) <- paste0("C", colnames(uddf), " ", ident_info_i[[1]][colnames(uddf)])
    str(uddf)
  }
  pdf(paste0(result_id, "_upset.pdf"), onefile = FALSE)
  print(UpSetR::upset(data = uddf, sets = colnames(uddf)))
  dev.off()
}

mygenes$gene <- gsub("'", "", mygenes$gene_name)
tvar <- mygenes$Dpct > .0; table(tvar)
mygenes <- mygenes[tvar, ]
tvar <- mygenes$avg_logFC > .25; table(tvar)
mygenes <- mygenes[tvar, ]
if(grepl("uniq", result_id)){
  tvar <- !grepl("&", mygenes$sCluster); table(tvar)
  mygenes <- mygenes[tvar, ]
}
tvar <- if(length(ident_info_i) > 0) names(ident_info_i[[1]]) else gtools::mixedsort(unique(mygenes$cluster))
mygenes$cluster <- factor(mygenes$cluster, tvar)
mygenes <- mygenes[order(mygenes$cluster), ]
# mygenes$nclust <- stringr::str_count(mygenes$sCluster, "&")
topgenes <- get_top_n(x = mygenes, n = ntop)
fname <- paste0(result_id, ifelse(nrow(topgenes) != nrow(mygenes), paste0("_top", ntop), ""))
# topgenes <- get_top_n(x = topgenes, n = ntop, orderby = 'nclust'); fname <- paste0(fname, "_nclustOrder")
topgenes <- topgenes[!duplicated(topgenes$gene), ]
genes <- gsub("'", "", topgenes$gene_name)
genes <- show_found(genes, rownames(eval(parse(text = myobject))), v = TRUE)
genesl <- make_list(x = topgenes, colname = "cluster", col_objects = "gene_name")

annor = reshape2::melt(lapply(genesl, gsub, pattern = "'", replacement = ""))
annor <- data.frame(
  RowGroup = factor(annor$L1, names(genesl)),
  row.names = as.character(annor$value))
write.csv(annor, file = paste0(fname, ".csv"))

source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
# pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
png(paste0(fname, ".png"), width = 1500, height = 1701, res = 250)
custom_heatmap(
  object = eval(parse(text = myobject)),
  rnames = genes,
  orderby = master_column,
  use_mean = master_column,
  sample_it = c(cname = master_column, maxln = "-1000"),
  scale_row = TRUE,
  categorical_col = c(master_column),
  feature_order = TRUE,
  couls = ident_info_i$colours,
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annor, annotation_names_row = FALSE, annotation_names_col = FALSE
)
graphics.off()

### X.X CITE-Seq normalise ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "cite_seq_"
csname = "/mnt/BioAdHoc/Groups/vd-vijay/ndu/Old_Backup/CITE_SEQ/results/008_10_025_4U_1010_CITE_RA_RPI3_S3_filter.csv"

# centered log-ratio (CLR) normalization - for each feature!
clr_function <- function(x) {
  apply(X = x, MARGIN = ifelse(nrow(x) < ncol(x), 1, 2), FUN = function(x){
    return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
  })
}

cite_dat <- readfile(csname, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

## Operations
mylib_n <- sub(".*\\-(.*)", "\\1", rownames(sc_cd4resting@meta.data[grep("003_10x_025_4U", sc_cd4resting@meta.data$origlib)[1], ]))
if(!grepl(paste0("-", mylib_n, "$"), rownames(cite_dat)[2])) rownames(cite_dat) <- paste0(rownames(cite_dat), "-", mylib_n)
cite_dat_clr <- clr_function(cite_dat)
colnames(cite_dat_clr) <- gsub("_$", "", gsub("_{1,}", "_", gsub(" |\\(|\\)", "_", colnames(cite_dat_clr))))
write.csv(cite_dat_clr, file = paste0(result_id, "normalized.csv"))

### CITE-seq data presence ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "cite_seq_"
csname = paste0(result_id, "normalized.csv")
dims = c(redu[[1]], sc_cd4resting_clust)
mycols = c("#fff5bf", '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')

cite_dat <- readfile(csname, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
sc_cd4resting@meta.data$CITEseq <- ifelse(rownames(sc_cd4resting@meta.data) %in% rownames(cite_dat), "CITEseq", "NonCITEseq")

tvar <- reshape2::melt(table(sc_cd4resting@meta.data[, c(sc_cd4resting_clust, "orig.donor", "CITEseq")]))
head(tvar)
write.csv(tvar, file = paste0(result_id, "cluster_donor.csv"), row.names = FALSE)

p <- DimPlot(
  object = sc_cd4resting,
  cols = c("CITEseq" = "black", "NonCITEseq" = "#f5f5f5"),
  reduction = names(redu),
  group.by = "CITEseq",
  pt.size = 0.5
) + labs(x = redu[[1]][1], y = redu[[1]][2])

fname <- paste0(result_id, names(redu))
pdf(paste0(fname, ".pdf"))
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"))
print(plot_blank(p))
graphics.off()

cellsinter <- show_found(rownames(sc_cd4resting@meta.data), rownames(cite_dat), v = TRUE)
cite_df <- data.frame(mat_names(rownames(sc_cd4resting@meta.data), c(dims, colnames(cite_dat))), check.names = FALSE)
cite_df[, dims] <- FetchData(sc_cd4resting, vars = dims)#, cells = cellsinter)
cite_df[cellsinter, colnames(cite_dat)] <- cite_dat[cellsinter, ] # cbind(cite_df, cite_dat_clr[cellsinter, ])
head(cite_df[complete.cases(cite_df), ])

# for(g in colnames(cite_df)[-c(1:3)]){
for(g in c("CD103_Integrin_E", "CD69")){
  cat(g, '\n')
  fname <- paste0(result_id, g)
  p <- ggplot(cite_df, aes_string(x = dims[1], y = dims[2], color = g)) +
    geom_point(size = 0.5) +
    scale_colour_gradientn(colours = mycols, na.value="#f5f5f5")
  pdf(paste0(fname, ".pdf"))
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"))
  print(plot_blank(p))
  dev.off()
}

### Per donor/cluster proportion tables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "cellcount_"
dcolumn = 'orig.donor'
subdiv = "orig.asthma"
corder = as.character(0:8)
roder = NULL

tvar1 = list(
  list(
    add_data = "cite_seq_normalized.csv",
    genes = c("ITGAE"),
    prefix = "cd103_resting",
    sselect = list(c('tag_CD103_Integrin_E', 'CD103_Integrin_E+'))
  ),
  list(
    add_data = "cite_seq_normalized.csv",
    genes = c("ITGAE"),
    prefix = "cd69_resting",
    sselect = list(c('tag_CD69', 'CD69+'))
  )
)
tvar2 <- lapply(c("GZMB", "CXCR4", "HOPX", "FKBP5", "IFNG"), function(x){
  list(
    add_data = NULL,
    genes = c(x),
    prefix = paste0(x, "_resting"),
    sselect = features_get_tag(x)
  )
})
fig_parms_l <- c(tvar1, tvar2)

for(fig_parms in fig_parms_l){
  mygenes <- show_found(fig_parms$genes, rownames(sc_cd4resting), v = TRUE)
  markersdf <- features_add_tag(
    lgenes = mygenes,
    annot = sc_cd4resting@meta.data,
    mat = sc_cd4resting@assays$RNA@data,
    verbose = TRUE
  )

  sselect_i <- fig_parms$sselect
  if(!is.null(fig_parms$add_data)){
    tvar <- readfile(fig_parms$add_data, row.names = 1, stringsAsFactors = FALSE)
    icells <- intersect(rownames(tvar), rownames(sc_cd4resting@meta.data))
    tvar <- tvar[icells, ]
    if(grepl("^tag_", fig_parms$sselect[[1]][1])){
      tvar <- features_add_tag(
        lgenes = gsub("tag_", "", fig_parms$sselect[[1]][1]),
        annot = sc_cd4resting@meta.data[icells, ],
        mat = t(tvar),
        verbose = TRUE
      )
      sselect <- lapply(fig_parms$sselect, casefold, upper = TRUE)
      sselect_i[[1]][1] <- gsub("TAG_", "tag_", sselect_i[[1]][1])
    }
    markersdf <- joindf(markersdf, tvar)
  }
  tailmat(markersdf)
  compiled_tab <- joindf(sc_cd4resting@meta.data, markersdf)
  tailmat(compiled_tab)

  scells <- filters_subset_df(sselect_i, compiled_tab, op = "or", v = TRUE)
  annot <- remove.factors(compiled_tab[scells, ])
  table(annot[, c(dcolumn, subdiv)])
  tvar <- make_list(annot, colname = dcolumn, col_objects = subdiv)
  tvar <- reshape2::melt(lapply(tvar, function(x) paste0(unique(x), collapse = "-") ))
  dgroup <- as.character(tvar[, 1])
  names(dgroup) <- as.character(tvar[, 2])

  ddfprops <- table(annot[, c(dcolumn, sc_cd4resting_clust)])
  # Proportions heatmap
  matpct <- as.matrix(as.data.frame.matrix(prop.table(ddfprops, margin = 1)))
  matpct <- round(matpct, 2)
  matpct <- if(!is.null(roder)){
    sufix = "_supervised"
    matpct[roder, corder]
  }else{
    sufix = "_unsupervised"
    matpct[, corder]
  }
  topz <- c(0, .7)
  matpct[matpct > topz[2]] <- topz[2]; matpct[matpct < topz[1]] <- topz[1];
  palettebreaks <- seq(from = 5, to = 70, by = 5)
  annoc <- data.frame(
    Group = dgroup, row.names = names(dgroup)
  )
  annoc <- annoc[rownames(matpct), sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
  anncolist <- lapply(annoc, function(x) v2cols(x, NULL, v = TRUE) )
  mypalette <- colorRampPalette(colors = c("white", "red"), space = 'Lab')
  fname <- paste0(result_id, "heatmap_proportions_", fig_parms$prefix, sufix)
  pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
  x <- NMF::aheatmap(
    x = matpct, annRow = annoc, annColors = anncolist,
    scale = 'none', Rowv = ifelse(is.character(roder), NA, TRUE), Colv = NA,
    col = mypalette(length(palettebreaks) - 1)
  )
  graphics.off()
  write.csv(ddfprops[rev(x$rowInd), ], file = paste0(fname, ".csv"))
}

### Per donor stats bulk/single-cell ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "expression_tab_multiseq"
genes = c(
  "GZMB", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "CCL20", "IL21", "IL17A",
  "IL17F", "TNF", "TNFSF14", "TNFSF12", "TNFSF10", "AREG", "IL2", "IL3", "CSF2",
  "IFNG", "IL4", "IL5", "IL13", "IL22", "IL10", "TGFB1", "IL6"
)
dconfigs = list(
  list(
    edata = "bulk_tpm", mdata = "bulk_metadata",
    filter = c("Cell_type", "CD4-TRM"),
    columns = c("Study_ID", "Cell_type")
  )
)
tmp <- features_parse_ensembl(rownames(bulk_tpm))
tvar <- show_found(genes, tmp, v = TRUE); tvar <- tmp %in% genes
tmp <- filters_subset_df(list(c("Cell_type", "CD4-TRM"), c("Stim", "Stim")), bulk_metadata, v = TRUE)
bulk_tpm_i = bulk_tpm[tvar, tmp]
tvar <- gsub("^([0-9]{1,}_NIH[A-Z]{1,}_[0-9]{1,}).*", "\\1", colnames(bulk_tpm_i))
colnames(bulk_tpm_i)[tvar %in% tvar[duplicated(tvar)]] # duplicated if leading num is ignored
colnames(bulk_tpm_i) <- tvar
rownames(bulk_tpm_i) <- features_parse_ensembl(rownames(bulk_tpm_i))
summarydf0 <- reshape2::melt(bulk_tpm_i)
colnames(summarydf0) <- c("Gene", "bulk_donor", "bulk_stim_tpm")
summarydf0$Donor <- gsub("^[0-9]{1,}_", "", summarydf0$bulk_donor)
summarydf0 <- summarydf0[, c("Gene", "Donor", "bulk_donor", "bulk_stim_tpm")]
head(summarydf0)

tvar <- show_found(genes,rownames(sc_cd4stim), v = TRUE)
sc_cd4stim_i <- stats_summary_table(
  mat = expm1(sc_cd4stim@assays$RNA@data)[genes, ] * 100,
  group = make_list(sc_cd4stim@meta.data, colname = "orig.donor", grouping = TRUE),
  moments = c("p", "pmn"),
  verbose = TRUE
)
summarydf1 <- reshape2::melt(cbind(Gene = rownames(sc_cd4stim_i), sc_cd4stim_i))
summarydf1$type <- gsub(".*_", "", summarydf1$variable)
summarydf1$variable <- gsub("_percent.*|_mean.*", "", summarydf1$variable)
summarydf1 <- reshape2::dcast(summarydf1, Gene + variable ~ type, value.var = "value")
colnames(summarydf1) <- c("Gene", "Donor", "sc_stim_meanOfPositive", "sc_stim_percentage")
head(summarydf1)

# Merging
bulk_ids <- paste0(summarydf0$Gene, summarydf0$Donor)
sc_ids <- paste0(summarydf1$Gene, summarydf1$Donor)
if(1){
  cat("Duplicates in bulk?", sum(duplicated(bulk_ids)), "\n") # technical replicates?
  bulk_ids[duplicated(bulk_ids)] <- paste0(summarydf0$Gene, summarydf0$bulk_donor)[duplicated(bulk_ids)]
  cat("After using the leading num of the tech. dup?", sum(duplicated(bulk_ids)), "\n")
  cat("Duplicates in single-cell?", sum(duplicated(sc_ids)), "\n")
  cat("All bulk donors in single-cell?", all(summarydf0$Donor %in% summarydf1$Donor), "\n")
  cat("All single-cell donors in bulk?", all(summarydf1$Donor %in% summarydf0$Donor), "\n")
}
rownames(summarydf1) <- sc_ids
rownames(summarydf0) <- bulk_ids
summarydf <- joindf(summarydf1, summarydf0)
summarydf <- summarydf[, unique(c(colnames(summarydf0), colnames(summarydf1)))]
summarydf <- data.table::rbindlist(list(summarydf, summarydf0[!rownames(summarydf0) %in% rownames(summarydf), ]), fill = TRUE)
head(summarydf)
summarydf[summarydf$Donor %in% "NIHMA_249", ]
write.csv(summarydf, file = paste0(result_id, ".csv"), row.names = FALSE)

### Per donor/cluster expression ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
result_id = "expression_tab_"
identity_column <- c("orig.donor", sc_cd4resting_clust)
genes = c("GZMB")
identity_column <- c("orig.donor")
genes = c(
  "IL21", "IL6", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "GZMB", "TNF",
  "TNFSF14", "TNFSF12", "TNFSF10", "AREG", "IL2", "IL3", "CSF2", "CCL20",
  "IFNG", "IL4", "IL5", "IL13", "IL17A", "IL17F", "IL22", "IL10", "TGFB1"
)

mygenes <- show_found(genes, rownames(sc_cd4resting), v = TRUE)
sc_cd4resting$Identity <- ident_combine(sc_cd4resting@meta.data, identity_column, sep = "~")
summarydf <- t(stats_summary_table(
  mat = expm1(sc_cd4resting@assays$RNA@data) * 100,
  group = make_list(sc_cd4resting@meta.data, colname = "Identity", grouping = TRUE),
  rnames = mygenes,
  moments = c("mn", "p"),
  sep_str = "~",
  verbose = TRUE
))
tvar <- t(sapply(rownames(summarydf), function(x) unlist(strsplit(x, "~")) ))
master_tab <- data.frame(tvar)
master_tab <- cbind(master_tab, round(summarydf, 3))
fname <- paste0(result_id, ifelse(length(mygenes) > 3, length(mygenes), paste0(mygenes, collapse = "-")))
head(master_tab)
write.csv(master_tab, file = paste0(fname, ".csv"), row.names = FALSE)

### Gene set lists from bulk specificity ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "signatures_bulk_"
resf = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/summary/resting_celltypes_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"

res <- readfile(resf, row.names = 1, stringsAsFactors = FALSE)
t(res[grep("STAT3|CALM3", rownames(res)), grep("group|_mean", colnames(res))])
signature_bulk_list <- make_list(x = res, colname = "group", col_objects = "gene_name")
signature_bulk_list <- lapply(signature_bulk_list, function(x) features_parse_ensembl(sub("'", "", x)) )
str(signature_bulk_list)
signature_bulk_df <- vlist2df(signature_bulk_list[-length(signature_bulk_list)])
signature_bulk_df[is.na(signature_bulk_df)] <- ""
head(signature_bulk_df)
write.csv(signature_bulk_df, file = paste0(result_id, basename(dirname(resf)), ".csv"), row.names = FALSE)

### Gene set lists from bulk ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "signatures_bulk_"
degfilt = c(pv = 0.05, fc = 1)
fnames = list(
  list(
    name = "resting_TRMDPvsTEFF",
    file = "/home/ciro/large/asthma_airways/results/deseq2/airways/comprs/resting_trmness/CD4-TRMvsCD4-Teff/results_CD4-TRMvsCD4-Teff_deseq2.csv"
  ),
  list(
    name = "resting_TRMavsTEFF",
    file = "/home/ciro/large/asthma_airways/results/deseq2/airways/comprs/resting_cd4_merge/TRMavsnonTRM/results_TRMavsnonTRM_deseq2.csv"
  )
)

ldegs_final <- lapply(
  X = fnames,
  FUN = function(fname){
  res <- read.csv(fname$file, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  dfname <- paste0(result_id, fname$name, "_padj", degfilt["pv"], "_fc", degfilt["fc"], ".csv")
  grps <- unlist(strsplit(basename(dirname(fname$file)), "vs"))
  ldegs <- list(
    getDEGenes(res, pv = degfilt['pv'], fc = degfilt['fc'], upreg = TRUE, further = paste0(grps[2], "_meanTPM>5"), v = TRUE),
    getDEGenes(res, pv = degfilt['pv'], fc = degfilt['fc'], upreg = FALSE, further = paste0(grps[1], "_meanTPM>5"), v = TRUE)
  ); names(ldegs) <- grps
  ldegs10 <- list(
    getDEGenes(res, pv = degfilt['pv'], fc = degfilt['fc'], upreg = TRUE, further = paste0(grps[2], "_meanTPM>10"), v = TRUE),
    getDEGenes(res, pv = degfilt['pv'], fc = degfilt['fc'], upreg = FALSE, further = paste0(grps[1], "_meanTPM>10"), v = TRUE)
  ); names(ldegs10) <- paste0(gsub("\\-", "", grps), "_", fname$name)

  ddf <- data.frame(log2(res[, paste0(grps, "_meanTPM")] + 1))
  ddf2 <- data.frame(log2(res[unlist(ldegs), paste0(grps, "_meanTPM")] + 1))
  ddf3 <- data.frame(log2(res[unlist(ldegs10), paste0(grps, "_meanTPM")] + 1))
  p <- ggplot(ddf, aes_string(x = "CD4.TRM_meanTPM", y = "CD4.Teff_meanTPM")) + geom_point(color = "#BEBEBE") +
    geom_hline(yintercept = log2(11)) + geom_vline(xintercept = log2(11)) +
    geom_hline(yintercept = log2(6), linetype = "dashed") + geom_vline(xintercept = log2(6), linetype = "dashed") +
    geom_point(data = ddf2, color = "red") + geom_point(data = ddf3, color = "blue", shape = 3)
    print(dfname)
  # pdf(sub("csv$", "pdf", dfname));
  # print(p);
  # dev.off()

  return(ldegs10)
})

ldegs_final <- unlist(ldegs_final, recursive = FALSE)
ldegs_final <- lapply(ldegs_final, function(x) features_parse_ensembl(x) )
sapply(ldegs_final, length)
degs_df <- vlist2df(ldegs_final)
degs_df[is.na(degs_df)] <- ""
headtail(degs_df)
dfname <- paste0(result_id, length(ldegs_final), "sets_padj", degfilt["pv"], "_fc", degfilt["fc"], ".csv")
write.csv(degs_df, file = dfname, quote = FALSE, row.names = FALSE)
system(paste("head", dfname))

# Removing genes from other signatures
load('/home/ciro/scripts/handy_functions/data/signatures_vijaylab.rdata')
str(signatures_vijaylab)
genes2remove <- signatures_vijaylab[grep("treg|tox", names(signatures_vijaylab))]; str(genes2remove)
genes2remove <- unique(unlist(genes2remove))
lapply(X = ldegs_final, FUN = function(x) x[x %in% genes2remove] )
ldegs_final_clean <- lapply(X = ldegs_final, FUN = function(x) x[!x %in% genes2remove] )
# ldegs_final_clean <- ldegs_final
# tvar <- ldegs_final_clean$CD4TRM_resting_TRMDPvsTEFF
# ldegs_final_clean$CD4TRM_resting_TRMDPvsTEFF <- tvar[!tvar %in% genes2remove]
str(ldegs_final); str(ldegs_final_clean)
degs_df_clean <- vlist2df(ldegs_final_clean)
head(degs_df_clean)
colnames(degs_df_clean) <- paste0(colnames(degs_df_clean), "clean")
write.csv(degs_df_clean, file = sub(".csv$", "_curated.csv", dfname), quote = FALSE, row.names = FALSE)
system(paste0("head ", sub(".csv$", "_curated.csv", dfname)))

### Signature and GSEA on resting cluster 0n1 single-cell ### %%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/R/clustering_utilities.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
result_id = "signatures_sc0n1_funcpathways_"
signatures_files = c(
  "/home/ciro/asthma_airways/info/signatures_funcpathways.csv"
)

globalsign = list()
signatures_files <- signatures_files[file.exists(signatures_files)]
tvar <- lapply(signatures_files, readfile, stringsAsFactors = FALSE)
tvar <- unlist(lapply(tvar, as.list), recursive = FALSE)
tvar <- sapply(tvar, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
tvar <- tvar[sapply(tvar, length) > 0]
str(tvar)
signatures_list <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
# signatures_list <- signatures_list[!grepl("tcr", names(signatures_list))]
str(signatures_list)
signatures_subset <- clean_feature_list(
  mat = sc_cd4resting_0n1@assays$RNA@data, features = signatures_list, filterby = "p~2", verbose = TRUE
  # ,return_stats =TRUE
)
lapply(names(signatures_subset), function(x){
  head(signatures_list[[x]][!signatures_list[[x]] %in% signatures_subset[[x]]])
})
# head(signatures_subset$stats[c("FAM153B", "FXYD7", "CCL4L2"), , drop = FALSE])

ddf <- vlist2df_diff(
  x = signatures_subset,
  y = signatures_list,
  delim = "Not found or < 2 pct"
)
str(headtail(ddf))
tvar <- which(ddf == "Not found or < 2 pct", arr.ind = TRUE)
ddf[unique(tvar[,1]), unique(tvar[,2])]
write.csv(ddf, file = paste0(result_id, "table.csv"))

source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
sc_cd4resting_0n1 <- signature_scoring(
  object = sc_cd4resting_0n1,
  prefix = paste0(result_id, "modulescore/"),
  lsignatures = signatures_subset,
  confounders = sc_cd4resting_0n1_clust,
  verbose = TRUE
)

### Signature and GSEA on resting single-cell ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/R/clustering_utilities.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
myobject = "sc_cd4resting"
signatures_filter = ""
result_id = "signatures_sc_"
signatures_files = c(
  "/home/ciro/asthma_airways/info/signatures_frombulk.csv"
)
result_id = "signatures_sc_13sets_"
signatures_files = c(
  "/home/ciro/asthma_airways/info/signatures_scresting.csv"
)
result_id = "signatures_sc_funcpathways_"
signatures_files = c(
  "/home/ciro/asthma_airways/info/signatures_funcpathways.csv"
)
result_id = "signatures_sc_ayearoftrmness_"
signatures_files = c(
  "/home/ciro/asthma_airways/info/trmness_combinations.csv"
)

signatures_filter = "Hombrink_Teichmann_2|Cytotoxicity|TREG|All"
result_id = "signatures_sc_stim_ayearoftrmness_"
myobject = "sc_cd4stim"
sc_cd4stim$cluster <- sc_cd4stim@meta.data[, sc_cd4stim_clust]

globalsign = list()
signatures_files <- signatures_files[file.exists(signatures_files)]
tvar <- lapply(signatures_files, readfile, stringsAsFactors = FALSE)
tvar <- unlist(lapply(tvar, as.list), recursive = FALSE)
tvar <- sapply(tvar, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
tvar <- tvar[sapply(tvar, length) > 0]
str(tvar)
signatures_list <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
signatures_list <- signatures_list[grepl(signatures_filter, names(signatures_list))]
str(signatures_list)
if("All_1" %in% names(signatures_list) && signatures_filter == ""){
  sig2elim <- unlist(signatures_list[c("TREG", "Cytotoxicity")], use.names = FALSE)
  signatures_list2 <- lapply(head(signatures_list, -2), function(x) x[!x %in% sig2elim] )
  names(signatures_list2) <- paste0(names(signatures_list2), "_filt")
  signatures_list <- c(signatures_list, signatures_list2)
}
eval(parse(text = paste0(
  "signatures_subset <- clean_feature_list(
  mat = ", myobject, "@assays$RNA@data, features = signatures_list,
  filterby = 'p~2', verbose = TRUE)"
)))
lapply(names(signatures_subset), function(x){
  head(signatures_list[[x]][!signatures_list[[x]] %in% signatures_subset[[x]]])
})
# head(signatures_subset$stats[c("FAM153B", "FXYD7", "CCL4L2"), , drop = FALSE])

ddf <- vlist2df_diff(
  x = signatures_subset,
  y = signatures_list,
  delim = "Not found or < 2 pct"
)
str(headtail(ddf))
tvar <- which(ddf == "Not found or < 2 pct", arr.ind = TRUE)
ddf[unique(tvar[,1]), unique(tvar[,2])]
write.csv(ddf, file = paste0(result_id, "table.csv"))

source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
eval(parse(text = paste0(myobject, " <- signature_scoring(
  object = ", myobject, ",
  prefix = paste0(result_id, 'modulescore/'),
  lsignatures = signatures_subset,
  confounders = c(", myobject, "_clust, 'cluster'),
  verbose = TRUE
)")))

gsea_results <- gsea_matrix(
  mat = expm1(eval(parse(text = paste0(myobject, "@assays$RNA@data")))),
  groups = eval(parse(text = paste0(myobject, "_clust"))),
  metadata = eval(parse(text = paste0(myobject, "@meta.data"))),
  metric = "Signal2Noise",
  gsea_list = signatures_list,
  method = "fgsea",
  path = paste0(result_id, "gsea/"),
  plot_it = TRUE,
  # exp_thr = 0.1,
  classical_plot = !TRUE,
  verbose = TRUE
)
for(i in c(0.05, 0.2)){
  for(j in c(1.5, 1)){
    x <- gsea_plot_summary(
      tests_list = gsea_results,
      path = paste0(result_id, "gsea/"),
      padjthr = i,
      nesthr = j,
      axes = list(
        columns = eval(parse(text = paste0("levels(", myobject, "@meta.data$cluster"))),
        rows = names(signatures_list)
      )
    )
  }
}

df_signatures <- read.csv(paste0(result_id, "modulescore/signatures.csv"), row.names = 1)
eval(parse(text = paste0(
  "mdata = FetchData(", myobject, ", vars = c(redu[[1]], colnames(", myobject, "@meta.data)))"
)))
df_signatures = joindf(df_signatures, mdata)
eval(parse(text = paste0(
  "df_signatures$Cluster = ", myobject, "@meta.data[rownames(df_signatures), 'cluster']"
)))
colours_i <- v2cols(unique(eval(parse(text = paste0(myobject, "@meta.data$cluster")))), colours)
png(paste0(result_id, "heatmap.png"), width = 1600, height = 1700, res = 250)
custom_heatmap(
  object = eval(parse(text = myobject)),
  annoc = df_signatures,
  rnames = signatures_subset[["Hombrink_Teichmann_2"]],
  # use_mean = "Cluster",
  orderby = "HOMBRINK_TEICHMANN_2.Score",
  scale_row = TRUE,
  categorical_col = c("Cluster", "orig.disease"),
  feature_order = "hclust",
  couls = colours_i,
  hcouls = c('yellow', 'black', 'blue'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE, fontsize = 7
)
graphics.off()

scales::percent(c(table(df_signatures$HOMBRINK_TEICHMANN_2.Score > 0) / nrow(df_signatures)))
var2plot = c("HOMBRINK_TEICHMANN_2.Score", "TREG.Score")
meansd <- sapply(
  X = setNames(nm = var2plot),
  FUN = function(x){
    c(
      mean = mean(df_signatures[, x]),
      mean_sd = mean(df_signatures[, x]) + sd(df_signatures[, x]),
      mean_0.5sd = mean(df_signatures[, x]) + (sd(df_signatures[, x]) / 2), # half of standard deviation
      q75 = unname(quantile(df_signatures[, x], probs = .75))
    )
})

for(i in seq(nrow(meansd))){
  p  <- ggplot(
    data = df_signatures,
    mapping = aes_string(x = var2plot[1], y = var2plot[2])
  ) + geom_point(size = 0.5)
  p <- plots_add_quadrants(p, limits = as.list(meansd[i, ]), type = "percent", color = "red")
  p <- p + labs(title = NULL, subtitle = NULL) + SetLegendPointsGG()
  sufix <- paste0(c(var2plot, rownames(meansd[i, , drop = FALSE])), collapse = "_")
  cat(sufix, "\n")
  fname <- paste0(result_id, "modulescore/signatures_scatter_", sufix)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(p)
  dev.off()
}

# Check the distribution
score = df_signatures$HOMBRINK_TEICHMANN_2.Score
p <- ggplot(df_signatures, aes(x=HOMBRINK_TEICHMANN_2.Score)) +
  geom_density(fill = "steelblue3", color="steelblue4") +
  geom_vline(xintercept = mean(score), color = "#BEBEBE", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = mean(score)+(sd(score)/4), y = 2.3, label = "Mean\n50%") +
  annotate("text", x = 0+(sd(score)/3), y = 2.3, label = paste0("+Score\n", scales::percent(mean(score > 0)))) +
  labs(x = NULL, y = "Density", title = "TRM module score distribution")
fname = paste0(result_id, "density")
pdf(paste0(fname, ".pdf"), width = 7, height = 5)
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"), width = 7, height = 5)
p1 <- p; p1$layers <- p1$layers[-2]
print(plot_blank(p1))
dev.off()

# Only the violin
tvar <- "HOMBRINK_TEICHMANN_2\\.|Cluster"
df_plot <- df_signatures[, grep(tvar, colnames(df_signatures)), drop = FALSE]
colnames(df_plot) <- gsub(".Score", "", colnames(df_plot))
df_plot <- reshape2::melt(data = df_plot)
if(length(table(df_plot$variable))>1)
  df_plot$variable <- factor(gsub("_.*", "", df_plot$variable), c("TRM", "NONTRM"))
# df_plot$Cluster <- factor(df_plot$Cluster, sc_cd4resting_ident$order)
head(df_plot)
df_plot$Group <- ident_combine(df_plot, c("Cluster", "variable"))
d2show <- df_plot[sample_even(df_plot, "Group", maxln = -100, v = T), ]

p <- ggplot(df_plot, aes_string(x = "Group", y = "value", fill = "variable", color = "variable")) +
  geom_jitter(data = d2show, size = 0.5, position = position_jitter(0.2), color = 'black', show.legend = FALSE)+
  geom_violin(trim=FALSE, scale = 'width', alpha = 0.8)+
  theme(axis.text.x=element_text(angle=45, vjust = 0.5))
if(nlevels(p$data$variable) > 1){
  p <- p + scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1",guide=FALSE)
}else{
  p <- p + scale_fill_manual(values="steelblue3") + scale_color_manual(values="steelblue4",guide=FALSE)
}
p <- p + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4)
p <- p + labs(x = NULL, y = "Score", fill = "Signature", color = NULL)
fname = paste0(result_id, "violin")
pdf(paste0(fname, ".pdf"), width = 10, height = 5)
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"), width = 8, height = 5)
print(plot_blank(p))
dev.off()

df_signatures$Signature <- df_signatures$HOMBRINK_TEICHMANN_2.Score
p <- ggplot(
  data = df_signatures,
  mapping = aes(x = UMAP_1, y = UMAP_2, color = Signature)
) + geom_point(size = 0.5) +
  scale_colour_gradientn(
    colours = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")
  )
fname = paste0(result_id, "scatter")
pdf(paste0(fname, ".pdf"))
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"))
print(plot_blank(p))
dev.off()

df_signatures$Signature <- ifelse(df_signatures$HOMBRINK_TEICHMANN_2.Score > 0, "Positive", "None")
p <- ggplot(
  data = df_signatures,
  mapping = aes(x = UMAP_1, y = UMAP_2, color = Signature)
) + geom_point(size = 0.2) +
  scale_colour_manual(values = c(Positive = "red", None = "yellow"))
fname = paste0(result_id, "scatter2colours")
pdf(paste0(fname, ".pdf"))
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"))
print(plot_blank(p))
dev.off()

library(dplyr)
summary_donor_df = df_signatures %>%
  group_by(orig.donor) %>%
  summarise(
    ScoreMean = round(mean(HOMBRINK_TEICHMANN_2.Score), 4),
    ScoreMedian = round(median(HOMBRINK_TEICHMANN_2.Score), 4),
    ScorePercentOver0 = scales::percent(mean(HOMBRINK_TEICHMANN_2.Score > 0), 4)
  )
write.csv(summary_donor_df, file = paste0(result_id, "donor_stats.csv"), row.names = FALSE)

### GSEA on bulk ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/R/clustering_utilities.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
fannog = "/home/ciro/asthma_airways/info/deprecated/gene_annotation_GRCh37.csv"
axes_order = list(
  columns = c("CD4-TRM", "CD4-TRM-like", "CD4-Teff", "CD4-Treg"),
  rows = c("cd4_trm_vieira", "cd4_cd103p_oja1", "cd4_tmc_vieira", "cd4_blood_oja2", "cd4_tfh_locci", "cd4_treg_schmiedel")
)

result_id = "signatures_bulkresting_"
sselect = list(c('Cell_type', '-CD4-TFH', '-CD8-TRM'), c('Stim', 'Unstim'))

result_id = "signatures_bulkstim_"
sselect = list(c('Cell_type', '-CD4-TFH', '-CD8-TRM'), c('Stim', 'Stim'))

annog <- read.csv(fannog, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
thesegenes <- features_parse_ensembl(rownames(annog[annog$gene_type_4 == "protein_coding", ]))

tmp_tab <- read.csv("/home/ciro/asthma_airways/info/signatures_2020_12_03.csv", stringsAsFactors = FALSE)
mysignatures <- lapply(X = tmp_tab, FUN = function(x) x[x!=""] )
mysignatures <- mysignatures[!grepl("tfh", names(mysignatures))]
str(mysignatures)

ddf <- vlist2df_diff(
  x = lapply(mysignatures, function(x) x[x %in% thesegenes] ),
  y = mysignatures,
  delim = "Not found"
)
str(headtail(ddf))
tvar <- which(ddf == "Not found", arr.ind = TRUE)
ddf[unique(c(tvar[,1], nrow(ddf))), unique(c(tvar[,2], ncol(ddf)))]
write.csv(ddf, file = paste0(result_id, "table.csv"))

# file.remove(list.files(path = result_id, pattern="a1", full.names=T))
gsea_results <- gsea_matrix(
  exp_thr = 0,
  mat = bulk_tpm[features_parse_ensembl(rownames(bulk_tpm)) %in% thesegenes, ], # count matrix or a table of metrics per comparison as columns
  groups = "Cell_type", # just the column and it'll do vs REST
  metadata = bulk_metadata[filters_subset_df(sselect, bulk_metadata, v = T), ],
  metric = "Signal2Noise",
  gsea_list = mysignatures,
  method = "fgsea",
  path = paste0(result_id, "gsea/"),
  plot_it = TRUE,
  classical_plot = !TRUE,
  verbose = TRUE
)
x <- gsea_plot_summary(
  tests_list = gsea_results,
  path = paste0(result_id, "gsea/"),
  padjthr = 0.05,
  nesthr = 1,
  axes = axes_order
)
x <- gsea_plot_summary(
  tests_list = gsea_results,
  path = paste0(result_id, "gsea/"),
  padjthr = 0.05,
  nesthr = 1.5,
  axes = axes_order
)

### Specific volcanos ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/volcano.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/devel/overlap.R")

padjthr = 0.05; fcthr = 1
result_id = "volcano_bulk/"; dir.create(result_id)
top_mean = 9; top_lfc = 6
further_thr = "BmeanTPM>5"
volc_configs = list(
  list(
    trimmer = "c15",
    showgenes = c("ACACB", "ALOX12P2", "CCDC61", "MCAM", "SNORA31", "CREM", "RGS1", "TNFAIP3", "PDE4B", "JUNB", "CXCR4", "FOSL2", "CD69", "C1QB"),
    resname = "resting_trmlike",
    resfile = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/resting_disease/CD4-TRM-like_SAvsCD4-TRM-like_MA/results_CD4-TRM-like_SAvsCD4-TRM-like_MA_deseq2.csv"
  ),
  list(
    trimmer = "c15",
    showgenes = c("CCL4L2", "FZD6", "HOPX", "NLRC5", "ZNF683", "TRAF1", "MYO7A", "DUOX1", "CREM", "DUSP1", "DUSP4", "RGS1", "TNFAIP3", "TOX", "PDE4B", "CXCR4", "CD28", "CTSH"),
    resname = "resting_trm",
    resfile = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/resting_disease/CD4-TRM_SAvsCD4-TRM_MA/results_CD4-TRM_SAvsCD4-TRM_MA_deseq2.csv"
  ),
  list(
    trimmer = "c15",
    showgenes = c("CCL16", "CTSF", "DDIT4", "EGR1", "EGR2", "JUNB", "KLF6", "TXNIP", "GZMA", "PLAC8", "XCL1", "BCAR3", "SOX4", "TXN"),
    resname = "stim_trmlike",
    resfile = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/stim_disease/CD4-TRM-like_SAvsCD4-TRM-like_MA/results_CD4-TRM-like_SAvsCD4-TRM-like_MA_deseq2.csv"
  ),
  list(
    trimmer = "c15",
    showgenes = c("AMICA1", "BIRC5", "CCL4", "CCL4L1", "CCL4L2", "CTSW", "FKBP5", "IL12A", "IL24", "ITGAX", "TNFSF11", "ZBTB16", "IL17F", "CTSF",
    "CTSL", "ADAM12", "IL20RB", "IL31", "IL37", "MMP19", "PLAUR", "XCL1"),
    resname = "stim_trm",
    resfile = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/stim_disease/CD4-TRM_SAvsCD4-TRM_MA/results_CD4-TRM_SAvsCD4-TRM_MA_deseq2.csv"
  ),
)

padjthr = 0.05; fcthr = 0.25
result_id = "volcano_sc/"; dir.create(result_id)
top_mean = NULL; top_lfc = NULL
further_thr = NULL
volc_configs = list(
  list(
    trimmer = "c400",
    showgenes = NULL,
    resname = "clusters_cov_sex_0vs1", size_feature = "pcts",
    resfile = "/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_cov_sex/0vs1/results_0vs1_mastlog2cpm.csv"
  ),
  list(
    trimmer = "c75",
    fc_lim = 2,
    showgenes = c("ITGAE", "ITGA1", "GZMB", "GZMH", "PRF1", "ZNF683", "TRAF3IP3", "IFITM1", "IFNG", "FKBP5", "CXCR3", "TNFSF14", "CREM", "TNFAIP3", "CXCR4", "METRNL", "DUSP1", "DUSP4"),
    resname = "clust0n1_SA_MalevsSA_Female", size_feature = "pcts",
    resfile = "/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clust0n1/SA_MalevsSA_Female/results_SA_MalevsSA_Female_mastlog2cpm.csv"
  )
)

for(config_i in volc_configs[2]){
  cat(basename(config_i$resname), "\n")
  res <- readfile(config_i$resfile, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  dtype <- sub("Bmean", "", colnames(res)[grepl("Bmean", colnames(res))])
  grps <- unlist(strsplit(gsub("results_|_mast.*|_dese.*", "", basename(config_i$resfile)), "vs"))
  tvar <- paste0(grps, "_mean", dtype)
  res_filt <- data.frame(res[!is.na(res$log2FoldChange), ], stringsAsFactors = F, check.names = FALSE)
  res_filt$padj[res_filt$padj>1] <- 1
  rownames(res_filt) <- features_parse_ensembl(rownames(res_filt))
  res_filt$gene <- features_parse_ensembl(gsub("'", "", res_filt$gene))
  res_filt$Mean <- round(log2(ifelse(res_filt$group == grps[1], res_filt[, tvar[1]], res_filt[, tvar[2]]) + 1), 1)
  if(!is.null(top_mean)) res_filt$Mean[which(res_filt$Mean > top_mean)] <- top_mean
  res_filt$pcts <- res_filt$pct_diff
  genes2plot <- mysignames <- getDEGenes(res_filt, pv = padjthr, fc = fcthr, further = further_thr, v = TRUE)
  if(is.null(further_thr)){
    tvar <- rownames(res_filt)[res_filt[, grep("minExp", colnames(res_filt), value = TRUE)]]
    genes2plot <- genes2plot[genes2plot %in% tvar]
  }
  res_filt$degs <- "Not_significant"
  res_filt[res_filt$gene %in% genes2plot, ]$degs <- "DEG"
  res_filt[!res_filt$gene %in% genes2plot, ]$Mean <- NA
  res_filt[!res_filt$gene %in% genes2plot, ]$pcts <- 0
  tvar <- cosmy(genes2plot, patties = "^rps|^rpl|^mt-|rp[0-9]{1,}-|^linc")
  showgenes <- if(is.null(config_i$showgenes) || any(config_i$showgenes %in% "add")){
    c(config_i$showgenes, bordering(res_filt[tvar, ], cnames = "log2FoldChange", n = 10))
  }else if(is.null(config_i$showgenes)){ FALSE }else{ config_i$showgenes }
  showgenes <- show_found(showgenes, rownames(res_filt), v = TRUE)
  tvar <- -log10(res_filt$padj); tvar[is.infinite(tvar)] <- max(tvar[is.finite(tvar)]); summary(tvar)
  trimit <- ifelse(isTRUE(max(tvar) < as.numeric(gsub("c", "", config_i$trimmer))), paste0("c", round(max(tvar))), config_i$trimmer)
  if(!is.null(config_i$fc_lim)){
    res_filt$log2FoldChange[res_filt$log2FoldChange>config_i$fc_lim] <- config_i$fc_lim
    res_filt$log2FoldChange[res_filt$log2FoldChange<(-config_i$fc_lim)] <- -config_i$fc_lim
  }
  void <- volplot(
    res_filt,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = 'padj',
    lfctype = 'log2FoldChange',
    col_feature = "Mean",
    size_feature = config_i$size_feature,
    gene_name = 'gene',
    check_genes = list(text = showgenes),
    titl = paste0("'", paste0(gsub("_", " ", grps), collapse = "' vs '"), "'"),
    return_plot = TRUE,
    clipp = trimit,
    v = TRUE
  ) + labs(
    size = "Delta %", color = "Group mean",
    title = NULL, subtitle = paste("DEGs:", length(genes2plot))
  )
    # scale_size(range = c(0, 7), limits = c(0, 100)) +
  if(!is.null(top_lfc)) void <- void + scale_x_continuous(n.breaks = 7, limits = c(-top_lfc, top_lfc))
  fname <- paste0(result_id, config_i$resname, "_fc", fcthr, make.names(further_thr), "_trim", trimit)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(void); if(is.null(config_i$size_feature)) update_geom_defaults("point",list(size=3))
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
  print(plot_blank(void))
  dev.off()
}

### Heatmap sc cluster - disease ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/devel/overlap.R")
result_id = "heatmap_sc_"
sig_thresh = list(pv = 0.05, fc = 0.5, fu = "pct25exp == TRUE")

color_by_cols = c(sc_cd4resting_clust, "orig.asthma")
res_file = c(
  SA_0vsMA_0 = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease_clusters/SA_0vsMA_0/results_SA_0vsMA_0_mastlog2cpm.csv",
  SA_1vsMA_1 = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease_clusters/SA_1vsMA_1/results_SA_1vsMA_1_mastlog2cpm.csv"
)
sselect = list(c(sc_cd4resting_clust, '0', '1'))

color_by_cols = c("orig.Sex", "orig.asthma", sc_cd4resting_clust)
res_file = c(
  "0_Male_SAvs0_Male_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/0_Male_SAvs0_Male_MA/results_0_Male_SAvs0_Male_MA_mastlog2cpm.csv",
  "0_Female_SAvs0_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/0_Female_SAvs0_Female_MA/results_0_Female_SAvs0_Female_MA_mastlog2cpm.csv"
)
sselect = list(c(sc_cd4resting_clust, '0'))

color_by_cols = c("orig.Sex", "orig.asthma", sc_cd4resting_clust)
res_file = c(
  "1_Male_SAvs1_Male_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/1_Male_SAvs1_Male_MA/results_1_Male_SAvs1_Male_MA_mastlog2cpm.csv",
  "1_Female_SAvs1_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/1_Female_SAvs1_Female_MA/results_1_Female_SAvs1_Female_MA_mastlog2cpm.csv"
)
sselect = list(c(sc_cd4resting_clust, '1'))

dir.create("clusters_female_disease")
result_id = "clusters_female_disease/heatmap_sc_"
color_by_cols = c("orig.Sex", "orig.asthma", sc_cd4resting_clust)
res_file = c(
  "0_Female_SAvs0_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/0_Female_SAvs0_Female_MA/results_0_Female_SAvs0_Female_MA_mastlog2cpm.csv",
  "1_Female_SAvs1_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/1_Female_SAvs1_Female_MA/results_1_Female_SAvs1_Female_MA_mastlog2cpm.csv"
)
sselect = list(c(sc_cd4resting_clust, '1', '0'), c("orig.Sex", "Female"))

fname <- paste0(result_id, paste0(basename(dirname(res_file)), collapse = "-"))
sc_cd4resting$orig.Sex = factor(sc_cd4resting$orig.Sex, c("Male", "Female"))
table(sc_cd4resting$orig.Sex)
sc_cd4resting$Identity = ident_combine(sc_cd4resting@meta.data, rev(color_by_cols), "_")
table(sc_cd4resting$Identity)
table(droplevels(sc_cd4resting@meta.data[ssamples, ]$Identity))
res_list = lapply(X = res_file, FUN = function(x) {
  y <- readfile(x, row.names = 1, stringsAsFactors = FALSE)
  y$pct25exp = y[, grepl("minExp", colnames(y))]; y
})
genes_list0 = unlist(lapply(X = res_list, FUN = function(x){
  list(
    up = getDEGenes(x = x, pv = sig_thresh$pv, fc = sig_thresh$fc, upreg = TRUE, further = sig_thresh$fu, verbose = !TRUE),
    down = getDEGenes(x = x, pv = sig_thresh$pv, fc = sig_thresh$fc, upreg = FALSE, further = sig_thresh$fu, verbose = !TRUE)
  )
}), recursive = FALSE)
names(genes_list0) <- gsub("(.*)vs.*\\.up", "\\1", names(genes_list0))
names(genes_list0) <- gsub(".*vs(.*)\\.down", "\\1", names(genes_list0))
str(genes_list0)
genes_list <- overlap_list(genes_list0); genes_list <- genes_list[sapply(genes_list, length) > 2]
tvar <- c("MA_0", "SA_0", "MA_1", "SA_1", paste0(rep(c(0:1), each=4), "_", c("Male_MA", "Male_SA", "Female_MA", "Female_SA")))
tvar <- tvar[tvar %in% names(genes_list)]
genes_list <- c(genes_list[tvar], genes_list[!names(genes_list) %in% tvar])
str(genes_list)
annor = reshape2::melt(genes_list)
annor <- data.frame(Genes = factor(annor$L1, names(genes_list)), row.names = as.character(annor$value))
write.csv(annor, file = paste0(fname, ".csv"))
ssamples <- filters_subset_df(sselect, sc_cd4resting@meta.data, v = T)
# ssamples <- sample_even(sc_cd4resting@meta.data[ssamples, ], cname = "Identity", maxln = -1000, v = TRUE)
tvar <- unique(unlist(lapply(sc_cd4resting@meta.data[ssamples, c(color_by_cols, "Identity")], as.character)))
colours_i <- v2cols(unique(c(tvar, levels(annor$Genes))), colours)
unique_name = function(x) unname(sapply(x, function(y) paste0(sort(unlist(strsplit(as.character(y), "_"))), collapse = "_") ))
colours_i = sapply(setNames(nm = names(colours_i)), function(x){
  tail(unname(colours_i[unique_name(names(colours_i)) == unique_name(x)]), 1)
})

edata <- expm1(as.matrix(sc_cd4resting@assays$RNA@data))
png(paste0(fname, ".png"), width = 1600, height = 1700, res = 250)
custom_heatmap(
  object = edata[, ssamples],
  annoc = sc_cd4resting@meta.data[ssamples, ],
  rnames = unlist(genes_list),
  use_mean = "Identity",
  orderby = "Identity", # color_by_cols
  scale_row = TRUE,
  categorical_col = c(color_by_cols, "Identity"),
  feature_order = FALSE,
  couls = colours_i,
  hcouls = c('yellow', 'black', 'blue'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annor,
  # gaps_row = unname(cumsum(sapply(genes_list, length))),
  # gaps_col = unname(cumsum(table(sc_cd4resting@meta.data[ssamples, "Identity"])[levels(sc_cd4resting@meta.data$Identity)]))
  annotation_names_col = FALSE, annotation_names_row = FALSE,
  do_log = TRUE
)
graphics.off()

### Heatmap sc resting clusters_sex_disease ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "clusters_female_disease/"
dir.create(result_id)
sig_thresh = list(pv = 0.05, fc = 0.25)
res_list_f = list.files(
  path = "/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease",
  pattern = "results_.*csv", full.names = TRUE, recursive = TRUE
)
res_list_f = grep(res_list_f, pattern = "Female", value = TRUE)
names(res_list_f) <- basename(dirname(res_list_f))

mdata = sc_cd4resting@meta.data
grouping_vars = c("cluster", "orig.Sex", "orig.asthma")
mdata$Group = ident_combine(mdata, grouping_vars, "_")
sselect = c("Group", unique(unlist(strsplit(names(res_list_f), "vs"))))
scells = filters_subset_df(sselect, mdata, v = T)
scells = sample_even(mdata[scells, ], cname = "Group", maxln = -500, v = T)
mdata = mdata[scells, ]
res_list = lapply(X = res_list_f, FUN = function(x){
  y <- readfile(x, stringsAsFactors = FALSE); rownames(y) <- y$gene
  y
})
source("/home/ciro/scripts/functions/group_specificity.R")
res_list_sp = g_sp(
  comps_stats = res_list,
  padjthr = sig_thresh$pv,
  fcthr = sig_thresh$fc,
  expr_mat = as.matrix(sc_cd4resting@assays$RNA@data[, rownames(mdata)]),
  annotation = mdata,
  cname = "Group",
  verbose = TRUE
)
mygenes = make_list(
  res_list_sp$summary[!is.na(res_list_sp$summary$group), ],
  "group", "rownams"
)
mdata$Group = droplevels(mdata$Group)
mdata$cluster = droplevels(mdata$cluster)
mygenes = mygenes[c(
  levels(mdata$Group)[levels(mdata$Group) %in% names(mygenes)],
  names(mygenes)[!names(mygenes) %in% levels(mdata$Group)]
)]

annor = reshape2::melt(mygenes)
annor <- data.frame(Genes = factor(annor$L1, names(mygenes)), row.names = as.character(annor$value))
tvar <- unique(unlist(lapply(mdata[, c("Group", grouping_vars)], as.character)))
colours_i <- v2cols(unique(c(tvar, levels(annor$Genes))), colours)

source("/home/ciro/scripts/handy_functions/devel/plots.R")
fname <- paste0(result_id, "fc", sig_thresh$fc, "_padj", sig_thresh$pv, "_")
png(paste0(fname, "heatmap.png"), width = 1500, height = 1700, res = 250)
custom_heatmap(
  object = sc_cd4resting@assays$RNA@data,
  annoc = mdata,
  rnames = unlist(mygenes),
  cnames = scells,
  orderby = "Group",,
  use_mean = "Group",
  scale_row = TRUE,
  categorical_col = c("Group", grouping_vars),
  # feature_order = "hclust",
  couls = colours_i,
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  gaps_row = unname(cumsum(sapply(mygenes, length))),
  annotation_row = annor
)
graphics.off()
write.csv(annor, file = paste0(fname, ".csv"))

### Heatmap bulk 4 cell types ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
result_id = "heatmap_bulk_"
color_by_cols = "aCell_type"
order_celltype = c("CD4-TRM", "CD4-TRM-like", "CD4-Teff", "CD4-Treg")

res_file = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/summary/resting4celltypes_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"
sselect = list(c('Cell_type', '-CD8-TRM', '-CD4-TFH'), c('Stim', 'Unstim'))
res_file = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/summary/stim4celltypes_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"
sselect = list(c('Cell_type', '-CD8-TRM', '-CD4-TFH'), c('Stim', 'Stim'))
res_file = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/summary/resting_trmness_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"
sselect = list(c('Cell_type', 'CD4-TRM', 'CD4-TRM-like', 'CD4-Teff'), c('Stim', 'Unstim'))
res_file = "/home/ciro/large/asthma_airways/results/deseq2/airways_platex/summary/stim_trmness_padj0.05_FC1/Cell_type_SummaryDEGsTable_suas.csv"
sselect = list(c('Cell_type', 'CD4-TRM', 'CD4-TRM-like', 'CD4-Teff'), c('Stim', 'Stim'))
color_by_cols = c("aCell_type", "Disease")

res = readfile(res_file, row.names = 1, stringsAsFactors = FALSE)
genes_list <- make_list(res[!is.na(res$group), ], colname = "group", col_objects = "gene_name")
genes_list <- lapply(X = genes_list, FUN = sub, pattern = "'", replacement = "")
tvar <- genes_list[grepl("n", names(genes_list))]
genes_list <- c(
  genes_list[!grepl("n", names(genes_list))],
  tvar[names(sort(sapply(tvar, length), dec = T))]
)

if(grepl("4celltypes", res_file)){
  genes_list <- c(
    genes_list[names(genes_list) %in% "CD4-TRMnCD4-TRM-like"],
    genes_list[!names(genes_list) %in% "CD4-TRMnCD4-TRM-like"]
  )
  color_by_cols = c("mCell_type", "aCell_type")
  bulk_metadata$mCell_type <- factor(bulk_metadata$Merged_Cell_type, c("CD4-TRM", "CD4-nonTRM", "CD4-Treg"))
}
ssamples <- filters_subset_df(sselect, bulk_metadata, v = T)
tvar <- as.character(bulk_metadata[ssamples, ]$Cell_type)
bulk_metadata$aCell_type <- factor(bulk_metadata$Cell_type, order_celltype[order_celltype %in% tvar])
annor = reshape2::melt(genes_list)
annor <- data.frame(Genes = factor(annor$L1, names(genes_list)), row.names = as.character(annor$value))

tvar <- unique(unlist(lapply(bulk_metadata[ssamples, color_by_cols], as.character)))
colours_i <- v2cols(unique(c(tvar, levels(annor$Genes))), colours)
# cat(paste0(names(colours_i), ",", colours_i, collapse = "\n"), "\n")
# tvar <- paste0(sapply(1:length(colours_i), function(x) paste0(names(colours_i)[x], "' = '", colours_i[x]) ), collapse = "', '")
# tvar <- paste0("colours_i = c('", tvar, "')")
fname <- paste0(result_id, basename(dirname(res_file)))
# pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
png(paste0(fname, ".png"), width = 1500, height = 1700, res = 250)
custom_heatmap(
  object = bulk_tpm,
  annoc = bulk_metadata[ssamples, ],
  rnames = unlist(genes_list),
  orderby = color_by_cols,
  scale_row = TRUE,
  categorical_col = color_by_cols,
  feature_order = FALSE,
  couls = colours_i,
  hcouls = c('yellow', 'black', 'blue'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annor,
  gaps_row = unname(cumsum(sapply(genes_list, length))),
  gaps_col = unname(cumsum(table(bulk_metadata[ssamples, "aCell_type"])[levels(bulk_metadata$aCell_type)]))
)
graphics.off()

### Heatmap sc stim disease ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "stim_disease_degs/"
if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
sig_thresh = list(pv = 0.05, fc = 0.25)
res_list_f = list.files(
  path = "/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/stim_disease",
  pattern = "results_.*csv", full.names = TRUE, recursive = TRUE
)

mdata = sc_cd4merged@meta.data
mdata$Group = ident_combine(mdata, c("orig.stim", "orig.asthma"), "_")
scells = sample_even(mdata, cname = "Group", maxln = -500, v = T)
mdata = mdata[scells, ]
names(res_list_f) <- basename(dirname(res_list_f))
res_list = lapply(X = res_list_f, FUN = function(x){
  y <- readfile(x, stringsAsFactors = FALSE); rownames(y) <- y$gene
  y
})
source("/home/ciro/scripts/functions/group_specificity.R")
res_list_sp = g_sp(
  comps_stats = res_list,
  padjthr = sig_thresh$pv,
  fcthr = sig_thresh$fc,
  expr_mat = as.matrix(sc_cd4merged@assays$RNA@data[, rownames(mdata)]),
  annotation = mdata,
  cname = "Group",
  verbose = TRUE,
  path_plot = result_id
)

### Violin plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")

myobject = "sc_cd4resting_0n1"; group_i = sc_cd4resting_0n1_clust; sselect = NULL
result_id = "violin_trajectory_derived/"; dir.create(result_id)
genes = c(
  "CCR6", "DUSP1", "FOS", "FOSB", "CREM", "DUSP2", "CCNH", "TNFAIP3", "FRMD4B",
  "MAF", "PDE4B", "CXCR4", "DUSP4", "GZMK", "CTSH", "FKBP5", "PDE4D", "CEBPD",
  "CHN1", "GZMA", "CCR4", "NR4A2", "CNOT6L", "SYTL3", "ZNF331", "CD28", "NR3C1",
  "STAT4", "PDE4A", "PRDM1", "RGS1", "ZNF683", "XCL2", "CCL5", "GZMH", "SOCS3",
  "BATF", "BCL3", "NKG7", "ITGAE", "CXCR6", "CCL4", "HOPX", "PLEKHF1", "IFNG",
  "GZMB", "ITGA1", "XCL1", "IL10RA", "TNFSF10", "ZFAS1", "BCL2", "TXNIP", "CXCR3"
)

myobject = "sc_cd4resting"; group_i = "cluster"; sselect = list(c("cluster", "0", "1"))
result_id = "violin_clust0n1/"; dir.create(result_id)
result_id = "violin_clust0n1_chilli/"; dir.create(result_id)
genes = c(
  "TBX21", "IFNG", "HOPX", "GZMB", "GZMH", "PRF1", "ZNF683", "CXCR6", "CCL4",
  "TNFSF14", "FKBP5", "DDIT4", "TRADD", "CREM", "TNFAIP3", "GATA3", "STAT4",
  "CCR4", "CCR6", "CXCR4", "CTLA4", "ICOS", "MAF", "DUSP4", "DUSP2", "TNFSF10"
)

sc_cd4resting$tmp <- factor(sc_cd4resting$cluster, rev(levels(sc_cd4resting$cluster)))
sc_cd4resting$cluster.disease = ident_combine(sc_cd4resting[[]], c("tmp", "orig.asthma"))
table(sc_cd4resting$cluster.disease)
myobject = "sc_cd4resting"; group_i = "cluster.disease"; sselect = list(c("cluster", "0", "1"))
result_id = "violin_clust0n1disease/chilli_"; dir.create("violin_clust0n1disease")
result_id = "violin_clust0n1disease_female/chilli_"; dir.create("violin_clust0n1disease_female")
sselect = list(c("cluster", "0", "1"), c("orig.Sex", "Female"))
genes = c(
  "IFITM1", "CXCR3", "TRAF3IP3", "TNFSF14", "FKBP5", "CREM", "TNFAIP3", "METRNL", "DUSP4"
)

mygenes <- show_found(genes, rownames(eval(parse(text = myobject))), verbose = TRUE)
tvar <- c(group_i, "cluster", "orig.donor", "orig.asthma", sapply(sselect, head, 1), mygenes)
pdata <- FetchData(
  object = eval(parse(text = myobject)),
  vars = unique(unlist(tvar))
)
if(!is.null(sselect)){
  pdata <- pdata[filters_subset_df(sselect, pdata, v = TRUE), ]
  pdata[, group_i] <- droplevels(pdata[, group_i])
}
for(g in mygenes){
  cat(g, "\n")
  p <- violin(
    dat = pdata,
    xax = group_i,
    yax = g,
    dots = grepl("dots", result_id),
    colour_by = "pct",
    chilli = grepl("chilli", result_id)
  ) + labs(y = "Seurat Normalised")

  # fname <- paste0(result_id, g)
  fname <- paste0(result_id, "boxless_", g); p <- plot_rm_layer(p, "Box")
  pdf(paste0(fname, ".pdf"))
  print(p)
  graphics.off()
  pdf(paste0(fname, "_blank.pdf"))
  print(plot_blank(p))
  graphics.off()
}

source("/home/ciro/scripts/handy_functions/devel/plots_dotplot.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R")
for(scale_mean_q in c(FALSE, TRUE)){
  p <- dot_plot(
    edata = as.matrix(t(pdata[, -c(1:3)])),
    mdata = pdata[, 1:3],
    columns = "cluster",
    scale_mean = scale_mean_q,
    features = mygenes,
    cols_limits = if(scale_mean_q) c(0, 0.6),
    values_norm = if(scale_mean_q) c else function(x) log2(x+1),
    clust_rows = !TRUE, cols = c('#fffef0', '#ff0000'), values_exprthr = 0
  )
  fname <- paste0(result_id, "_dotplot", ifelse(scale_mean_q, "_mean_scaled", ""))
  # if(!scale_mean_q) next
  # p <- DotPlot(
  #   object = sc_cd4resting[, rownames(pdata)],
  #   features = mygenes,
  #   group.by = "cluster",
  #   # cols = c('#fff4ba', '#ff0000'), col.min = -1.5, col.max = 1.5,
  #   cols = c('#fffef0', '#ff0000'), col.min = 0, col.max = 1.5,
  #   dot.min = 0.1
  # ) + coord_flip() +
  #   theme(
  #     axis.text.y = element_text(size = 13, face = "bold.italic"),
  #     axis.ticks.x = element_blank()
  #   ) + labs(y = NULL, x = NULL)
  # p$guides$colour$title <- "Z-scored\nAverage\nExpression"
  # p$guides$size$title <- "Percent\nExpressed"
  # fname <- paste0(result_id, "_seurat_dotplot", ifelse(scale_mean_q, "_mean_scaled", ""))
  cat(fname, "\n")
  pdf(paste0(fname, ".pdf"), width = 6)
  print(p)
  graphics.off()
  pdf(paste0(fname, "_blank.pdf"), width = 5)
  print(plot_blank(p))
  graphics.off()
}

pdf(paste0(result_id, "_heatmap.pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = as.matrix(t(pdata[, -c(1:3)])),
  annoc = pdata,
  rnames = mygenes,
  orderby = "cluster",,
  use_mean = c("cluster", "orig.donor"),
  scale_row = TRUE,
  categorical_col = c("cluster", "orig.asthma")[1],
  feature_order = "hclust",
  couls = colours,
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols_override = TRUE
)
graphics.off()

### Scatter pairs bulk ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
library(GGally)
dir.create("coexpression_pairs")
fconfigs = list(
  list(
    result_id = "coexpression_pairs/bulk",
    edata = "bulk_tpm",
    mdata = "bulk_metadata",
    features = c("IFNG", "IL4", "IL5", "IL13", "IL17A"),
    center_names = c("T[H]*1", "T[H]*2", "T[H]*2", "T[H]*2", "T[H]*17"),
    color_by = c("Cell_type", "Stim"), # c("Sex", "Disease")
    sselect = list(c('Cell_type', '-CD4-TFH', '-CD8-TRM')),
    blimit = list(log2(10 + 1), log2(10 + 1)), size = 15
  ),
  list(
    result_id = "coexpression_pairs/sc",
    edata = "sc_cd4stim@assays$RNA@data",
    mdata = "sc_cd4stim@meta.data",
    features = c("IFNG", "TNF", "XCL1", "IL17A", "IL21", "IL23A", "CCL20", "GZMB", "GZMA", "GZMH", "CCL3", "CCL4", "IL4", "IL5", "IL13"),
    center_names = c(paste0("X", 1:15)),
    color_by = c("cd4stim", "cd4stim"),
    sselect = "sample_even",
    blimit = list(log(0 + 1), log(0 + 1)), size = 30
  )
)

for(fconfig in fconfigs){}
allgenes <- rownames(eval(parse(text = fconfig$edata)))
mygenes <- show_found(fconfig$features, features_parse_ensembl(allgenes), v = TRUE)
mygenes <- allgenes[features_parse_ensembl(allgenes) %in% fconfig$features]
ssamples <- if("sample_even" %in% fconfig$sselect){
  sample_even(eval(parse(text = fconfig$mdata)), maxln = 0.2, v = T)
}else{
  filters_subset_df(fconfig$sselect, eval(parse(text = fconfig$mdata)), v = T)
}

df2plot = t(as.matrix(eval(parse(text = fconfig$edata))[mygenes, ssamples]))
if(fconfig$edata == "bulk_tpm") df2plot = log2(df2plot + 1)
df2plot =  data.frame(df2plot)
colnames(df2plot) <- features_parse_ensembl(colnames(df2plot))
df2plot <- df2plot[, fconfig$features]
if(any(fconfig$color_by %in% colnames(eval(parse(text = fconfig$mdata)))))
  df2plot = cbind(eval(parse(text = fconfig$mdata))[ssamples, fconfig$color_by, drop = F], df2plot)
head(df2plot)
pp <- list()
for(mycol in unique(fconfig$color_by)){
  pp[[mycol]] <- ggpairs(
    data = df2plot,
    mapping = if(mycol %in% df2plot) ggplot2::aes_string(colour=mycol) ,
    columns = fconfig$features,
    upper = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na")
  ) + theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "#f5f5f5")
  )
}
pp_blank <- pp
for(i in 1:pp[[1]]$nrow) { for(j in 1:pp[[1]]$ncol){
  mycol <- fconfig$color_by[ifelse(i > j, 1, 2)]
  p <- pp[[mycol]][i,j]
  if(i != j && !is.null(fconfig$blimit))
    p <- plots_add_quadrants(p, limits = fconfig$blimit, type = "percent", center = TRUE)
  if(any(mycol %in% colnames(df2plot)))
    pp[[1]][i,j] <- p + scale_color_manual(values = v2cols(df2plot[, mycol], colours))
  if(i == j){
    pp[[1]][i,j] <- ggplot() + theme_void() +
      annotate("text", x = 4, y = 25, size = 8, label = parse(text=center_names[i]))
  }
  pp_blank[[1]][i,j] <- plot_blank(pp[[1]][i,j]) + if(i == j) theme_nothing()
}}

fname <- paste0(c(fconfig$result_id, fconfig$color_by), collapse = "_")
pdf(paste0(fname, ".pdf"), width = fconfig$size, height = fconfig$size)
print(pp[[1]])
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = fconfig$size, height = fconfig$size)
print(plot_blank(pp_blank[[1]]))
graphics.off()

### Disease genes 0n1 subclustering ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_thresh = list(pv = 0.05, fc = 0.25)

result_id = "subclust0n1_disease_specific_"
result_id = "subclust0n1_disease_specific_fc0.5_"; sig_thresh$fc = 0.5
res_list_f = list.files(
  path = "/home/ciro/large/asthma_airways/results/scdgea/cluster0n1_20p/comprs/cluster_disease",
  pattern = "results_.*csv", full.names = TRUE, recursive = TRUE
)
mcols = c(sc_cd4resting_0n1_clust)
mdata = sc_cd4resting_0n1@meta.data
order_celltype = NULL

sex = "Female"
sex = "Male"
result_id = paste0("cluster0n1_diseasex_", sex, "_")
res_list_f = list.files(
  path = "~/ad_hoc/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease",
  pattern = "results_.*csv", full.names = TRUE, recursive = TRUE
); res_list_f <- grep(sex, res_list_f, value = TRUE)
mcols = c(sc_cd4resting_clust, "orig.Sex", "orig.asthma")
scells = filters_subset_df(list(c("orig.Sex", sex), c(sc_cd4resting_clust, "0", "1")), sc_cd4resting@meta.data, v = TRUE)
mdata = sc_cd4resting@meta.data[scells, ]
order_celltype = paste0(rep(c(0:1), each = 2), "_", sex, "_", rep(c("MA", "SA"), 2))

names(res_list_f) <- basename(dirname(res_list_f))
res_list = lapply(X = res_list_f, FUN = function(x){
  y <- readfile(x, stringsAsFactors = FALSE); rownames(y) <- y$gene
  y
})
mdata$Group = ident_combine(mdata, mcols, "_")
table(mdata$Group)

source("/home/ciro/scripts/functions/group_specificity.R")
res_list_ma = g_sp(
  comps_stats = res_list[!names(res_list_f) %in% c("0vs1", "0~1vs2~3~4~5")],
  padjthr = sig_thresh$pv,
  fcthr = sig_thresh$fc,
  expr_mat = as.matrix(sc_cd4resting_0n1@assays$RNA@data),
  annotation = mdata,
  cname = "Group",
  verbose = TRUE
)

res_list_u <- res_list_ma_u <- make_list(
  x = res_list_ma$summary[!is.na(res_list_ma$summary$group), ],
  colname = "group", col_objects = "rownams"
)
if("0~1vs2~3~4~5" %in% names(res_list_f)){
  source("/home/ciro/scripts/handy_functions/devel/filters.R")
  res_list_u0 = unlist(lapply(
    X = res_list[1:2],
    FUN = function(x){
      list(
        up = getDEGenes(x = x, pv = sig_thresh$pv, fc = sig_thresh$fc, upreg = TRUE, further = NULL, verbose = !TRUE),
        down = getDEGenes(x = x, pv = sig_thresh$pv, fc = sig_thresh$fc, upreg = FALSE, further = NULL, verbose = !TRUE)
      )
    }
  ), recursive = FALSE)
  names(res_list_u0) <- gsub("vs.*up", "", names(res_list_u0))
  names(res_list_u0) <- gsub(".*vs(.*).down", "\\1", names(res_list_u0))
  str(res_list_u0)
  res_list_ma_u <- res_list_ma_u[!grepl("n", names(res_list_ma_u))]
  res_list_u <- c(res_list_u0, res_list_ma_u)
}

genes_list <- overlap_calc(res_list_u, sep = "&")
genes_list <- genes_list[sapply(genes_list, length) > 0]
genes_list <- genes_list[!grepl("&", names(genes_list))]
str(genes_list)
table_supp = reshape2::melt(genes_list); colnames(table_supp) <- c("Gene", "Group")
write.csv(table_supp, file = paste0(result_id, "table.csv"), row.names = FALSE)

if(is.null(order_celltype))
 order_celltype = grep("~|&", names(genes_list), value = TRUE, inver = TRUE)
tGroup <- if("0~1vs2~3~4~5" %in% names(res_list)){
  tvar = table(mdata[, c(mcols[1], "orig.asthma")])
  tvar = setNames(colnames(tvar)[apply(X = tvar, MARGIN = 1, FUN = which.max)], rownames(tvar))
  mdata$Disease = tvar[as.character(mdata[, mcols[1]])]; mcols <- unique(c(mcols, "Disease"))
  order_celltype[order_celltype %in% tGroup]
}else{
  tvar <- order_celltype[order_celltype %in% names(genes_list)]
  genes_list <- c(genes_list[tvar], genes_list[!names(genes_list) %in% tvar])
  order_celltype
}
mdata$tGroup = factor(as.character(mdata$Group), tGroup)
table(mdata$tGroup, useNA = 'always')
ssamples <- sample_even(mdata, cname = "Group", maxln = -100, v = T)
annor = reshape2::melt(genes_list)
tvar <- if(!"0~1vs2~3~4~5" %in% names(res_list)){
  factor(annor$L1, names(genes_list))
}else{ annor$L1 }
annor <- data.frame(Genes = tvar, row.names = as.character(annor$value))

source("/home/ciro/scripts/handy_functions/devel/plots.R")
fname <- paste0(result_id)
png(paste0(fname, "heatmap.png"), width = 1500, height = 1700, res = 250)
custom_heatmap(
  object = as.matrix(sc_cd4resting_0n1@assays$RNA@data)[, ssamples],
  annoc = mdata[ssamples, ],
  rnames = sub("'", "", unlist(genes_list)),
  orderby = "tGroup",
  use_mean = "tGroup",
  scale_row = TRUE,
  categorical_col = mcols[2],
  feature_order = FALSE,
  couls = v2cols(c(unique(unlist(mdata[, mcols])), levels(annor$Genes)), colours),
  hcouls = c('yellow', 'black', 'blue'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  gaps_row = unname(cumsum(sapply(genes_list, length))),
  annotation_row = annor#,
  # gaps_col = unname(cumsum(table(mdata[ssamples, mcols]))[order_celltype])
)
graphics.off()

### Heatmap sc - spec genes  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
dir.create("heatmaps_sc")

result_id = "heatmaps_sc/stim_disease"
myobject = "sc_cd4merged"
color_by_cols = c("Identity", "Disease", "Stim")
genes = c(
  "AMICA1", "AREG", "BATF", "CCL20", "CCL3", "CCL4", "CCR5", "CCR6", "CD28",
  "CD40LG", "CD55", "CD69", "CEBPB", "CREBBP", "CREM", "CSF2", "CTLA4", "CXCR6",
  "DDIT4", "DUSP1", "DUSP2", "DUSP4", "FABP5", "FASLG", "FKBP5", "FOS", "GPR65",
  "GZMA", "GZMB", "GZMH", "HOPX", "ICOS", "IFITM1", "IFNG", "IL10", "IL10RA",
  "IL13", "IL16", "IL17A", "IL2", "IL21", "IL21R", "IL23A", "IL23R", "IL26",
  "IL2RA", "IL32", "IL6R", "IL6ST", "ITGAE", "NFKB1", "NR4A2", "PRDM1", "PTGER4",
  "RGS1", "RHOB", "RHOG", "RHOH", "SAMSN1", "SOCS3", "SOD1", "STAT1", "STAT3",
  "STAT4", "STAT5A", "STAT5B", "TBX21", "TCF7", "TGFB1", "TNF", "TNFAIP3",
  "TNFRSF10B", "TNFRSF12A", "TNFRSF14", "TNFRSF18", "TNFSF10", "TNFSF12",
  "TNFSF14", "TNFSF8", "TNFSF9", "XCL1", "XCL2", "ZNF683"
)
sselect = NULL

result_id = "heatmaps_sc/disease_in_stim"
myobject = "sc_cd4stim"
color_by_cols = c("Disease")
genes = c(
  "ADAM17", "AHR", "ANXA1", "AREG", "BATF", "BCL10", "BCL2A1", "BCL3", "BCL6",
  "BCL7B", "BIRC2", "BIRC3", "BTG1", "BTG2", "BTG3", "CASP3", "CCL20", "CCL3",
  "CCL3L1", "CCL4", "CCL4L2", "CD28", "CD40LG", "CD44", "CD55", "CD69", "CD83",
  "CD97", "CEBPB", "CREBBP", "CREM", "CSF1", "CSF2", "CTLA4", "DDIT4", "DUSP1",
  "DUSP10", "DUSP14", "DUSP16", "DUSP18", "DUSP2", "DUSP4", "DUSP5", "DUSP6",
  "DUSP8", "EGR2", "EGR3", "EGR4", "FABP5", "FAM46C", "FAS", "FASLG", "FKBP4",
  "FOS", "FOSB", "FRMD4B", "GPR65", "GZMB", "HSPA1A", "HSPA1B", "ICOS", "ID2",
  "IER2", "IER3", "IER5", "IFNG", "IL10", "IL13", "IL17A", "IL2", "IL21", "IL21R",
  "IL23A", "IL23R", "IL26", "IL2RA", "IL6R", "IL6ST", "IRF4", "JUNB", "NFKB1",
  "NFKB2", "NR3C1", "NR4A1", "NR4A2", "NR4A3", "PLAUR", "PRDM1", "PTGER4", "RGS1",
  "RHOB", "RHOG", "RHOH", "SAMSN1", "SOCS3", "SOCS4", "SOD1", "STAT3", "STAT4",
  "STAT5A", "STAT5B", "TBX21", "TCF7", "TGFB1", "TNF", "TNFAIP3", "TNFRSF10B",
  "TNFRSF12A", "TNFRSF18", "TNFSF14", "TNFSF8", "TNFSF9", "TRAF1", "TRAF3",
  "TRAF4", "XCL1", "XCL2"
)
sselect = NULL

result_id = "heatmaps_sc/donor_in_stim"
myobject = "sc_cd4stim"
color_by_cols = c("orig.donor", "Disease")
genes = c(
  "AREG", "CCL20", "CCL3", "CCL4", "CEBPB", "CREBBP", "CREM", "CSF2", "CTLA4",
  "DDIT4", "DUSP1", "DUSP2", "DUSP4", "FABP5", "FASLG", "GZMB", "ICOS", "IFNG",
  "IL10", "IL13", "IL17A", "IL2", "IL21", "IL21R", "IL23A", "IL23R", "IL26",
  "IL2RA", "IL6R", "IL6ST", "NFKB1", "NR4A2", "PRDM1", "PTGER4", "RGS1", "RHOB",
  "RHOG", "RHOH", "SAMSN1", "SOCS3", "SOD1", "STAT3", "STAT4", "STAT5A", "STAT5B",
  "TBX21", "TCF7", "TGFB1", "TNF", "TNFAIP3", "TNFRSF10B", "TNFRSF12A", "TNFRSF18",
  "TNFSF14", "TNFSF8", "TNFSF9", "XCL1", "XCL2", "FKBP5", "GZMA", "GZMH", "HOPX",
  "IL10RA", "ITGAE", "STAT1", "TNFRSF14", "TNFSF10", "TNFSF12", "ZNF683"
)
sselect = NULL

eval(parse(text = paste0(
  myobject, "$Disease <- factor(", myobject, "$orig.asthma, c('MA', 'SA'))"
)))
eval(parse(text = paste0(
  myobject, "$Stim <- factor(stringr::str_to_title(gsub('_.*', '', ",
  myobject, "$orig.ident)), c('Resting', 'Stim'))"
)))
eval(parse(text = paste0(
  myobject, "$Identity <- ident_combine(", myobject, "@meta.data, c('Stim', 'Disease'))"
)))
eval(parse(text = paste0("table(", myobject, "[['Identity']])")))
eval(parse(text = paste0(
  "ssamples = filters_subset_df(sselect, ", myobject, "@meta.data, v = T)"
)))
eval(parse(text = paste0(
  "ssamples = sample_even(", myobject, "@meta.data, '", color_by_cols[1], "', -1000, v = T)"
)))
eval(parse(text = paste0("table(", myobject, "[['", color_by_cols[1], "']])")))
mygenes = eval(parse(text = paste0("show_found(genes, rownames(", myobject, "), v = T)")))
tvar <- eval(parse(text = paste0(
  "unique(unlist(lapply(", myobject, "@meta.data[ssamples, c('",
  paste0(color_by_cols, collapse = "', '"), "')], as.character)))"
)))
colours_i <- v2cols(unique(c(tvar)), colours)

for(i in list(FALSE, 'hclust')){
  fname <- paste0(result_id, ifelse(is.character(i), paste0("_", i), ""))
  png(paste0(fname, ".png"), width = 2000, height = 2500, res = 300)
  custom_heatmap(
    object = eval(parse(text = myobject)),
    annoc = eval(parse(text = paste0(myobject, "@meta.data[ssamples, ]"))),
    rnames = mygenes,
    use_mean = unique(c(color_by_cols[1], "orig.donor")),
    orderby = c(ifelse(grepl("donor", color_by_cols[1]), "Disease", color_by_cols[1]), "hc"),
    scale_row = TRUE,
    categorical_col = color_by_cols,
    feature_order = i,
    couls = colours_i,
    hcouls = c('yellow', 'black', 'blue'),
    topz = 2,
    verbose = TRUE,
    type = "pheat",
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_cols = FALSE, fontsize = 7
  )
  graphics.off()
}

### Heatmap sc - Special one with donors ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
genes = c(
  "GZMB", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "CCL20", "IL21", "IL17A",
  "IL17F", "TNF", "TNFSF14", "TNFSF12", "TNFSF10", "AREG", "IL2", "IL3", "CSF2",
  "IFNG", "IL4", "IL5", "IL13", "IL22", "IL10", "TGFB1", "IL6"
)
fconfigs = list(
  list(
    result_id = "heatmaps_sc_donor_specific",
    edata = "sc_cd4stim@assays$RNA@data",
    mdata = "sc_cd4stim@meta.data",
    annot_cols = c("orig.donor", "orig.asthma", "orig.Sex"),
    order_col = "orig.asthma", sselect = NULL
  ),
  list(
    result_id = "heatmaps_bulk_donor_specific",
    edata = "bulk_tpm",
    mdata = "bulk_metadata",
    annot_cols = c("Study_ID", "Disease", "Sex"),
    order_col = "Disease", sselect = list(c("Stim", "Stim"), c("Cell_type", "CD4-TRM"))
  )
)

for(fconfig in fconfigs){
  cat(fconfig$result_id, "\n")
  if(!dir.exists(fconfig$result_id) && grepl("\\/$", fconfig$result_id))
    dir.create(fconfig$result_id)

  ssamples = filters_subset_df(fconfig$sselect, eval(parse(text = fconfig$mdata)), v = TRUE)
  if(grepl("_sc_", fconfig$result_id)){
    mygenes = show_found(genes, rownames(eval(parse(text = fconfig$edata))), v = FALSE)
    edata_i = eval(parse(text = paste0("as.matrix(expm1(", fconfig$edata, "[mygenes, ]))")))
  }else{
    edata_i = eval(parse(text = fconfig$edata))
    rownames(edata_i) = features_parse_ensembl(rownames(edata_i))
  }
  mygenes = show_found(genes, rownames(edata_i), v = TRUE)
  mdata_i = eval(parse(text = paste0(fconfig$mdata, "[ssamples, ]")))

  genes_stats = stats_summary_table(
    mat = edata_i[mygenes, ssamples],
    groups = make_list(mdata_i, fconfig$annot_cols[1], grouping = TRUE),
    moments = c("mn", "p"),
    expr_cutoff = ifelse(grepl("_sc_", fconfig$result_id), 0, 10),
    verbose = TRUE
  )

  for(do_scale in c("", "_zscored")){
    for(i in c(1,3)){
      moments_i = moments[i]
      mat2plot <- genes_stats[, grepl(moments_i, colnames(genes_stats))]
      mat2plot = if(moments_i[[1]] == "_mean") log2(mat2plot + 1) else as.matrix(mat2plot)
      colnames(mat2plot) <- sub(moments_i, "", colnames(mat2plot))
      annoc = summarise_table(mdata_i[, fconfig$annot_cols], fconfig$annot_cols[1])[, -1]
      # data.table::setorderv(x = annoc, cols = fconfig$order_col)
      # new_order = rownames(annoc)
      new_order = unlist(lapply(levels(annoc[, fconfig$order_col]), function(i){
        hc = hclust(dist(t(mat2plot[, rownames(annoc)[annoc[, fconfig$order_col] == i]])))
        hc$labels[hc$order]
      }))
      annoc = annoc[new_order, ]
      print(lapply(annoc, table))
      mat2plot = mat2plot[, new_order]
      fname <- paste0(fconfig$result_id, "/values", moments_i, ".csv")
      write.csv(mat2plot, file = fname)

      if(do_scale != ""){
        mat2plot <- as.matrix(t(scale(t(mat2plot))))
        topz <- max(c(min(abs(c(range(mat2plot), 2))), 1))
        mat2plot[mat2plot > topz] <- topz; mat2plot[mat2plot < (-topz)] <- -topz;
      }
      hcolours <- grDevices::colorRampPalette(rev(c('yellow', 'black', 'blue')))(256)
      source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
      fname <- paste0(fconfig$result_id, "/values", moments_i, do_scale)
      cat(fname, "\n")
      pdf(paste0(fname, ".pdf"), height = 9, width = 8)
      void <- try(pheatmap(
        mat = mat2plot, color = hcolours, scale = "none", show_colnames = TRUE,
        annotation_col = annoc, cluster_cols = FALSE, border_color = NA,
        annotation_colors = lapply(annoc, v2cols, sour = colours),
        gaps_col = table(annoc[, fconfig$order_col])[1], angle_col = 90, annotation_names_col = FALSE
      ))
      if(class(void) == "try-error"){
        pheatmap(
          mat = mat2plot, color = hcolours, scale = "none", show_colnames = TRUE,
          annotation_col = annoc, cluster_cols = FALSE, border_color = NA, cluster_rows = FALSE,
          annotation_colors = lapply(annoc, v2cols, sour = colours),
          gaps_col = table(annoc[, fconfig$order_col])[1], angle_col = 90, annotation_names_col = FALSE
        )
      }
      dev.off()
    }
  }
}

### UMAP per disease and sex ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "disease_sex/resting_umap_"
oname = "resting"
result_id = "disease_sex/stim_umap_"
oname = "stim"
dir.create(dirname(result_id))
identity_column = c("orig.Sex", "orig.asthma")
eval(parse(text = paste0(
  'sc_cd4', oname, '$Identity <- do.call(paste, c(sc_cd4', oname, '@meta.data[, identity_column], sep = "_"))'
)))
eval(parse(text = paste0(
  'scells <- unname(sample_even(annot = sc_cd4', oname, '@meta.data, cname = "Identity", v = TRUE))'
)))
# for(i in unique(sc_cd4resting@meta.data[, "Identity"])){
  # eval(parse(text = paste0(
  #   'scells_i <- filters_subset_df(c("Identity", i), sc_cd4', oname, '@meta.data[scells, ], v = TRUE)'
  # )))
i <- "panels"
p <- LabelClusters(
  plot = DimPlot(
  object = eval(parse(text = paste0("sc_cd4", oname, "[, unlist(scells)]"))),
  cells = unlist(scells),
  reduction = "umap",
  group.by = eval(parse(text = paste0("sc_cd4", oname, "_clust"))),
  split.by = "Identity", ncol = 2,
  cols = if(oname == "resting") sc_cd4resting_ident$colours
),
  id = eval(parse(text = paste0("sc_cd4", oname, "_clust"))),
  repel = TRUE, parse = TRUE
) + theme(legend.position = "none")
# p <- p + facet_wrap(facets = ~Identity)

fname <- paste0(result_id, i)
pdf(paste0(fname, ".pdf"), 10, 10)
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"), 10, 10)
print(plot_blank(p))
dev.off()
# }

### Monocle 3 UMAP per disease ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "trajectory_plots/"
dir.create(result_id)
cds_f = "/home/ciro/ad_hoc/asthma_pjs/results/trajectory/resting_cd4_15p_sng_nomacro_x7n9_20PCs/object.rdata"
cds <- readfile(cds_f)

p <- plot_cells(
  cds = cds,
  reduction_method = "UMAP",
  color_cells_by = sc_cd4resting_clust,
  group_cells_by = 'cluster',
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE, # black circles
  label_roots = TRUE, # white circles
  label_leaves = TRUE, # gray circles
  graph_label_size = 4
) + theme_cowplot() + theme(legend.title = element_blank()) +
  scale_color_manual(values = v2cols(cds@colData[, sc_cd4resting_clust], sc_cd4resting_ident$colours))
pdf(paste0(result_id, "umap_clusters.pdf"), width = 8, height = 8)
print(p)
dev.off()
source("/home/ciro/scripts/handy_functions/devel/plots.R")
p2 <- plot_rm_layer(p, "text")
p2$layers <- p2$layers[1:2]
pdf(paste0(result_id, "umap_clusters_blank.pdf"), width = 8, height = 8)
print(plot_blank(p2))
dev.off()

### Crater single-cell ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/mnt/BioHome/ciro/scripts/handy_functions/devel/plots_crater.R')
source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
dir.create("crater_plots_sc/")
dir.create("crater_plots_sc_additional/")
fig_configs = list(
  trm_vs_teff = list(
    fnames = c(
      "Resting_0_VS_1" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters/0vs1/results_0vs1_mastlog2cpm.csv',
      "Resting_Male_VS_Female" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex/MalevsFemale/results_MalevsFemale_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = list(c(sc_cd4resting_clust, "0", "1"))
  ),
  trm_vs_teff = list(
    fnames = c(
      "Resting_0_VS_1" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters/0vs1/results_0vs1_mastlog2cpm.csv',
      "Resting_SA_VS_MA" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease/SAvsMA/results_SAvsMA_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = list(c(sc_cd4resting_clust, "0", "1"))
  ),
  trm_vs_teff = list(
    fnames = c(
      "Resting_Male_VS_Female" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex/MalevsFemale/results_MalevsFemale_mastlog2cpm.csv',
      "Resting_SA_VS_MA" = '/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease/SAvsMA/results_SAvsMA_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = list(c(sc_cd4resting_clust, "0", "1"))
  ),
  activation_vs_disease = list(
    fnames = c(
      "All_Stim_VS_Unstim" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/activation_cov_sex/STvsUS/results_STvsUS_mastlog2cpm.csv',
      "All_SA_VS_MA" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/disease_cov_sex/SAvsMA/results_SAvsMA_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = NULL
  ),
  activation_sex_mild = list(
    fnames = c(
      "MA_ST_VS_MA_US" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/mild/MA_STvsMA_US/results_MA_STvsMA_US_mastlog2cpm.csv',
      "MA_Male_VS_MA_Female" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/mild/MA_MalevsMA_Female/results_MA_MalevsMA_Female_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = NULL
  ),
  activation_sex_severe = list(
    fnames = c(
      "SA_ST_VS_SA_US" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/severe/SA_STvsSA_US/results_SA_STvsSA_US_mastlog2cpm.csv',
      "SA_Male_VS_SA_Female" = '/home/ciro/large/asthma_airways/results/scdgea/airways_cd4/comprs/severe/SA_MalevsSA_Female/results_SA_MalevsSA_Female_mastlog2cpm.csv'
    ), result_id = "crater_plots_sc/", selectss = NULL
  )
)

fig_configs = list(
  female_disease0_vs_disease1 = list(
    fnames = c(
      "C0_Female_SAvs0_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/0_Female_SAvs0_Female_MA/results_0_Female_SAvs0_Female_MA_mastlog2cpm.csv",
      "C1_Female_SAvs1_Female_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/1_Female_SAvs1_Female_MA/results_1_Female_SAvs1_Female_MA_mastlog2cpm.csv"
    ), columns = c("orig.asthma", "orig.Sex", sc_cd4resting_clust, "combn"), result_id = "crater_plots_sc_additional/",
    selectss = list(c(sc_cd4resting_clust, "0", "1"), c("orig.Sex", "Female")),
    highlight_genes = c(
      "GZMB", "GZMH", "PRF1", "ZNF683", "TRAF3IP3", "IFITM1", "IFNG", "FKBP5",
      "CXCR3", "TNFSF14", "CREM", "TNFAIP3", "CXCR4", "METRNL", "DUSP1", "DUSP4")
  ),
  male_disease0_vs_disease1 = list(
    fnames = c(
      "C0_Male_SAvs0_Male_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/0_Male_SAvs0_Male_MA/results_0_Male_SAvs0_Male_MA_mastlog2cpm.csv",
      "C1_Male_SAvs1_Male_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex_disease/1_Male_SAvs1_Male_MA/results_1_Male_SAvs1_Male_MA_mastlog2cpm.csv"
    ), columns = c("orig.asthma", "orig.Sex", sc_cd4resting_clust, "combn"), result_id = "crater_plots_sc_additional/",
    selectss = list(c(sc_cd4resting_clust, "0", "1"), c("orig.Sex", "Male")),
    highlight_genes = c(
      "GZMB", "GZMH", "PRF1", "ZNF683", "TRAF3IP3", "IFITM1", "IFNG", "FKBP5",
      "CXCR3", "TNFSF14", "CREM", "TNFAIP3", "CXCR4", "METRNL", "DUSP1", "DUSP4")
  )
)

fig_configs = list(
  list(
    fnames = c(
      "C0_MalevsC0_Female" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex/0_Malevs0_Female/results_0_Malevs0_Female_mastlog2cpm.csv",
      "C0_SAvsC0_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease_clusters/SA_0vsMA_0/results_SA_0vsMA_0_mastlog2cpm.csv"
    ), columns = c("orig.asthma", "combn"), result_id = "crater_plots_sc_additional/",
    selectss = list(c(sc_cd4resting_clust, "0")),
    highlight_genes = c(
      "GZMB", "GZMH", "PRF1", "ZNF683", "TRAF3IP3", "IFITM1", "IFNG", "FKBP5",
      "CXCR3", "TNFSF14", "CREM", "TNFAIP3", "CXCR4", "METRNL", "DUSP1", "DUSP4")
  ),
  list(
    fnames = c(
      "C1_MalevsC1_Female" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_sex/1_Malevs1_Female/results_1_Malevs1_Female_mastlog2cpm.csv",
      "C1_SAvsC1_MA" = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/disease_clusters/SA_1vsMA_1/results_SA_1vsMA_1_mastlog2cpm.csv"
    ), columns = c("orig.asthma", "combn"), result_id = "crater_plots_sc_additional/",
    selectss = list(c(sc_cd4resting_clust, "0")),
    highlight_genes = c(
      "GZMB", "GZMH", "PRF1", "ZNF683", "TRAF3IP3", "IFITM1", "IFNG", "FKBP5",
      "CXCR3", "TNFSF14", "CREM", "TNFAIP3", "CXCR4", "METRNL", "DUSP1", "DUSP4")
  )
)

sc_cd4resting@meta.data$combn = ident_combine(sc_cd4resting@meta.data, rev(c("orig.asthma", "orig.Sex", sc_cd4resting_clust)), "_")
table(sc_cd4resting@meta.data$combn)
degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))
degfilt = list(mean = list("<0", NA), min_padj = list(">0.01", 1))
sig_thresh_fc = c(0.25, 0.5)

edata = as.matrix(expm1(sc_cd4resting@assays$RNA@data))
mygenes = grep("XIST|RPS4Y1", rownames(edata), value = TRUE, invert = TRUE)
for(fig_config in fig_configs){
  for(fc in as.character(sig_thresh_fc)){
    void <- crater_plot(
      tests_list = fig_config$fnames,
      edataf = edata,
      annotf = sc_cd4resting@meta.data,
      sample_filter = fig_config$selectss,
      gene_filter = degfilt,
      feature_subset = mygenes,
      topgenes = fig_config$highlight_genes,
      lfcthresh = fc,
      column4stats = fig_config$columns,
      outputname = fig_config$result_id,
      plot_interactive = TRUE,
      plot_squared = TRUE,
      v = TRUE
    )
  }
}

fnames = c(
  "clusters_female_disease/heatmap_sc_0_Female_SAvs0_Female_MA-1_Female_SAvs1_Female_MA.csv",
  "crater_plots_sc/FC0.25_c0_female_savs0_female_ma_against_c1_female_savs1_female_ma_stats.csv"
)
df_list = lapply(fnames, readfile, row.names = 1, check.names = FALSE)
colnames(df_list[[1]])[1] = "Group"
colnames(df_list[[2]])[2:3] = paste0(colnames(df_list[[2]])[2:3], ".fc")
colnames(df_list[[2]]) = gsub("^mean$", "Global_mean", colnames(df_list[[2]]))
df_list[[2]]$log2_mean = NULL
df_list[[2]]$filters = NULL
df_list[[2]]$lfcthresh = NULL
str(df_list)
dfm = joindf(df_list[[2]], df_list[[1]])[, unlist(sapply(df_list, colnames))]
str(dfm)
write.csv(dfm, file = gsub(".csv$", "_stats.csv", fnames[1]), row.names = FALSE)
system(paste("head", gsub(".csv$", "_stats.csv", fnames[1])))

library(openxlsx)
source("/home/ciro/scripts/handy_functions/devel/supp_table.R")
dgeaheaders = list(
  none = c(gene_name = "TEXT"),
  "MAST - Likelihood ratio test" = c(Group = "TEXT", "\\.fc$" = "NUMBER", "padj" = "SCIENTIFIC", "signif" = "NUMBER"),
  "Mean expression (Seurat Normalised)" = c("mean" = "NUMBER"),
  "Fraction of expressing cells (>0)" = c("_percen" = "NUMBER")
)
supptable <- supp_table(
  mytables = dfm,
  headers = dgeaheaders,
  title_name = "Differential gene expression analysis. Severe vs Mild in Female subjects.",
  filename = gsub(".csv$", "_stats", fnames[1]),
  verbose = 1
)

### Crater bulk ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('/mnt/BioHome/ciro/scripts/handy_functions/devel/plots_crater.R')
source("/home/ciro/scripts/handy_functions/devel/plots.R")
## Input
fig_configs = list(
  trmdp_activation_vs_disease = list(
    fnames = c(
      TRMDP_SA_VS_MA = '/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/disease/CD4-TRM_SAvsCD4-TRM_MA/results_CD4-TRM_SAvsCD4-TRM_MA_deseq2.csv',
      TRMDP_STIM_VS_RESTING = '/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/activation/CD4-TRM_StimvsCD4-TRM_Unstim/results_CD4-TRM_StimvsCD4-TRM_Unstim_deseq2.csv'
    ),
    selectss = list(c("Cell_type", "CD4-TRM"))
  ),
  trmsp_activation_vs_disease = list(
    fnames = c(
      TRMSP_SA_VS_MA = '/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/disease/CD4-TRM-like_SAvsCD4-TRM-like_MA/results_CD4-TRM-like_SAvsCD4-TRM-like_MA_deseq2.csv',
      TRMSP_STIM_VS_RESTING = '/home/ciro/large/asthma_airways/results/deseq2/airways_platex/comprs/activation/CD4-TRM-like_StimvsCD4-TRM-like_Unstim/results_CD4-TRM-like_StimvsCD4-TRM-like_Unstim_deseq2.csv'
    ),
    selectss = list(c("Cell_type", "CD4-TRM-like"))
  )
)
result_id = "crater_plots_bulk/"
degfilt = list(mean = list("<10", NA), min_padj = list(">0.05", 1))
fannog = "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv"

## Reading
annog <- readfile(fannog, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
annog$ensembl_id <- rownames(annog)
rownames(annog) <- paste0(annog$ensembl_id, "_", annog$gene_name)

## Results
dir.create(result_id)
for(fig_config in fig_configs){
  void <- crater_plot(
    tests_list = fig_config$fnames,
    edataf = bulk_tpm,
    annotf = bulk_metadata,
    sample_filter = fig_config$selectss,
    gene_filter = degfilt,
    feature_subset = rownames(annog[annog$gene_type_4 == "protein_coding", ]),
    topgenes = 10,
    lfcthresh = "0.25",
    outputname = paste0(result_id, "prot_coding_"),
    plot_interactive = TRUE,
    # plot_squared = TRUE,
    v = TRUE
  )
}

### Scatter: co-expression ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
myobject = "sc_cd4stim"
cut_offs = c(0, log(1+1))
cut_offs = log(2+1)
cut_offs = log(10+1)
cut_offs = log(0+1)
## ---------------------------------------------------
yaxes = c("GZMB", "HOMBRINK_TEICHMANN_2.Score")
result_id = "coexpression/"; dir.create("coexpression")
result_id = "coexpression/disease_"
genes = c(
  "IL21", "IL6", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "TNF", "TNFSF14",
  "TNFSF12", "TNFSF10", "AREG", "IL2", "IL3", "CSF2", "CCL20", "IFNG", "IL4",
  "IL5", "IL13", "IL17A", "IL17F", "IL22", "IL10", "TGFB1"
)

## ---------------------------------------------------
yaxes = c("IFNG", "IL17A")
result_id = "coexpression_specific/disease_"; dir.create("coexpression_specific")
genes = c("IL17A", "IL13", "IL13")

## ---------------------------------------------------
yaxes = c("GZMB")
result_id = "coexpression_donor/"; dir.create("coexpression_donor")
genes = c(
  "IL4", "IL5", "IL13", "IL17A", "IL21", "TNF", "CCL3", "CCL4", "CCL5", "CCL20", "CSF2"
)

mygenes = show_found(genes, rownames(eval(parse(text = myobject))), v = TRUE)
pdata <- FetchData(
  object = eval(parse(text = myobject)),
  vars = unique(c("cluster", "orig.donor", "orig.asthma", yaxes, mygenes))
)
df_signatures <- read.csv("signatures_sc_stim_ayearoftrmness_modulescore/signatures.csv", row.names = 1)
pdata <- joindf(pdata, df_signatures)
facet_var = if(grepl("disease", result_id)){
  # pdata <- pdata[sample_even(pdata, "orig.asthma", v = TRUE), ]
  "orig.asthma"
}else if(grepl("donor", result_id)){
  pdata <- pdata[table(pdata$orig.donor)[pdata$orig.donor] > 50, ]
  "orig.donor"
}else{ "orig.ident" }
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
axes <- data.frame(
  rep(mygenes, length(yaxes)),
  rep(yaxes, each = length(mygenes)), stringsAsFactors = FALSE
)
tvar <- apply(axes, 1, function(x) paste(unique(sort(x)), collapse = "_") )
axes <- axes[!duplicated(tvar) & stringr::str_count(tvar, "_") > 0, ]
groups <- names(table(pdata[, facet_var]))
mytab = suppressMessages(reshape2::melt(lapply(
  X = 1:nrow(axes),
  FUN = function(i){
  pdata_i = pdata[, c(axes[i, 1], axes[i, 2])]
  colnames(pdata_i) <- c("x", "y")
  lapply(
    X = setNames(nm = groups),
    FUN = function(x){
      y <- stat_quadrant_fun(
        pdata_i[pdata[, facet_var] == x, ],
        xintercept = cut_offs[1], yintercept = cut_offs[1],
        type = "percent", loc_type = "component"
      )
      y[, 1] <- paste0(axes[i, 1], y[, 1])
      y[, 2] <- paste0(axes[i, 2], y[, 2]); y
    }
  )
})))
str(mytab)
fname <- paste0(result_id, "table_threshold", expm1(cut_offs[1]), ".csv")
# head(reshape2::dcast(mytab[, -ncol(mytab)], x + y ~ L2, value.var = "count"))
write.csv(mytab[, -ncol(mytab)], file = fname, row.names = FALSE)
max_lev = max(pdata[, mygenes])
for(j in cut_offs){
  cat("Threshold:", j, "\n")
  for(i in 1:nrow(axes)){
    cat(axes[i, 1], "vs", axes[i, 2], "\n")
    possitive = rowSums(pdata[, c(axes[i, 1], axes[i, 2])] > 0) > 1; pdata$Density <- 0
    pdata[possitive, ]$Density <- MASS_kde2d(
      x = pdata[possitive, axes[i, 1]],
      y = pdata[possitive, axes[i, 2]]
    )
    p <- ggplot(
      data = pdata,
      mapping = aes_string(x = axes[i, 1], y = axes[i, 2], color = "Density")
    ) + geom_point(size = 0.5) +
      geom_density2d(
        data = pdata[possitive, ], colour = "#c0c5ce"
      ) + viridis::scale_color_viridis(option = "magma")
    if(grepl("disease", result_id)) p <- p + facet_wrap(~orig.asthma)
    if(grepl("donor", result_id)) p <- p + facet_wrap(~orig.donor) + xlim(-0.5, max_lev)
    p <- plots_add_quadrants(
      p, limits = list(j, ifelse(j > 0 && grepl("HOMBRINK", axes[i, 2]), 0.1, j)),
      type = "percent", label.size = 2)
    fname <- paste0(result_id, axes[i, 2], "_", axes[i, 1], "_threshold", expm1(j))
    widths = if(grepl("disease", result_id)) 12 else if(grepl("donor", result_id)) 18 else 7
    heights = if(grepl("donor", result_id)) 16 else 7
    pdf(paste0(fname, ".pdf"), width = widths, height = heights)
    print(p)
    graphics.off()
    pdf(paste0(fname, "_blank.pdf"), width = widths, height = heights)
    print(plot_blank(p, "ext|ine"))
    graphics.off()
  }
}

### Scatter: markers on UMAP ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/plots.R")
dir.create("markers")
result_id = "markers/"
result_id = "markers/disease_"
disease_params = function(x){
  x + facet_wrap(~orig.asthma) +
    theme(
      strip.background = element_rect(fill = NA),
      strip.text = element_text(face = "bold")
    )
}
genes = c(
  "IL4", "IL5", "IL13", "IL17A", "IL17F", "IL22", "TNF", "IFNG", "CCL3", "CCL4",
  "CCL5", "CCL20", "IL21", "IL10", "TNFSF14", "AREG", "CSF2", "TGFB1", "GZMB"
)
myobject = "sc_cd4stim"
scells = sample_even(eval(parse(text = myobject))[[]], "orig.asthma", v = TRUE)

mygenes = show_found(genes, y = rownames(eval(parse(text = myobject))), verbose = TRUE)
ddf = FetchData(eval(parse(text = myobject)), vars = c(redu[[1]], mygenes, colnames(eval(parse(text = myobject))[[]])))
cols = couls_opt$red_gradient[[3]]
for(g in mygenes){
  cat(g, "\n")
  p <- ggplot(
    data = ddf[scells, ],
    mapping = aes_string(x = redu[[1]][1], y = redu[[1]][2], color = g)
  ) + geom_point(size = 0.2) + theme_cowplot() +
    labs(color = "Seurat\nNormalized") +
    scale_colour_gradientn(colours = cols)
  if(grepl("disease_$", result_id)) p <- disease_params(p)

  fname <- paste0(result_id, myobject, "_umap_", g)
  pdf(paste0(fname, ".pdf"), width = ifelse(grepl("disease", result_id), 14, 7))
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = ifelse(grepl("disease", result_id), 14, 7))
  print(plot_blank(p))
  dev.off()
}

### Module score ### -----------------------------------------------------------
eval(parse(text = paste0(myobject, "= Seurat::AddModuleScore(",
  "object = ", myobject, ",",
  "features = list(custom_signature = mygenes))"
)))
ddf = FetchData(eval(parse(text = myobject)), vars = c(redu[[1]], colnames(eval(parse(text = myobject))[[]])))
p <- ggplot(
  data = ddf[scells, ],
  mapping = aes_string(x = redu[[1]][1], y = redu[[1]][2], color = "Cluster1")
) + geom_point(size = 0.2) + theme_cowplot() +
  labs(color = "Module\nScore") +
  scale_colour_gradientn(colours = cols)
if(grepl("disease_$", result_id)) p <- disease_params(p)

fname <- paste0(result_id, myobject, "_umap_module_score")
pdf(paste0(fname, ".pdf"), width = ifelse(grepl("disease", result_id), 14, 7))
print(p)
dev.off()

### Stim vs Unstim heatmap ### -------------------------------------------------
if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)

do_raster = TRUE
sc_cd4merged <- ScaleData(object = sc_cd4merged, features = mygenes)
sc_cd4merged$Identity <- ident_combine(sc_cd4merged@meta.data, c("orig.stim", "orig.Sex", "orig.asthma"))
table(sc_cd4merged$Identity)

p <- DoHeatmap(
  object = sc_cd4merged,
  features = genes,
  cells = sample_even(sc_cd4merged@meta.data, "orig.stim", maxln = -3000),
  group.by = "Identity",
  raster = do_raster
)

fname <- paste0(result_id, "merged_heatmap", ifelse(do_raster, "", "_noraster"))
pdf(paste0(fname, ".pdf"))
print(p)
dev.off()

### Merging some comparisons ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setwd("/home/ciro/ad_hoc/asthma_airways/results/scdgea/airways_cd4/comprs")
fnames0 = c(
  in_all_stim_vs_unstim = "activation/STvsUS/results_STvsUS_mastlog2cpm.csv",
  in_ma_stim_vs_unstim = "stim_disease/ST_MAvsUS_MA/results_ST_MAvsUS_MA_mastlog2cpm.csv",
  in_sa_stim_vs_unstim = "stim_disease/ST_SAvsUS_SA/results_ST_SAvsUS_SA_mastlog2cpm.csv",
  in_stim_sa_vs_ma = "stim_disease/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv"
  # in_all_sa_vs_ma = "disease/SAvsMA/results_SAvsMA_mastlog2cpm.csv"
)
fnames <- unlist(lapply(names(fnames0), function(x){
  y <- fnames0[[x]]
  z <- sub("\\/.*", "_cov_sex/", y); w <- paste0(basename(dirname(y)), "/", basename(y))
  setNames(c(y, paste0(z, w)), c(x, paste0(x, "_sexcov")))
}))
file.exists(fnames)
res_list = lapply(X = fnames, FUN = readfile, row.names = 1, stringsAsFactors = FALSE)
str(res_list[1:2])
res_list_stats = lapply(X = names(res_list), FUN = function(x){
  y <- res_list[[x]][, grepl("log2F|padj|group|percentage", colnames(res_list[[x]]))]
  colnames(y) <- paste0(x, "_", colnames(y)); y
})
str(res_list_stats[1:3])
tvar <- unique(unlist(lapply(res_list, rownames)))
res_summ = data.frame(gene_name = paste0("'", tvar), row.names = tvar)
for(i in 1:length(res_list_stats)){
  res_summ <- joindf(res_summ, res_list_stats[[i]])
}
dim(res_summ)
headmat(res_summ)
write.csv(res_summ, file = "stim_disease_comparisons.csv", row.names = FALSE)

### Checking MAST LFC calculation ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res_file="/home/ciro/ad_hoc/asthma_airways/results/scdgea/airways_cd4/comprs/stim_disease/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv"
res_file="/home/ciro/ad_hoc/asthma_airways/results/scdgea/stim_cd4_15p_sng_nomacro9n10/comprs/stim_disease_cov_sex/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv"
res_file="/home/ciro/ad_hoc/asthma_airways/results/scdgea/airways_cd4/comprs/stim_disease/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv"
res_file="/home/ciro/ad_hoc/asthma_airways/results/scdgea/airways_cd4/comprs/stim_disease_cov_sex/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv"
res_file="/home/ciro/ad_hoc/asthma_airways/results/scdgea/airways_cd4/comprs/disease/SAvsMA/results_SAvsMA_mastlog2cpm.csv"
grep -E "padj|IL21|CCL3" ${res_file} | cut -d, -f3,6,8,9
# https://rdrr.io/bioc/MAST/man/logFC.html
# u(contrast1)v(contrast1)-u(contrast0)v(contrast0)
# u(x)= crossprod(x, coefC)
# v(x)= 1/(1+exp(-crossprod(coefD, x)))
# MAST::logFC(summaryCond, contrast0 = "SA")

load("~/ad_hoc/asthma_airways/results/scdgea/stim_cd4_15p_sng_nomacro9n10/comprs/disease/SAvsMA/a1_rdata/summ_SAvsMA.RData")
load("~/ad_hoc/asthma_airways/results/scdgea/stim_cd4_15p_sng_nomacro9n10/comprs/stim_disease_cov_sex/ST_SAvsST_MA/a1_rdata/summ_ST_SAvsST_MA.RData")
res_summ = summaryCond[[1]]
head(res_summ)
res_summ[res_summ$primerid=="IL21", ]

load("~/ad_hoc/asthma_airways/results/scdgea/stim_cd4_15p_sng_nomacro9n10/comprs/disease/SAvsMA/a1_rdata/regression_SAvsMA.RData")
zlm_res
x = "IL21"; cont = "conditionSA"
some_genes = c("IL21", "IFNG", "IL5")
logFC(zlm_res[some_genes, ])
u_fun  = function(x, cont){
  crossprod(zlm_res@vcovC[cont, , x], zlm_res@coefC[x, ])
}
v_fun  = function(x, cont){
  1/(1+exp(-crossprod(zlm_res@coefD[x, ], zlm_res@vcovD[, cont, x])))
}
(u_fun("IL21", "(Intercept)") * v_fun("IL21", "(Intercept)")) - (u_fun("IL21", "conditionSA") * v_fun("IL21", "conditionSA"))
zz <- MAST::zlm(~condition+cngeneson, zlm_res@sca[some_genes,])
MAST::logFC(zz)
MAST::getLogFC(zz)


### Donors summary ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/clustering/R/subject_summary.R")
dir.create("donor_summary")
fconfigs = list(
  list(
    mdata = c(
      Unfilt = "/home/ciro/large/asthma_pjs/results/demuxlet/sever_asthma/resting_cd4.rdata",
      resting_cd4_15p = "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p/clustering/zetInfo/metadata_19PCs_30Ks_0.06667JD.RData",
      Filt = "/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/metadata_20PCs_30Ks_0.06667JD.RData"
    ), cluster = sc_cd4resting_clust
  ),
  list(
    mdata = c(
      Unfilt = "/home/ciro/large/asthma_pjs/results/demuxlet/sever_asthma/stim_cd4.rdata",
      xCD8 = "/home/ciro/large/asthma_pjs/results/clust_seurat/stim_cd4_15p/clustering/zetInfo/metadata_17PCs_30Ks_0.06667JD.RData",
      Filt = "/home/ciro/large/asthma_pjs/results/clust_seurat/stim_cd4_15p_sng_nomacro9n10/clustering/zetInfo/metadata_17PCs_30Ks_0.06667JD.RData"
    ), cluster = sc_cd4stim_clust
  )
)

for(fconfig in fconfigs){
  ddf <- subject_summary(
    metadata_list = fconfig$mdata,
    subject_id = "orig.donor",
    append_cols = list(tables = 1, names = c("orig.asthma", "origlib", "orig.stim", "orig.cell_type", "orig.Sex")),
    append_num = list(tables = seq(fconfig$mdata), names = c(fconfig$cluster, "orig.class"))
  )

  colnames(ddf) <- gsub("orig|orig\\.", "", colnames(ddf))
  tvar <- paste0(paste0(names(fconfig$mdata)[-length(fconfig$mdata)], "\\.[0-9]{1,}$"), collapse = "|")
  ddf <- data.frame(ddf)[ddf$subject_name1 != "NA", !grepl(tvar, colnames(ddf))]
  str(ddf)
  fname <- paste0("donor_summary/", basename(gsub("clustering/ze.*", "", fconfig$mdata[length(fconfig$mdata)])), ".csv")
  cat(fname, "\n")
  write.csv(ddf, file = fname, quote = FALSE, row.names = FALSE)
}
# tvar = "/home/ciro/large/asthma_pjs/results/clust_seurat/stim_cd3_degs/filtering/cd4_stimfromovCD4ST_aggr_mapped_and_stimfrom006_cells_noCD4p.csv"
# kk = readfile("/home/ciro/large/asthma_pjs/results/clust_seurat/stim_cd3_degs/filtering/cd4_stimfromovCD4ST_aggr_mapped_and_stimfrom006_cells_noCD8p.rdata")
# stim_cells = readfile(tvar, stringsAsFactors = FALSE)[, 2]
# table(kk$origlib[rownames(kk) %in% stim_cells])

# less /home/ciro/ad_hoc/hayley/raw/sever_asthma3/COUNTS_hg19/006_10x_018_GS_SH_3S_18/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l
# CLUSTRES="stim_cd4_15p/zetInfo/logs/seuratExploration_stim_cd4_15p_FPTM_17PCs.log
# stim_cd4_15p_sng_nomacro9n10/zetInfo/logs/seuratExploration_stim_cd4_15p_sng_nomacro9n10_FPTM_17PCs.log
# resting_cd4_15p/zetInfo/logs/seuratExploration_resting_cd4_15p_FPT_19PCs.log
# resting_cd4_15p_sng_nomacro_x7n9/zetInfo/logs/seuratExploration_resting_cd4_15p_sng_nomacro_x7n9_FPT_20PCs.log"
# for CLUSTRE in ${CLUSTRES[@]}; do
#   echo ${CLUSTRE/\/*/}
#   grep -E ".filt.cells =|subs.names =" /home/ciro/ad_hoc/asthma_airways/results/clust_seurat/${CLUSTRE}
# done
# grep -E "_FILT_CELLS" /home/ciro/asthma_airways/scripts/jobConstSeu_setting_*_cd4.sh

### QC metrics ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "qc_plots/"; dir.create(result_id)
fconfigs = list(
  list(
    object = "sc_cd4resting",
    elbow = 20
  ),
  list(
    object = "sc_cd4stim",
    elbow = 17
  )
)
qc_vars = grep("_RNA|percent.mt", colnames(sc_cd4resting[[]]), value = TRUE)
for(fconfig in fconfigs){
  cat("- Elbow\n")
  pelbow <- ElbowPlot(
    object = eval(parse(text = fconfig$object)),
    ndims = length(eval(parse(text = fconfig$object))@reductions$pca@stdev)
  ) + geom_vline(xintercept = fconfig$elbow, linetype = "dotted", color = "gray") +
    geom_point(aes(x = dims, y = stdev), size = 3)
  fname <- paste0(result_id, fconfig$object, "_", fconfig$elbow, "pcs")
  pdf(paste0(fname, '_elbow.pdf'), width = 7, height = 5); print(pelbow); dev.off()
  pdf(paste0(fname, '_elbow_blank.pdf'), width = 7, height = 5)
  print(plot_blank(pelbow, "ext|ine")); dev.off()

  cat("- HVGs\n")
  top_n <- head(x = VariableFeatures(object = eval(parse(text = fconfig$object))), 30)
  p1 <- VariableFeaturePlot(eval(parse(text = fconfig$object)))
  phvg <- LabelPoints(plot = p1, points = top_n, repel = TRUE, size = 5) +
    theme(legend.position = c(0, 75), legend.direction = "vertical")
  pdf(paste0(fname, '_hvg.pdf'), width = 7, height = 5.5); print(phvg); dev.off()
  pdf(paste0(fname, '_hvg_blank.pdf'), width = 7, height = 5.5)
  print(plot_blank(phvg)); dev.off()

  ddf = FetchData(
    object = eval(parse(text = fconfig$object)),
    vars = c(qc_vars, eval(parse(text = paste0(fconfig$object, "_clust"))))
  )
  couls_i = eval(parse(text = paste0(fconfig$object, "_ident$colours")))
  if(is.null(couls_i)){
    couls_i <- gg_color_hue(levels(ddf[, eval(parse(text = paste0(fconfig$object, "_clust")))]))
  }
  cat("- Violins\n")
  for(qc_i in qc_vars){
    cat(" *", qc_i, "\n")
    # ddf$tmp <- if(qc_i == "percent.mt") ddf[, qc_i] else log2(ddf[, qc_i] + 1)
    pviolin <- violin(
      ddf,
      eval(parse(text = paste0(fconfig$object, "_clust"))),
      qc_i
    ) + theme(legend.position = "none") +
      scale_fill_manual(values = couls_i) +
      scale_color_manual(values = couls_i)
    pdf(paste0(fname, '_violin_', qc_i, '.pdf'), width = 7, height = 5)
    print(pviolin); dev.off()
    pdf(paste0(fname, '_violin_', qc_i, '_blank.pdf'), width = 7, height = 5)
    print(plot_blank(pviolin)); dev.off()
  }
}

### Specific GSEAs on resting single-cell ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cp /Users/ciro/Documents/liai/asthma/airways/instructions/GSEA_supp_table_8_22JUNE21.csv /Volumes/ciro/asthma_airways/info/
# clean_csv /home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
comparisons_df = data.frame(
  group1 = c(levels(sc_cd4resting@meta.data[, sc_cd4resting_clust]), "0"),
  group2 = c(rep("REST", nlevels(sc_cd4resting@meta.data[, sc_cd4resting_clust])), "1"),
  column = sc_cd4resting_clust, stringsAsFactors = FALSE
)
fconfigs = list(
  list(
    result_id = "gsea_spec_sc_cd4resting/",
    object = "sc_cd4resting",
    file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
    comparisons = comparisons_df
  )
)

for (fconfig in fconfigs) {
  globalsign = list()
  signatures_files <- fconfig$file[file.exists(fconfig$file)]
  tvar <- lapply(signatures_files, readfile, stringsAsFactors = FALSE)
  tvar <- unlist(lapply(tvar, as.list), recursive = FALSE)
  tvar <- sapply(tvar, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
  tvar <- tvar[sapply(tvar, length) > 0]; str(tvar)
  signatures_list <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
  str(signatures_list)
  tmp <- unname(gsub(" .*", "", sapply(signatures_list, head, 1)))
  tvar <- paste0(gsub("\\.{2,}", ".", gsub("signature|genes|Bulk.", "", names(signatures_list))), tmp)
  names(signatures_list) <- tvar
  signatures_subset <- lapply(signatures_list, function(x){
    x[x %in% rownames(eval(parse(text = fconfig$object)))]
  })

  ddf <- vlist2df_diff(
    x = signatures_subset,
    y = signatures_list,
    delim = "Not found or < 2 pct"
  )
  str(headtail(ddf))
  tvar <- which(ddf == "Not found or < 2 pct", arr.ind = TRUE)
  ddf[unique(tvar[,1]), unique(tvar[,2])]
  if(!dir.exists(fconfig$result_id) && grepl("\\/$", fconfig$result_id))
    dir.create(fconfig$result_id)
  write.csv(ddf, file = paste0(fconfig$result_id, "features_lists.csv"))

  gsea_results <- gsea_matrix(
    mat = expm1(eval(parse(text = paste0(fconfig$object, "@assays$RNA@data")))),
    groups = fconfig$comparisons,
    metadata = eval(parse(text = paste0(fconfig$object, "@meta.data"))),
    metric = "Signal2Noise",
    gsea_list = signatures_list,
    method = "fgsea",
    path = fconfig$result_id,
    plot_it = TRUE,
    verbose = TRUE
  )
  for(i in c(0.05, 0.2)){
    for(j in c(1.5, 1)){
      x <- gsea_plot_summary(
        tests_list = gsea_results,
        path = fconfig$result_id,
        padjthr = i,
        nesthr = j
      )
    }
  }
}

### Per cluster radar plot ### -------
source("/home/ciro/scripts/handy_functions/devel/plots.R")
result_id = "gsea_spec_sc_cd4resting/specific/"
if(!dir.exists(result_id) && grepl("\\/$", result_id)) dir.create(result_id)
gsea_results = readRDS("gsea_spec_sc_cd4resting/tests_RNA_snn_res.0.4_0_18lists.rds")

gsea_summary_list <- gsea_summary(path = dirname(result_id))
mysum_radar <- data.frame(gsea_summary_list[[1]], check.names = FALSE)
mysum_radar <- mysum_radar[, colnames(mysum_radar) != "0vs1"]

ddf = reshape2::melt(cbind(mysum_radar, rownames(mysum_radar)))
colnames(ddf) <- c("Pathway", "Cluster", "NES")
ddf$NES[ddf$NES < 0] <- 0
# ddf$Pathway <- stringr::str_wrap(ddf$Pathway, width = 10)
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto(
    "CordRadar", CoordPolar, theta = theta, r = r, start = start,
    direction = sign(direction),
    is_linear = function(coord) TRUE)
}
p <- ggplot(ddf, aes(x = Cluster, y = NES, color = Pathway, group = Pathway)) +
  geom_line(color = "#BEBEBE") + geom_point(size = 4, ) +
  ylim(c(-0.5, max(ddf$NES))) + coord_radar() +
  facet_wrap(~ Pathway) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    panel.grid.major.x = element_line(size = 0.1, colour = "#ebebeb"),
    panel.grid.major.y = element_line(size = 0.1, colour = "#ebebeb"),
    panel.border = element_rect(fill = NA)
  )
pdf(paste0(result_id, "_radar_facet.pdf"), width = 14, height = 14)
print(p)
dev.off()
pdf(paste0(result_id, "_radar_facet_blank.pdf"), width = 14, height = 14)
print(plot_blank(p))
dev.off()

mysum_radar[mysum_radar < 0] <- NA
mysum_radar[is.na(mysum_radar)] <- 0
max_point_all = max(mysum_radar)
for (i in rownames(mysum_radar)) {
  cat(i, "\n")
  p1 <- make_radar(
    mysum_radar[i, ], mid_point = 1, max_point = max_point_all
  ) + theme(legend.position = "none")
  pdf(paste0(result_id, "radar_", i, ".pdf"), width = 10, height = 10)
  print(p1)
  dev.off()
  pdf(paste0(result_id, "radar_", i, "_blank.pdf"), width = 10, height = 10)
  print(plot_blank(p1))
  dev.off()
}
