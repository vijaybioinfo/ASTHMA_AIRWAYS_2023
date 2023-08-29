#!/usr/bin/R

###############
# Global data #
###############

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-02-22
# ---

# This script loads all the files used in a global way

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }
cat(red_bold("------------------ Setting global parameters -----------------\n"))

### Objects names ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat(cyan("------------------ Output location, objects, and file names\n"))
fig_dir = '/home/ciro/ad_hoc/asthma_airways/results/figures_cd4'
colours_f = "/home/ciro/asthma_airways/info/airways_global_colours.csv"
redu = list(umap = c("UMAP_1", "UMAP_2"))
sc_cd4resting_clust = "RNA_snn_res.0.4"
sc_cd4stim_clust = "RNA_snn_res.0.4"
global_objects_f = c(
  bulk_metadata = "/home/ciro/ad_hoc/asthma_airways/raw/mapping_2020_08_27/metadata_filtered_airways_varsubset_mtype.csv",
  bulk_raw = "/home/ciro/ad_hoc/asthma_airways/raw/mapping_2020_08_27/raw.Rdata",
  bulk_tpm = "/home/ciro/ad_hoc/asthma_airways/raw/batchcor/combat_corrected_22401genes_325samples_Sex-plateCorrected.Rdata",
  sc_cd4resting = "/home/ciro/ad_hoc/asthma_airways/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData",
  sc_cd4stim = "/home/ciro/ad_hoc/asthma_airways/results/clust_seurat/stim_cd4_15p_sng_nomacro9n10/clustering/zetInfo/clustCells17PCs_30Ks_0.06667JD.RData"
)
if(!exists("include")) include = names(global_objects_f)
if("all" %in% include) include = names(global_objects_f)
object_names = ls(pattern = "_f$|^redu|_clust")

### Environment/packages ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create(fig_dir, showWarnings = FALSE); setwd(fig_dir)
cat("Working at (fig_dir):", fig_dir, "\n")

cat(cyan("------------------ Packages and functions\n"))
packages_funcs = c(
  "Seurat", "ggplot2", "cowplot", "patchwork",
  "/home/ciro/scripts/handy_functions/devel/file_reading.R",
  "/home/ciro/scripts/handy_functions/devel/filters.R",
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  "/home/ciro/scripts/handy_functions/devel/plots.R")
loaded <- lapply(
  X = packages_funcs,
  FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
}); theme_set(theme_cowplot())

### Global objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat(cyan("------------------ Loading objects\n"))
colours <- readfile(colours_f, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
for(i in names(global_objects_f)){
  command <- paste0(i, " <- readfile(global_objects_f['", i, "'], stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)")
  if(!exists(i) && isTRUE(i %in% include)){
    cat(command, "\n"); object_names = unique(c(object_names, i))
    eval(parse(text = command))
  }
}

sc_cd4resting_ident = list(
  celltype_subset = c(
    "0" = "TRMDP", "1" = "TRMSP", "2" = "TCM", "3" = "Treg", "4" = "TFH",
    "5" = "ThIFNr", "6" = "Activated", "7" = "CellCycle", "8" = "CTLs"
  ),
  colours = c("0" = "#99141B", "1" = "#FFA303", "2" = "#FFE610", "3" = "#21B84C",
    "4" = "#14B9FA", "5" ="#A6B531", "6" = "#000000", "7" = "#A3A3A3", "8" = "#C944C8"),
  order = as.character(c(0:5, 7:8))
)
if(exists("sc_cd4resting")){
  source("/home/ciro/scripts/clustering/R/utilities.R")
  sc_cd4resting_ident = ident_order(sc_cd4resting_ident)
  sc_cd4resting = ident_set_names(
    object = sc_cd4resting,
    ident = sc_cd4resting_ident,
    cluster_column = sc_cd4resting_clust
  )
  sc_cd4resting_ident <- ident_colours(
    sc_cd4resting_ident, mdata = sc_cd4resting@meta.data)
  sc_cd4resting = sc_cd4resting[, rownames(sc_cd4resting@meta.data)]
}
sc_cd4stim_ident = list()

print(format(
  object_names,
  justify = "centre", quote = FALSE)
)
# mdata = FetchData(sc_cd4resting, vars = c(colnames(sc_cd4resting@meta.data), unlist(redu)))
# saveRDS(mdata, file = "metadata.rds"); rm(mdata)
