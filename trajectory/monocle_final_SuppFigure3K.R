#!/usr/bin/R

#############
# Monocle 3 #
#############

#install.packages("R.utils")

args = R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
args[[1]]<-"cd4_resting2"
suppressPackageStartupMessages({
  # BiocManager::install(c('batchelor','terra')) 'leidenbase', 'rsample', 'spdep', 'terra'
  library(monocle3)
})

get_earliest_principal_node <- function(
  cds,
  time_bin = 1,
  revert = FALSE
){
  time_bin_i <- names(which.max(table(cds[[time_bin]])))
  cell_ids <- which(colData(cds)[, time_bin] == time_bin_i)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  tvar <- if(revert){
    which.min(table(closest_vertex[cell_ids,]))
  }else{ which.max(table(closest_vertex[cell_ids,])) }
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(tvar))]
}

redu = c("UMAP", "tSNE", "PCA", "LSI", "Aligned")[1]
source("/Users/fcastaneda/Desktop/Trajectory/file_reading.R")
source("/Users/fcastaneda/Desktop/Trajectory/filters.R")
source("/Users/fcastaneda/Desktop/Trajectory/plots.R")
library(ggplot2)
cat("Dataset:", args[[1]], "\n");
if(isTRUE(grepl("cd4_resting2", args[[1]], ignore.case = TRUE))){
  edataf = setNames(nm = "/Users/fcastaneda/Desktop/Trajectory/clustCells20PCs_30Ks_0.06667JD.RData")
  npcs = 20; cclust = "RNA_snn_res.0.4"; root_clust = c("2")
}
hvgf = NULL; selectss = list(c(cclust, "0", "1", "2")) ; sufix = args[[1]]

out_dir = paste0("/Users/fcastaneda/Desktop/Trajectory/",
  sub("\\.rds$", "", sub('_object_lock_mean.*_pct', '_pct', basename(edataf))), sufix, "_only_TRMs_3")

out_dir = gsub("\\.rds", "", out_dir)
dir.create(out_dir, recursive = TRUE); setwd(out_dir)
cat("Working at:", getwd(), "\n")


## Preprocessing

cds_f = list.files(pattern = "object.r")
if(!isTRUE(file.exists(cds_f))){
  cat("Running analysis\n")
  edata <- readfile(edataf)
  scells = filters_subset_df(selectss, edata@meta.data, v = TRUE)
  edata <- edata[, scells]
  if(!is.null(hvgf)){
    hvgdat <- read.csv(hvgf, stringsAsFactors = FALSE, row.names = 1)
    sum(hvgdat[, 'variable']); tvar <- rownames(hvgdat[hvgdat[, 'variable'], ])
    cat(commas(tvar), '\n')
    Seurat::VariableFeatures(edata) <- tvar
  }

  ## The first step in working with Monocle 3 is to load up your data into Monocle 3's main class, cell_data_set:
  gene_annotation <- edata@assays$RNA@meta.features
  gene_annotation <- cbind(gene_short_name = rownames(gene_annotation), gene_annotation)
  cds <- new_cell_data_set(
    expression_data = edata@assays$RNA@counts,
    cell_metadata = edata@meta.data,
    gene_metadata = gene_annotation
  )

  ## Step 1: Normalize and pre-process the data
  cds <- preprocess_cds(
    cds = cds,
    method = "PCA",
    num_dim = npcs,
    norm_method = "log",
    use_genes = Seurat::VariableFeatures(edata),
    pseudo_count = NULL,
    scaling = TRUE,
    verbose = TRUE
  )

  #cds <- align_cds(cds, residual_model_formula_str = "~nCount_RNA+percent.mt", alignment_group = NULL)

    # residual_model_formula_str = "~nCount_RNA+percent.mt",
    # alignment_group = NULL,


  p <- plot_pc_variance_explained(cds)
  pdf("1_variance_explained_elbow.pdf")
  print(p)
  dev.off()
  # Taking objects PCA embeddings
  # evolqg::PCAsimilarity(
  #   cds@reducedDims@listData[["PCA"]], edata@reductions[["pca"]]@cell.embeddings
  # )
  #cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings
  cds@int_colData$reducedDims[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings
  p <- plot_pc_variance_explained(cds)
  pdf("1_variance_explained_elbow_seurat.pdf")
  print(p)
  dev.off()
  #cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings[, 1:npcs]
  cds@int_colData$reducedDims[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings[, 1:npcs]


  ## Step 2: Remove batch effects with cell alignment
  # cds <- align_cds(cds, alignment_group = "batch")
  cds <- align_cds(
    cds = cds,
    preprocess_method = "PCA",
    residual_model_formula_str = "~nCount_RNA+percent.mt",
    verbose = TRUE
  )

  ## Step 3: Reduce the dimensions using UMAP
  # Taking objects UMAP embeddings

  cds@int_colData$reducedDims[["UMAP"]] <- edata@reductions[["umap"]]@cell.embeddings
  #cds@reducedDims@listData[["UMAP"]] <- edata@reductions[["umap"]]@cell.embeddings

  # cds <- reduce_dimension(
  #   cds = cds,
  #   max_components = 2,
  #   reduction_method = redu,
  #   preprocess_method = NULL,
  #   umap.metric = "cosine",
  #   umap.min_dist = 0.1,
  #   umap.n_neighbors = 15L,
  #   umap.fast_sgd = TRUE,
  #   umap.nn_method = "annoy",
  #   cores = 1,
  #   verbose = TRUE
  # )

  ## Step 4: Cluster the cells
  # Monocle is able to learn when cells should be placed in the same trajectory as
  # opposed to separate trajectories through its clustering procedure. Recall that
  # we run cluster_cells(), each cell is assigned not only to a cluster but also to
  # a partition. When you are learning trajectories, each partition will eventually
  # become a separate trajectory
  # Seurat::FindClusters is: 1 = original Louvain algorithm
  cds <- cluster_cells(
    cds = cds,
    reduction_method = redu,
    k = 180,
    cluster_method = c("leiden", "louvain")[2],
    num_iter = 2,
    partition_qval = 0.05,
    weight = FALSE,
    resolution = NULL,
    random_seed = NULL,
    verbose = TRUE
  )

  ## Step 5: Learn a graph
  cds <- learn_graph(
    cds = cds,
    use_partition = TRUE,
    close_loop = TRUE,
    learn_graph_control = list(minimal_branch_len=30),
    verbose = TRUE
  ) #here 1:27 pm

  p <- plot_cells(
    cds = cds,
    color_cells_by = "partition",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )

  pdf("2_partitions.pdf")
  print(p)
  dev.off()

  table(cds@clusters@listData$UMAP$partitions)

  if(!'partitions_bk' %in%  names(cds@clusters@listData$UMAP))
    cds@clusters@listData$UMAP$partitions_bk <- cds@clusters@listData$UMAP$partitions


  #why merging all in 1 partition??!!!!
 #tvar <- factor(rep('1', ncol(cds))); names(tvar) <- names(cds@clusters@listData$UMAP$partitions)
  # cds@clusters@listData$UMAP$partitions <- tvar

  ## Step 5: Learn a graph - merge all partitions
   #cds <- learn_graph(
   # cds = cds,
   # use_partition = TRUE,
   # close_loop = TRUE,
   #  learn_graph_control = NULL,
    # verbose = TRUE
   #)
  # p <- plot_cells(
  #   cds = cds,
   #  color_cells_by = "partition",
   #  label_cell_groups = FALSE,
  #   label_leaves = FALSE,
   #  label_branch_points = FALSE
   #)

  pdf("2_partitions2.pdf")
  print(p)
dev.off()
  table(cds@clusters@listData$UMAP$partitions) #1:42 pm

  pdf("3_select_root.pdf", width = 9, height = 7)
  p <- plot_cells(
    cds = cds,
    reduction_method = redu,
    color_cells_by = cclust,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 1.5
  )
  print(p)
dev.off()
}else{
  cat("Object file found:", cds_f, "\n")
  #cds = readfile(cds_f)
  cds=readfile("/Users/fcastaneda/Desktop/Trajectory/clustCells20PCs_30Ks_0.06667JD.RDatacd4_resting2_only_TRMs_3/object.rds")
}

## Step 6: Order cells
if(is.null(root_clust)) root_clust = names(which.max(table(cds[[cclust]])))

root_cell <- rownames(colData(cds)[as.character(colData(cds)[[cclust]]) == root_clust, ])[1]

cds <- order_cells(
  cds = cds,
  reduction_method = "UMAP",
  # root_pr_nodes = get_earliest_principal_node(cds, time_bin = cclust),
  root_cells = root_cell,
  verbose = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  trajectory_graph_color="#e02914",
  trajectory_graph_segment_size=1.8,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5, cell_size = 1, alpha=0.7
)
pdf("4_pseudotime.pdf")
print(p)
dev.off()
pdf("4_pseudotime_blank.pdf")
print(plot_blank(p))
dev.off()

pdf("4_pseudotime_ONLY.pdf")
p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE,
  # trajectory_graph_color="wheat2", # #e02914 white
  # trajectory_graph_segment_size= 4,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 4, cell_size = 1, alpha=1
)
p$layers[[1]]$aes_params$colour <- "white"

print(p)
dev.off()

for(cname in c(grep("RNA_snn_res", colnames(colData(cds)), value = TRUE))){
  p <- plot_cells(
    cds = cds,
    reduction_method = redu,
    color_cells_by = cname,
    group_cells_by = 'cluster',
    label_cell_groups = !TRUE,
    label_groups_by_cluster = !TRUE,
    label_branch_points = TRUE, # black circles
    label_roots = TRUE, # white circles
    label_leaves = !TRUE, # gray circles
    graph_label_size = 3,
    cell_size = 1
  )
  pdf(paste0("5_", casefold(redu), "_", cname, ".pdf"), width = 7.5)
  print(p)
  dev.off()
}
graphics.off()


cname <-  "RNA_snn_res.0.4"
#pdf(paste0("5_", casefold(redu), "_", cname, "_COLORS_UMAP_ALPHA_mod.pdf"), width = 7.5)
  colors_sc_cd4resting_ident<-c('0'="#99141B", '1'="#FFA303", '2'="#FFE610", '3'="#21B84C",'4'="#14B9FA", '5'="#A6B531",'8'="#C944C8", '7'="#A3A3A3", '6'="#A3A3A3")
options(ggplot2.discrete.colour= colors_sc_cd4resting_ident)

p <- plot_cells(
  cds = cds,
  reduction_method = redu,
  color_cells_by = cname,
  group_cells_by = 'cluster',
  label_cell_groups = FALSE,
  label_groups_by_cluster = !TRUE,
  label_branch_points = FALSE, # black circles
  label_roots = TRUE, # white circles
  label_leaves = FALSE, # gray circles
  cell_size = 3.2,
  graph_label_size = 4,
  trajectory_graph_segment_size= 1,
  alpha=0.05,
  cell_stroke = 0
  ) #+ geom_point(colour="blue")
#+ scale_color_manual(values = colors_sc_cd4resting_ident, name = "cluster")
#+ scale_colour_manual(values = colors_sc_cd4resting_ident, name = "cluster")
p$layers[[1]]$aes_params$colour <- 'white'

pdf(paste0("5_", casefold(redu), "_", cname, "_COLORS_UMAP_ALPHA_mod.pdf"), width = 7.5)
print(p)
dev.off()


p <- plot_cells(
  cds = cds,
  reduction_method = redu,
  color_cells_by = cclust,
  group_cells_by = 'cluster',
  label_cell_groups = !TRUE,
  label_groups_by_cluster = !TRUE,
  label_branch_points = !TRUE, # black circles
  label_roots = !TRUE, # white circles
  label_leaves = !TRUE, # gray circles
  cell_size = 1
)
pdf(paste0("6_", casefold(redu), "_", cclust, ".pdf"), width = 7.5)
print(p)
dev.off()
pdf(paste0("6_", casefold(redu), "_", cclust, "_blank.pdf"))
print(plot_blank(p, "ext|nstanc"))
dev.off()

saveRDS(cds, file = "object.rds")
cat("Finished\n"); #q(save = "no")
