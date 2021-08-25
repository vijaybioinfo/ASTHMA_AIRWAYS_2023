#!/usr/bin/R

# sinew::makeOxygen(dimred_hvg)

#' @title Dimentional reduction with Highly Variable Features
#' @description Explore counts data with Dimentionality reduction based on
#' highly variable features.
#' @param mdata Metadata.
#' @param edata_counts Matrix of expression.
#' @param edata_norm Normalized data, Default: NULL.
#' @param pc_n Number of Principal Components to use, Default: 2:4.
#' @param perplex_or_neigh t-SNE's Perplexity or UMAP's n_neighbors,
#' Default: c(7, 10, 12).
#' @param showgenes Features or genes to plot, Default: NULL.
#' @param cnames Columns to plot from mdata, Default: NULL.
#' @param output_dir Folder to put the results in, Default: NULL.
#' @param samples_filter Samples to use, Default: NULL.
#' @param features_filter Features or genes to use, Default: NULL.
#' @param pt_size Size of the points on the scatter plot, Default: 3.
#' @param verbose Show progress, Default: TRUE.
#' @param ... Extra parameters for getMostVariableGenes and custom_heatmap.
#' @return NULL
#' @examples
#' \dontrun{
#' if(interactive()){
#'   dimred_hvg(mdata = metadata, edata_counts = counts)
#' }
#' @seealso
#'  \code{\link[stats]{dist}}
#'  \code{\link[Rtsne]{Rtsne}}
#'  \code{\link[uwot]{umap}}
#' @rdname dimred_hvg
#' @export
#' @importFrom stats dist
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap

dimred_hvg <- function(
  mdata,
  edata_counts,
  edata_norm = NULL,
  pc_n = 2:4,
  perplex_or_neigh = c(7, 10, 12),
  showgenes = NULL,
  cnames = NULL,
  output_dir = NULL,
  samples_filter = NULL,
  features_filter = NULL,
  pt_size = 3,
  verbose = TRUE,
  ...
) {
  if(is.null(output_dir)) output_dir <- getwd()
  if(!dir.exists(output_dir)) dir.create(output_dir)
  if(is.null(samples_filter)) samples_filter <- rownames(mdata)
  if(is.null(features_filter)) features_filter <- rownames(edata_counts)
  mdata <- mdata[samples_filter, ]
  edata_counts <- edata_counts[features_filter, samples_filter]
  if(!is.null(edata_norm)) edata_norm <- edata_norm[features_filter, samples_filter]

  fname <- paste0(output_dir, '/_highly_variable_genes', list(...)[["top_n"]], '.pdf')
  if(!is.file.finished(fname)) pdf(fname)
  var_df <- getMostVariableGenes(
    counts = edata_counts,
    plot = !is.file.finished(fname),
    verbose = verbose,
    ...
  )
  graphics.off()
  if(is.null(edata_norm)) edata_norm <- var_df$edata
  keep_genes <- names(var_df$means)[var_df$top_n]
  edata_ss <- log2(edata_norm[keep_genes, ] + 1)
  summary(rowMeans(edata_ss[keep_genes, ]))
  tmp <- if(!is.null(showgenes)){
    grep(paste0(showgenes, collapse = "|"), rownames(edata_norm), value = TRUE)
  }else{ set.seed(27); sample(keep_genes, 3) }
  cnames <- unique(c(cnames, tmp))


  # Clustering
  df <- stats::dist(t(scale(edata_ss))) # t(scale(edata_ss)) same result
  cluster_prefix = "hcluster"
  if(verbose){ print(summary(c(abs(df)))); str(df) }
  sub_grps <- nclust_optimal(
    mat = df, prefix = paste0(output_dir, '/_', cluster_prefix, length(keep_genes)),
    verbose = verbose)
  mdata <- mdata[, grep(cluster_prefix, colnames(mdata), invert = TRUE)]
  for(i in names(sub_grps)){
    clusters <- sub_grps[[i]]
    k_clust = gsub("N", "", i)
    mdata$hcluster = unname(clusters[rownames(mdata)])
    mdata$hcluster[!is.na(mdata$hcluster)] <- paste0(cluster_prefix, "_", mdata$hcluster[!is.na(mdata$hcluster)])
    cnames <- unique(c(cnames, paste0(cluster_prefix, k_clust)))
    colnames(mdata) <- gsub("hcluster$", paste0(cluster_prefix, k_clust), colnames(mdata))
  }

  pca_res <- prcomp(edata_ss, center = TRUE)
  pca_df <- data.frame(pca_res$rotation, stringsAsFactors = FALSE, check.names = FALSE)
  mdata$Sex_Disease <- paste0(mdata$Sex, "_", mdata$Disease)

  fname <- paste0(output_dir, '/_elbow.pdf')
  if(!is.file.finished(fname)){
    pdf(fname)
    plot(1:length(pca_res$sdev), pca_res$sdev)
    dev.off()
  }

  fname <- paste0(output_dir, '/_heatmap.pdf')
  if(!is.file.finished(fname)){
    cnames_i <- cnames[cnames %in% colnames(mdata)]
    pdf(fname, width = 15, height = 15, onefile = FALSE)
    custom_heatmap(
      object = edata_norm[keep_genes, ],
      annoc = mdata[, cnames_i, drop = FALSE],
      # use_mean = cnames_i[1],
      feature_order = "pca",
      verbose = verbose,
      ...
    )
    graphics.off()
  }

  for(pc in pc_n){
    for(px in perplex_or_neigh){
      set.seed(27)
      drdata <- if(grepl("tsne", output_dir)){
        dim_axes = c("t-SNE 1", "t-SNE 2")
        Rtsne::Rtsne(dist(pca_df[, 1:pc]), perplexity = px, pca = FALSE, theta = 0, is_distance = TRUE, max_iter = 1000)$Y
      }else{
        dim_axes = c("UMAP 1", "UMAP 2")
        uwot::umap(t(edata_ss), pca = pc, n_neighbors = px)
      }
      ddf <- data.frame(x.var = drdata[, 1], y.var = drdata[, 2])
      rownames(ddf) <- rownames(pca_df)
      for(cname in cnames){
        fname <- paste0(output_dir, '/hv', length(keep_genes), 'genes_', paste0(pc, "PCs_", px), 'PX_', cname, '_size', pt_size, '.pdf')
        if(verbose) cat(fname, '\n'); if(is.file.finished(fname)) next
        ddf$col <- if(cname %in% colnames(mdata)) factormix(mdata[, cname]) else "NONE"
        ddf$col <- if(cname %in% rownames(edata_norm)) log2(c(edata_norm[cname, rownames(ddf)]) + 1) else ddf$col
        scolor <- if(is.numeric(ddf$col)){
          scale_colour_gradientn(colours = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000'))
        }else{ scale_colour_manual(values = v2cols(levels(ddf$col), list(...)[["couls"]])) }

        gg.d <- ggplot(ddf, aes(x = x.var, y = y.var)) +
            geom_point(size = pt_size, shape = 19, alpha = 0.7, aes(colour = col)) + scolor +
            labs(
              title = gsub("_|.pdf", " ", basename(fname)),
              x = dim_axes[1], y = dim_axes[2], colour = "Group"
            )
        pdf(fname, width = 8, height = 7)
        print(gg.d)
        dev.off()
      }
    }
  }; return(invisible(x = NULL))
}
