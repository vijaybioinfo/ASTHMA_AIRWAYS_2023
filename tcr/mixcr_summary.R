#!/usr/bin/R

#################
# MiXCR results #
#################

# This script creates a summary file from a series of MiXCR outputs
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
}); theme_set(theme_cowplot())
source("/home/ciro/scripts/tracer/plot-tracer-graphviz-general.R")
source('/home/ciro/scripts/handy_functions/devel/utilities.R')
source('/home/ciro/scripts/handy_functions/devel/filters.R')
system("ls -loh")

colsname <- "/home/ciro/asthma_pjs/info/airways_global_colours.csv"
grcols <- read.csv(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
colnames(grcols) <- "colors"

mixcr_output = "asthma_airways_cd4"
metadataf = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/raw/mapping_2020_08_27/metadata_filtered_airways_varsubset.csv"
var_names = c("Plate", "Stim", "Disease", "Study_ID", "Cell_type", "Severity")
sselect = list(c("Cell_type", "-CD4-TFH", "-CD8-TRM"))
clone_info = c("cloneCount", "cloneFraction", "clonalSequence", "aaSeqCDR3", "bestVHit", "bestJHit", "bestDHit")
tcrtype = c("A", "B")
heatmap_cols = c("Cell_type", "Stim", "Disease", "Study_ID")
box_cols = c("Disease", "value", "variable", "Cell_type")

setwd("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/mixcr_rnseq")
system("ls -loh")
clone_files <- list.files(path = mixcr_output, pattern = "clones.txt", full.names = TRUE, recursive = TRUE)
cat(length(clone_files), "clonotype files\n")
cat(head(clone_files), sep = "\n")

clone_definition = c("clonalSequence", "bestVHit", "bestJHit", "bestDHit")
clone_nested = c("cloneFraction", "cloneCount", "clonalSequence")
dtype = "all"
exchains = list(TRA = "TRD", TRD = "TRA")
clones_list <- lapply(
  X = clone_files,
  FUN = function(x){
    y <- data.frame(data.table::fread(x), stringsAsFactors = FALSE)[, clone_info]
    if(nrow(y) > 0){
      y <- cbind(sample_name = basename(dirname(x)), y)
      y$sample_name <- as.character(y$sample_name)
    }
    y$key <- make.names(do.call(paste0, c(y[, clone_definition, drop = FALSE])))
    chain_files <- paste0(dirname(x), "/clones_TR", tcrtype, ".txt")
    names(chain_files) <- paste0("TR", tcrtype)
    chain_freq <- lapply(
      X = chain_files,
      FUN = function(chain){
        z <- data.frame(data.table::fread(chain), stringsAsFactors = FALSE)[, clone_info]
        rownames(z) <- make.names(do.call(paste0, c(z[, clone_definition])))
        # all(z$cloneFraction == z$cloneCount / sum(z$cloneCount))
        if(grepl("TRA|TRD", chain)){
          exchain <- exchains[which(sapply(names(exchains), function(x) grepl(x, basename(chain)) ))]
          z <- z[!grepl(exchain[[1]], z$bestJHit), ]
          z <- z[!grepl(exchain[[1]], z$bestDHit), ]
          z$cloneFraction = z$cloneCount / sum(z$cloneCount)
        }
        z[, clone_nested, drop = FALSE]
    })
    if(length(chain_freq) == 0) return(y)
    for(i in 1:length(chain_freq)){
      chain_freq_i <- chain_freq[[i]]
      colnames(chain_freq_i) <- paste0("TR", tcrtype[i], ".", colnames(chain_freq_i))
      y <- cbind(y, chain_freq_i[y$key, ])
    }
    y <- y[, !colnames(y) %in% "key"]
    y
  }
)
chains_found = unique(unlist(lapply(clones_list, function(x) gsub("(...).*", "\\1", x[, 'bestVHit']) )))
chains_found; chains_found = paste0("TR", tcrtype)
clones_df <- make_clone_pairs_rbind(l = clones_list, verbose = TRUE)

besties = grep(pattern = "best", x = colnames(clones_df), value = TRUE)
chain_segments = reshape2::melt(lapply(
  X = setNames(chains_found, gsub("TR", "", chains_found)),
  FUN = function(x){
    k <- lapply(
      X = setNames(besties, gsub("best|Hit", "", besties)),
      FUN = function(y){
        z <- clones_df[!is.na(clones_df[, paste0(x, ".clonalSequence")]), ]
        unique(gsub("(...).*", "\\1", z[!is.na(z[, y]), y]))
    })#; unlist(k)
})); colnames(chain_segments) <- c("Gene", "Segment", "TCR")
chain_segments[, 3:1]
unique(gsub("(...).*", "\\1", clones_df[!is.na(clones_df[, "TRA.clonalSequence"]), "bestDHit"]))

tvar <- grep(pattern = "cloneFraction", x = colnames(clones_df), value = TRUE)
for(i in unique(clones_df[, 1])[1:2]) # per sample
  for(j in tvar) cat(j, "=>", sum(clones_df[clones_df[, 1] == i, j], na.rm = TRUE), "\n")

metadata <- read.csv(metadataf, row.names = 1, stringsAsFactors = FALSE)
clinic_data = read.csv('/home/ciro/asthma_pjs/info/airways_metadata_donors_processed.csv', stringsAsFactors = FALSE, row.names = 1)
head(clinic_data[, 1:3])
str(metadata)
all(metadata$Study_ID %in% rownames(clinic_data))
metadata$Severity <- clinic_data[metadata$Study_ID, "Adapted.Asthma.Severity.Score"]
# clones_annotation <- cbind(clones_df, metadata[clones_df$sample_name, var_names])
clones_annotation <- cbind(metadata[clones_df$sample_name, var_names], clones_df)
clones_annotation <- clones_annotation[filters_subset_df(sselect, clones_annotation, verbose = TRUE), ]
rownames(clones_annotation) <- NULL
head(clones_annotation)
tvar <- grep(pattern = "bestVHit", x = colnames(clones_annotation), value = TRUE)
head(clones_annotation[!is.na(clones_annotation[, tvar[1]]) & !is.na(clones_annotation[, tail(tvar,1)]), ])
fname <- paste0(mixcr_output, "_", dtype, "_clones.txt")
write.table(x = clones_annotation, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
system(paste0("head ", fname))

tvar <- grep("sample_name|\\.clonalSequence|\\.cloneFraction|\\.cloneCount", colnames(clones_df))
clones_annotation_ext = reshape2::melt(clones_df[, tvar])
clones_annotation_ext <- clones_annotation_ext[!is.na(clones_annotation_ext$value), ]
# colnames(clones_annotation_ext) <- c("sample_name", "chain", "cloneFraction")
clones_annotation_ext$chain <- gsub("\\..*", "", clones_annotation_ext$variable)
clones_annotation_ext$variable <- gsub(".*\\.", "", clones_annotation_ext$variable)
clones_annotation_ext$clonalSequence <- clones_annotation_ext$TRA.clonalSequence
tvar <- is.na(clones_annotation_ext$clonalSequence)
clones_annotation_ext$clonalSequence[tvar] <- clones_annotation_ext$TRB.clonalSequence[tvar]
clones_annotation_ext <- clones_annotation_ext[, !grepl("\\.clonalSe", colnames(clones_annotation_ext))]
clones_annotation_ext1 = reshape2::dcast(
  clones_annotation_ext, sample_name + chain + clonalSequence ~ variable)
clones_annotation_ext1 <- cbind(clones_annotation_ext1, metadata[clones_annotation_ext$sample_name, var_names])
rownames(clones_annotation_ext1) <- NULL

id_cols = c("sample_name", "chain")
id_cols = c("Study_ID", "Cell_type", "chain")
clones_annotation_ext1$id = do.call("paste", c(clones_annotation_ext1[, id_cols], sep = "_"))
ddf = clones_annotation_ext1 %>% group_by(id) %>% # dplyr::summarize(
  mutate(
    frac_gt5pct = sum(cloneFraction > 0.05),
    top_clone = head(clonalSequence[which(cloneCount == max(cloneCount))], 1),
    top_clone_freq = max(cloneFraction)
  )
unique(ddf[ddf$id == tail(ddf$id, 1), "top_clone_freq"])
ddf <- ddf[!duplicated(ddf$id), ]; ddf <- ddf[, !grepl("^clon", colnames(ddf))]
ddf$id = NULL
head(as.data.frame(ddf))
fname <- paste0("summary_", mixcr_output, "_", paste0(id_cols, collapse = "_"), ".txt")
write.table(x = ddf, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)

### Plots ### ------------------------------------------------------------------
chain_filters = c(10, 0.0025)
filtered <- sapply(
  X = paste0("TR", tcrtype),
  FUN = function(chain){
    chain_count <- grep(paste0(chain, ".*Count"), colnames(clones_annotation), value = TRUE)
    chain_fracs <- grep(paste0(chain, ".*Fra"), colnames(clones_annotation), value = TRUE)
    print(c(chain_count, chain_fracs))
    df_i <- clones_annotation[, c(chain_count, chain_fracs)]
    df_i <- df_i[complete.cases(df_i), ]
    tvar <- which(df_i[, chain_count] > chain_filters[1] & df_i[, chain_fracs] > chain_filters[2])
    print(length(tvar) / nrow(df_i))
    rownames(df_i[tvar, ])
})
str(filtered); sum(sapply(filtered, length))
str(sapply(filtered, unique)); sum(sapply(filtered, function(x) length(unique(x))))
length(unique(unlist(filtered)))
chain_cols <- grep("TR.*Fra|TR.*Coun", colnames(clones_annotation), value = TRUE)
clones_annotation_u <- clones_annotation[unlist(filtered, use.names= FALSE), ]
nrow(clones_annotation_u)/nrow(clones_annotation)
head(clones_annotation_u[, chain_cols])
head(clones_annotation_u[rowSums(clones_annotation_u[, chain_cols]>0, na.rm = T) > 1, chain_cols])
summary(clones_annotation_u[, chain_cols])

dir.create(paste0(mixcr_output, "_heatmap/"))
for(chain_col in chain_cols){
  vGenes <- reshape::cast(
    data = clones_annotation_u[!is.na(clones_annotation_u[, chain_col]), ],
    bestVHit ~ sample_name,
    value = chain_col, sum
  )
  rownames(vGenes) = as.character(vGenes$bestVHit)
  vGenes$bestVHit = NULL; if(grepl("ount", chain_col)) vGenes <- log2(vGenes + 1)
  annoc = metadata[colnames(vGenes), heatmap_cols]
  data.table::setorderv(x = annoc)
  ssamples = rownames(annoc)#[order(annoc[, heatmap_cols[1]]), ])
  annoc = annoc[ssamples, ]
  vGenes <- vGenes[, ssamples]
  mat2plot <- as.matrix(t(scale(t(vGenes)))); colnames(mat2plot) <- colnames(vGenes)
  topz <- max(c(min(abs(c(range(mat2plot), 2))), 1))
  mat2plot[mat2plot > topz] <- topz; mat2plot[mat2plot < (-topz)] <- -topz;
  RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
  fname <- paste0(mixcr_output, "_heatmap/", dtype, "_", chain_col, "")
  pdf(paste0(fname, ".pdf"), height = 10, width = 9)
  pheatmap::pheatmap(
    mat = mat2plot, color = RedBlue, scale = "none", show_colnames = FALSE,
    annotation_col = annoc, cluster_cols = FALSE,
    annotation_colors = lapply(annoc, v2cols, sour = grcols)
  )
  dev.off()
}

chain_cols <- grep("TR.*Fra", colnames(clones_annotation), value = TRUE)
for(chain in paste0("TR", tcrtype)){
  print(chain)
  print(sum(!is.na(clones_annotation_u[, paste0(chain, ".cloneFraction")])))
  print(sum(grepl(chain, clones_annotation_u$bestVHit) | (grepl(chain, clones_annotation_u$bestJHit) | grepl(chain, clones_annotation_u$bestDHit))))
}
var2plot <- c("sample_name", var_names)
df2plot <- reshape2::melt(clones_annotation_u[, c(var2plot, chain_cols)], id.vars = var2plot)
df2plot <- df2plot[!is.na(df2plot$value), ]
str(df2plot)
df2plot$variable <- gsub("\\..*", "", as.character(df2plot$variable))
table(df2plot$variable)
summary(df2plot$value)

p <- ggplot(
  data = df2plot,
  mapping = aes_string(x = box_cols[1], y = box_cols[2], color = box_cols[1], fill = box_cols[1])
) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), color = "black", show.legend = FALSE) +
  geom_boxplot(width=0.1, fill = "white", color = "white", alpha = 0.25, outlier.shape = NA, show.legend = FALSE) +
  scale_color_brewer(palette="Paired") +
  labs(x = NULL, y = expression(Log[2]*"(Fraction)"), color = "Chain") +
  scale_y_continuous(labels = scales::percent, trans = "log2")
if(!is.na(box_cols[3])){
  myformula = paste0(box_cols[3:length(box_cols)], collapse = "+")
  p <- p + facet_wrap(facets = as.formula(paste0("~", myformula)))
}
tmp <- box_cols[box_cols %in% colnames(clones_annotation_u)]
if(length(tmp)){
  p <- p + scale_color_manual(values = v2cols(df2plot[, tmp[1]], grcols))
  p <- p + scale_fill_manual(values = v2cols(df2plot[, tmp[1]], grcols))
  p <- p + stat_compare_means(method = "wilcox.test")
}

mychains <- paste0(gsub("TR", "", unique(df2plot$variable)), collapse = "-")
fname <- paste0(c(paste0(mixcr_output, "_", dtype, "_", gsub(".*\\.", "", chain_cols[1])),
  mychains, tmp), collapse = "-")
pdf(paste0(fname, ".pdf"), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
print(plot_blank(p))
dev.off()

library(circlize)
# I aim to capture the chains with the highest sharing and samples with the gratest
# amount of chains
# Create adjancy matrix
subset_i <- dtype
clones_annotation_u_i <- clones_annotation_u
subset_i <- "TRG"
subset_i <- "TRD"
clones_annotation_u_i <- clones_annotation_u[!is.na(clones_annotation_u[, paste0(subset_i, ".cloneFraction")]), ]

min_samples = 2 # clonal
min_chains = 0.1
#   s1 s2 s3
# c1 1  0  0
# c2 1  0  1
# c3 1  1  1
# c4 1  1  0
# c5 1  0  1

df2chord <- t(table(clones_annotation_u_i[, c("sample_name", "clonalSequence")]))
head(sort(rowSums(df2chord), decr = T))
summary(rowSums(df2chord))
summary(colSums(df2chord))
write.csv(df2chord, file = paste0(mixcr_output, "_", subset_i[1], "_chord.csv"))

df2chord <- df2chord[rowSums(df2chord>0) >= min_samples, ]
df2chord <- df2chord[, colSums(df2chord>0) >= (nrow(df2chord) * min_chains)]
# df2chord <- df2chord[rowSums(df2chord) >= min_samples, ]
df2chord <- as.matrix(as.data.frame.matrix(df2chord))
str(df2chord)
hc <- hclust(dist(t(df2chord)))
df2chord <- df2chord[, hc$labels[hc$order]]
hc <- hclust(dist(df2chord))
df2chord <- df2chord[hc$labels[hc$order], ]

samples <- metadata[colnames(df2chord), "Disease"]
names(samples) <- colnames(df2chord)
chains <- rownames(df2chord)

fname <- paste0(mixcr_output, "_", subset_i[1], "_chord_minS", min_samples, "_pctC", min_chains)
write.csv(df2chord, file = paste0(fname, ".csv"))

df2chord_mat <- df2chord
colnames(df2chord_mat) <- paste0("S", 1:ncol(df2chord_mat))
rownames(df2chord_mat) <- paste0("C", 1:nrow(df2chord_mat))
grid.col = v2cols(samples, grcols)[samples]
names(grid.col) = colnames(df2chord_mat)
pdf(paste0(fname, ".pdf"))
chordDiagram(df2chord_mat, grid.col = grid.col)
dev.off()

df2chord_df = data.frame(
  from = rep(colnames(df2chord), each = nrow(df2chord)),
  to = rep(rownames(df2chord), times = ncol(df2chord)),
  stringsAsFactors = FALSE
)
df2chord_df_id <- do.call(paste, c(df2chord_df[, 1:2], sep = "_"))
df2chord_df <- df2chord_df[complete.cases(df2chord_df), ]
df2chord_df$to <- paste0("C", as.numeric(factor(df2chord_df$to, chains)))
# df2chord_df$from <- factor(df2chord_df$from, names(samples)) # df2chord_df <- df2chord_df[order(df2chord_df$from), ]
df2chord_df$from <- paste0("S", as.numeric(factor(df2chord_df$from, names(samples))))
df2chord_df <- df2chord_df[, 2:1]
colnames(df2chord_df) <- c("from", "to")
head(df2chord_df)

for(plotthis in c("cloneCount", "cloneFraction")){
  sample_chain <- clones_annotation_u_i[, plotthis]
  names(sample_chain) <- do.call(paste, c(clones_annotation_u_i[, c("sample_name", "clonalSequence")], sep = "_"))
  df2chord_df$value <- sample_chain[df2chord_df_id]
  df2chord_df$value[is.na(df2chord_df$value)] <- 0
  pdf(paste0(fname, "_", plotthis, ".pdf"))
  chordDiagram(df2chord_df, grid.col = grid.col)
  dev.off()
  list.files(pattern = "chord")
}
