#!/usr/bin/R

### conda activate tracer ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Requirements
suppressPackageStartupMessages({
  library(data.table)
  library(RColorBrewer)
  library(ggplot2)
  library(cowplot)
  library(eulerr)
  library(UpSetR)
  library(VennDiagram)
  library(grDevices)
})
theme_set(theme_cowplot())
## Python 3 libraries
## argparse
## graphviz

source("/home/ciro/scripts/tracer/plot-tracer-graphviz-general.R")
source('/home/ciro/scripts/handy_functions/devel/filters.R')
source('/home/ciro/scripts/handy_functions/devel/utilities.R')
source('/home/ciro/scripts/handy_functions/devel/plots.R')

## python script
py_path = "/home/amadrigal/SingleCellRNAseq/Simon/guo-2018-reanalysis/tracer/2019-05-13-analysis/plot-tracer-graph-general.py"

# Getting colours
colsname <- "/home/ciro/asthma_pjs/info/airways_global_colours.csv"
grcols <- read.csv(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
colnames(grcols) <- "colors"

## Setting path
group.clonotype = "Study_ID" # define expanded for each of these
annotf_all = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/mixcr_rnseq/asthma_airways_cd4_A-B_clones.txt"
outdir = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/mixcr_rnseq/clonotype_networks_cd4/"

dir.create(outdir); setwd(outdir); system("ls -loh")
## Reading Tracer annotation
## Generating annotation with TCR columns
input_annot <- read.table(annotf_all, sep = "\t", stringsAsFactors = FALSE, header = 1)
rownames(input_annot) <- NULL
undefined_by_marker = ""
selected <- which(!input_annot$sample_name %in% undefined_by_marker)
input_annot <- input_annot[selected, ]
if(is.null(input_annot$cloneCount)){
  tvar <- grep("cloneCount", colnames(input_annot))
  input_annot$cloneCount <- rowMeans(input_annot[, tvar])
}

# QC plot
chain = "TRB"
tvar <- paste0(chain, c(".cloneCount", ".cloneFraction"))
chain_filters = c(10, 0.0025); names(chain_filters) <- tvar
source("/home/ciro/scripts/tracer/plot-tracer-graphviz-general.R")
qc <- clones_qc(
  dat = input_annot[!is.na(input_annot[, tvar[1]]), ],
  xax = tvar[1],
  yax = tvar[2],
  threshold = chain_filters
)
p <- qc$plot
p <- p + scale_x_continuous(trans = "log2", n.breaks = 10) +
  scale_y_continuous(labels = scales::percent, trans = "log2", n.breaks = 10) +
  labs(
    x = expression(Log[2]*"(Clone count)"), y = expression(Log[2]*"(Fraction)")
  )
pdf(paste0("qc_", chain, ".pdf"))
print(p)
dev.off()
pdf(paste0("qc_", chain, "_blank.pdf"), width = 6.6)
print(plot_blank(p))
dev.off()

## Vector of variables to plot in each tracer map
vars.look = c("Stim", "Cell_type", "Study_ID")
groups.look = c("Study_ID", "set")
input_annot$set <- "all_cells"
tvar <- v2cols(names(table(input_annot$Study_ID)))
grcols[names(tvar), 1] <- unname(tvar)
sapply(input_annot[, vars.look, drop = FALSE], table)#names([[3]]) %in% rownames(grcols)

# Both chains must be shared
sharing_definition = "both"
# # Alpha sharing are enough to be clonally related
# sharing_definition = "a"
# Beta sharing are enough to be clonally related
sharing_definition = "b"

# Selecting cells
sselect = list(name = "celltypes4mild", subset = list(c("Disease", "MA"), c("Cell_type", "CD4-Teff", "CD4-TRM", "CD4-TRM-like", "CD4-Treg")))
sselect = list(name = "celltypes4severe", subset = list(c("Disease", "SA"), c("Cell_type", "CD4-Teff", "CD4-TRM", "CD4-TRM-like", "CD4-Treg")))
new_selections = list(c("no_column", "jaja"))
scells1 <- filters_subset_df(
  x = sselect$subset,
  df = input_annot[qc$filtered, ],
  v = TRUE
)
new_selections = lapply(X = unique(input_annot[scells1, ]$Study_ID), FUN = function(x) c("Study_ID", x) )

tvar <- paste0(gsub(paste0(chain, ".clone"), "", names(chain_filters)), chain_filters)
sufix = paste0(c(sharing_definition, group.clonotype, sselect$name, tvar), collapse = "_")
new_root = paste0(sufix, "_per_patient")
if(new_selections[[1]][1] != "no_column") dir.create(new_root)
for(new_selection in new_selections){
  further_select = if(new_selection[[1]] != "no_column") new_selection
  scells <- filters_subset_df(
    x = c(sselect$subset, if(!is.null(further_select)) list(further_select) ),
    df = input_annot[qc$filtered, ],
    v = TRUE
  )

  suffix = if(!is.null(further_select)){
    paste0(new_root, "/", further_select[2])
  }else{
    sufix
  }
  annotf <- paste0(suffix, ".rdata")
  tcr_annot <- tcr_def_expansion(
    annot = input_annot[scells, ],
    chain_1 = "TRA.clonalSequence", chain_2 = "TRB.clonalSequence",
    within_subset = group.clonotype,
    sharing_chain = sharing_definition,
    exclusive = TRUE
  )
  # table(tcr_annot[, c('Freq.TRA.and.TRB', 'GreaterTHAN1')], useNA = 'ifany')
  cat("File:", annotf, ifelse(file.exists(annotf), "exists", "writting"), "\n");
  # save(tcr_annot, file = annotf)

  clone_count_thr = 1
  shared_chain = c(both = "Name.TRA.and.TRB", a = "Name.TRA", b = "Name.TRB")[[sharing_definition]]
  clonecount = c(both = "cloneCount", a = "TRA.cloneCount", b = "TRB.cloneCount")[[sharing_definition]]
  tcr_annot <- tcr_annot[which(tcr_annot[, clonecount] > clone_count_thr), ]

  # Expanded
  clonal_definition = 1
  dinamic_outdir = "./expanded_"
  tcr_annot <- tcr_annot[which(tcr_annot[, "GreaterTHAN1"] > clonal_definition), ]
  if(nrow(tcr_annot) == 0) next
  # # Expanded + Non-expanded
  # clonal_definition = 0
  # dinamic_outdir = "./everything_"

  if(is.null(further_select)) write.csv(tcr_annot, file = paste0(suffix, "_summary.csv"))

  # Venn Diagram
  source('/home/ciro/scripts/handy_functions/devel/overlap.R')
  summ_treat <- make_list(
    x = tcr_annot,
    colname = "Cell_type",
    col_objects = shared_chain
  )
  euler_i <- overlap_list(summ_treat, sep = "&")
  str(summ_treat)
  euler_i <- euler_i[sapply(euler_i, length) > 0]
  euler_i <- sapply(euler_i, length)

  fit2 <- euler(euler_i)
  pdf(paste0(suffix, "_euler.pdf"))
  print(plot(fit2, quantities = TRUE, lty = 1:3, labels = list(font = 4), fill = grcols[names(summ_treat), ]))
  dev.off()

  uddf <- as.data.frame.matrix(table(
    tcr_annot[, c(shared_chain, "Cell_type")]
  ))
  uddf[uddf > 0] <- 1
  pdf(paste0(suffix, "_upset.pdf"), onefile = FALSE)
  print(upset(data = uddf, sets = colnames(uddf), keep.order = TRUE, sets.bar.color = v2cols(names(summ_treat), grcols)))
  dev.off()

  temp <- venn.diagram(summ_treat, fill = v2cols(names(summ_treat), grcols), filename = NULL)
  pdf(paste0(suffix, "_venn.pdf"))
  grid.draw(temp)
  dev.off()
  kk <- file.remove(list.files(pattern = "VennD"))
}

dinamic_dir_i <- paste0(dinamic_outdir, suffix, "/")
dir.create(dinamic_dir_i)
system("ls -loh")
# system("rm */first_pair*")

for(subdivision in groups.look){
  tcr_make_dot(
    annot.all = tcr_annot,
    work.path = dinamic_dir_i,
    var.clonotype = subdivision,
    clonal.definition = clonal_definition,
    name.TRA.and.TRB = shared_chain,
    vars.look = vars.look,
    df.colors = grcols,
    clonal.expansion = "GreaterTHAN1"
  )
  ## Generating PDFs from dot files
  for( dir.work in names(sort(table(tcr_annot[, subdivision]))) ){
      for( var.look in vars.look[!vars.look %in% subdivision] ){
          dot_path<- paste0( dinamic_dir_i, dir.work, "/", var.look, ".dot")
          cat("Path is:", dot_path, "\n")
          if( file.exists(dot_path) ){
              cat("Plotting\n")
              system(paste0( "python ", py_path, " --path_dot ", dot_path))
          }else{
              cat("Path does not exist\n")
          }
      }
      cat("\n")
  }
}
