#!/usr/bin/R

#####################
# MiXCR exploration #
#####################

# This script creates plots exploring our bulk TCR data

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})
theme_set(theme_cowplot())
source('/home/ciro/scripts/handy_functions/devel/filters.R')
source('/home/ciro/scripts/handy_functions/devel/utilities.R')

colours_f = "/home/ciro/asthma_pjs/info/airways_global_colours.csv"
colours <- read.csv(colours_f, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

setwd("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/a1_final_figures_cd4")
result_id = "tcr_resting_"
sselect = list(c("Stim", "Unstim"), c("Cell_type", "CD4-Teff", "CD4-TRM", "CD4-TRM-like"))
result_id = "tcr_stim_"
sselect = list(c("Stim", "Stim"), c("Cell_type", "CD4-Teff", "CD4-TRM", "CD4-TRM-like"))
clones_annotation_f = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/mixcr_rnseq/asthma_airways_cd4_all_clones.txt"
clones_annotation_f = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/mixcr_rnseq/asthma_airways_cd4_A-B_clones.txt"
clone_name = c("Study_ID", "TRA.clonalSequence", "TRB.clonalSequence")[-1]
yax = "cloneFraction"
yax = "log2(cloneCount)"
clone_count_thr = 1

clones_annotation <- read.table(clones_annotation_f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(clones_annotation)

if(is.null(clones_annotation$cloneCount)){
  tvar <- grep("cloneCount", colnames(clones_annotation))
  clones_annotation$cloneCount <- rowMeans(clones_annotation[, tvar])
}
clones_annotation$cloneName <- do.call(paste, c(clones_annotation[, clone_name], sep = "_"))
scells <- filters_subset_df(
  x = sselect,
  df = clones_annotation[clones_annotation$cloneCount > clone_count_thr, ],
  v = TRUE
)

p <- ggplot(
    data = clones_annotation[scells, ],
    mapping = aes_string(x = "Cell_type", y = yax, fill = "Cell_type")
  ) +
  facet_wrap(~Disease, nrow = 1) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  geom_line(mapping = aes(group = cloneName), color = "blue") +
  geom_jitter(width = 0.3, size = 0.4, color = "#595959", alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(values = v2cols(clones_annotation[scells, "Cell_type"], colours)) +
  labs(
    x = NULL, y = expression("Log"[2]*"(MiXCR clone count)"),
    caption = "Blue lines: shared clonal sequence",
    subtitle = paste("Clones with a clone count >", clone_count_thr)
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf(paste0(result_id, yax, ".pdf"))
print(p)
dev.off()
