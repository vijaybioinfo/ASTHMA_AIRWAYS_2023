#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-07-05
# ---

# This script will create plots for any figure from our cd4 data
# Each figure's code will then probably be moved to its correct file

### Environment/packages and global objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# include = "sc_cd4resting"
source("/home/ciro/asthma_airways/scripts/final_figures_cd4/global.R")
sc_cd4resting@meta.data$sex_disease = factor(ident_combine(
  sc_cd4resting@meta.data, c("orig.Sex", "orig.asthma")),
  c("Male-MA", "Male-SA", "Female-MA", "Female-SA"))
tvar <- sc_cd4resting@assays$RNA@counts["GZMB", ] > 0
sc_cd4resting@meta.data$tag_GZMB = ifelse(tvar, "GZMBp", "GZMBn")
tvar <- sc_cd4stim@assays$RNA@counts["GZMB", ] > 0
sc_cd4stim@meta.data$tag_GZMB = ifelse(tvar, "GZMBp", "GZMBn")

{ cat(red_bold("#### QC metrics #### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("qc", showWarnings = FALSE)
  qc_metrics = grep("_RNA|mt", colnames(sc_cd4stim@meta.data), value = TRUE)
  pp <- lapply(qc_metrics, function(x){
    violin(sc_cd4stim@meta.data, "orig.stim", yax = x) +
      theme(legend.position = "none") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 7))
  })
  pdf("qc/stim_metrics.pdf", width = 14, height = 7)
  print(cowplot::plot_grid(plotlist = pp, ncol = 3))
  graphics.off()
  pp = lapply(pp, function(x) x + theme_axes() )
  pdf("qc/stim_metrics_blank.pdf", width = 14, height = 7)
  print(plot_blank(cowplot::plot_grid(plotlist = pp, ncol = 3)))
  graphics.off()
}

{ cat(red_bold("#### Main dim. red. #### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  result_id = "dim_reduction/"; dir.create("dim_reduction", showWarnings = FALSE)
  source("/home/ciro/scripts/clustering/R/plotting.R")
  fconfigs = list(
    list(object = "sc_cd4resting",
      vars2col = c("cluster.celltype_subset"), plot_main = TRUE,
      facet = c("sex_disease"), redu = redu[[1]],
      filters = NULL, sample_even = TRUE, per_facet = TRUE,
      cols = sc_cd4resting_ident$col)
  )

  for (fconfig in fconfigs) {
    scells = filters_subset_df(fconfig$filters,
      eval(parse(text = fconfig$object))@meta.data, verbose = TRUE)
    ddf <- Seurat::FetchData(eval(parse(text = fconfig$object)),
      vars = unique(c(fconfig$redu, fconfig$vars2col, fconfig$facet)))
    ddf <- ddf[scells, ]
    for(i in fconfig$vars2col){
      fname0 <- paste0(result_id, fconfig$object, "_", gsub("orig\\.", "", i))
      cat(cyan("-", fname0, "\n"))
      ddf$tmp <- droplevels(factor(ddf[, i]))
      if(isTRUE(fconfig$plot_main)){
        p <- Seurat::LabelClusters(
          plot = Seurat::DimPlot(
          object = eval(parse(text = fconfig$object))[, scells],
          cols = v2cols(levels(ddf$tmp), fconfig$cols),
          reduction = "umap", group.by = i
        ), id = i, repel = TRUE, parse = TRUE) + theme(legend.position = "none")
        pdf(paste0(fname0, ".pdf"))
        print(p); graphics.off()
        pdf(paste0(fname0, "_blank.pdf"))
        print(plot_blank(p)); graphics.off()
      }
      for(j in fconfig$facet){
        cat("  *", j, "\n")
        ddf_j <- ddf[!is.na(ddf[, j]), ]
        if(isTRUE(fconfig$sample_even))
          ddf_j <- ddf_j[sample_even(ddf_j, j), ]
        ddf_j[, j] <- droplevels(factor(ddf_j[, j]))
        pp <- plot_grids(ddf_j, x = fconfig$redu[1], y = fconfig$redu[2],
          color = i, facet = j, colours = fconfig$cols)
        fname1 <- paste0(fname0, "_", gsub("orig\\.", "", j))
        sizes = ifelse(make_grid(nlevels(ddf_j[, j])) == 1, 9, 18)
        pdf(paste0(fname1, ".pdf"), width = sizes[1], height = sizes[2])
        print(pp); graphics.off()
        pdf(paste0(fname1, "_blank.pdf"), width = sizes[1], height = sizes[2])
        print(plot_blank(pp)); graphics.off()

        if(isTRUE(fconfig$per_facet)){
          for(k in levels(ddf_j[, j])){
            cat("   %", k, "\n")
            scells = filters_subset_df(c(j, k), ddf_j, verbose = TRUE)
            p <- Seurat::LabelClusters(
              plot = Seurat::DimPlot(
              object = eval(parse(text = fconfig$object))[, scells],
              cols = v2cols(levels(ddf$tmp), fconfig$cols),
              reduction = "umap", group.by = i
            ), id = i, repel = TRUE, parse = TRUE) + theme(legend.position = "none")
            fname <- paste0(fname1, "_", k)
            pdf(paste0(fname, ".pdf"))
            print(p); graphics.off()
            pdf(paste0(fname, "_blank.pdf"))
            print(plot_blank(p)); graphics.off()
          }
        }
      }
    }
  }
}

{ cat(red_bold("### Markers: violins/dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("violins", showWarnings = FALSE)
  dir.create("dim_reduction_markers", showWarnings = FALSE)

  # This is missing the previous steps [to generate this results] in this script
  cite_f = paste0(dirname(getwd()), "/a1_final_figures_cd4/cite_seq_normalized.csv")
  cite_df <- readfile(cite_f, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(cite_df) = paste0("cite_", gsub("_Integrin_E", "", colnames(cite_df)))
  sc_cd4resting@meta.data = joindf(sc_cd4resting@meta.data, cite_df)
  signatures_f <- paste0(dirname(getwd()), "/final_figures_cd4/signatures_sc_ayearoftrmness_modulescore/signatures.csv")
  signatures_df <- read.csv(signatures_f, row.names = 1)
  sc_cd4resting@meta.data = joindf(sc_cd4resting@meta.data, signatures_df)
  signatures_f <- paste0(dirname(getwd()), "/final_figures_cd4/signatures_sc_funcpathways_modulescore/signatures.csv")
  signatures_df <- read.csv(signatures_f, row.names = 1)
  sc_cd4resting@meta.data = joindf(sc_cd4resting@meta.data, signatures_df)
  signatures_f <- paste0(dirname(getwd()), "/final_figures_cd4/signatures_sc_stim_ayearoftrmness_modulescore/signatures.csv")
  signatures_df <- read.csv(signatures_f, row.names = 1)
  sc_cd4stim@meta.data = joindf(sc_cd4stim@meta.data, signatures_df)
  fconfigs = list(
    f1g = list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting/",
      features = c("ITGAE", "CD69", "HOPX", "ZNF683"), sufix = "C0n1n2_",
      axis_x = list(name = "cluster", order = c("0", "1", "2"))),
    f2d = list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting_f3bf/",
      features = c("GZMA", "GZMB", "GZMH", "FASLG", "HOPX", "ZNF683",
        "HLA-DRB1", "HLA-DRB5", "HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-DQB1",
        "IFNG", "TNF", "TNFSF14", "CCL4", "CKLF", "CCL5"
      ), axis_x = list(name = "cluster", order = c("1", "0"))),
    list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting_cl0_sex_disease/",
      features = c("HLA-DRB1", "HLA-DRB5", "HLA-DRA", "HLA-DPB1", "HLA-DPA1",
        "HLA-DQB1", "IFNG", "TNF", "TNFSF14", "CCL4", "CKLF", "CCL5"),
      axis_x = list(name = c("cluster", "orig.Sex", "orig.asthma"),
        order = c("0_Female_MA", "0_Female_SA", "0_Male_MA", "0_Male_SA"))),
    list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting_cl0to5_disease/",
      features = c("CREM", "DUSP1", "DUSP2", "DUSP4", "TNFAIP3", "FKBP5", "DDIT4"),
      axis_x = list(name = c("cluster", "orig.asthma"),
      order = paste0(0:5, rep(c("_MA", "_SA"), each = 6)))),
    list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting_cl4_disease/",
      features = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TNFRSF18"),
      axis_x = list(name = c("cluster", "orig.asthma"),
        order = c("4_MA", "4_SA"))),
    list(object = "sc_cd4resting", result_id = "violins/sc_cd4resting_cl3_disease/",
      features = c("JUNB", "JUN", "FOS", "FOSB"),
      axis_x = list(name = c("cluster", "orig.asthma"), order = c("3_MA", "3_SA"))),
    list(object = "sc_cd4merged", result_id = "violins/sc_cd4merged_gzmb_cells/",
      features = c("CCL3", "CCL4", "IL13", "IL4", "IFNG", "IL17A", "IL21", "TNF",
        "CCL5", "IL23A", "CSF2", "IL32", "CCL20", "IL5", "XCL1", "XCL2", "TNFSF14",
        "TNFSF12", "TNFSF9"), axis_x = list(name = "gzmb_stim_disease",
      order = c("GZMBp-US-MA", "GZMBp-US-SA", "GZMBp-ST-MA", "GZMBp-ST-SA")))
  )
  fconfigs = c(fconfigs,
    lapply(X = fconfigs[1:3], FUN = function(x){
       x$features = c("cite_CD103", "cite_CD69"); x$result_id = "violins/sc_cd4resting_cite/"; x
    }))
  tmp = list(list("_TCM", "2", genes = c("S1PR1", "CCR7", "TCF7", "ICAM2", "SELL")),
    list("_TREG", "3", genes = c("FOXP3", "IL2RA", "IKZF2")),
    list("_TFH", "4", genes = c("MAF", "PDCD1", "CXCL13")),
    list("_TIFNR", "5", genes = c("IFIT1", "IFIT3", "OAS1")),
    list("_TCYCLE", "7", genes = c("TOP2A", "MKI67")),
    list("_TCTL", "8", genes = c("GZMA", "GZMB", "GZMH", "GNLY", "PRF1")))
  fconfigs = c(fconfigs, lapply(X = tmp, FUN = function(x){
    list(object = "sc_cd4resting", result_id = paste0("violins/sc_cd4resting", x[[1]], "/"),
      axis_x = list(name = "cluster", order = c(x[[2]], "REST")), features = x[[3]])
  }))
  fconfigs = list(
    list(object = "sc_cd4resting", result_id = paste0("dim_reduction_markers/sc_cd4resting/"),
      features = c("HOMBRINK_TEICHMANN_2_FILT.Score")),
    list(object = "sc_cd4resting", result_id = paste0("dim_reduction_markers/sc_cd4resting_sex_disease/"),
      features = c("GR_138G_GCGS.Score"), facet = "sex_disease"),
    list(object = "sc_cd4stim", result_id = paste0("dim_reduction_markers/sc_cd4stim_disease/"),
      facet = "orig.asthma", features = c("IL23A", "IL32", "XCL1", "XCL2",
        "TNFSF12", "TNFSF9", "HOMBRINK_TEICHMANN_2_FILT.Score")),
    list(object = "sc_cd4stim", result_id = paste0("violins/sc_cd4stim_disease/"),
      axis_x = list(name = "orig.asthma", order = c("MA", "SA")),
      features = c("IL23A", "IL32", "XCL1", "XCL2", "TNFSF12", "TNFSF9", "HOMBRINK_TEICHMANN_2_FILT.Score"))
  )
  # couls_opt_i = unlist(couls_opt, recursive = FALSE)
  # names(couls_opt_i) <- paste0(gsub("_gradient", "", names(couls_opt_i)), "_")
  # couls_opt_i <- couls_opt_i[grepl("brew", names(couls_opt_i))]
  fconfigs = c(fconfigs,
    lapply(X = fconfigs[sapply(fconfigs, function(x) grepl("violin", x$result_id) )],
      FUN = function(x){ x$sufix <- paste0(x$sufix, "chilli_"); x }))
  # [sapply(fconfigs, function(x) grepl("cite", x$result_id) )]
  for (fconfig in fconfigs) {
    if(!exists("sc_cd4merged") && fconfig$object == "sc_cd4merged"){
      sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
      sc_cd4merged@meta.data$gzmb_stim_disease = ident_combine(
        sc_cd4merged@meta.data, c("tag_GZMB", "orig.stim", "orig.asthma"))
    }
    features_df = features_find(
      features = fconfig$features,
      universe = c(rownames(eval(parse(text = fconfig$object))),
        colnames(eval(parse(text = fconfig$object))@meta.data)),
      verbose = TRUE)
    tvar <- !dir.exists(fconfig$result_id) && grepl("\\/$", fconfig$result_id)
    if(tvar) dir.create(fconfig$result_id)
    fname <- paste0(fconfig$result_id, "_features.csv")
    cat(cyan(c("- ", fconfig$result_id, fconfig$sufix, "\n")), sep = "")
    write.csv(features_df, file = fname, row.names = FALSE)

    myfeatures <- unique(features_df[features_df$found != "", ]$found)
    ddf <- Seurat::FetchData(eval(parse(text = fconfig$object)),
      vars = unique(c(redu[[1]], myfeatures, fconfig$axis_x$name, fconfig$facet)))
    ddf = fig_set_identities(ddf, fconfig$axis_x, verbose = TRUE)
    ddf <- ddf[!is.na(ddf$Identity), ] # table(ddf[, c("cluster", "Identity")])
    if(isTRUE(fconfig$sample_even) && !is.null(fconfig$facet))
      ddf <- ddf[sample_even(ddf, fconfig$facet, verbose = TRUE), ]
    # for(grad_i in names(couls_opt_i)){
    for(g in myfeatures){
      cat(" *", g, "\n");
      # fname <- paste0(fconfig$result_id, "zcolortest_", grad_i, fconfig$sufix, g)
      fname <- paste0(fconfig$result_id, fconfig$sufix, g)
      # if(file.exists(paste0(fname, ".pdf"))) next
      sizes = c(7, 7)
      if(grepl("violin", fname)){
        p = violin(ddf[!is.na(ddf[, g]), ], "Identity", g, colour_by = "pct",
          # couls = couls_opt_i[[grad_i]],
          chilli = grepl("chilli", fname)
        ) + labs(x = NULL, y = "Seurat Normalised") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        # p <- plot_rm_layer(p, "Box")
      }else{
        tmp = grepl("score", fname, ignore.case = TRUE)
        ddf$Feature = if(tmp) ifelse(ddf[, g] < 0, 0, ddf[, g]) else df[, g]
        tvar <- ifelse(tmp, "Module\nScore", "Seurat\nNormalized")
        p <- ggplot(data = ddf[!is.na(ddf$Feature), ],
          mapping = aes_string(x = redu[[1]][1], y = redu[[1]][2], color = "Feature")
        ) + geom_point(size = 0.2) +
          labs(color = tvar) +
          scale_colour_gradientn(colours = couls_opt$red_gradient$strong_white)
        if(!is.null(fconfig$facet)){
          p <- plot_facet_wrap(p, fconfig$facet)
          sizes = plot_size(make_grid(nlevels(factor(ddf[, fconfig$facet]))))
        }
      }
      pdf(paste0(fname, ".pdf"), width = sizes[1], height = sizes[2])
      print(p)
      graphics.off()
      pdf(paste0(fname, "_blank.pdf"), width = sizes[1], height = sizes[2])
      print(plot_blank(p))
      graphics.off()
    }#}
  }
}

{ cat(red_bold("### Crater plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/ciro/scripts/handy_functions/devel/plots_crater.R")
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("craters", showWarnings = FALSE)
  dgea_dir = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/cluster"
  edata = as.matrix(expm1(sc_cd4resting@assays$RNA@data))
  mygenes = grep("XIST|RPS4Y1", rownames(edata), value = TRUE, invert = TRUE)
  mdata = sc_cd4resting@meta.data
  fconfigs = lapply(X = levels(sc_cd4resting$cluster),
    FUN = function(x){
      list(fnames = c(
        "Male_MA_vs_SA" = paste0(dgea_dir, x, "/Male_SAvsMale_MA/results_Male_SAvsMale_MA_mastlog2cpm.csv"),
        "Female_MA_vs_SA" = paste0(dgea_dir, x, "/Female_SAvsFemale_MA/results_Female_SAvsFemale_MA_mastlog2cpm.csv")
      ), result_id = paste0("craters/colsizeaxes_cluster", x, "_"), selectss = list(c(sc_cd4resting_clust, x)),
        columns = c("sex_disease"), limits_col = c(0, 9), limits_size = c(0, 120), plot_squared = 6.5,
        highlight_genes = c("CREM", "DUSP1", "DUSP2", "DUSP4", "TNFAIP3", "DDIT4", "FKBP5"))
  })
  # if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
  # edata = expm1(sc_cd4merged@assays$RNA@data) # MEMORY intensive
  # mygenes = grep("XIST|RPS4Y1", rownames(edata), value = TRUE, invert = TRUE)
  # mdata = sc_cd4merged@meta.data
  # rm(sc_cd4merged); gc(); edata = as.matrix(edata)
  fconfigs = list(list(
      fnames = c(
        "Male_MA_vs_SA" = "/home/ciro/ad_hoc/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex_disease/SA_MalevsMA_Male/results_SA_MalevsMA_Male_mastlog2cpm.csv",
        "Female_MA_vs_SA" = "/home/ciro/ad_hoc/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex_disease/SA_FemalevsMA_Female/results_SA_FemalevsMA_Female_mastlog2cpm.csv"
      ), result_id = "craters/resting_", columns = c("sex_disease"), plot_squared = TRUE,
      highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL")
  ))
  degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))
  # sapply(fconfigs, function(x) file.exists(x$fnames) )

  for(fconfig in fconfigs){
    cat(cyan(c("- ", fconfig$result_id, "\n")), sep = "")
    for(fc in as.character(c(0.25))){
      cat(" * FC:", fc, "\n");
      void <- crater_plot(
        tests_list = fconfig$fnames,
        edataf = edata,
        annotf = mdata,
        sample_filter = fconfig$selectss,
        gene_filter = degfilt,
        feature_subset = mygenes,
        topgenes = c("top10", fconfig$highlight_genes),
        lfcthresh = fc,
        # column4stats = fconfig$columns,
        outputname = fconfig$result_id,
        plot_interactive = TRUE,
        plot_squared = fconfig$plot_squared,
        limits_col = fconfig$limits_col,
        limits_size = fconfig$limits_size
      )
    }
  }
}

{ cat(red_bold("### GSEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("gsea", showWarnings = FALSE)
  tvar <- as.character(sc_cd4resting@meta.data$cluster)
  sc_cd4resting@meta.data$TRMclusters = ifelse(tvar %in% c("0", "1"), tvar, "nonTRM")
  sc_cd4resting@meta.data$TRM0n1 = ifelse(tvar %in% c("0", "1"), "TRM", "nonTRM")
  sc_cd4resting@meta.data$cluster_disease = as.character(ident_combine(
    sc_cd4resting@meta.data, c("cluster", "orig.asthma"), "_"))
  tvar <- sc_cd4resting@meta.data$cluster_disease
  tvar <- ifelse(sc_cd4resting@meta.data$TRMclusters %in% c("0", "1"), tvar, "nonTRM")
  sc_cd4resting@meta.data$trm_disease = tvar
  sc_cd4resting@meta.data$trm_ma = ifelse(!tvar %in% c("0_MA", "nonTRM"), "TRMnon0_MA", tvar)
  sc_cd4resting@meta.data$trm_sa = ifelse(!tvar %in% c("0_SA", "nonTRM"), "TRMnon0_SA", tvar)
  sc_cd4resting@meta.data$cluster_sex = as.character(ident_combine(
    sc_cd4resting@meta.data, c("cluster", "orig.Sex"), "_"))
  sc_cd4resting@meta.data$cluster_sex_disease = as.character(ident_combine(
    sc_cd4resting@meta.data, c("cluster", "orig.Sex", "orig.asthma"), "_"))
  table(sc_cd4resting@meta.data[, c("sex_disease", "orig.Sex")])
  fconfigs = list(
    list(
      result_id = "gsea/sc_cd4resting_clusters/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = "cluster",
      list_in = "TRM.si|TCM.c|TREG.s|TFH.s|Type.I|Cell.cyc|Cytot.*genes$",
      list_order = c(
        "TRM.Oja", "TCM.circulating.Oja", "TREG.DICE", "TFH.Locci",
        "Type.I.and.II.IFN.signaling.Seumois", "Cell.cycle.Best", "Cytotoxic.Patil"
      ))
  )
  fconfigs = list(
    list(
      result_id = "gsea/sc_cd4resting_sex_disease/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("Male-SA", "Female-SA"), group2 = c("Male-MA", "Female-MA"),
        column = "sex_disease", stringsAsFactors = FALSE
      ), list_in = "cAMP|Glucocorticoid"),
    list(
      result_id = "gsea/sc_cd4resting_TRMclusters/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("0", "1", "TRM"), group2 = c("nonTRM", "nonTRM", "nonTRM"),
        column = c("TRMclusters", "TRMclusters", "TRM0n1"), stringsAsFactors = FALSE
      ), list_in = "TRM.Oja"),
    list(
      result_id = "gsea/sc_cd4resting_trm_disease/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("0_MA", "0_SA", "0_SA"),
        group2 = c("TRMnon0_MA", "TRMnon0_SA", "0_MA"),
        column = c("trm_ma", "trm_sa", "trm_disease"), stringsAsFactors = FALSE
      ), list_in = "TCR"),
    list(
      result_id = "gsea/sc_cd4resting_cluster_sex/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("0_Male"), group2 = c("0_Female"),
        column = c("cluster_sex"), stringsAsFactors = FALSE
      ), list_in = "TCR"),
    list(
      result_id = "gsea/sc_cd4resting_cluster_sex_disease/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("0_Male_SA"), group2 = c("0_Female_SA"),
        column = c("cluster_sex_disease"), stringsAsFactors = FALSE
      ), list_in = "Patil"),
    list(
      result_id = "gsea/sc_cd4resting_cluster_disease/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("0_MA", "0_SA"), group2 = c("1_MA", "1_SA"),
        column = c("cluster_disease"), stringsAsFactors = FALSE
      ), list_in = "TCR"),
    list(
      result_id = "gsea/sc_cd4resting_disease/", object = "sc_cd4resting",
      file = "/home/ciro/asthma_airways/info/GSEA_supp_table_8_22JUNE21.csv",
      comparisons = data.frame(
        group1 = c("SA"), group2 = c("MA"),
        column = "orig.asthma", stringsAsFactors = FALSE
      ), list_in = "cAMP|Glucocorticoid")
  )

  for (fconfig in fconfigs) {
    signatures_files <- fconfig$file[file.exists(fconfig$file)]
    signatures_list = gsea_process_list(
      lapply(signatures_files, readfile, stringsAsFactors = FALSE),
      include = fconfig$list_in, exclude = fconfig$list_out)
    if(any(grepl("8_22JUNE21", signatures_files))){
      tmp <- unname(gsub(" .*", "", sapply(signatures_list, head, 1)))
      tvar <- paste0(gsub("\\.{2,}", ".", gsub("signature|genes|Bulk.", "", names(signatures_list))), tmp)
      names(signatures_list) <- tvar
    }
    str(signatures_list)
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
    write.csv(ddf, file = paste0(fconfig$result_id, "_features_lists.csv"))

    gsea_results <- gsea_matrix(
      mat = expm1(eval(parse(text = paste0(fconfig$object, "@assays$RNA@data")))),
      groups = fconfig$comparisons,
      metadata = eval(parse(text = paste0(fconfig$object, "@meta.data"))),
      gsea_list = signatures_list,
      path = fconfig$result_id,
      verbose = TRUE
    )
    for(i in c(0.05)){
      for(j in c(1.5)){
        x <- try(gsea_plot_summary(
          tests_list = gsea_results,
          path = fconfig$result_id,
          padjthr = i, nesthr = j,
          axes = list(rows = fconfig$list_order), axes_cluster = FALSE
        ), silent = TRUE)
      }
    }
  }
}

{ cat(red_bold("### Markers: contour ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("scatter_contour", showWarnings = FALSE)
  fconfigs = list(
    list(object = "sc_cd4resting", result_id = "scatter_contour/sc_cd4resting_GZMBp/",
      features = list(x = c("CCL3", "IFNG"), y = c("IFNG", "IL13", "IL17A", "TNF")),
      filters = c("tag_GZMB", "GZMBp")),
    list(object = "sc_cd4stim", result_id = "scatter_contour/sc_cd4stim_GZMBp/",
      features = list(x = c("CCL3", "IFNG"), y = c("IFNG", "IL13", "IL17A", "TNF")),
      filters = c("tag_GZMB", "GZMBp"))
  )

  for (fconfig in fconfigs) {
    cat(cyan(c("- ", fconfig$result_id, fconfig$sufix, "\n")), sep = "")
    tvar <- !dir.exists(fconfig$result_id) && grepl("\\/$", fconfig$result_id)
    if(tvar) dir.create(fconfig$result_id)
    scells = filters_subset_df(fconfig$filters,
      eval(parse(text = fconfig$object))@meta.data, verbose = TRUE)
    ddf <- Seurat::FetchData(eval(parse(text = fconfig$object)),
      vars = unique(c(fconfig$redu, unlist(fconfig$features), fconfig$facet)))
    ddf <- ddf[scells, ]
    if(isTRUE(fconfig$sample_even) && !is.null(fconfig$facet))
      ddf <- ddf[sample_even(ddf, fconfig$facet, verbose = TRUE), ]
    myfeatures = data.frame(t(data.frame(
      x = rep(fconfig$features$x, length(fconfig$features$y)),
      y = rep(fconfig$features$y, each = length(fconfig$features$x))
    )), stringsAsFactors = FALSE)
    mytab = scatter_summary_df(ddf, myfeatures)
    fname0 = paste0(fconfig$result_id, fconfig$sufix)
    fname <- paste0(fname0, "_features_stats.csv")
    write.csv(mytab, file = fname, row.names = FALSE)
    for(g in myfeatures){
      cat(" *", g, "\n");
      fname <- paste0(fname0, paste0(g, collapse = "_"))
      if(file.exists(paste0(fname, ".pdf"))) next
      sizes = c(7, 7)
      p <- scatter_contour(ddf, g[1], g[2]) +
        viridis::scale_color_viridis(option = "magma")
      if(!is.null(fconfig$facet)){
        p <- plot_facet_wrap(p, fconfig$facet); sizes[1] = 14
        tmp = nlevels(factor(ddf[, fconfig$facet]))
        sizes = ifelse(make_grid(tmp) == 1, 7, 14)
      }
      pdf(paste0(fname, ".pdf"), width = sizes[1], height = sizes[2])
      print(p)
      graphics.off()
      pdf(paste0(fname, "_blank.pdf"), width = sizes[1], height = sizes[2])
      print(plot_blank(p))
      graphics.off()
    }
  }
}

{ cat(red_bold("### Metadata for HLA comparisons ### %%%%%%%%%%%%%%%%%%%%%%\n"))
  mdata = sc_cd4resting@meta.data
  scells = filters_subset_df(c("cluster", "0", "1"), mdata, v = T)
  mdata = mdata[scells, grepl("asthma$|cluster$|Sex", colnames(mdata))]
  sfeatures = c("HLA-DRB1", "HLA-DRB5", "HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-DQB1")
  tmp = t(as.matrix(sc_cd4resting@assays$RNA@counts[sfeatures, scells]))
  hla_cells = rowSums(tmp > 0) > 0
  mdata$HLA = ifelse(hla_cells[rownames(mdata)], "HLAp", "HLAn")
  saveRDS(mdata, file = "metadata_hla.rds")
}

{ cat(red_bold("### Summary stats for stim_disease ### %%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  dir.create("stats_summary", showWarnings = FALSE)
  if(!exists("sc_cd4merged")){
    sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
    sc_cd4merged@meta.data$stim_disease = as.character(ident_combine(
      sc_cd4merged@meta.data, c("orig.stim", "orig.asthma"), sep = "_"))
    sc_cd4merged@meta.data$gzmb_stim_disease = ident_combine(
      sc_cd4merged@meta.data, c("tag_GZMB", "orig.stim", "orig.asthma"))
  }
  fconfigs = list(
    list(object = "sc_cd4merged", result_id = "stats_summary/sc_cd4merged_stim_disease.csv",
      features = c("IL21", "CCL4"), axis_x = "stim_disease"),
    list(object = "sc_cd4merged", result_id = "violins/sc_cd4merged_gzmb_cells/_stats.csv",
      features = c("CCL3", "CCL4", "IL13", "IL4", "IFNG", "IL17A", "IL21", "TNF",
        "CCL5", "IL23A", "CSF2", "IL32", "CCL20", "IL5", "XCL1", "XCL2", "TNFSF14",
        "TNFSF12", "TNFSF9"), axis_x = "gzmb_stim_disease")
  )

  for(fconfig in fconfigs){
    features = show_found(fconfig$features, rownames(sc_cd4merged), verbose = TRUE)
    ddf <- stats_summary_table(
      mat = expm1(sc_cd4merged@assays$RNA@data[features, ]),
      group = setNames(sc_cd4merged@meta.data[, fconfig$axis_x], rownames(sc_cd4merged@meta.data)),
      moments = c("p", "mn", "md"),
      verbose = TRUE
    ); write.csv(ddf, file = fconfig$result_id)
  }

}

{ cat(red_bold("### Radar plots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  dir.create("radars", showWarnings = FALSE)
  gsea_results = readRDS("../final_figures_cd4/gsea_spec_sc_cd4resting/tests_RNA_snn_res.0.4_0_18lists.rds")
  gsea_summary_list <- gsea_summary(gsea_results, path = dirname(result_id), cache = FALSE)
  mysum_radar <- gsea_summary_list[[1]]
  ddf <- mysum_radar[1, colnames(mysum_radar) != "0vs1", drop = FALSE]
  ddf[is.na(ddf)] <- 0

  p = plot_radar(ddf, cols = "red") + #facet_wrap(~feature) +
    theme(legend.position = "none")
  pdf(paste0("radars/apoptosis.pdf"))
  print(p)
  dev.off()
  # p$layers <- p$layers[-5] # only when threshold given
  pdf(paste0("radars/apoptosis_blank.pdf"))
  print(plot_blank(p))
  dev.off()
}
