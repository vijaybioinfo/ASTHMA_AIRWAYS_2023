## Authors: Francisco Emmanuel Castaneda Castro & Ciro Ramirez-Suastegui 
### This code will do the specified figures on the CD4 BAL paper, bulk and scRNA seq
########################################
####     Rebuttal CD4_BAL          #####
########################################
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
library(GEOquery)
# library(hdf5r)
library(Seurat)
library(dplyr)

resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # filters_columns is.file.finished show_commas (cluster_reports)
  "/home/ciro/scripts/clustering/R/plotting.R", # cluster_reports
  "/home/ciro/scripts/clustering/R/utilities.R", # get_top_n
  "/home/ciro/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/handy_functions/devel/plots.R", # plot_pct getlegend mytheme
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
  "/home/fcastaneda/bin/clustering/R/plotting.R" #plots of clustering pipeline
)
for(i in resources){ source(i) }
source("/home/ciro/scripts/figease/figease.R")
#source("/home/ciro/scripts/figease/figease.R")
library(Seurat)
library(ggplot2)
library(tidyverse)

  outdir<-"/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal_Figures_Final"

{ cat(redb("### Global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  #######CD4_BAL

  source("/home/ciro/asthma_airways/scripts/final_figures_cd4/global.R")

  dir.create(outdir)
  setwd(outdir)

  sc_cd4resting@meta.data = joindf(sc_cd4resting@meta.data,
      as.data.frame(sc_cd4resting@reductions$umap@cell.embeddings))

  sc_cd4stim@meta.data = joindf(sc_cd4stim@meta.data,
      as.data.frame(sc_cd4stim@reductions$umap@cell.embeddings))

  dir.create("resting")
  setwd("resting")

  #definition of traetment for dgea taking Tratmen as a covariate
    a<-readRDS("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/Herrera_GSENNNNNN_annotation.rds")
    a$Treatment2<-ifelse(a$Treatment == "", "None", a$Treatment)
    table(a$Treatment, a$Treatment2)
    saveRDS(a, "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/Herrera_GSENNNNNN_annotation_to_dgea.rds")

    ## To DGE Treatment vs Responders
    sc_cd4resting$Treatment2<-ifelse(sc_cd4resting$Treatment == "", "None", sc_cd4resting$Treatment)
    table(sc_cd4resting$Treatment, sc_cd4resting$Treatment2)

    dgea_treatment<-c("NIHW_00060"='Treatment_Not_responders' ,  "NIHW_00352"='Treatment_Responders' ,  "NIHW_00011"='Treatment_Responders' ,  "NIHW_00367"='No_treatment' ,  "NIHW_00334"='Treatment_Responders' , "NIHW_00386"='Treatment_Not_responders		' , "NIHW_00070"='No_treatment' ,  "NIHW_00150"='No_treatment' ,  "NIHW_00304"='Treatment_Responders' ,  "NIHW_00361"='Treatment_Responders' ,  "NIHW_00356"='No_treatment' , "NIHW_00362"='Treatment_Responders' ,  "NIHW_00419"='Treatment_Not_responders' ,  "NIHW_00364"='No_treatment' ,  "NIHW_00398"='Treatment_Responders' ,  "NIHW_00046"='No_treatment')

    sc_cd4resting$to_dgea_Treatment<-dgea_treatment[sc_cd4resting$orig.donor]
    table(sc_cd4resting$orig.donor, sc_cd4resting$to_dgea_Treatment)

    saveRDS(sc_cd4resting@meta.data, "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/Herrera_GSENNNNNN_annotation_to_dgeaTreatment_Responders.rds")

}

### Figure 1

{ cat(redb("### Dotplot CHECK - Fig 1D ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
      dir.create("dotplots")
      source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/figease.R") #no margin circle color
      # source("/home/ciro/sc ipts/figease/figease.R")
      sc_cd4resting$orig.sex_asthma<-paste0(sc_cd4resting$orig.Sex, "_", sc_cd4resting$orig.asthma)

      fconfigs = list(
        list(result_id = "./dotplots/", sufix = "19_Main_figure_19aprl2023",
        object = "sc_cd4resting",
          edata = "sc_cd4resting@assays$RNA@data", metadata = "sc_cd4resting@meta.data",
          size = c(7,10),
          axis_x = list(
          col = "celltype_subset", order = c("TRMDP", "TRMSP", "TCM", "Treg", "TFH", "ThIFNr", "CellCycle", "CTLs")
        # col = "celltype", order = c("Basal", "Club", "Ciliated", "Fibroblast", "Ionocyte", "Tcell", "Bcell", "DC", "Neutrophils", "Mast", "Cycling")
          ),
          features =  rev(c("ITGAE","CD69","ITGA1","HOPX","AMICA1","ZNF683","S1PR1","CCR7","TCF7","ICAM2","SELL","FOXP3","IL2RA","IKZF2","MAF","PDCD1","CXCL13","IFIT1","IFIT3","OAS1","TOP2A","MKI67","GZMA","GZMB","GZMH","GNLY","PRF1"))
        )
      )

      source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/figease.R") #no margin circle color
      pp_curtains = fig_plot_curtain(fconfigs[1], dot.scale = 10, verbose = 2)
      pp_curtains = fig_plot_curtain(fconfigs[1], verbose = 2,return_plot=TRUE, col.min = 0)
      pp_curtains = fig_plot_curtain(fconfigs[3], verbose = 2,return_plot=TRUE, col.min = 0)

### Figure 2

{ cat(redb("### Histograms - Fig 2B ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
      library(reshape2)
      dir.create("histograms_donors")

      colors_sc_cd4resting_ident<-c("#99141B", "#FFA303", "#FFE610", "#21B84C", "#14B9FA", "#A6B531", "#A3A3A3", "#C944C8")
      names(colors_sc_cd4resting_ident)<-names(sc_cd4resting_ident$colours[1:8])

      sc_cd4resting$order_donors<-paste0(sc_cd4resting$orig.donor, "_", sc_cd4resting$orig.Sex, "_", sc_cd4resting$orig.asthma)

      table_prop_donor<-as.data.frame.matrix(table(sc_cd4resting$order_donors, sc_cd4resting$celltype_subset))

      patab<-table_prop_donor/rowSums(table_prop_donor)

      patab$donor<-rownames(patab)

      a<-melt(patab, id.vars="donor")
      write.csv(patab, "./histograms_donors/proportions_per_cluser_percentage.csv")
      write.csv(table_prop_donor, "./histograms_donors/proportions_per_cluser_raw_counts.csv")

        order1<-c("TRMDP", "TRMSP", "TCM", "Treg", "TFH", "ThIFNr", "CellCycle", "CTLs")
        # barpl_2$ARTE<-gsub("HC0[0-9]_", "", barpl_2$Groups)
        # barpl_2$Donor<-gsub("_[A-Z]{3}", "", barpl_2$Groups)

        order_donor<-c("NIHMA_246_Female_MA", "NIHMA_254_Female_MA",  "NIHMA_244_Female_MA", "NIHW_00356_Female_SA", "NIHW_00150_Female_SA", "NIHW_00361_Female_SA", "NIHW_00334_Female_SA", "NIHW_00398_Female_SA", "NIHW_00352_Female_SA", "NIHW_00362_Female_SA", "NIHW_00011_Female_SA", "NIHMA_245_Male_MA", "NIHMA_228_Male_MA", "NIHMA_206_Male_MA",  "NIHMA_226_Male_MA", "NIHMA_225_Male_MA", "NIHMA_249_Male_MA", "NIHW_00367_Male_SA", "NIHW_00364_Male_SA", "NIHW_00419_Male_SA", "NIHW_00386_Male_SA", "NIHW_00070_Male_SA", "NIHW_00060_Male_SA", "NIHW_00046_Male_SA", "NIHW_00304_Male_SA")

        p<-ggplot(a, aes(x=factor(donor, order_donor) , y=value, fill=variable)) +
        geom_bar(position="stack", stat="identity") +
        scale_fill_manual(values=colors_sc_cd4resting_ident) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")


           pdf("./histograms_donors/histograms_clusters_donors_order_TRM.pdf", width=6)
                  print(p)
            dev.off()


    }

### Figure 3

{ cat(redb("### Blanks volcano TRMDP vs TRMSP Sex covariate -Fig3A ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("volcano")
  setwd("volcano")
  dir.create("CD103pvsCD103n_sexcov_plot/")
  #I ran the documented line a

  # nd change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_cov_sex/0vs1/results_0vs1_mastlog2cpm.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "0"
    group2 = "1"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0(x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

    cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  # dir.create("volcano")

  library("plotly")
  pdf("CD103pvsCD103n_sexcov_plot/CD103pvsCD103n_sexcov_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-3.2,3.2),
      ylims=c(0,320),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

  p_html <- volplot(
    resdf,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = "padj",
    lfctype = "log2FoldChange",
    col_feature = means,
    size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
    gene_name = "gene",
    group = "degs",
    check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))), #Delate GZMA??
    return_plot = TRUE,
    clipp = 4,
    xlims=c(-2.5,2.5),
    ylims=c(0,320),
    verbose = verbose ,
    interact= TRUE
  )

  htmlwidgets::saveWidget(
          as_widget(p_html), 'CD103pvsCD103n_sexcov_plot/CD103pvsCD103n_sexcov_blank.html')

  pdf("CD103pvsCD103n_sexcov_plot/CD103pvsCD103n_sexcov.pdf") #Change name for CD8
    volplot_wnames(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(showgenes)),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-3.2,3.2),
      ylims=c(0,320),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

}

setwd(paste0(outdir, "/resting"))

{ cat(redb("### Heatmap pseudo-bulk for TRMSP vs TRMDP Male and Female - Fig3B ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  source("/mnt/BioHome/fcastaneda/bin/dgea/R/group_specific.R")

  library(tidyverse)
     fnames <- rev(list.files(
      path = "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/group_specific/info/sex", pattern = "_mastlog2cpm_results.csv", full.names = TRUE))

      names(fnames) <- rev(c("Female_TRMDPvsFemale_TRMSP", "Male_TRMDPvsMale_TRMSP")) # needs to be named

     results = lapply(X = fnames, read.csv, row.names = 1, check.names = FALSE)

   # Filter unwanted genes
     sc_cd4resting$cluster_sex_disease<-paste0(sc_cd4resting$celltype_subset, "_", sc_cd4resting$orig.asthma, "_", sc_cd4resting$orig.Sex)

     sc_cd4resting$Sex_cluster<-paste0(sc_cd4resting$orig.Sex, "_", sc_cd4resting$celltype_subset)

     table(sc_cd4resting$Sex_cluster)

     sc_cd4resting$keep<-ifelse(sc_cd4resting$celltype_subset %in% c("TRMSP", "TRMDP"), "keep", "no_k")

     table(sc_cd4resting$keep, sc_cd4resting$celltype_subset)

    sc_cd4resting_2<-SplitObject(sc_cd4resting, "keep")
    sc_cd4resting_keep<-sc_cd4resting_2$keep

      setwd(paste0(outdir, "/resting/"))
      library(data.table)
     mdata = sc_cd4resting_keep@meta.data
     edata =  expm1(sc_cd4resting_keep@assays$RNA@data)


        results_list = results
        path = paste0(outdir, "/resting/SEX_23May")
        colours = NULL
        return_report = TRUE
        optimize = FALSE
        verbose = TRUE
        # to calculate the stats
        mdata = mdata
        edata = edata
        hname = "Sex_cluster"
        redundant = TRUE
        sep="&"
        fcthr=0.25 #,
        padjthr = 0.05
        vs = "vs"

        if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%% GSpec report %%%%%%%%%%%%%%%%%%%%\n")
        if(!is.null(path) && verbose) cat("Report at:\n", path, "\n")
        the_report = list()
        results_group_specific = group_specific_features(
          results = results_list,  redundant = TRUE,  sep="&", fcthr=0.25, padjthr = 0.05, verbose = TRUE)

            table(is.na(results_group_specific$group_specific$'log2FoldChange(Male_TRMDPvsMale_TRMSP)'))

            table(is.na(results_group_specific$group_specific$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)'))

            results_group_specific$group_specific$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)'[is.na(results_group_specific$group_specific$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)')] <- 0

            table(is.na(results_group_specific$group_specific$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)'))

          # check _summary
        # if(isTRUE(list(...)$check)) return(results_group_specific)
        if(verbose > 1) str(results_group_specific)
        summ_df = results_group_specific$group_specific
        if(is.null(results_group_specific$stats)){
          if (verbose) cat("Summary statistics fetching or calculation\n")
          stat_df = data.frame(
            row.names = unique(unlist(lapply(results_list, rownames))))
          stat_df$gene_name = rownames(stat_df)
          tmp <- paste0(levels(summ_df$group),
            rep(moments[c("mn", "p")], each = nlevels(summ_df$group)))
          for(i in tmp){
            for(j in 1:length(results_list)){
              y <- which(grepl(paste0("^", i), colnames(results_list[[j]])))
              if(any(grepl(paste0("^", i), colnames(stat_df)))) next
              if(length(y))
                stat_df = joindf(stat_df, results_list[[j]][, y, drop = FALSE])
            }
          }
          if(!is.null(mdata) && !is.null(edata)){
            results_group_specific$stats = stats_summary_table(
              mat = edata,
              groups = setNames(mdata[, hname], rownames(mdata)),
              rname = rownames(results_group_specific[[1]]),
              moments = c("mn", "p"),
              verbose = verbose
            )
          }else if(ncol(stat_df) > 1){
            if(verbose) cat("------- Stat summary -------\nTaken from results\n")
            if(verbose) warning("Some measurements may be missing from comparisons\n")
            colnames(stat_df) <- gsub(".eurat.*|CPM.*", "", colnames(stat_df))
            results_group_specific$stats = stat_df[, -1, drop = FALSE]; rm(stat_df)
          }
        }

        if(!is.null(results_group_specific$stats)){
          summ_df = joindf(summ_df, results_group_specific$stats)
        }
        # summ_df <- summ_df[order(summ_df$group), ]
        # groups = levels(summ_df$group)[levels(summ_df$group) %in% colnames(summ_df)]
        # tvar <- colnames(summ_df)[!colnames(summ_df) %in% groups]
        # summ_df <- summ_df[, c(tvar, groups)]
        the_report$summary = summ_df
        # if(grepl("summary", return_report)) return(the_report)
        # if(!isTRUE(return_report))
          write.csv(summ_df, file = paste0(path, "_summary_2.csv"), row.names = FALSE)

        if(!is.null(results_group_specific$stats)){
          if(verbose) cat("--- Heatmap\n")
          # annor = summ_df[summ_df[, "group"] != "NDE", "group", drop = FALSE]
          # annor[, "group"] <- droplevels(annor[, "group"])

          annor_comb = summ_df[summ_df[, "combined"] != "NDE", c("combined", "log2FoldChange(Male_TRMDPvsMale_TRMSP)", "log2FoldChange(Female_TRMDPvsFemale_TRMSP)"), drop = FALSE]

          annor_comb$avg_FCs<-rowMeans(annor_comb[,c( 'log2FoldChange(Male_TRMDPvsMale_TRMSP)', 'log2FoldChange(Female_TRMDPvsFemale_TRMSP)')])

          table(is.na(annor_comb$avg_FCs))
          # annor_comb$

          annor_comb[, "combined"] <- droplevels(factor(annor_comb[, "combined"]))
          table(annor_comb[, "combined"])

          annor_comb<-annor_comb[order(annor_comb$avg_FCs),]

          # annor_comb$shared_level <- NULL
          # annor_comb$shared_level2 <- NULL
          annor_comb$'log2FoldChange(Male_TRMDPvsMale_TRMSP)' <- NULL
          annor_comb$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)' <- NULL
          annor_comb$avg_FCs <- NULL

          columns <- grep("_percent", colnames(summ_df), value = TRUE)
          if(length(columns) == 0){ str(summ_df); stop("No stats columns found") }
          tvar <- apply( # keeping features with a % difference > 25
            X = summ_df[rownames(annor_comb), columns],
            MARGIN = 1, FUN = function(x) abs(diff(range(x)))
          )
          labels_row_i <- rownames(annor_comb)
           labels_row_i[!rownames(annor_comb) %in% names(tvar[tvar > 20])] <- ""
          the_report$heatmap = list()
          for(i in c("mn", "p")){
             i<-"mn"
            moments_i = moments[i]
            if(verbose) cat("    *", moments_i, "\n")
            anno_cnames = grepl(paste0(moments_i, "$"), colnames(summ_df))
            mat2plot <- as.matrix(summ_df[rev(rownames(annor_comb)), anno_cnames])
            if(moments_i[[1]] == "_mean") mat2plot <- log2(mat2plot + 1)
            if(1){ # row-wise scale and get the range
              mat2plot <- t(scale(t(mat2plot)))
              topz <- max(c(min(abs(c(range(mat2plot), 2))), 1))
            }else{ topz <- 75 } # this was for percentages...
            mat2plot[mat2plot > topz] <- topz; mat2plot[mat2plot < (-topz)] <- -topz;
            colnames(mat2plot) <- sub(moments_i, "", colnames(mat2plot))
            tvar <- intersect(levels(annor_comb[, "combined"]), colnames(mat2plot))
            mat2plot <- mat2plot[, tvar]


            annoc = data.table(
              Group = c(colnames(mat2plot)[c(2,4,1,3)]))

            annoc<-data.frame(annoc[, c("Sex", "cluster") := tstrsplit(annoc$Group, "_", fixed=TRUE)])

            rownames(annoc) <- c(colnames(mat2plot)[c(2,4,1,3)])
            annoc<-annoc[,c("cluster", "Sex")]

            colors_i = c(as.list(annor_comb), list(Group = colnames(mat2plot)))
            tmp <- v2cols(unname(unlist(colors_i)), colours)
            colors_i <- lapply(colors_i, v2cols, sour = tmp)
            colors_i$cluster<-c('TRMDP'="#99141B", 'TRMSP'="#FFA303")
            # colors_i$disease<-c('MA'="#160cc9", 'SA'="#f55236")
            colors_i$sex<-c('Female'="#d3a0eb", 'Male'="#86bf73")

            # if(!isTRUE(return_report)){
              fname <- paste0(path, "_heatmap", moments_i)
              # if(optimize){
              #   png(paste0(fname, ".png"), width = 1500, height = 1700, res = 250)
              # }else{
            # pdf(paste0(fname, "_legend.pdf"), width = 7, height = 12) #}
            pdf(paste0(fname, "_corrected2_legend.pdf"), width = 7, height = 12)
            # pdf(paste0(fname, "_corrected2.pdf"), width = 7, height = 12)

            # }

            cols <- grDevices::colorRampPalette(rev(c('yellow', 'black', 'blue')))(256)

            # annor = summ_df[summ_df[, "group"] != "NDE", "group", drop = FALSE]
            # annor[, "group"] <- droplevels(annor[, "group"])

            head(mat2plot[,annoc$Group])
            mat2plot<-mat2plot[,rownames(annoc)]
            head(mat2plot)

            pheatmap(
              mat = mat2plot,
              color = cols, border_color = NA, scale = "none",
              cluster_cols = FALSE, cluster_rows = FALSE,
              annotation_row = annor_comb, annotation_col = annoc,
              annotation_colors = colors_i,
              fontsize_row = 10/(4*log10(nrow(mat2plot))),
              annotation_names_col = FALSE, annotation_names_row = FALSE,
              show_colnames = FALSE, show_rownames = TRUE,
              gaps_col = cumsum(table(annoc$Group)),
               # gaps_row = cumsum(table(annor_comb[, "combined"])),
              labels_row = labels_row_i, annotation_legend= TRUE #FALSE
            ); graphics.off()
          }
          if(grepl(paste0("heatmap"), return_report)) return(the_report)
        }
}

{ cat(redb("### GSVA per donor list of genes -Fig3I & Fig 4D&E  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # dir.create("gsva")
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary, gsea_process_list")
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R")

      library(scGSVA)
      library(ggpubr)

      gdf <- readfile("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/info/signatures_table_2.csv", stringsAsFactors = FALSE)

      slist = gsea_process_list(lapply(gdf, function(x) x[-1] ))

      stim_resting_signatures_Sara <- readfile("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/info/signature_list_Sara_CD4_Rebbutal_March_stim.csv", stringsAsFactors = FALSE)

      slist2<-slist[c("X...Genes.upregulated.in.cytotoxic.CD4.T.cells", "TCR.signalling.pathway", "TCR.signaling.genes", "cAMP.signaling.genes")]

      slist2$TH2_Seumois <- stim_resting_signatures_Sara$TH2_Seumois
      slist2$TCR.signalling.pathway <- NULL

      sc_signatures_lists_no_sc<-lapply(slist2, function(x){
        # x<-sc_signatures_lists_without_ribosomals[1]
        genes_2<-as.vector(x)
        genes_2_li<-genes_2[genes_2 %in% rownames(sc_cd4resting@assays$RNA@counts) ]
        return(genes_2_li)
      })


      signature_list_dgea_to_gsva<-melt(sc_signatures_lists_no_sc)

      colnames(signature_list_dgea_to_gsva)<-c("GeneID", "PATH")
      signature_list_dgea_to_gsva$Annot<-signature_list_dgea_to_gsva$PATH

        res<-scgsva(sc_cd4resting,signature_list_dgea_to_gsva, useTerm =FALSE)

        sc_cd4resting@meta.data <- joindf(sc_cd4resting@meta.data, as.data.frame(res@gsva))


        gsva_list2<-sc_cd4resting@meta.data %>% group_by(orig.donor, celltype_subset) %>% summarize(across(TCR.signaling.genes:cAMP.signaling.genes, mean))

        sexes<-distinct(sc_cd4resting@meta.data, orig.donor, orig.Sex, orig.asthma)

        gsva_list_3<-merge(gsva_list2, sexes, by="orig.donor")

        dir.create("gsva")
        write.csv(gsva_list_3, "./gsva/gsva_per_donor_per_cluster_per_donor.csv")

}

{ cat(redb("### Shannon Index -Fig3K %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  ### Shannon Index
  setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting")
  dir.create("shannon_index")
  dir.create("shannon_index/raw")

  bulkTCR<-read.csv("/mnt/BioAdHoc/Groups/vd-vijay/ciro/asthma_pjs/results/mixcr_rnseq/asthma_airways_cd4_A-B_clones.txt", sep="\t")
  bulkTCR$Cell_type<-gsub("-", "_", bulkTCR$Cell_type)
  unstim_bulktcr<- bulkTCR %>% filter(Stim == "Unstim" & !is.na(TRB.clonalSequence))

  #make a table with #count freq cdr3nt cdr3aa v d j why the # at the ??

  lapply(unique(unstim_bulktcr$Cell_type), function(celltype){
    cat(celltype, "\n")
    # celltype<-unique(unstim_bulktcr$Cell_type[1])
    direc<-paste0("shannon_index/raw/", celltype)
    dir.create(direc)
    subset<-unstim_bulktcr %>% filter(unstim_bulktcr$Cell_type == celltype)
    lapply(unique(subset$Study_ID), function(don){
      cat(don, "\n")
        # don<-unique(subset$Study_ID)[1]
        subset_donor<-subset %>% filter(Study_ID == don)
        subset_donor$proportion_clone<-subset_donor$TRB.cloneCount/sum(subset_donor$TRB.cloneCount)

        info_clonotypes<-data.frame(
          count= subset_donor$TRB.cloneCount,
          freq= subset_donor$proportion_clone,
          cdr3nt= subset_donor$TRB.clonalSequence,
          cdr3aa= subset_donor$TRB.aaSeqCDR3,
          v= ifelse(!grepl("TRBV", subset_donor$TRB.bestVHit), NA, subset_donor$TRB.bestVHit),
          d= ifelse(!grepl("TRBJ", subset_donor$TRB.bestJHit), NA, subset_donor$TRB.bestJHit),
          j= ifelse(!grepl("TRBD", subset_donor$TRB.bestDHit), NA, subset_donor$TRB.bestDHit)
        )
        colnames(info_clonotypes)[1]<-"#count"
      write.table(info_clonotypes, file=paste0(direc,"/", don, "_", celltype, "_unstim_cd4.tsv"), sep="\t", row.names = FALSE, quote=FALSE)
    })
    # #file.name      sample.id
    files_to_vdjtools<-data.frame(
      file.name=paste0(getwd(), "/shannon_index/raw/", celltype, "/",list.files(direc)),
      sample.id=gsub("_unstim_cd4.tsv", "",  list.files(direc)))

    colnames(files_to_vdjtools)[1]<-"#file.name"
    write.table(files_to_vdjtools, file=paste0(direc, "_file_to_vdjtools.tsv"), sep="\t", row.names = FALSE, quote=FALSE)
  })
}

{ cat(redb("### Upset plots -Fig3L & Supl Fig5C ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
    bulkTCR<-read.csv("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/sharing_bulk_combinations/all_table_info2.csv")
    dir.create("upset_plots")
    dir.create("upset_plots/resting")

    #to do the compraision donor-wise
    length(unique(bulkTCR$TRB.aaSeqCDR3))
    bulkTCR$TRB.aaSeqCDR3_donor<-paste0(bulkTCR$TRB.aaSeqCDR3, "_", bulkTCR$Study_ID)
    length(unique(bulkTCR$TRB.aaSeqCDR3_donor))

    length(unique(bulkTCR$TRB.aaSeqCDR3_donor)) - length(unique(bulkTCR$TRB.aaSeqCDR3)) #public clonotypes???

    ## All
    sharing_celltypes<-as.data.frame.matrix(table(bulkTCR$TRB.aaSeqCDR3_donor, bulkTCR$Celltype))

    comb<-colnames(sharing_celltypes)
    range(sharing_celltypes)
    comb_logical<-ifelse(sharing_celltypes>=1, 1, 0)
    sharing_celltypes<-as.data.frame.matrix(comb_logical)
    all(comb_logical == sharing_celltypes)
    range(sharing_celltypes)


    sharing_celltypes$TRB.aaSeqCDR3_donor <-rownames(sharing_celltypes)

    clone_sizes<-bulkTCR %>% group_by(TRB.aaSeqCDR3_donor) %>% dplyr::summarize(
      clone_size=sum(TRB.cloneCount),
      donors = paste(unique(unlist(str_split(Study_ID, ";"))), collapse=";"),
      asthma= unique(Disease)
    )

    clone_sizes$donor_num<-str_count(clone_sizes$donors, ";")+1
    table(str_count(clone_sizes$donors, ";"))

    shar_cellt<-merge(sharing_celltypes, clone_sizes, by="TRB.aaSeqCDR3_donor")

    library(ggplot2)
    library(ComplexUpset)

    donors_name<-unique(shar_cellt$donors)
    donors_color<-ifelse(grepl("NIHMA", donors_name), "#160cc9", "#f55236")
    names(donors_color)<-donors_name

    sharing_celltypes<-as.data.frame.matrix(table(bulkTCR$TRB.aaSeqCDR3_donor[bulkTCR$Celltype %in% c("CD4.TRM", "CD4.TRM.like", "CD4.Teff")], bulkTCR$Celltype[bulkTCR$Celltype %in% c("CD4.TRM", "CD4.TRM.like", "CD4.Teff")]))

    comb<-colnames(sharing_celltypes)
    range(sharing_celltypes)
    comb_logical<-ifelse(sharing_celltypes>=1, 1, 0)
    sharing_celltypes<-as.data.frame.matrix(comb_logical)
    all(comb_logical == sharing_celltypes)
    range(sharing_celltypes)


    sharing_celltypes$TRB.aaSeqCDR3_donor <-rownames(sharing_celltypes)

    clone_sizes<-bulkTCR %>% group_by(TRB.aaSeqCDR3_donor) %>% dplyr::summarize(
      clone_size=sum(TRB.cloneCount),
      donors = paste(unique(unlist(str_split(Study_ID, ";"))), collapse=";"),
      asthma= unique(Disease)
    )

    clone_sizes$donor_num<-str_count(clone_sizes$donors, ";")+1
    table(str_count(clone_sizes$donors, ";"))

    shar_cellt<-merge(sharing_celltypes, clone_sizes, by="TRB.aaSeqCDR3_donor")

    library(ggplot2)
    library(ComplexUpset)

    donors_name<-unique(shar_cellt$donors)
    donors_color<-ifelse(grepl("NIHMA", donors_name), "#160cc9", "#f55236")
    names(donors_color)<-donors_name

    # c('MA'="#160cc9", 'SA'="#f55236")
    library(reshape2)
    shar_cellt

    ##All clonotypes
    # don <- "NIHMA_244"

    table_upset<-lapply(unique(shar_cellt$donors), function(don){
    cat("   ", don, "\n")

      list_clust_clon<-lapply(c("CD4.TRM", "CD4.TRM.like", "CD4.Teff"), function(celltype){
        cat("       ", celltype, "\n")
        # a<-"CD4.TRM"
        clonotypes<-unique(shar_cellt$TRB.aaSeqCDR3_donor[shar_cellt[,celltype]  == "1" & shar_cellt$donors == don])   #& fgal$clon.size.tag>1
        clonotypes<-clonotypes[!is.na(clonotypes)]
        return(clonotypes)
        }
      )

      # num_unique_b_sequences<-length(unique(unlist(list_clust_clon)))

      names(list_clust_clon)<-c("CD4.TRM", "CD4.TRM.like", "CD4.Teff")

      combs_num <- Reduce(c,lapply(2:length(list_clust_clon),
                  function(x) combn(1:length(list_clust_clon),x,simplify=FALSE)))

      combs<-lapply(combs_num, function(x) Reduce(intersect,list_clust_clon[x]))

      names(combs)<-unlist(lapply(combs_num, function(x){
        paste(names(list_clust_clon)[unlist(x)], collapse=".vs.")
      }))

      combs<-lapply(names(combs), function(aaah){
          # aaah<-names(combs)[1]
        if(str_count(aaah, "vs") > 1) return(unlist(combs[aaah]))

        aver<-combs[[aaah]][!combs[[aaah]] %in% combs$CD4.TRM.vs.CD4.TRM.like.vs.CD4.Teff]
        # names(aver)<-aaah
        return(aver)
      })

      names(combs)<-unlist(lapply(combs_num, function(x){
        paste(names(list_clust_clon)[unlist(x)], collapse=".vs.")
      }))



      combinations_length<-lapply(combs, length)

      # names(combinations_length)<-unlist(lapply(combs_num, function(x){
      #   paste(names(list_clust_clon)[unlist(x)], collapse=".vs.")
      # }))

      combin_unique<-lapply(list_clust_clon, function(differ){

        return(setdiff(differ, unique(unlist(combs))))
      })

      names(combin_unique)<-names(list_clust_clon)
      combinations_unique_length<-lapply(combin_unique, length)
      # names(combinations_unique_length)<-names(combin_unique)

      unlist(list(combinations_length, combinations_unique_length))

      a<-melt(unlist(list(combinations_length, combinations_unique_length)))
      colnames(a)<-don
      return(t(a))
    })

    upset_table_final<-do.call(rbind, table_upset)

    colSums(upset_table_final)
    write.csv(upset_table_final, "./upset_plots/resting/upset_table_TRM.csv")

    pdf(paste0("./upset_plots/resting/upset_comb_TRMs.pdf"))
    upset(shar_cellt,comb, n_intersections=20, width_ratio=0.2, height_ratio=0.5,
    annotations = list(
      donor_num=(
        ggplot(mapping=aes(fill=donors, colour="black"))
        + geom_bar(stat='count', position='fill', inherit.aes = TRUE)
        + scale_colour_manual(values=c('black' = "black"))
        + scale_fill_manual(values=donors_color)
        # + scale_y_continuous(labels=scales::percent_format())
        # + ylab("Asthma") + theme(legend.position="none")
        + ylab("Asthma") + theme(legend.position="none")
                #, geom_text(aes(y=clone_size, label=clone_size), vjust=-90, size=1.5)
                ),
      Clone_size=(
        ggplot(mapping=aes(y=clone_size))
          + geom_boxplot(na.rm=TRUE)
        ) + ylab(" Clone Size")
        ),
        queries=list(
        upset_query(
            intersect=c('CD4.TRM'),
            color='#99141B',
            fill='#99141B',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
        upset_query(
            intersect=c('CD4.TRM.like'),
            color='#FFA303',
            fill='#FFA303',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
        upset_query(
            intersect=c('CD4.Teff'),
            color='#FFE610',
            fill='#FFE610',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
          upset_query(set='CD4.TRM', fill='#99141B'),
          upset_query(set='CD4.TRM.like', fill='#FFA303'),
          upset_query(set='CD4.Teff', fill='#FFE610')
      ) #,
      # themes=upset_default_themes(
      #       axis.ticks.x=element_blank(),
      #       axis.text.x=element_blank(),
      #       axis.ticks.y=element_blank(),
      #       axis.text.y=element_blank()
        #     ),
        # intersections_matrix=theme(
        #     axis.ticks.x=element_blank(),
        #     axis.text.x=element_blank(),
        #     axis.text.y=element_blank(),
        #     axis.line.x = element_blank(),
        #     axis.line.y = element_blank()
        #     )
      # ),
      # base_annotations= list(
      # ' '=intersection_size(counts=TRUE)
      )
    # )
    dev.off()

    pdf(paste0("./upset_plots/resting/upset_comb_TRMs_blank.pdf"))
    upset(shar_cellt,comb, n_intersections=20, width_ratio=0.2, height_ratio=0.5,
    annotations = list(
      donor_num=(
        ggplot(mapping=aes(fill=donors, colour="black"))
        + geom_bar(stat='count', position='fill', inherit.aes = TRUE)
        + scale_colour_manual(values=c('black' = "black"))
        + scale_fill_manual(values=donors_color)
        # + scale_y_continuous(labels=scales::percent_format())
        # + ylab("Asthma") + theme(legend.position="none")
        + ylab("") + theme(legend.position="none")
                #, geom_text(aes(y=clone_size, label=clone_size), vjust=-90, size=1.5)
                ),
      Clone_size=(
        ggplot(mapping=aes(y=clone_size))
          + geom_boxplot(na.rm=TRUE)
        ) + ylab("")
        ),
        queries=list(
        upset_query(
            intersect=c('CD4.TRM'),
            color='#99141B',
            fill='#99141B',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
        upset_query(
            intersect=c('CD4.TRM.like'),
            color='#FFA303',
            fill='#FFA303',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
        upset_query(
            intersect=c('CD4.Teff'),
            color='#FFE610',
            fill='#FFE610',
            only_components=c('intersections_matrix', 'Intersection size')
        ),
          upset_query(set='CD4.TRM', fill='#99141B'),
          upset_query(set='CD4.TRM.like', fill='#FFA303'),
          upset_query(set='CD4.Teff', fill='#FFE610')
      ),
      themes=upset_default_themes( # upset_default_themes
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text=element_blank()
        #     ),
        # intersections_matrix=theme(
        #     axis.ticks.x=element_blank(),
        #     axis.text.x=element_blank(),
        #     axis.text.y=element_blank(),
        #     axis.line.x = element_blank(),
        #     axis.line.y = element_blank()
        #     )
      ),
      base_annotations= list(
      'Intersection size'= intersection_size(counts=TRUE, text=list(size=0.001))
      )
      )
    # )
    dev.off()

}

### Figure 3 in ASTHMA_AIRWAYS_2021/trajectory/monocle_final.R
### Supp Figure 3K in ASTHMA_AIRWAYS_2021/trajectory/monocle_final.R

### Figure 4

{ cat(redb("### Dotplot Resting -Fig4B ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  sc_cd4resting$cluster_disease<-paste0(sc_cd4resting$celltype_subset, "_", sc_cd4resting$orig.asthma)
  table(sc_cd4resting$cluster_disease)

  sc_cd4resting$cluster_disease2 <- sc_cd4resting$cluster_disease
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TRMDP_MA"] <- "A_TRMDP_MA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TRMSP_MA"] <- "B_TRMSP_MA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TCM_MA"] <- "C_TCM_MA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "Treg_MA"] <- "D_Treg_MA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TFH_MA"] <- "E_TFH_MA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "ThIFNr_MA"] <- "F_ThIFNr_MA"

  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TRMDP_SA"] <- "G_TRMDP_SA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TRMSP_SA"] <- "H_TRMSP_SA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TCM_SA"] <- "I_TCM_SA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "Treg_SA"] <- "J_Treg_SA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "TFH_SA"] <- "K_TFH_SA"
  sc_cd4resting$cluster_disease[sc_cd4resting$cluster_disease == "ThIFNr_SA"] <- "L_ThIFNr_SA"
  table(sc_cd4resting$cluster_disease2, sc_cd4resting$cluster_disease)

  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R")

  fconfigs = list(
    list(result_id = "./dotplots/", sufix = "figure4b_zscore", #Changing scale_mean = FALSE in /home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R to use the real Mean or Z-score
    object = "sc_cd4resting",
      edata = "sc_cd4resting@assays$RNA@data", metadata = "sc_cd4resting@meta.data[sc_cd4resting@meta.data$celltype_subset %in% c('TRMDP', 'TRMSP', 'TCM', 'Treg', 'TFH', 'ThIFNr'), ]",
      size = c(5.5,4),
      col=c('white', 'white', '#ed1602'),
      scale_mean=TRUE,
      size_scale=15,
      axis_x = list(
        col = "cluster_disease"
      ),
      features =  c("CREM", "DUSP1", "DUSP2", "DUSP4", "TNFAIP3", "FKBP5")
    )
  )

  pp_curtains = fig_plot_curtain(fconfigs, use_zscore_seurat = FALSE, verbose = 2, clust_cols=FALSE, scale_fun='size', size_limits=c(0,100), size_scale=10, cols_limits=c(-1, 1.5))


  # saveRDS(sc_cd4stim@meta.data, "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/sc_cd4stim_annotation_GZMB_tag.rds")


}

{ cat(red_bold("### Crater plots -Fig4A & Suppl Fig 4A  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # source("/home/ciro/scripts/handy_functions/devel/plots_crater.R")
  # source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")

  # source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater.R") # /home/ciro/scripts/handy_functions/devel/plots_crater.R
  # source('/mnt/BioHome/ciro/scripts/handy_functions/devel/plots_crater.R')
  source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/plots_crater.R") #Bigger point
  source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")

  outdir_crater<-paste0(outdir, "resting/craters_plot")
  dir.create(outdir_crater, showWarnings = FALSE)
  setwd(outdir_crater)

  dgea_dir = "/home/ciro/ad_hoc/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/cluster"
  edata = as.matrix(expm1(sc_cd4resting@assays$RNA@data))
  mygenes = grep("XIST|RPS4Y1", rownames(edata), value = TRUE, invert = TRUE)
  mdata = sc_cd4resting@meta.data
  fconfigs = lapply(X = levels(sc_cd4resting$cluster),
    FUN = function(x){
      list(fnames = c(
        "Male_MA_vs_SA" = paste0(dgea_dir, x, "/Male_SAvsMale_MA/results_Male_SAvsMale_MA_mastlog2cpm.csv"),
        "Female_MA_vs_SA" = paste0(dgea_dir, x, "/Female_SAvsFemale_MA/results_Female_SAvsFemale_MA_mastlog2cpm.csv")
      ), result_id = paste0("./colsize_cluster", x, "_"), selectss = list(c(sc_cd4resting_clust, x)),
        columns = c("sex_disease"), #limits_col = c(0, 9), limits_size = c(0, 120),
        # plot_squared = 6.5,
        highlight_genes = c("CREM", "DUSP1", "DUSP2", "DUSP4", "TNFAIP3", "DDIT4", "FKBP5", "GZMB", "GZMH", "ZNF683", "HOPX", "CCL4", "CTSW"))
  })
  # if(!exists("sc_cd4merged")) sc_cd4merged = merge(sc_cd4resting, sc_cd4stim)
  # edata = expm1(sc_cd4merged@assays$RNA@data) # MEMORY intensive
  # mygenes = grep("XIST|RPS4Y1", rownames(edata), value = TRUE, invert = TRUE)
  # mdata = sc_cd4merged@meta.data
  # rm(sc_cd4merged); gc(); edata = as.matrix(edata)

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
        feature_subset = c(mygenes),
        topgenes = c(fconfig$highlight_genes),
        lfcthresh = fc,
        # check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMB", "GZMH", "ZNF683", "HOPX", "CCL4", "CTSW"))),
        # column4stats = fconfig$columns,
        outputname = fconfig$result_id,
        plot_interactive = TRUE,
        plot_squared = fconfig$plot_squared,
        limits_col = fconfig$limits_col,
        limits_size = fconfig$limits_size
      )
    }
  }

  fconfigs = list(list(
      fnames = c(
        "Male_MA_vs_SA" = "/home/ciro/ad_hoc/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex_disease/SA_MalevsMA_Male/results_SA_MalevsMA_Male_mastlog2cpm.csv",
        "Female_MA_vs_SA" = "/home/ciro/ad_hoc/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/sex_disease/SA_FemalevsMA_Female/results_SA_FemalevsMA_Female_mastlog2cpm.csv"
      ), result_id = "./resting_", columns = c("sex_disease"), plot_squared = TRUE
      highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL")
  ))

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
        feature_subset = c(mygenes),
        topgenes = c(fconfig$highlight_genes),
        lfcthresh = fc,
        # check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMB", "GZMH", "ZNF683", "HOPX", "CCL4", "CTSW"))),
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

#### Figure 5
dir.create(paste0(outdir, "/stim"))

setwd(paste0(outdir, "/stim"))

{ cat(redb("### Crater plot stim Male & Female MAvsSA - Fig5A ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

        library(tidyverse)
        # dir.create("./crater_disease")
        dir.create("./crater_disease_good")


        ##### Crater plots
        source("/home/kmlanderos/scripts/handy_functions/devel/plots_crater.R") # /home/ciro/scripts/handy_functions/devel/plots_crater.R
        source('/mnt/BioHome/ciro/scripts/handy_functions/devel/plots_crater.R')
        source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
        # dir.create("craters", showWarnings = FALSE)

        edata =  expm1(sc_cd4stim@assays$RNA@data)
        mygenes = grep("XIST|RPS4Y1|^RP11|^RP", rownames(edata), value = TRUE, invert = TRUE) # Filter unwanted genes
        mdata = sc_cd4stim@meta.data

          f1 <- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female_good/Female_SAvsMA/SAvsMA/mastlog2cpm_results.csv")

          f1 <- f1 %>% mutate(padj = case_when(
            padj < 1*10**-100 ~ 1*10**-100,
            TRUE ~ padj
          ))

          write.csv(f1, "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female_good/Female_SAvsMA/SAvsMA/.mastlog2cpm_results.csv", row.names = F)

          f2 <- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female_good/Male_SAvsMA/SAvsMA/mastlog2cpm_results.csv")
          f2 <- f2 %>% mutate(padj = case_when(
            padj < 1*10**-100 ~ 1*10**-100,
            TRUE ~ padj
          ))

          write.csv(f2, "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female/Female_SAvsMA/.mastlog2cpm_results.csv", row.names = F)


          # NOTE: Change Version if necessary
          fconfigs = list(list(
              fnames = c(
                "Male_SAvsMA" = "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female/Female_SAvsMA/.mastlog2cpm_results.csv",
                "Female_SAvsMA" = "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/stim/dgea/MAvsSA_Male_and_Female_good/Female_SAvsMA/SAvsMA/.mastlog2cpm_results.csv"
              ), result_id = "./crater_disease_good/Female_Male_Ma_vsSA", #plot_squared = TRUE, #, columns = c("sex_disease")
               # highlight_genes = c("ALOX5AP", "CKLF", "BATF", "IFNG", "NR3C1", "METRNL"),
              selectss = list(c('orig.Sex', 'Male', 'Female'), c('orig.asthma','MA', 'SA')) # brain, lung; # CD4, CD8
            )
          )

          degfilt = list(mean = list("<0", NA), min_padj = list(">0.05", 1))

          for(fconfig in fconfigs){
            cat(c("- ", fconfig$result_id, "\n"), sep = "")
            for(fc in as.character(c(0.25))){
              cat(" * FC:", fc, "\n");
              void <- crater_plot(
                tests_list = fconfig$fnames,
                edataf = edata,
                annotf = mdata,
                sample_filter = fconfig$selectss,
                gene_filter = degfilt,
                feature_subset = mygenes,
                # topgenes = c("top10", fconfig$highlight_genes),
                lfcthresh = fc,
                # column4stats = fconfig$columns,
                outputname = fconfig$result_id,
                plot_interactive = TRUE,
                plot_squared = fconfig$plot_squared,
                limits_col = fconfig$limits_col,
                limits_size = fconfig$limits_size,
                verbose = TRUE
                # return_out = TRUE
              )
            }
        }
}

cat(redb("### Dotplot UNSTIM & STIM -Fig5B  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")) {

  a<-load("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/integration_CD4/integrated_cd4/zetInfo_integrated_30ANC_24PC_obj.RData")

  setwd(outdir)

  dir.create("resting_and_stim")
  setwd("resting_and_stim")
  dir.create("dotplot")

  source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/lung08/scripts/plots_dotplot.R")

  immune.combined$orig.stim_disease <- paste0(immune.combined$orig.stim, "_", immune.combined$orig.asthma)

  immune.combined$orig.stim_disease2<-ifelse(grepl("US", immune.combined$orig.stim_disease), paste0("A_", immune.combined$orig.stim_disease), immune.combined$orig.stim_disease)

  table(immune.combined$orig.stim_disease, immune.combined$orig.stim_disease2)

  # source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R")

  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R")

  fconfigs = list(
    list(result_id = "./dotplot/", sufix = "genes_zscore_good_June20", #Changing scale_mean = FALSE in /home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R to use the real Mean or Z-score
    object = "immune.combined",
      edata = "immune.combined@assays$RNA@data", metadata = "immune.combined@meta.data",
      size = c(4,8.5),
      col=c('#f7f7f0', '#f5f5b8', '#e3e312', '#f23513', '#a81103'),
      scale_mean=TRUE,
      size_scale=15,
      axis_x = list(
    col = "orig.stim_disease2"
    # col = "celltype", order = c("Basal", "Club", "Ciliated", "Fibroblast", "Ionocyte", "Tcell", "Bcell", "DC", "Neutrophils", "Mast", "Cycling")
      ),
      features = c("ITGAE", "ITGA1", "GZMB", "FASLG", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "IFNG", "TNF", "IL4", "IL5", "IL13", "CSF2", "IL17A", "IL17F", "CCL20", "IL21", "TNFSF14", "AREG", "TGFB1")
    )
  )

  pp_curtains = fig_plot_curtain(fconfigs, use_zscore_seurat = FALSE, verbose = 2, clust_cols=FALSE, scale_fun='size', size_limits=c(0,100), size_scale=10)

}

setwd(paste0(outdir, "/stim"))

{ cat(redb("### Blanks volcano GZMBpvsGZMBn -Fig 5C ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/dgea")
  #I ran the documented line and change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/dgea/MAvsSA_Male_and_Female_good_sex_cov/GZMBpvsGZMn/GZMBpvsGZMBn/mastlog2cpm_results.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "GZMBp"
    group2 = "GZMBn"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

     cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  dir.create("volcano")
  pdf("volcano/GZMBpvsGZMMn_sex_cov_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-2.5,2.5),
      ylims=c(0,320),
      verbose = verbose
      # interact= TRUE
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

  p_html <- volplot(
    resdf,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = "padj",
    lfctype = "log2FoldChange",
    col_feature = means,
    size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
    gene_name = "gene",
    group = "degs",
    check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
    return_plot = TRUE,
    clipp = 4,
    xlims=c(-2.5,2.5),
    ylims=c(0,320),
    verbose = verbose ,
    interact= TRUE
  )

  htmlwidgets::saveWidget(
          as_widget(p_html), 'volcano/GZMBpvsGZMMn_sex_cov_blank.html')


  # datavis <- GetAssayData(gpr25_exp1_day30)
  # annot <- gpr25_exp1_day30@meta.data
  # datatype <- "SeuratNormalized"
  # verbose <- TRUE
  #
  # aver<-stats_summary_table(
  #   mat = datavis,
  #   groups = make_list(x = annot, colname = 'origlib', grouping = TRUE),
  #   rnames = rownames(datavis),
  #   datatype = datatype,
  #   verbose = verbose
  # )

  # write.csv(aver, "./volcano/mean_and_percetage_per_gene.csv")

}

cat(redb("### Dotplot STIM - Fig5D ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")) {
  ### Dotplots

  dir.create("./dotplots")


  # tvar <- sc_cd4resting@assays$RNA@counts["GZMB", ] > 0
  # sc_cd4resting@meta.data$tag_GZMB = ifelse(tvar, "GZMBp", "GZMBn")
  tvar <- sc_cd4stim@assays$RNA@counts["GZMB", ] > 0
  sc_cd4stim@meta.data$tag_GZMB = ifelse(tvar, "GZMBp", "GZMBn")

  # source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/lung08/scripts/plots_dotplot.R"

  sc_cd4stim$gzmb_seps<-ifelse(sc_cd4stim$tag_GZMB == "GZMBp", paste0(sc_cd4stim$orig.asthma, "_", sc_cd4stim$tag_GZMB), sc_cd4stim$tag_GZMB)

  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R")

  fconfigs = list(
    list(result_id = "./dotplots/", sufix = "GZMBp_disease_GZMBn_zscore_20Jun", #Changing scale_mean = FALSE in /home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R to use the real Mean or Z-score
    object = "sc_cd4stim",
      edata = "sc_cd4stim@assays$RNA@data", metadata = "sc_cd4stim@meta.data",
      size = c(4,5.5),
      col=c('#f7f7f0', '#f5f5b8', '#e3e312', '#f23513', '#a81103'),
      scale_mean=TRUE,
      size_scale=15,
      axis_x = list(
    col = "gzmb_seps"
    # col = "celltype", order = c("Basal", "Club", "Ciliated", "Fibroblast", "Ionocyte", "Tcell", "Bcell", "DC", "Neutrophils", "Mast", "Cycling")
  ),
      features =  c("GZMA", "GZMH", "CCL3", "CCL4", "CCL5", "IFNG", "TNF", "CSF2", "CCL20", "IL21", "TNFSF14", "TGFB1")
      #"CCL3", "CCL4", "CCL5", "CCL20", "GZMA", "GZMH", "IL21", "XCL1", "TNFSF14", "CSF2"

      #c("ITGAE", "ITGA1", "GZMB", "FASLG", "CCL3", "CCL4", "CCL5", "CCL20", "XCL1", "XCL2", "IFNG", "TNF", "IL4", "IL5", "IL13", "IL17A", "IL17F", "IL21", "IL10", "CSF2", "TNFSF14",  "AREG", "TGFB1")
    )
  )

  #  c("GZMB", "CCL3", "CCL4", "CCL5", "FASLG", "IFNG", "GZMA", "GZMH", "TNFSF14", "METRNL", "CSF2", "IL21", "CCL20", "XCL1", "XCL2", "TNF", "PRF1", "TGFB1")

  pp_curtains = fig_plot_curtain(fconfigs[1], use_zscore_seurat = FALSE, verbose = 2, clust_cols=FALSE, scale_fun='size', size_limits=c(0,100), size_scale=10)

}

{ cat(red_bold("### Coexpression Scatter Stim - Fig5E ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  dir.create("./scatter_coexpression")
  dir.create("./scatter_coexpression/14_June/")

  # dir.create("figures_vijay_grant/scatter_contour/sc_biop_cd8_AREGvsOthers")
  fconfigs = list(
    list(result_id = "./scatter_coexpression/14_June/",
      edata = "sc_cd4stim@assays$RNA@data", metadata = "sc_cd4stim@meta.data[sc_cd4stim@meta.data$tag_GZMB == 'GZMBp', ]",
      axis_x = "RNA_snn_res.0.4", plot_main = FALSE,
      features = list(x = c("CCL3"), y = c("CCL4", "CCL5", "IL21", "IFNG", "TNF", "TGFB1"))
    )
  )

  plot_add_quadrants<- function(plot, limits = list(0.5, 0.5), ...) {
    plot +
      geom_vline(xintercept = limits[[1]], linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = limits[[2]], linetype = "dashed", color = "gray50") +
      stat_quadrant(xintercept = limits[[1]], yintercept = limits[[2]], ...)
  }

  pp_contour = fig_plot_contour(fconfigs[1],
    theme_extra = function(x){
      plot_add_quadrants(x +  viridis::scale_color_viridis(option = "magma"), type = "percent") + geom_point(size=3)
    }, return_plot = TRUE
  )

}

### Supplementary Figures
setwd(paste0(outdir, "/resting"))
{ cat(redb("### Blanks volcano TRMDP vs TRMSP Treatment covariate -Supp Fig3A ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting")
  #I ran the documented line and change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/TRMDSvsTRMSP_treatment_cov/TRMDPvsTRMSP/TRMDPvsTRMSP/mastlog2cpm_results.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "TRMDP"
    group2 = "TRMSP"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

     cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  dir.create("volcano")
  pdf("volcano/TRMDPvsTRMSP_Treatment_cov_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-3,3),
      ylims=c(0,320),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()


  pdf("volcano/TRMDPvsTRMSP_Treatment_cov.pdf") #Change name for CD8
    volplot_wnames(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA", "IFNG", "TNFSF14", "FKBP5", "DDIT4", "CTSW"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-3,3),
      ylims=c(0,320),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
    dev.off()


}

{ cat(redb("### Correlation DGEA TRMDP vs TRMSP sexcov and treatment covariate -Sup Fig3B & Sup Fig3C ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  library(data.table)
  library(dplyr)

  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting")
  dir.create("dgeas_comp")

  treat_cov<-read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/TRMDSvsTRMSP_treatment_cov/TRMDPvsTRMSP/TRMDPvsTRMSP/mastlog2cpm_results.csv") %>% data.table()

  sex_cov<-read.csv("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_pjs/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_cov_sex/0vs1/results_0vs1_mastlog2cpm.csv") %>% data.table()

  treat_cov2<-treat_cov[abs(log2FoldChange) >= 0.25 & padj < 0.05, .(gene=gene, log2FC_treat=log2FoldChange)]
  sex_cov2<-sex_cov[abs(log2FoldChange) >= 0.25 & padj < 0.05, .(gene=gene, log2FC_sex=log2FoldChange)]

  # dgeas_correlation<-merge(treat_cov2, sex_cov2, by="gene")
  dgeas_correlation2<-merge(treat_cov2, sex_cov2, by="gene")


  library(ggvenn)

  dgeas_intersect<-list(
    treat_cov2 = treat_cov2$gene,
    sex_cov2 = sex_cov2$gene
  )

  pdf("./dgeas_comp/TRMDPvsTRMDP_covSex_and_covTreat_check.pdf")
  ggvenn(dgeas_intersect)
  dev.off()

  library(ggpubr)
  pdf("./dgeas_comp/correlation_TRMDPvsTRMDP_covSex_and_covTreat_check.pdf")
  ggscatter(dgeas_correlation2, x = "log2FC_sex", y = "log2FC_treat",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "log2FC TRMDP vs TRMSP sex covariate", ylab = "log2FC TRMDP vs TRMSP Treatment as covariate")
  dev.off()


  sex_cov23<-sex_cov[abs(log2FoldChange) >= 0.25 & padj < 0.05 & minExp0in3 == "TRUE", .(gene=gene, log2FC_sex=log2FoldChange)]
  treat_cov23<-treat_cov[abs(log2FoldChange) >= 0.25 & padj < 0.05 & minExp0in4933samples == "TRUE", .(gene=gene, log2FC_treat=log2FoldChange)]

  dgeas_correlation_23<-merge(treat_cov23, sex_cov23, by="gene")
  library(ggvenn)

  dgeas_intersect_23<-list(
    treat_cov2 = treat_cov23$gene,
    sex_cov2 = sex_cov23$gene
  )

  #Sup Fig4B
  pdf("./dgeas_comp/TRMDPvsTRMDP_covSex_and_covTreat_26_Jun.pdf")
  ggvenn(dgeas_intersect_23)
  dev.off()

  #Sup Fig4C
  library(ggpubr)
  pdf("./dgeas_comp/correlation_TRMDPvsTRMDP_covSex_and_covTreat_26_Jun.pdf")
  ggscatter(dgeas_correlation_23, x = "log2FC_sex", y = "log2FC_treat",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "log2FC TRMDP vs TRMSP sex covariate", ylab = "log2FC TRMDP vs TRMSP Treatment as covariate")
  dev.off()



  length_higher<-length(as.character(intersect(dgeas_intersect[[1]], dgeas_intersect[[2]])))
  intersection_table

  intersection_dgea <- as.character(intersect(dgeas_intersect[[1]], dgeas_intersect[[2]]))

  treat_cov_unique <- as.character(dgeas_intersect[[1]][!dgeas_intersect[[1]] %in% dgeas_intersect[[2]]])

  sex_cov_unique <- as.character(dgeas_intersect[[2]][!dgeas_intersect[[2]] %in% dgeas_intersect[[1]]])

  length(treat_cov_unique) <- length_higher
  length(sex_cov_unique) <- length_higher

  intersec_table<-cbind(intersection_dgea, treat_cov_unique,sex_cov_unique)

  write.csv(intersec_table, "./dgeas_comp/intersection_table_check.csv")

  # list.rbind(.intersection_table)
}

{ cat(redb("### Blanks volcano TRMDP biologics_vs_nonBiologics -SupFig 3F ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting")
  #I ran the documented line and change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/Biologics_vs_NoBiologics_and_Corticoisteroid/Biologics_vs_NoBiologics_TRMDP/biologicsvsno_biologics/mastlog2cpm_results.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "biologics"
    group2 = "no_biologics"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

     cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  dir.create("volcano")
  pdf("volcano/TRMDP_biologics_vs_nonBiologics_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-1.55,1.55),
      ylims=c(0,65),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()


  pdf("volcano/TRMDP_biologics_vs_nonBiologics.pdf") #Change name for CD8
    volplot_wnames(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(showgenes)),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-1.55,1.55),
      ylims=c(0,65),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
    dev.off()


  # datavis <- GetAssayData(gpr25_exp1_day30)
  # annot <- gpr25_exp1_day30@meta.data
  # datatype <- "SeuratNormalized"
  # verbose <- TRUE
  #
  # aver<-stats_summary_table(
  #   mat = datavis,
  #   groups = make_list(x = annot, colname = 'origlib', grouping = TRUE),
  #   rnames = rownames(datavis),
  #   datatype = datatype,
  #   verbose = verbose
  # )

  # write.csv(aver, "./volcano/mean_and_percetage_per_gene.csv")

}

{ cat(redb("### Blanks volcano TRMDP OCS_vs_noOCS -SupFig 3E### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting")
  #I ran the documented line and change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/mnt//bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/Biologics_vs_NoBiologics_and_Corticoisteroid/OCS_vs_No_OCS_TRMDP/OCSvsnOCS/mastlog2cpm_results.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "OCS"
    group2 = "nOCS"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

     cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  dir.create("volcano")
  pdf("volcano/TRMDP_OCS_vs_nonOCS_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-2.5,2.5),
      ylims=c(0,90),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()


  pdf("volcano/TRMDP_OCS_vs_nonOCS.pdf") #Change name for CD8
    volplot_wnames(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(showgenes)),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-2.5,2.5),
      ylims=c(0,65),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
    dev.off()


  # datavis <- GetAssayData(gpr25_exp1_day30)
  # annot <- gpr25_exp1_day30@meta.data
  # datatype <- "SeuratNormalized"
  # verbose <- TRUE
  #
  # aver<-stats_summary_table(
  #   mat = datavis,
  #   groups = make_list(x = annot, colname = 'origlib', grouping = TRUE),
  #   rnames = rownames(datavis),
  #   datatype = datatype,
  #   verbose = verbose
  # )

  # write.csv(aver, "./volcano/mean_and_percetage_per_gene.csv")

}

{ cat(redb("### Public datasets GZMB+ Healthy - Sup Fig 3I ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  #Integration in /Volumes/bioadhoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/Rebuttal_CD4/scripts/integration/integration_clustering.R
  dir.create("./tables_clustering")

  tcells_combined<-readRDS("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/Rebuttal_CD4/results/integration/hvg1500_pc20_4studies_genesAirways/object.rds")
  redu = list(umap = c("UMAP_1", "UMAP_2"))
  # colsname <- "/home/ciro/scripts/handy_functions/data/colours.csv"
  # colors_df <- read.csv(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
  colors_df<-vector()

  tcells_combined@meta.data = joindf(
    tcells_combined@meta.data,
    as.data.frame(tcells_combined@reductions$umap@cell.embeddings))

  #Changing name of study to a smaller one, just to representation in the plot
    table(tcells_combined$orig.study)
    tcells_combined$orig.study2<-tcells_combined$orig.study
    #orig.study2<-ifelse(tcells_combined$orig.celltype=="CD8", "_CD8", "")

    tcells_combined$orig.study<-gsub("NA", "", paste(tcells_combined$orig.study2, tcells_combined$study, sep=""))
    table(tcells_combined$orig.study)

    tcells_combined@meta.data$orig.study_disease<-gsub("_NA", "", paste(tcells_combined$orig.study, tcells_combined$orig.disease_severity, sep="_"))

    table(tcells_combined@meta.data$orig.study, tcells_combined@meta.data$orig.study_disease)
    table(tcells_combined@meta.data$orig.study_disease, tcells_combined@meta.data$orig.disease_severity)

    unique(tcells_combined@meta.data$orig.study[is.na(tcells_combined@meta.data$orig.disease_severity)]) #Which cells doesn't have mild or severe asthma (The healthy ones)

    tcells_combined@meta.data$orig.study_disease2<-"study_disease2"
    tcells_combined$orig.study_disease2[tcells_combined$orig.study_disease=="liao_naturemedicine_2020"]<- "liao_healthy"
    tcells_combined$orig.study_disease2[tcells_combined$orig.study_disease=="Grant_nature_2021"]<- "grant_healthy"
    tcells_combined$orig.study_disease2[tcells_combined$orig.study_disease=="herrera_mag_2022_MA"]<- "herrera_MA"
    tcells_combined$orig.study_disease2[tcells_combined$orig.study_disease=="herrera_mag_2022_SA"]<- "herrera_SA"
    tcells_combined$orig.study_disease2[tcells_combined$orig.study_disease=="morse_ERJ_2019"]<- "morse_healthy"


    table(tcells_combined@meta.data$orig.study, tcells_combined@meta.data$orig.study_disease2)
    table(tcells_combined@meta.data$orig.study_disease2, tcells_combined@meta.data$orig.disease_severity)

    tcells_combined@meta.data$orig.disease_severity_healthy<-gsub("Healthy_NA", "Healthy", paste(tcells_combined$orig.disease, tcells_combined$orig.disease_severity, sep="_"))

    table(tcells_combined$orig.celltype)
    tcells_combined$orig.celltype[grepl("CD4", tcells_combined$orig.celltype)]<-"CD4_Healthy"
    table(tcells_combined$orig.celltype)


  #Idents(tcells_combined) <- tcells_combined$orig.Asthma_Severity2
  tcells_combined_disease<-SplitObject(tcells_combined, split.by = "orig.disease")
  tcells_combined_healthy<-tcells_combined_disease$Healthy

  #For healthy libraries.. I know, i know this could be a`function jijiji
  tvar<-c("GZMB", "ITGAE")
  tmp<-add_gene_tag(tvar, tcells_combined_healthy@meta.data, as.matrix(tcells_combined_healthy@assays$RNA@data[tvar,]), thresh = 0.01)

  tcells_combined_healthy@meta.data$GZMB_expr<-tmp$tag_GZMB
  tcells_combined_healthy@meta.data$ITGAE_expr<-tmp$tag_ITGAE

  mdata<-tcells_combined_healthy@meta.data
  gz<-table(mdata$orig.subject, mdata$GZMB_expr=="GZMB+")[,2]
  itgea<-table(mdata$orig.subject, mdata$ITGAE_expr=="ITGAE+")[,2]

  gzpos_itgea_pos<-table(mdata$orig.subject, mdata$GZMB_expr=="GZMB+" & mdata$ITGAE_expr=="ITGAE+")[,2]

  table(is.na(mdata$orig.subject)) #There are no NA's

  total<-as.matrix(table(mdata$orig.subject))[,1]

  gene_expr_donor_count<-data.frame(gz,itgea,gzpos_itgea_pos, total)

  colnames(gene_expr_donor_count)<-c("GZMB+",  "ITGAE+", "GZMB+ITGAE+", "Total_cells")
  write.table(gene_expr_donor_count,"./tables_clustering/genes_counts_perdonor_healthy_integrated.csv", sep=",", col.names = NA, row.names=TRUE) #row.names=TRUE ?? what for? doen't change something

  gene_expr_donor_perc<-gene_expr_donor_count/gene_expr_donor_count$Total_cells #I doesn't round the number but i did it in the blood tables /mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/SaHe_Blood/scripts/figures/init.R section  Gene tables (GZMB+ PRF1+ and GZMB+PRF1+)
  #gene_expr_donor_perc$Total_cells<-total
  write.table(gene_expr_donor_perc,"./tables_clustering/genes_percentage_perdonor_healthy_integrated.csv", sep=",", col.names = NA)

  #For BAL .. I know, i know this could be a function 

    tcells_combined_asthma<-tcells_combined_disease$Asthma
    tmp<-add_gene_tag(tvar, tcells_combined_asthma@meta.data, as.matrix(tcells_combined_asthma@assays$RNA@data[tvar,]), thresh = 0.01)

    tcells_combined_asthma@meta.data$GZMB_expr<-tmp$tag_GZMB
    tcells_combined_asthma@meta.data$ITGAE_expr<-tmp$tag_ITGAE

    mdata<-tcells_combined_asthma@meta.data
    gz<-table(mdata$orig.subject, mdata$GZMB_expr=="GZMB+")[,2]
    itgea<-table(mdata$orig.subject, mdata$ITGAE_expr=="ITGAE+")[,2]

    gzpos_itgea_pos<-table(mdata$orig.subject, mdata$GZMB_expr=="GZMB+" & mdata$ITGAE_expr=="ITGAE+")[,2]

    table(is.na(mdata$orig.subject)) #There are no NA's. All should be FALSE, table function works because is a TRUE/FALSE vector

    total<-as.matrix(table(mdata$orig.subject))[,1]

    gene_expr_donor_count<-data.frame(gz,itgea,gzpos_itgea_pos, total)

    colnames(gene_expr_donor_count)<-c("GZMB+",  "ITGAE+", "GZMB+ITGAE+", "Total_cells")
    write.table(gene_expr_donor_count,"./tables_clustering/genes_counts_perdonor_Asthma_integrated.csv", sep=",", col.names = NA, row.names=TRUE) #row.names=TRUE ?? what for? doen't change something

    gene_expr_donor_perc<-gene_expr_donor_count/gene_expr_donor_count$Total_cells #I doesn't round the number but i did it in the blood tables /mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/SaHe_Blood/scripts/figures/init.R section  Gene tables (GZMB+ PRF1+ and GZMB+PRF1+)
    #gene_expr_donor_perc$Total_cells<-total
    write.table(gene_expr_donor_perc,"./tables_clustering/genes_percentage_perdonor_Asthma_integrated.csv", sep=",", col.names = NA)
}

{ cat(red_bold("### Kara Mould in healthy GZMB expression -Sup Fig4I ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

  dir.create("tables_clustering/Mould")
  setwd("tables_clustering/Mould")
  a<-readRDS("/home/fcastaneda/fcastaneda/rnaseq-sc-standar/Rebuttal_CD4/raw/Mould/seurat.rds")

  tvar <- c("CD4", "CD8B", "GZMB")
  #tmp2 <- add_gene_tag(tvar, liao_sep@meta.data, as.matrix(liao_sep@assays$RNA@data[tvar, ]))
  tmp <- add_gene_tag(tvar, a@meta.data, as.matrix(a@assays$RNA@data[tvar, ]), thresh=0.01)
  all(rownames(a@meta.data) == rownames(tmp))

  a$tag_GZMB<-tmp$tag_GZMB
  a$tag_CD4<-tmp$tag_GZMB
  a$tag_CD8<-tmp$tag_GZMB

  a$CD4_CD8<-paste(tmp$tag_CD4,tmp$tag_CD8, sep="")
  a$keep<-ifelse(a$CD4_CD8 == "CD4+CD8B-", "CD4", "Tcell_noCD4")
  a$cluster_cd4<-paste0(a$overall_cluster, "_", a$keep)
  mould<-SplitObject(a, "cluster_cd4")
  tcells_cd4_mould<-mould$'3_CD4'
  table(tcells_cd4_mould$orig.ident, tcells_cd4_mould$tag_GZMB)/rowSums(table(tcells_cd4_mould$orig.ident, tcells_cd4_mould$tag_GZMB))

  prop_GZMB_in_healthy_mould<-table(tcells_cd4_mould$orig.ident, tcells_cd4_mould$tag_GZMB)

  make_prop<-function(dcounts, name){
     write.table(dcounts, paste0("./",name,"_counts.csv") , col.names=NA, sep=",")
     write.table(dcounts/rowSums(dcounts), paste0("./", name,"_proportions.csv"), col.names=NA, sep=",")
   }

  make_prop(prop_GZMB_in_healthy_mould, "GZMBp_CD4Tcells_Mould")
}

{ cat(red_bold("### GSEA -SupFig 4B & SupFig 4C GSEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
    source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
    source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
    dir.create("gsea")

    gdf <- readfile("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/info/signatures_table_2.csv", stringsAsFactors = FALSE)

    slist = gsea_process_list(lapply(gdf, function(x) x[-1] ))
    slist2<-slist[c(1,11)]

    fconfigs = list(
      list(
        result_id = "./gsea/sc_cd4resting_TRMDP_comp_3F/",
        edata = "expm1(sc_cd4resting@assays$RNA@data)",
        metadata = "sc_cd4resting@meta.data[sc_cd4resting$celltype_subset=='TRMDP',]",
        lists = slist2, comparisons = data.frame(
          group1 = c("SA", "Male"), group2 = c("MA", "Female"),
          column = c("orig.asthma", "orig.Sex"), stringsAsFactors = FALSE)
        ),
      list(
        result_id = "./gsea/sc_cd4resting_TRMDP_severe_comp_3F/",
        edata = "expm1(sc_cd4resting@assays$RNA@data)",
        metadata = "sc_cd4resting@meta.data[sc_cd4resting$celltype_subset=='TRMDP' & sc_cd4resting$orig.asthma == 'SA',]",
        lists = slist2, comparisons = data.frame(
          group1 = c("Male"), group2 = c("Female"),
          column = c("orig.Sex"), stringsAsFactors = FALSE)
          )
      )

      pp_gsea = fig_gsea(fconfigs)

  }

{ cat(redb("### Dotplot STIM -SupFig 5B  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  setwd(paste0(outdir, "/resting_and_stim/"))

  # a<-load("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/integration_CD4/integrated_cd4/zetInfo_integrated_30ANC_24PC_obj.RData")

  # setwd("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/")
  dir.create("dotplot")

  source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/lung08/scripts/plots_dotplot.R")

  immune.combined$orig.stim_disease <- paste0(immune.combined$orig.stim, "_", immune.combined$orig.asthma)

  immune.combined$orig.stim_disease2<-ifelse(grepl("US", immune.combined$orig.stim_disease), paste0("A_", immune.combined$orig.stim_disease), immune.combined$orig.stim_disease)

  table(immune.combined$orig.stim_disease, immune.combined$orig.stim_disease2)

  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R")

  fconfigs = list(
    list(result_id = "./dotplot/", sufix = "genes", #Changing scale_mean = FALSE in /home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/source.R to use the real Mean or Z-score
    object = "immune.combined",
      edata = "immune.combined@assays$RNA@data", metadata = "immune.combined@meta.data",
      size = c(5,3),
      col=c('#f7f7f0', '#f5f5b8', '#e3e312', '#f23513', '#a81103'),
      scale_mean=FALSE,
      size_scale=15,
      axis_x = list(
    col = "orig.stim_disease2"
    # col = "celltype", order = c("Basal", "Club", "Ciliated", "Fibroblast", "Ionocyte", "Tcell", "Bcell", "DC", "Neutrophils", "Mast", "Cycling")
      ),
      features =  c("GZMB","GZMA","GZMH","GZMK","GZMM","GNLY")
    )
  )

  pp_curtains = fig_plot_curtain(fconfigs, use_zscore_seurat = FALSE, verbose = 2, clust_cols=FALSE, scale_fun='size', size_limits=c(0,100), size_scale=10)

}

setwd(paste0(outdir, "/stim"))
{ cat(red_bold("### GZMB stats -SupFig 5C ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  gzmb_stat<-as.data.frame.matrix(t(table(sc_cd4stim$tag_GZMB, sc_cd4stim$orig.donor))/colSums(table(sc_cd4stim$tag_GZMB, sc_cd4stim$orig.donor)))

  sexes<-distinct(sc_cd4stim@meta.data, orig.donor, orig.Sex)

  gzmb_stat$orig.donor <- rownames(gzmb_stat)
  gzmb_stats2<-merge(gzmb_stat, sexes, "orig.donor")

  write.csv(gzmb_stats2, "./gzmb_stats.csv")
}

{ cat(redb("### Blanks volcano  GZMBp SA vs MA sexcov -SupFig 5D ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/dgea")
  #I ran the documented line and change the name of the output file
  # dir.create("volcano")
  results<- read.csv("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/dgea/GZMBp_diasease_comp_sexcov_plot/GZMBp_SAvsMA/SAvsMA/mastlog2cpm_results.csv")

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "SA"
    group2 = "MA"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

    cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  dir.create("volcano")
  pdf("volcano/GZMBp_SAvsMA_sex_cov_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-2,2),
      ylims=c(0,140),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

  p_html <- volplot(
    resdf,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = "padj",
    lfctype = "log2FoldChange",
    col_feature = means,
    size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
    gene_name = "gene",
    group = "degs",
    check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))), #Delate GZMA??
    return_plot = TRUE,
    clipp = 4,
    xlims=c(-2,2),
    ylims=c(0,140),
    verbose = verbose ,
    interact= TRUE
  )

  htmlwidgets::saveWidget(
          as_widget(p_html), 'volcano/GZMBp_SAvsMA_sex_cov_blank.html')


  # datavis <- GetAssayData(gpr25_exp1_day30)
  # annot <- gpr25_exp1_day30@meta.data
  # datatype <- "SeuratNormalized"
  # verbose <- TRUE
  #
  # aver<-stats_summary_table(
  #   mat = datavis,
  #   groups = make_list(x = annot, colname = 'origlib', grouping = TRUE),
  #   rnames = rownames(datavis),
  #   datatype = datatype,
  #   verbose = verbose
  # )

  # write.csv(aver, "./volcano/mean_and_percetage_per_gene.csv")

}

{ cat(redb("### BULK Blanks volcano CD4 TRM Stim vs Unstim -SupFig 5E ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/dgea")
  dir.create("CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_plot/")
  #I ran the documented line and change the name of the output file... what?? June 20 2023
  # dir.create("volcano")
  results<-NULL
  resdf<-NULL

  results<- read.csv("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_airways/results/deseq2/airways_platex/comprs/activation/CD4-TRM_StimvsCD4-TRM_Unstim/results_CD4-TRM_StimvsCD4-TRM_Unstim_deseq2.csv")

  results$group<-gsub("CD4-TRM", "CD4.TRM", results$group)

  #results<-read.csv("/home/fcastaneda/vd-vijay/cramirez/asthma_biopsy/results/dgea/p5pe_clean1/epithelial_covsex/NIHWvsNIHMA/mastlog2cpm_results.csv")

    group1 = "CD4.TRM_Stim"
    group2 = "CD4.TRM_Unstim"
    padjthr = 0.05
    fcthr = 0.25
    cols_filt = "minExp"
    output = "./"
    verbose = TRUE
    return_report = FALSE
    ##############
    pseuc = 1
    dtype = "CST"
    # volcano parameters
    showgenes = NULL
    # Rebbuttal parameters
    cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
    means = NULL
    if(is.null(cols_filt)) cols_filt = "pattern123"
    the_report = list()
    if(is.null(names(group1))) names(group1) <- group1
    if(is.null(names(group2))) names(group2) <- group2

    if(verbose) cat("---------------------- Filtering results ---------------\n")
    resdf <- data.frame(
      results[which(results$padj <= 1), ],
      stringsAsFactors = FALSE, check.names = FALSE
    ) # order should ideally be by padj all the way through
    resdf <- resdf[order(resdf$padj), ]
    if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
    resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
    rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
    tvar <- list(
      which(resdf[, "log2FoldChange"] <= -fcthr),
      which(resdf[, "log2FoldChange"] >= fcthr)
    )
    if(!"group" %in% colnames(resdf)){
      if(verbose) cat("Addding 'group' column\n")
      resdf$group = NA; resdf$group[tvar[[1]]] <- group1
      resdf$group[tvar[[2]]] <- group2
    }else{
      tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
    }
    group_m <- sapply(c(names(group1), names(group2)), function(x){
      grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
    })
    if(length(group_m) == 2){
      if(verbose){
        cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
        str(group_m, vec.len = 10, no.list = 4)
      }
      resdf$Mean <- 0; means = "Mean"
      resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
      resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
    }
    genes2plot <- mysignames <- getDEGenes(
      resdf, pv = padjthr, fc = fcthr,
      gene_name = "gene", further = NULL, verbose = verbose
    )
    if(any(grep(cols_filt, colnames(resdf)))){
      tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
      genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
      if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
    }
    resdf$degs <- "Not_significant"
    resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
    if(length(genes2plot) == 0)
      genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
    if(length(genes2plot) == 0)
      genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
    if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
    resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
    if("pct_diff" %in% colnames(resdf))
      resdf$pct_diff <- abs(resdf$pct_diff)
      resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
    if(is.null(showgenes)){
      showgenes <- unique(c(
        bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
        head(resdf[genes2plot, "gene"], 10)))
    }

    cat("---------------------- Volcano -------------------------\n")

    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/scripts/volplot.R") # volplot
    # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
    source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
    source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
    source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

    # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

    # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
    #ephitelials
  # dir.create("volcano")

  library("plotly")
  pdf("CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_plot/CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_blank.pdf") #Change name for CD8
    volplot(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-15,15),
      ylims=c(0,300),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

  p_html <- volplot(
    resdf,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = "padj",
    lfctype = "log2FoldChange",
    col_feature = means,
    size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
    gene_name = "gene",
    group = "degs",
    check_genes = list(text = features_parse_ensembl(c(showgenes, "GZMA"))), #Delate GZMA??
    return_plot = TRUE,
    clipp = 4,
    xlims=c(-15,15),
    ylims=c(0,300),
    verbose = verbose ,
    interact= TRUE
  )

  htmlwidgets::saveWidget(
          as_widget(p_html), 'CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_plot/CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_blank.html')

  pdf("CD4-TRM_StimvsCD4-TRM_Unstim_deseq2_plot/CD4-TRM_StimvsCD4-TRM_Unstim_deseq2.pdf") #Change name for CD8
    volplot_wnames(
      resdf,
      pvalth = padjthr,
      lfcth = fcthr,
      pvaltype = "padj",
      lfctype = "log2FoldChange",
      col_feature = means,
      size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
      gene_name = "gene",
      group = "degs",
      check_genes = list(text = features_parse_ensembl(showgenes)),
      return_plot = TRUE,
      clipp = 4,
      xlims=c(-15,15),
      ylims=c(0,300),
      verbose = verbose
    ) + labs(
      size = "Delta %", color = paste0("Mean (", dtype, ")"),
      title = paste(group2, "(-) vs ", group1, "(+)")
    )
  dev.off()

}

### Suplementary table
{ cat(redb("### Crater table and Heatmap  Female_SAvsFemale_MA & Male_MAvsMale_SA --Check tables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  source("/mnt/BioHome/fcastaneda/bin/dgea/R/group_specific.R") #aquiii22

  library(tidyverse)
     fnames <- rev(list.files(
      path = "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/tables/raw", pattern = "_mastlog2cpm_results.csv", full.names = TRUE))

      names(fnames) <- rev(c("Female_SAvsFemale_MA", "Male_MAvsMale_SA")) # needs to be named

     results = lapply(X = fnames, read.csv, row.names = 1, check.names = FALSE)

   # Filter unwanted genes
    #  sc_cd4resting$cluster_sex_disease<-paste0(sc_cd4resting$celltype_subset, "_", sc_cd4resting$orig.asthma, "_", sc_cd4resting$orig.Sex)
    #
    #  sc_cd4resting$Sex_cluster<-paste0(sc_cd4resting$orig.Sex, "_", sc_cd4resting$celltype_subset)
    #
    #  table(sc_cd4resting$Sex_cluster)
    #
    #  sc_cd4resting$keep<-ifelse(sc_cd4resting$celltype_subset %in% c("TRMSP", "TRMDP"), "keep", "no_k")
    #
    #  table(sc_cd4resting$keep, sc_cd4resting$celltype_subset)
    #
    # sc_cd4resting_2<-SplitObject(sc_cd4resting, "keep")
    # sc_cd4resting_keep<-sc_cd4resting_2$keep

      setwd("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/tables/")
      library(data.table)

    sc_cd4stim$orig.sex_asthma<-paste0(sc_cd4stim$orig.Sex, "_", sc_cd4stim$orig.asthma)

     mdata = sc_cd4stim@meta.data
     edata =  expm1(sc_cd4stim@assays$RNA@data)
    #
    gs_report <- group_specific_report(
       results_list = results,
      path = "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/new_figures_paper/stim/tables/",
      mdata = mdata,
      edata = edata,
      hname = "orig.sex_asthma",
      redundant = TRUE,
      sep="&",
      fcthr=0.25 #,
      # check = TRUE
      )

}



#### 

{ cat(redb("### Heatmap pseudo-bulk for TRMSP vs TRMDP Male and Female - Fig3B ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

  source("/mnt/BioHome/fcastaneda/bin/dgea/R/group_specific.R")

  library(tidyverse)
     # fnames <- rev(list.files(
     #  path = "/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/group_specific/info/sex", pattern = "_mastlog2cpm_results.csv", full.names = TRUE))

      fnames<-list(
        "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/Treatment_Responders_TRMDP/Treatment_Responders_TRMDP/No_treatmentvsTreatment_Not_responders/mastlog2cpm_results.csv",

        "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/Treatment_Responders_TRMDP/Treatment_Responders_TRMDP/No_treatmentvsTreatment_Responders/mastlog2cpm_results.csv",

        "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/CD8andCD4_BalBiopsy/results/CD4_Rebuttal/resting/dgea/Treatment_Responders_TRMDP/Treatment_Responders_TRMDP/Treatment_RespondersvsTreatment_Not_responders/mastlog2cpm_results.csv"
      )

      names(fnames) <- c("No_treatmentvsTreatment_Not_responders", "No_treatmentvsTreatment_Responders", "Treatment_RespondersvsTreatment_Not_responders") # needs to be named

     results = lapply(X = fnames, read.csv, row.names = 1, check.names = FALSE)

   # # Filter unwanted genes
   #   sc_cd4resting$cluster_sex_disease<-paste0(sc_cd4resting$celltype_subset, "_", sc_cd4resting$orig.asthma, "_", sc_cd4resting$orig.Sex)
   #
   #   sc_cd4resting$Sex_cluster<-paste0(sc_cd4resting$orig.Sex, "_", sc_cd4resting$celltype_subset)
   #
   #   table(sc_cd4resting$Sex_cluster)



     sc_cd4resting$keep<-ifelse(sc_cd4resting$celltype_subset %in% c("TRMDP") & !is.na(sc_cd4resting$to_dgea_Treatment), "keep", "no_k")

     table(sc_cd4resting$keep, sc_cd4resting$to_dgea_Treatment, useNA="always")

    sc_cd4resting_2<-SplitObject(sc_cd4resting, "keep")
    sc_cd4resting_keep<-sc_cd4resting_2$keep

    table(sc_cd4resting_keep$keep, sc_cd4resting_keep$to_dgea_Treatment, useNA="always")



      setwd(paste0(outdir, "/resting/"))
      library(data.table)
     mdata = sc_cd4resting_keep@meta.data
     edata =  expm1(sc_cd4resting_keep@assays$RNA@data)


        results_list = results
        path = paste0(outdir, "/resting/Treatment_Responders_27July")
        colours = NULL
        return_report = TRUE
        optimize = FALSE
        verbose = TRUE
        # to calculate the stats
        mdata = mdata
        edata = edata
        hname = "to_dgea_Treatment"
        redundant = TRUE
        sep="&"
        fcthr=0.25 #,
        padjthr = 0.05
        vs = "vs"

        if(verbose) cat("\n%%%%%%%%%%%%%%%%%%%% GSpec report %%%%%%%%%%%%%%%%%%%%\n")
        if(!is.null(path) && verbose) cat("Report at:\n", path, "\n")
        the_report = list()
        results_group_specific = group_specific_features(
          results = results_list,  redundant = TRUE,  sep="&", fcthr=0.25, padjthr = 0.05, verbose = TRUE)

            table(is.na(results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Not_responders)'))

            table(is.na(results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Responders)'))

            table(is.na(results_group_specific$group_specific$'log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)'))

            #
            results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Responders)'[is.na(results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Responders)')] <- 0
            #
            results_group_specific$group_specific$'log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)'[is.na(results_group_specific$group_specific$'log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)')] <- 0

            table(is.na(results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Not_responders)'))

            table(is.na(results_group_specific$group_specific$'log2FoldChange(No_treatmentvsTreatment_Responders)'))

            table(is.na(results_group_specific$group_specific$'log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)')) ##ASK SARA

            # table(is.na(results_group_specific$group_specific$'log2FoldChange(Female_TRMDPvsFemale_TRMSP)'))  What for??

          # check _summary
        # if(isTRUE(list(...)$check)) return(results_group_specific)
        if(verbose > 1) str(results_group_specific)
        summ_df = results_group_specific$group_specific
        if(is.null(results_group_specific$stats)){
          if (verbose) cat("Summary statistics fetching or calculation\n")
          stat_df = data.frame(
            row.names = unique(unlist(lapply(results_list, rownames))))
          stat_df$gene_name = rownames(stat_df)
          tmp <- paste0(levels(summ_df$group),
            rep(moments[c("mn", "p")], each = nlevels(summ_df$group)))
          for(i in tmp){
            for(j in 1:length(results_list)){
              y <- which(grepl(paste0("^", i), colnames(results_list[[j]])))
              if(any(grepl(paste0("^", i), colnames(stat_df)))) next
              if(length(y))
                stat_df = joindf(stat_df, results_list[[j]][, y, drop = FALSE])
            }
          }
          if(!is.null(mdata) && !is.null(edata)){
            results_group_specific$stats = stats_summary_table(
              mat = edata,
              groups = setNames(mdata[, hname], rownames(mdata)),
              rname = rownames(results_group_specific[[1]]),
              moments = c("mn", "p"),
              verbose = verbose
            )
          }else if(ncol(stat_df) > 1){
            if(verbose) cat("------- Stat summary -------\nTaken from results\n")
            if(verbose) warning("Some measurements may be missing from comparisons\n")
            colnames(stat_df) <- gsub(".eurat.*|CPM.*", "", colnames(stat_df))
            results_group_specific$stats = stat_df[, -1, drop = FALSE]; rm(stat_df)
          }
        }

        if(!is.null(results_group_specific$stats)){
          summ_df = joindf(summ_df, results_group_specific$stats)
        }
        # summ_df <- summ_df[order(summ_df$group), ]
        # groups = levels(summ_df$group)[levels(summ_df$group) %in% colnames(summ_df)]
        # tvar <- colnames(summ_df)[!colnames(summ_df) %in% groups]
        # summ_df <- summ_df[, c(tvar, groups)]
        the_report$summary = summ_df
        # if(grepl("summary", return_report)) return(the_report)
        # if(!isTRUE(return_report))
          write.csv(summ_df, file = paste0(path, "_summary_Treatment_Responders.csv"), row.names = FALSE)

        if(!is.null(results_group_specific$stats)){
          if(verbose) cat("--- Heatmap\n")
          # annor = summ_df[summ_df[, "group"] != "NDE", "group", drop = FALSE]
          # annor[, "group"] <- droplevels(annor[, "group"])

          annor_comb = summ_df[summ_df[, "combined"] != "NDE", c("combined", "log2FoldChange(No_treatmentvsTreatment_Not_responders)", "log2FoldChange(No_treatmentvsTreatment_Responders)", "log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)"), drop = FALSE]

          annor_comb$avg_FCs<-rowMeans(annor_comb[,c( "log2FoldChange(No_treatmentvsTreatment_Not_responders)", "log2FoldChange(No_treatmentvsTreatment_Responders)", "log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)")])

          annor_comb[, "combined"] <- droplevels(factor(annor_comb[, "combined"]))
          table(annor_comb[, "combined"])

          annor_comb$shared_level<-str_count(annor_comb$combined, "&") + 1
          annor_comb$shared_level<-ifelse(annor_comb$combined %in% c("Treatment_Responders", "Treatment_Responders&Treatment_Responders"), 0.5, annor_comb$shared_level)

          annor_comb$shared_level<-ifelse(annor_comb$combined %in% c("Treatment_Not_responders", "Treatment_Not_responders&Treatment_Not_responders"), 0.1, annor_comb$shared_level)

          annor_comb$shared_level<-ifelse(annor_comb$combined %in% c("No_treatment", "No_treatment&No_treatment"), 1, annor_comb$shared_level)

          table(annor_comb$combined, annor_comb$shared_level)

          # str_count(annor_comb$combined, "&") + 1

          table(is.na(annor_comb$avg_FCs))
          # annor_comb$

          annor_comb<-annor_comb[order(annor_comb$shared_level, annor_comb$avg_FCs),]

          annor_comb$shared_level <- NULL
          # annor_comb$shared_level2 <- NULL
          annor_comb$'log2FoldChange(No_treatmentvsTreatment_Not_responders)' <- NULL
          annor_comb$'log2FoldChange(No_treatmentvsTreatment_Responders)' <- NULL
          annor_comb$'log2FoldChange(Treatment_RespondersvsTreatment_Not_responders)' <- NULL
          annor_comb$avg_FCs <- NULL

          columns <- grep("_percent", colnames(summ_df), value = TRUE)
          if(length(columns) == 0){ str(summ_df); stop("No stats columns found") }
          tvar <- apply( # keeping features with a % difference > 25
            X = summ_df[rownames(annor_comb), columns],
            MARGIN = 1, FUN = function(x) abs(diff(range(x)))
          )
          labels_row_i <- rownames(annor_comb)
           labels_row_i[!rownames(annor_comb) %in% names(tvar[tvar > 20])] <- ""
          the_report$heatmap = list()
          for(i in c("mn", "p")){
             i<-"mn"
            moments_i = moments[i]
            if(verbose) cat("    *", moments_i, "\n")
            anno_cnames = grepl(paste0(moments_i, "$"), colnames(summ_df))
            mat2plot <- as.matrix(summ_df[rev(rownames(annor_comb)), anno_cnames])
            if(moments_i[[1]] == "_mean") mat2plot <- log2(mat2plot + 1)
            if(1){ # row-wise scale and get the range
              mat2plot <- t(scale(t(mat2plot)))
              topz <- max(c(min(abs(c(range(mat2plot), 2))), 1))
            }else{ topz <- 75 } # this was for percentages...
            mat2plot[mat2plot > topz] <- topz; mat2plot[mat2plot < (-topz)] <- -topz;
            colnames(mat2plot) <- sub(moments_i, "", colnames(mat2plot))
            tvar <- intersect(levels(annor_comb[, "combined"]), colnames(mat2plot))
            mat2plot <- mat2plot[, tvar]


            annoc = data.table(
              Group = c(colnames(mat2plot)[c(1,3,2)]))

            # annoc<-data.frame(annoc[, c("Sex", "cluster") := tstrsplit(annoc$Group, "_", fixed=TRUE)])

            # rownames(annoc) <- c(colnames(mat2plot)[c(2,4,1,3)])
            # annoc<-annoc[,c("cluster", "Sex")]

            rownames(annoc) <- annoc$Group
            annoc$Group2<-annoc$Group

            colors_i = c(as.list(annor_comb), list(Group = colnames(mat2plot)))
            tmp <- v2cols(unname(unlist(colors_i)), colours)
            colors_i <- lapply(colors_i, v2cols, sour = tmp)
            # colors_i$cluster<-c('TRMDP'="#99141B", 'TRMSP'="#FFA303")
            # colors_i$disease<-c('MA'="#160cc9", 'SA'="#f55236")
            # colors_i$sex<-c('Female'="#d3a0eb", 'Male'="#86bf73")

            # if(!isTRUE(return_report)){
              fname <- paste0(path, "_heatmap", moments_i)
              # if(optimize){
              #   png(paste0(fname, ".png"), width = 1500, height = 1700, res = 250)
              # }else{
            # pdf(paste0(fname, "_legend.pdf"), width = 7, height = 12) #}
            pdf(paste0(fname, "Treatment_responders_corrected2_legend.pdf"), width = 7, height = 12)
            # pdf(paste0(fname, "_corrected2.pdf"), width = 7, height = 12)

            # }

            cols <- grDevices::colorRampPalette(rev(c('yellow', 'black', 'blue')))(256)

            # annor = summ_df[summ_df[, "group"] != "NDE", "group", drop = FALSE]
            # annor[, "group"] <- droplevels(annor[, "group"])

            head(mat2plot[,annoc$Group])
            mat2plot<-mat2plot[rownames(annor_comb),rownames(annoc)]

            head(mat2plot)

            pheatmap(
              mat = mat2plot,
              color = cols, border_color = NA, scale = "none",
              cluster_cols = FALSE, cluster_rows = FALSE,
              annotation_row = annor_comb, #annotation_col = annoc,
              annotation_colors = colors_i,
              fontsize_row = 10/(4*log10(nrow(mat2plot))),
              annotation_names_col = FALSE, annotation_names_row = FALSE,
              show_colnames = FALSE, show_rownames = TRUE,
              gaps_col = cumsum(table(annoc$Group)),
               # gaps_row = cumsum(table(annor_comb[, "combined"])),
              labels_row = labels_row_i, annotation_legend= FALSE #FALSE
            ); graphics.off()
              write.csv(annor_comb, "order_genes_Treatment_responders_corrected2_legend.csv")

          }
          if(grepl(paste0("heatmap"), return_report)) return(the_report)
        }
}
