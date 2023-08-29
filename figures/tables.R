#!/usr/bin/R

#####################
# Tables compendium #
#####################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-05-24
# ---

# This script will format supplementary tables.

source("/home/ciro/scripts/handy_functions/devel/supp_table.R")
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
source("/home/ciro/scripts/handy_functions/devel/utilities.R")
source("/home/ciro/scripts/handy_functions/devel/filters.R")
fig_dir = '/home/ciro/large/asthma_airways/results/figures_cd4'
dir.create(fig_dir, showWarnings = FALSE); setwd(fig_dir)
cat("Working at (fig_dir):", fig_dir, "\n")
system("ls -loh")
dir.create("supplementary_tables", showWarnings = FALSE)

{ cat("### Single-cell QC ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  # Rscript /home/ciro/scripts/cellranger_wrappeR/summary.R -i \
  #   /mnt/BioAdHoc/Groups/vd-vijay/cramirez/hayley/raw/sever_asthma3/COUNTS_hg19
  # grep -E "libraries_summary|Q30|Mapped" ~/*/scripts/*/*
  cr_path = "/home/ciro/ad_hoc/hayley/raw/sever_asthma3"
  outs = c("CD4US_filt_mapped", "ovCD4ST_aggr_mapped")
  library_mdata_exclude = "CD8|library_type|aggr"

  aggr_file = paste0(cr_path, "/COUNTS_hg19/AGGR/", outs, "/outs/aggregation.csv")
  library_mdata_f = "/home/ciro/asthma_airways/info/metadata_library.csv"
  library_mdata = read.csv(library_mdata_f, stringsAsFactors = FALSE, row.names = 1)
  aggr_libs = unlist(lapply(
    X = aggr_file,
    FUN = function(x) read.csv(x, stringsAsFactors = FALSE)[, 1]
  ))
  library_stats_f = paste0(cr_path, "/COUNTS_hg19/SingleCell3v3_libraries_summary.csv")
  library_stats = read.csv(library_stats_f, stringsAsFactors = FALSE, check.names = FALSE)

  merges = c("_pre_normalization_raw_reads_per_filtered_bc",
    "_pre_normalization_cmb_reads_per_filtered_bc")
  aggr_report_list = lapply(
    X = paste0(dirname(aggr_file), "/summary.json"),
    FUN = function(fname){
      aggr_report = rjson::fromJSON(paste(readLines(fname, warn=F), collapse=""))
      for(i in merges) aggr_report[[i]] = unlist(aggr_report[paste0(aggr_report$batches, i)])
      tvar <- sapply(aggr_report, length) # added more vectors of libraries
      setNames(list(
        data.frame(aggr_report[tvar == max(tvar)]),
        aggr_report[grep("multi_transcrip|total_reads", names(aggr_report))]
      ), c("table", "aggr"))
  }); names(aggr_report_list) = basename(gsub("outs.*", "", aggr_file))
  aggr_report = data.frame(data.table::rbindlist(lapply(aggr_report_list, "[[", 1)))
  aggr_report$frac_reads_kept <- scales::percent(aggr_report$frac_reads_kept)
  rownames(aggr_report) <- aggr_report[, 1]
  aggr_report <- aggr_report[, -c(1, 3)] # don't need No. of cells per library
  colnames(aggr_report) <- paste("Aggregation", c(
    "Fraction of Reads Kept", "Pre-Normalization Total Reads per Cell",
    "Pre-Normalization Confidently Mapped Barcoded Reads per Cell"
  ))

  if(1){
    cat("All libraries in stats?", all(aggr_libs %in% library_stats$Library), "\n")
    cat("All libraries in aggr?", all(aggr_libs %in% rownames(aggr_report)), "\n")
    cat("All libraries in metadata?", all(aggr_libs %in% rownames(library_mdata)), "\n")
  }

  lib_stats = library_stats[library_stats$Library %in% aggr_libs, ]
  tvar <- !grepl(library_mdata_exclude, colnames(library_mdata))
  lib_stats <- cbind(
    Library = lib_stats$Library,
    library_mdata[lib_stats$Library, tvar],
    lib_stats[, -1],
    aggr_report[lib_stats$Library, ]
  )
  lib_stats$Library <- as.character(lib_stats$Library)
  str(lib_stats)
  lib_stats[, (ncol(lib_stats) - 1):ncol(lib_stats)]

  st1 = supp_table(
    mytables = list("Gene Expression." = lib_stats),
    headers = list(
      none = c("Library" = "TEXT"),
      Metadata = c(setNames("TEXT", paste0(colnames(library_mdata)[1:5], collapse = "|"))),
      "Cells stats" = c(
        "Number|per|Detected" = "COMMA", Barcodes = "PERCENTAGE",
        Fraction = "PERCENTAGE", exclude = "Aggregation"),
      "Reads Mapped (hg19-3.0.0)" = c("Reads Mapped" = "PERCENTAGE"),
      "Quality" = c("Q30" = "PERCENTAGE", "Saturation" = "PERCENTAGE"),
      "Aggregation" = c("^aggr" = "TEXT", "Kept" = "PERCENTAGE", "Normalization" = "COMMA")
    ), title_name = "Single-cell sequencing quality control.",
    rename_columns = "Reads Mapped |Aggregation |aggr.", body_start_n = 68
  )
  fname = "supplementary_tables/single-cell_qc.xlsx" # openxlsx::getSheetNames(fname)

  # Add further aggregation info
  tmp = lapply(aggr_report_list, function(x){
    tvar <- lapply(x[2], function(y) grepl("total_reads", names(y)) )
    x[[2]][which(tvar[[1]])]
  }); tmp = reshape2::melt(tmp)[, 3:1]
  colnames(tmp) <- tail(colnames(lib_stats), 3)
  tmp[, 1] <- gsub("_$|^_", "", gsub("_{2,}", "_", gsub("ov|aggr|filt|mapped", "", tmp[, 1])))
  tmp[, 2] <- stringr::str_to_title(gsub("_", " ", tmp[, 2]))
  openxlsx::writeData(wb = st1, sheet = st1$sheet_names[1], x = tmp,
    startCol = ncol(lib_stats) - (ncol(tmp) + 1),
    startRow = nrow(lib_stats) + 5, colNames = FALSE)
  openxlsx::setColWidths(wb = st1, sheet = st1$sheet_names[1], cols = 1, widths = 25)
  # formats = list(
  #   list(format = "TEXT", header = ncol(lib_stats) - (ncol(tmp) + 1)),
  #   list(format = "TEXT", header = ncol(lib_stats) - (ncol(tmp))),
  #   list(format = "COMMA", header = ncol(lib_stats) - (ncol(tmp) - 1))
  # )
  # st1 = .body_style_fun(wfile = st1, sheet = st1$sheet_names[1],
  #   row_start = nrow(lib_stats) + 5, length_i = nrow(tmp),
  #   formats = formats, verbose = 2)
  # addStyle(wb = st1, sheet = st1$sheet_names[1], style = .body_section_style(),
  #   rows = (nrow(lib_stats) + 5):(nrow(lib_stats) + 5 + nrow(tmp)),
  #   cols = (ncol(lib_stats) - ncol(tmp)):ncol(lib_stats),
  #   gridExpand = TRUE, stack = TRUE)
  openxlsx::saveWorkbook(wb = st1, file = fname, overwrite = TRUE)
}

{ cat("### Clinical correlation ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cor_f = c(
    "/home/ciro/large/asthma_pjs/results/associations/cd4_male/summary_table_20792tests.csv",
    "/home/ciro/large/asthma_pjs/results/associations/cd4_female/summary_table_20048tests.csv")
  names(cor_f) <- gsub("cd4_", "", basename(dirname(cor_f)))
  axis_y = c(
    "Adapted_Asthma_Severity_Score",
    "ACQ_6_Score",
    "GINA_Management_Step",
    "No_of_OCS_past_12_mths_at_NIH_EOSA_enrolment",
    "No_of_Hospital_Admissions_past_12_mths",
    "Pre_BD_FEV1_pred", "Post_BD_FEV1_pred",
    "Pre_BD_FVC_pred", "Post_BD_FVC_pred",
    "Pre_BD_FEV1_FVC_ratio", "Post_BD_FEV1_FVC_ratio",
    "FeNO_ppb_at_study_enrolment",
    "Blood_Eosinophils_Count_10_9_L_at_Bronchoscopy",
    "Blood_Eosinophils_Count_Current_10_9_L_at_enrolment_3_mths",
    "Age_at_Bronchoscopy_yrs",
    "Age_of_Asthma_Diagnosis_yrs",
    "BMI_kg_m2"
  )
  axis_x = paste0("resting_CD4_CL", c(0:5, 7:8))
  # head(grep("resting_", cor_df$group2, value = TRUE))
  cor_df = data.table::rbindlist(lapply(X = names(cor_f), FUN = function(x){
    y <- readfile(cor_f[[x]], stringsAsFactors = FALSE);
    y$padj = p.adjust(y$pval, method = "BH"); y$Sex = x; y
  }))
  tvar <- (cor_df$group1 %in% c(axis_y, axis_x) | grepl("resting_CD4", cor_df$group1)) &
    (cor_df$group2 %in% c(axis_y, axis_x) | grepl("resting_CD4", cor_df$group2))
  cor_df_subset = data.frame(cor_df[which(tvar), ], stringsAsFactors = FALSE)
  cor_df_subset$test = NULL; #cor_df_subset$Type = "Correlation"
  cor_df_subset$Name = NULL
  cor_df_subset$Sex = stringr::str_to_sentence(cor_df_subset$Sex)
  cor_df_subset$method = stringr::str_to_sentence(cor_df_subset$method)
  tvar <- c("Group_1", "Group_2", "Size", "P-value", "Adjusted P-value", "Method")
  colnames(cor_df_subset)[1:6] <- tvar
  # names(table(cor_df_subset$group1)) %in% axis_y
  # names(table(cor_df_subset$group2)) %in% axis_y; table(cor_df_subset$Sex)
  tmp = c("Pre_BD_FEV1_FVC_ratio", "Post_BD_FEV1_FVC_ratio")
  tvar <- cor_df_subset$Group_1 %in% tmp | cor_df_subset$Group_2 %in% tmp
  cor_df_subset$Size[tvar] <- cor_df_subset$Size[tvar] * -1

  supp_table(
    mytables = list("Associations." = cor_df_subset),
    headers = list(
      Groups = c(Sex = "TEXT", Group = "TEXT"),
      Tests = c(Method = "TEXT", "Size|value" = "NUMBER")
    ), title_name = "Clinical features associated with CD4 T cell subsets.",
    filename = "supplementary_tables/associations"
  )
}

{ cat("### Single-cell resting DGEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  marknames = paste0(
    "/home/ciro/large/asthma_airways/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/",
    "markers/20PCs_RNA_snn_res.0.4_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
  )
  mygenes <- read.csv(marknames, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  dgea_dir = "/home/ciro/large/asthma_airways/results/scdgea"
  trm_dgea_files = setNames(c(
    paste0(dgea_dir, "/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters/0vs1/results_0vs1_mastlog2cpm.csv"),
    paste0(dgea_dir, "/airways_cd4/comprs/activation/STvsUS/results_STvsUS_mastlog2cpm.csv"),
    # paste0(dgea_dir, "/stim_cd4_15p_sng_nomacro9n10/comprs/stim_disease_cov_sex/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv")
    paste0(dgea_dir, "/airways_cd4/comprs/stim_disease/ST_SAvsST_MA/results_ST_SAvsST_MA_mastlog2cpm.csv")
  ), c("TRM103+ vs TRM103-", "Activation", "Activation disease"))
  trm_dgea_list = lapply(X = trm_dgea_files, FUN = function(x){
    trm_dgea <- read.csv(x, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
    trm_dgea$gene = gsub("'", "", trm_dgea$gene_name)
    trm_dgea
  })
  sexes = c("Female", "Male")
  clusters = rep(0:1, each = 2)
  trm_dgea0n1_files = setNames(paste0(dgea_dir, "/resting_cd4_15p_sng_nomacro_x7n9/comprs/cluster",
    clusters, "/", sexes, "_SAvs", sexes, "_MA/results_", sexes, "_SAvs", sexes,
    "_MA_mastlog2cpm.csv"), paste0(sexes, "_", clusters))
  trm_dgea0n1_list = lapply(X = names(trm_dgea0n1_files), FUN = function(x){
    fname = trm_dgea0n1_files[[x]]
    trm_dgea <- read.csv(fname, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
    trm_dgea$gene = gsub("'", "", trm_dgea$gene_name)
    colnames(trm_dgea) <- gsub("Female_|Male_", "", colnames(trm_dgea))
    colnames(trm_dgea)[-ncol(trm_dgea)] <- paste0(x, "_", colnames(trm_dgea)[-ncol(trm_dgea)])
    trm_dgea
  })
  trm_dgea0n1 = trm_dgea0n1_list[[1]]
  for(i in 2:length(trm_dgea0n1_list))
    trm_dgea0n1 = dplyr::left_join(trm_dgea0n1, trm_dgea0n1_list[[i]], by = "gene")
  supp_table(
    mytables = c(list(
      `Resting CD4` = mygenes,
      `Sex-specific differences` = trm_dgea0n1), trm_dgea_list),
    headers = list(
      none = c("gene$" = "TEXT"),
      "Significance" = c("^cluster|group$" = "TEXT", "log.*F" = "NUMBER",
        "p_val$|pvalue$" = "SCIENTIFIC", "adj$" = "SCIENTIFIC", "pct\\." = "NUMBER",
        "sCluster" = "TEXT"),
      "Mean expression (CPM)" = c("meanC" = "NUMBER"),
      "Fraction of expressing cells (CPM > 0)" = c("exprFrac|_percentage" = "NUMBER"),
      "Mean of expressing cells (CPM > 0)" = c("exprMean" = "NUMBER")
    ),
    title_name = "Single-cell differential gene expression analysis.",
    filename = "supplementary_tables/single-cell_resting_markers",
    verbose = 1
  )
}

{ cat("### Per cluster sex-disease DGEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dgea_f = list.files(path = "craters", pattern = "stat", full.names = TRUE)
  names(dgea_f) = gsub("_.*", "", basename(dgea_f))
  dgea_f <- dgea_f[file.exists(dgea_f)]; tvar <- grepl("resting", dgea_f)
  dgea_f <- c(dgea_f[tvar], dgea_f[!tvar])
  names(dgea_f) <- gsub("cluster", "", names(dgea_f))
  dgea_list = lapply(X = dgea_f, FUN = function(x){
    y <- readfile(x, stringsAsFactors = FALSE, row.names = 1)
    y$filters <- !y$filters; colnames(y) <- gsub("filters", "DEG", colnames(y))
    colnames(y)[2:3] <- paste0(colnames(y)[2:3], ".log2FC"); y
  })
  names(dgea_list[[1]]) <- stringr::str_replace_all(
    names(dgea_list[[1]]), c("log2FC" = "L2FC", "\\.padj" = ".adjp"))
  dgea_df = data.frame(row.names = unique(unlist(lapply(dgea_list, rownames))))
  dgea_df = joindf(dgea_df, dgea_list[[1]])
  dgea_df = dgea_list[[1]]
  colnames(dgea_df) <- paste0(names(dgea_list)[1], ".", colnames(dgea_df))
  for (i in c(2:length(dgea_list))){
    tmp = dgea_list[[i]]
    colnames(tmp) <- paste0(names(dgea_list)[i], ".", colnames(tmp))
    dgea_df = joindf(dgea_df, tmp); rm(tmp)
  }
  dgea_df$gene = rownames(dgea_df)
  tvar <- !grepl("gene_name|group|lfcth|signi|min_|log2_|\\.mean$", colnames(dgea_df))
  tmp <- which(rowSums(dgea_df[, grepl("DEG", colnames(dgea_df))], na.rm = TRUE) > 0)
  dgea_df <- dgea_df[tmp, tvar]
  names(dgea_df) <- gsub("resting.DEG", "resting.DEGs", names(dgea_df))
  str(dgea_df, list.len = 20)

  supp_table(
    mytables = list(`Per cluster sex-disease` = dgea_df),
    headers = list(
      none = c("gene$" = "TEXT"),
      "All resting cells" = c("L2FC" = "NUMBER", "adjp$" = "SCIENTIFIC", DEGs = "TEXT"),
      "Significant" = c("^cluster|group|DEG$" = "TEXT"),
      "Significance" = c("log.*F" = "NUMBER", "adj$" = "SCIENTIFIC"),
      "Mean expression" = c("_mean" = "NUMBER"),
      "Fraction of expressing cells" = c("_exprFrac|_percentage" = "NUMBER")
    ), body_start_n = 26,
    title_name = "Single-cell differential gene expression analysis.",
    filename = "supplementary_tables/single-cell_per_cluster_sex-disease_all"
  )
}

{ cat("### Per cluster disease DGEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  dgea_f = list.files(
    "/home/ciro/large/asthma_airways/results/scdgea/resting_cd4_15p_sng_nomacro_x7n9/comprs/clusters_disease",
    pattern = "results.*.csv", full.names = TRUE, recursive = TRUE
  ); names(dgea_f) = basename(dirname(dgea_f))
  names(dgea_f) <- gsub("_.*", "", names(dgea_f))
  dgea_list = lapply(X = dgea_f, FUN = function(x){
    y <- readfile(x, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    y$DEG <- rownames(y) %in% filters_thresholds(y, fc = 0.05, verbose = TRUE); y
  })
  dgea_df = data.frame(row.names = unique(unlist(lapply(dgea_list, rownames))))
  dgea_df = joindf(dgea_df, dgea_list[[1]])
  colnames(dgea_df) <- gsub(paste0(names(dgea_list)[1], "_"), "", colnames(dgea_df))
  colnames(dgea_df) <- paste0(names(dgea_list)[1], ".", colnames(dgea_df))
  for (i in c(2:length(dgea_list))){
    tmp = dgea_list[[i]]
    colnames(tmp) <- gsub(paste0(names(dgea_list)[i], "_"), "", colnames(tmp))
    colnames(tmp) <- paste0(names(dgea_list)[i], ".", colnames(tmp))
    dgea_df = joindf(dgea_df, tmp); rm(tmp)
  }
  dgea_df$gene = rownames(dgea_df)
  tvar <- !grepl("gene_name|group|pcor|pvalue|diff|min", colnames(dgea_df))
  tmp <- which(rowSums(dgea_df[, grepl("DEG", colnames(dgea_df))], na.rm = TRUE) > 0)
  dgea_df <- dgea_df[tmp, tvar]
  str(dgea_df, list.len = 20)
  head(sort(table(dgea_df$gene)))

  supp_table(
    mytables = list(`Per cluster sex-disease` = dgea_df),
    headers = list(
      none = c("gene$" = "TEXT"),
      "Significant" = c("^cluster|group|DEG$" = "TEXT"),
      "Significance" = c("log.*F" = "NUMBER", "adj$" = "SCIENTIFIC"),
      "Mean expression" = c("_mean" = "NUMBER"),
      "Fraction of expressing cells" = c("_exprFrac|_percentage" = "NUMBER")
    ), body_start_n = 26,
    title_name = "Single-cell differential gene expression analysis.",
    filename = "supplementary_tables/single-cell_per_cluster_disease"
  )
}
