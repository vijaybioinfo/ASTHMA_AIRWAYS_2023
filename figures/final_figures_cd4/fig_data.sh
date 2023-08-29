#!/usr/bin/bash

#########################
# Figure data gathering #
#########################

# This script will place all the necessary files for a figures folder

OUTDIR=/home/ciro/large/asthma_pjs/results/a1_final_figures_cd4
FIGDATA=/home/ciro/asthma_pjs/scripts/final_figures_cd4/fig_data.csv
UNLINK1ST=TRUE

. /home/ciro/scripts/handy_functions/devel/fig_data.sh

# cd /home/ciro/large/asthma_pjs/results/a1_final_figures_cd4
# ### In R ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## Gene sets
# load("/home/ciro/scripts/handy_functions/data/signatures_vijaylab.rdata")
# load("/home/ciro/scripts/handy_functions/data/signatures_cd8covid.rdata")
# ssigna_patterns = "tfh|cytotox|th17|interf|chmied|cycl|saragreg|bulk|th1|th2|tcm|acti|hombr|clar"
# signames <- c(
#   "/home/ciro/asthma_pjs/info/deprecated/signatures_hdmplus.csv",
#   "/home/ciro/covid19/info/tcell_activation_signatures.csv",
#   "/home/ciro/asthma_pjs/info/airways_bulk_resting_trmness_TRMDPvsTEFF_padj0.05_fc1_curated.csv"
# )
# tvar <- lapply(signames, read.csv, stringsAsFactors = FALSE)
# tvar <- unlist(lapply(tvar, as.list), recursive = FALSE)
# globalsign <- signatures_vijaylab
# globalsign <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
# globalsign <- sapply(globalsign, function(x) x[!(x == "" | is.na(x))] )
# globalsign <- c(globalsign, signatures_cd8covid[!names(signatures_cd8covid) %in% names(globalsign)])
# str(globalsign)
# tvar <- list(
#   trm_signature_saragreg = c("ITGAE", "ITGA1", "CXCR6", "AMICA1", "ALOX5AP", "CLU", "MYO7A"),
#   nontrm_signature_saragreg = c("S1PR1", "SELL", "KLF2", "KLF3", "TCF7", "MAL")
# )
# globalsign <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
# tvar <- grep(ssigna_patterns, names(globalsign), value = TRUE)
# global_signatures <- globalsign[tvar]
# str(global_signatures)
# save(global_signatures, file = "../data/global_signatures.rdata")
