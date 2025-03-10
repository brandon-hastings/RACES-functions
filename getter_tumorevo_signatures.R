


# quality control the input parameters
#-------------------------------------------------------------------------------
qc_args <- function(name, coverage, purity, tools) {
  if (!is.numeric(name) || name %% 1 != 0) {
    stop("Error: 'name' must be integer") # name : integer
  }
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!(tools %in% c("SigProfiler", "SparseSignatures"))) {
    stop("Error: 'type' must be one of 'SigProfiler' and 'SparseSignatures'")
  }
  return(TRUE)
}


mutational_signature_getter <- function(MAIN_PATH, SPN_ID, coverage, purity, tools) {
  
  qc_args(SPN_ID, coverage, purity, tools)
  
  path <- paste0(MAIN_PATH, "SPN0", SPN_ID, "/tumourevo/", coverage, "x_", purity, "p/results_tumourevo_mseq/MSeq/signature_deconvolution/")
  
  object_list <- list()
  for (tool in tools) {
    object_list[[tool]] <- mutational_signature_getter_aux(path, tool)
  }
  return(object_list)
}



# input
  # path: absolute path
  # tool: c('SigProfiler', 'SparseSignatures')
mutational_signature_getter_aux <- function(path, tool) {
  
  # SigProfiler
  if (tool == "SigProfiler") {
    file_paths <- c(
      COSMIC_exp = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt"), 
      COSMIC_sig = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt"), 
      DENOVO_exp = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"), 
      DENOVO_sig = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt")
    )
    output <- lapply(file_paths, read.delim)
    names(output) <- basename(file_paths)
    return(output)
  # SparseSignatures
  } else if (tool == "SparseSignatures") {
    output <- readRDS(paste0(path, "SparseSignatures/MSeq_nmf_Lasso_out.rds"))
    return(output)
  } else {
    stop("Invalid!")
  }
}





