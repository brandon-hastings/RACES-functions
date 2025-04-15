

#===============================================================================
# SigProfiler getter
#===============================================================================
sigprofiler_getter <- function(dataset, context) {
  
  if (context == "SBS") {
    ctx = "SBS96"
  } else if (context == "DBS") {
    ctx = "DBS78"
  } else {
    stop("Wrong context!")
  }
  
  p1 <- paste0("reuslts_tumorevo_", tolower(dataset), "/")
  
  p2 <- paste0(p1, dataset, "/signature_deconvolution/SigProfiler/results/", ctx, "/Suggested_Solution/")
  
  output <- list(
    COSMIC_EXPOSURE_PATH=paste0(p1, "COSMIC_", ctx, "_Decomposed_Solution/Activities/COSMIC_", ctx, "_Activities.txt"), 
    COSMIC_SIGNATURES_PATH=paste0(p1, "COSMIC_", ctx, "_Decomposed_Solution/Signatures/COSMIC_", ctx, "_Signatures.txt"), 
    DENOVO_EXPOSURE_PATH=paste0(p1, ctx, "_De-Novo_Solution/Activities/", ctx, "_De-Novo_Activities_refit.txt"), 
    DENOVO_SIGNATURE_PATH=paste0(p1, ctx, "_De-Novo_Solution/Signatures/", ctx, "_De-Novo_Signatures.txt")
  )
  return(output)
}


#===============================================================================
# SparseSignatures getter
#===============================================================================
sparsesignatures_getter <- function(dataset) {
  
  p1 <- paste0("reuslts_tumorevo_", tolower(dataset), "/")
  
  p2 <- paste0(p1, dataset, "/signature_deconvolution/SparseSignatures/")
  
  output <- list(
    CV=paste0(p1, "MSeq_cv_means_mse.rds"), 
    NMF=paste0(p1, "MSeq_nmf_Lasso_out.rds")
  )
  return(output)
}


#===============================================================================
# signature deconvolution getter
#===============================================================================
signature_getter <- function(tool, dataset, context) {
  
  if (tolower(tool) == "sparsesignature") {
    output <- sparsesignatures_getter(dataset)
  } else if (tolower(tool) == "sigprofiler") {
    output <- sigprofiler_getter(dataset, context)
  } else {
    stop("wrong tool name specified!")
  }
  return(output)
}




