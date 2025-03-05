

rm(list = ls()) # clears objects from the workspace
library(dplyr)
library(rRACES)



mutational_signature_getter <- function(path, tool) {
  if (tool == "SigProfiler") {
    # SigProfiler
    file_paths <- c(
      COSMIC_exp = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt"), 
      COSMIC_sig = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt"), 
      DENOVO_exp = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"), 
      DENOVO_sig = paste0(path, "SigProfiler/results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt")
    )
    output <- lapply(file_paths, read.delim)
    names(output) <- basename(file_paths)
    return(output)
  } else if (tool == "SparseSignatures") {
    # SparseSignatures
    output <- readRDS(paste0(path, "SparseSignatures/MSeq_nmf_Lasso_out.rds"))
    return(output)
  } else {
    stop("Invalid!")
  }
}

PATH <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPN_ID <- 1
coverage <- 100
purity <- 0.6
variant_caller <- "XXX"
pipeline <- "signature_deconvolution"
tool <- "SigProfiler"


MAIN_PATH <- paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", "tumourevo/", 
                    coverage, "x_", purity, "p/", variant_caller, "_ascat/", 
                    pipeline, "/")


file_paths <- filepath_getter(MAIN_PATH, tool)




  













tools <- list(variant_caller=c("A"), 
              signature_deconvolution=c("SparseSignatures", "SigProfiler"), 
              subclonal_deconvolution=c("viber", "ctree", "pyclone", "mobster"))

p <- "signature_deconvolution"
for (tool in tools[[p]]) {
  print(tool)
}

tool <- tools[[p]][1]

paste0(PATH, "SPN0", SPN_ID, "/tumourevo/", coverage, "x_", purity, "p/", "variant_caller_ascat", "/", tool, "Suggested_Solution")


tool <- "SigProfiler"
context <- 96
signature <- "COSMIC"

sig_type <- ifelse(signature == "COSMIC", "COSMIC_SBS96_Decomposed_Solution", "SBS96_De-Novo_Solution")

ms_path <- paste0("MSeq/signature_deconvolution/", tool, "/results/SBS", context, "/Suggested_Solution/", sig_type)
exp <- paste0(ms_path, "/Activities/COSMIC_SBS96_Activities.txt")
sigs <- paste0(ms_path, "/Signatures")














