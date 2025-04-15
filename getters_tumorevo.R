
source("getters_tumorevo_signatures.R")
source("getters_tumorevo_subclone.R")

tumorevo_getter <- function(PATH, SPN, coverage, purity, subclone_tool, signature_tool, dataset, context, patient_id=NULL, sample_id=NULL) {
  
  main_path <- paste0(PATH, SPN, "/tumourevo/", coverage, "x_", purity, "/")
  
  a <- signature_getter(signature_tool, dataset, context)
  b <- subclone_getter(subclone_tool, dataset, patient_id, sample_id)
  aa <- paste0(main_path, a)
  bb <- paste0(main_path, b)
  cc <- c(aa, bb)
  return(cc)
}

