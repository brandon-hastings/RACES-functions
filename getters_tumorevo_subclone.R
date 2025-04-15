
#===============================================================================
# mobster getter
#===============================================================================
mobster_getter <- function(dataset, patient_id = NULL, sample_id = NULL) {
  
  p1 <- paste0("reuslts_tumorevo_", tolower(dataset), "/subclonal_deconvolution/mobster/", dataset, "/")
  
  if ( is.null(patient_id) & is.null(sample_id) ) {
    files <- list.files(p1, pattern = "*_mobsterh_st_best_fit.rds", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else if ( !is.null(patient_id) & is.null(sample_id) ) {
    p2 <- paste0(p1, patient_id, "/")
    files <- list.files(p2, pattern = "*_mobsterh_st_best_fit.rds", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else if ( !is.null(patient_id) & !is.null(sample_id) ) {
    p2 <- paste0(p1, patient_id, "/")
    files <- list.files(p2, pattern = paste0("*_", sample_id, "_mobsterh_st_best_fit.rds"), full.names = TRUE, recursive = TRUE) # List all RDS files
  } else {
    stop("something wrong!")
  }
  return(files)
}


#===============================================================================
# pyclonevi getter
#===============================================================================
pyclonevi_getter <- function(dataset, patient_id = NULL) {
  
  p1 <- paste0("reuslts_tumorevo_", tolower(dataset), "/subclonal_deconvolution/pyclonevi/", dataset, "/")
  
  if ( is.null(patient_id) ) {
    files <- list.files(p1, pattern = "*_remove_tail_all_best_fit.txt", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else if ( !is.null(patient_id) ) {
    p2 <- paste0(p1, patient_id, "/")
    files <- list.files(p2, pattern = "*_remove_tail_all_best_fit.txt", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else {
    stop("something wrong!")
  }
  return(files)
}

#===============================================================================
# viber getter
#===============================================================================
viber_getter <- function(dataset, patient_id = NULL) {
  
  p1 <- paste0("reuslts_tumorevo_", tolower(dataset), "/subclonal_deconvolution/viber/", dataset, "/")
  
  if ( is.null(patient_id) ) {
    files <- list.files(p1, pattern = "*_remove_tail_all_viber_best_st_fit.rds", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else if ( !is.null(patient_id) ) {
    p2 <- paste0(p1, patient_id, "/")
    files <- list.files(p2, pattern = "*_remove_tail_all_viber_best_st_fit.rds", full.names = TRUE, recursive = TRUE) # List all RDS files
  } else {
    stop("something wrong!")
  }
  return(files)
}


#===============================================================================
# subclonal deconvolution getter
#===============================================================================
subclone_getter <- function(tool, dataset, patient_id=NULL, sample_id=NULL) {
  
  if (tool=="mobster") {
    output <- mobster_getter(dataset, patient_id, sample_id)
  }
  else if (tool=="pyclonevi") {
    output <- pyclonevi_getter(dataset, patient_id)
  }
  else if (tool=="viber") {
    output <- viber_getter(dataset, patient_id)
  }
  else {
    stop("wrong tool specified!")
  }
  return(output)
}







