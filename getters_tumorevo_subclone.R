




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
  if (!(tools %in% c("mobster", "pyclonevi", "viber"))) {
    stop("Error: 'type' must be one of 'mobster', 'pyclonevi' and 'viber'")
  }
  return(TRUE)
}


getter_mobster <- function(directory, patient_id = NULL, sample_id = NULL) {
  # List all RDS files in the directory
  files <- list.files(directory, pattern = "*.rds", full.names = TRUE, recursive = TRUE)
  
  # Construct regex pattern dynamically
  patient_pattern <- if (!is.null(patient_id)) paste0("MSeq_", patient_id, "_") else "MSeq_.*_"
  sample_pattern  <- if (!is.null(sample_id)) paste0(sample_id, "_") else ".*_"
  
  # Complete pattern
  pattern <- paste0(patient_pattern, ".*_", sample_pattern, "mobsterh_st_best_fit\\.rds$")
  
  # Filter files matching the pattern
  matching_files <- files[grepl(pattern, files)]
  
  return(matching_files)
}


getter_pyclonevi <- function(directory, patient_id = NULL) {
  # List all RDS files in the directory
  files <- list.files(directory, pattern = "*.txt", full.names = TRUE, recursive = TRUE)
  
  # Construct regex pattern dynamically
  patient_pattern <- if (!is.null(patient_id)) paste0("MSeq_", patient_id, "_") else "MSeq_.*_"
  #sample_pattern  <- if (!is.null(sample_id)) paste0(sample_id, "_") else ".*_"
  
  # Complete pattern
  pattern <- paste0(patient_pattern, ".*_", "remove_tail_all_best_fit\\.txt$")
  
  # Filter files matching the pattern
  matching_files <- files[grepl(pattern, files)]
  
  return(matching_files)
}


getter_viber <- function(directory, patient_id = NULL) {
  # List all RDS files in the directory
  files <- list.files(directory, pattern = "*.rds", full.names = TRUE, recursive = TRUE)
  
  # Construct regex pattern dynamically
  patient_pattern <- if (!is.null(patient_id)) paste0("MSeq_", patient_id, "_") else "MSeq_.*_"
  #sample_pattern  <- if (!is.null(sample_id)) paste0(sample_id, "_") else ".*_"
  
  # Complete pattern
  pattern <- paste0(patient_pattern, "remove_tail_all_viber_best_st_fit\\.rds$")
  
  # Filter files matching the pattern
  matching_files <- files[grepl(pattern, files)]
  
  return(matching_files)
}



# return the file object given the file path
#-------------------------------------------------------------------------------
load_multiple_rds <- function(file_paths) {
  ext <- tolower( tail(strsplit(basename(file_paths), "\\.")[[1]], 1) )
  if (ext != "rds") {
    stop("Invalid file type!")
  }
  data_list <- lapply(file_paths, readRDS)
  names(data_list) <- basename(file_paths)  # Use file names as list names
  return(data_list)
}

load_multiple_txt <- function(file_paths) {
  ext <- tolower( tail(strsplit(basename(file_paths), "\\.")[[1]], 1) )
  if (ext != "txt") {
    stop("Invalid file type!")
  }
  data_list <- lapply(file_paths, read.delim)
  names(data_list) <- basename(file_paths)  # Use file names as list names
  return(data_list)
}



getters_tumorevo_subclone <- function(MAIN_PATH, SPN_ID, coverage, purity, patient_id = NULL, sample_id = NULL, tool) {
  
  path <- paste0(MAIN_PATH, "SPN0", SPN_ID, "/tumorevo/", coverage, "x/", purity, "p/", "results_tumourevo_mseq/subclonal_deconvolution/")
  
  if (tool=="mobster") {
    files_list <- getter_mobster(path, patient_id, sample_id)
    if (length(files_list) != 0) {
      return(load_multiple_rds(files_list))
    }
  }
  if (tool=="pyclonevi") {
    files_list <- getter_pyclonevi(path, patient_id)
    if (length(files_list) != 0) {
      return(load_multiple_txt(files_list))
    }
  }
  if (tool=="viber") {
    files_list <- getter_viber(path, patient_id)
    if (length(files_list) != 0) {
      return(load_multiple_rds(files_list))
    }
  }
  else {
    stop("Error in 'getters_tumorevo_subclone' function!")
  }
  
}





