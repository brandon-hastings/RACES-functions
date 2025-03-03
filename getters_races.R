

#rm(list = ls()) # clears objects from the workspace
library(dplyr)
library(rRACES)


#===============================================================================
# AUXILIARY FUNCTIONS
#===============================================================================

# quality control the input parameters
#-------------------------------------------------------------------------------
qc_args <- function(name, coverage, purity, timepoint, sample_id, type) {
  if (!is.numeric(name) || name %% 1 != 0) {
    stop("Error: 'name' must be integer") # name : integer
  }
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!is.null(timepoint) && (!is.numeric(timepoint) || timepoint %% 1 != 0)) {
    stop("Error: 'timepoint' must be integer") # timepoint : integer
  }
  if (!is.null(sample_id) && (!is.numeric(sample_id) || sample_id %% 1 != 0)) {
    stop("Error: 'sample_id' must be integer") # sample_id : integer
  }
  if (!(type %in% c("snv", "cna", "phylo"))) {
    stop("Error: 'type' must be one of 'SNV' and 'CNA'")
  }
  return(TRUE)
}

# generate the pattern for string matching for each required scenario
#-------------------------------------------------------------------------------
gen_ptr <- function(SPN_ID, timepoint, sample_id, type) {
  
  if (type == "cna") {
    timepoint_ph <- ifelse(test = is.null(timepoint), yes = "\\d+", no = timepoint)
    sample_id_ph <- ifelse(test = is.null(sample_id), yes = "\\d+", no = sample_id)
    ptr <- paste0("^SPN0", SPN_ID, "_", timepoint_ph, "\\.", sample_id_ph, "_cna\\.rds$")
    return(ptr)
  }
  
  else if (type == "snv") {
    ptr <- paste0("SPN0", SPN_ID, "/races/purity_", purity, "/seq_results_muts_merged_coverage_", coverage, "x.rds")
    return(ptr)
  }
  
  else if (type == "phylo") {
    ptr <- paste0("SPN0", SPN_ID, "/races/phylo_forest.sff")
    return(ptr)
  }
  
  else {
    stop("Error: invalid!")
  }
}

# return list of absolute file paths as a result of string matching
#-------------------------------------------------------------------------------
search_files <- function(directory, ptr) {
  files <- list.files(path = directory, pattern = ptr, full.names = TRUE, recursive = TRUE)
  return(files)
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



#===============================================================================
# MAIN GETTER FUNCTION
#===============================================================================
# input:
  # PATH      : (string) - abosolute path
  # SPN_ID    : (integer) - SPN ID
  # coverage  : (integer)
  # purity    : (float)
  # timepoint : (integer)
  # sample_id : (integer)
  # type      : ("cna", "snv", "phylo")
# output:
  # list of file object
getter_races <- function(PATH = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", 
                         SPN_ID = 1, 
                         coverage = 100, 
                         purity = 0.6, 
                         timepoint = NULL, 
                         sample_id = NULL, 
                         type = "cna") {
  
  type <- tolwer(type)
  
  tryCatch(
    {
      qc_args(name=SPN_ID, coverage=coverage, purity=purity, timepoint=timepoint, sample_id=sample_id, type=type)
      ptr <- gen_ptr(SPN_ID, timepoint, sample_id, type)
      #message(ptr)
      files_list <- search_files(directory=PATH, ptr)
      
      if (type=="phylo") {
        return(rRACES::load_phylogenetic_forest(files_list))
      }
      
      if ( length(files_list) != 0 ) {
        return(load_multiple_rds(files_list))
      }
      
      # if nothing find
      else {
        message("NOT FOUND!")
        return(NULL)
      }
      
    },
    error = function(cond) {
      message(conditionMessage(cond))
    }, 
    warning = function(cond) {
      message(conditionMessage(cond))
    }, 
    finally = {
      #message("JOB DONE!")
    }
  )
}





