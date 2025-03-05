library(data.table)
library(tidyverse)
source("getters_sarek.R")

# parse args. Make more specific for optional args (coverage, purity, etc.)
parse_args <- function() {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    if (dir.exists(args[1])) {
      basedir <- args[1]
      SPN_number <- args[2]
    } else if (dir.exists(args[2])) {
      basedir <- args[2]
      SPN_number <- args[1]
    } else {
      stop("could not parse arguments or basedir does not exist")
    }
    else {
      stop("Requires 2 arguments, the basedir and the SPN number")
    }
  }
  return(list(basedir, SPN_number))
}

arguments <- parse_args
sarek_output <- batch_retrieve_sarek_variant_caller(basedir = arguments[1],
                                                    SPN_id_number = arguments[2])
# 
check_vector_convert <- function(argument) {
  if (!is.vector(argument)) {
    ret_val <- c(argument)
  } else {
    ret_val <- argument
  }
  return(ret_val)
}


tumoruevo_patient_sample_naming <- function(samplesheet) {
  # for ascat and strelka
  # take sample dir naming from row names
  # split on _vs_
  # first one will be SPN01_{PATIENT}.{SAMPLE}
}

format_samplesheet <- function(samplesheet, cna_caller, cancer_type, normalID) {
  renamed_samplesheet <- samplesheet %>% 
    rename_at(vars(contains("vcf")), ~ "vcf") %>% 
    rename_at(vars(contains("tbi")), ~ "tbi") %>% 
    rename_at(vars(contains("segments")), ~ "cna_segments") %>% 
    rename_at(vars(contains("purity_ploidy")), ~ "cna_extra")
  renamed_samplesheet$caller <- c(rep(cna_caller, nrow(renamed_samplesheet)))
  renamed_samplesheet$cancer_type <- c(rep(cancer_type, nrow(renamed_samplesheet)))
  renamed_samplesheet$normal_sample <- c(rep(normalID, nrow(renamed_samplesheet)))
  return(renamed_samplesheet)
  
}


# dataframe where caller is the column name, sample is row name,
# file is cell value
callers_to_data_frame <- function(caller_list) {
  callers_datatables <- lapply(caller_list, rbindlist(fill = TRUE))
  # make sure rownames match
  raw_samplesheet <- merge(callers_datatables[1], callers_datatables[2],
        by = 'row.names',
        all = TRUE)
  return(raw_samplesheet)
}

filter_callers <- function(samplesheet) {
  ascat_files <- c("segments", "purity_ploidy")
  mutect_files <- c("vcf", "tbi")
  strelka_files <- c("snvs_vcf", "snvs_tbi")
  # combine above lists here. Easier to maintain if they are separate
  search_columns <- c(strelka_files, mutect_files, ascat_files)
  filtered_samplesheet <- samplesheet[ ,colnames(samplesheet) %in% search_columns, with=FALSE]
  return(filtered_samplesheet)
}


# can do single or batch processing
create_tumourevo_samplesheet <- function(basedir,
                                         SPN_id_number,
                                         coverage,
                                         purity,
                                         snp_caller,
                                         cna_caller="ascat",
                                         cancer_type="PANCANCER",
                                         normalID="normal_sample") {
  # hacky input conversion for batch processing single coverage/purity combo
  coverage <- check_vector_convert(coverage)
  purity <- check_vector_convert(purity)
  caller <- c(snp_caller, cna_caller)
  
  getter <- batch_retrieve_sarek_variant_caller(basedir = basedir,
                                                SPN_id_number = SPN_id_number,
                                                coverage = coverage,
                                                purity = purity,
                                                variant_callers = caller,
                                                normalID = normalID)
  # process from batch structure, even if only one coverage/purity is given
  # make a dataframe for each coverage/purity combination
  for (c_p in getter) {
    # combine snp and cna caller on samples, return all files as columns
    raw_samplesheet <- callers_to_data_frame(c_p)
    # filter to needed columns
    caller_filtered_samplesheet <- filter_callers(raw_samplesheet)
    # construct patient/sample naming
    named_samplesheet <- tumoruevo_patient_sample_naming(caller_filtered_samplesheet)
    # final formatting/column ordering
    tumourevo_samplesheet <- format_samplesheet(named_samplesheet,
                                                cna_caller,
                                                cancer_type,
                                                normalID)
    # save samplesheet as csv
    filename <- "some file path name"
    write.csv(tumourevo_samplesheet, file = filename)
    # check exist and return message
    if (file.exists(filename)) {
      print(paste0("samplesheet created at ", filename))
    }
  }
}
