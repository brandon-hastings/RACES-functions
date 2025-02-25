
# Returns a multi-layered list structure, where coverage is the first layer,
# purity is the second, name of the variant caller is the third, and sample ID
# is the last, holding variant call files for each sample/caller/purity/coverage
# combination.
# Coverage, purity, and variant caller can be supplied as parameters in the
# function call. Sample IDs are generated from the file structure.

# ---Work in progress---
parse_sarek_variant_caller_dir <- function(basedir,
                                       SPN_id_number,
                                       coverage=c(50, 100, 150, 200),
                                       purity=c(0.3, 0.6, 0.9),
                                       variant_callers=c("mutect2", "strelka")) {
  sarek_path <- file.path(as.character(basedir),
                              paste("SPN", as.character(SPN_id_number)),
                              "sarek")
  # can be improved here. apply function? access through not a list?
  coverages <- list(rep(NA, length(coverage)))
  for (c in 1:length(coverage))  {
    purities <- list(rep(NA, length(purity)))
    for (p in 1:length(purity)) {
      coverage_purity_combination <- paste(as.character(coverage[c]),
                                           "x_",
                                           as.character(purity[p]),
                                           "p")
      callers <- list(rep(NA, length(variant_callers)))
      for (vc in 1:length(variant_callers)) {
        vc_path <- file.path(sarek_path,
                             coverage_purity_combination,
                             "variant_calling",
                             variant_callers[vc])
        sample_dirs <- list.dir(path = vc_path)
        samples <- c(rep(NA, length(sample_dirs)))
        for (d in 1:length(sample_dirs)) {
          files_list <- list.files(sample_dirs[d])
          samples[d] <- files_list
        }
        names(samples) <- lapply(sample_dirs, basename)
        callers[vc] <- samples
      }
      names(callers) <- variant_callers
      purities[p] <- callers
    }
    names(purities) <- as.character(purity)
    coverages[c] <- purities
  }
  names(coverages) <- as.character(coverage)
  
  return(coverages)
}

# return list of files from sarek for a given caller/purity/coverage combination
get_sarek_variant_called_files <- function(basedir,
                                           SPN_id_number,
                                           coverage,
                                           purity,
                                           variant_caller,
                                           sampleID,
                                           normalID="normal_sample",
                                           patientID=NULL) {
  # input checking
  accepted_callers <- c("mutect2", "strelka")
  if (!(variant_caller %in% accepted_callers)) {
    stop("Variant caller not supported. Or check spelling of variant_caller arg.")
  }
  if (is.null(patientID) & variant_caller == "mutect2") {
    stop("getting mutect2 requires patientID")
  }
  
  # get files
  if (variant_caller == "strelka") {
    sample_naming <- paste(sampleID, "_vs_", normalID)
    output_files <- list.files(file.path(basedir, SPN_id_number, "sarek",
                                  paste(as.character(coverage),
                                        "x_",
                                        as.character(purity),
                                        "p"),
                                  "variant_calling",
                                  variant_caller,
                                  sample_naming))
  } else if (variant_caller == "mutect2") {
    sample_naming <- paste(sampleID, paste(sampleID,
                                           "_vs_",
                                           normalID),
                           patientID, sep = ",")
    output_files <- list.files(file.path(basedir, SPN_id_number, "sarek",
                                  paste(as.character(coverage),
                                        "x_",
                                        as.character(purity),
                                        "p"),
                                  "variant_calling",
                                  variant_caller,
                                  sample_naming
                                  ))
  }
  return(output_files)
}

# take a list of filenames from strelka or mutect and structure them into
# a named list
parse_sarek_variant_called_files <- function(list_of_output_files) {
  # check which caller we are looking at
  if (grepl("mutect", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "vcf" 
      } else if (endsWith(list_of_output_files[i], "vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "tbi" 
      } else if (endsWith(list_of_output_files[i], "vcf.gz.stats")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "stats" 
      } else if (endsWith(list_of_output_files[i], "contamination.table")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "contamination_table" 
      } else if (endsWith(list_of_output_files[i], "segmentation.table")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "segmentation_table" 
      } else if (endsWith(list_of_output_files[i], "artifactprior.tar.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "artifactprior" 
      } else if (endsWith(list_of_output_files[i], "pileups.table")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "pileups_table" 
      } else if (endsWith(list_of_output_files[i], "filtered.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "filtered_vcf"
      } else if (endsWith(list_of_output_files[i], "filtered.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "filtered_tbi" 
      } else if (endsWith(list_of_output_files[i], "filteringStats.tsv")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "filtered_stats"
      }
    }
  } else if (grepl("strelka", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "indels.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "indels_vcf" 
      } else if (endsWith(list_of_output_files[i], "indels.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "indels_tbi" 
      } else if (endsWith(list_of_output_files[i], "snvs.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "snvs_vcf" 
      } else if (endsWith(list_of_output_files[i], "snvs.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "snvs_tbi" 
      }
    }
  }
  return(named_files)
}