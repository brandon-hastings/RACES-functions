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
  accepted_callers <- c("mutect2", "strelka", "ascat")
  if (!(variant_caller %in% accepted_callers)) {
    stop("Variant caller not supported. Or check spelling of variant_caller arg.")
  }
  if (is.null(patientID) & variant_caller == "mutect2") {
    stop("getting mutect2 requires patientID")
  }
  
  # get files
  if (variant_caller %in% c("strelka", "ascat")) {
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


# take a list of filenames from strelka, mutect2, or ascat and structure them 
# into a named list
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
  } else if (grepl("ascat", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "tumour.ASCATprofile.png")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "cn_profile_image" 
      } else if (endsWith(list_of_output_files[i], "tumour.ASPCF.png")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "cn_segmentation_image" 
      } else if (endsWith(list_of_output_files[i], "tumour.sunrise.png")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "sunrise_image" 
      } else if (endsWith(list_of_output_files[i], "cnvs.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "cnvs" 
      } else if (endsWith(list_of_output_files[i], "purityploidy.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "purity_ploidy" 
      } else if (endsWith(list_of_output_files[i], "segments.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "segments" 
      } else if (endsWith(list_of_output_files[i], "tumour_tumourBAF.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "tumour_baf" 
      } else if (endsWith(list_of_output_files[i], "tumour_normalBAF.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "normal_baf"
      } else if (endsWith(list_of_output_files[i], "tumour_tumourLogR.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "tumour_logR" 
      } else if (endsWith(list_of_output_files[i], "tumour_normalLogR.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "normal_logR"
      }
    }
  }
  return(named_files)
}


# Returns a multi-layered list structure, where variant caller is the first layer,
# coverage is the second, purity is the third, and sample ID
# is the last, holding variant call files for each caller/coverage/purity/sample
# combination.
# Coverage, purity, and variant caller can be supplied as parameters in the
# function call. Sample IDs are generated from the file structure.
batch_retrieve_sarek_variant_caller <- function(basedir,
                                                SPN_id_number,
                                                coverage=c(50, 100, 150, 200),
                                                purity=c(0.3, 0.6, 0.9),
                                                variant_callers=c("mutect2",
                                                                  "strelka",
                                                                  "ascat"),
                                                normalID="normal_sample") {
  SPN_ID <- paste("SPN", as.character(SPN_id_number))
  sarek_path <- file.path(as.character(basedir),
                          paste("SPN", as.character(SPN_id_number)),
                          "sarek")
  
  callers <- list(rep(NA, length(variant_callers)))
  for (vc in 1:length(variant_callers)) {
    coverages <- list(rep(NA, length(coverage)))
    for (c in 1:length(coverage))  {
      purities <- list(rep(NA, length(purity)))
      for (p in 1:length(purity)) {
        coverage_purity_combination <- paste(as.character(coverage[c]),
                                             "x_",
                                             as.character(purity[p]),
                                             "p")
        vc_path <- file.path(sarek_path,
                             coverage_purity_combination,
                             "variant_calling",
                             variant_callers[vc])
        sample_dirs <- list.dir(path = vc_path)
        samples <- c(rep(NA, length(sample_dirs)))
        for (s in 1:length(sample_dirs)) {
          if (variant_callers[vc] == "mutect2") {
            mutect_naming <- strsplit(basename(sample_dirs[s]), ",")[[1]]
            sampleID <- mutect_naming[1]
            patientID <- mutect_naming[3]
          } else {
            sampleID <- sample_dirs[s]
          }
          files_list <- get_sarek_variant_called_files(basedir = basedir,
                                                       SPN_id_number = SPN_ID,
                                                       coverage = coverage[c],
                                                       purity = purity[p],
                                                       variant_caller = variant_callers[vc],
                                                       sampleID = sampleID,
                                                       normalID = normalID,
                                                       patientID = patientID)
          parsed_files <- parse_sarek_variant_called_files(files_list)
          patientID <- NULL
          samples[s] <- parsed_files
        }
        names(samples) <- lapply(sample_dirs, basename)
        purities[p] <- samples
      }
      names(purities) <- lapply(purity, as.character)
      coverages[c] <- purities
    }
    names(coverages) <- lapply(coverage, as.character)
    callers[vc] <- coverages
  }
  names(callers) <- variant_callers
  return(callers)
}
