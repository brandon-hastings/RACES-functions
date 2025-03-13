library(rlist)
# return list of files from sarek for a given caller/purity/coverage combination
get_sarek_variant_called_files <- function(basedir,
                                           coverage,
                                           purity,
                                           variant_caller,
                                           sampleID,
                                           normalID="normal_sample",
                                           patientID=NULL) {
  # testing
  # input checking
  accepted_callers <- c("mutect2", "strelka", "ascat", "freebayes")
  if (!(variant_caller %in% accepted_callers)) {
    stop("Variant caller not supported. Or check spelling of variant_caller arg.")
  }
  
  # get files
  if (variant_caller %in% c("strelka", "ascat", "freebayes")) {
    sample_naming <- paste0(sampleID, "_vs_", normalID)
    path_to_files <- file.path(basedir, "sarek",
                               paste0(as.character(coverage),
                                     "x_",
                                     as.character(purity),
                                     "p"),
                               "variant_calling",
                               variant_caller,
                               sample_naming)
    output_files <- list.files(path_to_files, full.names = TRUE)
  } else if (variant_caller == "mutect2") {
    # sample_naming <- paste(sampleID, paste(sampleID,
    #                                        "_vs_",
    #                                        normalID),
    #                        patientID, sep = ",")
    sample_naming <- sampleID
    path_to_files <- file.path(basedir, "sarek",
                                paste0(as.character(coverage),
                                      "x_",
                                      as.character(purity),
                                      "p"),
                                "variant_calling",
                                variant_caller,
                                sample_naming
                                )
    output_files <- list.files(path_to_files, full.names = TRUE)
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
      if (endsWith(list_of_output_files[i], "filtered.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "vcf" 
      } else if (endsWith(list_of_output_files[i], "filtered.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "tbi" 
      }
    }
  } else if (grepl("strelka", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "snvs.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "snvs_vcf" 
      } else if (endsWith(list_of_output_files[i], "snvs.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "snvs_tbi" 
      }
    }
  } else if (grepl("freebayes", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "freebayes.vcf.gz")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "vcf" 
      } else if (endsWith(list_of_output_files[i], "freebayes.vcf.gz.tbi")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "tbi" 
      }
    }
  } else if (grepl("ascat", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- list(rep(NA, length(list_of_output_files)))
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "purityploidy.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "purity_ploidy" 
      } else if (endsWith(list_of_output_files[i], "segments.txt")) {
        named_files[i] <- list_of_output_files[i]
        names(named_files)[i] <- "segments"
      }
    }
  }
  named_files <- named_files[!is.na(names(named_files))]
  return(named_files)
}

# identify normal sample in strelka and drop it
drop_normal <- function(samples_list,
                        normal="normal_sample",
                        return_normal=FALSE) {
  index <- NULL
  for (i in 1:length(samples_list)) {
    if (basename(samples_list[i]) == normal) {
      index <- i
    } else {
      next
    }
  }
  if (is.null(index)) {
    return(samples_list)
  } else if (return_normal == TRUE) {
    return(index)
  } else {
    return(samples_list[-index])
  }
}
# Returns a multi-layered list structure, where coverage/purity is the first layer,
# caller is the second, and sample ID is the last, holding variant call files 
# for each caller/coverage/purity/sample combination.
# Coverage, purity, and variant caller can be supplied as parameters in the
# function call. Sample IDs are generated from the file structure.
batch_retrieve_sarek_variant_caller <- function(basedir,
                                                SPN_ID,
                                                coverage=c(50, 100, 150, 200),
                                                purity=c(0.3, 0.6, 0.9),
                                                variant_callers=c("mutect2",
                                                                  "strelka",
                                                                  "ascat"),
                                                normalID="normal_sample",
                                                patientID=NULL) {
  
  coverage_purity_combo <- list()
  for (c in 1:length(coverage))  {
    purities <- list(rep(NA, length(purity)))
    for (p in 1:length(purity)) {
      coverage_purity_combination <- paste0(as.character(coverage[c]),
                                           "x_",
                                           as.character(purity[p]),
                                           "p")
      callers <- list()
        for (vc in 1:length(variant_callers)) {
          vc_path <- file.path(basedir,
                               "sarek",
                               coverage_purity_combination,
                               "variant_calling",
                               variant_callers[vc])
          sample_dirs <- list.dirs(path = vc_path)[-1]
          # TODO: solve for mutect but can't just move it into if block
          sample_dirs <- drop_normal(sample_dirs)
          samples <- list()
          if (variant_callers[vc] != "mutect2") {
            for (s in 1:length(sample_dirs)) {
              sampleID <- strsplit(basename(sample_dirs[s]), "_vs_")[[1]][1]
              files_list <- get_sarek_variant_called_files(basedir = basedir,
                                                           coverage = coverage[c],
                                                           purity = purity[p],
                                                           variant_caller = variant_callers[vc],
                                                           sampleID = sampleID,
                                                           normalID = normalID,
                                                           patientID = patientID)
              parsed_files <- parse_sarek_variant_called_files(files_list)
              patientID <- NULL
              samples <- append(samples, list(parsed_files))
            }
            names(samples) <- lapply(sample_dirs, basename)
          }
          else if (variant_callers[vc] == "mutect2") {
            # same as SPN ID
            joint_called_dir_name <- basename(file.path(basedir))
            # get the joint called files once
            files_list <- get_sarek_variant_called_files(basedir = basedir,
                                                                   coverage = coverage[c],
                                                                   purity = purity[p],
                                                                   variant_caller = variant_callers[vc],
                                                                   sampleID = joint_called_dir_name,
                                                                   normalID = normalID,
                                                                   patientID = patientID)
            # drop the joint called folder from samples
            sample_dirs <- drop_normal(sample_dirs, normal=joint_called_dir_name)
            for (s in 1:length(sample_dirs)) {
              parsed_files <- parse_sarek_variant_called_files(files_list)
              samples <- append(samples, list(parsed_files))
              sampleID <- basename(sample_dirs[s])
              # SPN01_SPN01_1.3
            }
            print(lapply(sample_dirs, basename))
            names(samples) <- lapply(sample_dirs, basename)
            # names(samples) <- lapply(lapply(sample_dirs, basename), function(x)paste0(joint_called_dir_name, "_", x))
          }
        # names(samples) <- lapply(sample_dirs, basename)
        callers <- append(callers, list(samples))
        }
      names(callers) <- variant_callers
      coverage_purity_combo <- append(coverage_purity_combo, list(callers))
      names(coverage_purity_combo)[c*p] <- coverage_purity_combination
    }
  }
  return(coverage_purity_combo)
}
