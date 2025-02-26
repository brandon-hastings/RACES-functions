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

sarek_output <- batch_retrieve_sarek_variant_caller()

# take in list of parameters? IDK yet. Probably bad idea
format_tumourevo_samplesheet <- function(params_list) {
  tumourevo_samplesheet <- NULL
  write.csv(tumourevo_samplesheet)
  return(tumourevo_samplesheet)
}
