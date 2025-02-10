# master script for analysis of rRACES runs
# if slurm is needed, don't know how complex the other functions are
library(rslurm)

# results 
combine_results <- function(results_lists){
  # combine results function in progress
  warning("combine_results function was run but is not complete")
  
  for (i in results_lists){
    # pass
  }
}
# assuming that each function has different args, input each as list
launch_scripts <- function(cna_args, driver_args, sbatch_options = NULL){
  # basic argument checking
  stopifnot(is.character(results_dir))
  stopifnot(is.list(cna_args))
  stopifnot(is.list(driver_args))
  
  if (is.null(sbatch_options)){
    # generic slurm settings for each job
    sopt <- list(time = '1:00:00', account = 'cdslab', partition = 'THIN',
                 mem = '8G')
  } else if (is.list(sbatch_options)){
    # check that correct names are present in supplied list
    if ("time" %in% names(sbatch_options) && "account" %in% names(sbatch_options) && "partition" %in% names(sbatch_options) && "mem" %in% names(sbatch_options)){
      sopt <- sbatch_options
    } else {
      warning("Names supplied in sbatch_options list argument are incorrect or missing")
    }
  } else {
    warning("sbatch_options argument was supplied, but must be a list")
  }

  # QUESTION are the outputs of these functions the same structure?
  # example of launching each job as an individual slurm job
  cna_job <- slurm_call(CNA_function, params = cna_args,
                        jobname = 'CNA_val', submit = TRUE, slurm_options = sopt)
  drivers_job <- slurm_call(driver_function, params = driver_args,
                               jobname = 'driver_val', submit = TRUE, slurm_options = sopt)
  
  # collect results for all jobs in list of lists and clean up temp files
  submitted_jobs <- c(cna_job, drivers_job)
  collected_results <- list()
  # avoid bottleneck by constantly evaluating this and not waiting for slurm_out
  while (length(submitted_jobs) != length(collected_results)) {
    # check every 5 minutes until all jobs are completed
    Sys.sleep(60*5)
    for (i in submitted_jobs){
      # QUESTION what happens to get_slurm_out when a job fails?
      collected_results[[length(collected_results)+1]] <- get_slurm_out(i, wait = FALSE)
      cleanup_files(i)
    }
  }
  
  # combine results from each job into a uniform structure
  combine_results(collected_results)
}

# create lists of parameters for each tool and run main function
cna_args <- list()
driver_args <- list()
launch_scripts(cna_args = cna_args, driver_args = driver_args)
