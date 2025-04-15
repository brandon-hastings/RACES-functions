
source("getters_races.R")
source("getter_tumorevo_signatures.R")
source("getters_tumorevo_subclone.R")

wrapper <- function(MAIN_PATH = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", 
                    SPN_ID = 1, 
                    coverage = 100, 
                    purity = 0.6, 
                    timepoint = NULL, 
                    sample_id = NULL, 
                    type = "cna") {
  
  message("hello!")
  
}





#MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type
list(races=c(MAIN_PATH, 
              SPN_ID, 
              coverage, 
              purity, 
              timepoint, 
              sample_id, 
              type), 
     tumorevo=0, 
     sarek=0)


# getters_races.R
races <- getter_races(MAIN_PATH, 
                      SPN_ID, 
                      coverage, 
                      purity, 
                      timepoint, 
                      sample_id, 
                      type)

# getter_tumorevo_signatures.R
signature <- mutational_signature_getter(MAIN_PATH, 
                                         SPN_ID, 
                                         coverage, 
                                         purity, 
                                         tools, 
                                         dataset, 
                                         context)

# getters_tumorevo_subclone.R
sunclone <- getters_tumorevo_subclone(MAIN_PATH, 
                                      SPN_ID, 
                                      coverage, 
                                      purity, 
                                      patient_id = NULL, 
                                      sample_id = NULL, 
                                      tool)
  
  
  
  
  
  
  
  
  
  
  