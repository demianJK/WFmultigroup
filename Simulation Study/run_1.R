# multiple scripts in individual sessions are run for using multiple cores
# copy "run_1.R" times the number of cores you want to use and change the number in the name and in the respective script (it's used for the random seed below)
script <- 1

# set paths
pathS <- "scripts/functions/" # path to source local functions 
pathF <- "data/fit/" # path to save fit 

# source local functions
source(paste0(pathS, "generateData.R")) # simulate data
source(paste0(pathS, "analyseData.R")) # analyse with WFmultigroup in lavaan
source(paste0(pathS, "projectInfo.R"))

# (install and) load packages
packages <- c("lavaan", # (version 0.6-17)
              "tidyr", # (version 1.3.1)
              "dplyr") # (version 1.1.4)
newPackages <-
  packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# get simulation design
pI <- projectInfo()
SimConds <- readRDS("objects/SimConds.rds")
nconds <- nrow(SimConds)

# use multiple cores through initializing multiple R sessions
ncores <- 20
nrep <- pI$nrep
scriptReps <- seq(script, nrep, ncores)
set.seed(as.numeric(paste0(script, nrep, ncores, collapse="")))

for (rep in scriptReps){
  for (cond in 1:nconds){
    
    n <- SimConds$n[cond]
    g <- SimConds$g[cond]
    p <- SimConds$p[cond]
    var_B <- SimConds$var_B[cond]
    VR <- SimConds$VR[cond]
    
    ## generate data
    data_LF <- generateData(var_B, VR, n, g, p)

    ## reformat to WF
    data_WF <- pivot_wider(data_LF, names_from = "i", values_from = paste0("x", 1:p), names_sep = ".")
    
    ## analyse data
    fit <- analyseData(data_WF, VR, n, g, p)
    saveRDS(fit, file = paste0(pathF, "fit_C", cond, "_R", rep, ".rds"))
    
  }
  gc()
  print(rep)
}

# for reproducibility
if (script == 1){
  # hardware and software details
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  # create an .R file with all of the function definitions in the current search space
  dump(lsf.str(), file="sourcedFunctions.R")
}
