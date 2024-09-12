projectInfo <- function(){
  
  out <- list(
  approaches = c("WFmultigroup"),
  nrep = 1000,
  # population
  VR = c(2, 5),
  var_B = c(0.05, 0.25),
  # sample
  p = c(2),
  n = c(2, 10, 30), 
  g = c(200, 500, 1000)
  )

  if ( !file.exists("objects/SimConds.rds") ){ # create Sim Table if it does not exist
    final <- expand.grid(out[3:length(out)])
    final$var_W_1 <- 1 - final$var_B # within-cluster variance in group 1
    final$var_W_2 <- final$var_W_1/final$VR # within-cluster variance in group 2
    final <- final[order(final$VR),] # sort with increasing p
    final$cond <-  1:nrow(final)
    final <- final[, c("cond", "VR", "var_B", "var_W_1", "var_W_2", "p", "n", "g")]
    saveRDS(final, file = "objects/SimConds.rds" ) 
  }
  
  return(out)
  
}


