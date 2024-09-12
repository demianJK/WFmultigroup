# this function simulates data in LF separately for the between- and within-cluster level
# where the latter is discerned in two groups (k=2)

generateData <- function(var_B, VR, n, g, p){ 

  ## variances
  
  var_W_1 <- 1-var_B
  var_W_2 <- var_W_1/VR
  
  var_w_1_All <- c()
  var_w_2_All <- c()
  var_b_All <- c()
  for (i in 1:p) {
    var_w_1_All[i] <- paste0("x", i, "~~", var_W_1, "*", "x", i)
    var_w_2_All[i] <- paste0("x", i, "~~", var_W_2, "*", "x", i)
    var_b_All[i] <- paste0("x", i, "~~", var_B, "*", "x", i)
  }
  var_w_1_All <- paste(var_w_1_All, collapse = "; ")
  var_w_2_All <- paste(var_w_2_All, collapse = "; ")
  var_b_All <- paste(var_b_All, collapse = "; ")
  
  ## covariances 
  
  # cov depends on corr (fix) and var (varying)
  cor_B <- 0.3
  cor_W <- 0.3
  # cov_xy  = corr_xy * var_x/y (because sd/var same for every variable, e.g., x and y)
  cov_w_1 <- cor_W * var_W_1
  cov_w_2 <- cor_W * var_W_2
  cov_b <- cor_B * var_B
  covs_w_1 <- c()
  covs_w_2 <- c()
  covs_b <- c()
  count = 0
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j & j > i) {
        count <- count + 1
        covs_w_1[count] <- paste("x", i, "~~", cov_w_1, "*", "x", j, sep = "")
        covs_w_2[count] <- paste("x", i, "~~", cov_w_2, "*", "x", j, sep = "")
        covs_b[count] <- paste("x", i, "~~", cov_b, "*", "x", j, sep = "")
      }
    }
  }
  covs_w_1 <- paste(covs_w_1, collapse = "; ")
  covs_w_2 <- paste(covs_w_2, collapse = "; ")
  covs_b <- paste(covs_b, collapse = "; ")
  
  # put variances and covariances together (means are 0 per default)
  popModel_W_1 <- paste(var_w_1_All, covs_w_1, sep = ";")
  popModel_W_2 <- paste(var_w_2_All, covs_w_2, sep = ";")
  popModel_B <- paste(var_b_All, covs_b, sep = ";")
  
  # generate data for each level
  sample_B <- simulateData(popModel_B, sample.nobs = g, model.type = "lavaan")
  sample_W_1 <- simulateData(popModel_W_1, sample.nobs = n*g/2, model.type = "lavaan")
  sample_W_2 <- simulateData(popModel_W_2, sample.nobs = n*g/2, model.type = "lavaan")
  
  # merge data 
  count <- 0 # unique index  
  clusters <- rep(1:g, each = n)
  data <- as.data.frame(matrix(NA, ncol=p, nrow=g*n))
  for (k in 1:2) { # merge the sampled data from both levels
    sample_W <- get(paste0("sample_W_", k)) 
    for (i in 1:(n*g/2)){
      count <- count + 1
      j <- clusters[count]
      data[count,] <- sample_W[i,] + sample_B[j,]
    }
  }
  colnames(data) <- paste0("x", 1:p) # unfortunately deleted in loop before
  # coding variables
  data$j <- as.factor(clusters) # clusters 
  data$i <- rep(1:n, g) # unit in group
  data$k <- rep(c(1, 2), each=n*g/2) # grouping variable
  
  return(data)
}