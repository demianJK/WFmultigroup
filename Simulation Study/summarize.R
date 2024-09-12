### summarize results of Monte Carlo study

# (install and) load packages
packages <- c("lavaan", # (version 0.6-17)
              "tidyr", # (version 1.3.1)
              "dplyr", # (version 1.1.4)
              "DescTools" # (version 0.99.50)
)
newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# set paths and load data
pathO <- "objects/" # path for aggregated infos
pathF <- "data/fit/" 
source("scripts/functions/projectInfo.R") # infos about MC study
pI <- projectInfo()

SimConds <- readRDS(paste0(pathO, "SimConds.rds"))
ncond <- nrow(SimConds)
nrep <- pI$nrep
approaches <- pI$approaches
napproach <- length(approaches)

# fixed simulation factors
p <- 2
k <- 2
pc <- p * (p + 1) / 2 # number of unique covariance matrix elements (in LF)
c <- pc-p # numb of unique covariances


## Data -> ListStats ######################################################################################################################################################
# put all relevant info from single rds files lavaan objects into one list

critsList <- c("conv", # Did the model converge?
               "time", # How long did it take till the model converged?
               "W_1", # within-cluster model parameters
               "W_2", 
               "B", # between-cluster model parameters
               "se_W_1", # is there at least one standard error missing?
               "se_W_2",
               "se_B",
               "coverage_W_1", # TRUE/FALSE whether confidence interval encompasses population parameter
               "coverage_W_2",
               "coverage_B",
               "errorMessages"
) 

# create list to save results of each replication
ListStats <- setNames(as.list(approaches), approaches)
ListStats <- lapply(ListStats, function(x){
  setNames(vector("list", length(critsList)), critsList)
})
for (i in 1:napproach){ # initializing necessary
  for (j in 1:length(critsList)){
    ListStats[[i]][[j]] <- setNames(vector("list", ncond), 1:ncond)
    for (k in 1:ncond){
      ListStats[[i]][[j]][[k]] <- setNames(vector("list", nrep), 1:nrep)
    }
  }
}

for (rep in 1:nrep){    
  for (cond in 1:ncond){
    
    ## sample characteristics
    n <- SimConds$n[cond]
    g <- SimConds$g[cond]
    
    ## population parameters
    # all
    VR <- SimConds$VR[cond]
    popVar_B <- SimConds$var_B[cond]
    popVar_W_1 <- SimConds$var_W_1[cond]
    cor <- 0.3
    popCov_B <- cor * popVar_B
    popCov_W_1 <- cor * popVar_W_1
    popParams_W_1 <- c(rep(popVar_W_1, p), rep(popCov_W_1, c)) 
    popParams_B <- c(rep(popVar_B, p), rep(popCov_B, c)) 
    popVar_W_2 <- SimConds$var_W_2[cond]
    popCov_W_2 <- cor * popVar_W_2
    popParams_W_2 <- c(rep(popVar_W_2, p), rep(popCov_W_2, c)) 
    popParams_B <- c(rep(popVar_B, p), rep(popCov_B, c))
    
    # get model fits
    fitList <- readRDS(paste0(pathF, "fit_C", cond, "_R", rep, ".rds"))

    for (i in 1:napproach){
      fit <- fitList[[i]]
      
      ListStats[[approaches[i]]][["conv"]][[cond]][[rep]] <- FALSE
      ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["W_1"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["W_2"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["se_W_1"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["se_W_2"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["se_B"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["coverage_W_1"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["coverage_W_2"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["coverage_B"]][[cond]][[rep]] <- NA
      
      if(!is.atomic(fit) && fit@Fit@converged){ # if exists and optimizer says converged
        ListStats[[approaches[i]]][["conv"]][[cond]][[rep]] <- fit@Fit@converged
        ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- unname(fit@timing[["total"]]) # in s
          
        ## parameter estimates
          
          # between
          allB <- which(grepl("p", parameterEstimates(fit)$label, fixed = TRUE))
          id_x1_b <- allB[1]# which(parameterEstimates(fit)$label == ".p15.")[1]
          id_x2_b <- allB[2]# which(parameterEstimates(fit)$label == ".p16.")[1]
          id_x12_b <- allB[3]# which(parameterEstimates(fit)$label == ".p17.")[1]
          id_b <- c(id_x1_b, id_x2_b, id_x12_b)
          
          # within (group 1)
          id_x1_w_g1 <- which(parameterEstimates(fit)$label == "Vx1_w_g1")[1]
          id_x2_w_g1 <- which(parameterEstimates(fit)$label == "Vx2_w_g1")[1]
          id_x12_w_g1 <- which(parameterEstimates(fit)$label == "Cx12_w_g1")[1]
          id_w_g1 <- c(id_x1_w_g1, id_x2_w_g1, id_x12_w_g1)
          
          # within (group 2)
          id_x1_w_g2 <- which(parameterEstimates(fit)$label == "Vx1_w_g2")[1]
          id_x2_w_g2 <- which(parameterEstimates(fit)$label == "Vx2_w_g2")[1]
          id_x12_w_g2 <- which(parameterEstimates(fit)$label == "Cx12_w_g2")[1]
          id_w_g2 <- c(id_x1_w_g2, id_x2_w_g2, id_x12_w_g2)
          
          ListStats[[approaches[i]]][["W_1"]][[cond]][[rep]] <- parameterEstimates(fit)[id_w_g1, "est"]
          ListStats[[approaches[i]]][["W_2"]][[cond]][[rep]] <- parameterEstimates(fit)[id_w_g2, "est"]
          ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- parameterEstimates(fit)[id_b, "est"]
          
          
          ## se and coverage
          
          se_W_1 <- !anyNA(as.matrix(parameterEstimates(fit)[id_w_g1, "se"]))
          ListStats[[approaches[i]]][["se_W_1"]][[cond]][[rep]] <- se_W_1
          if (se_W_1){
            ListStats[[approaches[i]]][["coverage_W_1"]][[cond]][[rep]] <- popParams_W_1 %()% as.matrix(parameterEstimates(fit)[id_w_g1, c("ci.lower", "ci.upper")])
          }
          se_W_2 <- !anyNA(as.matrix(parameterEstimates(fit)[id_w_g2, "se"]))
          ListStats[[approaches[i]]][["se_W_2"]][[cond]][[rep]] <- se_W_2
          if (se_W_2){
            ListStats[[approaches[i]]][["coverage_W_2"]][[cond]][[rep]] <- popParams_W_2 %()% as.matrix(parameterEstimates(fit)[id_w_g2, c("ci.lower", "ci.upper")])
          }
          
          se_B <- !anyNA(as.matrix(parameterEstimates(fit)[id_b, "se"]))
          ListStats[[approaches[i]]][["se_B"]][[cond]][[rep]] <- se_B
          if (se_B){
            ListStats[[approaches[i]]][["coverage_B"]][[cond]][[rep]] <- popParams_B %()% as.matrix(parameterEstimates(fit)[id_b, c("ci.lower", "ci.upper")])
          }
          
          
      } else if (!is.null(attr(fit, "condition"))) { # save error note
        ListStats[[approaches[i]]][["errorMessages"]][[cond]][[rep]] <- attr(fit, "condition")[[1]] # only error, not function call
      }
      
    }

    gc()
  }
  print(rep)
}

saveRDS(ListStats, file = paste0(pathO, "ListStats.rds"))



### ListStats --> TableStats  ######################################################################################################################################################
# summarize info from list into table

estim <- c("RMSE", "Bias")
critsTable <- c("conv", 
                "time", "timesd",
                
                ## within
                # variance (heterogenous)
                paste(estim, "var", "W_1", sep="_"),
                paste(paste0("rel", estim), "var", "W_1", sep="_"),
                paste(estim, "var", "W_2", sep="_"),
                paste(paste0("rel", estim), "var", "W_2", sep="_"),
                # variance (total)
                paste(estim, "var", "W", sep="_"),
                paste(paste0("rel", estim), "var", "W", sep="_"),
                # covariance (heterogenous)
                paste(estim, "cov", "W_1", sep="_"),
                paste(paste0("rel", estim), "cov", "W_1", sep="_"),
                paste(estim, "cov", "W_2", sep="_"),
                paste(paste0("rel", estim), "cov", "W_2", sep="_"),
                # covariance (total)
                paste(estim, "cov", "W", sep="_"),
                paste(paste0("rel", estim), "cov", "W", sep="_"),
                # overall (heterogenous)
                paste(estim, "W_1", sep="_"),
                paste(paste0("rel", estim), "W_1", sep="_"),
                paste(estim, "W_2", sep="_"),
                paste(paste0("rel", estim), "W_2", sep="_"),
                # overall (total)
                paste(estim, "W", sep="_"),
                paste(paste0("rel", estim), "W", sep="_"),
                
                # negative variances
                "negVar_W", "negVar_W_1", "negVar_W_2",
                
                ## between
                # variance
                paste(estim, "var", "B", sep="_"),
                paste(paste0("rel", estim), "var", "B", sep="_"),
                # covariance
                paste(estim, "cov", "B", sep="_"),
                paste(paste0("rel", estim), "cov", "B", sep="_"),
                # overall
                paste(estim, "B", sep="_"),
                paste(paste0("rel", estim), "B", sep="_"),
                
                # negative variances
                "negVar_B",
                
                # se and coverage
                "se_W",
                "se_W_1",
                "se_W_2",
                "se_B",
                "coverage_W",
                "coverage_W_1",
                "coverage_W_2",
                "coverage_B",
                
                # Intraclass Correlation
                paste(estim, "ICC", sep="_"),
                paste(paste0("rel", estim), "ICC", sep="_"),
                paste(estim, "ICC_1", sep="_"),
                paste(paste0("rel", estim), "ICC_1", sep="_"),
                paste(estim, "ICC_2", sep="_"),
                paste(paste0("rel", estim), "ICC_2", sep="_")
)

colTable <- c()
for (i in 1:napproach){
  colTable <- append(colTable, paste(critsTable, approaches[i], sep="_"))
}
tmp <- as.data.frame(matrix(NA, ncol=length(colTable), nrow=ncond))
colnames(tmp) <- colTable
TableStats <- SimConds
TableStats <- cbind(TableStats, tmp)

# note that the relative parameters (e.g., relRMSE) are not computed in the order according to the formula, 
# but the order is changed on mathematical basis to allow for more efficient code

for (cond in 1:ncond){
  
  n <- SimConds$n[cond]
  g <- SimConds$g[cond]
  
  ## population parameters
  # all
  VR <- SimConds$VR[cond]
  popVar_B <- SimConds$var_B[cond]
  popVar_W_1 <- SimConds$var_W_1[cond]
  cor <- 0.3
  popCov_B <- cor * popVar_B
  popCov_W_1 <- cor * popVar_W_1
  popICC_1 <- popVar_B/(popVar_B+popVar_W_1)
  popVar_W_2 <- SimConds$var_W_2[cond]
  popVar_W_mean <- mean(c(popVar_W_1, popVar_W_2)) 
  popCov_W_2 <- cor * popVar_W_2
  popCov_W_mean <- mean(c(popCov_W_1, popCov_W_2))
  popICC_2 <-  popVar_B/(popVar_B+popVar_W_2) 
  
  
  for (i in 1:napproach){ 
    
    # convergence
    tmpList <- unname(unlist(ListStats[[approaches[i]]][["conv"]][[cond]]))
    tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
    tmp <- table(tmpList)
    TableStats[cond, paste("conv", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
    
    if (TableStats[cond, paste("conv", approaches[i], sep="_")] > 0){
      
      # computation time
      tmp <- unname(unlist(ListStats[[approaches[i]]][["time"]][[cond]])) 
      TableStats[cond, paste("time", approaches[i], sep="_")] <- round(mean(tmp, na.rm=TRUE), 2) 
      TableStats[cond, paste("timesd", approaches[i], sep="_")] <- round(sd(tmp, na.rm=TRUE), 2) 
      
      ##### Parameter Estimates

      ### between (same in VR=1 and VR>1)
      B <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["B"]][[cond]])))
      
      # variances
      var_B <- B[,1:p] 
      
      # negative variances 
      tv <- c()
      for (rep in 1:nrep){
        tv[rep] <- any(var_B[rep,] < 0) 
      }
      tv <- factor(tv, levels=c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_B", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      RMSE_var_B <- mean( sqrt(colMeans(sweep(var_B, 2, popVar_B)^2)) )
      TableStats[cond, paste("RMSE_var_B", approaches[i], sep="_")] <- round(RMSE_var_B, 2) 
      relRMSE_var_B <- mean( (sqrt(colMeans(sweep(var_B, 2, popVar_B)^2))/popVar_B)  ) * 100
      TableStats[cond, paste("relRMSE_var_B", approaches[i], sep="_")] <- round( relRMSE_var_B, 2) 
      
      Bias_var_B <- mean( colMeans(sweep(var_B, 2, popVar_B)) )
      TableStats[cond, paste("Bias_var_B", approaches[i], sep="_")] <- round(Bias_var_B, 2)
      relBias_var_B <- mean( (colMeans(sweep(var_B, 2, popVar_B))/popVar_B)  ) * 100
      TableStats[cond, paste("relBias_var_B", approaches[i], sep="_")] <- round( relBias_var_B, 2)
      
      ## covariance
      cov_B <- as.data.frame(B[,(p+1):(p+c)])
      
      RMSE_cov_B <- mean( sqrt(colMeans(sweep(cov_B, 2, popCov_B)^2)) )
      TableStats[cond, paste("RMSE_cov_B", approaches[i], sep="_")] <- round(RMSE_cov_B, 2) 
      relRMSE_cov_B <- mean( (sqrt(colMeans(sweep(cov_B, 2, popCov_B)^2))/popCov_B)  ) * 100
      TableStats[cond, paste("relRMSE_cov_B", approaches[i], sep="_")] <- round( relRMSE_cov_B, 2) 
      
      Bias_cov_B <- mean( colMeans(sweep(cov_B, 2, popCov_B)) )
      TableStats[cond, paste("Bias_cov_B", approaches[i], sep="_")] <- round(Bias_cov_B, 2)
      relBias_cov_B <- mean( (colMeans(sweep(cov_B, 2, popCov_B))/popCov_B)  ) * 100
      TableStats[cond, paste("relBias_cov_B", approaches[i], sep="_")] <- round( relBias_cov_B, 2)
      
      # overall (variance and covariance)
      TableStats[cond, paste("RMSE_B", approaches[i], sep="_")] <- round( weighted.mean( c(RMSE_var_B, RMSE_cov_B), c(p, c) ), 2)
      TableStats[cond, paste("Bias_B", approaches[i], sep="_")] <- round( weighted.mean( c(Bias_var_B, Bias_cov_B), c(p, c) ), 2)
      TableStats[cond, paste("relRMSE_B", approaches[i], sep="_")] <- round( weighted.mean( c(relRMSE_var_B, relRMSE_cov_B), c(p, c) ), 2)
      TableStats[cond, paste("relBias_B", approaches[i], sep="_")] <- round( weighted.mean( c(relBias_var_B, relBias_cov_B), c(p, c) ), 2)
      
      # se and coverage
      se_B <- unname(unlist(ListStats[[approaches[i]]][["se_B"]][[cond]]))
      TableStats[cond, paste("se_B", approaches[i], sep = "_")] <- round(sum(se_B) /
                                                                             length(se_B) * 100, 2)
      coverage_B <- unname(unlist(ListStats[[approaches[i]]][["coverage_B"]][[cond]]))
      TableStats[cond, paste("coverage_B", approaches[i], sep="_")] <- round(sum(coverage_B)/length(coverage_B) * 100, 2)
      
      ## group 1
      
      W_1 <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["W_1"]][[cond]])))
      
      # extract all p variances
      var_W_1 <- W_1[, 1:p]
      
      # check for negative variances
      tv <- c()
      for (rep in 1:nrep) {
        tv[rep] <- any(var_W_1[rep, ] < 0)
      }
      tv <- factor(tv, levels = c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_W_1", approaches[i], sep = "_")] <- tmp[1] / (tmp[1] + tmp[2]) * 100
      
      RMSE_var_W_1 <- mean(sqrt(colMeans(sweep(
        var_W_1, 2, popVar_W_1
      ) ^ 2)))
      TableStats[cond, paste("RMSE_var_W_1", approaches[i], sep = "_")] <- round(RMSE_var_W_1, 2)
      relRMSE_var_W_1 <- mean((sqrt(colMeans(
        sweep(var_W_1, 2, popVar_W_1) ^ 2
      )) / popVar_W_1)) * 100
      TableStats[cond, paste("relRMSE_var_W_1", approaches[i], sep = "_")] <- round(relRMSE_var_W_1, 2)
      
      Bias_var_W_1 <- mean(colMeans(sweep(var_W_1, 2, popVar_W_1)))
      TableStats[cond, paste("Bias_var_W_1", approaches[i], sep = "_")] <- round(Bias_var_W_1, 2)
      relBias_var_W_1 <- mean((colMeans(sweep(
        var_W_1, 2, popVar_W_1
      )) / popVar_W_1)) * 100
      TableStats[cond, paste("relBias_var_W_1", approaches[i], sep = "_")] <- round(relBias_var_W_1, 2)
      
      # extract all c covariances
      cov_W_1 <- as.data.frame(W_1[, (p + 1):(p + c)])
      
      RMSE_cov_W_1 <- mean(sqrt(colMeans(sweep(
        cov_W_1, 2, popCov_W_1
      ) ^ 2)))
      TableStats[cond, paste("RMSE_cov_W_1", approaches[i], sep = "_")] <- round(RMSE_cov_W_1, 2)
      relRMSE_cov_W_1 <- mean((sqrt(colMeans(
        sweep(cov_W_1, 2, popCov_W_1) ^ 2
      )) / popCov_W_1)) * 100
      TableStats[cond, paste("relRMSE_cov_W_1", approaches[i], sep = "_")] <- round(relRMSE_cov_W_1, 2)
      
      Bias_cov_W_1 <- mean(colMeans(sweep(cov_W_1, 2, popCov_W_1)))
      TableStats[cond, paste("Bias_cov_W_1", approaches[i], sep = "_")] <- round(Bias_cov_W_1, 2)
      relBias_cov_W_1 <- mean((colMeans(sweep(
        cov_W_1, 2, popCov_W_1
      )) / popCov_W_1)) * 100
      TableStats[cond, paste("relBias_cov_W_1", approaches[i], sep = "_")] <- round(relBias_cov_W_1, 2)
      
      # overall (variance and covariance)
      RMSE_W_1 <- weighted.mean(c(RMSE_var_W_1, RMSE_cov_W_1), c(p, c))
      TableStats[cond, paste("RMSE_W_1", approaches[i], sep = "_")] <- round(RMSE_W_1, 2)
      Bias_W_1 <- weighted.mean(c(Bias_var_W_1, Bias_cov_W_1), c(p, c))
      TableStats[cond, paste("Bias_W_1", approaches[i], sep = "_")] <- round(Bias_W_1, 2)
      relRMSE_W_1 <- weighted.mean(c(relRMSE_var_W_1, relRMSE_cov_W_1), c(p, c))
      TableStats[cond, paste("relRMSE_W_1", approaches[i], sep = "_")] <- round(relRMSE_W_1, 2)
      relBias_W_1 <- weighted.mean(c(relBias_var_W_1, relBias_cov_W_1), c(p, c))
      TableStats[cond, paste("relBias_W_1", approaches[i], sep = "_")] <- round(relBias_W_1, 2)
      
      # se and coverage
      se_W_1 <- unname(unlist(ListStats[[approaches[i]]][["se_W_1"]][[cond]]))
      TableStats[cond, paste("se_W_1", approaches[i], sep = "_")] <- round(sum(se_W_1) /
                                                                             length(se_W_1) * 100, 2)
      coverage_W_1 <- unname(unlist(ListStats[[approaches[i]]][["coverage_W_1"]][[cond]]))
      TableStats[cond, paste("coverage_W_1", approaches[i], sep = "_")] <- round(sum(coverage_W_1) /
                                                                                   length(coverage_W_1) * 100, 2)
      
      # ICC
      ICC_1 <- var_B / (var_W_1 + var_B)
      RMSE_ICC_1 <- mean(sqrt(colMeans((ICC_1 - popICC_1) ^ 2)))
      TableStats[cond, paste("RMSE_ICC_1", approaches[i], sep = "_")] <- round(RMSE_ICC_1, 2)
      relRMSE_ICC_1 <- (RMSE_ICC_1 / popICC_1) * 100
      TableStats[cond, paste("relRMSE_ICC_1", approaches[i], sep = "_")] <- round(relRMSE_ICC_1, 2)
      Bias_ICC_1 <- mean(colMeans(ICC_1 - popICC_1))
      TableStats[cond, paste("Bias_ICC_1", approaches[i], sep = "_")] <- round(Bias_ICC_1, 2)
      relBias_ICC_1 <- (Bias_ICC_1 / popICC_1) * 100
      TableStats[cond, paste("relBias_ICC_1", approaches[i], sep = "_")] <- round(relBias_ICC_1, 2)
      
      ## group 2
      W_2 <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["W_2"]][[cond]])))
      
      # extract all p variances
      var_W_2 <- W_2[, 1:p]
      
      # check for negative variances
      tv <- c()
      for (rep in 1:nrep) {
        tv[rep] <- any(var_W_2[rep, ] < 0)
      }
      tv <- factor(tv, levels = c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_W_2", approaches[i], sep = "_")] <- tmp[1] / (tmp[1] + tmp[2]) * 100
      
      RMSE_var_W_2 <- mean(sqrt(colMeans(sweep(
        var_W_2, 2, popVar_W_2
      ) ^ 2)))
      TableStats[cond, paste("RMSE_var_W_2", approaches[i], sep = "_")] <- round(RMSE_var_W_2, 2)
      relRMSE_var_W_2 <- mean((sqrt(colMeans(
        sweep(var_W_2, 2, popVar_W_2) ^ 2
      )) / popVar_W_2)) * 100
      TableStats[cond, paste("relRMSE_var_W_2", approaches[i], sep = "_")] <- round(relRMSE_var_W_2, 2)
      
      Bias_var_W_2 <- mean(colMeans(sweep(var_W_2, 2, popVar_W_2)))
      TableStats[cond, paste("Bias_var_W_2", approaches[i], sep = "_")] <- round(Bias_var_W_2, 2)
      relBias_var_W_2 <- mean((colMeans(sweep(
        var_W_2, 2, popVar_W_2
      )) / popVar_W_2)) * 100
      TableStats[cond, paste("relBias_var_W_2", approaches[i], sep = "_")] <- round(relBias_var_W_2, 2)
      
      # extract all c covariances
      cov_W_2 <- as.data.frame(W_2[, (p + 1):(p + c)])
      
      RMSE_cov_W_2 <- mean(sqrt(colMeans(sweep(
        cov_W_2, 2, popCov_W_2
      ) ^ 2)))
      TableStats[cond, paste("RMSE_cov_W_2", approaches[i], sep = "_")] <- round(RMSE_cov_W_2, 2)
      relRMSE_cov_W_2 <- mean((sqrt(colMeans(
        sweep(cov_W_2, 2, popCov_W_2) ^ 2
      )) / popCov_W_2)) * 100
      TableStats[cond, paste("relRMSE_cov_W_2", approaches[i], sep = "_")] <- round(relRMSE_cov_W_2, 2)
      
      Bias_cov_W_2 <- mean(colMeans(sweep(cov_W_2, 2, popCov_W_2)))
      TableStats[cond, paste("Bias_cov_W_2", approaches[i], sep = "_")] <- round(Bias_cov_W_2, 2)
      relBias_cov_W_2 <- mean((colMeans(sweep(
        cov_W_2, 2, popCov_W_2
      )) / popCov_W_2)) * 100
      TableStats[cond, paste("relBias_cov_W_2", approaches[i], sep = "_")] <- round(relBias_cov_W_2, 2)
      
      # overall (variance and covariance)
      RMSE_W_2 <- weighted.mean(c(RMSE_var_W_2, RMSE_cov_W_2), c(p, c))
      TableStats[cond, paste("RMSE_W_2", approaches[i], sep = "_")] <- round(RMSE_W_2, 2)
      Bias_W_2 <- weighted.mean(c(Bias_var_W_2, Bias_cov_W_2), c(p, c))
      TableStats[cond, paste("Bias_W_2", approaches[i], sep = "_")] <- round(Bias_W_2, 2)
      relRMSE_W_2 <- weighted.mean(c(relRMSE_var_W_2, relRMSE_cov_W_2), c(p, c))
      TableStats[cond, paste("relRMSE_W_2", approaches[i], sep = "_")] <- round(relRMSE_W_2, 2)
      relBias_W_2 <- weighted.mean(c(relBias_var_W_2, relBias_cov_W_2), c(p, c))
      TableStats[cond, paste("relBias_W_2", approaches[i], sep = "_")] <- round(relBias_W_2, 2)
      
      # se and coverage
      se_W_2 <- unname(unlist(ListStats[[approaches[i]]][["se_W_2"]][[cond]]))
      TableStats[cond, paste("se_W_2", approaches[i], sep = "_")] <- round(sum(se_W_2) /
                                                                             length(se_W_2) * 100, 2)
      coverage_W_2 <- unname(unlist(ListStats[[approaches[i]]][["coverage_W_2"]][[cond]]))
      TableStats[cond, paste("coverage_W_2", approaches[i], sep = "_")] <- round(sum(coverage_W_2) /
                                                                                   length(coverage_W_2) * 100, 2)
      
      # ICC
      ICC_2 <- var_B / (var_W_2 + var_B)
      RMSE_ICC_2 <- mean(sqrt(colMeans((ICC_2 - popICC_2) ^ 2)))
      TableStats[cond, paste("RMSE_ICC_2", approaches[i], sep = "_")] <- round(RMSE_ICC_2, 2)
      relRMSE_ICC_2 <- (RMSE_ICC_2 / popICC_2) * 100
      TableStats[cond, paste("relRMSE_ICC_2", approaches[i], sep = "_")] <- round(relRMSE_ICC_2, 2)
      Bias_ICC_2 <- mean(colMeans(ICC_2 - popICC_2))
      TableStats[cond, paste("Bias_ICC_2", approaches[i], sep = "_")] <- round(Bias_ICC_2, 2)
      relBias_ICC_2 <- (Bias_ICC_2 / popICC_2) * 100
      TableStats[cond, paste("relBias_ICC_2", approaches[i], sep = "_")] <- round(relBias_ICC_2, 2)
      
      ## overall for both heteroscedastic groups
      
      # variances
      TableStats[cond, paste("RMSE_var_W", approaches[i], sep = "_")] <- round(mean(c(RMSE_var_W_1, RMSE_var_W_2)), 2)
      TableStats[cond, paste("relRMSE_var_W", approaches[i], sep = "_")] <- round(mean(c(relRMSE_var_W_1, relRMSE_var_W_2)), 2)
      TableStats[cond, paste("Bias_var_W", approaches[i], sep = "_")] <- round(mean(c(Bias_var_W_1, Bias_var_W_2)), 2)
      TableStats[cond, paste("relBias_var_W", approaches[i], sep = "_")] <- round(mean(c(relBias_var_W_1, relBias_var_W_2)), 2)
      
      TableStats[cond, paste("negVar_W", approaches[i], sep = "_")] <- round(mean(c(TableStats[cond, paste("negVar_W_1", approaches[i], sep =
                                                                                                             "_")], TableStats[cond, paste("negVar_W_2", approaches[i], sep = "_")])), 2)
      
      
      # covariances
      TableStats[cond, paste("RMSE_cov_W", approaches[i], sep = "_")] <- round(mean(c(RMSE_cov_W_1, RMSE_cov_W_2)), 2)
      TableStats[cond, paste("relRMSE_cov_W", approaches[i], sep = "_")] <- round(mean(c(relRMSE_cov_W_1, relRMSE_cov_W_2)), 2)
      TableStats[cond, paste("Bias_cov_W", approaches[i], sep = "_")] <- round(mean(c(Bias_cov_W_1, Bias_cov_W_2)), 2)
      TableStats[cond, paste("relBias_cov_W", approaches[i], sep = "_")] <- round(mean(c(relBias_cov_W_1, relBias_cov_W_2)), 2)
      
      # overall (variance and covariance)
      TableStats[cond, paste("RMSE_W", approaches[i], sep = "_")] <- round(mean(c(RMSE_W_1, RMSE_W_2)), 2)
      TableStats[cond, paste("relRMSE_W", approaches[i], sep = "_")] <- round(mean(c(relRMSE_W_1, relRMSE_W_2)), 2)
      TableStats[cond, paste("Bias_W", approaches[i], sep = "_")] <- round(mean(c(Bias_W_1, Bias_W_2)), 2)
      TableStats[cond, paste("relBias_W", approaches[i], sep = "_")] <- round(mean(c(relBias_W_1, relBias_W_2)), 2)
      
      # coverage
      TableStats[cond, paste("se_W", approaches[i], sep = "_")] <- round(mean(c(TableStats[cond, paste("se_W_1", approaches[i], sep =
                                                                                                         "_")], TableStats[cond, paste("se_W_2", approaches[i], sep = "_")])), 2)
      TableStats[cond, paste("coverage_W", approaches[i], sep = "_")] <- round(mean(c(TableStats[cond, paste("coverage_W_1", approaches[i], sep =
                                                                                                               "_")], TableStats[cond, paste("coverage_W_2", approaches[i], sep = "_")])), 2)
      
      # ICC
      TableStats[cond, paste("RMSE_ICC", approaches[i], sep = "_")] <- round(mean(c(RMSE_ICC_1, RMSE_ICC_2)), 2)
      TableStats[cond, paste("relRMSE_ICC", approaches[i], sep = "_")] <- round(mean(c(relRMSE_ICC_1, relRMSE_ICC_2)), 2)
      TableStats[cond, paste("Bias_ICC", approaches[i], sep = "_")] <- round(mean(c(Bias_ICC_1, Bias_ICC_2)), 2)
      TableStats[cond, paste("relBias_ICC", approaches[i], sep = "_")] <- round(mean(c(relBias_ICC_1, relBias_ICC_2)), 2)
      
    }
  }
}

saveRDS(TableStats, file = paste0(pathO, "TableStats.rds"))


### TableStats --> dat  ######################################################################################################################################################
# rearrange sumamrized table into data frame that can be used for figures

dat <- data.frame(cond = rep(1:ncond, napproach),
                  VR = rep(TableStats$VR, napproach),
                  var_B = rep(TableStats$var_B, napproach),
                  var_W_1 = rep(TableStats$var_W_1, napproach),
                  var_W_2 = rep(TableStats$var_W_2, napproach),
                  #p = rep(TableStats$p, napproach),
                  n = rep(TableStats$n, napproach),
                  g = rep(TableStats$g, napproach),
                  conv = unlist(TableStats[, grepl("conv", colnames(TableStats))], use.names=FALSE),
                  time = c(unlist(TableStats[, grepl("time_", colnames(TableStats))], use.names=FALSE)),
                  # variance
                  RMSE_var_B = unlist(TableStats[, paste0("RMSE_var_B_", approaches)], use.names=FALSE),
                  RMSE_var_W = unlist(TableStats[, paste0("RMSE_var_W_", approaches)], use.names=FALSE),
                  RMSE_var_W_1 = unlist(TableStats[, paste0("RMSE_var_W_1_", approaches)], use.names=FALSE),
                  RMSE_var_W_2 = unlist(TableStats[, paste0("RMSE_var_W_2_", approaches)], use.names=FALSE),
                  Bias_var_B = unlist(TableStats[, paste0("Bias_var_B_", approaches)], use.names=FALSE),
                  Bias_var_W = unlist(TableStats[, paste0("Bias_var_W_", approaches)], use.names=FALSE),
                  Bias_var_W_1 = unlist(TableStats[, paste0("Bias_var_W_1_", approaches)], use.names=FALSE),
                  Bias_var_W_2 = unlist(TableStats[, paste0("Bias_var_W_2_", approaches)], use.names=FALSE),
                  relRMSE_var_B = unlist(TableStats[, paste0("relRMSE_var_B_", approaches)], use.names=FALSE),
                  relRMSE_var_W = unlist(TableStats[, paste0("relRMSE_var_W_", approaches)], use.names=FALSE),
                  relRMSE_var_W_1 = unlist(TableStats[, paste0("relRMSE_var_W_1_", approaches)], use.names=FALSE),
                  relRMSE_var_W_2 = unlist(TableStats[, paste0("relRMSE_var_W_2_", approaches)], use.names=FALSE),
                  relBias_var_B = unlist(TableStats[, paste0("relBias_var_B_", approaches)], use.names=FALSE),
                  relBias_var_W = unlist(TableStats[, paste0("relBias_var_W_", approaches)], use.names=FALSE),
                  relBias_var_W_1 = unlist(TableStats[, paste0("relBias_var_W_1_", approaches)], use.names=FALSE),
                  relBias_var_W_2 = unlist(TableStats[, paste0("relBias_var_W_2_", approaches)], use.names=FALSE),
                  negVar_B = unlist(TableStats[, paste0("negVar_B_", approaches)], use.names=FALSE),
                  negVar_W = unlist(TableStats[, paste0("negVar_W_", approaches)], use.names=FALSE),
                  negVar_W_1 = unlist(TableStats[, paste0("negVar_W_1_", approaches)], use.names=FALSE),
                  negVar_W_2 = unlist(TableStats[, paste0("negVar_W_2_", approaches)], use.names=FALSE),
                  # covariance
                  RMSE_cov_B = unlist(TableStats[, paste0("RMSE_cov_B_", approaches)], use.names=FALSE),
                  RMSE_cov_W = unlist(TableStats[, paste0("RMSE_cov_W_", approaches)], use.names=FALSE),
                  RMSE_cov_W_1 = unlist(TableStats[, paste0("RMSE_cov_W_1_", approaches)], use.names=FALSE),
                  RMSE_cov_W_2 = unlist(TableStats[, paste0("RMSE_cov_W_2_", approaches)], use.names=FALSE),
                  Bias_cov_B = unlist(TableStats[, paste0("Bias_cov_B_", approaches)], use.names=FALSE),
                  Bias_cov_W = unlist(TableStats[, paste0("Bias_cov_W_", approaches)], use.names=FALSE),
                  Bias_cov_W_1 = unlist(TableStats[, paste0("Bias_cov_W_1_", approaches)], use.names=FALSE),
                  Bias_cov_W_2 = unlist(TableStats[, paste0("Bias_cov_W_2_", approaches)], use.names=FALSE),
                  relRMSE_cov_B = unlist(TableStats[, paste0("relRMSE_cov_B_", approaches)], use.names=FALSE),
                  relRMSE_cov_W = unlist(TableStats[, paste0("relRMSE_cov_W_", approaches)], use.names=FALSE),
                  relRMSE_cov_W_1 = unlist(TableStats[, paste0("relRMSE_cov_W_1_", approaches)], use.names=FALSE),
                  relRMSE_cov_W_2 = unlist(TableStats[, paste0("relRMSE_cov_W_2_", approaches)], use.names=FALSE),
                  relBias_cov_B = unlist(TableStats[, paste0("relBias_cov_B_", approaches)], use.names=FALSE),
                  relBias_cov_W = unlist(TableStats[, paste0("relBias_cov_W_", approaches)], use.names=FALSE),
                  relBias_cov_W_1 = unlist(TableStats[, paste0("relBias_cov_W_1_", approaches)], use.names=FALSE),
                  relBias_cov_W_2 = unlist(TableStats[, paste0("relBias_cov_W_2_", approaches)], use.names=FALSE),
                  # level
                  RMSE_B = unlist(TableStats[, paste0("RMSE_B_", approaches)], use.names=FALSE),
                  RMSE_W = unlist(TableStats[, paste0("RMSE_W_", approaches)], use.names=FALSE),
                  RMSE_W_1 = unlist(TableStats[, paste0("RMSE_W_1_", approaches)], use.names=FALSE),
                  RMSE_W_2 = unlist(TableStats[, paste0("RMSE_W_2_", approaches)], use.names=FALSE),
                  Bias_B = unlist(TableStats[, paste0("Bias_B_", approaches)], use.names=FALSE),
                  Bias_W = unlist(TableStats[, paste0("Bias_W_", approaches)], use.names=FALSE),
                  Bias_W_1 = unlist(TableStats[, paste0("Bias_W_1_", approaches)], use.names=FALSE),
                  Bias_W_2 = unlist(TableStats[, paste0("Bias_W_2_", approaches)], use.names=FALSE),
                  relRMSE_B = unlist(TableStats[, paste0("relRMSE_B_", approaches)], use.names=FALSE),
                  relRMSE_W = unlist(TableStats[, paste0("relRMSE_W_", approaches)], use.names=FALSE),
                  relRMSE_W_1 = unlist(TableStats[, paste0("relRMSE_W_1_", approaches)], use.names=FALSE),
                  relRMSE_W_2 = unlist(TableStats[, paste0("relRMSE_W_2_", approaches)], use.names=FALSE),
                  relBias_B = unlist(TableStats[, paste0("relBias_B_", approaches)], use.names=FALSE),
                  relBias_W = unlist(TableStats[, paste0("relBias_W_", approaches)], use.names=FALSE),
                  relBias_W_1 = unlist(TableStats[, paste0("relBias_W_1_", approaches)], use.names=FALSE),
                  relBias_W_2 = unlist(TableStats[, paste0("relBias_W_2_", approaches)], use.names=FALSE),
                  # se and coverage
                  se_W = unlist(TableStats[, paste0("se_W_", approaches)], use.names=FALSE),
                  se_W_1 = unlist(TableStats[, paste0("se_W_1_", approaches)], use.names=FALSE),
                  se_W_2 = unlist(TableStats[, paste0("se_W_2_", approaches)], use.names=FALSE),
                  se_B = unlist(TableStats[, paste0("se_B_", approaches)], use.names=FALSE),
                  coverage_W = unlist(TableStats[, paste0("coverage_W_", approaches)], use.names=FALSE),
                  coverage_W_1 = unlist(TableStats[, paste0("coverage_W_1_", approaches)], use.names=FALSE),
                  coverage_W_2 = unlist(TableStats[, paste0("coverage_W_2_", approaches)], use.names=FALSE),
                  coverage_B = unlist(TableStats[, paste0("coverage_B_", approaches)], use.names=FALSE),
                  # ICC
                  RMSE_ICC = unlist(TableStats[, paste0("RMSE_ICC_", approaches)], use.names=FALSE),
                  RMSE_ICC_1 = unlist(TableStats[, paste0("RMSE_ICC_1_", approaches)], use.names=FALSE),
                  RMSE_ICC_2 = unlist(TableStats[, paste0("RMSE_ICC_2_", approaches)], use.names=FALSE),
                  relRMSE_ICC = unlist(TableStats[, paste0("relRMSE_ICC_", approaches)], use.names=FALSE),
                  relRMSE_ICC_1 = unlist(TableStats[, paste0("relRMSE_ICC_1_", approaches)], use.names=FALSE),
                  relRMSE_ICC_2 = unlist(TableStats[, paste0("relRMSE_ICC_2_", approaches)], use.names=FALSE),
                  Bias_ICC = unlist(TableStats[, paste0("Bias_ICC_", approaches)], use.names=FALSE),
                  Bias_ICC_1 = unlist(TableStats[, paste0("Bias_ICC_1_", approaches)], use.names=FALSE),
                  Bias_ICC_2 = unlist(TableStats[, paste0("Bias_ICC_2_", approaches)], use.names=FALSE),
                  relBias_ICC = unlist(TableStats[, paste0("relBias_ICC_", approaches)], use.names=FALSE),
                  relBias_ICC_1 = unlist(TableStats[, paste0("relBias_ICC_1_", approaches)], use.names=FALSE),
                  relBias_ICC_2 = unlist(TableStats[, paste0("relBias_ICC_2_", approaches)], use.names=FALSE)
                  )


saveRDS(dat, paste0(pathO, "dat.rds"))

