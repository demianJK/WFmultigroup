# this function analyses data in theWFmultigroup approach with an intercept-only model (variances, covariances, means)
# suited only for two groups (k=2)

analyseData <- function(data_WF, VR, n, g, p){
  
  ## 1) set model specification (different for within and between)
  
  ## within (p * n)
  
  # variances (p * n)
  tmp2 <- c()
  resid_w <- c() # residual variances (p * n - equality among n blocks)
  tmp3 <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp2[i] <- paste0("x", j, ".", i)
      tmp3[i] <- paste0(tmp2[i], "~~c(Vx", j, "_w_g1, Vx", j, "_w_g2)*", tmp2[i])  ## changed
    } 
    resid_w[j] <- paste(tmp3, collapse="; ")
  }
  resid_w <- paste(resid_w, collapse="; ")
  
  # covariances (p-pc=c*n with c-wise equality constraints)
  resid_cov <- c() # manifest correlations (p * n - equality among n blocks)
  count <- 0
  for (i in 1:n){ # n-unit-wise ordering
    for(j in 1:p){   
      for(m in 1:p){
        if(j != m & m > j){
          count <- count + 1
          resid_cov[count] <- paste0("x", j, ".", i, "~~c(Cx", j, m, "_w_g1, Cx", j, m, "_w_g2)*", "x", m, ".", i)
        } 
      } 
    }
  }
  resid_cov <- paste(resid_cov, collapse="; ")
  
  # "If variables are both at the within- and the between-level, the intercepts at the within-level should be fixed to zero." (p.701f, Barendse & Rosseel, 2020)
  fac_int_w <- c()
  tmp <- c()
  count <- 0
  for (j in 1:p){
    for (i in 1:n){
      count <- count + 1
      tmp[count] <- paste0("x", j, ".", i, "~0*1")
    }
  }
  fac_int_w <- paste(tmp, collapse = "; ")
  
  model_WFmultigroup_W <- paste(resid_w, resid_cov, fac_int_w, sep = "; ")
  
  ## between (p)
  
  fac_b <- c() # latent factor = random intercepts (p - loadings fixed to 1)
  tmp <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp[i] <- paste0("1*x", j, ".", i)
    }
    fac_b[j] <- paste0("fx", j, "=~", paste(tmp, collapse="+"))
  }
  fac_b <- paste(fac_b, collapse="; ")
  
  # variances
  fac_var_b <- c() # factor variance (p)
  fac_int_b <- c() # to estimate means
  for (j in 1:p){
    fac_var_b[j] <- paste0("fx", j, "~~fx", j)
    fac_int_b[j] <- paste0("fx", j, "~1")
  }
  fac_var_b <- paste(fac_var_b, collapse="; ")
  fac_int_b <- paste(fac_int_b, collapse="; ")
  
  # covariances
  fac_cov_b <- c() # correlations between factors
  count <- 0
  for(j in 1:p){   
    for(m in 1:p){
      if(j != m & m > j){
        count <- count + 1
        fac_cov_b[count] <- paste0("fx", j, "~~", "fx", m)
      }
    } 
  }
  fac_cov_b <- paste(fac_cov_b, collapse = "; ")
  
  model_WFmultigroup_B <- paste(fac_b, fac_var_b, fac_cov_b, fac_int_b, sep="; ")
  
  model_WFmultigroup <- paste(model_WFmultigroup_W, model_WFmultigroup_B, sep="; ")
  

  fit_WFmultigroup <- try(sem(model = model_WFmultigroup, # Barendse & Rossel (2020) use sem() as well (see Appendix A, p.718)
                    #likelihood = "normal", # S, SE, and test statistics based on /(g)
                    data = data_WF,
                    group="k",
                    group.equal = c("lv.variances", "lv.covariances")),
                silent = TRUE)
  
  
  
  return(list(fit_WFmultigroup=fit_WFmultigroup) ) 
}

# Barendse, M. T., & Rosseel, Y. (2020). Multilevel Modeling in the ‘Wide Format’ Approach with Discrete Data: A Solution for Small Cluster Sizes. Structural Equation Modeling: A Multidisciplinary Journal, 27(5), 696–721. https://doi.org/10.1080/10705511.2019.1689366
