#### Example Data and Models for Fig. 1
# Note that code that works only for n=2, p=2, and k=2. 
# You may change g, and some other characteristics, however.

library("lavaan")
library("tidyr")

## population characteristics
var_B <- 0.05
VR <- 2 # set larger than 1 if you want to model heteroscedasticity
cor <- 0.3 # correlation is the same at both levels
# don't change:
var_W_1 <- 1-var_B
var_W_2 <- var_W_1/VR
# ICC_1 <- var_B / (var_B + var_W_1)
# ICC_2 <- var_B / (var_B + var_W_2)
# meanICC <- mean(c(ICC_1, ICC_2))
cov_B <- var_B * cor
cov_W_1 <- var_W_1 * cor
cov_W_2 <- var_W_2 * cor

## sample characteristics
g <- 1000
n <- 2
p <- 2 
# p*n < p (cols < rows of WF) has to be satisfied

### data generation
set.seed(7742)

# generate data for each level 
popModel_B <- paste0("x1~~", var_B, "*x1; ", "x2~~", var_B, "*x2; ", "x1~~", cov_B, "*x2") 
popModel_W_1 <- paste0("x1~~", var_W_1, "*x1; ", "x2~~", var_W_1, "*x2; ", "x1~~", cov_W_1, "*x2")  
popModel_W_2 <- paste0("x1~~", var_W_2, "*x1; ", "x2~~", var_W_2, "*x2; ", "x1~~", cov_W_2, "*x2")  
sample_B <- simulateData(popModel_B, sample.nobs = g, model.type = "lavaan")
sample_W_1 <- simulateData(popModel_W_1, sample.nobs = n*(g/2), model.type = "lavaan")
sample_W_2 <- simulateData(popModel_W_2, sample.nobs = n*(g/2), model.type = "lavaan")
  
# merge data 
count <- 0 # unique index  
cluster <- rep(1:g, each = n)
data <- as.data.frame(matrix(NA, ncol=p, nrow=g*n))
for (k in 1:2) { # merge the sampled data from both levels
  sample_W <- get(paste0("sample_W_", k)) 
  for (i in 1:(n*g/2)){
    count <- count + 1
    j <- cluster[count]
    data[count,] <- sample_W[i,] + sample_B[j,]
  }
}
colnames(data) <- paste0("x", 1:p) 
# coding variables
data$unit <- rep(1:n, g) # unit in group (i)
data$cluster <- as.factor(cluster) # cluster (j)
data$k <- rep(c(1, 2), each = n * g / 2) # grouping variable
LF <- cbind(data[, (p + 1):(p + 3)], data[, 1:p]) # rearrange columns and rename
round(LF[,4:5], 0)
WF <- pivot_wider(LF, names_from = "unit", values_from = paste0("x", 1:p), names_sep = ".")
round(WF[,3:6], 0)

model_WFmultigroup <-
  "x1.1~~c(Vx1_w_g1, Vx1_w_g2)*x1.1; x1.2~~c(Vx1_w_g1, Vx1_w_g2)*x1.2;
  x2.1~~c(Vx2_w_g1, Vx2_w_g2)*x2.1; x2.2~~c(Vx2_w_g1, Vx2_w_g2)*x2.2;
  x1.1~~c(Cx12_w_g1, Cx12_w_g2)*x2.1; x1.2~~c(Cx12_w_g1, Cx12_w_g2)*x2.2;
  x1.1~0*1; x1.2~0*1; x2.1~0*1; x2.2~0*1;
  fx1=~1*x1.1+1*x1.2; fx2=~1*x2.1+1*x2.2; fx1~~fx1; fx2~~fx2; fx1~~fx2;
  fx1~1; fx2~1"

fit_WFmultigroup <- sem(
  model = model_WFmultigroup,
  data = WF,
  group = "k",
  group.equal = c("lv.variances", "lv.covariances")
)

# summary(fit_WFmultigroup)
S_WF_1 <- lavInspect(fit_WFmultigroup, "sampstat")$`1`$cov
S_WF_2 <- lavInspect(fit_WFmultigroup, "sampstat")$`2`$cov
round(S_WF_1, 2)
round(S_WF_2, 2)

## the other data matrices and sample covariance matrices 
round(LF[,4:(4+p-1)], 0) 
S_LF <- cov(LF[,4:(4+p-1)])
round(S_LF, 2)

round(WF[,3:(3+(p*n)-1)], 0)
S_WF <- cov(WF[,3:(3+(p*n)-1)])*(g-1)/g # "biased" estimator with g in denominator (i.e., no Bessel's correction, g-1) as in lavaan
round(S_WF, 2)
