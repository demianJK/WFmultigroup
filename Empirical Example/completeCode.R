###### (0) Prerequisites

## load required packages
library("dplyr") # select and filter data (version 1.1.4)
library("foreign") # read SPSS (version 0.8-87)
library("ggplot2") # figures (version 3.5.1)
library("huxtable") # APA table (version 5.5.6)
library("lavaan") # ML MG SEM (version 0.6-18)
# Note that this CRAN version of lavaan does not yield the same results in the homogeneous model in the "genuine" ML MG SEM approach
# as the WFmultigroup approach does. However, the most recent version on Github (0.6-19.2187) does so.
# install.packages("devtools")
# library("devtools")
# install_github("yrosseel/lavaan")
library("lme4") # logistic regression of missingness (version 1.1-35.5)
library("mice") # multiple imputation (version 3.16.0)
library("naniar") # MCAR test (version 1.1.0)
library("patchwork") # combining ggplots by + (version 1.2.0)
library("psych") # descriptive stats (version 2.4.6.26)
library("tidyr") # reformating (version 1.3.1)

## load data

# Go to https://www.oecd.org/pisa/data/2022database/ 
# Navigate to SPSS (TM) Data Files (compressed) >>> Student Questionnaire data file and download the file
PISA <- read.spss("../CY08MSP_STU_QQQ.SAV", to.data.frame=TRUE, use.value.labels = FALSE) # otherwise numerical vectors might be handled as factors
# the data frame is in LF (i.e., each row corresponds to a student)

# If you don't want to run the multiple imputation, simply load the final data frame and continue in line 409.
PISA_short_balanced_imp <- read.csv(file = "../PISA_short_balanced_imp.csv")



##### (1) Data Subsetting

## select relevant variables
PISA_short <- select(PISA, 
                     CNTSTUID, # unique student ID (level-1)
                     CNTSCHID, # school (level-2)
                     CNT, # CNT (group)
                     CREATAS, # Creative Activities at school
                     GROSAGR # Growth Mindset
) 
# PISA_short is "LF unbalanced"

## select relevant cases (Albania and Ireland) of between-cluster variable country
PISA_short <- filter(PISA_short, CNT == "ALB" | CNT == "IRL")



##### (2) Inspecting the Data I: Data Structure and Data Types

## inspect data structure and data types
str(PISA_short)

# is it not necessary to factorise the discrete ID indicators CNTSTUID and CNTSCHID...

## ... but we recode the grouping variable for the figures
PISA_short$CNT <- ifelse(PISA_short$CNT == "IRL", yes="Ireland", no="Albania")
# (we do not factorise bc otherwise we would introduce problems with data subsetting and multiple imputation later on)



##### (3) Inspecting the Data II: Unbalanced Cluster Sizes

## get information on the selected subsample
N <- nrow(PISA_short)
schools <- unique(PISA_short$CNTSCHID)
g <- length(schools) 
n <- as.vector(table(PISA_short$CNTSCHID))
n_mean <- mean(n)
n_min <- min(n)
n_max <- max(n)

country <- c()
for (j in 1:g){
  country[j] <-  unique(PISA_short$CNT[PISA_short$CNTSCHID == schools[j]]) 
}

nData <- data.frame(country = country,
                    school = schools, 
                    n = n)

a <- ggplot(nData, aes(x=g, fill=country)) + 
  geom_bar(data = transform(nData, country = NULL), fill = "grey85") +
  geom_bar(show.legend = FALSE) + facet_grid(. ~ country) +
  scale_y_continuous(name="Frequency", expand=c(0,0)) +
  scale_x_discrete(name="Number of Schools (g)",) +
  scale_fill_manual(values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")
# table(country)

b <- ggplot(nData, aes(x=n, fill=country)) + 
  geom_histogram(data = transform(nData, country = NULL), fill = "grey85", binwidth=1) +
  geom_histogram(binwidth=1, show.legend = FALSE) + facet_grid(country ~ .) +
  scale_y_continuous(name="Number of Schools (g)", expand=c(0,0)) +
  scale_x_continuous(name="School Size (n)", expand=c(0.01,0.01), limits=c(0, NA)) +
  scale_fill_manual(values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

a + b # Fig.3

# N=11.698 with g=444 and the distribution of cluster sizes (n) differs fairly.
# country-wise:
table(PISA_short$CNT) # N
table(nData$country) # g



##### (4) Inspecting the Data III: Distribution of Variables

## Raw Data

# univariate
a <- ggplot(PISA_short, aes(x=CREATAS, fill=CNT)) + 
  geom_histogram(show.legend = FALSE, position = "identity", alpha=0.5) + 
  scale_x_continuous(name="Creative Activities at School (CREATAS)", expand=c(0,0), limits=c(-6, 6)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 2000)) +
  scale_fill_manual(values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="Raw Data", subtitle="A")

b <- ggplot(PISA_short, aes(x=GROSAGR, fill=CNT)) + 
  geom_histogram(show.legend = FALSE, position = "identity", alpha=0.5) + 
  scale_x_continuous(name="Growth Mindset (GROSAGR)", expand=c(0,0), limits=c(-6, 6)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 2000)) +
  scale_fill_manual(values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(subtitle="B")

# bivariate
c <- ggplot(PISA_short, aes(x=CREATAS, y=GROSAGR, col=CNT)) + 
  geom_point(show.legend = FALSE, alpha=0.3) + 
  scale_x_continuous(name="Creative Activities at School (CREATAS)", expand=c(0, 0), limits=c(-6.5, 6.5)) +
  scale_y_continuous(name="Growth Mindset (GROSAGR)", expand=c(0,0), limits=c(-6, 6)) +
  scale_color_manual(values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(subtitle="C")

## Cluster means

# estimate cluster means and create data frame
CREATAS_cluster_means <- aggregate(PISA_short$CREATAS, list(PISA_short$CNTSCHID), FUN=mean, na.rm=TRUE, na.action=NULL)
GROSAGR_cluster_means <- aggregate(PISA_short$GROSAGR, list(PISA_short$CNTSCHID), FUN=mean, na.rm=TRUE, na.action=NULL)

PISA_short <- PISA_short[order(PISA_short$CNTSCHID),]

j <- c()
country <- c()
for (i in 1:nrow(PISA_short)){
  tmp_j <- PISA_short$CNTSCHID[i]
  if (i==1){
    country <- append(country, PISA_short$CNT[i])
    j <- append(j, tmp_j)
  } else {
    if (tmp_j > tail(j, n=1)){
      country <- append(country, PISA_short$CNT[i])
      j <- append(j, tmp_j )
    }
  }
}

PISA_short_cluster_means <- data.frame(j=1:444, country=country, CREATAS=CREATAS_cluster_means$x, GROSAGR=GROSAGR_cluster_means$x)

# univariate
d <- ggplot(PISA_short_cluster_means, aes(x=CREATAS, fill=country)) + 
  geom_histogram(show.legend = FALSE, position = "identity", alpha=0.5) + 
  scale_x_continuous(name="Creative Activities at School (CREATAS)", expand=c(0, 0), 
                     limits=c(-6, 6)) +
  scale_y_continuous(name="Frequency", limits=c(0, 150), expand=c(0, 0),) +
  scale_fill_manual(name="Country", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(title="Cluster Means", subtitle="D")


e <- ggplot(PISA_short_cluster_means, aes(x=GROSAGR, fill=country)) + 
  geom_histogram(show.legend = FALSE, position = "identity", alpha=0.5) + 
  scale_x_continuous(name="Growth Mindset (GROSAGR)", expand=c(0,0), 
                     limits=c(-6, 6)) +
  scale_y_continuous(name="Frequency", limits=c(0, 150), expand=c(0, 0),) +
  scale_fill_manual(name="Country", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(subtitle="E")

# bivariate
f <- ggplot(PISA_short_cluster_means, aes(x=CREATAS, y=GROSAGR, col=country)) +
  geom_point(alpha=0.3) +
  scale_x_continuous(name="Creative Activities at School (CREATAS)", expand=c(0, 0),
                     limits=c(-6, 6)
  ) +
  scale_y_continuous(name="Growth Mindset (GROSAGR)", expand=c(0,0),
                     limits=c(-6, 6)
  ) +
  scale_color_manual(name="Country", values=c("#002654", "#ffce00")) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))  +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(subtitle="F")

a + b + c + d + e + f + plot_layout(nrow=2, guides='collect') & theme(text = element_text("serif"), legend.position = "bottom") # Fig.4



##### (5) Inspecting the Data IV: Missing Data

## What is the proportion of missingness?
vis_miss(PISA_short)
# 28% of CREATAS and 20% of GROSAGR missing
# for each country:
table(is.na(PISA_short$CREATAS), PISA_short$CNT)#/nrow(PISA_short) 
table(is.na(PISA_short$GROSAGR), PISA_short$CNT)#/nrow(PISA_short)
# numbers from footnote Fig.3

## Is the missingness systematical?
# MCAR: missings are completely independent of other variables and the missing value itself
# MAR: missings are dependent on other variables but not on the missing itself
# MNAR: missings are independent of the other variables but they are not random

## Let's check the missing patterns (= co-occurence of missings in multiple variables).

## (a) descriptive 
# by figure with percentages
md.pattern(PISA_short, rotate.names = TRUE) # note this function is from package mice but mcar_test is from package naniar

# rows: missing patterns
# numbers to left: cases for each missing pattern
# number to right: number of missings in missing pattern
# numbers at bottom: number of missing cases for each variable (column) --> absolute numbers we got in figure before

# 4 patterns
# most often all variables existent (1. row), 
8137 / (8137 + 1182 + 312 + 2067) # approx. 70% cases without any missings, thus, 30% of cases with at least one missing!
# then one missing in CREATAS (2. row), 
(1182) / (8137 + 1182 + 312 + 2067) # approx 10% of only missing CREATAS
# then missings in CREATAS and GROSAGR (4. row), 
(2067) / (8137 + 1182 + 312 + 2067) # approx 18% of missing CREATAS and GROSAGR
# Note 10% + 18% add up to the 28% missing cases reported for CREATAS before
# then one missing in GROSAGR (3. row)
(312) / (8137 + 1182 + 312 + 2067) # approx 3% of only missing GROSAGR


## (b) inferential
# by using Little's (1988) test that compares patterns of missingness
# H0: MCAR
# H1: not MCAR
# Note CNT and CNTSCHID are perfectly correlated and can thus not be used in the same test bc of multicollinearity (i.e., singularity)
# we drop CNT
mcar_test(PISA_short[, c("CNTSCHID", "CREATAS", "GROSAGR")]) 
# test is significant, thus evidence that MCAR does not hold

## explore MAR assumption

# create missing data indicators (missing=1, existent=0)
PISA_short$missing_CREATAS <- ifelse(is.na(PISA_short$CREATAS), yes=1, no=0) 
PISA_short$missing_GROSAGR <- ifelse(is.na(PISA_short$GROSAGR), yes=1, no=0)

## (a) descriptive 
# by correlation table
cor_data <- PISA_short
cor_data$CNT <- ifelse(cor_data$CNT == "Albania", yes=1, no=0) # recode to numeric bc character does not work
cor <- cor(cor_data, use = "pairwise.complete.obs")
cor[upper.tri(cor)] <- NA
print(round(cor, 2), na.print="") 

# missingness has large correlation with country (0.393 and 0.431)
# missingness has large correlation with cluster (-0.393 and -0.431)
# contingency of missingness (or presence) of both variables is quite large (0.667), we see this in the missing patterns
# together, this suggest a design effect (i.e., questionnaires not administered in certain clusters in countries)
# missingness has small correlation with other variable (-0.065 and 0.120)
# most importantly, country has moderate to large correlation with the other variable (0.425 and -0.235)

## (b) inferential 
# by fitting logistic mixed-effects models to predict missingness
# Note that a variable and their missingness indicator cannot be used in the same model because of multicollinearity (e.g. GROSAGR and missing_GROSAGR).
# Thus, we consider one model for each.

# CREATAS
model_CREATAS <- glmer(missing_CREATAS ~ CNT * GROSAGR + (1 | CNTSCHID), family = binomial, data = PISA_short)
summary(model_CREATAS)
# CNT and GROSAGR predict NA in CREATAS
model_CREATAS_mi <- glmer(missing_CREATAS ~ CNT * missing_GROSAGR + (1 | CNTSCHID), family = binomial, data = PISA_short)
summary(model_CREATAS_mi)
# CNT, NA in GROSAGR, and their interaction predict NA in CREATAS

# GROSAGR
model_GROSAGR <- glmer(missing_GROSAGR ~ CNT * CREATAS + (1 | CNTSCHID), family = binomial, data = PISA_short)
summary(model_GROSAGR)
# CNT predicts NA in CREATAS
model_GROSAGR_mi <- glmer(missing_GROSAGR ~ CNT * missing_CREATAS + (1 | CNTSCHID), family = binomial, data = PISA_short)
summary(model_GROSAGR_mi)
# CNT, NA in CREATAS, and their interaction predict NA in GROSAGR

# evidence for MAR: missingness can be predicted by other variables (or missingness of other variables) in data and country
# thus imputation is warranted, but first we inspect another source of missingness and estimation problems



##### (6) Reformating I: Balanced Cluster Sizes in LF
# necessary for imputing unbalanced data, and to reformat to WF later

## create new data frame with balanced number of students 
PISA_short_balanced <- data.frame(
  j = rep(1:g, each=n_max),
  i = rep(1:n_max, times=g),
  CNTSCHID = rep(NA, n_max*g) , # incomplete
  CNTSTUID = rep(NA, n_max*g), # incomplete
  CNT = rep(NA, n_max*g),
  CREATAS = rep(NA, n_max*g),
  GROSAGR = rep(NA, n_max*g),
  missing_CREATAS = rep(1, n_max*g),
  missing_GROSAGR = rep(1, n_max*g)
)

# sort data by school
PISA_short <- PISA_short[with(PISA_short, order(CNTSCHID)), ]

# fill in existing data
for (j in 1:g) {
  school <- unique(PISA_short$CNTSCHID)[j]
  students <- filter(PISA_short, CNTSCHID == school)$CNTSTUID
  nSchool <- length(students)
  PISA_short_balanced$CNTSCHID[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- school
  PISA_short_balanced$CNTSTUID[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- students
  PISA_short_balanced$CNT[((j - 1) * n_max + 1):((j - 1) * n_max + n_max)] <- unique(PISA_short$CNT[which(PISA_short$CNTSCHID == school)]) 
  PISA_short_balanced$CREATAS[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- PISA_short$CREATAS[which(PISA_short$CNTSCHID == school)]
  PISA_short_balanced$GROSAGR[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- PISA_short$GROSAGR[which(PISA_short$CNTSCHID == school)]
  PISA_short_balanced$missing_CREATAS[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- PISA_short$missing_CREATAS[which(PISA_short$CNTSCHID == school)]
  PISA_short_balanced$missing_GROSAGR[((j - 1) * n_max + 1):((j - 1) * n_max + nSchool)] <- PISA_short$missing_GROSAGR[which(PISA_short$CNTSCHID == school)]
}

# Now N=n_max*g = 19980 level-1 units.
# Final subsample per country: g*n_max = N 
table(nData$country)*n_max 

# "genuine" missings and unbalanced data
table(PISA_short_balanced$missing_CREATAS, PISA_short_balanced$CNT) 
table(PISA_short_balanced$missing_GROSAGR, PISA_short_balanced$CNT) 
# numbers from footnote Fig.4



##### (7) Multiple Imputation
# in LF and country-wise

# set imputation method for CREATAS and GROSAGR
meth <- mice(PISA_short_balanced, maxit = 0)$method
meth["CNTSCHID"] <- "" 
meth[c("CREATAS", "GROSAGR")] <- "2l.pan"  # homogeneous variances in each group (i.e., country) assumed 

# create imputation models for CREATAS and GROSAGR
pred <- make.predictorMatrix(PISA_short_balanced)
pred[ , "j"] <- -2  # Set cluster variable
pred[c("j", "i", "CNTSCHID", "CNTSTUID", "CNT", "missing_CREATAS", "missing_GROSAGR"), ] <- 0  # no models for these variables
pred[ , c("i", "CNTSCHID", "CNTSTUID", "CNT", "missing_CREATAS", "missing_GROSAGR") ] <- 0  # not used as predictors ###### no CNT

# impute
imp_Albania <- mice(filter(PISA_short_balanced, CNT == "Albania"), predictorMatrix = pred, method = meth, seed = 123)
imp_Ireland <- mice(filter(PISA_short_balanced, CNT == "Ireland"), predictorMatrix = pred, method = meth, seed = 123)

# inspect single imputed data sets
stripplot(imp_Albania, CREATAS, pch = 19, xlab = "Imputation number")
stripplot(imp_Ireland, CREATAS, pch = 19, xlab = "Imputation number")
stripplot(imp_Albania, GROSAGR, pch = 19, xlab = "Imputation number")
stripplot(imp_Ireland, GROSAGR, pch = 19, xlab = "Imputation number")
# Because the imputed data sets appear quite similar, we will combine them instead of estimating models for each
# data set and pooled the results. 

# compare descriptive stats of existent and imputed data (Tab.1)
ex_Alb <- describe(select(PISA_short_balanced[PISA_short_balanced$CNT == "Albania",], CREATAS, GROSAGR))
imp_Alb <- describe(select(complete(imp_Albania), CREATAS, GROSAGR))
ex_Ire <- describe(select(PISA_short_balanced[PISA_short_balanced$CNT == "Ireland",], CREATAS, GROSAGR))
imp_Ire <- describe(select(complete(imp_Ireland), CREATAS, GROSAGR))
# for both countries, mean and sd are quite similar in the existent and imputed data

# combine imputed data sets of both groups (i.e., countries)
PISA_short_balanced_imp <- rbind(complete(imp_Albania), complete(imp_Ireland))

# plot imputed data (Fig.5)
ggplot(PISA_short_balanced_imp, aes(x=CREATAS, y=GROSAGR, col=CNT)) + 
  geom_point(data = transform(PISA_short_balanced_imp, CNT = NULL), col="grey85", alpha=0.5) +
  geom_point(show.legend = FALSE, alpha=0.3) + 
  facet_grid(CNT ~ missing_CREATAS, margins=TRUE, # adds an additional facet for all levels combined
             labeller=as_labeller(c('0'="Existent", '1'="Imputed", '(all)'="All", 'Albania'="Albania", 'Ireland'="Ireland"))) +
  scale_x_continuous(name="Creative Activities at School (CREATAS)", expand=c(0, 0), limits=c(-6.5, 6.5)) +
  scale_y_continuous(name="Growth Mindset (GROSAGR)", expand=c(0.05,0.05), limits=c(-6, 6)) +
  scale_color_manual(values=c("#002654", "#ffce00", "black")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))



##### (8) Reformating II: Format LF to WF
# where each row corresponds to a school
PISA_short_balanced_imp_WF <- select(PISA_short_balanced_imp, -c("CNTSCHID", "CNTSTUID", "missing_CREATAS", "missing_GROSAGR")) # drop variables, otherwise formating faulty
PISA_short_balanced_imp_WF <- pivot_wider(PISA_short_balanced_imp_WF, names_from = i, values_from = c("CREATAS", "GROSAGR"), names_sep = ".") 



##### (9) Model Estimation 

## Homogeneity/Heterogeneity is set differently for both levels:
# Level-1: in model syntax (by using same or different parameter labels)
# Level-2: with function parameter “group.equal = c("lv.variances", "lv.covariances")” (by setting it or leaving it out)
# Thus, the model syntax below is the same for all models.
# (Note that the variances at each level in the homogeneous models equal the pooled variances in the heterogeneous models.)

varNames <- c("CREATAS", "GROSAGR") # variable names in vector required for loop
p <- length(varNames)

## within (p*n)

# means set to 0
means_w <- c()
tmp <- c()
count <- 0
for (j in 1:p){
  for (i in 1:n_max){
    count <- count + 1
    tmp[count] <- paste0(varNames[j], ".", i, "~0*1")
  }
}
means_w <- paste(tmp, collapse = "; ")


## between (p)

# factor loadings
fac_load_b <- c() 
tmp <- c()
for (j in 1:p){
  for (i in 1:n_max){
    tmp[i] <- paste0("1*", varNames[j], ".", i)
  }
  fac_load_b[j] <- paste0("f", varNames[j], " =~", paste(tmp, collapse="+"))
}
fac_load_b <- paste(fac_load_b, collapse="; ")

# variances and means
fac_var_b <- c() 
fac_int_b <- c()
for (j in 1:p){
  fac_var_b[j] <- paste0("f", varNames[j], "~~f", varNames[j])
  fac_int_b[j] <- paste0("f", varNames[j], "~1")
}
fac_var_b <- paste(fac_var_b, collapse="; ")
fac_int_b <- paste(fac_int_b, collapse="; ")

# covariances
fac_cov_b <- c() 
count <- 0
for(j in 1:p){   
  for(m in 1:p){
    if(j != m & m > j){
      count <- count + 1
      fac_cov_b[count] <- paste0("f", varNames[j], "~~", "f", varNames[m])
    }
  } 
}
fac_cov_b <- paste(fac_cov_b, collapse = "; ")

model_WF_B <- paste(fac_load_b, fac_var_b, fac_cov_b, fac_int_b, sep="; ")


### Model with Homogeneous Within- and Between-Cluster (Co)variances

## within (p*n)

# variances 
tmp2 <- c()
resid_var_w_homo <- c() 
tmp3 <- c()
for (j in 1:p){
  for (i in 1:n_max){
    tmp2[i] <- paste0(varNames[j], ".", i)
    tmp3[i] <- paste0(tmp2[i], "~~c(", varNames[j], "_both, ", varNames[j], "_both)*", tmp2[i]) # same label for parameter ACROSS groups
  }
  resid_var_w_homo[j] <- paste(tmp3, collapse="; ")
}
resid_var_w_homo <- paste(resid_var_w_homo, collapse="; ")

# covariances 
resid_cov_w_homo <- c() 
count <- 0
for (i in 1:n_max){ 
  for(j in 1:p){   
    for(m in 1:p){
      if(j != m & m > j){
        count <- count + 1
        resid_cov_w_homo[count] <- paste0(varNames[j], ".", i, "~~c(", varNames[j], "_", varNames[m], "_both, ", varNames[j], "_", varNames[m], "_both)*", varNames[m], ".", i) # same label for parameter ACROSS groups
      } 
    } 
  }
}
resid_cov_w_homo <- paste(resid_cov_w_homo, collapse="; ")

model_WF_W_homo <- paste(resid_var_w_homo, resid_cov_w_homo, means_w, sep = "; ")

model_WFmultigroup_homo <- paste(model_WF_W_homo, model_WF_B, sep="; ")

fit_WFmultigroup_homo <- sem(model = model_WFmultigroup_homo,
                             data = PISA_short_balanced_imp_WF,
                             group="CNT",
                             group.equal = c("lv.variances", "lv.covariances"))
summary(fit_WFmultigroup_homo)


## the "genuine" lavaan ML MG SEM returns very similar estimates in the most recent version on Github (0.6-19.2186).
# Note that here only the function parameter "group.equal" controls whether a fully homogeneous/heterogeneous model is estimated.

model_MLMGSEM <- c(
  "
  Group: 1
  Level: 1
  CREATAS ~~ CREATAS
  GROSAGR ~~ GROSAGR
  CREATAS ~~ GROSAGR
  Level: 2
  CREATAS ~~ CREATAS
  GROSAGR ~~ GROSAGR
  CREATAS ~~ GROSAGR
  
  Group: 2
  Level: 1
  CREATAS ~~ CREATAS
  GROSAGR ~~ GROSAGR
  CREATAS ~~ GROSAGR
  Level: 2
  CREATAS ~~ CREATAS
  GROSAGR ~~ GROSAGR
  CREATAS ~~ GROSAGR
  "
) # Note that the same model syntax is used for the fully heterogeneous model later.
# Alternatively, one could use parameter labels to denote homogeneous (i.e., same label in both groups) or heterogeneous (i.e., differen labels in both groups) parameters 
# (just as in the WFmultigroup approach; see also the models that are heterogeneous at one level in the genuine ML MG SEM approach).

fit_MLMGSEM_homo <- sem(model = model_MLMGSEM,
                      data = PISA_short_balanced_imp, # data in LF!
                      cluster="j", 
                      group="CNT",
                      group.equal = c("residuals", "residual.covariances") # homogeneous
)
summary(fit_MLMGSEM_homo)



### Model with Heterogeneous Within-Cluster (Co)variances

## within (p*n)

# variances (with n-wise equality constraints)
tmp2 <- c() 
tmp3 <- c() 
resid_var_w_hetero <- c() 
for (j in 1:p){
  for (i in 1:n_max){
    tmp2[i] <- paste0(varNames[j], ".", i)
    tmp3[i] <- paste0(tmp2[i], "~~c(", varNames[j], "_albania, ", varNames[j], "_ireland)*", tmp2[i]) # same label for parameter WITHIN groups
  }
  resid_var_w_hetero[j] <- paste(tmp3, collapse="; ")
}
resid_var_w_hetero <- paste(resid_var_w_hetero, collapse="; ")

# covariances (with n-wise equality constraints)
resid_cov_w_hetero <- c() 
count <- 0
for (i in 1:n_max){ 
  for(j in 1:p){   
    for(m in 1:p){
      if(j != m & m > j){
        count <- count + 1
        resid_cov_w_hetero[count] <- paste0(varNames[j], ".", i, "~~c(", varNames[j], "_", varNames[m], "_albania, ", varNames[j], "_", varNames[m], "_ireland)*", varNames[m], ".", i) # same label for parameter WITHIN groups
      } 
    } 
  }
}
resid_cov_w_hetero <- paste(resid_cov_w_hetero, collapse="; ")

model_WF_W_hetero <- paste(resid_var_w_hetero, resid_cov_w_hetero, means_w, sep = "; ")

model_WFmultigroup_hetero_W <- paste(model_WF_W_hetero, model_WF_B, sep="; ")

fit_WFmultigroup_hetero_W <- sem(model = model_WFmultigroup_hetero_W,
                               data = PISA_short_balanced_imp_WF,
                               group="CNT",
                               group.equal = c("lv.variances", "lv.covariances")) 
summary(fit_WFmultigroup_hetero_W)


## Here you can use the CRAN version of lavaan (0.6-18) with its "genuine" ML MG SEM which yields very similar estimates.

model_MLMGSEM_hetero_W <- c(
  "
  Group: 1
  Level: 1
  CREATAS ~~ CREATAS_albania*CREATAS
  GROSAGR ~~ GROSAGR_albania*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_albania*GROSAGR
  Level: 2
  CREATAS ~~ CREATAS_both*CREATAS
  GROSAGR ~~ GROSAGR_both*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_both*GROSAGR
  
  Group: 2
  Level: 1
  CREATAS ~~ CREATAS_ireland*CREATAS
  GROSAGR ~~ GROSAGR_ireland*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_ireland*GROSAGR
  Level: 2
  CREATAS ~~ CREATAS_both*CREATAS
  GROSAGR ~~ GROSAGR_both*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_both*GROSAGR
  "
)

fit_MLMGSEM_hetero_W <- sem(model = model_MLMGSEM_hetero_W,
                        data = PISA_short_balanced_imp, 
                        cluster="j", 
                        group="CNT"
                        )
summary(fit_MLMGSEM_hetero_W)



### Model with Heterogeneous Between-Cluster (Co)variances

model_WFmultigroup_hetero_B <- paste(model_WF_W_homo, model_WF_B, sep="; ")

fit_WFmultigroup_hetero_B <- sem(model = model_WFmultigroup_hetero_B,
                                 data = PISA_short_balanced_imp_WF,
                                 group="CNT"#,
                                 #group.equal = c("lv.variances", "lv.covariances")
                                 ) 
summary(fit_WFmultigroup_hetero_B)


## Here you have to use the most recent version on Github (0.6-19.2186) again with its "genuine" ML MG SEM which yields very similar estimates.

model_MLMGSEM_hetero_B <- c(
  "
  Group: 1
  Level: 1
  CREATAS ~~ CREATAS_both*CREATAS
  GROSAGR ~~ GROSAGR_both*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_both*GROSAGR
  Level: 2
  CREATAS ~~ CREATAS_albania*CREATAS
  GROSAGR ~~ GROSAGR_albania*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_albania*GROSAGR
  
  Group: 2
  Level: 1
  CREATAS ~~ CREATAS_both*CREATAS
  GROSAGR ~~ GROSAGR_both*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_both*GROSAGR
  Level: 2
  CREATAS ~~ CREATAS_ireland*CREATAS
  GROSAGR ~~ GROSAGR_ireland*GROSAGR
  CREATAS ~~ CREATAS_GROSAGR_ireland*GROSAGR
  "
)

fit_MLMGSEM_hetero_B <- sem(model = model_MLMGSEM_hetero_B,
                            data = PISA_short_balanced_imp, 
                            cluster="j", 
                            group="CNT"
)
summary(fit_MLMGSEM_hetero_B)



### Model with Heterogeneous Within and Between-Cluster (Co)variances

model_WFmultigroup_hetero_WB <- paste(model_WF_W_hetero, model_WF_B, sep="; ")

fit_WFmultigroup_hetero_WB <- sem(model = model_WFmultigroup_hetero_WB,
                                 data = PISA_short_balanced_imp_WF,
                                 group="CNT"#,
                                 #group.equal = c("lv.variances", "lv.covariances")
) 
summary(fit_WFmultigroup_hetero_WB)

## the "genuine" lavaan ML MG SEM returns very similar estimates

fit_MLMGSEM_hetero_WB <- sem(model = model_MLMGSEM,
                        data = PISA_short_balanced_imp, 
                        cluster="j", 
                        group="CNT"#,
                        #group.equal = c("residuals", "residual.covariances") 
) 
summary(fit_MLMGSEM_hetero_WB)



### Model Comparisons

anova(fit_WFmultigroup_homo, fit_WFmultigroup_hetero_W)
anova(fit_WFmultigroup_homo, fit_WFmultigroup_hetero_B)
anova(fit_WFmultigroup_homo, fit_WFmultigroup_hetero_WB)
anova(fit_WFmultigroup_hetero_B, fit_WFmultigroup_hetero_WB)
# the most complex model, that has heterogeneous within- and between-cluster (co)variances, fits the data best