### Histograms of raw data and cluster means help to decide whether within- and between-cluster variances might be heterogeneous.
# Here we generated (full) homogeneous, heterogeneous level-1, heterogeneous level-2, and (full) heterogeneous level-1 and level-2 data
# to check how these unfold in the plotted raw data and cluster means.
# Fig. A1 in Appendix

# load packages
library(lavaan)
library(ggplot2)
library(patchwork)

## set population characteristics

p <- 2 # nb of variables
VR <- 2 # variance ratio
k <- 2 # nb of groups
var_B_1 <- 0.10
cor_B <- 0.3
cor_W <- 0.3

# the rest is determined by these variables

var_B_2 <- var_B_1 / VR
var_W_1 <- 1 - var_B_1
var_W_2 <- var_W_1 / VR
cov_B_1 <- var_B_1 * cor_B
cov_B_2 <- var_B_2 * cor_B
cov_W_1 <- var_W_1 * cor_W
cov_W_2 <- var_W_2 * cor_W


## set sample characteristics

g <- 1000 # nb of clusters
n <- 30 # cluster size
# both have to be divisible by k=2


## variances

vars_w_1 <- c()
vars_w_2 <- c()
vars_b_1 <- c()
vars_b_2 <- c()
for (i in 1:p) {
  vars_w_1[i] <- paste0("x", i, "~~", var_W_1, "*", "x", i)
  vars_w_2[i] <- paste0("x", i, "~~", var_W_2, "*", "x", i)
  vars_b_1[i] <- paste0("x", i, "~~", var_B_1, "*", "x", i)
  vars_b_2[i] <- paste0("x", i, "~~", var_B_2, "*", "x", i)
}
vars_w_1 <- paste(vars_w_1, collapse = "; ")
vars_w_2 <- paste(vars_w_2, collapse = "; ")
vars_b_1 <- paste(vars_b_1, collapse = "; ")
vars_b_2 <- paste(vars_b_2, collapse = "; ")


## covariances

covs_w_1 <- c()
covs_w_2 <- c()
covs_b_1 <- c()
covs_b_2 <- c()
count = 0
for (i in 1:p) {
  for (j in 1:p) {
    if (i != j & j > i) {
      count <- count + 1
      covs_w_1[count] <- paste("x", i, "~~", cov_W_1, "*", "x", j, sep = "")
      covs_w_2[count] <- paste("x", i, "~~", cov_W_2, "*", "x", j, sep = "")
      covs_b_1[count] <- paste("x", i, "~~", cov_B_1, "*", "x", j, sep = "")
      covs_b_2[count] <- paste("x", i, "~~", cov_B_2, "*", "x", j, sep = "")
    }
  }
}
covs_w_1 <- paste(covs_w_1, collapse = "; ")
covs_w_2 <- paste(covs_w_2, collapse = "; ")
covs_b_1 <- paste(covs_b_1, collapse = "; ")
covs_b_2 <- paste(covs_b_2, collapse = "; ")

# put variances and covariances together (means are 0 per default)
popModel_W_1 <- paste(vars_w_1, covs_w_1, sep = ";")
popModel_W_2 <- paste(vars_w_2, covs_w_2, sep = ";")
popModel_B_1 <- paste(vars_b_1, covs_b_1, sep = ";")
popModel_B_2 <- paste(vars_b_2, covs_b_2, sep = ";")


### generate data
set.seed(1678)

## (I) homo

# sample data
sample_B_1 <- simulateData(popModel_B_1, sample.nobs = g, model.type = "lavaan")
sample_W_1 <- simulateData(popModel_W_1,
                           sample.nobs = n * g,
                           model.type = "lavaan")

# merge data
count <- 0
clusters <- rep(1:g, each = n)
data_homo <- as.data.frame(matrix(NA, ncol = p, nrow = g * n))
for (i in 1:(n * g)) {
  count <- count + 1
  j <- clusters[count]
  data_homo[count, ] <- sample_W_1[i, ] + sample_B_1[j, ]
}


## (II) hetero-1

sample_B_1 <- simulateData(popModel_B_1, sample.nobs = g, model.type = "lavaan")
sample_W_1 <- simulateData(popModel_W_1,
                           sample.nobs = n * g / k,
                           model.type = "lavaan")
sample_W_2 <- simulateData(popModel_W_2,
                           sample.nobs = n * g / k,
                           model.type = "lavaan")

# merge data
count <- 0
clusters <- rep(1:g, each = n)
data_hetero_1 <- as.data.frame(matrix(NA, ncol = p, nrow = g * n))
for (l in 1:k) {
  sample_W <- get(paste0("sample_W_", l))
  for (i in 1:(n * g / k)) {
    count <- count + 1
    j <- clusters[count]
    data_hetero_1[count, ] <- sample_W[i, ] + sample_B_1[l, ]
  }
}


## (III) hetero-2

sample_B_1 <- simulateData(popModel_B_1,
                           sample.nobs = g / k,
                           model.type = "lavaan")
sample_B_2 <- simulateData(popModel_B_1,
                           sample.nobs = g / k,
                           model.type = "lavaan")
sample_W_1 <- simulateData(popModel_W_1,
                           sample.nobs = n * g,
                           model.type = "lavaan")

# merge data
count <- 0
clusters <- rep(1:g, each = n)
data_hetero_2 <- as.data.frame(matrix(NA, ncol = p, nrow = g * n))
for (l in 1:k) {
  # merge the sampled data from both levels
  sample_B <- get(paste0("sample_B_", l))
  for (i in 1:(n * g / k)) {
    count <- count + 1
    j <- clusters[count] - g/k * (l-1)
    data_hetero_2[count, ] <- sample_W_1[count, ] + sample_B[j, ]
  }
}


## (IV) hetero-12

sample_B_1 <- simulateData(popModel_B_1,
                           sample.nobs = g / k,
                           model.type = "lavaan")
sample_B_2 <- simulateData(popModel_B_1,
                           sample.nobs = g / k,
                           model.type = "lavaan")
sample_W_1 <- simulateData(popModel_W_1,
                           sample.nobs = n * g / k,
                           model.type = "lavaan")
sample_W_2 <- simulateData(popModel_W_2,
                           sample.nobs = n * g / k,
                           model.type = "lavaan")

# merge data
count <- 0
clusters <- rep(1:g, each = n)
data_hetero_12 <- as.data.frame(matrix(NA, ncol = p, nrow = g * n))
for (l in 1:k) {
  # merge the sampled data from both levels
  sample_B <- get(paste0("sample_B_", l))
  sample_W <- get(paste0("sample_W_", l))
  for (i in 1:(n * g / k)) {
    count <- count + 1
    j <- clusters[count]
    data_hetero_12[count, ] <- sample_W[i, ] + sample_B[l, ]
  }
}

## create combined data frame

data <- data.frame(
  j = as.factor(clusters),
  i = rep(1:n, g),
  k = as.factor(rep(1:k, each = n * g / k)),
  x1_homo = data_homo$V1,
  x2_homo = data_homo$V2,
  x1_hetero_1 = data_hetero_1$V1,
  x2_hetero_1 = data_hetero_1$V2,
  x1_hetero_2 = data_hetero_2$V1,
  x2_hetero_2 = data_hetero_2$V2,
  x1_hetero_12 = data_hetero_12$V1,
  x2_hetero_12 = data_hetero_12$V2
)

## create data frame with cluster means

cluster_means <- aggregate(list(data$x1_homo, data$x2_homo, 
                                data$x1_hetero_1, data$x2_hetero_1, 
                                data$x1_hetero_2, data$x2_hetero_2,
                                data$x1_hetero_12, data$x2_hetero_12),
                           list(data$j), FUN=mean)

colnames(cluster_means) <- c("j", "x1_homo", "x2_homo", "x1_hetero_1", "x2_hetero_1", "x1_hetero_2", "x2_hetero_2", "x1_hetero_12", "x2_hetero_12") 
cluster_means$k <- as.factor(rep(1:k, each=(g/k)))


### plot raw data (-1)

a1 <- ggplot(data, aes(x=x1_homo, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-5, 5)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 3000)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))  + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position="none") + 
  labs(title="Homogeneous")

b1 <- ggplot(data, aes(x=x1_hetero_1, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-5, 5)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 3000)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position="none") + 
  labs(title="Heterogeneous Level-1")

c1 <- ggplot(data, aes(x=x1_hetero_2, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-5, 5)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 3000)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position="none") + 
  labs(title="Heterogeneous Level-2")

d1 <- ggplot(data, aes(x=x1_hetero_12, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-5, 5)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 3000)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position="none") + 
  labs(title="Heterogeneous Level-1 and Level-2")



### plot cluster means (-2)

a2 <- ggplot(cluster_means, aes(x=x1_homo, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-2, 2)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 200)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5))  + 
  theme(legend.position="none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

b2 <- ggplot(cluster_means, aes(x=x1_hetero_1, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-2, 2)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 200)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  theme(legend.position="none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

c2 <- ggplot(cluster_means, aes(x=x1_hetero_2, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-2, 2)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 200)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  theme(legend.position="none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

d2 <- ggplot(cluster_means, aes(x=x1_hetero_12, fill=k)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_x_continuous(name="", expand=c(0.05, 0.05), limits=c(-2, 2)) +
  scale_y_continuous(name="Frequency", expand=c(0, 0), limits=c(0, 200)) +
  scale_fill_manual(name="Group", values=c("#002654", "#ffce00")) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) + 
  theme(legend.position="none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

a1 + b1 + c1 + d1 + a2 + b2 + c2 + d2 + plot_layout(nrow=2, guides='collect') & theme(text = element_text("serif"), legend.position = "bottom") 

ggsave(
  filename = "FigA1.jpeg",
  device = "jpeg",
  plot = last_plot(),
  width = 1400*3, 
  height = 500*3,
  units = "px",
  dpi = 300
)
