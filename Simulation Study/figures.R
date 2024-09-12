## Code to create Fig. 2

# load packages
library(ggplot2)
library(dplyr) 
library(patchwork)
library(cowplot) # get_legend

# set paths and load data
pathO <- "objects/"
pathF <- "figures/"
dat <- readRDS(paste0(pathO, "dat.rds"))

# data preparation

dat$n_fac <- factor(dat$n, ordered=TRUE, levels=unique(dat$n))
dat$VR_fac <- factor(dat$VR, levels = unique(dat$VR), labels = paste0("VR==", unique(dat$VR)))
dat$var_B_fac <- factor(dat$var_B, levels = unique(dat$var_B), labels = paste0("sigma[B]^2==", unique(dat$var_B)))
dat$meanICC_fac <- apply(X=dat, MARGIN=1, FUN=function(x){paste0("bar(rho)==", round((dat$var_B/(dat$var_B+(dat$var_W_1+dat$var_W_2)/2) ), 2))})[,1]
dat$meanTheta_B_fac <- apply(X=dat, MARGIN=1, FUN=function(x){paste0("bar(Theta)[B]==", round((2*dat$var_B+dat$var_B*0.3)/3, 2))})[,1]
dat$meanTheta_W_fac <- apply(X=dat, MARGIN=1, FUN=function(x){paste0("bar(Theta)[W]==", round(((2*dat$var_W_1+dat$var_W_1*0.3)/3 + (2*dat$var_W_2+dat$var_W_2*0.3)/3)/2, 2))})[,1]

expdD <- 0.08 # expand down
expdU <- 0.75 # expand up when strong effects with larger g (x-axis) (larger because of label with parameter values)
dist <- 40 # distance between all three values of n at a given g


### Fig. 2: Estimation Accuracy of Between-Cluster, Within-Cluster, and ICC Parameter Estimates 

# between
a <- ggplot(dat, aes(x=g, y=relRMSE_B, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanTheta_B_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative RMSE (%)", expand=expand_scale(mult = c(expdD, expdD)), limits=c(0, NA)) +
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  scale_shape_manual(name="Approach", values=c(17, 16), guide = "none") +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") ))  +
  labs(title="Between-Cluster")

b <- ggplot(dat, aes(x=g, y=relBias_B, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  geom_hline(yintercept=0, size=0.5, col="darkgrey", lty="longdash") +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanTheta_B_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative Bias (%)", expand=expand_scale(mult = c(expdD, expdU)), limits=c(NA, NA)) +
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  scale_shape_manual(name="Approach", values=c(17, 16), guide = "none") +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5))  +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") ))

# within
c <- ggplot(dat, aes(x=g, y=relRMSE_W, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanTheta_W_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative RMSE (%)", expand=expand_scale(mult = c(expdD, expdU)), limits=c(0, NA)) +
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  scale_shape_manual(name="Approach", values=c(17, 16), guide = "none") +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") )) +
  labs(title="Within-Cluster")

d <- ggplot(dat, aes(x=g, y=relBias_W, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  geom_hline(yintercept=0, size=0.5, col="darkgrey", lty="longdash") +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanTheta_W_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative Bias (%)", expand=expand_scale(mult = c(expdD, expdD)), limits=c(NA, NA)) +
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  scale_shape_manual(name="Approach", values=c(17, 16)) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5))  +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") ))

# ICC
e <- ggplot(dat, aes(x=g, y=relRMSE_ICC, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanICC_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative RMSE (%)", expand=expand_scale(mult = c(expdD, expdU)), limits=c(0, NA)) +
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  scale_shape_manual(name="Approach", values=c(17, 16)) +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") ))  +
  labs(title="ICC")

f <- ggplot(dat, aes(x=g, y=relBias_ICC, col=n_fac)) + 
  facet_grid(var_B_fac ~ VR_fac, labeller=label_parsed) +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  geom_hline(yintercept=0, size=0.5, col="darkgrey", lty="longdash") +
  stat_summary(data=dat, fun="mean", alpha=0.6) +
  geom_label(aes(x=Inf, y=Inf, label=meanICC_fac, family="serif"), col="black", alpha=0.05, parse=TRUE, hjust=1, vjust=1) +
  scale_y_continuous(name="Relative Bias (%)", expand=expand_scale(mult = c(expdD, expdU))) + 
  scale_x_continuous(name="Number of Clusters (g)", breaks=unique(dat$g)) +
  scale_color_discrete(name="Cluster Size (n)") +
  theme_minimal() +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, margin=margin(10,0,0,0)),
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5))  +
  theme(legend.position="none") + 
  guides(col = guide_legend(override.aes = list(linetype = "blank") ))


a + b + c + d + e + f + plot_layout(nrow=3, ncol=2, guides='collect') & theme(text = element_text("serif"), legend.position = "bottom") 


ggsave(
  filename = "fig2.jpeg",
  path = pathF,
  device = "jpeg",
  plot = last_plot(),
  width = 800*3, # values used for png multiplied with 3 (bc for png, dpi=97)
  height = 1000*3,
  units = "px",
  dpi = 300
)

