###################################################################
# R codes for making source figures
# 07: making figures
###################################################################
rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(R.matlab)

source("../../../00_utils.R")

library(tidyverse)
library(ggplot2)
library(reshape2)
library(colorspace)
library(ggpubr)

my_cleaner_theme = theme(
  legend.position  = "none",        # place the legend under the plot
  panel.border     = element_blank(), # remove the border
  panel.background = element_blank(), # remove the background
  panel.grid.major = element_blank(), # remove part of grid
  panel.grid.minor = element_blank(), # remove other part of grid
  axis.text.x      = element_blank(), # remove "1 2 3 ..." from x-axis
  axis.ticks.x     = element_blank(), # remove x-axis tick marks
  axis.title.x     = element_blank(), # remove "Var1"
  axis.text.y      = element_blank(), # remove "1 2 3 ..." from y-axis
  axis.ticks.y     = element_blank(), # remove y-axis tick marks
  axis.title.y     = element_blank()  # remove "Var2"
)


##########################################################################################
# simulate data from noisyICA setting
set.seed(42)
simData = SimFMRI123(var.inactive = 0,noisyICA = TRUE, snr=0.4,nTR = 50,phi = 0.47) 
xmat = simData$X
smat = simData$S
mmat = simData$Ms

p1=matrix(smat[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p2=matrix(smat[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p3=matrix(smat[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

################################################################################
# PCA+SparseICA
a=BIC_sparseICA(xData = xmat, n.comp = 3, whiten = "eigenvec", method = "C",
                eps = 1e-6,maxit = 1000, BIC_plot = F)
best_nu=a$best_nu

my_sparseICA = sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", restarts = 40, 
                         nu = best_nu,eps = 1e-6,maxit = 1000)
a=matchICA(my_sparseICA$estS,smat)

p4=matrix(a[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p5=matrix(a[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p6=matrix(a[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)


################################################################################
my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = 1e-6,maxit = 1000,restarts = 40)
b=matchICA(my_infomaxICA$S,smat)

p7=matrix(b[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p8=matrix(b[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p9=matrix(b[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

##########################################################################################
## fast ICA
my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = 1000,tol=1e-6,restarts = 40)
c=matchICA(my_fastICA$S,smat)

p10=matrix(c[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p11=matrix(c[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p12=matrix(c[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

# SICA-EBM
dat = readMat("test_EBM_123.mat")
temp = cbind(smat,matrix(rnorm(47*1089),nrow = 1089))
my_S_EBM = matchICA(t(dat$myS),temp)
d = my_S_EBM[,1:3]

p13=matrix(d[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p14=matrix(d[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p15=matrix(d[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

# sparsefastica
dat = readMat("test_sparsefastica_123.mat")
e = matchICA(t(dat$myS),smat)

p16=matrix(e[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p17=matrix(e[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

p18=matrix(e[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

ggarrange(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,p10,p11,p12,
          p13,p14,p15,p16,p17,p18,
          labels = c("A","","","B","","",
                     "C","","","D","","",
                     "E","","","F","",""),
          ncol = 6, nrow = 3)
ggsave("../../Figures/sim123_sources.png",width = 18, height = 9)
