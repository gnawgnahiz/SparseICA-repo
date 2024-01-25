###################################################################
# R codes for comparing methods for selecting the number of ICs
# 02: make figures for estimation under mis-specification
###################################################################
rm(list = ls())

library(fastICA)
library(SparseICA)
library(ggplot2)
library(reshape2)
library(colorspace)
library(ggpubr)

source("../../../00_utils.R")

###################################################################
# make plots for mis-specifying the number of ICs

# set themes
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


set.seed(42)
simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=0.4,nTR = 50,phi = 0.47)
xmat = simData$X
smat = simData$S
mmat = simData$Ms

###################################################################
# plots for Sparse ICA
###################################################################
# truth
psmat1 = matrix(smat[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
psmat2 = matrix(smat[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
psmat3 = matrix(smat[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)


# number of components correctly specified
select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", method = "C",
                                eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
my_nu = select_sparseICA$best_nu

my_sparseICA1 = sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", restarts = 40,
                          nu = my_nu,eps = 1e-6, maxit = 1000)

a=matchICA(my_sparseICA1$estS,smat)

pQ3_1 = matrix(a[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ3_2 = matrix(a[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ3_3 = matrix(a[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)


# number of components under-specified: 2
select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", method = "C",
                                eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
my_nu = select_sparseICA$best_nu

my_sparseICA2 = sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", restarts = 40,
                          nu = my_nu,eps = 1e-6, maxit = 1000)

b=matchICA(cbind(my_sparseICA2$estS,rnorm(1089,0,0.01)),smat)

pQ2_1 = matrix(b[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ2_2 = matrix(b[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ2_3 = matrix(b[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

# number of components over-specified: 4
select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", method = "C",
                                eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
my_nu = select_sparseICA$best_nu

my_sparseICA3 = sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", restarts = 40,
                          nu = my_nu,eps = 1e-6, maxit = 1000)
smat3 = cbind(smat,rnorm(1089,0,0.01))
c=matchICA(my_sparseICA3$estS,smat3)

pQ4_1 = matrix(c[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE) +
  xlab("IC1")
pQ4_2 = matrix(c[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ4_3 = matrix(c[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ4_4 = matrix(c[,4],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

pempty = ggplot() + my_cleaner_theme

ggarrange(psmat1,psmat2,psmat3,pempty,pQ2_2,pQ2_3,pempty,
          pQ3_1,pQ3_2,pQ3_3,pQ4_1,pQ4_2,pQ4_3,pQ4_4,
          labels = c("Truth","","","Q=2","","","",
                     "Q=3","","","Q=4","","",""),
          ncol = 7, nrow = 2,label.y = 1.013)
ggsave("Supp_Figure12.png")


###################################################################
# plots for Fast ICA
###################################################################
# truth
psmat1 = matrix(smat[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
psmat2 = matrix(smat[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
psmat3 = matrix(smat[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)


# number of components correctly specified
my_fastICA1 = fastICA(xmat,n.comp = 3,maxit = 1000,tol=1e-6,method = "C")

a=matchICA(my_fastICA1$S,smat)

pQ3_1 = matrix(a[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ3_2 = matrix(a[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ3_3 = matrix(a[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)


# number of components under-specified: 2
my_fastICA2 = fastICA(xmat,n.comp = 2,maxit = 1000,tol=1e-6,method = "C")

b=matchICA(cbind(my_fastICA2$S,rnorm(1089,0,0.01)),smat)

pQ2_1 = matrix(b[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ2_2 = matrix(b[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ2_3 = matrix(b[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

# number of components over-specified: 4
my_fastICA3 = fastICA(xmat,n.comp = 4,maxit = 1000,tol=1e-6,method = "C")
smat3 = cbind(smat,rnorm(1089,0,0.01))
c=matchICA(my_fastICA3$S,smat3)

pQ4_1 = matrix(c[,1],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE) +
  xlab("IC1")
pQ4_2 = matrix(c[,2],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ4_3 = matrix(c[,3],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)
pQ4_4 = matrix(c[,4],33) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  my_cleaner_theme +
  scale_fill_continuous_sequential(palette = "YlOrRd",rev=FALSE)

pempty = ggplot() + my_cleaner_theme

ggarrange(psmat1,psmat2,psmat3,pempty,pQ2_2,pQ2_3,pempty,
          pQ3_1,pQ3_2,pQ3_3,pQ4_1,pQ4_2,pQ4_3,pQ4_4,
          labels = c("Truth","","","Q=2","","","",
                     "Q=3","","","Q=4","","",""),
          ncol = 7, nrow = 2,label.y = 1.013)
ggsave("Supp_Figure13.png")

