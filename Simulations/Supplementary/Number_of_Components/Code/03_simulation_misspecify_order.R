###################################################################
# R codes for comparing methods for selecting the number of ICs
# 03: conduct sim123 simulations
###################################################################
rm(list = ls())

library(steadyICA)
library(fastICA)
library(irlba)
library(SparseICA)

source("../../../00_utils.R")

#################################################################################################
# For three given SNR, compare PMSE of S and M
my_SNR = c(0.4,1.5,3)
# Number of replications for each variable type
rep = 100
# Level for variance of inactive points in components
var_level = 0
# set a seed
seed1 = 412036
seed2 = 365298
seed3 = 753169
# number of time points in each simulation
my_nTR = 50
# AR structure in each simulation
my_phi = 0.47
# convergence matric
my_eps = 1e-6
# maximum number of iterations
my_maxit = 1000
# whether it's noisy ICA model
my_noisyICA = TRUE 


#######################################################################################################
# Low SNR
S_PMSE1 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep),
                     fast2=numeric(rep),fast3=numeric(rep),fast4=numeric(rep))
M_PMSE1 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep),
                     fast2=numeric(rep),fast3=numeric(rep),fast4=numeric(rep))

set.seed(seed1)
for (i in 1:rep) {
  # set.seed(rep + (seed-1)*i)
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[1],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms

  # standardized xmat and PCA
  xmat_center = scale(xmat,center = T,scale = F)

  ###############################################################################
  # sparse ICA single start Q=3
  select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", method = "C",
                                  eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE1$sparse3[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE1$sparse3[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=3
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE1$fast3[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$fast3[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  
  ###############################################################################
  # sparse ica single start Q=2
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", method = "C",
                                 eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu

  my_sparseICA = sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)

  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE1$sparse2[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE1$sparse2[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=2
  my_fastICA = fastICA(xmat,n.comp = 2,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE1$fast2[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$fast2[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  ###############################################################################
  # sparse ica single start Q=4
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", method = "C",
                                  eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE1$sparse4[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE1$sparse4[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=4
  my_fastICA = fastICA(xmat,n.comp = 4,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE1$fast4[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$fast4[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  cat(i," finished!\n")
}


#######################################################################################################
# Medium SNR
S_PMSE2 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep))
M_PMSE2 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep))

set.seed(seed2)
for (i in 1:rep) {
  # set.seed(rep + (seed-1)*i)
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[2],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # standardized xmat and PCA
  xmat_center = scale(xmat,center = T,scale = F)
  
  ###############################################################################
  # sparse ICA single start Q=3
  select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", method = "C",
                                  eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE2$sparse3[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE2$sparse3[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=3
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE2$fast3[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$fast3[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  
  ###############################################################################
  # sparse ica single start Q=2
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", method = "C",
                                   eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE2$sparse2[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE2$sparse2[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=2
  my_fastICA = fastICA(xmat,n.comp = 2,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE2$fast2[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$fast2[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  ###############################################################################
  # sparse ica single start Q=4
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", method = "C",
                                   eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE2$sparse4[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE2$sparse4[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=4
  my_fastICA = fastICA(xmat,n.comp = 4,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE2$fast4[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$fast4[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  cat(i," finished!\n")
}


#######################################################################################################
# high SNR
S_PMSE3 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep))
M_PMSE3 = data.frame(sparse2=numeric(rep),sparse3=numeric(rep),sparse4=numeric(rep))

set.seed(seed3)
for (i in 1:rep) {
  # set.seed(rep + (seed-1)*i)
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[3],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # standardized xmat and PCA
  xmat_center = scale(xmat,center = T,scale = F)
  
  ###############################################################################
  # sparse ICA single start Q=3
  select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", method = "C",
                                  eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE3$sparse3[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE3$sparse3[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=3
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE3$fast3[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$fast3[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  
  ###############################################################################
  # sparse ica single start Q=2
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", method = "C",
                                   eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 2,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE3$sparse2[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE3$sparse2[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=2
  my_fastICA = fastICA(xmat,n.comp = 2,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE3$fast2[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$fast2[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  ###############################################################################
  # sparse ica single start Q=4
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", method = "C",
                                   eps = 1e-6,maxit = 1000, BIC_plot = F,verbose=F,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 4,whiten = "eigenvec", restarts = 1,
                           nu = my_nu,eps = 1e-6, maxit = 1000)
  
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE3$sparse4[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE3$sparse4[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  
  # fast ICA single start Q=4
  my_fastICA = fastICA(xmat,n.comp = 4,maxit = 1000,tol=1e-6,method = "C")
  my_M_fast = est.M.ols(my_fastICA$S,xmat_center)
  S_PMSE3$fast4[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$fast4[i] = sqrt(frobICA(M1=my_M_fast,M2=mmat,standardize = TRUE))
  
  cat(i," finished!\n")
}

save(S_PMSE1,S_PMSE2,S_PMSE3,M_PMSE1,M_PMSE2,M_PMSE3,file = "03_PRMSE_misQ.RData")
