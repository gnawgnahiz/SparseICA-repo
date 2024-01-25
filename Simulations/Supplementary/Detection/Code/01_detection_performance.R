###################################################################
# R codes for evaluating the detection performance of Sparse ICA
###################################################################
rm(list = ls())

library(steadyICA)
library(fastICA)
library(SparseICA)
library(mltools)
library(MLmetrics)

source("00_utils.R")

##############################################################################################
# sim123 low SNR
# repeat 100 times
rep = 100
measure_lowSNR = data.frame(AUC = numeric(rep), Matthew=numeric(rep), F1score=numeric(rep))
set.seed(285323)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=0.4, nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  nu_selection = BIC_sparseICA(xmat,3)
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,col.stand = T,row.stand = F,iter.stand=0,
                           U.list=NULL,whiten = "eigenvec", restarts = 1,lambda = sqrt(2)/2,
                           nu = nu_selection$best_nu,eps = 1e-6, maxit = 1000, converge_plot = F,verbose = T)
  
  S_true01 = as.numeric(smat)
  S_true01[S_true01!=0]=1
  
  S_sparse = matchICA(my_sparseICA$estS,smat)
  S_est01 = as.numeric(S_sparse)
  S_est01[S_est01!=0]=1
  
  measure_lowSNR$AUC[i] = AUC(S_est01,S_true01)
  measure_lowSNR$Matthew[i] = mcc(S_est01,S_true01)
  measure_lowSNR$F1score[i] = F1_Score(S_true01,S_est01)
}


mean(measure_lowSNR$AUC)
sd(measure_lowSNR$AUC)
mean(measure_lowSNR$Matthew)
sd(measure_lowSNR$Matthew)
mean(measure_lowSNR$F1score)
sd(measure_lowSNR$F1score)
boxplot(measure_lowSNR$AUC,measure_lowSNR$Matthew,measure_lowSNR$F1score)


##############################################################################################
# sim123 medium SNR
# repeat 100 times
rep = 100
measure_mediumSNR = data.frame(AUC = numeric(rep), Matthew=numeric(rep), F1score=numeric(rep))
set.seed(75468)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=1.5, nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  nu_selection = BIC_sparseICA(xmat,3)
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,col.stand = T,row.stand = F,iter.stand=0,
                           U.list=NULL,whiten = "eigenvec", restarts = 1,lambda = sqrt(2)/2,
                           nu = nu_selection$best_nu,eps = 1e-6, maxit = 1000, converge_plot = F,verbose = T)
  
  S_true01 = as.numeric(smat)
  S_true01[S_true01!=0]=1
  
  S_sparse = matchICA(my_sparseICA$estS,smat)
  S_est01 = as.numeric(S_sparse)
  S_est01[S_est01!=0]=1
  
  measure_mediumSNR$AUC[i] = AUC(S_est01,S_true01)
  measure_mediumSNR$Matthew[i] = mcc(S_est01,S_true01)
  measure_mediumSNR$F1score[i] = F1_Score(S_true01,S_est01)
}


mean(measure_mediumSNR$AUC)
sd(measure_mediumSNR$AUC)
mean(measure_mediumSNR$Matthew)
sd(measure_mediumSNR$Matthew)
mean(measure_mediumSNR$F1score)
sd(measure_mediumSNR$F1score)
boxplot(measure_mediumSNR$AUC,measure_mediumSNR$Matthew,measure_mediumSNR$F1score)


##############################################################################################
# sim123 high SNR
# repeat 100 times
rep = 100
measure_highSNR = data.frame(AUC = numeric(rep), Matthew=numeric(rep), F1score=numeric(rep))
set.seed(197216)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=3, nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  nu_selection = BIC_sparseICA(xmat,3)
  
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,col.stand = T,row.stand = F,iter.stand=0,
                           U.list=NULL,whiten = "eigenvec", restarts = 1,lambda = sqrt(2)/2,
                           nu = nu_selection$best_nu,eps = 1e-6, maxit = 1000, converge_plot = F,verbose = T)
  
  S_true01 = as.numeric(smat)
  S_true01[S_true01!=0]=1
  
  S_sparse = matchICA(my_sparseICA$estS,smat)
  S_est01 = as.numeric(S_sparse)
  S_est01[S_est01!=0]=1
  
  measure_highSNR$AUC[i] = AUC(S_est01,S_true01)
  measure_highSNR$Matthew[i] = mcc(S_est01,S_true01)
  measure_highSNR$F1score[i] = F1_Score(S_true01,S_est01)
}


mean(measure_highSNR$AUC)
sd(measure_highSNR$AUC)
mean(measure_highSNR$Matthew)
sd(measure_highSNR$Matthew)
mean(measure_highSNR$F1score)
sd(measure_highSNR$F1score)
boxplot(measure_highSNR$AUC,measure_highSNR$Matthew,measure_highSNR$F1score)


##############################################################################################
# in realistic simulation
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench') 
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 

load("../Data/group_sparseICA.RData")
load("../Data/template_xifti.RData")
load("../Data/template_timecourses2.RData")


##############################################################################################
# realistic low SNR
# repeat 100 times
rep = 100
nu_list = seq(0.1,3,0.1)
measure_highdim = data.frame(AUC = numeric(rep), Matthew=numeric(rep), F1score=numeric(rep))
set.seed(197216)
for (i in 1:rep) {
  simData = SimFMRIreal(snr=0.4, noisyICA = T, nTR=120, phi=0.47, 
                        var.inactive = 0, components = c(17,6,23), 
                        true_data = dat_iter5, template_xifti = template_image,
                        template_time = sample_time)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", method = "C",
                                  lambda = sqrt(2)/2, eps = 1e-6,maxit = 500, 
                                  BIC_plot = F,nu_list = nu_list)
  my_nu = select_sparseICA$best_nu
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,col.stand = T,row.stand = F,iter.stand=0,method = "C",
                           U.list=NULL,whiten = "eigenvec", restarts = 1,lambda = sqrt(2)/2, 
                           nu = my_nu,eps = 1e-6,maxit = 500, converge_plot = F,verbose = F)
  
  S_true01 = as.numeric(smat)
  S_true01[S_true01!=0]=1
  
  S_sparse = matchICA(my_sparseICA$estS,smat)
  S_est01 = as.numeric(S_sparse)
  S_est01[S_est01!=0]=1
  
  measure_highdim$AUC[i] = AUC(S_est01,S_true01)
  measure_highdim$Matthew[i] = mcc(S_est01,S_true01)
  measure_highdim$F1score[i] = F1_Score(S_true01,S_est01)
  
  cat("#####################The ",i," th run is finished !###########################\n")
  print(measure_highdim[i,])
}

mean(measure_highdim$AUC)
sd(measure_highdim$AUC)
mean(measure_highdim$Matthew)
sd(measure_highdim$Matthew)
mean(measure_highdim$F1score)
sd(measure_highdim$F1score)
boxplot(measure_highdim$AUC,measure_highdim$Matthew,measure_highdim$F1score)

summary(measure_highdim$Matthew)

save(measure_lowSNR,measure_mediumSNR,measure_highSNR,measure_highdim,file = "../Results/AUC_MCC_F1_all.RData")
