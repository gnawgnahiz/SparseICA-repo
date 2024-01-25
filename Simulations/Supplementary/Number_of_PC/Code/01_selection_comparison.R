###################################################################
# R codes for comparing methods for selecting the number of ICs
###################################################################
rm(list = ls())

# library(steadyICA)
# library(fastICA)
# library(irlba)
# #library(R.matlab)
# library(SparseICA)
library(pesel)
#library(fMRItools)

source("00_utils.R")
source("minka.R")

# formula for MDL method
mdl = function(N,K,M){
  num = prod(eigen_values[(N+1):K])^(1/(K-N))
  dem = sum(eigen_values[(N+1):K])/(K-N)
  likelihood = log(num/dem)
  res = -M*(K-N)*likelihood+log(M)*(1+N*K+(N-1)/2)/2
  return(res)
}

rep = 100
###################################################################
# single-subject sim123 simulation, low SNR

order_lowSNR = data.frame(minka=numeric(rep),
                          pesel=numeric(rep),
                          mdl=numeric(rep))
set.seed(25976)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=0.4,nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # Minka method
  eigen_values = eigen(as.matrix(cov(scale(xmat))))$val
  order_lowSNR$minka[i] = minka2001(lambda = eigen_values, M = 50, BIC=T,tau = 0.0001)
  
  # PESEL method
  my_pesel = pesel(X=xmat,npc.min = 0,npc.max = 50,
                 prior = NULL,scale = F,
                 method = "homogenous",
                 asymptotics = NULL)
  order_lowSNR$pesel[i] = my_pesel$nPCs
  
  # MDL method
  # M: number of voxels
  M = 1089
  # K: number of time points
  K = 50
  # N: number of PCs
  my_mdl = data.frame(N = 1:K,
                      MDL = rep(NA,K))
  for (j in 1:K) {
    my_mdl$MDL[j]=mdl(N=my_mdl$N[j],K,M)
  }
  order_lowSNR$mdl = my_mdl$N[which.min(my_mdl$MDL)]
}
median(order_lowSNR$minka)
median(order_lowSNR$pesel)
median(order_lowSNR$mdl)


###################################################################
# single-subject sim123 simulation, medium SNR

order_mediumSNR = data.frame(minka=numeric(rep),
                          pesel=numeric(rep),
                          mdl=numeric(rep))
set.seed(653248)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=1.5,nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # Minka method
  eigen_values = eigen(as.matrix(cov(scale(xmat))))$val
  order_mediumSNR$minka[i] = minka2001(lambda = eigen_values, M = 50, BIC=T,tau = 0.0001)
  
  # PESEL method
  my_pesel = pesel(X=xmat,npc.min = 0,npc.max = 50,
                   prior = NULL,scale = F,
                   method = "homogenous",
                   asymptotics = NULL)
  order_mediumSNR$pesel[i] = my_pesel$nPCs
  
  # MDL method
  # M: number of voxels
  M = 1089
  # K: number of time points
  K = 50
  # N: number of PCs
  my_mdl = data.frame(N = 1:K,
                      MDL = rep(NA,K))
  for (j in 1:K) {
    my_mdl$MDL[j]=mdl(N=my_mdl$N[j],K,M)
  }
  order_mediumSNR$mdl = my_mdl$N[which.min(my_mdl$MDL)]
}
median(order_mediumSNR$minka)
median(order_mediumSNR$pesel)
median(order_mediumSNR$mdl)


###################################################################
# single-subject sim123 simulation, high SNR

order_highSNR = data.frame(minka=numeric(rep),
                             pesel=numeric(rep),
                             mdl=numeric(rep))
set.seed(7856324)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=3,nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # Minka method
  eigen_values = eigen(as.matrix(cov(scale(xmat))))$val
  order_highSNR$minka[i] = minka2001(lambda = eigen_values, M = 50, BIC=T,tau = 0.0001)
  
  # PESEL method
  my_pesel = pesel(X=xmat,npc.min = 0,npc.max = 50,
                   prior = NULL,scale = F,
                   method = "homogenous",
                   asymptotics = NULL)
  order_highSNR$pesel[i] = my_pesel$nPCs
  
  # MDL method
  # M: number of voxels
  M = 1089
  # K: number of time points
  K = 50
  # N: number of PCs
  my_mdl = data.frame(N = 1:K,
                      MDL = rep(NA,K))
  for (j in 1:K) {
    my_mdl$MDL[j]=mdl(N=my_mdl$N[j],K,M)
  }
  order_highSNR$mdl = my_mdl$N[which.min(my_mdl$MDL)]
}
median(order_highSNR$minka)
median(order_highSNR$pesel)
median(order_highSNR$mdl)

###################################################################
# single-subject sim123 simulation, super high SNR
order_superhighSNR = data.frame(minka=numeric(rep),
                           pesel=numeric(rep),
                           mdl=numeric(rep))
set.seed(852146)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=10,nTR = 50,phi = 0.47)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # Minka method
  eigen_values = eigen(as.matrix(cov(scale(xmat))))$val
  order_superhighSNR$minka[i] = minka2001(lambda = eigen_values, M = 50, BIC=T,tau = 0.0001)
  
  # PESEL method
  my_pesel = pesel(X=xmat,npc.min = 0,npc.max = 50,
                   prior = NULL,scale = F,
                   method = "homogenous",
                   asymptotics = NULL)
  order_superhighSNR$pesel[i] = my_pesel$nPCs
  
  # MDL method
  # M: number of voxels
  M = 1089
  # K: number of time points
  K = 50
  # N: number of PCs
  my_mdl = data.frame(N = 1:K,
                      MDL = rep(NA,K))
  for (j in 1:K) {
    my_mdl$MDL[j]=mdl(N=my_mdl$N[j],K,M)
  }
  order_superhighSNR$mdl = my_mdl$N[which.min(my_mdl$MDL)]
}
median(order_superhighSNR$minka)
median(order_superhighSNR$pesel)
median(order_superhighSNR$mdl)

###################################################################
# high-dimensional simulation, low SNR
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench') 
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 

load("group_sparseICA.RData")
load("template_xifti.RData")
load("template_timecourses2.RData")

order_highdim = data.frame(minka=numeric(rep),
                           pesel=numeric(rep),
                           mdl=numeric(rep))
set.seed(658943)
for (i in 1:rep) {
  simData = SimFMRIreal(snr=0.4, noisyICA = T, nTR=120, phi=0.47, 
                        var.inactive = 0, components = c(17,6,23), 
                        true_data = dat_iter5, template_xifti = template_image,
                        template_time = sample_time)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  
  # Minka method
  eigen_values = eigen(as.matrix(cov(scale(xmat))))$val
  order_highdim$minka[i] = minka2001(lambda = eigen_values, M = 120, BIC=T,tau = 0.0001)
  
  # PESEL method
  my_pesel = pesel(X=xmat,npc.min = 0,npc.max = 120,
                   prior = NULL,scale = F,
                   method = "homogenous",
                   asymptotics = NULL)
  order_highdim$pesel[i] = my_pesel$nPCs
  
  # MDL method
  # M: number of voxels
  M = 59412
  # K: number of time points
  K = 120
  # N: number of PCs
  my_mdl = data.frame(N = 1:K,
                      MDL = rep(NA,K))
  for (j in 1:K) {
    my_mdl$MDL[j]=mdl(N=my_mdl$N[j],K,M)
  }
  order_highdim$mdl = my_mdl$N[which.min(my_mdl$MDL)]
  
  cat(i," finished!\n")
}
median(order_highdim$minka)
median(order_highdim$pesel)
median(order_highdim$mdl)

save(order_lowSNR,order_mediumSNR,order_highSNR,order_superhighSNR,order_highdim,file = "01_selection_comparison.RData")
