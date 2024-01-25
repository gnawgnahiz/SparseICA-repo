###################################################################
# R codes for simulation setting 2 (High-Dimensional) in single-subject simulations
# 04: generate simulated high-dim data, truth is nonsparse
###################################################################
rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(R.matlab)
library(ciftiTools)
source("00_utils.R")

# ciftiTools package requires users using the following codes to point to the path to workbench
# For Mac
ciftiTools.setOption('wb_path', '/Applications/workbench')
# For Windows
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
#ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench') 


###################################################################
# Create fastICA version true components for high-dimensional simulations
###################################################################
load("../../Data/PC30_subPC85.RData")
my_fastICA = fastICA_restarts(PC30,n.comp = 30,restarts = 40,method = "C",maxit = 500,tol = 1e-6,verbose = T)
a=signchange(t(my_fastICA$S))
my_fastICA$estS_sign=t(a$S)

cifti_image=read_cifti("../../Data/sub-29287_ses-1_task-rest_acq-rc8chan_run-1_space-fsLR_den-91k_bold.dtseries.nii",brainstructures = "all")
cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = 30)

view_cifti_surface(xifti = select_xifti(cifti_image,13),zlim = c(-2.43,2.82))
view_cifti_surface(xifti = select_xifti(cifti_image,23),zlim = c(-2.43,2.82))
view_cifti_surface(xifti = select_xifti(cifti_image,30),zlim = c(-2.43,2.82))

dat_fastICA = my_fastICA
save(dat_fastICA,file = "../../Data/group_fastICA.RData")

###################################################################
# Load true components data
###################################################################
load("../../Data/group_fastICA.RData")
load("../../Data/template_xifti.RData")
load("../../Data/template_timecourses2.RData")


###################################################################
# General simulation settings
###################################################################
# Simulate in low SNR level
my_SNR = 0.4
# Number of replications
rep = 100
# Level for variance of inactive points in components
var_level = 0.01
# Random seeds
seed1 = 159963
seed2 = 753741
seed3 = 852367
seed4 = 951753
seed5 = 364192
# number of time points in each simulation
my_nTR = 120
# AR structure in each simulation
my_phi = 0.47
# Convergence criterion
my_eps = 1e-6
# Maximum number of iterations
my_maxit = 1000
# Whether it's noisy ICA
my_noisyICA = TRUE
# group components to be used
my_components = c(13,23,30)
# the list of candidate nu
nu_list = seq(0.1,3,0.1)


#######################################################################################################
#generate data for matlab
#split into 5 parts

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed1)
for (i in 1:20) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i]] = simData$X
  smat_list[[i]] = simData$S
  mmat_list[[i]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_nonsparse_1_20.RData")
writeMat("../../Data/highdim_nonsparse_1_20.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed2)
for (i in 21:40) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-20]] = simData$X
  smat_list[[i-20]] = simData$S
  mmat_list[[i-20]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_nonsparse_21_40.RData")
writeMat("../../Data/highdim_nonsparse_21_40.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed3)
for (i in 41:60) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-40]] = simData$X
  smat_list[[i-40]] = simData$S
  mmat_list[[i-40]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_nonsparse_41_60.RData")
writeMat("../../Data/highdim_nonsparse_41_60.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed4)
for (i in 61:80) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-60]] = simData$X
  smat_list[[i-60]] = simData$S
  mmat_list[[i-60]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_nonsparse_61_80.RData")
writeMat("../../Data/highdim_nonsparse_61_80.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed5)
for (i in 81:100) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-80]] = simData$X
  smat_list[[i-80]] = simData$S
  mmat_list[[i-80]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_nonsparse_81_100.RData")
writeMat("../../Data/highdim_nonsparse_81_100.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)
