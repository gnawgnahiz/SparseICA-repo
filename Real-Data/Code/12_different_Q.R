###################################################################
# R codes for generating sparse ICs with Q = 10, 20, 30, 40, 50
# Note: This script was run on a personal computer
###################################################################
rm(list = ls())

library(steadyICA)
library(fastICA)
library(irlba)
library(SparseICA)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench') 
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
#ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench') 

###############################################################################################
# Q=10
load("../Data/PC10_subPC85.RData")
best_nu = 2
Q = 10
my_sparseICA = sparseICA(xData = PC10, n.comp = Q, restarts = 40,irlba = T, whiten = "none",
                         nu = best_nu,eps = 1e-4,maxit = 500, converge_plot = F)
# sign change
a=signchange(t(my_sparseICA$estS))
my_sparseICA$estS_sign=t(a$S)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_sparseICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/sparse_Q10_nu2")
#save(my_sparseICA,file = "../Results/group_sparseICA_Q10_nu1.RData")

my_fastICA = fastICA(PC10, n.comp = Q,method = "C",maxit = 500,tol = 1e-6,verbose = T)
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA$estS_sign)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/fast_Q10")

###############################################################################################
# Q=20
load("../Data/PC20_subPC85.RData")
best_nu = 2
Q = 20
my_sparseICA = sparseICA(xData = PC10, n.comp = Q, restarts = 40,irlba = T, whiten = "none",
                         nu = best_nu,eps = 1e-4,maxit = 500, converge_plot = F)
# sign change
a=signchange(t(my_sparseICA$estS))
my_sparseICA$estS_sign=t(a$S)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_sparseICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/sparse_Q20_nu2")
#save(my_sparseICA,file = "../Results/group_sparseICA_Q10_nu1.RData")

my_fastICA = fastICA(PC20, n.comp = Q,method = "C",maxit = 500,tol = 1e-6,verbose = T)
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA$estS_sign)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/fast_Q20")


###############################################################################################
# Q=30
load("../Data/PC30_subPC85.RData")
best_nu = 2
Q = 30
my_sparseICA = sparseICA(xData = PC10, n.comp = Q, restarts = 40,irlba = T, whiten = "none",
                         nu = best_nu,eps = 1e-4,maxit = 500, converge_plot = F)
# sign change
a=signchange(t(my_sparseICA$estS))
my_sparseICA$estS_sign=t(a$S)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_sparseICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/sparse_Q30_nu2")
#save(my_sparseICA,file = "../Results/group_sparseICA_Q10_nu1.RData")

my_fastICA = fastICA(PC30, n.comp = Q,method = "C",maxit = 500,tol = 1e-6,verbose = T)
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA$estS_sign)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/fast_Q30")


###############################################################################################
# Q=40
load("../Data/PC40_subPC85.RData")
best_nu = 2
Q = 40
my_sparseICA = sparseICA(xData = PC10, n.comp = Q, restarts = 40,irlba = T, whiten = "none",
                         nu = best_nu,eps = 1e-4,maxit = 500, converge_plot = F)
# sign change
a=signchange(t(my_sparseICA$estS))
my_sparseICA$estS_sign=t(a$S)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_sparseICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/sparse_Q40_nu2")
#save(my_sparseICA,file = "../Results/group_sparseICA_Q10_nu1.RData")

my_fastICA = fastICA(PC40, n.comp = Q,method = "C",maxit = 500,tol = 1e-6,verbose = T)
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA$estS_sign)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/fast_Q40")


###############################################################################################
# Q=50
load("../Data/PC50_subPC85.RData")
best_nu = 2
Q = 50
my_sparseICA = sparseICA(xData = PC10, n.comp = Q, restarts = 40,irlba = T, whiten = "none",
                         nu = best_nu,eps = 1e-4,maxit = 500, converge_plot = F)
# sign change
a=signchange(t(my_sparseICA$estS))
my_sparseICA$estS_sign=t(a$S)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_sparseICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/sparse_Q50_nu2")
#save(my_sparseICA,file = "../Results/group_sparseICA_Q10_nu1.RData")

my_fastICA = fastICA(PC50, n.comp = Q,method = "C",maxit = 500,tol = 1e-6,verbose = T)
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA$estS_sign)

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")

cifti_image$data$cortex_left=my_fastICA$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = Q)

write_cifti(cifti_image,"../Results/fast_Q50")
