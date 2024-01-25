###################################################################
# R codes for real data application
# 03: Select the tuning parameter in group Sparse ICA
# Note: This script was run on a personal computer.
###################################################################

rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(ciftiTools)
#ciftiTools.setOption('wb_path', '/usr/local/workbench')
ciftiTools.setOption('wb_path', '/Applications/workbench') 
# ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
# ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench')
source("../../Simulations/00_utils.R")

load("../Data/PC30_subPC85.RData")
cat("PC30 loaded!\n")

nu_selection = BIC_sparseICA(xData = PC30, n.comp = 30,whiten = "none",method = "C", 
                              irlba = TRUE,eps = 1e-6,maxit = 500,BIC_plot = T,nu_list = seq(0.1,4,0.05))

save(nu_selection,file = "../Results/nu_selection_groupPC_BIC.RData")


