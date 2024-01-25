###################################################################
# R codes for making BIC selection plots
# 01: making BIC plots
###################################################################

rm(list = ls())

library(SparseICA)

source("../../../00_utils.R")

#################################################################################################
# make BIC plot for single-subject spatio-temporal simulations
set.seed(2123)
simData = SimFMRI123(var.inactive = 0,noisyICA = T, snr=0.4,nTR = 50,phi = 0.47)
xmat = simData$X
smat = simData$S
mmat = simData$Ms

pdf("../Figures/Supp_Figure15.pdf",width = 8,height = 4)
select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec", eps = 1e-6,
                                     maxit = 500,BIC_plot = T,verbose=F,nu_list = seq(0.1,4,0.1))
my_nu = select_sparseICA$best_nu
dev.off()

###########################################################################################
# make BIC plot for single-subject high-dimensional simulations
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench') 
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
#ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench') 

# load true components data
load("../Data/group_SparseICA.RData")
load("../Data/template_xifti.RData")
load("../Data/template_timecourses2.RData")

my_SNR = 0.4
rep = 100
var_level = 0
my_nTR = 120
my_phi = 0.47
my_eps = 1e-6
my_maxit = 1000
my_noisyICA = TRUE
my_components = c(17,6,23)

set.seed(2123)
simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, 
                      var.inactive = var_level, components = my_components, 
                      true_data = dat_iter5,template_xifti = template_image,template_time = sample_time)
xmat = simData$X
smat = simData$S
mmat = simData$Ms

nu_list = seq(0.1,4,0.1)
pdf("../Figures/Supp_Figure16.pdf",width = 8,height = 4)
select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec",restarts = 30, 
                                eps = my_eps,maxit = my_maxit, BIC_plot = F,nu_list = nu_list)
my_nu = select_sparseICA$best_nu
dev.off()

##########################################################################################
# make plot for BIC selection in group simulations
source("00_utils_group.R")

simData = list()
for (i in 1:nsub){
  simData[[i]] = SimFMRI.ngca(snr = SNR, noisyICA=noisyICA,
                              gamma.rate = gamma.rate, gamma.shape = gamma.shape, 
                              FWHM = FWHM,var.inactive = var.inactive)
}

pdf("../Figures/Supp_Figure17.pdf",width = 8,height = 4)
bic = BIC_group_sparse_ICA(simData,n.group = 3,BIC_plot  = T)
dev.off()

###########################################################################################
# make BIC plot for real data application of ABIDE

load("../Data/PC30_subPC85.RData")
cat("PC30 loaded!\n")

nu_selection2 = BIC_sparseICA(xData = PC30, n.comp = 30,whiten = "none", restarts = 40,
                              eps = 1e-6,maxit = 500,BIC_plot = F,nu_list = seq(0.1,4,0.05))

pdf("../Figures/Supp_Figure18.pdf",width = 8,height = 4)
plot(nu_selection2$nu_list,nu_selection2$BIC, main = "BIC Plot", xlab = "nu", ylab = "BIC", type = "l")
abline(v=nu_selection2$best_nu,col="red")
dev.off()
