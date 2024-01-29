###################################################################
# R codes for MDS analysis of group ICA
# Note: This script was run on a personal computer
###################################################################
rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(ciftiTools)

ciftiTools.setOption('wb_path', '/Applications/workbench') 
# ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
# ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench')

load("../Data/PC30_subPC85.RData")

################################################################################
# visualize multiple initializations

# generate the random restarts matrices for all methods
my_W_init = gen.inits(p=30,d=30,runs = 40,orth.method="svd")

################################################################################
# sparse ICA
best_nu=2
sparseICA_restarts=list()
for (k in 1:40) {
  my_sparseICA=sparseICA(xData = PC30, n.comp = 30, restarts = 1, irlba = T, U.list = list(my_W_init[[k]]),
                         whiten = "none",nu = best_nu,eps = 1e-5,maxit = 500, converge_plot = F)
  sparseICA_restarts[[k]]=my_sparseICA$estS
}

# distance of S
dis_S_sparse = matrix(nrow = 40,ncol = 40)
for (i in 1:40) {
  for (j in 1:40) {
    dis_S_sparse[i,j]=sqrt(frobICA(S1=sparseICA_restarts[[i]],S2=sparseICA_restarts[[j]],standardize = TRUE))
    cat(i," th init and ",j," th init finish!\n")
  }
}
mds_S_sparse=cmdscale(dis_S_sparse,k=2)
save(dis_S_sparse,mds_S_sparse,file = "../Results/MDS_S_sparse.RData")

################################################################################
# infomax ICA
infomaxICA_restarts=list()
for (k in 1:40) {
  my_infomaxICA=infomaxICA(X=PC30,n.comp = 30,whiten = FALSE,eps = 1e-5,maxit = 500,W.list = list(my_W_init[[k]]),verbose = T)
  infomaxICA_restarts[[k]]=my_infomaxICA$S
}

# distance of S
dis_S_infomax = matrix(nrow = 40,ncol = 40)
for (i in 1:40) {
  for (j in 1:40) {
    dis_S_infomax[i,j]=sqrt(frobICA(S1=infomaxICA_restarts[[i]],S2=infomaxICA_restarts[[j]],standardize = TRUE))
    cat(i," th init and ",j," th init finish!\n")
  }
}
mds_S_infomax=cmdscale(dis_S_infomax,k=2)

save(dis_S_infomax,mds_S_infomax,file = "../Results/MDS_S_infomax.RData")


################################################################################
# fast ICA
fastICA_restarts=list()
for (k in 1:40) {
  my_fastICA=fastICA(PC30,n.comp = 30,maxit = 500,tol=1e-5,method = "C",w.init = my_W_init[[k]])
  fastICA_restarts[[k]]=my_fastICA$S
}

# distance of S
dis_S_fast = matrix(nrow = 40,ncol = 40)
for (i in 1:40) {
  for (j in 1:40) {
    dis_S_fast[i,j]=sqrt(frobICA(S1=fastICA_restarts[[i]],S2=fastICA_restarts[[j]],standardize = TRUE))
    cat(i," th init and ",j," th init finish!\n")
  }
}
mds_S_fast=cmdscale(dis_S_fast,k=2)

save(dis_S_fast,mds_S_fast,file = "../Results/MDS_S_fast.RData")

load("../Results/MDS_S_sparse.RData")
load("../Results/MDS_S_infomax.RData")
load("../Results/MDS_S_fast.RData")

pdf("../Figures/Supp_Figure14.pdf",width = 9,height = 3)
par(mfrow=c(1,3))
plot(mds_S_sparse[,1],mds_S_sparse[,2],xlab = "x",ylab = "y",main = "Sparse ICA",xlim = c(-0.35,0.35),ylim = c(-0.35,0.35))
plot(mds_S_infomax[,1],mds_S_infomax[,2],xlab = "x",ylab = "y",main = "Infomax ICA",xlim = c(-0.35,0.35),ylim = c(-0.35,0.35))
plot(mds_S_fast[,1],mds_S_fast[,2],xlab = "x",ylab = "y",main = "Fast ICA",xlim = c(-0.35,0.35),ylim = c(-0.35,0.35))
par(mfrow=c(1,1))
dev.off()
