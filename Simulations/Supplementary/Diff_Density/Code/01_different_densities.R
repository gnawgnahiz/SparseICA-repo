###################################################################
# R codes for evaluating the performance of Sparse ICA at different source densities
###################################################################
rm(list = ls())

library(SparseICA)
library(ProDenICA)
library(steadyICA)
library(fastICA)

source("00_utils.R")

################################################################################################
### Can use letters a-r below for dist
jordan_simu = function(p=2,N=1024,dist="a"){
  # generate mixing matrix
  A0 = mixmat(p)
  # generate source matrix
  s = rjordan(dist,N)
  if (p>1){
    for (i in 1:(p-1)) {
      s = cbind(s,rjordan(dist,N))
    }
  }
  s = scale(s)
  # generate x
  x = s %*% A0
  return(list(x=as.matrix(x),s=as.matrix(s),m=as.matrix(A0)))
}


########################################################
# 8 example densities:
# 1st row: super-gaussian distributions: laplace, exponential, t with df=3, t with df=5
# 2ed row: sub-gaussian distributions

dist_list = c("a","b","d","e",
              "c","f","k","n")
rep = 100
dat_sparse = matrix(nrow = length(dist_list),ncol = rep)
rownames(dat_sparse) = dist_list
dat_fast = matrix(nrow = length(dist_list),ncol = rep)
rownames(dat_fast) = dist_list
dat_info = matrix(nrow = length(dist_list),ncol = rep)
rownames(dat_info) = dist_list

set.seed(2580)
for (i in 1:length(dist_list)) {
  for (j in 1:rep) {
    dat = jordan_simu(p=2,N=1024,dist = dist_list[i])
    
    my_sparse = sparseICA_Rcpp(dat$x,n.comp = 2,nu=1,restarts.pbyd = 1)
    estM_sparse = est.M.ols(my_sparse$estS,dat$x)
    dat_sparse[i,j] =frobICA(M1=estM_sparse,M2=dat$m,standardize = T)
    
    my_fast = fastICA(dat$x,2)
    dat_fast[i,j] = frobICA(M1=my_fast$A,M2=dat$m,standardize = T)
    
    my_infomax = infomaxICA(dat$x,2)
    dat_info[i,j] =frobICA(M1=my_infomax$M,M2=dat$m,standardize = T)
  }
}
save(dat_sparse,dat_fast,dat_info,file = "../Results/PRMSE_jorson_8version.RData")

#load("../Results/PRMSE_jorson.RData")
xpos = seq(-5, 5, by = 0.01) 
pdf(file="../Figure/Supp_Figure1.pdf",width = 12,height = 12)
par(mfrow=c(4,4))
plot(djordan("a",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="a")
temp = djordan("b",xpos)
plot(temp[temp!=0],type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="b")
plot(djordan("d",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="c")
plot(djordan("e",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="d")
boxplot(dat_sparse[1,],dat_fast[1,],dat_info[1,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[2,],dat_fast[2,],dat_info[2,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[3,],dat_fast[3,],dat_info[3,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[4,],dat_fast[4,],dat_info[4,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
plot(djordan("c",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="e")
plot(djordan("f",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="f")
plot(djordan("k",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="g")
plot(djordan("n",xpos),type = "l",col = "darkred",lwd=2,xaxt="n",yaxt="n",xlab = "",ylab = "",main="h")
boxplot(dat_sparse[5,],dat_fast[5,],dat_info[5,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[6,],dat_fast[6,],dat_info[6,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[7,],dat_fast[7,],dat_info[7,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
boxplot(dat_sparse[8,],dat_fast[8,],dat_info[8,],names=c("SparseICA","FastICA","InfomaxICA"),ylab="PRMSE")
par(mfrow=c(1,1))
dev.off()
