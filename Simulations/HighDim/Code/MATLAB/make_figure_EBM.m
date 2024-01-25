% test the sparseICA-EBM function

%sparse  ICA  sim 2
%addpath("/Users/ben/Desktop/work2/DICA/2stageNEW/MIToolbox")
%addpath("/Users/zihang/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Ben_RA/SparseICA/Simulations/SparseICA_EBM");
load('simreal.mat');

tStart = tic; 
[W,totalIterSparse,Cost,independence,sparsity] = ICA_EBM_Sparse(xmat',0.01,0.01);
myS = W * xmat';
tEnd = toc(tStart);
save('test_EBM_real.mat','myS','W','xmat','tEnd');








