% run the sparseICA-EBM function on the simulated real images

%sparse  ICA  sim

load('../../Data/highdim_sparse_1_20.mat');

for i = 1:20
    xmat = xmat_matlab_lowSNR(:,(1+120*(i-1)):(120+120*(i-1)));
    tStart = tic; 
    [W,totalIterSparse,Cost,independence,sparsity] = ICA_EBM_Sparse(xmat',0.01,0.01);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/SICA_EBM/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','xmat','tEnd');
end




