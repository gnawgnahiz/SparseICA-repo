%run the SparseICA-EBM function on the simulated 123 images

%sparse  ICA  sim

load('../../Data/sim123_nonsparse_lowSNR.mat');

for i = 1:100
    xmat = xmat_matlab_lowSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    tStart = tic; 
    [W,totalIterSparse,Cost,independence,sparsity] = ICA_EBM_Sparse(xmat',0.01,0.01);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/nonsparse/EBM/low/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','xmat','tEnd');
end

load('../../Data/sim123_nonsparse_mediumSNR.mat');

for i = 1:100
    xmat = xmat_matlab_mediumSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    tStart = tic; 
    [W,totalIterSparse,Cost,independence,sparsity] = ICA_EBM_Sparse(xmat',0.01,0.01);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/nonsparse/EBM/medium/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','xmat','tEnd');
end

load('../../Data/sim123_nonsparse_highSNR.mat');

for i = 1:100
    xmat = xmat_matlab_highSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    tStart = tic; 
    [W,totalIterSparse,Cost,independence,sparsity] = ICA_EBM_Sparse(xmat',0.01,0.01);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/nonsparse/EBM/high/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','xmat','tEnd');
end

