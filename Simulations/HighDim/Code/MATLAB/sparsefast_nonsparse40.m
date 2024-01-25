% test the sparsefastica function

%sparse fastICA  sim

approach = 'symm';
numOfIC = 3;
finetune = 'off';
myy = 1;
stabilization = 'on';
epsilon1 = 0.001;
epsilon2 = 0.000001;
maxNumIterations = 1000;
maxFinetune = 100;
initState = 'rand';
guess = 0;
sampleSize = 1;
displayMode = 'off';
displayInterval = 1;
s_verbose = 'on';
sigma = 1;
rou = 0.965;
func = 'tanh';
round_2nd = 20;

% low SNR situation
load('../../Data/highdim_nonsparse_21_40.mat')

for i = 1:20
    xmat = xmat_matlab_lowSNR(:,(1+120*(i-1)):(120+120*(i-1)));
    [E, D] = pcamat(xmat',1,3);
    [nv, wm, dwm] = whitenv(xmat', E, D);

    X = nv;
    whiteningMatrix = wm;
    dewhiteningMatrix = dwm;

    tStart = tic; 
    [A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, finetune, myy, stabilization, epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, guess, sampleSize, displayMode, displayInterval, s_verbose, sigma, rou, func, round_2nd);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/nonsparse/sparsefast/estS_' num2str(i+20) '.mat' ];
    save(filenm,'myS','W','A','tEnd','xmat');
end








