% test the sparsefastica function

%sparse fastICA  sim 2
%addpath("/Users/ben/Desktop/work2/DICA/2stageNEW/MIToolbox")
%addpath("/Users/zihang/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Ben_RA/SparseICA/Simulations/sparsefastica")
%load('sim123.mat')
load('simreal.mat')

xmatt=xmat';
[E, D] = pcamat(xmatt,1,3);
[nv, wm, dwm] = whitenv(xmatt, E, D);

X = nv;
whiteningMatrix = wm;
dewhiteningMatrix = dwm;
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

tStart = tic; 
[A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, finetune, myy, stabilization, epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, guess, sampleSize, displayMode, displayInterval, s_verbose, sigma, rou, func, round_2nd);

%[A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, 3,finetune,myy,stabilization);
%[icasig, A, W] = fastica(xmatt,'approach', 'symm', 'g', 'tanh','lastEig', 3,'numOfIC', 3);
tEnd = toc(tStart);
myS = W * xmatt;
save('test_sparsefastica_real.mat','myS','W','A','xmat');








