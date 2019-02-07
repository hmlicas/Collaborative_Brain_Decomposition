clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-data/vol_sbjLst.txt';
maskFile = 'path-to-data/mask.nii.gz';
prepDataFile = 'path-to-output-vol/prepData.mat';
outDir = 'path-to-output-vol/init';
spaR = 1;
vxI = 0;
ard = 1;
iterNum = 1000;
K = 17;
tNum = 118;
alpha = 2;
beta = 10;
resId = 'fmri_vol';

deployFuncInit_vol(sbjListFile,maskFile,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId);
