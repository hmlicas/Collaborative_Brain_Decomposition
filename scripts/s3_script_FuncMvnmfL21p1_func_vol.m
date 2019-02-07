clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-data/vol_sbjLst.txt';
maskFile = 'path-to-data/mask.nii.gz';
prepDataFile = 'path-to-output-vol/prepData.mat';
outDir = 'path-to-output-vol/res';

resId = 'fmri_vol';
initName = 'path-to-output/robustInit/init.mat';
K = 17;
alphaS21 = 2;
alphaL = 10;
vxI = 0;
spaR = 1;
ard = 1;
eta = 1;
iterNum = 30;
calcGrp = 1;
parforOn = 0;

deployFuncMvnmfL21p1_func_vol(sbjListFile,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);
