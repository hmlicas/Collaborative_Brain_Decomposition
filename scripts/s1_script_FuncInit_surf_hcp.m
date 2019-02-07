clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-hcp-data/hcp_sbjLst.txt';
wbPath = 'path-to-workbench/wb_command.exe';
surfL = 'path-to-hcp-data/100307/MNINonLinear/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii';
surfR = 'path-to-hcp-data/100307/MNINonLinear/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii';
surfML = 'path-to-hcp-data/100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii';
surfMR = 'path-to-hcp-data/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii';
prepDataFile = 'path-to-output-hcp/prepData.mat';
outDir = 'path-to-output-hcp/init';
spaR = 1;
vxI = 0;
ard = 0;
iterNum = 1000;
K = 17;
tNum = 1200;
alpha = 2;
beta = 10;
resId = 'hcp';

deployFuncInit_surf_hcp(sbjListFile,wbPath,surfL,surfR,surfML,surfMR,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId); %,svFileL,svFileR);
