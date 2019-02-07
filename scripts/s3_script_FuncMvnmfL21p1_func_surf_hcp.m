clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-hcp-data/hcp_sbjLst.txt';
wbPath = 'path-to-workbench/wb_command.exe';
prepDataFile = 'path-to-output-hcp/prepData.mat';
outDir = 'path-to-output-hcp/res';

resId = 'hcp';
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

deployFuncMvnmfL21p1_func_surf_hcp(sbjListFile,wbPath,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);
